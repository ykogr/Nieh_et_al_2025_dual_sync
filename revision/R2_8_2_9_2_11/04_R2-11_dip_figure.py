#!/usr/bin/env python3
"""
04_R2-11_dip_figure.py
Supplementary figure for R2-11: dequantized (jittered) distributions of the
synchronization time lags (dt) with Hartigan's dip test results.

The figure shows, for each tolerance delta = 0.1, 0.2, 0.3 s, density
histograms of one representative jitter replicate (uniform jitter on
(-0.05, +0.05) s added to the raw quantized dt values; same seed family as
03_R2-11_dip_test.py) for the visible and invisible conditions, annotated
with the median dip statistic and median p-value across the jitter
replicates taken from dip_test_jitter_summary.csv (output of
03_R2-11_dip_test.py). Reading the statistics from that file (instead of
recomputing) guarantees that the figure annotations match the response
letter and dip_test_summary.txt exactly.

Requirements
------------
  conda activate nieh_iscience_2026   # numpy, pandas, matplotlib

Usage
-----
  python 04_R2-11_dip_figure.py --dtcsv dt.csv \
      --jitter-summary dip_test_out/dip_test_jitter_summary.csv \
      --outdir dip_figure_out [--seed 20260815] [--bin-width 0.01] \
      [--descending-x]

Outputs (written to --outdir)
-----------------------------
  fig_S_dip_dequantized.pdf : vector figure for the supplement
  fig_S_dip_dequantized.png : 300-dpi raster preview
  dip_figure_stats.txt      : the annotated statistics, for cross-checking
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESOLUTION = 0.1  # annotation grid (s)
COL_VIS = "#E2711D"   # orange, matches Figure 2
COL_INVIS = "#2C7FB8" # blue,   matches Figure 2
DELTAS = ["0.1", "0.2", "0.3"]


def find_col(cols, delta, cond):
    """Find the dt.csv column for a given delta ('0.1') and condition."""
    for c in cols:
        name = str(c)
        if not name.startswith(delta):
            continue
        if cond == "invis" and "invis" in name:
            return c
        if cond == "vis" and "vis" in name and "invis" not in name:
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dtcsv", required=True)
    ap.add_argument("--jitter-summary", required=True,
                    help="dip_test_jitter_summary.csv from "
                         "03_R2-11_dip_test.py")
    ap.add_argument("--outdir", default="dip_figure_out")
    ap.add_argument("--seed", type=int, default=20260815,
                    help="seed for the representative jitter replicate")
    ap.add_argument("--bin-width", type=float, default=0.01)
    ap.add_argument("--descending-x", action="store_true",
                    help="orient the x-axis positive-to-negative, as in "
                         "main Figure 2 (default: conventional ascending)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.dtcsv)
    summ = pd.read_csv(args.jitter_summary).set_index("column")
    rng = np.random.default_rng(args.seed)
    half = RESOLUTION / 2.0

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.2), sharey=False)
    stats_lines = ["R2-11 supplementary figure: annotated statistics",
                   f"(from {args.jitter_summary}; jitter U(-{half}, +{half})"
                   f" s; representative replicate seed {args.seed})", ""]

    for k, (ax, delta) in enumerate(zip(axes, DELTAS)):
        dmax = float(delta)
        edges = np.arange(-(dmax + half), dmax + half + 1e-9, args.bin_width)
        for cond, color, label in [("vis", COL_VIS, "vis"),
                                   ("invis", COL_INVIS, "invis")]:
            col = find_col(df.columns, delta, cond)
            if col is None:
                sys.exit(f"Column for delta={delta}, {cond} not found "
                         f"in {args.dtcsv}")
            if col not in summ.index:
                sys.exit(f"Column '{col}' not found in {args.jitter_summary}")
            x = pd.to_numeric(df[col], errors="coerce").dropna().values
            xj = x + rng.uniform(-half, half, size=len(x))
            ax.hist(xj, bins=edges, density=True, histtype="stepfilled",
                    color=color, alpha=0.45, edgecolor=color, linewidth=1.0,
                    label=label)
            s = summ.loc[col]
            stats_lines.append(
                f"{col}: n = {int(s['n'])}, median Dip = "
                f"{s['dip_median']:.4f}, median p = {s['p_median']:.2f}, "
                f"prop(p < 0.05) = {s['prop_p_below_0.05']:.4f}")

        s_vis = summ.loc[find_col(df.columns, delta, "vis")]
        s_inv = summ.loc[find_col(df.columns, delta, "invis")]
        anno = (f"vis:    Dip = {s_vis['dip_median']:.3f}, "
                f"p = {s_vis['p_median']:.2f}\n"
                f"invis: Dip = {s_inv['dip_median']:.3f}, "
                f"p = {s_inv['p_median']:.2f}")
        ax.text(0.03, 0.97, anno, transform=ax.transAxes, fontsize=8.5,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                          edgecolor="#B0B0B0", linewidth=0.8))

        ax.set_title(rf"$\delta$ = {delta}", fontsize=12)
        ax.set_xlabel(r"Time lag $\Delta t$ [s]", fontsize=11)
        if k == 0:
            ax.set_ylabel("Density", fontsize=11)
        ax.text(-0.08, 1.12, f"({chr(97 + k)})", transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="top")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(-(dmax + half), dmax + half)
        # extra headroom so the annotation box and legend sit clear of the bars
        ax.autoscale(axis="y")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.40)
        if args.descending_x:
            ax.invert_xaxis()

    axes[-1].legend(loc="upper right", frameon=False, fontsize=10,
                    title=None)
    fig.suptitle("")
    fig.tight_layout()

    pdf = os.path.join(args.outdir, "fig_S_dip_dequantized.pdf")
    png = os.path.join(args.outdir, "fig_S_dip_dequantized.png")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    txt = "\n".join(stats_lines)
    with open(os.path.join(args.outdir, "dip_figure_stats.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    print(f"\nwritten: {pdf}\nwritten: {png}")


if __name__ == "__main__":
    main()
