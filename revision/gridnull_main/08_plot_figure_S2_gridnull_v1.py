#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
08_plot_figure_S2_gridnull_v1.py

Publication-quality segmentation figure (Figure S2) from the grid-null
timelines produced by 03_Micro_Analysis_NOLPM_gridnull.py.

This is the grid-null replacement for 08_plot_figure_S1.py (which read
group_mean_zS.csv from the legacy 01_segment_boundaries.py pipeline).
It reproduces the same visual design (Okabe-Ito palette, thin per-condition
group-mean curves, thick pooled curve, dashed adopted period boundaries,
shaded P1/P2/P3 bands) but computes the curves from the grid-null
timelines_all_master.csv files, using exactly the same pooling rule as
08_changepoint_gridnull.py (a time point is kept if at least --min_frac of
sources have a valid value; pooled mean over all pair columns of both
conditions).

If --group_curve is given (changepoint_group_curve.csv from
08_changepoint_gridnull.py), the pooled curve computed here is verified
against it (max |diff| must be < 1e-9 on the shared time points); the script
aborts if they disagree, so the figure is guaranteed to show the same curve
that the change-point estimate (c1 = 8.5 s, c2 = 20.0 s at delta = 0.1) was
computed from.

Usage (local, from ROOT/gridnull_main, conda env nieh_iscience_2026):
  python 08_plot_figure_S2_gridnull_v1.py \
      --timelines "out_grid_visible_d0.1/timelines_all_master.csv,out_grid_invisible_d0.1/timelines_all_master.csv" \
      --labels "visible,invisible" \
      --group_curve out_changepoint_grid/changepoint_group_curve.csv \
      --boundaries 20 60 --delta 0.1 --outdir fig_S2_gridnull

Options:
  --timelines CSV1,CSV2   comma-separated timelines_all_master.csv paths
                          (one per condition, same order as --labels)
  --labels A,B            condition labels (default: visible,invisible)
  --metric NAME           column prefix (default: z_S_t)
  --min_frac F            min fraction of valid sources per time point
                          (default: 0.5; must match 08_changepoint_gridnull.py)
  --group_curve CSV       optional changepoint_group_curve.csv for verification
  --boundaries B1 B2      adopted period boundaries in seconds (default: 20 60)
  --tmax T                right edge of the x-axis in seconds (default: 180)
  --sem                   add +/- SEM bands around the condition curves
  --delta LABEL           delta label shown in the corner (default: none)
  --formats ...           output formats (default: pdf png)
  --outdir DIR            output directory (default: fig_S2_gridnull)

Outputs:
  <outdir>/figure_S2_gridnull.<fmt>   the figure
  <outdir>/figure_S2_source_data.csv  t_center, per-condition mean (+SEM),
                                      pooled mean, n_valid (source data)
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# colorblind-safe palette (Okabe-Ito), identical to 08_plot_figure_S1.py
COLORS = {"visible": "#0072B2", "invisible": "#D55E00"}
FALLBACK_COLORS = ["#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
POOLED_COLOR = "#000000"
PERIOD_SHADES = ["#F2F2F2", "#FFFFFF", "#F2F2F2"]


def load_condition(path, metric):
    df = pd.read_csv(path)
    cols = [c for c in df.columns if c.startswith(metric + "__")]
    if not cols:
        raise SystemExit(f"no {metric}__* columns in {path}")
    return df[["t_center"] + cols].set_index("t_center")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timelines", type=str, required=True)
    ap.add_argument("--labels", type=str, default="visible,invisible")
    ap.add_argument("--metric", type=str, default="z_S_t")
    ap.add_argument("--min_frac", type=float, default=0.5)
    ap.add_argument("--group_curve", type=str, default=None)
    ap.add_argument("--boundaries", type=float, nargs=2, default=[20.0, 60.0])
    ap.add_argument("--tmax", type=float, default=180.0)
    ap.add_argument("--sem", action="store_true")
    ap.add_argument("--delta", type=str, default=None)
    ap.add_argument("--formats", nargs="+", default=["pdf", "png"])
    ap.add_argument("--outdir", type=str, default="fig_S2_gridnull")
    args = ap.parse_args()

    paths = [p.strip() for p in args.timelines.split(",") if p.strip()]
    labels = [s.strip() for s in args.labels.split(",") if s.strip()]
    if len(paths) != len(labels):
        raise SystemExit("--timelines and --labels must have the same length")

    conds = {lab: load_condition(p, args.metric) for lab, p in zip(labels, paths)}

    # pooled wide table across all conditions (same construction as
    # 08_changepoint_gridnull.py: outer join on t_center, all pair columns)
    frames = []
    for k, (lab, sub) in enumerate(conds.items()):
        s = sub.copy()
        s.columns = [f"src{k}__{c}" for c in s.columns]
        frames.append(s)
    wide = pd.concat(frames, axis=1).sort_index()
    n_valid = wide.notna().sum(axis=1)
    coverage = n_valid / wide.shape[1]
    keep = coverage >= args.min_frac
    pooled = wide.mean(axis=1, skipna=True)

    # optional verification against the change-point pipeline's own curve
    if args.group_curve:
        gc = pd.read_csv(args.group_curve).set_index("t_center")
        common = wide.index.intersection(gc.index)
        if len(common) == 0:
            raise SystemExit("group_curve has no overlapping t_center values")
        diff = np.nanmax(np.abs(pooled.loc[common].values
                                - gc.loc[common, "mean_zS"].values))
        if not (diff < 1e-9):
            raise SystemExit(
                f"pooled curve does not match {args.group_curve} "
                f"(max |diff| = {diff:.3e}); check --min_frac/--metric/inputs")
        print(f"[check] pooled curve matches {args.group_curve} "
              f"(max |diff| = {diff:.3e} over {len(common)} time points)")

    os.makedirs(args.outdir, exist_ok=True)

    # source data CSV
    out = pd.DataFrame({"t_center": wide.index, "n_valid": n_valid.values,
                        "coverage": coverage.values,
                        "mean_zS_pooled": pooled.values,
                        "kept_for_changepoint": keep.values})
    for lab, sub in conds.items():
        out[f"mean_zS_{lab}"] = sub.mean(axis=1, skipna=True).reindex(wide.index).values
        out[f"sem_zS_{lab}"] = (sub.std(axis=1, skipna=True)
                                / np.sqrt(sub.notna().sum(axis=1))
                                ).reindex(wide.index).values
    src_csv = os.path.join(args.outdir, "figure_S2_source_data.csv")
    out.to_csv(src_csv, index=False)

    # figure
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    b1, b2 = sorted(args.boundaries)
    edges = [0.0, b1, b2, args.tmax]
    names = ["P1", "P2", "P3"]
    for i in range(3):
        ax.axvspan(edges[i], edges[i + 1], color=PERIOD_SHADES[i], zorder=0)
        ax.text((edges[i] + edges[i + 1]) / 2.0, 0.97, names[i],
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=10, color="#666666")
    fb = iter(FALLBACK_COLORS)
    for lab, sub in conds.items():
        m = sub.mean(axis=1, skipna=True)
        col = COLORS.get(lab, next(fb))
        tt, mm = m.index.values, m.values
        ax.plot(tt, mm, color=col, lw=1.2, label=lab, zorder=3)
        if args.sem:
            sem = (sub.std(axis=1, skipna=True)
                   / np.sqrt(sub.notna().sum(axis=1))).values
            ax.fill_between(tt, mm - sem, mm + sem, color=col, alpha=0.18,
                            lw=0, zorder=2)
    tk, yk = wide.index.values[keep.values], pooled.values[keep.values]
    ax.plot(tk, yk, color=POOLED_COLOR, lw=2.2,
            label="pooled (change-point input)", zorder=4)
    for b in (b1, b2):
        ax.axvline(b, ls="--", color="#444444", lw=1.0, zorder=5)
    ax.axhline(0.0, color="#999999", lw=0.6, zorder=1)
    if args.delta:
        ax.text(0.99, 0.02, f"$\\delta$ = {args.delta} s",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                color="#444444")
    ax.set_xlim(0.0, args.tmax)
    ax.set_xlabel("Time in session (s)")
    ax.set_ylabel("Group-mean $z_S(t)$ (grid null)")
    ax.legend(frameon=False, fontsize=9, loc="upper right",
              bbox_to_anchor=(1.0, 0.90))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    for fmt in args.formats:
        fp = os.path.join(args.outdir, f"figure_S2_gridnull.{fmt}")
        fig.savefig(fp, dpi=300 if fmt == "png" else None)
        print(f"[out] {fp}")
    print(f"[out] {src_csv}")


if __name__ == "__main__":
    sys.exit(main())
