#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
10_plot_figure2_timelags.py  (trialfix_kit v4.4.1)

Redraw manuscript Figure 2 (distribution of synchronization time lags dt for
delta = 0.1, 0.2, 0.3 s; visible vs invisible) from dt_trialfix.csv, the wide
CSV written by make_dt_csv.py (columns "<delta>_vis", "<delta>_invis";
dt = t_A - t_B of every NLOPM-matched event, already rounded to the 0.1-s grid).

Layout follows the 2026-08-23 manuscript figure: one panel per delta, grouped
bars (proportion of matched events at each admissible lag), x axis from +delta
down to -delta in 0.1-s steps, visible = orange, invisible = blue, legend to the
right of the last panel, panel letters (a)-(c) and "delta = ..." above each panel.

Usage (called by run_gridnull_main_trialfix.sh, STAGE=fig2 / all / downstream):
    python 10_plot_figure2_timelags.py --dtcsv dt_trialfix.csv --outdir fig2_timelags
Outputs:
    fig2_timelags/timelags.pdf / .png           -> Overleaf figs/timelags.pdf
    fig2_timelags/timelags_source_data.csv      proportions and counts per bar
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL = {"vis": "#E36C09", "invis": "#0070C0"}
LABEL = {"vis": "vis", "invis": "invis"}


def lag_grid(delta):
    k = int(round(delta / 0.1))
    return [round(0.1 * i, 1) for i in range(k, -k - 1, -1)]        # +delta ... 0 ... -delta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dtcsv", default="dt_trialfix.csv")
    ap.add_argument("--deltas", default="0.1,0.2,0.3")
    ap.add_argument("--outdir", default="fig2_timelags")
    ap.add_argument("--ymax", type=float, default=0.6, help="upper y limit (proportion); 0 = automatic")
    args = ap.parse_args()

    deltas = [float(x) for x in args.deltas.split(",")]
    df = pd.read_csv(args.dtcsv)
    os.makedirs(args.outdir, exist_ok=True)

    rows = []
    fig, axes = plt.subplots(1, len(deltas), figsize=(2.6 * len(deltas) + 1.2, 2.3), sharey=True)
    axes = np.atleast_1d(axes)
    letters = "abcdefg"
    for ax, delta, letter in zip(axes, deltas, letters):
        lags = lag_grid(delta)
        x = np.arange(len(lags))
        w = 0.38
        for j, cond in enumerate(("vis", "invis")):
            col = f"{delta:g}_{cond}"
            if col not in df.columns:
                raise SystemExit(f"[ERROR] column {col} not in {args.dtcsv} (columns: {list(df.columns)})")
            dt = np.round(df[col].dropna().to_numpy(float), 1)
            n = len(dt)
            counts = np.array([np.sum(np.isclose(dt, L, atol=1e-6)) for L in lags])
            prop = counts / n if n else np.zeros(len(lags))
            ax.bar(x + (j - 0.5) * w, prop, width=w, color=COL[cond], label=LABEL[cond])
            for L, c, p in zip(lags, counts, prop):
                rows.append(dict(delta=delta, condition=cond, lag_s=L, n=int(c), n_total=int(n), proportion=round(float(p), 4)))
            outside = int(np.sum(np.abs(dt) > delta + 1e-6))
            if outside:
                print(f"[WARN] delta={delta:g} {cond}: {outside} events with |dt| > delta (not plotted)")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{L:g}" for L in lags], fontsize=8)
        ax.set_xlabel(r"Time lag $\Delta t$ [s]", fontsize=9)
        ax.set_title(rf"$\delta$ = {delta:g}", fontsize=9)
        ax.text(-0.12, 1.12, f"({letter})", transform=ax.transAxes, fontsize=10, va="top", ha="left")
        ax.tick_params(axis="y", labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("Proportion", fontsize=9)
    if args.ymax > 0:
        axes[0].set_ylim(0, args.ymax)
    axes[-1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)
    fig.tight_layout()

    for fmt in ("pdf", "png"):
        fp = os.path.join(args.outdir, f"timelags.{fmt}")
        fig.savefig(fp, dpi=300 if fmt == "png" else None, bbox_inches="tight")
        print(f"Saved: {fp}")
    src = os.path.join(args.outdir, "timelags_source_data.csv")
    pd.DataFrame(rows).to_csv(src, index=False)
    print(f"Saved: {src}")
    for delta in deltas:
        sub = [r for r in rows if r["delta"] == delta]
        for cond in ("vis", "invis"):
            s = [r for r in sub if r["condition"] == cond]
            print(f"delta={delta:g} {cond:5s} n={s[0]['n_total']:5d}  " +
                  "  ".join(f"{r['lag_s']:+.1f}:{r['proportion']:.3f}" for r in s))


if __name__ == "__main__":
    main()
