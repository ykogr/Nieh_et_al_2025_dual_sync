#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coverage_table.py -- count how many pairs cover each time bin (t_center).

When the time axis is video time, each pair's series starts at the whole second after its
actual onset (3-20 s), so the leading bins have group means determined by only a few pairs.
Produces a table for deciding where to start drawing Fig 4 / Fig S2, and a coverage plot.

Usage:
  python coverage_table.py --timelines out_grid_visible_d0.1/timelines_all_master.csv,out_grid_invisible_d0.1/timelines_all_master.csv \
      --labels visible,invisible --min_pairs 15 --outdir coverage_grid
Outputs:
  <outdir>/coverage_by_bin.csv      t_center, n_<label>..., n_total
  <outdir>/coverage_summary.txt     first/last t_center with at least min_pairs, etc.
  <outdir>/coverage.png
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timelines", required=True, help="comma-separated wide/master csv files")
    ap.add_argument("--labels", default="visible,invisible")
    ap.add_argument("--min_pairs", type=int, default=15,
                    help="recommended lower bound (bins below this count are dropped from the figure or annotated)")
    ap.add_argument("--outdir", default="coverage_grid")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    files = args.timelines.split(",")
    labels = args.labels.split(",")
    assert len(files) == len(labels), "--timelines and --labels have different counts"

    out = None
    lines = []
    for f, lab in zip(files, labels):
        w = pd.read_csv(f)
        cols = [c for c in w.columns if c.startswith("z_S_t__")]
        cov = pd.DataFrame({"t_center": w["t_center"].round(1), f"n_{lab}": w[cols].notna().sum(axis=1).values})
        out = cov if out is None else out.merge(cov, on="t_center", how="outer")
        ok = cov[cov[f"n_{lab}"] >= args.min_pairs]
        lines.append(f"{lab}: {len(cols)} pairs, bins {cov.t_center.min():.1f}-{cov.t_center.max():.1f} s; "
                     f">= {args.min_pairs} pairs from {ok.t_center.min():.1f} to {ok.t_center.max():.1f} s; "
                     f"full coverage ({len(cols)} pairs) from "
                     f"{cov[cov[f'n_{lab}'] == len(cols)].t_center.min():.1f} s")
    out = out.sort_values("t_center").fillna(0)
    ncols = [c for c in out.columns if c.startswith("n_")]
    out["n_total"] = out[ncols].sum(axis=1).astype(int)
    for c in ncols:
        out[c] = out[c].astype(int)
    out.to_csv(os.path.join(args.outdir, "coverage_by_bin.csv"), index=False)

    fig, ax = plt.subplots(figsize=(8, 3))
    for c in ncols:
        ax.plot(out["t_center"], out[c], label=c[2:])
    ax.axhline(args.min_pairs, color="grey", ls="--", lw=0.8, label=f"min_pairs = {args.min_pairs}")
    ax.set_xlabel("Time in session (s)")
    ax.set_ylabel("pairs covering the bin")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "coverage.png"), dpi=150)

    txt = "\n".join(lines) + "\n"
    with open(os.path.join(args.outdir, "coverage_summary.txt"), "w") as fh:
        fh.write(txt)
    print(txt, end="")
    print(f"[out] {args.outdir}/coverage_by_bin.csv, coverage_summary.txt, coverage.png")


if __name__ == "__main__":
    main()
