#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
09_plot_figure4ad_gridnull_v2.py

Regenerate Figure 4 panels (a)-(d) from the grid-null timelines produced by
03_Micro_Analysis_NOLPM_gridnull.py (timelines_all_master.csv).

Panel layout matches the current manuscript figure (p. 7 of the 2026-08-23
revision PDF):
  (a) visible   z_S(t)        (orange)
  (b) invisible z_S(t)        (blue)
  (c) visible   S_excess(t) %  (orange)
  (d) invisible S_excess(t) %  (blue)
Lines: group mean across pairs; shaded bands: 95% CI of the group mean
(t-distribution, mean +/- t_{0.975, n-1} * SEM, n = pairs with a valid value
at that time point).

The script also writes peaks.txt reporting, for each panel, the peak of the
group-mean curve and the 95% CI band at the peak time point (the quantities
quoted in the Figure 4 caption and Results text, e.g. "3.80 [3.55-4.08]").
Panel (e) (Stan beta(t)) is NOT drawn here; it depends on the delta = 0.1
re-fit and will be added by a separate step after that fit is accepted.

Usage (local, from ROOT/gridnull_main, conda env nieh_iscience_2026):
  python 09_plot_figure4ad_gridnull_v2.py \
      --timelines "out_grid_visible_d0.1/timelines_all_master.csv,out_grid_invisible_d0.1/timelines_all_master.csv" \
      --labels "visible,invisible" --outdir fig4ad_gridnull

Options:
  --timelines CSV1,CSV2  one timelines_all_master.csv per condition
  --labels A,B           condition labels (default: visible,invisible)
  --colors C1,C2         line colors per condition
                         (default: tab orange / tab blue, as in the current
                         manuscript Figure 4; note Figure S1/S2 use the
                         Okabe-Ito palette with the opposite assignment)
  --tmax T               x-axis right edge in seconds (default: 180)
  --sharey row|none      row: panels (a,b) and (c,d) share one y-range each
                         (default, v2.1); none: independent y-axes (v2)
  --ylim_min_pairs N     when computing the shared y-range, use only time
                         points with at least N pairs in BOTH conditions
                         (default 0 = all points). The curves themselves are
                         always drawn in full; this only decides how much of
                         the wide early CI band (few pairs) sets the scale.
                         N = 15 matches coverage_table.py --min_pairs 15.
  --formats ...          output formats (default: pdf png)
  --outdir DIR           output directory (default: fig4ad_gridnull)

Outputs:
  <outdir>/figure4ad_gridnull.<fmt>      the 2x2 figure
  <outdir>/figure4ad_source_data.csv     per-condition mean / CI curves
  <outdir>/figure4ad_peaks.txt           peak values with 95% CI
  <outdir>/figure4ad_ylim.txt            y-range used for each row (v2.1)

v2 (2026-09-02): timelines_all_master.csv can contain the same table
stacked more than once (it concatenates all *__wide.csv in the folder).
v1 assumed unique t_center rows, so each curve was drawn repeatedly with
return strokes and the source CSV was inflated. v2 verifies duplicates are
identical and collapses them. Mean/CI/peak values were not affected.

v2.1 (2026-09-10, trialfix_kit v4.3.3): the two panels of each row now share
the y-axis (a/b: z_S(t); c/d: S_excess(t)), so that visible and invisible
curves can be compared at a glance. The common range is the union of the two
conditions' CI bands with a 5% margin (--sharey row, default); --sharey none
restores the independent axes of v2. Mean/CI/peak values and the source-data
CSV are unchanged; the chosen y-ranges are written to figure4ad_ylim.txt.
Optionally --ylim_min_pairs N restricts the range computation to time points
covered by >= N pairs (the first seconds have very few pairs and a very wide
CI band that would otherwise dominate the scale).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def dedupe_t_center(df, path):
    """Collapse exact duplicate t_center blocks in timelines_all_master.csv.

    The master file is built by vertically stacking every *__wide.csv in the
    output folder, so it can contain the same 174-row table more than once.
    Duplicates are only removed after verifying that all rows sharing a
    t_center are identical in every column (NaN == NaN); any conflict aborts.
    """
    if not df["t_center"].duplicated().any():
        print(f"[check] {path}: no duplicate t_center rows")
        return df
    nun = df.groupby("t_center").nunique(dropna=False)
    bad = [(t, c) for c in nun.columns for t in nun.index[nun[c] > 1]]
    if bad:
        raise SystemExit(
            f"[ERROR] {path}: rows sharing a t_center DIFFER "
            f"(e.g. t={bad[0][0]}, column {bad[0][1]}); refusing to guess. "
            f"Inspect the *__wide.csv files in that folder.")
    n0 = len(df)
    df = df.drop_duplicates(subset=["t_center"]).sort_values("t_center")
    print(f"[check] {path}: {n0} rows contain exact duplicate blocks; "
          f"collapsed to {len(df)} unique t_center rows "
          f"(verified identical in all columns before collapsing)")
    return df



def load_metric(path, metric):
    df = pd.read_csv(path)
    df = dedupe_t_center(df, path)
    cols = [c for c in df.columns if c.startswith(metric + "__")]
    if not cols:
        raise SystemExit(f"no {metric}__* columns in {path}")
    return df[["t_center"] + cols].set_index("t_center")


def mean_ci(sub):
    """Group mean and 95% t-CI of the mean across pair columns."""
    m = sub.mean(axis=1, skipna=True)
    n = sub.notna().sum(axis=1)
    sem = sub.std(axis=1, skipna=True) / np.sqrt(n)
    tcrit = pd.Series(np.where(n > 1, stats.t.ppf(0.975, np.maximum(n - 1, 1)),
                               np.nan), index=sub.index)
    return m, m - tcrit * sem, m + tcrit * sem, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timelines", type=str, required=True)
    ap.add_argument("--labels", type=str, default="visible,invisible")
    ap.add_argument("--colors", type=str, default="#ff7f0e,#1f77b4")
    ap.add_argument("--tmax", type=float, default=180.0)
    ap.add_argument("--sharey", type=str, default="row", choices=["row", "none"])
    ap.add_argument("--ylim_min_pairs", type=int, default=0)
    ap.add_argument("--formats", nargs="+", default=["pdf", "png"])
    ap.add_argument("--outdir", type=str, default="fig4ad_gridnull")
    args = ap.parse_args()

    paths = [p.strip() for p in args.timelines.split(",") if p.strip()]
    labels = [s.strip() for s in args.labels.split(",") if s.strip()]
    colors = [c.strip() for c in args.colors.split(",") if c.strip()]
    if not (len(paths) == len(labels) == len(colors) == 2):
        raise SystemExit("--timelines/--labels/--colors must each list 2 items")

    os.makedirs(args.outdir, exist_ok=True)
    metrics = [("z_S_t", "$z$-value", 1.0), ("S_excess", "$S_{excess}$ (%)", 100.0)]
    short = {"visible": "vis", "invisible": "invis"}

    fig, axes = plt.subplots(2, 2, figsize=(8.4, 5.6), sharex=True,
                             sharey="row" if args.sharey == "row" else False)
    panel = iter("abcd")
    peak_lines, src = [], None
    ylim_lines = []
    for r, (metric, ylab, scale) in enumerate(metrics):
        for k, (lab, path, col) in enumerate(zip(labels, paths, colors)):
            sub = load_metric(path, metric)
            m, lo, hi, n = mean_ci(sub)
            m, lo, hi = m * scale, lo * scale, hi * scale
            ax = axes[r, k]
            ax.fill_between(m.index.values, lo.values, hi.values, color=col,
                            alpha=0.25, lw=0)
            ax.plot(m.index.values, m.values, color=col, lw=1.0)
            ax.axhline(0.0, color="#999999", lw=0.5)
            tag = short.get(lab, lab)
            name = "$z_S(t)$" if metric == "z_S_t" else "$S_{excess}(t)$"
            ax.set_title(f"({next(panel)}) {tag} {name}", fontsize=10, loc="left")
            ax.set_xlim(0, args.tmax)
            if k == 0:
                ax.set_ylabel(ylab)
            if r == 1:
                ax.set_xlabel("Time [s]")
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            # peak of the group-mean curve and CI band at that time point
            i = int(np.nanargmax(m.values))
            peak_lines.append(
                f"{metric} {lab}: peak mean = {m.values[i]:.2f} "
                f"[{lo.values[i]:.2f}, {hi.values[i]:.2f}] "
                f"at t_center = {m.index.values[i]:g} s (n = {int(n.values[i])} pairs)")
            # source data
            d = pd.DataFrame({"t_center": m.index.values,
                              f"{metric}_{lab}_mean": m.values,
                              f"{metric}_{lab}_ci_lo": lo.values,
                              f"{metric}_{lab}_ci_hi": hi.values,
                              f"{metric}_{lab}_n": n.values}).set_index("t_center")
            src = d if src is None else src.join(d, how="outer")
        # v2.1: common y-range per row = union of both conditions' CI bands (+5% margin)
        if args.sharey == "row":
            ncols = [c for c in src.columns if c.startswith(metric) and c.endswith("_n")]
            ok = (src[ncols].fillna(0) >= args.ylim_min_pairs).all(axis=1)
            lo_all = np.nanmin(np.concatenate([src.loc[ok, c].values for c in src.columns
                                               if c.startswith(metric) and c.endswith("_ci_lo")]))
            hi_all = np.nanmax(np.concatenate([src.loc[ok, c].values for c in src.columns
                                               if c.startswith(metric) and c.endswith("_ci_hi")]))
            pad = 0.05 * (hi_all - lo_all)
            axes[r, 0].set_ylim(lo_all - pad, hi_all + pad)   # propagates to axes[r, 1]
            ylim_lines.append(
                f"{metric}: shared y-range = [{lo_all - pad:.3f}, {hi_all + pad:.3f}] "
                f"(union of 95% CI bands, {labels[0]} and {labels[1]}, +5% margin; "
                f"computed from {int(ok.sum())} of {len(ok)} time points with >= "
                f"{args.ylim_min_pairs} pairs in both conditions)")
        else:
            for k in range(2):
                ylim_lines.append(f"{metric} {labels[k]}: independent y-range = "
                                  f"[{axes[r, k].get_ylim()[0]:.3f}, {axes[r, k].get_ylim()[1]:.3f}]")

    fig.tight_layout()
    for fmt in args.formats:
        fp = os.path.join(args.outdir, f"figure4ad_gridnull.{fmt}")
        fig.savefig(fp, dpi=300 if fmt == "png" else None)
        print(f"[out] {fp}")
    src_csv = os.path.join(args.outdir, "figure4ad_source_data.csv")
    src.reset_index().to_csv(src_csv, index=False)
    print(f"[out] {src_csv}")
    pk = os.path.join(args.outdir, "figure4ad_peaks.txt")
    with open(pk, "w") as f:
        f.write("Peak of group-mean curve with 95% CI band at the peak time point\n")
        f.write("(S_excess values are proportions x100, i.e. percent)\n\n")
        f.write("\n".join(peak_lines) + "\n")
    print(f"[out] {pk}")
    print("\n".join(peak_lines))
    yl = os.path.join(args.outdir, "figure4ad_ylim.txt")
    with open(yl, "w") as f:
        f.write(f"--sharey {args.sharey}  --ylim_min_pairs {args.ylim_min_pairs}\n"
                + "\n".join(ylim_lines) + "\n")
    print(f"[out] {yl}")
    print("\n".join(ylim_lines))


if __name__ == "__main__":
    sys.exit(main())
