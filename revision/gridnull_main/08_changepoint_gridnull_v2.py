#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_changepoint_gridnull_v2.py  (2026-08-30)

Data-driven segmentation of the session from the GRID-NULL window-level
z_S(t) curves (STAR Methods, "Period segmentation and period-resolved
analyses").

Input: one or more timelines_all_master.csv files produced by
03_Micro_Analysis_NOLPM_gridnull.py (wide format; columns z_S_t__<pair_id>
keyed by t_center), e.g. the visible and invisible main runs at delta = 0.1.

Procedure (fully deterministic, no dependencies beyond numpy/pandas):
  1. Pool all z_S_t series across the supplied files (all pairs, both
     conditions).
  2. Retain time points where at least --min_frac of the series have a valid
     (non-missing) window value.
  3. Average across series at each retained time point (group-average curve).
  4. Exact least-squares K-change-point segmentation (--n_cp, default 2):
     boundaries minimizing the within-segment sum of squared deviations from
     segment means, subject to a minimum segment length (--min_seg seconds).
     Solved by dynamic programming over segment start indices, which returns
     the same global optimum as exhaustive search (asserted for K = 2).
  5. Leave-one-series-out stability check: re-estimate the boundaries
     dropping each series in turn and report the ranges.

Usage:
  python 08_changepoint_gridnull_v2.py \
    --timelines "out_grid_visible_d0.1/timelines_all_master.csv,out_grid_invisible_d0.1/timelines_all_master.csv" \
    --outdir out_changepoint_grid
  # 3 change points:
  python 08_changepoint_gridnull_v2.py \
    --timelines "out_grid_visible_d0.1/timelines_all_master.csv,out_grid_invisible_d0.1/timelines_all_master.csv" \
    --n_cp 3 --outdir out_changepoint_grid_k3

Outputs:
  <outdir>/changepoint_group_curve.csv  t_center, n_valid, coverage, mean_zS
  <outdir>/changepoint_result.txt       boundaries + stability check

v2 (2026-09-02): collapse exact duplicate t_center blocks in
timelines_all_master.csv before segmentation (earlier runs treated the
stacked duplicates as extra time points: "retained 348 / 348", which
shifted c1 slightly, e.g. delta=0.1 K=2: 8.5 -> 9.0 s; c2=20.0 unchanged).
"""

import argparse
import os

import numpy as np
import pandas as pd


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



def best_two_changepoints(y: np.ndarray, t: np.ndarray, min_seg: float):
    """Exhaustive SSE minimization over all (c1, c2) with segment length >= min_seg.
    Boundaries are BETWEEN samples: segment1 = [0:i), segment2 = [i:j), segment3 = [j:n).
    Returns (t_c1, t_c2, sse). Boundary time = midpoint between adjacent t_centers."""
    n = len(y)
    cs = np.concatenate([[0.0], np.cumsum(y)])
    cs2 = np.concatenate([[0.0], np.cumsum(y ** 2)])

    def sse(i, j):  # SSE of y[i:j]
        m = j - i
        s = cs[j] - cs[i]
        s2 = cs2[j] - cs2[i]
        return s2 - s * s / m

    best = (None, None, np.inf)
    for i in range(1, n - 1):
        if t[i - 1] - t[0] < min_seg:
            continue
        # (unchanged exhaustive search continues below)
        if t[-1] - t[i] < 2 * min_seg:
            break
        for j in range(i + 1, n):
            if t[j - 1] - t[i] < min_seg:
                continue
            if t[-1] - t[j] < min_seg:
                break
            v = sse(0, i) + sse(i, j) + sse(j, n)
            if v < best[2]:
                best = (i, j, v)
    i, j, v = best
    b1 = 0.5 * (t[i - 1] + t[i])
    b2 = 0.5 * (t[j - 1] + t[j])
    return b1, b2, v


def best_changepoints(y: np.ndarray, t: np.ndarray, min_seg: float, n_cp: int):
    """Exact least-squares segmentation with n_cp change points (n_cp + 1
    segments), each segment [a:b) satisfying t[b-1] - t[a] >= min_seg
    (same convention as best_two_changepoints). Dynamic programming over
    segment boundaries; returns (list of boundary times, sse). The global
    optimum is identical to exhaustive search."""
    n = len(y)
    K = n_cp + 1
    cs = np.concatenate([[0.0], np.cumsum(y)])
    cs2 = np.concatenate([[0.0], np.cumsum(y ** 2)])

    def sse(i, j):  # SSE of y[i:j]
        m = j - i
        s = cs[j] - cs[i]
        s2 = cs2[j] - cs2[i]
        return s2 - s * s / m

    INF = np.inf
    dp = np.full((K + 1, n + 1), INF)
    arg = np.full((K + 1, n + 1), -1, dtype=int)
    dp[0][0] = 0.0
    for k in range(1, K + 1):
        for b in range(1, n + 1):
            for a in range(0, b):
                if dp[k - 1][a] == INF:
                    continue
                if t[b - 1] - t[a] < min_seg:
                    break  # larger a only shortens the segment
                v = dp[k - 1][a] + sse(a, b)
                if v < dp[k][b]:
                    dp[k][b] = v
                    arg[k][b] = a
    if dp[K][n] == INF:
        raise SystemExit("[ERROR] no feasible segmentation under min_seg")
    cuts = []
    b = n
    for k in range(K, 0, -1):
        a = arg[k][b]
        if k > 1:
            cuts.append(a)
        b = a
    cuts = sorted(cuts)
    bounds = [0.5 * (t[c - 1] + t[c]) for c in cuts]
    return bounds, dp[K][n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timelines", type=str, required=True,
                    help="comma-separated timelines_all_master.csv paths")
    ap.add_argument("--metric", type=str, default="z_S_t")
    ap.add_argument("--min_frac", type=float, default=0.5,
                    help="min fraction of series with valid values at a time point")
    ap.add_argument("--min_seg", type=float, default=5.0,
                    help="min segment length (s)")
    ap.add_argument("--n_cp", type=int, default=2,
                    help="number of change points (segments = n_cp + 1)")
    ap.add_argument("--outdir", type=str, default="out_changepoint_grid")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    frames = []
    for k, path in enumerate([p.strip() for p in args.timelines.split(",") if p.strip()]):
        df = pd.read_csv(path)
        cols = [c for c in df.columns if c.startswith(args.metric + "__")]
        if not cols:
            raise SystemExit(f"no {args.metric}__* columns in {path}")
        sub = df[["t_center"] + cols].copy()
        sub.columns = ["t_center"] + [f"src{k}__{c}" for c in cols]
        sub = dedupe_t_center(sub, path)
        frames.append(sub.set_index("t_center"))
        print(f"[IN] {path}: {len(cols)} series, {len(sub)} time points")
    wide = pd.concat(frames, axis=1, join="outer").sort_index()

    n_series = wide.shape[1]
    n_valid = wide.notna().sum(axis=1)
    coverage = n_valid / n_series
    keep = coverage >= args.min_frac
    mean_z = wide.mean(axis=1, skipna=True)

    curve = pd.DataFrame({
        "t_center": wide.index, "n_valid": n_valid.values,
        "coverage": coverage.values, "mean_zS": mean_z.values,
        "retained": keep.values.astype(int),
    })
    curve.to_csv(os.path.join(args.outdir, "changepoint_group_curve.csv"), index=False)

    t = wide.index.values[keep.values]
    y = mean_z.values[keep.values]
    bounds, v = best_changepoints(y, t, args.min_seg, args.n_cp)
    if args.n_cp == 2:
        # cross-check: DP must reproduce the original exhaustive search
        eb1, eb2, ev = best_two_changepoints(y, t, args.min_seg)
        assert abs(bounds[0] - eb1) < 1e-9 and abs(bounds[1] - eb2) < 1e-9 \
            and abs(v - ev) < 1e-6, "DP/exhaustive mismatch"

    btxt = ", ".join(f"c{k+1} = {b:.1f} s" for k, b in enumerate(bounds))
    lines = [
        "=" * 68,
        "Change-point segmentation of the group-average z_S(t) curve",
        f"(grid-aligned null; {n_series} series pooled across "
        f"{len(frames)} input files; min coverage = {args.min_frac:g}; "
        f"min segment = {args.min_seg:g} s; n_cp = {args.n_cp})",
        "=" * 68,
        f"retained time points: {int(keep.sum())} / {len(keep)} "
        f"(t = {t[0]:.1f} .. {t[-1]:.1f} s)",
        f"estimated boundaries: {btxt}  (SSE = {v:.4f})",
        "",
        "Leave-one-series-out stability:",
    ]
    all_b = [[] for _ in range(args.n_cp)]
    for c in wide.columns:
        sub = wide.drop(columns=[c])
        nv = sub.notna().sum(axis=1) / sub.shape[1]
        kp = (nv >= args.min_frac).values
        tt = sub.index.values[kp]
        yy = sub.mean(axis=1, skipna=True).values[kp]
        xb, _ = best_changepoints(yy, tt, args.min_seg, args.n_cp)
        for k in range(args.n_cp):
            all_b[k].append(xb[k])
    for k in range(args.n_cp):
        lines.append(f"  c{k+1} range: {min(all_b[k]):.1f} .. {max(all_b[k]):.1f} s "
                     f"(median {np.median(all_b[k]):.1f})")
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(args.outdir, "changepoint_result.txt"), "w") as f:
        f.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
