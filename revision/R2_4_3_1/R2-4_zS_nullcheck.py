#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R2-4 / discretization check: quantify how the circular-shift null
convention affects the period-resolved zS on the ORIGINAL annotations.

Background
----------
The published pipeline (03_Micro_Analysis_NOLPM.py) draws circular-shift
offsets as continuous uniform random numbers (u = rng.uniform(0, T)) and
uses the annotation times as-is.  The annotations lie on a 0.1-s grid.
For grid-aligned data, distances between the two participants' events
are multiples of 0.1 s, so at delta = 0.1 s a pair of grid-aligned
series can match at inter-event distances {-0.1, 0, +0.1} (three grid
cells), whereas a continuously shifted (off-grid) null series can only
match within a 0.2-s-wide continuous window (two grid cells' worth of
probability mass).  The expected match count of the observed data is
therefore ~3:2 higher than the null's even in the complete absence of
synchrony, which biases zS upward (the bias shrinks as delta grows:
5:4 at delta = 0.2, 7:6 at 0.3, 11:10 at 0.5, 21:20 at 1.0).

This script recomputes period-resolved zS for ALL pair videos in the
repository under BOTH null conventions:

  null = "continuous" : offsets uniform in (0, T)   [published pipeline]
  null = "grid"       : offsets = integer multiples of the 0.1-s grid
                        [unbiased for grid-aligned data]

and reports, per delta x period x condition:
  - mean zS under each null and their difference (the artifact size),
  - per-video-period "prevalence" counts (one-sided p < 0.05),
  - the visible-invisible contrast in mean zS under each null.

Everything is computed from the original annotations only; no
re-annotation data is involved.

Inputs
------
--annotated_data : root of the original annotation repository
                   (visible_pair_Human/{A,B}, invisible_pair_Human/{A,B});
                   ALL pair videos with both A and B files are used.

Optional
--------
--deltas     0.1,0.2,0.3,0.5,1.0   (full manuscript grid; default)
--n_perm     2000    permutations per video x period x delta x null.
                     (2000 keeps the total runtime ~1 h; increase for
                     final citable numbers if desired, e.g. 5000.)
--periods    P1:0-20,P2:20-60,P3:60-180
--seed       20260827
--out_prefix R2-4_nullcheck

Outputs
-------
<out_prefix>_per_video.csv : video x condition x period x delta x null:
                             n_a, n_b, n_match, S, zS, p
<out_prefix>_summary.txt   : artifact-size table + prevalence table +
                             condition-contrast table

Run (Mac / zsh, conda env nieh_iscience_2026, from
iscience_revision_analysis/R2_4_3_1):

    python R2-4_zS_nullcheck.py \
        --annotated_data ../annotation_raw_data \
        --out_prefix R2-4_nullcheck

Approximate runtime with defaults: ~1 hour for ~60 pair videos
(progress printed per video).  Quick smoke test:
"--deltas 0.1 --n_perm 200 --out_prefix R2-4_nullcheck_smoke" (minutes).
"""

import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd


GRID = 0.1


# ---------------------------------------------------------------- I/O
def read_event_series(path, grid=GRID):
    df = pd.read_csv(path, dtype=str)
    chew_cols = [c for c in df.columns if c.lower().startswith("chew")]
    if not chew_cols:
        raise ValueError(f"No Chew columns found in {path}")
    vals = [pd.to_numeric(df[c], errors="coerce") for c in chew_cols]
    t = pd.concat(vals).dropna().to_numpy(dtype=float)
    if grid and grid > 0:
        t = np.round(t / grid) * grid
    t = t[t >= 0]
    t = np.unique(t)
    return t


def strip_participant_suffix(fname):
    base = fname[:-4] if fname.lower().endswith(".csv") else fname
    if base.endswith(("_A", "_B")):
        base = base[:-2]
    return base


def collect_all_pair_videos(annotated_data):
    layout = [("visible", "visible_pair_Human"),
              ("invisible", "invisible_pair_Human")]
    out = []
    for cond, folder in layout:
        a_dir = os.path.join(annotated_data, folder, "A")
        b_dir = os.path.join(annotated_data, folder, "B")
        if not (os.path.isdir(a_dir) and os.path.isdir(b_dir)):
            sys.exit(f"ERROR: folder not found: {a_dir} or {b_dir}")
        b_map = {strip_participant_suffix(os.path.basename(p)): p
                 for p in glob.glob(os.path.join(b_dir, "*.csv"))}
        for p in sorted(glob.glob(os.path.join(a_dir, "*.csv"))):
            v = strip_participant_suffix(os.path.basename(p))
            if v in b_map:
                out.append(dict(video=v, condition=cond,
                                left=p, right=b_map[v]))
    if not out:
        sys.exit("ERROR: no pair videos found.")
    return out


# ----------------------------------------------------------------- NLOPM
# [FLOAT-TOL PATCH 2026-09-07] tolerance for comparisons against delta.
# Event times are on a 0.1-s grid; in binary floating point the difference of
# two grid values is not exactly a multiple of 0.1 (e.g. 0.8 - 0.7 =
# 0.10000000000000009 > 0.1), so events exactly delta apart were matched or
# not depending on their float representation. 1e-6 s is far below the grid.
FLOAT_TOL = float(os.environ.get("NLOPM_FLOAT_TOL", "1e-6"))



def nlopm_matches(A, B, delta):
    n, m = len(A), len(B)
    i, j = 0, 0
    n_match = 0
    while i < n and j < m:
        a = A[i]
        while j < m and B[j] < a - delta - FLOAT_TOL:  # [FLOAT-TOL PATCH]
            j += 1
        if j >= m:
            break
        k = j
        cands = []
        while k < m and B[k] <= a + delta + FLOAT_TOL:  # [FLOAT-TOL PATCH]
            cands.append(k)
            k += 1
        if not cands:
            i += 1
            continue
        best = min(cands, key=lambda kk: (round(abs(B[kk] - a), 6), B[kk]))  # [FLOAT-TOL PATCH]
        if i + 1 < n and abs(A[i + 1] - B[best]) < abs(a - B[best]) - FLOAT_TOL:  # [FLOAT-TOL PATCH]
            earlier = [kk for kk in cands if kk < best]
            if earlier:
                best2 = min(earlier, key=lambda kk: (round(abs(B[kk] - a), 6), B[kk]))  # [FLOAT-TOL PATCH]
                n_match += 1
                j = best2 + 1
            i += 1
            continue
        n_match += 1
        j = best + 1
        i += 1
    return n_match


def zs_period(t_left, t_right, p0, p1, delta, n_perm, rng,
              null_shift, grid=GRID):
    if len(t_left) == 0 or len(t_right) == 0:
        return dict(n_a=0, n_b=0, n_match=0, S=np.nan, zS=np.nan, p=np.nan)
    tmin = max(t_left.min(), t_right.min(), p0)
    tmax = min(t_left.max(), t_right.max(), p1)
    if tmax <= tmin:
        return dict(n_a=0, n_b=0, n_match=0, S=np.nan, zS=np.nan, p=np.nan)
    A = t_left[(t_left >= tmin) & (t_left <= tmax)]
    B = t_right[(t_right >= tmin) & (t_right <= tmax)]
    if len(A) == 0 or len(B) == 0:
        return dict(n_a=len(A), n_b=len(B), n_match=0,
                    S=np.nan, zS=np.nan, p=np.nan)
    T = tmax - tmin
    n_match = nlopm_matches(A, B, delta)
    s_obs = 2.0 * n_match / (len(A) + len(B))
    n_steps = max(int(round(T / grid)), 2)
    s_null = np.empty(n_perm)
    for r in range(n_perm):
        if null_shift == "grid":
            u = rng.integers(1, n_steps) * grid
            Bs = tmin + (B - tmin + u) % T
            Bs = np.round(Bs / grid) * grid   # clean float error only
            Bs = np.unique(Bs)
        else:
            u = rng.uniform(0.0, T)
            Bs = np.sort(tmin + (B - tmin + u) % T)
        s_null[r] = 2.0 * nlopm_matches(A, Bs, delta) \
            / (len(A) + len(B))
    sd = s_null.std(ddof=1)
    zs = (s_obs - s_null.mean()) / sd if sd > 0 else np.nan
    p = (np.sum(s_null >= s_obs) + 1) / (n_perm + 1)
    return dict(n_a=len(A), n_b=len(B), n_match=n_match,
                S=s_obs, zS=zs, p=p)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description="zS null-convention comparison on original "
                    "annotations (R2-4 discretization check)")
    ap.add_argument("--annotated_data", required=True)
    ap.add_argument("--deltas", default="0.1,0.2,0.3,0.5,1.0")
    ap.add_argument("--n_perm", type=int, default=2000)
    ap.add_argument("--periods", default="P1:0-20,P2:20-60,P3:60-180")
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--out_prefix", default="R2-4_nullcheck")
    args = ap.parse_args()

    deltas = [float(x) for x in args.deltas.split(",")]
    periods = []
    for tok in args.periods.split(","):
        name, rng_ = tok.split(":")
        a, b = rng_.split("-")
        periods.append((name, float(a), float(b)))

    videos = collect_all_pair_videos(args.annotated_data)
    n_v = sum(1 for v in videos if v["condition"] == "visible")
    n_i = sum(1 for v in videos if v["condition"] == "invisible")
    print(f"Found {len(videos)} pair videos "
          f"({n_v} visible, {n_i} invisible)")

    rows = []
    for vi, v in enumerate(videos):
        tl = read_event_series(v["left"])
        tr = read_event_series(v["right"])
        for di, delta in enumerate(deltas):
            for ni, null_shift in enumerate(("continuous", "grid")):
                child = np.random.SeedSequence([args.seed, vi, di, ni])
                rng = np.random.default_rng(child)
                for name, p0, p1 in periods:
                    res = zs_period(tl, tr, p0, p1, delta,
                                    args.n_perm, rng, null_shift)
                    rows.append(dict(video=v["video"],
                                     condition=v["condition"],
                                     period=name, delta=delta,
                                     null=null_shift, **res))
        print(f"  [{vi + 1}/{len(videos)}] {v['video']} "
              f"({v['condition']}) done")

    df = pd.DataFrame(rows)
    per_video_csv = f"{args.out_prefix}_per_video.csv"
    df.to_csv(per_video_csv, index=False)

    lines = []
    w = lines.append
    w("R2-4 zS null-convention comparison (original annotations only)")
    w("=" * 66)
    w(f"Pair videos    : {len(videos)} ({n_v} visible, {n_i} invisible)")
    w(f"deltas         : {deltas}")
    w(f"n_perm         : {args.n_perm}   seed: {args.seed}")
    w("null=continuous: offsets ~ U(0,T)  [published pipeline]")
    w("null=grid      : offsets = k * 0.1 s  [grid-aligned, unbiased]")
    w("")
    w("Theoretical inflation of expected matches (continuous null),")
    w("grid-aligned data: 3:2 @0.1, 5:4 @0.2, 7:6 @0.3, 11:10 @0.5, "
      "21:20 @1.0")
    w("")

    for name, _, _ in periods:
        w(f"===== {name} =====")
        w(f"{'delta':>6} {'null':>11} | {'mean zS V':>10} {'mean zS I':>10}"
          f" {'V-I':>7} | {'prev V':>7} {'prev I':>7}")
        for delta in deltas:
            for null_shift in ("continuous", "grid"):
                d = df[(df["period"] == name) & (df["delta"] == delta)
                       & (df["null"] == null_shift)]
                dv = d[d["condition"] == "visible"]
                di_ = d[d["condition"] == "invisible"]
                mv, mi = dv["zS"].mean(), di_["zS"].mean()
                pv = int((dv["p"] < 0.05).sum())
                pi_ = int((di_["p"] < 0.05).sum())
                w(f"{delta:>6} {null_shift:>11} | {mv:>+10.3f} {mi:>+10.3f}"
                  f" {mv - mi:>+7.3f} | {pv:>4}/{len(dv):<2} "
                  f"{pi_:>4}/{len(di_):<2}")
            # artifact size for this delta
            dc = df[(df["period"] == name) & (df["delta"] == delta)]
            piv = dc.pivot_table(index=["video", "condition"],
                                 columns="null", values="zS")
            if {"continuous", "grid"} <= set(piv.columns):
                diff = (piv["continuous"] - piv["grid"]).dropna()
                w(f"{'':>6} {'artifact':>11} | mean zS(cont - grid) = "
                  f"{diff.mean():+.3f} (SD {diff.std(ddof=1):.3f}, "
                  f"n = {len(diff)})")
        w("")

    summary_txt = f"{args.out_prefix}_summary.txt"
    with open(summary_txt, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote: {per_video_csv}\nWrote: {summary_txt}")


if __name__ == "__main__":
    main()
