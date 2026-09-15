#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R2-4: Jitter-injection simulation — does condition-independent timing
noise attenuate zS toward zero (rather than create spurious condition
differences)?

Rationale
---------
Both reliability studies (intra-rater and inter-rater) converged on a
median absolute timing discrepancy of 0.20 s between annotation passes
(SD of signed discrepancies 0.235 s intra / 0.249 s inter).  This
simulation injects timing noise of exactly this magnitude into the
ORIGINAL coder's annotations of ALL pair videos in the repository and
recomputes the period-resolved standardized synchrony score zS
(NLOPM reimplementation, same as R2-4_zS_reproduction.py):

  Arm 1 "real"    : original annotations, unmodified (baseline).
  Arm 2 "jitter"  : every event time of BOTH participants independently
                    perturbed by N(0, jitter_sd), re-rounded to the
                    0.1-s annotation grid (--reps replicates).
  Arm 3 "null"    : participant B's series first circularly shifted by
                    a random offset (destroying any true synchrony),
                    then jittered as in Arm 2.  Demonstrates that
                    timing noise does not CREATE synchrony or condition
                    differences where none exist.

Expected pattern if the attenuation claim is correct:
  - Arm 2: mean zS shrinks toward 0 in BOTH conditions by a similar
    factor; the visible-invisible contrast shrinks but does not
    systematically reverse sign.
  - Arm 3: mean zS ~ 0 in both conditions; contrast ~ 0.

Inputs
------
--annotated_data : root of the original annotation repository
                   (contains visible_pair_Human/{A,B},
                    invisible_pair_Human/{A,B}).  ALL pair videos with
                   both A and B files are used.

Optional
--------
--jitter_sd  0.24   SD of injected Gaussian timing noise in seconds
                    (observed SD of signed coder discrepancies:
                     0.235 intra / 0.249 inter).
--miss_rate  0.0    additionally delete this fraction of events at
                    random in Arm 2/3 (e.g. 0.14 = inter-rater
                    detection-threshold difference).  Default off:
                    the primary claim concerns timing noise.
--delta      0.1    tolerance window for zS (the contested tolerance).
--reps       10     jitter replicates per video and arm.
--n_perm     1000   circular-shift permutations per zS.
--periods    P1:0-20,P2:20-60,P3:60-180
--seed       20260827
--out_prefix R2-4_jitter_sim

Outputs
-------
<out_prefix>_per_video.csv : video x condition x arm x rep x period:
                             n_a, n_b, n_match, S, zS, p
<out_prefix>_summary.txt   : noise-magnitude check + per-period
                             attenuation table + null-arm table

Run (Mac / zsh, conda env nieh_iscience_2026, from iscience_revision_analysis/R2_4_3_1):

    python R2-4_jitter_simulation.py \
        --annotated_data ../annotation_raw_data \
        --out_prefix R2-4_jitter_sim

Approximate runtime with defaults: 30-60 min for ~60 pair videos
(progress is printed per video).  For a quick smoke test add
"--reps 2 --n_perm 200 --out_prefix R2-4_jitter_smoke".
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
    """All pair videos in the repository with both A and B files."""
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
              null_shift="grid", grid=GRID):
    """Period-restricted zS.

    null_shift = "grid"       : circular-shift offsets drawn as integer
                                multiples of the annotation grid, so the
                                shifted (null) series stays grid-aligned
                                exactly like the observed series.
    null_shift = "continuous" : offsets drawn uniformly (as in
                                R2-4_zS_reproduction.py).  NOTE: for
                                grid-aligned data and delta equal to the
                                grid this null is biased LOW, which
                                inflates zS of any grid-aligned series
                                (approximately 3:2 in expected matches at
                                delta = 0.1 with a 0.1-s grid).  Kept for
                                comparison with the original pipeline."""
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


# ------------------------------------------------------- noise injection
def jitter_series(t, jitter_sd, miss_rate, rng, grid=GRID):
    """Independent Gaussian timing noise per event, re-rounded to the
    annotation grid; optional random event deletion."""
    if miss_rate > 0:
        keep = rng.uniform(size=len(t)) >= miss_rate
        t = t[keep]
    tj = t + rng.normal(0.0, jitter_sd, size=len(t))
    if grid and grid > 0:
        tj = np.round(tj / grid) * grid
    tj = tj[tj >= 0]
    return np.unique(tj)


def circular_shift_series(t, rng, grid=GRID):
    """Destroy true synchrony: circularly shift the whole series by a
    random offset within its own span (offset drawn from the central
    60% of the span to stay far from the identity).

    The offset is drawn as an integer multiple of the annotation grid
    so that the shifted series remains grid-aligned, exactly like the
    real annotations.  (A continuous offset would move events off the
    0.1-s grid and systematically alter coincidence counts at delta =
    grid, biasing zS in the null arm.)"""
    if len(t) < 2:
        return t
    t0, t1 = t.min(), t.max()
    span = t1 - t0
    if span <= 0:
        return t
    n_steps = int(round(span / grid))
    if n_steps < 2:
        return t
    k = rng.integers(int(0.2 * n_steps), max(int(0.8 * n_steps),
                                             int(0.2 * n_steps) + 1))
    u = k * grid
    ts = t0 + (t - t0 + u) % span
    ts = np.round(ts / grid) * grid          # clean float error only
    return np.unique(ts[ts >= 0])


# -------------------------------------------------- noise-magnitude check
def matched_abs_dt(orig, jit, tol=0.5):
    """Greedy one-to-one nearest matching (same spirit as the agreement
    script) to measure the realized timing discrepancy of the injected
    noise, comparable to the observed coder discrepancy."""
    dts = []
    j = 0
    used = np.zeros(len(jit), bool)
    for a in orig:
        lo = np.searchsorted(jit, a - tol)
        hi = np.searchsorted(jit, a + tol, side="right")
        best, bd = -1, None
        for k in range(lo, hi):
            if used[k]:
                continue
            d = abs(jit[k] - a)
            if bd is None or d < bd:
                best, bd = k, d
        if best >= 0:
            used[best] = True
            dts.append(bd)
    return np.array(dts)


def f1_at(orig, jit, tol):
    if len(orig) == 0 or len(jit) == 0:
        return np.nan
    m = len(matched_abs_dt(orig, jit, tol))
    prec = m / len(jit)
    rec = m / len(orig)
    return 2 * prec * rec / (prec + rec) if prec + rec > 0 else 0.0


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description="Jitter-injection attenuation simulation (R2-4)")
    ap.add_argument("--annotated_data", required=True)
    ap.add_argument("--jitter_sd", type=float, default=0.24)
    ap.add_argument("--miss_rate", type=float, default=0.0)
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--n_perm", type=int, default=1000)
    ap.add_argument("--periods", default="P1:0-20,P2:20-60,P3:60-180")
    ap.add_argument("--null_shift", choices=["grid", "continuous"],
                    default="grid",
                    help="circular-shift null: grid-aligned offsets "
                         "(unbiased for grid data; default) or "
                         "continuous offsets (as in the original "
                         "reimplementation; biased at delta = grid)")
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--out_prefix", default="R2-4_jitter_sim")
    args = ap.parse_args()

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

    master = np.random.SeedSequence(args.seed)
    rows = []
    noise_dts, noise_f1_01, noise_f1_02 = [], [], []

    for vi, v in enumerate(videos):
        tl = read_event_series(v["left"])
        tr = read_event_series(v["right"])
        child = np.random.SeedSequence([args.seed, vi])
        rng = np.random.default_rng(child)

        # Arm 1: real
        for name, p0, p1 in periods:
            res = zs_period(tl, tr, p0, p1, args.delta, args.n_perm, rng,
                            null_shift=args.null_shift)
            rows.append(dict(video=v["video"], condition=v["condition"],
                             arm="real", rep=0, period=name, **res))

        for rep in range(1, args.reps + 1):
            # Arm 2: jitter
            tlj = jitter_series(tl, args.jitter_sd, args.miss_rate, rng)
            trj = jitter_series(tr, args.jitter_sd, args.miss_rate, rng)
            if rep == 1:                      # noise-magnitude check
                for o, jt in ((tl, tlj), (tr, trj)):
                    d = matched_abs_dt(o, jt, tol=0.5)
                    noise_dts.append(d)
                    noise_f1_01.append(f1_at(o, jt, 0.1))
                    noise_f1_02.append(f1_at(o, jt, 0.2))
            for name, p0, p1 in periods:
                res = zs_period(tlj, trj, p0, p1, args.delta,
                                args.n_perm, rng,
                                null_shift=args.null_shift)
                rows.append(dict(video=v["video"],
                                 condition=v["condition"],
                                 arm="jitter", rep=rep, period=name, **res))
            # Arm 3: null (shift B, then jitter both)
            trn = circular_shift_series(tr, rng)
            tln_j = jitter_series(tl, args.jitter_sd, args.miss_rate, rng)
            trn_j = jitter_series(trn, args.jitter_sd, args.miss_rate, rng)
            for name, p0, p1 in periods:
                res = zs_period(tln_j, trn_j, p0, p1, args.delta,
                                args.n_perm, rng,
                                null_shift=args.null_shift)
                rows.append(dict(video=v["video"],
                                 condition=v["condition"],
                                 arm="null", rep=rep, period=name, **res))
        print(f"  [{vi + 1}/{len(videos)}] {v['video']} "
              f"({v['condition']}) done")

    df = pd.DataFrame(rows)
    per_video_csv = f"{args.out_prefix}_per_video.csv"
    df.to_csv(per_video_csv, index=False)

    # ------------------------------------------------------------ summary
    dt_all = np.concatenate(noise_dts) if noise_dts else np.array([])
    lines = []
    w = lines.append
    w("R2-4 jitter-injection attenuation simulation")
    w("=" * 60)
    w(f"Pair videos            : {len(videos)} "
      f"({n_v} visible, {n_i} invisible)")
    w(f"delta (zS tolerance)   : {args.delta}")
    w(f"jitter SD              : {args.jitter_sd} s "
      f"(observed coder SD: 0.235 intra / 0.249 inter)")
    w(f"miss rate              : {args.miss_rate}")
    w(f"replicates x n_perm    : {args.reps} x {args.n_perm}")
    w(f"null shift             : {args.null_shift}")
    w(f"seed                   : {args.seed}")
    w("")
    w("Noise-magnitude check (rep 1, injected vs original):")
    if len(dt_all):
        w(f"  median |dt| = {np.median(dt_all):.3f} s "
          f"(observed coder value: 0.200 s)")
        w(f"  mean F1 @0.1 = {np.nanmean(noise_f1_01):.3f} "
          f"(observed: 0.536 intra / 0.450 inter)")
        w(f"  mean F1 @0.2 = {np.nanmean(noise_f1_02):.3f} "
          f"(observed: 0.784 intra / 0.694 inter)")
    w("")

    for name, _, _ in periods:
        d = df[df["period"] == name]
        real = d[d["arm"] == "real"]
        rv = real[real["condition"] == "visible"]["zS"].mean()
        ri = real[real["condition"] == "invisible"]["zS"].mean()
        w(f"--- {name} ---")
        w(f"  real    : mean zS  visible {rv:+.3f}   invisible {ri:+.3f}"
          f"   contrast (V-I) {rv - ri:+.3f}")
        for arm in ("jitter", "null"):
            a = d[d["arm"] == arm]
            if a.empty:
                continue
            cv, ci, cc = [], [], []
            for rep, g in a.groupby("rep"):
                mv = g[g["condition"] == "visible"]["zS"].mean()
                mi = g[g["condition"] == "invisible"]["zS"].mean()
                cv.append(mv); ci.append(mi); cc.append(mv - mi)
            cv, ci, cc = map(np.asarray, (cv, ci, cc))
            w(f"  {arm:<7} : mean zS  visible {cv.mean():+.3f}   "
              f"invisible {ci.mean():+.3f}   contrast {cc.mean():+.3f} "
              f"(SD over reps {cc.std(ddof=1):.3f})")
            if arm == "jitter":
                att_v = cv.mean() / rv if rv else np.nan
                att_i = ci.mean() / ri if ri else np.nan
                same_sign = np.mean(np.sign(cc) == np.sign(rv - ri))
                w(f"            attenuation factor  visible "
                  f"{att_v:.2f}   invisible {att_i:.2f}")
                w(f"            replicates preserving contrast sign: "
                  f"{same_sign:.0%} ({int(same_sign * len(cc))}/{len(cc)})")
            if arm == "null":
                exceed = np.mean(np.abs(cc) >= abs(rv - ri))
                w(f"            replicates with |null contrast| >= "
                  f"|real contrast|: {exceed:.0%} "
                  f"({int(round(exceed * len(cc)))}/{len(cc)})")
        w("")

    summary_txt = f"{args.out_prefix}_summary.txt"
    with open(summary_txt, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote: {per_video_csv}\nWrote: {summary_txt}")


if __name__ == "__main__":
    main()
