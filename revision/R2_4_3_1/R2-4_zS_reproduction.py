#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R2-4_zS_reproduction.py
Coder-robustness check of the paper's key micro-scale statistic
(Nieh et al. 2026, iScience revision; reviewer comment R2-4).

For every PAIR video annotated by the second coder, this script computes
the period-resolved standardized synchrony score zS twice:
  (1) from the ORIGINAL coder's annotations of the two participants, and
  (2) from the SECOND coder's annotations of the same two participants,
following the STAR Methods pipeline:
  - NLOPM (nearest-with-look-ahead order-preserving matching) at
    tolerance delta (default 0.1 s),
  - synchronization rate S = 2 * n_match / (n_a + n_b),
  - null distribution from circular shifts of participant B's events
    within the period-restricted overlapping interval
    (default 10,000 permutations),
  - zS = (S_obs - mean(S_null)) / SD(S_null),
  - periods P1 (0-20 s), P2 (20-60 s), P3 (60-180 s).

NOTE: this is an independent reimplementation of the NLOPM matcher from
the STAR Methods description for the reliability check; at delta = 0.1 s
(mean inter-chew interval 0.70 s) tolerance windows almost always
contain at most one candidate, where all order-preserving nearest
matchers coincide.

Usage (from the folder that contains IRR_final):
  python R2-4_zS_reproduction.py --annotated_data ../annotation_raw_data \
      --irr_root IRR_final --out_prefix R2-4_zS

Outputs:
  <out_prefix>_per_video.csv : video x period x coder: n_a, n_b,
                               n_match, S, zS, Monte Carlo p
  <out_prefix>_summary.txt   : cross-coder agreement of zS per period
                               (ICC(A,1), Pearson r) + condition means
"""

import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd


# ---------------------------------------------------------------- I/O
def read_event_series(path, grid=0.1):
    df = pd.read_csv(path, dtype=str)
    chew_cols = [c for c in df.columns if c.lower().startswith("chew")]
    if not chew_cols:
        raise ValueError(f"No Chew columns found in {path}")
    vals = [pd.to_numeric(df[c], errors="coerce") for c in chew_cols]
    t = pd.concat(vals).dropna().to_numpy(dtype=float)
    if grid and grid > 0:
        t = np.round(t / grid) * grid
    t = t[t >= 0]
    t = np.unique(t)          # remove duplicates, sort (STAR Methods)
    return t


def strip_participant_suffix(fname):
    base = fname[:-4] if fname.lower().endswith(".csv") else fname
    if base.endswith(("_A", "_B")):
        base = base[:-2]
    return base


def find_original(dirpath, video):
    hits = [p for p in glob.glob(os.path.join(dirpath, "*.csv"))
            if strip_participant_suffix(os.path.basename(p)) == video]
    if len(hits) != 1:
        sys.exit(f"ERROR: expected exactly one original CSV for '{video}' "
                 f"in {dirpath}, found {len(hits)}")
    return hits[0]


def find_second_coder(dirpath):
    hits = sorted(glob.glob(os.path.join(dirpath, "*.csv")))
    if len(hits) != 1:
        sys.exit(f"ERROR: expected exactly one CSV in {dirpath}, "
                 f"found {len(hits)}")
    return hits[0]


def collect_pair_videos(annotated_data, irr_root):
    """Return list of dicts: video, condition, per-coder file paths for
    both participants (left = A, right = B)."""
    layout = [("visible", "visible_pair_Human"),
              ("invisible", "invisible_pair_Human")]
    video_condition = {}
    for cond, folder in layout:
        scan = os.path.join(annotated_data, folder, "A")
        if not os.path.isdir(scan):
            sys.exit(f"ERROR: folder not found: {scan}")
        for p in glob.glob(os.path.join(scan, "*.csv")):
            video_condition[strip_participant_suffix(
                os.path.basename(p))] = (cond, folder)
    out = []
    for v in sorted(os.listdir(irr_root)):
        vdir = os.path.join(irr_root, v)
        if not os.path.isdir(vdir) or v not in video_condition:
            continue
        if not (os.path.isdir(os.path.join(vdir, "A"))
                and os.path.isdir(os.path.join(vdir, "B"))):
            continue                      # solo video -> skip
        cond, folder = video_condition[v]
        out.append(dict(
            video=v, condition=cond,
            orig_left=find_original(
                os.path.join(annotated_data, folder, "A"), v),
            orig_right=find_original(
                os.path.join(annotated_data, folder, "B"), v),
            sec_left=find_second_coder(os.path.join(vdir, "A")),
            sec_right=find_second_coder(os.path.join(vdir, "B"))))
    if not out:
        sys.exit("ERROR: no pair videos found.")
    return out


def collect_pair_videos_repo(annotated_data, rerun_root):
    """Same as collect_pair_videos, but the second dataset is a
    re-annotation tree in REPOSITORY layout (visible_pair_Human/{A,B},
    invisible_pair_Human/{A,B}); e.g. the same coder's second pass."""
    layout = [("visible", "visible_pair_Human"),
              ("invisible", "invisible_pair_Human")]
    out = []
    for cond, folder in layout:
        scan = os.path.join(rerun_root, folder, "A")
        if not os.path.isdir(scan):
            continue
        for p in sorted(glob.glob(os.path.join(scan, "*.csv"))):
            v = strip_participant_suffix(os.path.basename(p))
            out.append(dict(
                video=v, condition=cond,
                orig_left=find_original(
                    os.path.join(annotated_data, folder, "A"), v),
                orig_right=find_original(
                    os.path.join(annotated_data, folder, "B"), v),
                sec_left=find_original(
                    os.path.join(rerun_root, folder, "A"), v),
                sec_right=find_original(
                    os.path.join(rerun_root, folder, "B"), v)))
    if not out:
        sys.exit(f"ERROR: no pair videos found under {rerun_root}.")
    return out


# ----------------------------------------------------------------- NLOPM
# [FLOAT-TOL PATCH 2026-09-07] tolerance for comparisons against delta.
# Event times are on a 0.1-s grid; in binary floating point the difference of
# two grid values is not exactly a multiple of 0.1 (e.g. 0.8 - 0.7 =
# 0.10000000000000009 > 0.1), so events exactly delta apart were matched or
# not depending on their float representation. 1e-6 s is far below the grid.
FLOAT_TOL = float(os.environ.get("NLOPM_FLOAT_TOL", "1e-6"))



def nlopm_matches(A, B, delta):
    """Nearest-with-look-ahead order-preserving one-to-one matching
    (STAR Methods reimplementation). A, B sorted 1-D arrays.
    Returns number of matches."""
    n, m = len(A), len(B)
    i, j = 0, 0
    n_match = 0
    while i < n and j < m:
        a = A[i]
        while j < m and B[j] < a - delta - FLOAT_TOL:  # [FLOAT-TOL PATCH]
            j += 1
        if j >= m:
            break
        # candidate events in B within [a - delta, a + delta]
        k = j
        cands = []
        while k < m and B[k] <= a + delta + FLOAT_TOL:  # [FLOAT-TOL PATCH]
            cands.append(k)
            k += 1
        if not cands:
            i += 1
            continue
        best = min(cands, key=lambda kk: (round(abs(B[kk] - a), 6), B[kk]))  # [FLOAT-TOL PATCH]
        # look-ahead: if the NEXT event in A is strictly closer to the
        # chosen candidate, reserve that candidate for the next event
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


def sync_rate(A, B, delta):
    na, nb = len(A), len(B)
    if na + nb == 0:
        return np.nan
    return 2.0 * nlopm_matches(A, B, delta) / (na + nb)


def zs_period(t_left, t_right, p0, p1, delta, n_perm, rng,
              null_shift="continuous", grid=0.1):  # [GRID-NULL PATCH]
    """Period-restricted zS: truncate both series to the overlapping
    interval intersected with [p0, p1]; circular-shift the right
    participant's events within that interval (STAR Methods)."""
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
    s_null = np.empty(n_perm)
    n_steps = max(int(round(T / grid)), 2)  # [GRID-NULL PATCH]
    for r in range(n_perm):
        if null_shift == "grid":  # [GRID-NULL PATCH] same construction as R2-4_zS_nullcheck.py
            u = rng.integers(1, n_steps) * grid
            Bs = tmin + (B - tmin + u) % T
            Bs = np.round(Bs / grid) * grid   # clean float error only
            Bs = np.unique(Bs)
        else:
            u = rng.uniform(0.0, T)
            Bs = np.sort(tmin + (B - tmin + u) % T)
        s_null[r] = 2.0 * nlopm_matches(A, Bs, delta) / (len(A) + len(B))
    sd = s_null.std(ddof=1)
    zs = (s_obs - s_null.mean()) / sd if sd > 0 else np.nan
    p = (np.sum(s_null >= s_obs) + 1) / (n_perm + 1)
    return dict(n_a=len(A), n_b=len(B), n_match=n_match,
                S=s_obs, zS=zs, p=p)


# ------------------------------------------------------------------- ICC
def icc_a1(x, y):
    d = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    d = d[np.all(np.isfinite(d), axis=1)]
    n, k = d.shape
    if n < 2:
        return np.nan
    mean_r = d.mean(axis=1)
    mean_c = d.mean(axis=0)
    grand = d.mean()
    msr = k * np.sum((mean_r - grand) ** 2) / (n - 1)
    msc = n * np.sum((mean_c - grand) ** 2) / (k - 1)
    sse = np.sum((d - mean_r[:, None] - mean_c[None, :] + grand) ** 2)
    mse = sse / ((n - 1) * (k - 1))
    denom = msr + (k - 1) * mse + k * (msc - mse) / n
    return (msr - mse) / denom if denom > 0 else np.nan


def pearson(x, y):
    d = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    d = d[np.all(np.isfinite(d), axis=1)]
    if len(d) < 2 or d[:, 0].std() == 0 or d[:, 1].std() == 0:
        return np.nan
    return float(np.corrcoef(d[:, 0], d[:, 1])[0, 1])


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description="Coder-robustness of period-resolved zS (R2-4)")
    ap.add_argument("--annotated_data", required=True)
    ap.add_argument("--irr_root",
                    help="second coder's tree (ICC_annotate layout)")
    ap.add_argument("--rerun_root",
                    help="re-annotation tree in repository layout "
                         "(e.g. same coder's second pass, ICC_result)")
    ap.add_argument("--out_prefix", default="R2-4_zS")
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--grid", type=float, default=0.1)
    ap.add_argument("--n_perm", type=int, default=10000)
    ap.add_argument("--periods", default="P1:0-20,P2:20-60,P3:60-180")
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--null", choices=["continuous", "grid"], default="continuous",
                    help="circular-shift offsets: continuous U(0,T) [original] or "
                         "integer multiples of the grid [letter, General note]")  # [GRID-NULL PATCH]
    args = ap.parse_args()

    periods = []
    for tok in args.periods.split(","):
        name, rng_s = tok.split(":")
        a, b = rng_s.split("-")
        periods.append((name, float(a), float(b)))

    if args.rerun_root and args.irr_root:
        sys.exit("Use either --irr_root or --rerun_root, not both")
    if args.rerun_root:
        videos = collect_pair_videos_repo(args.annotated_data,
                                          args.rerun_root)
    elif args.irr_root:
        videos = collect_pair_videos(args.annotated_data, args.irr_root)
    else:
        sys.exit("Provide --irr_root or --rerun_root")
    print(f"{len(videos)} pair video(s) found.")
    rng = np.random.default_rng(args.seed)

    rows = []
    for v in videos:
        series = dict(
            original=(read_event_series(v["orig_left"], args.grid),
                      read_event_series(v["orig_right"], args.grid)),
            second=(read_event_series(v["sec_left"], args.grid),
                    read_event_series(v["sec_right"], args.grid)))
        for coder, (tl, tr) in series.items():
            for name, p0, p1 in periods:
                res = zs_period(tl, tr, p0, p1, args.delta,
                                args.n_perm, rng,
                                null_shift=args.null, grid=args.grid)  # [GRID-NULL PATCH]
                rows.append(dict(video=v["video"], condition=v["condition"],
                                 coder=coder, period=name, **res))
                print(f"  {v['video']} {coder:8s} {name}: "
                      f"zS = {res['zS']:.2f}" if np.isfinite(res['zS'])
                      else f"  {v['video']} {coder:8s} {name}: zS = NA",
                      flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(f"{args.out_prefix}_per_video.csv", index=False)

    with open(f"{args.out_prefix}_summary.txt", "w") as fh:
        fh.write("=====================================================\n")
        fh.write("R2-4 coder-robustness of the period-resolved zS\n")
        fh.write("=====================================================\n\n")
        if args.rerun_root:
            fh.write("MODE: intra-rater  (coder 'second' = re-annotation "
                     f"tree '{args.rerun_root}', same coder)\n")
        else:
            fh.write("MODE: inter-rater  (coder 'second' = independent "
                     f"second coder tree '{args.irr_root}')\n")
        fh.write(f"Pair videos           : {len(videos)}\n")
        fh.write(f"delta (tolerance)     : {args.delta} s\n")
        fh.write(f"Permutations per cell : {args.n_perm}\n")
        fh.write(f"Null shift            : {args.null}"
                 + (f" (offsets = k * {args.grid} s, grid-aligned)" if args.null == "grid"
                    else " (offsets ~ U(0,T), original convention)") + "\n")  # [GRID-NULL PATCH]
        fh.write(f"Periods               : {args.periods}\n\n")

        for name, _, _ in periods:
            d = df[df["period"] == name].pivot_table(
                index=["video", "condition"], columns="coder", values="zS",
                aggfunc="first").reset_index()
            fh.write(f"--- {name} ---\n")
            fh.write(f"  ICC(A,1) original vs second coder = "
                     f"{icc_a1(d['original'], d['second']):.3f}   "
                     f"Pearson r = {pearson(d['original'], d['second']):.3f}"
                     "\n")
            for cond in ("visible", "invisible"):
                dc = d[d["condition"] == cond]
                fh.write(f"  {cond:9s} (n = {len(dc)}): mean zS original = "
                         f"{dc['original'].mean():.2f}, second = "
                         f"{dc['second'].mean():.2f}\n")
            dv = d[d["condition"] == "visible"]
            di = d[d["condition"] == "invisible"]
            fh.write(f"  visible - invisible difference: original = "
                     f"{dv['original'].mean() - di['original'].mean():+.2f},"
                     f" second coder = "
                     f"{dv['second'].mean() - di['second'].mean():+.2f}\n")
            for _, r in d.iterrows():
                fh.write(f"    {r['video']:22s} {r['condition']:9s} "
                         f"original = {r['original']:6.2f}   "
                         f"second = {r['second']:6.2f}\n")
            fh.write("\n")

        fh.write("Note: zS computed per STAR Methods (NLOPM matching,\n"
                 "S = 2*n_match/(n_a+n_b), period-restricted circular-shift\n"
                 "null); independent reimplementation for the reliability\n"
                 "check, not the original analysis code.\n")

    print(f"Saved: {args.out_prefix}_per_video.csv, "
          f"{args.out_prefix}_summary.txt")


if __name__ == "__main__":
    main()
