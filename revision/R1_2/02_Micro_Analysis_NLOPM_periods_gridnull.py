#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2 (R1-2): Period-restricted micro-scale synchrony analysis.

Extends the original Micro_Analysis_NOLPM.py so that all statistics are computed
SEPARATELY within user-defined periods (e.g., non-steady onset ~0-20 s vs each
subsequent steady-state period), instead of pooling the full 3-min session.

Per pair x period:
  - NLOPM one-to-one matching restricted to the period -> S, n_matched
  - circular-shift permutation WITHIN the period (shift range = period length)
      -> z_S, right-tail p (synchrony > chance) and left-tail p (synchrony < chance;
         relevant for the below-zero dip period)
  - optional sliding-window episode metrics within the period -> STR (normalized by
    period length) and Max_episode  [use with caution for short periods]
  - long-format Delta-t table of matched events (for Step 5 latency analyses)

Period convention:
  --periods takes boundaries RELATIVE TO OVERLAP START of each pair, e.g.
    --periods "P1=0:20,P2=20:60,P3=60:180"
  Labels are optional: "0:20,20:60,60:180" auto-labels P1..Pn.
  Upper bounds are clipped to the actual overlap end of each pair.

Outputs (in --outdir):
  chew_sync_summary_periods_{delta}_{tag}.csv  : one row per pair x period
  dt_long_periods_{delta}_{tag}.csv            : one row per matched event
  chew_sync_summary_periods_master.csv         : concat of all summary files in outdir
  dt_long_periods_master.csv                   : concat of all dt files in outdir

Usage example:
  python 02_Micro_Analysis_NLOPM_periods.py \
    --indir ".../Annotated_Data/visible_pair_Human" --tag visible \
    --periods "P1=0:20,P2=20:60,P3=60:180" \
    --delta 0.1 --perms 10000 --episodes --winperms 300 \
    --win 5 --step 1 --episode_z 1.64 --episode_minmatch 3 --episode_minsec 3.0 \
    --outdir out_periods
"""
import argparse
import glob
import os
import re as regex
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# ====================== [GRID-NULL PATCH 2026-08-28] ======================
# Annotated event times lie on a 0.1-s grid. Continuous circular-shift
# offsets (u ~ U(0,T)) move the shifted series off the grid, undercounting
# expected matches and inflating z_S (3:2 at delta=0.1, 5:4 at 0.2, ...).
# Corrected default: offsets are integer multiples of the grid, applied to
# the series AS ANALYZED (sub-grid alignment offsets, if any, are preserved).
# Set --null_shift continuous to reproduce the previous behaviour exactly.
ANNOTATION_GRID = 0.1
DEFAULT_NULL_SHIFT = "grid"

# [PERIOD-ORIGIN PATCH 2026-09-08] "absolute" (default): periods t0-t1 are video
# time, intersected with the pair's overlap interval; "relative": seconds since
# the first shared chew (overlap start), which is what the code did before.
PERIOD_ORIGIN = os.environ.get("NLOPM_PERIOD_ORIGIN", "absolute")

def draw_circ_shift_offset(rng, T, null_shift=None, grid=ANNOTATION_GRID):
    """Draw one circular-shift offset ('grid' or 'continuous')."""
    if null_shift is None:
        null_shift = DEFAULT_NULL_SHIFT
    if null_shift == "continuous":
        return float(rng.uniform(0.0, T))
    K = int(round(T / grid))
    if K < 2:
        return float(rng.uniform(0.0, T))
    return grid * float(rng.integers(1, K))
# ==========================================================================



# ------------------------ I/O utilities (as in the original script) ------------------------
def _read_csv_smart(path: str) -> pd.DataFrame:
    for enc in ("utf-8-sig", "utf-8", "cp932", "shift_jis", "latin1"):
        try:
            return pd.read_csv(path, header=0, encoding=enc)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, header=0)


def load_times_from_csv(path: str) -> np.ndarray:
    """
    [EVENT-COLUMN PATCH 2026-09-07]
    Read chewing-event times from one annotation CSV (one participant).
      - Long format with a "t" column: use that column (unchanged behaviour).
      - Wide format (Trial, Chew_1, Chew_2, ...): use ONLY the columns whose
        name starts with "chew". The previous implementation numericised
        EVERY cell (df.stack / to_numpy().ravel()), so the Trial index
        (1, 2, ..., k) leaked into the event series as spurious events at
        t = 1 s, 2 s, ..., k s.
    Keep non-negative values only; sort and drop duplicates (as before).
    """
    df = _read_csv_smart(path)
    if "t" in df.columns:
        s = pd.to_numeric(df["t"], errors="coerce")
    else:
        cols = [c for c in df.columns if str(c).strip().lower().startswith("chew")]
        if not cols:
            # fallback: drop obvious index-like columns, keep the rest
            drop = {"trial", "bite", "index", "unnamed: 0", ""}
            cols = [c for c in df.columns if str(c).strip().lower() not in drop]
        s = pd.to_numeric(pd.Series(df[cols].to_numpy().ravel()), errors="coerce")
    arr = s.dropna().astype(float).values
    arr = arr[arr >= 0]
    arr = np.unique(np.sort(arr))
    return arr


def restrict_overlap(a: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float, float]:
    if len(a) == 0 or len(b) == 0:
        return np.array([]), np.array([]), 0.0, 0.0
    start = max(a[0], b[0])
    end = min(a[-1], b[-1])
    if end <= start:
        return np.array([]), np.array([]), 0.0, 0.0
    return a[(a >= start) & (a <= end)], b[(b >= start) & (b <= end)], float(start), float(end)


# ------------------------ Core matching (identical to the original) ------------------------
# [FLOAT-TOL PATCH 2026-09-07] tolerance for comparisons against delta.
# Event times are on a 0.1-s grid; in binary floating point the difference of
# two grid values is not exactly a multiple of 0.1 (e.g. 0.8 - 0.7 =
# 0.10000000000000009 > 0.1), so events exactly delta apart were matched or
# not depending on their float representation. 1e-6 s is far below the grid.
FLOAT_TOL = float(os.environ.get("NLOPM_FLOAT_TOL", "1e-6"))



def NLOPM_match(a: np.ndarray, b: np.ndarray, delta: float):
    i = j = 0
    m, n = len(a), len(b)
    matched_pairs: List[Tuple[float, float]] = []
    dt_list: List[float] = []
    while i < m and j < n:
        while j < n and b[j] < a[i] - delta - FLOAT_TOL:  # [FLOAT-TOL PATCH]
            j += 1
        if j >= n:
            break
        while i < m and a[i] < b[j] - delta - FLOAT_TOL:  # [FLOAT-TOL PATCH]
            i += 1
        if i >= m:
            break
        j_star = j
        while (j_star + 1) < n \
                and (b[j_star + 1] <= a[i] + delta + FLOAT_TOL) \
                and (abs(a[i] - b[j_star + 1]) <= abs(a[i] - b[j_star]) + FLOAT_TOL):
            j_star += 1
        if abs(a[i] - b[j_star]) <= delta + FLOAT_TOL:  # [FLOAT-TOL PATCH]
            matched_pairs.append((float(a[i]), float(b[j_star])))
            dt_list.append(float(a[i] - b[j_star]))
            i += 1
            j = j_star + 1
        else:
            if a[i] < b[j]:
                i += 1
            else:
                j += 1
    return matched_pairs, dt_list


def sync_ratio(n_a: int, n_b: int, n_matched: int) -> float:
    denom = n_a + n_b
    return 0.0 if denom == 0 else (2.0 * n_matched) / denom


# ------------------------ Segment-level permutation ------------------------
def segment_perm_stats(a: np.ndarray, b: np.ndarray, seg_start: float, seg_end: float,
                       delta: float, n_perm: int, seed: int) -> Dict[str, float]:
    """Circular-shift permutation restricted to [seg_start, seg_end].
    Returns S_obs, z_S, right/left-tail p, null mean/sd. A is shifted, B fixed."""
    T = seg_end - seg_start
    M, _ = NLOPM_match(a, b, delta)
    S_obs = sync_ratio(len(a), len(b), len(M))
    out = {"S_obs": float(S_obs), "n_matched": int(len(M)),
           "z_S": np.nan, "p_right": np.nan, "p_left": np.nan,
           "null_mean": np.nan, "null_sd": np.nan}
    if T <= 0 or n_perm <= 0 or (len(a) == 0 or len(b) == 0):
        return out
    rng = np.random.default_rng(seed)
    null_stats = np.empty(n_perm, dtype=float)
    for p in range(n_perm):
        u = draw_circ_shift_offset(rng, T)  # [GRID-NULL PATCH]
        a_shift = a + u
        a_shift = np.where(a_shift > seg_end, a_shift - T, a_shift)
        if DEFAULT_NULL_SHIFT == "grid":
            a_shift = np.round(a_shift, 1)  # [GRID-NULL PATCH]
        a_shift.sort()
        M0, _ = NLOPM_match(a_shift, b, delta)
        null_stats[p] = sync_ratio(len(a_shift), len(b), len(M0))
    mu = float(np.mean(null_stats))
    sd = float(np.std(null_stats, ddof=1))
    out["null_mean"], out["null_sd"] = mu, sd
    out["z_S"] = (S_obs - mu) / sd if sd > 0 else np.nan
    out["p_right"] = float((np.sum(null_stats >= S_obs) + 1.0) / (n_perm + 1.0))
    out["p_left"] = float((np.sum(null_stats <= S_obs) + 1.0) / (n_perm + 1.0))
    return out


# ------------------------ Windowed episodes within a segment ------------------------
def episodes_in_segment(a: np.ndarray, b: np.ndarray, seg_start: float, seg_end: float,
                        delta: float, win: float, step: float, n_perm: int,
                        z_thr: float, min_match: int, min_duration: float,
                        seed: int) -> Dict[str, float]:
    """Sliding-window z_S(t) with a within-segment circular-shift null; episode metrics.
    STR is normalized by the SEGMENT length. Returns NaN when the segment is too short."""
    T = seg_end - seg_start
    if T < 2 * win or len(a) == 0 or len(b) == 0:
        return {"STR": np.nan, "Max_episode": np.nan, "n_windows": 0, "n_episodes": np.nan}
    rng = np.random.default_rng(seed)
    centers = np.arange(seg_start + win / 2.0, seg_end - win / 2.0 + 1e-9, step)

    def window_S(aa_full: np.ndarray) -> np.ndarray:
        s_vals = np.empty(len(centers))
        nm_vals = np.empty(len(centers), dtype=int)
        for ci, c in enumerate(centers):
            lo, hi = c - win / 2.0, c + win / 2.0
            aa = aa_full[(aa_full >= lo) & (aa_full <= hi)]
            bb = b[(b >= lo) & (b <= hi)]
            M0, _ = NLOPM_match(aa, bb, delta)
            s_vals[ci] = sync_ratio(len(aa), len(bb), len(M0))
            nm_vals[ci] = len(M0)
        return s_vals, nm_vals

    S_obs, nM_obs = window_S(a)
    null_mat = np.empty((len(centers), n_perm))
    for p in range(n_perm):
        u = draw_circ_shift_offset(rng, T)  # [GRID-NULL PATCH]
        a_shift = a + u
        a_shift = np.where(a_shift > seg_end, a_shift - T, a_shift)
        if DEFAULT_NULL_SHIFT == "grid":
            a_shift = np.round(a_shift, 1)  # [GRID-NULL PATCH]
        a_shift.sort()
        null_mat[:, p] = window_S(a_shift)[0]
    mu = np.nanmean(null_mat, axis=1)
    sd = np.nanstd(null_mat, axis=1, ddof=1)
    z = (S_obs - mu) / np.where(sd > 0, sd, np.nan)

    mask = (z >= z_thr) & (nM_obs >= min_match)
    episodes = []
    in_ep, ep_start = False, None
    for i, flag in enumerate(mask):
        if flag and not in_ep:
            in_ep, ep_start = True, centers[i]
        elif (not flag) and in_ep:
            in_ep = False
            if (centers[i] - ep_start) >= min_duration:
                episodes.append((ep_start, centers[i]))
    if in_ep and (centers[-1] - ep_start) >= min_duration:
        episodes.append((ep_start, centers[-1]))

    total = sum(e - s for s, e in episodes) if episodes else 0.0
    return {
        "STR": total / T,
        "Max_episode": max((e - s for s, e in episodes), default=0.0),
        "n_windows": int(len(centers)),
        "n_episodes": len(episodes),
    }


# ------------------------ Period parsing / pair discovery ------------------------
def parse_periods(spec: str) -> List[Tuple[str, float, float]]:
    out = []
    for k, item in enumerate(spec.split(","), 1):
        item = item.strip()
        if "=" in item:
            label, rng = item.split("=", 1)
        else:
            label, rng = f"P{k}", item
        t0, t1 = (float(x) for x in rng.split(":"))
        if t1 <= t0:
            raise SystemExit(f"[ERROR] bad period {item}")
        out.append((label.strip(), t0, t1))
    return out


_SUFFIX_RE = regex.compile(r'(?i)(?:[_-])([AB])\.csv$')


def auto_discover_pairs(indir: str):
    """Discover A/B pairs under two supported layouts:
      (i)  flat:      indir/<base>_A.csv + indir/<base>_B.csv
      (ii) subfolder: indir/A/<base>_A.csv + indir/B/<base>_B.csv
           (files without the _A/_B suffix inside A/ and B/ also work,
            as long as base names match between the two subfolders)
    """
    buckets: Dict[str, Dict[str, str]] = {}

    def add_file(path: str, role_hint: str = None):
        fn = os.path.basename(path)
        m = _SUFFIX_RE.search(fn)
        if m:
            role = m.group(1).upper()
            base = os.path.splitext(_SUFFIX_RE.sub(".csv", fn))[0]
        elif role_hint:
            role = role_hint
            base = os.path.splitext(fn)[0]
        else:
            return
        buckets.setdefault(base, {})[role] = path

    sub_a, sub_b = os.path.join(indir, "A"), os.path.join(indir, "B")
    if os.path.isdir(sub_a) and os.path.isdir(sub_b):
        for fn in os.listdir(sub_a):
            if fn.lower().endswith(".csv"):
                add_file(os.path.join(sub_a, fn), role_hint="A")
        for fn in os.listdir(sub_b):
            if fn.lower().endswith(".csv"):
                add_file(os.path.join(sub_b, fn), role_hint="B")
    else:
        for fn in os.listdir(indir):
            if fn.lower().endswith(".csv"):
                add_file(os.path.join(indir, fn))
    return [(b, r["A"], r["B"]) for b, r in sorted(buckets.items()) if "A" in r and "B" in r]


# ------------------------ Main ------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", type=str, required=True)
    ap.add_argument("--tag", type=str, default=None, help="condition tag (e.g., visible)")
    ap.add_argument("--periods", type=str, required=True,
                    help='e.g. "P1=0:20,P2=20:60,P3=60:180" (s, relative to overlap start)')
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--perms", type=int, default=10000, help="segment-level permutations")
    ap.add_argument("--min_events", type=int, default=5,
                    help="min events per participant per period; below -> sufficient=0")
    ap.add_argument("--episodes", action="store_true",
                    help="also compute within-period STR / Max_episode (slower)")
    ap.add_argument("--winperms", type=int, default=300, help="window-level permutations")
    ap.add_argument("--win", type=float, default=5.0)
    ap.add_argument("--step", type=float, default=1.0)
    ap.add_argument("--episode_z", type=float, default=1.64)
    ap.add_argument("--episode_minmatch", type=int, default=3)
    ap.add_argument("--episode_minsec", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null_shift", type=str, default="grid", choices=["grid", "continuous"],
                    help="[GRID-NULL PATCH] 'grid' = 0.1-s integer-multiple shifts (corrected, default); "
                         "'continuous' = U(0,T) (previous behaviour)")
    ap.add_argument("--outdir", type=str, default="out_periods")
    args = ap.parse_args()
    global DEFAULT_NULL_SHIFT
    DEFAULT_NULL_SHIFT = args.null_shift


    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or regex.sub(r'\W+', '_', os.path.basename(os.path.abspath(args.indir)))
    periods = parse_periods(args.periods)
    pairs = auto_discover_pairs(args.indir)
    if not pairs:
        raise SystemExit(f"[ERROR] no A/B pairs in {args.indir}")
    print(f"[INFO] {len(pairs)} pairs, periods: {periods}")

    summary_rows, dt_rows = [], []
    for pid, fa, fb in pairs:
        a = load_times_from_csv(fa)
        b = load_times_from_csv(fb)
        a2, b2, start, end = restrict_overlap(a, b)
        for label, t0, t1 in periods:
            if PERIOD_ORIGIN == 'absolute':  # [PERIOD-ORIGIN PATCH]
                seg_start = max(float(t0), start)
                seg_end = min(float(t1), end)
            else:
                seg_start = start + t0
                seg_end = min(start + t1, end)
            T_seg = seg_end - seg_start
            aa = a2[(a2 >= seg_start) & (a2 <= seg_end)]
            bb = b2[(b2 >= seg_start) & (b2 <= seg_end)]
            row = {"pair_id": pid, "condition": tag, "delta": args.delta,
                   "period": label, "t0_rel": t0, "t1_rel": t1, "T_seg": T_seg,
                   "n_A": len(aa), "n_B": len(bb),
                   "sufficient": int(min(len(aa), len(bb)) >= args.min_events
                                     and T_seg > 0)}
            if T_seg <= 0:
                print(f"[WARN] {pid} {label}: empty segment (overlap ends before t0)")
                summary_rows.append({**row, "n_matched": 0, "S": np.nan, "z_S": np.nan,
                                     "p_right": np.nan, "p_left": np.nan,
                                     "null_mean": np.nan, "null_sd": np.nan,
                                     "STR": np.nan, "Max_episode": np.nan})
                continue
            st = segment_perm_stats(aa, bb, seg_start, seg_end, args.delta,
                                    args.perms, seed=args.seed)
            row.update({"n_matched": st["n_matched"], "S": st["S_obs"], "z_S": st["z_S"],
                        "p_right": st["p_right"], "p_left": st["p_left"],
                        "null_mean": st["null_mean"], "null_sd": st["null_sd"]})
            if args.episodes:
                ep = episodes_in_segment(aa, bb, seg_start, seg_end, args.delta,
                                         args.win, args.step, args.winperms,
                                         args.episode_z, args.episode_minmatch,
                                         args.episode_minsec, seed=args.seed)
                row.update({"STR": ep["STR"], "Max_episode": ep["Max_episode"],
                            "n_windows": ep["n_windows"], "n_episodes": ep["n_episodes"]})
            else:
                row.update({"STR": np.nan, "Max_episode": np.nan})
            summary_rows.append(row)

            M, dts = NLOPM_match(aa, bb, args.delta)
            for (ta, tb), d in zip(M, dts):
                dt_rows.append({"pair_id": pid, "condition": tag, "delta": args.delta,
                                "period": label, "a_time": ta, "b_time": tb, "dt": d,
                                "a_time_rel": ta - start, "b_time_rel": tb - start})
        print(f"[INFO] done {pid}")

    dstr = str(args.delta)
    sum_path = os.path.join(args.outdir, f"chew_sync_summary_periods_{dstr}_{tag}.csv")
    dt_path = os.path.join(args.outdir, f"dt_long_periods_{dstr}_{tag}.csv")
    pd.DataFrame(summary_rows).to_csv(sum_path, index=False)
    pd.DataFrame(dt_rows).to_csv(dt_path, index=False)
    print(f"[DONE] {sum_path}\n[DONE] {dt_path}")

    # master merges (concat everything currently in outdir)
    for pattern, master in [("chew_sync_summary_periods_*.csv", "chew_sync_summary_periods_master.csv"),
                            ("dt_long_periods_*.csv", "dt_long_periods_master.csv")]:
        paths = [p for p in sorted(glob.glob(os.path.join(args.outdir, pattern)))
                 if not p.endswith(master)]
        if len(paths) >= 1:
            pd.concat([pd.read_csv(p) for p in paths if os.path.getsize(p) > 0],
                      ignore_index=True).to_csv(os.path.join(args.outdir, master), index=False)
            print(f"[DONE] {os.path.join(args.outdir, master)}")


if __name__ == "__main__":
    main()
