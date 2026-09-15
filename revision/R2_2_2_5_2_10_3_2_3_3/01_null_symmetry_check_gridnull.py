#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 1 of the R2-10/R2-5/R2-2 branch (R2-10): Empirical check that circular-shift nulls are equivalent
whether A or B is shifted.

Background:
  The manuscript states "Owing to symmetry, shifting A or B yields equivalent
  results". Reviewer 2 (R2-10) notes this holds only if the null distribution
  is symmetric and asks for empirical confirmation.

Method (per pair x period, same conventions as 02_Micro_Analysis_NLOPM_periods.py):
  - S_obs from NLOPM matching restricted to the period.
  - Null 1 ("shiftA"): circularly shift A within the period, B fixed
      (= implementation used in the R1-2 period pipeline).
  - Null 2 ("shiftB"): circularly shift B within the period, A fixed
      (= procedure as described in the manuscript Methods).
  - Null 3 ("shiftA_ref", optional): a second shiftA null with an independent
      seed. |z_shiftA - z_shiftA_ref| quantifies pure Monte Carlo error and
      serves as the benchmark against which |z_shiftA - z_shiftB| is judged.
  Each null: n_perm replicates -> null mean, SD, sample skewness, z_S,
  right/left-tail Monte Carlo p.

Outputs (in --outdir):
  null_symmetry_{delta}_{tag}.csv  : one row per pair x period
  null_symmetry_master.csv         : concat of all null_symmetry_*.csv in outdir
  null_symmetry_summary.txt        : agreement metrics computed from the master
      (Pearson r of z_S, Bland-Altman mean diff and 95% limits of agreement,
       max |dz|, max |dp|, null skewness summary, KS distance between nulls,
       and the Monte Carlo error benchmark) -- recomputed on every run, so run
      the script for both conditions into the SAME outdir, then read the txt.

Usage example (run once per condition, same outdir):
  python 01_null_symmetry_check.py \
    --indir ".../Annotated_Data/visible_pair_Human" --tag visible \
    --periods "P1=0:20,P2=20:60,P3=60:180" \
    --delta 0.1 --perms 10000 --seed 0 --outdir out_symmetry
  python 01_null_symmetry_check.py \
    --indir ".../Annotated_Data/invisible_pair_Human" --tag invisible \
    --periods "P1=0:20,P2=20:60,P3=60:180" \
    --delta 0.1 --perms 10000 --seed 0 --outdir out_symmetry
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



# ------------------------ I/O utilities (identical to 02) ------------------------
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


# ------------------------ Core matching (identical to 02) ------------------------
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


# ------------------------ Null generation ------------------------
def sample_skewness(x: np.ndarray) -> float:
    """Adjusted Fisher-Pearson sample skewness (as in scipy.stats.skew(bias=False))."""
    n = len(x)
    if n < 3:
        return np.nan
    m = np.mean(x)
    s = np.std(x, ddof=1)
    if s == 0:
        return np.nan
    g1 = np.mean(((x - m) / s) ** 3)
    return float(g1 * n * n / ((n - 1) * (n - 2)))


def ks_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Two-sample Kolmogorov-Smirnov statistic (descriptive; no p-value)."""
    x = np.sort(x)
    y = np.sort(y)
    grid = np.concatenate([x, y])
    cdf_x = np.searchsorted(x, grid, side="right") / len(x)
    cdf_y = np.searchsorted(y, grid, side="right") / len(y)
    return float(np.max(np.abs(cdf_x - cdf_y)))


def circular_shift_null(a: np.ndarray, b: np.ndarray, seg_start: float, seg_end: float,
                        delta: float, n_perm: int, seed: int, shift_side: str) -> np.ndarray:
    """Null S distribution: circularly shift ONE side within [seg_start, seg_end]."""
    T = seg_end - seg_start
    rng = np.random.default_rng(seed)
    null_stats = np.empty(n_perm, dtype=float)
    for p in range(n_perm):
        u = draw_circ_shift_offset(rng, T)  # [GRID-NULL PATCH]
        if shift_side == "A":
            a_shift = a + u
            a_shift = np.where(a_shift > seg_end, a_shift - T, a_shift)
            if DEFAULT_NULL_SHIFT == "grid":
                a_shift = np.round(a_shift, 1)  # [GRID-NULL PATCH]
            a_shift.sort()
            M0, _ = NLOPM_match(a_shift, b, delta)
            null_stats[p] = sync_ratio(len(a_shift), len(b), len(M0))
        else:
            b_shift = b + u
            b_shift = np.where(b_shift > seg_end, b_shift - T, b_shift)
            if DEFAULT_NULL_SHIFT == "grid":
                b_shift = np.round(b_shift, 1)  # [GRID-NULL PATCH]
            b_shift.sort()
            M0, _ = NLOPM_match(a, b_shift, delta)
            null_stats[p] = sync_ratio(len(a), len(b_shift), len(M0))
    return null_stats


def null_summary(S_obs: float, null_stats: np.ndarray, prefix: str) -> Dict[str, float]:
    n_perm = len(null_stats)
    mu = float(np.mean(null_stats))
    sd = float(np.std(null_stats, ddof=1))
    return {
        f"null_mean_{prefix}": mu,
        f"null_sd_{prefix}": sd,
        f"null_skew_{prefix}": sample_skewness(null_stats),
        f"z_S_{prefix}": (S_obs - mu) / sd if sd > 0 else np.nan,
        f"p_right_{prefix}": float((np.sum(null_stats >= S_obs) + 1.0) / (n_perm + 1.0)),
        f"p_left_{prefix}": float((np.sum(null_stats <= S_obs) + 1.0) / (n_perm + 1.0)),
    }


# ------------------------ Period parsing / pair discovery (identical to 02) ------------------------
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


# ------------------------ Agreement summary (from master) ------------------------
def write_summary(master_path: str, out_txt: str, mc_ref: bool) -> None:
    df = pd.read_csv(master_path)
    df = df[df["sufficient"] == 1].copy()
    df = df.dropna(subset=["z_S_shiftA", "z_S_shiftB"])
    lines = []
    lines.append("=" * 68)
    lines.append("R2-10 null-symmetry check: agreement between shift-A and shift-B")
    lines.append(f"cells (pair x period, sufficient data): {len(df)}")
    lines.append(f"conditions in master: {sorted(df['condition'].unique().tolist())}")
    lines.append(f"delta values in master: {sorted(df['delta'].unique().tolist())}")
    lines.append("=" * 68)
    for delta, d in df.groupby("delta"):
        zA, zB = d["z_S_shiftA"].values, d["z_S_shiftB"].values
        diff = zA - zB
        r = float(np.corrcoef(zA, zB)[0, 1]) if len(d) > 2 else np.nan
        lines.append(f"\n--- delta = {delta} (n cells = {len(d)}) ---")
        lines.append(f"Pearson r(z_shiftA, z_shiftB)      : {r:.4f}")
        lines.append(f"Bland-Altman mean diff (A - B)     : {np.mean(diff):+.4f}")
        lines.append(f"Bland-Altman SD of diff            : {np.std(diff, ddof=1):.4f}")
        lo = np.mean(diff) - 1.96 * np.std(diff, ddof=1)
        hi = np.mean(diff) + 1.96 * np.std(diff, ddof=1)
        lines.append(f"Bland-Altman 95% LoA               : [{lo:+.4f}, {hi:+.4f}]")
        lines.append(f"max |z_shiftA - z_shiftB|          : {np.max(np.abs(diff)):.4f}")
        dp = np.abs(d["p_right_shiftA"].values - d["p_right_shiftB"].values)
        lines.append(f"max |p_right_shiftA - p_right_shiftB|: {np.max(dp):.4f}")
        lines.append(f"mean KS distance between nulls     : {np.mean(d['ks_AB'].values):.4f}"
                     f"  (max {np.max(d['ks_AB'].values):.4f})")
        lines.append("null skewness  shiftA: mean {:+.3f}  range [{:+.3f}, {:+.3f}]".format(
            np.mean(d["null_skew_shiftA"]), np.min(d["null_skew_shiftA"]), np.max(d["null_skew_shiftA"])))
        lines.append("null skewness  shiftB: mean {:+.3f}  range [{:+.3f}, {:+.3f}]".format(
            np.mean(d["null_skew_shiftB"]), np.min(d["null_skew_shiftB"]), np.max(d["null_skew_shiftB"])))
        if mc_ref and "z_S_shiftA_ref" in d.columns and d["z_S_shiftA_ref"].notna().any():
            dref = np.abs(d["z_S_shiftA"].values - d["z_S_shiftA_ref"].values)
            lines.append("Monte Carlo error benchmark (same side, independent seed):")
            lines.append(f"  mean |z_shiftA - z_shiftA_ref|   : {np.mean(dref):.4f}  (max {np.max(dref):.4f})")
            lines.append(f"  mean |z_shiftA - z_shiftB|       : {np.mean(np.abs(diff)):.4f}"
                         f"  (max {np.max(np.abs(diff)):.4f})")
            lines.append("  -> A-vs-B discrepancy comparable to the same-side MC error")
            lines.append("     indicates empirical equivalence of shifting A or B.")
    txt = "\n".join(lines) + "\n"
    with open(out_txt, "w") as f:
        f.write(txt)
    print(txt)


# ------------------------ Main ------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", type=str, required=True)
    ap.add_argument("--tag", type=str, default=None, help="condition tag (e.g., visible)")
    ap.add_argument("--periods", type=str, required=True,
                    help='e.g. "P1=0:20,P2=20:60,P3=60:180" (s, relative to overlap start)')
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--min_events", type=int, default=5)
    ap.add_argument("--no_mc_ref", action="store_true",
                    help="skip the second shift-A null used as Monte Carlo error benchmark")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null_shift", type=str, default="grid", choices=["grid", "continuous"],
                    help="[GRID-NULL PATCH] 'grid' = 0.1-s integer-multiple shifts (corrected, default); "
                         "'continuous' = U(0,T) (previous behaviour)")
    ap.add_argument("--outdir", type=str, default="out_symmetry")
    args = ap.parse_args()
    global DEFAULT_NULL_SHIFT
    DEFAULT_NULL_SHIFT = args.null_shift


    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or regex.sub(r'\W+', '_', os.path.basename(os.path.abspath(args.indir)))
    periods = parse_periods(args.periods)
    pairs = auto_discover_pairs(args.indir)
    if not pairs:
        raise SystemExit(f"[ERROR] no A/B pairs in {args.indir}")
    mc_ref = not args.no_mc_ref
    print(f"[INFO] {len(pairs)} pairs, periods: {periods}, perms={args.perms}, mc_ref={mc_ref}")

    rows = []
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
                   "sufficient": int(min(len(aa), len(bb)) >= args.min_events and T_seg > 0)}
            if T_seg <= 0 or len(aa) == 0 or len(bb) == 0:
                print(f"[WARN] {pid} {label}: empty segment; skipped")
                rows.append(row)
                continue
            M, _ = NLOPM_match(aa, bb, args.delta)
            S_obs = sync_ratio(len(aa), len(bb), len(M))
            row.update({"n_matched": len(M), "S_obs": S_obs})
            # deterministic per-cell base seed (independent across pairs/periods/sides)
            base = abs(hash((pid, label, round(args.delta, 6), args.seed))) % (2 ** 31)
            null_A = circular_shift_null(aa, bb, seg_start, seg_end, args.delta,
                                         args.perms, seed=base + 1, shift_side="A")
            null_B = circular_shift_null(aa, bb, seg_start, seg_end, args.delta,
                                         args.perms, seed=base + 2, shift_side="B")
            row.update(null_summary(S_obs, null_A, "shiftA"))
            row.update(null_summary(S_obs, null_B, "shiftB"))
            row["ks_AB"] = ks_distance(null_A, null_B)
            if mc_ref:
                null_A2 = circular_shift_null(aa, bb, seg_start, seg_end, args.delta,
                                              args.perms, seed=base + 3, shift_side="A")
                row.update(null_summary(S_obs, null_A2, "shiftA_ref"))
            rows.append(row)
        print(f"[INFO] done {pid}")

    dstr = str(args.delta)
    out_csv = os.path.join(args.outdir, f"null_symmetry_{dstr}_{tag}.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"[DONE] {out_csv}")

    master = os.path.join(args.outdir, "null_symmetry_master.csv")
    paths = [p for p in sorted(glob.glob(os.path.join(args.outdir, "null_symmetry_*.csv")))
             if not p.endswith("null_symmetry_master.csv")]
    pd.concat([pd.read_csv(p) for p in paths if os.path.getsize(p) > 0],
              ignore_index=True).to_csv(master, index=False)
    print(f"[DONE] {master}")
    write_summary(master, os.path.join(args.outdir, "null_symmetry_summary.txt"), mc_ref)


if __name__ == "__main__":
    main()
