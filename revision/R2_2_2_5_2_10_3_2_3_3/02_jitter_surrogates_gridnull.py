#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2 of the R2-10/R2-5/R2-2 branch (R2-5): Stricter surrogate tests for micro-scale synchrony.

Background:
  Reviewer 2 (R2-5) notes that circular-shift surrogates preserve only the
  marginal event-rate structure of each series: a null rejection could in
  principle reflect slow covariation of chewing rates (rate matching) rather
  than fine-grained event-level coincidence. This script adds two stricter
  null models that PRESERVE slow rate covariation and therefore isolate
  fine-grained timing:

  Null models (per pair x period, same NLOPM conventions as
  02_Micro_Analysis_NLOPM_periods.py; A fixed, B surrogated):
    - "jitter" (primary): event-wise jitter b'_k = b_k + eta_k,
        eta_k ~ Uniform(-w, +w) i.i.d., reflected at the period boundaries.
        Destroys timing at scales <= w while preserving the slow rate
        envelope (w << period length). Primary w = 0.35 s (= half the mean
        inter-chew interval of 0.70 s); sensitivity w = 0.20, 0.50 s.
    - "isi" (secondary): within-period shuffle of inter-event intervals
        (first event anchored). Preserves the ISI distribution (rhythm)
        and event count; destroys ISI ordering and slow envelope alignment.
    - "circ" (reference): circular shift of B, as in the manuscript.
        Included so that all three nulls are computed with identical seeds,
        matching conventions, and output format in a single run.

  For each null: n_perm surrogates -> null mean/SD, z_S, right/left Monte
  Carlo p (add-one rule), as in the main pipeline.

Group-level summary (recomputed from the master on every run):
  per delta x null_type x w x condition x period:
    n pairs, mean z_S, SD, one-sample t vs 0 (two-sided p; pure-Python
    Student-t via regularized incomplete beta -- no scipy dependency),
    Holm correction across the 6 condition x period cells WITHIN each
    delta x null_type x w family (same policy as the R1-2 analyses),
    and descriptive prevalence: n pairs with p_right < .05 (uncorrected).

Outputs (in --outdir):
  jitter_null_{delta}_{tag}.csv   : one row per pair x period x null config
  jitter_null_master.csv          : concat of all jitter_null_*.csv in outdir
  jitter_null_group_summary.txt   : group-level table (from master)

Usage example (run once per condition, same outdir):
  python 02_jitter_surrogates.py \
    --indir ".../visible_pair_Human" --tag visible \
    --periods "P1=0:20,P2=20:60,P3=60:180" \
    --delta 0.1 --w "0.35,0.2,0.5" --nulls "jitter,isi,circ" \
    --perms 10000 --seed 0 --outdir out_jitter
  python 02_jitter_surrogates.py \
    --indir ".../invisible_pair_Human" --tag invisible \
    --periods "P1=0:20,P2=20:60,P3=60:180" \
    --delta 0.1 --w "0.35,0.2,0.5" --nulls "jitter,isi,circ" \
    --perms 10000 --seed 0 --outdir out_jitter
"""
import argparse
import glob
import math
import os
import re as regex
from typing import Dict, List, Optional, Tuple

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



# ------------------------ I/O utilities (identical to 01/02) ------------------------
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


# ------------------------ Core matching (identical to 01/02) ------------------------
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


# ------------------------ Null generators ------------------------
def jitter_surrogate_null(a: np.ndarray, b: np.ndarray, seg_start: float, seg_end: float,
                          delta: float, w: float, n_perm: int, seed: int) -> np.ndarray:
    """Event-wise jitter of B (A fixed), reflected at period boundaries."""
    rng = np.random.default_rng(seed)
    null_stats = np.empty(n_perm, dtype=float)
    for p in range(n_perm):
        if DEFAULT_NULL_SHIFT == "grid":
            K = int(np.floor(w / ANNOTATION_GRID + 1e-9))  # [GRID-NULL PATCH]
            b2 = b + ANNOTATION_GRID * rng.integers(-K, K + 1, size=len(b))
        else:
            b2 = b + rng.uniform(-w, w, size=len(b))
        # reflect at segment boundaries (w << T_seg, one reflection suffices)
        b2 = np.where(b2 < seg_start, 2.0 * seg_start - b2, b2)
        b2 = np.where(b2 > seg_end, 2.0 * seg_end - b2, b2)
        if DEFAULT_NULL_SHIFT == "grid":
            b2 = np.round(b2, 1)  # [GRID-NULL PATCH]
        b2.sort()
        M0, _ = NLOPM_match(a, b2, delta)
        null_stats[p] = sync_ratio(len(a), len(b2), len(M0))
    return null_stats


def isi_shuffle_null(a: np.ndarray, b: np.ndarray,
                     delta: float, n_perm: int, seed: int) -> Optional[np.ndarray]:
    """Within-period shuffle of B's inter-event intervals (first event anchored)."""
    if len(b) < 3:
        return None
    isi = np.diff(b)
    rng = np.random.default_rng(seed)
    null_stats = np.empty(n_perm, dtype=float)
    for p in range(n_perm):
        b2 = b[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(isi))])
        M0, _ = NLOPM_match(a, b2, delta)
        null_stats[p] = sync_ratio(len(a), len(b2), len(M0))
    return null_stats


def circular_shift_null(a: np.ndarray, b: np.ndarray, seg_start: float, seg_end: float,
                        delta: float, n_perm: int, seed: int) -> np.ndarray:
    """Circular shift of B within the period (manuscript procedure)."""
    T = seg_end - seg_start
    rng = np.random.default_rng(seed)
    null_stats = np.empty(n_perm, dtype=float)
    for p in range(n_perm):
        u = draw_circ_shift_offset(rng, T)  # [GRID-NULL PATCH]
        b2 = b + u
        b2 = np.where(b2 > seg_end, b2 - T, b2)
        if DEFAULT_NULL_SHIFT == "grid":
            b2 = np.round(b2, 1)  # [GRID-NULL PATCH]
        b2.sort()
        M0, _ = NLOPM_match(a, b2, delta)
        null_stats[p] = sync_ratio(len(a), len(b2), len(M0))
    return null_stats


def null_summary(S_obs: float, null_stats: np.ndarray) -> Dict[str, float]:
    n_perm = len(null_stats)
    mu = float(np.mean(null_stats))
    sd = float(np.std(null_stats, ddof=1))
    return {
        "null_mean": mu,
        "null_sd": sd,
        "z_S": (S_obs - mu) / sd if sd > 0 else np.nan,
        "p_right": float((np.sum(null_stats >= S_obs) + 1.0) / (n_perm + 1.0)),
        "p_left": float((np.sum(null_stats <= S_obs) + 1.0) / (n_perm + 1.0)),
    }


# ------------------------ Student-t p-value (pure Python, no scipy) ------------------------
def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function (Numerical Recipes)."""
    MAXIT, EPS, FPMIN = 300, 3e-12, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS:
            break
    return h


def betainc_reg(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_bt = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log(1.0 - x))
    bt = math.exp(ln_bt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_two_sided_p(t: float, df: float) -> float:
    """Two-sided p for Student t: p = I_x(df/2, 1/2), x = df / (df + t^2)."""
    if not np.isfinite(t) or df <= 0:
        return float("nan")
    return betainc_reg(df / 2.0, 0.5, df / (df + t * t))


def holm(pvals: List[float]) -> List[float]:
    """Holm step-down adjusted p-values (NaN-safe)."""
    idx = [i for i, p in enumerate(pvals) if np.isfinite(p)]
    m = len(idx)
    order = sorted(idx, key=lambda i: pvals[i])
    adj = [float("nan")] * len(pvals)
    running = 0.0
    for rank, i in enumerate(order):
        val = (m - rank) * pvals[i]
        running = max(running, val)
        adj[i] = min(1.0, running)
    return adj


# ------------------------ Period parsing / pair discovery (identical to 01/02) ------------------------
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


# ------------------------ Group-level summary (from master) ------------------------
def write_group_summary(master_path: str, out_txt: str) -> None:
    df = pd.read_csv(master_path)
    df = df[(df["sufficient"] == 1) & df["z_S"].notna()].copy()
    df["jitter_w"] = df["jitter_w"].fillna(-1.0)  # -1 = not applicable (isi/circ)
    lines = []
    lines.append("=" * 76)
    lines.append("R2-5 surrogate tests: group-level one-sample t of z_S vs 0")
    lines.append("(Holm across the 6 condition x period cells within each")
    lines.append(" delta x null_type x w family; prevalence = n pairs p_right < .05, unc.)")
    lines.append(f"conditions in master: {sorted(df['condition'].unique().tolist())}")
    lines.append(f"delta values in master: {sorted(df['delta'].unique().tolist())}")
    lines.append("=" * 76)
    for (delta, ntype, w), fam in df.groupby(["delta", "null_type", "jitter_w"]):
        wtxt = f", w={w:g}s" if w >= 0 else ""
        lines.append(f"\n--- delta = {delta}, null = {ntype}{wtxt} ---")
        lines.append(f"{'cond':<10} {'period':<6} {'n':>3} {'mean_zS':>8} {'sd':>7} "
                     f"{'t':>7} {'p':>9} {'p_holm':>9} {'prev':>6}")
        cells = []
        for (cond, per), cell in fam.groupby(["condition", "period"]):
            n = len(cell)
            mz = float(np.mean(cell["z_S"]))
            sd = float(np.std(cell["z_S"], ddof=1)) if n > 1 else float("nan")
            t = mz / (sd / math.sqrt(n)) if (n > 1 and sd > 0) else float("nan")
            p = t_two_sided_p(t, n - 1)
            prev = int(np.sum(cell["p_right"] < 0.05))
            cells.append({"cond": cond, "per": per, "n": n, "mz": mz, "sd": sd,
                          "t": t, "p": p, "prev": prev})
        adj = holm([c["p"] for c in cells])
        for c, ph in zip(cells, adj):
            lines.append(f"{c['cond']:<10} {c['per']:<6} {c['n']:>3} {c['mz']:>8.3f} "
                         f"{c['sd']:>7.3f} {c['t']:>7.2f} {c['p']:>9.4f} {ph:>9.4f} "
                         f"{c['prev']:>3d}/{c['n']:<2d}")
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
    ap.add_argument("--w", type=str, default="0.35,0.2,0.5",
                    help="comma-separated jitter half-widths in s (primary first)")
    ap.add_argument("--nulls", type=str, default="jitter,isi,circ",
                    help="comma-separated subset of {jitter,isi,circ}")
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--min_events", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null_shift", type=str, default="grid", choices=["grid", "continuous"],
                    help="[GRID-NULL PATCH] 'grid' = 0.1-s integer-multiple shifts (corrected, default); "
                         "'continuous' = U(0,T) (previous behaviour)")
    ap.add_argument("--outdir", type=str, default="out_jitter")
    args = ap.parse_args()
    global DEFAULT_NULL_SHIFT
    DEFAULT_NULL_SHIFT = args.null_shift


    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or regex.sub(r'\W+', '_', os.path.basename(os.path.abspath(args.indir)))
    periods = parse_periods(args.periods)
    w_list = [float(x) for x in args.w.split(",") if x.strip()]
    null_types = [x.strip() for x in args.nulls.split(",") if x.strip()]
    for nt in null_types:
        if nt not in ("jitter", "isi", "circ"):
            raise SystemExit(f"[ERROR] unknown null type: {nt}")
    pairs = auto_discover_pairs(args.indir)
    if not pairs:
        raise SystemExit(f"[ERROR] no A/B pairs in {args.indir}")
    print(f"[INFO] {len(pairs)} pairs, periods: {periods}, nulls={null_types}, "
          f"w={w_list}, perms={args.perms}")

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
            base_row = {"pair_id": pid, "condition": tag, "delta": args.delta,
                        "period": label, "t0_rel": t0, "t1_rel": t1, "T_seg": T_seg,
                        "n_A": len(aa), "n_B": len(bb),
                        "sufficient": int(min(len(aa), len(bb)) >= args.min_events and T_seg > 0)}
            if T_seg <= 0 or len(aa) == 0 or len(bb) == 0:
                print(f"[WARN] {pid} {label}: empty segment; skipped")
                rows.append(dict(base_row, null_type="none", jitter_w=np.nan))
                continue
            M, _ = NLOPM_match(aa, bb, args.delta)
            S_obs = sync_ratio(len(aa), len(bb), len(M))
            base = abs(hash((pid, label, round(args.delta, 6), args.seed))) % (2 ** 31)
            configs = []
            if "jitter" in null_types:
                configs += [("jitter", w) for w in w_list]
            if "isi" in null_types:
                configs.append(("isi", np.nan))
            if "circ" in null_types:
                configs.append(("circ", np.nan))
            for k, (ntype, w) in enumerate(configs):
                row = dict(base_row, null_type=ntype, jitter_w=w,
                           n_matched=len(M), S_obs=S_obs)
                seed_k = base + 101 + k
                if ntype == "jitter":
                    null_stats = jitter_surrogate_null(aa, bb, seg_start, seg_end,
                                                       args.delta, w, args.perms, seed_k)
                elif ntype == "isi":
                    null_stats = isi_shuffle_null(aa, bb, args.delta, args.perms, seed_k)
                else:
                    null_stats = circular_shift_null(aa, bb, seg_start, seg_end,
                                                     args.delta, args.perms, seed_k)
                if null_stats is None:
                    print(f"[WARN] {pid} {label} {ntype}: too few events; skipped")
                    rows.append(row)
                    continue
                row.update(null_summary(S_obs, null_stats))
                rows.append(row)
        print(f"[INFO] done {pid}")

    dstr = str(args.delta)
    out_csv = os.path.join(args.outdir, f"jitter_null_{dstr}_{tag}.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"[DONE] {out_csv}")

    master = os.path.join(args.outdir, "jitter_null_master.csv")
    paths = [p for p in sorted(glob.glob(os.path.join(args.outdir, "jitter_null_*.csv")))
             if not p.endswith("jitter_null_master.csv")]
    pd.concat([pd.read_csv(p) for p in paths if os.path.getsize(p) > 0],
              ignore_index=True).to_csv(master, index=False)
    print(f"[DONE] {master}")
    write_group_summary(master, os.path.join(args.outdir, "jitter_null_group_summary.txt"))


if __name__ == "__main__":
    main()
