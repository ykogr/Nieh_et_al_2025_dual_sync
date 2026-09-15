#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 3 of the R2-10/R2-5/R2-2 branch (R2-2): Formal comparison of zero-lag
alignment against a minimal reactive model.

Background:
  Reviewer 2 (R2-2) notes that near-zero-lag synchrony does not, by itself,
  distinguish prediction from very fast reactive coupling, and asks for a
  formal comparison against a minimal reactive model (e.g., a 150-250 ms
  shift). Because delta = 0.1 s limits the support of Delta-t to [-0.1, +0.1]
  and event times are annotated at 0.1 s resolution, a continuous
  distribution fit is inappropriate. Instead we use a SHIFTED-ALIGNMENT
  test:

  For each pair x period we compute z_S (vs. the within-period circular-
  shift null, identical conventions as 02_Micro_Analysis_NLOPM_periods.py)
  under 1 + 2*len(taus) alignments:
    - "zero"   : the observed alignment (tau = 0)
    - "B+tau"  : B' = B + tau  (re-alignment implied by "A reacts to B with
                 latency tau"), tau in {0.15, 0.20, 0.25} s by default
    - "A+tau"  : A' = A + tau  (re-alignment implied by "B reacts to A"),
                 same taus (leadership is not known a priori, so both
                 directions are tested)
  Prediction: if a minimal reactive model with latency tau* holds, z_S is
  maximized near the corresponding shifted alignment; if synchrony is
  predictive (zero-lag), z_S(zero) > max over all shifted alignments.

  Paired statistic per pair x period: d = z_S(zero) - max_{tau,dir} z_S.
  Note this is CONSERVATIVE for the zero-lag hypothesis: the max over 6
  noisy re-alignments is upward-biased, which favors the reactive side.
  Group inference: Wilcoxon signed-rank per condition x period cell
  (exact distribution when |d| has no ties, normal approximation with tie
  correction otherwise), Holm across the cells within each delta, plus a
  pooled test (descriptive) across cells.

  Shifted events that leave the period are dropped (no wrap); the loss is
  <= tau / T_seg (<= 1.25% for P1, <= 0.21% for P3), negligible.

Secondary (band-density excess) test, delta = 0.2 s only:
  If reactive following at 150-250 ms were present, matched latencies
  should concentrate in |Delta-t| in [0.15, 0.20] s (with 0.1 s annotation
  resolution, effectively |Delta-t| = 0.2). Statistic per pair x period:
    D_band = 2 * #{matched pairs with band_lo <= |Delta-t| <= band_hi} / (n_A + n_B),
  compared with the circular-shift null of the same statistic (z_band,
  right-tail Monte Carlo p = excess). Group: one-sample t of z_band per
  cell with Holm across cells, plus prevalence. Runs automatically when
  --delta 0.2 is given (auto mode); signed band counts (Delta-t = a - b > 0
  vs < 0) are also recorded descriptively.

Outputs (in --outdir):
  shift_test_{delta}_{tag}.csv  : one row per pair x period x alignment
  shift_test_master.csv         : concat of all shift_test_*.csv in outdir
  band_test_{delta}_{tag}.csv   : one row per pair x period (delta = 0.2 runs)
  band_test_master.csv          : concat of all band_test_*.csv in outdir
  shift_test_group_summary.txt  : group-level tables (both tests)

Usage (run once per condition x delta, same outdir):
  python 03_reactive_shift_test.py \
    --indir ".../visible_pair_Human" --tag visible \
    --periods "P1=0:20,P2=20:60,P3=60:180" \
    --delta 0.1 --taus "0.15,0.2,0.25" \
    --perms 10000 --seed 0 --outdir out_shift
  python 03_reactive_shift_test.py \
    --indir ".../invisible_pair_Human" --tag invisible \
    --periods "P1=0:20,P2=20:60,P3=60:180" \
    --delta 0.1 --taus "0.15,0.2,0.25" \
    --perms 10000 --seed 0 --outdir out_shift
  (then the same two commands with --delta 0.2; the band test runs
   automatically for delta = 0.2)
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


def band_stat(n_a: int, n_b: int, dts: List[float], lo: float, hi: float) -> Tuple[float, int, int, int]:
    """Band-density statistic and signed counts. eps guards 0.1 s quantization."""
    eps = 1e-9
    n_band = n_pos = n_neg = 0
    for dt in dts:
        if lo - eps <= abs(dt) <= hi + eps:
            n_band += 1
            if dt > 0:
                n_pos += 1
            elif dt < 0:
                n_neg += 1
    denom = n_a + n_b
    return (0.0 if denom == 0 else 2.0 * n_band / denom), n_band, n_pos, n_neg


# ------------------------ Null generator (identical scheme to 02) ------------------------
def circular_shift_null(a: np.ndarray, b: np.ndarray, seg_start: float, seg_end: float,
                        delta: float, n_perm: int, seed: int,
                        band: Optional[Tuple[float, float]] = None):
    """Circular shift of B within the period (A fixed). Returns null S stats;
    if band is given, also returns null band-density stats."""
    T = seg_end - seg_start
    rng = np.random.default_rng(seed)
    null_S = np.empty(n_perm, dtype=float)
    null_D = np.empty(n_perm, dtype=float) if band is not None else None
    for p in range(n_perm):
        u = draw_circ_shift_offset(rng, T)  # [GRID-NULL PATCH]
        b2 = b + u
        b2 = np.where(b2 > seg_end, b2 - T, b2)
        if DEFAULT_NULL_SHIFT == "grid":
            b2 = np.round(b2, 6)  # [GRID-NULL PATCH] kill float epsilon; 6 dp PRESERVES tau sub-grid offsets
        b2.sort()
        M0, dts0 = NLOPM_match(a, b2, delta)
        null_S[p] = sync_ratio(len(a), len(b2), len(M0))
        if band is not None:
            null_D[p], _, _, _ = band_stat(len(a), len(b2), dts0, band[0], band[1])
    return null_S, null_D


def null_summary(obs: float, null_stats: np.ndarray, prefix: str = "") -> Dict[str, float]:
    n_perm = len(null_stats)
    mu = float(np.mean(null_stats))
    sd = float(np.std(null_stats, ddof=1))
    return {
        f"{prefix}null_mean": mu,
        f"{prefix}null_sd": sd,
        f"{prefix}z": (obs - mu) / sd if sd > 0 else np.nan,
        f"{prefix}p_right": float((np.sum(null_stats >= obs) + 1.0) / (n_perm + 1.0)),
        f"{prefix}p_left": float((np.sum(null_stats <= obs) + 1.0) / (n_perm + 1.0)),
    }


# ------------------------ Student-t p-value (pure Python, identical to 02) ------------------------
def _betacf(a: float, b: float, x: float) -> float:
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
    if not np.isfinite(t) or df <= 0:
        return float("nan")
    return betainc_reg(df / 2.0, 0.5, df / (df + t * t))


def holm(pvals: List[float]) -> List[float]:
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


# ------------------------ Wilcoxon signed-rank (exact when tie-free) ------------------------
def wilcoxon_signed_rank(diffs: List[float]) -> Dict[str, float]:
    """Two-sided Wilcoxon signed-rank test. Zeros dropped. Exact DP
    distribution when |d| are tie-free (the usual case for continuous z_S);
    normal approximation with tie correction otherwise."""
    d = np.asarray([x for x in diffs if np.isfinite(x) and x != 0.0], dtype=float)
    n = len(d)
    out = {"n_eff": n, "W_pos": float("nan"), "p": float("nan"), "exact": np.nan}
    if n < 5:
        return out
    absd = np.abs(d)
    order = np.argsort(absd)
    ranks = np.empty(n, dtype=float)
    sorted_abs = absd[order]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_abs[j + 1] == sorted_abs[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0  # midranks
        i = j + 1
    W_pos = float(np.sum(ranks[d > 0]))
    out["W_pos"] = W_pos
    has_ties = len(np.unique(absd)) < n
    if not has_ties:
        # exact null distribution of W_pos via DP over ranks 1..n
        M = n * (n + 1) // 2
        counts = np.zeros(M + 1, dtype=float)
        counts[0] = 1.0
        for r in range(1, n + 1):
            counts[r:] += counts[:-r].copy()
        total = counts.sum()  # = 2^n
        w = int(round(W_pos))
        cdf_le = counts[:w + 1].sum() / total
        cdf_ge = counts[w:].sum() / total
        out["p"] = float(min(1.0, 2.0 * min(cdf_le, cdf_ge)))
        out["exact"] = 1
    else:
        mu = n * (n + 1) / 4.0
        # tie correction on the variance
        _, tie_counts = np.unique(absd, return_counts=True)
        tie_term = float(np.sum(tie_counts ** 3 - tie_counts))
        var = n * (n + 1) * (2 * n + 1) / 24.0 - tie_term / 48.0
        z = (W_pos - mu) / math.sqrt(var) if var > 0 else float("nan")
        out["p"] = float(math.erfc(abs(z) / math.sqrt(2.0))) if np.isfinite(z) else float("nan")
        out["exact"] = 0
    return out


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


# ------------------------ Group-level summary ------------------------
def _fmt(x, w=8, d=3):
    return f"{x:>{w}.{d}f}" if (x is not None and np.isfinite(x)) else " " * (w - 3) + "nan"


def write_group_summary(shift_master: str, band_master: Optional[str], out_txt: str) -> None:
    lines = []
    lines.append("=" * 78)
    lines.append("R2-2 shifted-alignment test: zero-lag vs minimal reactive re-alignments")
    lines.append("(d = z_S(zero) - max over {B+tau, A+tau} z_S; d > 0 favors zero-lag.")
    lines.append(" Wilcoxon signed-rank two-sided; Holm across cells within each delta.")
    lines.append(" Note: max over re-alignments is upward-biased -> conservative for zero-lag.)")
    lines.append("=" * 78)
    df = pd.read_csv(shift_master)
    df = df[(df["sufficient"] == 1) & df["z_S"].notna()].copy()
    for delta, fam in df.groupby("delta"):
        aligns = sorted(fam["alignment"].unique(),
                        key=lambda s: (s != "zero", s))
        lines.append(f"\n--- delta = {delta}: mean z_S per alignment ---")
        header = f"{'cond':<10} {'period':<6} {'n':>3}" + "".join(f" {a:>8}" for a in aligns)
        lines.append(header)
        wide = fam.pivot_table(index=["condition", "period", "pair_id"],
                               columns="alignment", values="z_S")
        for (cond, per), cell in wide.groupby(level=["condition", "period"]):
            vals = "".join(f" {_fmt(float(cell[a].mean()))}" if a in cell else " " * 9
                           for a in aligns)
            lines.append(f"{cond:<10} {per:<6} {len(cell):>3}{vals}")
        lines.append(f"\n--- delta = {delta}: paired test zero vs best reactive ---")
        lines.append(f"{'cond':<10} {'period':<6} {'n':>3} {'median_d':>9} {'mean_d':>8} "
                     f"{'n_d>0':>6} {'n_zero_best':>11} {'W':>7} {'p':>9} {'p_holm':>9} {'ex':>3}")
        cells = []
        for (cond, per), cell in wide.groupby(level=["condition", "period"]):
            if "zero" not in cell.columns:
                continue
            react_cols = [a for a in cell.columns if a != "zero"]
            sub = cell.dropna(subset=["zero"])
            react_max = sub[react_cols].max(axis=1)
            dvec = (sub["zero"] - react_max).dropna()
            wt = wilcoxon_signed_rank(list(dvec.values))
            cells.append({"cond": cond, "per": per, "n": len(dvec),
                          "med": float(np.median(dvec)) if len(dvec) else float("nan"),
                          "mean": float(np.mean(dvec)) if len(dvec) else float("nan"),
                          "npos": int(np.sum(dvec > 0)),
                          "nzb": int(np.sum(sub[["zero"] + react_cols].idxmax(axis=1) == "zero")),
                          "W": wt["W_pos"], "p": wt["p"], "ex": wt["exact"],
                          "d": list(dvec.values)})
        adj = holm([c["p"] for c in cells])
        for c, ph in zip(cells, adj):
            lines.append(f"{c['cond']:<10} {c['per']:<6} {c['n']:>3} {_fmt(c['med'], 9)} "
                         f"{_fmt(c['mean'])} {c['npos']:>4}/{c['n']:<2} {c['nzb']:>8}/{c['n']:<2} "
                         f"{_fmt(c['W'], 7, 1)} {_fmt(c['p'], 9, 4)} {_fmt(ph, 9, 4)} {c['ex']!s:>3}")
        pooled = [x for c in cells for x in c["d"]]
        wt = wilcoxon_signed_rank(pooled)
        lines.append(f"{'pooled':<10} {'all':<6} {len(pooled):>3} "
                     f"{_fmt(float(np.median(pooled)) if pooled else float('nan'), 9)} "
                     f"{_fmt(float(np.mean(pooled)) if pooled else float('nan'))} "
                     f"{int(np.sum(np.array(pooled) > 0)):>4}/{len(pooled):<2} {'-':>11} "
                     f"{_fmt(wt['W_pos'], 7, 1)} {_fmt(wt['p'], 9, 4)} {'-':>9} {wt['exact']!s:>3}")
    if band_master and os.path.exists(band_master):
        bdf = pd.read_csv(band_master)
        bdf = bdf[(bdf["sufficient"] == 1) & bdf["band_z"].notna()].copy()
        if len(bdf):
            lines.append("\n" + "=" * 78)
            lines.append("R2-2 band-density excess test (reactive-latency band, |dt| in [lo, hi])")
            lines.append("(one-sample t of z_band vs 0 per cell; Holm across cells within delta;")
            lines.append(" prev = n pairs with right-tail Monte Carlo p < .05, unc.;")
            lines.append(" n_pos/n_neg = mean signed band counts, dt = a - b)")
            lines.append("=" * 78)
            for (delta, lo, hi), fam in bdf.groupby(["delta", "band_lo", "band_hi"]):
                lines.append(f"\n--- delta = {delta}, band = [{lo:g}, {hi:g}] s ---")
                lines.append(f"{'cond':<10} {'period':<6} {'n':>3} {'mean_z':>8} {'sd':>7} "
                             f"{'t':>7} {'p':>9} {'p_holm':>9} {'prev':>6} "
                             f"{'mean_n_band':>11} {'n_pos':>6} {'n_neg':>6}")
                cells = []
                for (cond, per), cell in fam.groupby(["condition", "period"]):
                    n = len(cell)
                    mz = float(np.mean(cell["band_z"]))
                    sd = float(np.std(cell["band_z"], ddof=1)) if n > 1 else float("nan")
                    t = mz / (sd / math.sqrt(n)) if (n > 1 and sd > 0) else float("nan")
                    cells.append({"cond": cond, "per": per, "n": n, "mz": mz, "sd": sd, "t": t,
                                  "p": t_two_sided_p(t, n - 1),
                                  "prev": int(np.sum(cell["band_p_right"] < 0.05)),
                                  "nb": float(np.mean(cell["n_band"])),
                                  "npos": float(np.mean(cell["n_band_pos"])),
                                  "nneg": float(np.mean(cell["n_band_neg"]))})
                adj = holm([c["p"] for c in cells])
                for c, ph in zip(cells, adj):
                    lines.append(f"{c['cond']:<10} {c['per']:<6} {c['n']:>3} {_fmt(c['mz'])} "
                                 f"{_fmt(c['sd'], 7)} {_fmt(c['t'], 7, 2)} {_fmt(c['p'], 9, 4)} "
                                 f"{_fmt(ph, 9, 4)} {c['prev']:>3}/{c['n']:<2} {_fmt(c['nb'], 11, 1)} "
                                 f"{_fmt(c['npos'], 6, 1)} {_fmt(c['nneg'], 6, 1)}")
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
    ap.add_argument("--taus", type=str, default="0.15,0.2,0.25",
                    help="comma-separated reactive latencies in s")
    ap.add_argument("--band", type=str, default="0.15:0.2",
                    help="band for the density-excess test, 'lo:hi' in s")
    ap.add_argument("--do_band", type=str, default="auto", choices=["auto", "on", "off"],
                    help="band test: auto = only when delta == 0.2")
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--min_events", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null_shift", type=str, default="grid", choices=["grid", "continuous"],
                    help="[GRID-NULL PATCH] 'grid' = 0.1-s integer-multiple shifts (corrected, default); "
                         "'continuous' = U(0,T) (previous behaviour)")
    ap.add_argument("--outdir", type=str, default="out_shift")
    args = ap.parse_args()
    global DEFAULT_NULL_SHIFT
    DEFAULT_NULL_SHIFT = args.null_shift


    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or regex.sub(r'\W+', '_', os.path.basename(os.path.abspath(args.indir)))
    periods = parse_periods(args.periods)
    taus = [float(x) for x in args.taus.split(",") if x.strip()]
    band_lo, band_hi = (float(x) for x in args.band.split(":"))
    do_band = (args.do_band == "on") or (args.do_band == "auto" and abs(args.delta - 0.2) < 1e-9)
    pairs = auto_discover_pairs(args.indir)
    if not pairs:
        raise SystemExit(f"[ERROR] no A/B pairs in {args.indir}")
    alignments = [("zero", "none", 0.0)] \
        + [(f"B+{t:g}", "B", t) for t in taus] \
        + [(f"A+{t:g}", "A", t) for t in taus]
    print(f"[INFO] {len(pairs)} pairs, periods: {periods}, delta={args.delta}, "
          f"taus={taus}, band_test={do_band}, perms={args.perms}")

    shift_rows, band_rows = [], []
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
                        "sufficient": int(min(len(aa), len(bb)) >= args.min_events and T_seg > 0)}
            if T_seg <= 0 or len(aa) == 0 or len(bb) == 0:
                print(f"[WARN] {pid} {label}: empty segment; skipped")
                shift_rows.append(dict(base_row, alignment="zero", z_S=np.nan))
                continue
            base = abs(hash((pid, label, round(args.delta, 6), args.seed))) % (2 ** 31)
            for k, (aname, side, tau) in enumerate(alignments):
                if side == "B":
                    bb_al = bb + tau
                    bb_al = bb_al[(bb_al >= seg_start) & (bb_al <= seg_end)]
                    aa_al = aa
                elif side == "A":
                    aa_al = aa + tau
                    aa_al = aa_al[(aa_al >= seg_start) & (aa_al <= seg_end)]
                    bb_al = bb
                else:
                    aa_al, bb_al = aa, bb
                row = dict(base_row, alignment=aname, tau=tau,
                           n_A=len(aa_al), n_B=len(bb_al))
                if min(len(aa_al), len(bb_al)) < args.min_events:
                    row["z_S"] = np.nan
                    shift_rows.append(row)
                    continue
                M, _ = NLOPM_match(aa_al, bb_al, args.delta)
                S_obs = sync_ratio(len(aa_al), len(bb_al), len(M))
                null_S, _ = circular_shift_null(aa_al, bb_al, seg_start, seg_end,
                                                args.delta, args.perms, base + 201 + k)
                s = null_summary(S_obs, null_S)
                row.update(n_matched=len(M), S_obs=S_obs,
                           null_mean=s["null_mean"], null_sd=s["null_sd"],
                           z_S=s["z"], p_right=s["p_right"], p_left=s["p_left"])
                shift_rows.append(row)
            if do_band:
                M, dts = NLOPM_match(aa, bb, args.delta)
                D_obs, n_band, n_pos, n_neg = band_stat(len(aa), len(bb), dts, band_lo, band_hi)
                _, null_D = circular_shift_null(aa, bb, seg_start, seg_end,
                                                args.delta, args.perms, base + 901,
                                                band=(band_lo, band_hi))
                s = null_summary(D_obs, null_D, prefix="band_")
                band_rows.append(dict(base_row, band_lo=band_lo, band_hi=band_hi,
                                      n_A=len(aa), n_B=len(bb), n_matched=len(M),
                                      n_band=n_band, n_band_pos=n_pos, n_band_neg=n_neg,
                                      D_obs=D_obs, band_null_mean=s["band_null_mean"],
                                      band_null_sd=s["band_null_sd"], band_z=s["band_z"],
                                      band_p_right=s["band_p_right"],
                                      band_p_left=s["band_p_left"]))
        print(f"[INFO] done {pid}")

    dstr = str(args.delta)
    out_csv = os.path.join(args.outdir, f"shift_test_{dstr}_{tag}.csv")
    pd.DataFrame(shift_rows).to_csv(out_csv, index=False)
    print(f"[DONE] {out_csv}")
    if band_rows:
        bcsv = os.path.join(args.outdir, f"band_test_{dstr}_{tag}.csv")
        pd.DataFrame(band_rows).to_csv(bcsv, index=False)
        print(f"[DONE] {bcsv}")

    shift_master = os.path.join(args.outdir, "shift_test_master.csv")
    paths = [p for p in sorted(glob.glob(os.path.join(args.outdir, "shift_test_*.csv")))
             if not p.endswith("shift_test_master.csv")]
    pd.concat([pd.read_csv(p) for p in paths if os.path.getsize(p) > 0],
              ignore_index=True).to_csv(shift_master, index=False)
    band_master = os.path.join(args.outdir, "band_test_master.csv")
    bpaths = [p for p in sorted(glob.glob(os.path.join(args.outdir, "band_test_*.csv")))
              if not p.endswith("band_test_master.csv")]
    if bpaths:
        pd.concat([pd.read_csv(p) for p in bpaths if os.path.getsize(p) > 0],
                  ignore_index=True).to_csv(band_master, index=False)
    else:
        band_master = None
    write_group_summary(shift_master, band_master,
                        os.path.join(args.outdir, "shift_test_group_summary.txt"))


if __name__ == "__main__":
    main()
