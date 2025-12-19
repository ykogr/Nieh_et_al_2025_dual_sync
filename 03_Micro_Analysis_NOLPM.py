#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose:
  Evaluate micro-scale synchrony from chewing event series of person A and B.
  - Event level: one-to-one NLOPM matching to compute Δt = a_i - b_j (positive if A lags behind B).
  - Window level: observe synchrony ratio S per sliding window; build a null distribution via circular shift;
    then compute z_S(t) and S_excess(t).
  - Episodes: detect episode intervals using a z-threshold, a minimum number of matches per window,
    and a minimum duration (seconds).
  - Whole-interval stats: compute whole-interval z_S, STR (episode time ratio), and Max_episode.

Outputs:
  1) chew_sync_timeline_{pair_id}_{tag}.csv
       Long-format timeline per pair. Columns include:
         t_center, S, S_null_mean, S_null_sd, z_S_t, S_excess, nA, nB, nM, episode, dt
       If nM=0, one row is emitted with dt=NA to keep the window present.
  2) chew_sync_summary_{tag}.csv
       One row per pair summary (n_A, n_B, n_matched, z_S, p_S, STR, Max_episode, ...).
  3) timelines_all_{tag}__wide.csv
       Wide-merged per-window metrics across pairs keyed by t_center.
       Note: dt is not included.

Arguments:
  --indir   : Path to the input directory containing CSV files (e.g., /Annotated_Data/...)
  --outdir  : Path to the output directory (e.g., /out)

Usage example:
  python Micro_Analysis_NOLPM.py --indir "...Annotated_Data/visible_pair_Human" \
    --delta 0.1 --t 180 --perms 10000 --win 5 --step 1 --episode_z 1.64 \
    --episode_minmatch 3 --episode_minsec 3.0 --outdir out
"""
import argparse
import os
import re as regex
import glob
from typing import List, Tuple, Dict
import numpy as np
import pandas as pd


# ------------------------ I/O utilities ------------------------
def _read_csv_smart(path: str) -> pd.DataFrame:
    """Read a CSV with simple encoding auto-detection."""
    for enc in ("utf-8-sig", "utf-8", "cp932", "shift_jis", "latin1"):
        try:
            return pd.read_csv(path, header=0, encoding=enc)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, header=0)

def load_times_from_csv(path: str) -> np.ndarray:
    """
    Read an event-time column from one CSV (one participant) and return np.ndarray[float].
    - Preferred column name: "t" (seconds)
    - Keep non-negative values only; sort and drop duplicates
    """
    df = _read_csv_smart(path)
    if "t" in df.columns:
        s = pd.to_numeric(df["t"], errors="coerce")
    else:
        s = pd.to_numeric(df.stack(dropna=False), errors="coerce")
    arr = s.dropna().astype(float).values
    arr = arr[arr >= 0]
    arr = np.unique(np.sort(arr))
    return arr

# ------------------------ Time-domain alignment ------------------------
def restrict_overlap(a: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Restrict to the overlapping interval between a and b. All analyses are limited to this interval.
    Returns: a2,b2,start,end
    """
    if len(a) == 0 or len(b) == 0:
        return np.array([]), np.array([]), 0.0, 0.0
    start = max(a[0], b[0])
    end = min(a[-1], b[-1])
    if end <= start:
        return np.array([]), np.array([]), 0.0, 0.0
    a2 = a[(a >= start) & (a <= end)]
    b2 = b[(b >= start) & (b <= end)]
    return a2, b2, float(start), float(end)

# ------------------------ Core Matching ------------------------
def NLOPM_match(a: np.ndarray, b: np.ndarray, delta: float) -> Tuple[List[Tuple[float, float]], List[float]]:
    """
    One-to-one, order-preserving (non-crossing) NLOPM matching.
      - Within |a[i]-b[j]|<=δ, connect the closest pair
      - If A lags behind B, then Δt = a_i - b_j > 0
    Returns:
      matched_pairs: [(ai, bj), ...] (time-ordered)
      dt_list: [ai-bj, ...]
    """
    i = j = 0
    m, n = len(a), len(b)
    matched_pairs: List[Tuple[float, float]] = []
    dt_list: List[float] = []
    while i < m and j < n:
        while j < n and b[j] < a[i] - delta:
            j += 1
        if j >= n:
            break
        while i < m and a[i] < b[j] - delta:
            i += 1
        if i >= m:
            break
        j_star = j
        while (j_star + 1) < n \
                and (b[j_star + 1] <= a[i] + delta) \
                and (abs(a[i] - b[j_star + 1]) <= abs(a[i] - b[j_star])):
            j_star += 1
        if abs(a[i] - b[j_star]) <= delta:
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
    """Synchrony ratio S = 2 * n_matched / (n_a + n_b)"""
    denom = n_a + n_b
    return 0.0 if denom == 0 else (2.0 * n_matched) / denom

# ------------------------ Null distribution (circular shift) ------------------------
def circ_shift_permutation(a: np.ndarray, b: np.ndarray, start: float, end: float, delta: float,
                           n_perm: int = 1000, stat: str = "S", random_state: int = 0) -> Dict[str, float]:
    """
    Build a null distribution by shifting A by a uniform random offset within overlap length T (with wrap-around).
    stat = "S" or "median_abs_dt"
    Returns: dict(p_value, stat_obs, stat_null_mean, stat_null_std)
    """
    T = end - start
    rng = np.random.default_rng(random_state)
    M, dt = NLOPM_match(a, b, delta)
    if stat == "S":
        stat_obs = sync_ratio(len(a), len(b), len(M))
    elif stat == "median_abs_dt":
        stat_obs = np.median(np.abs(dt)) if len(dt) > 0 else np.inf
    else:
        raise ValueError("Unknown stat")
    if T <= 0:
        return {"p_value": np.nan, "stat_obs": stat_obs, "stat_null_mean": np.nan, "stat_null_std": np.nan}
    null_stats: List[float] = []
    for _ in range(n_perm):
        u = rng.uniform(0, T)
        a_shift = a + u
        a_shift = np.where(a_shift > (start + T), a_shift - T, a_shift)
        a_shift.sort()
        M0, dt0 = NLOPM_match(a_shift, b, delta)
        if stat == "S":
            s0 = sync_ratio(len(a_shift), len(b), len(M0))
        else:
            s0 = np.median(np.abs(dt0)) if len(dt0) > 0 else np.inf
        null_stats.append(s0)
    null_arr = np.array(null_stats, dtype=float)
    if stat == "S":
        p = (np.sum(null_arr >= stat_obs) + 1.0) / (len(null_arr) + 1.0)
    else:
        p = (np.sum(null_arr <= stat_obs) + 1.0) / (len(null_arr) + 1.0)
    return {
        "p_value": float(p),
        "stat_obs": float(stat_obs),
        "stat_null_mean": float(np.mean(null_arr)),
        "stat_null_std": float(np.std(null_arr, ddof=1))
    }

# ------------------------ Window-level z and long-format Δt table ------------------------
def sliding_window_with_dt(a: np.ndarray, b: np.ndarray, start: float, end: float,
                           delta: float, win: float, step: float,
                           n_perm: int = 300, seed: int = 0) -> pd.DataFrame:
    """
    Observe S per window, and estimate E[S] and SD[S] via circular-shift permutations.
    Also expand Δt of matched pairs in each window into a long-format "dt" column.
    Returns columns: t_center, S, S_null_mean, S_null_sd, z_S_t, S_excess, nA, nB, nM, dt
    """
    rng = np.random.default_rng(seed)

    centers = np.arange(start + win/2.0, end - win/2.0 + 1e-9, step)
    obs_rows = []
    dt_per_window = []
    for c in centers:
        lo, hi = c - win/2.0, c + win/2.0
        ia0 = np.searchsorted(a, lo, side="left"); ia1 = np.searchsorted(a, hi, side="right")
        ib0 = np.searchsorted(b, lo, side="left"); ib1 = np.searchsorted(b, hi, side="right")
        aa = a[ia0:ia1]; bb = b[ib0:ib1]
        M, dt = NLOPM_match(aa, bb, delta)
        nA, nB, nM = len(aa), len(bb), len(M)
        S = sync_ratio(nA, nB, nM)
        obs_rows.append((c, S, nA, nB, nM))
        dt_per_window.append(dt)

    obs_df = pd.DataFrame(obs_rows, columns=["t_center", "S", "nA", "nB", "nM"])

    # Null distribution (scan all windows with the same shift value)
    T = end - start
    null_mat = np.empty((len(centers), n_perm), dtype=float)
    for p in range(n_perm):
        u = rng.uniform(0, T)
        a_shift = a + u
        a_shift = np.where(a_shift > (start + T), a_shift - T, a_shift)
        a_shift.sort()
        s_vals = []
        for c in centers:
            lo, hi = c - win/2.0, c + win/2.0
            ia0 = np.searchsorted(a_shift, lo, side="left"); ia1 = np.searchsorted(a_shift, hi, side="right")
            ib0 = np.searchsorted(b, lo, side="left"); ib1 = np.searchsorted(b, hi, side="right")
            aa = a_shift[ia0:ia1]; bb = b[ib0:ib1]
            M0, _ = NLOPM_match(aa, bb, delta)
            s_vals.append(sync_ratio(len(aa), len(bb), len(M0)))
        null_mat[:, p] = np.array(s_vals, dtype=float)

    mu = np.nanmean(null_mat, axis=1)
    sd = np.nanstd(null_mat, axis=1, ddof=1)
    obs_df["S_null_mean"] = mu
    obs_df["S_null_sd"]   = sd
    obs_df["z_S_t"]       = (obs_df["S"].to_numpy() - mu) / np.where(sd > 0, sd, np.nan)
    obs_df["S_excess"]    = obs_df["S"].to_numpy() - mu

    # Expand Δt into long format (nM rows per window; if nM=0, emit one row with dt=NA)
    long_rows = []
    for (_, row), dt_list in zip(obs_df.iterrows(), dt_per_window):
        if len(dt_list) == 0:
            long_rows.append({
                "t_center": row["t_center"], "S": row["S"],
                "S_null_mean": row["S_null_mean"], "S_null_sd": row["S_null_sd"],
                "z_S_t": row["z_S_t"], "S_excess": row["S_excess"],
                "nA": row["nA"], "nB": row["nB"], "nM": row["nM"],
                "dt": np.nan
            })
        else:
            for d in dt_list:
                long_rows.append({
                    "t_center": row["t_center"], "S": row["S"],
                    "S_null_mean": row["S_null_mean"], "S_null_sd": row["S_null_sd"],
                    "z_S_t": row["z_S_t"], "S_excess": row["S_excess"],
                    "nA": row["nA"], "nB": row["nB"], "nM": row["nM"],
                    "dt": float(d)
                })
    return pd.DataFrame(long_rows)

# ------------------------ Episode detection ------------------------
def detect_episodes_adaptive(df: pd.DataFrame,
                             z_thr: float = 1.64,   # ≈ one-sided 95%
                             min_match: int = 3,
                             min_duration: float = 3.0,  # seconds
                             step: float = 1.0) -> Tuple[pd.DataFrame, List[Tuple[float, float]]]:
    """Define an episode as a consecutive run of windows satisfying z and nM conditions."""
    core = (df.drop_duplicates(subset=["t_center"]).sort_values("t_center"))
    mask = (core["z_S_t"] >= z_thr) & (core["nM"] >= min_match)
    episodes: List[Tuple[float, float]] = []
    in_ep = False
    start = None
    for i, flag in enumerate(mask.to_numpy()):
        if flag and not in_ep:
            in_ep = True; start = core["t_center"].iat[i]
        elif (not flag) and in_ep:
            in_ep = False; end = core["t_center"].iat[i]
            if (end - start) >= min_duration:
                episodes.append((start, end))
    if in_ep:
        end = core["t_center"].iat[-1]
        if (end - start) >= min_duration:
            episodes.append((start, end))
    df_out = df.copy(); df_out["episode"] = 0
    for s, e in episodes:
        df_out.loc[(df_out["t_center"] >= s) & (df_out["t_center"] <= e), "episode"] = 1
    return df_out, episodes

# ------------------------ Auto-discover input A/B pairs ------------------------
_SUFFIX_RE = regex.compile(r'(?i)(?:[_-])([AB])\.csv$')
def split_base_and_role(filename: str):
    m = _SUFFIX_RE.search(filename)
    if not m:
        return None, None
    role = m.group(1).upper()
    base = _SUFFIX_RE.sub(".csv", filename)
    base = os.path.splitext(base)[0]
    return base, role

def auto_discover_pairs(indir: str):
    files = [f for f in os.listdir(indir) if f.lower().endswith(".csv")]
    buckets = {}
    for fn in files:
        base_role = split_base_and_role(fn)
        if base_role == (None, None):
            continue
        base, role = base_role
        d = buckets.setdefault(base, {})
        d[role] = os.path.join(indir, fn)
    pairs = []
    for base, roles in sorted(buckets.items()):
        if "A" in roles and "B" in roles:
            pairs.append((base, roles["A"], roles["B"]))
    return pairs

# ------------------------ Single-pair analysis ------------------------
def analyze_pair(file_a: str, file_b: str,
                 delta: float, T_total: float, perms: int, seed: int,
                 win: float, step: float,
                 episode_z: float, episode_minmatch: int, episode_minsec: float,
                 out_prefix: str):
    a = load_times_from_csv(file_a)
    b = load_times_from_csv(file_b)
    a2, b2, start, end = restrict_overlap(a, b)
    if T_total is not None and T_total > 0:
        end = min(end, start + float(T_total))
        a2 = a2[(a2 >= start) & (a2 <= end)]
        b2 = b2[(b2 >= start) & (b2 <= end)]

    # Whole-interval z_S (right-tailed)
    M, dt_global = NLOPM_match(a2, b2, delta=delta)
    S_global = sync_ratio(len(a2), len(b2), len(M))
    perm_S = circ_shift_permutation(a2, b2, start, end, delta=delta,
                                    n_perm=perms, stat="S", random_state=seed)
    null_mean_S = perm_S["stat_null_mean"]
    null_std_S  = perm_S["stat_null_std"]
    z_S_global  = (S_global - null_mean_S) / null_std_S if (null_std_S and null_std_S > 0) else np.nan

    # Window-level z and long-format Δt table
    df_win_long = sliding_window_with_dt(a2, b2, start, end, delta, win, step, n_perm=1000, seed=seed)
    df_win_long, episodes = detect_episodes_adaptive(
        df_win_long, z_thr=episode_z, min_match=episode_minmatch, min_duration=episode_minsec, step=step
    )

    # STR / Max_episode
    T_overlap  = float(end - start)
    total_sync = float(sum((e - s) for (s, e) in episodes)) if episodes else 0.0
    STR        = (total_sync / T_overlap) if T_overlap > 0 else np.nan
    Max_episode = float(max((e - s) for (s, e) in episodes)) if episodes else 0.0

    # Output: long-format timeline (includes dt)
    csv_path = out_prefix + "_timeline.csv"
    df_win_long.to_csv(csv_path, index=False)

    return {
        "n_A": int(len(a2)), "n_B": int(len(b2)),
        "n_matched": int(len(M)),
        "z_S": float(z_S_global),
        "p_S": float(perm_S["p_value"]),
        "STR": float(STR),
        "Max_episode": float(Max_episode),
        "T_overlap": float(T_overlap),
        "timeline_csv": csv_path,
    }

# ------------------------ Merging / summary outputs ------------------------
def merge_timelines(rows: list, outdir: str, tag: str, filename: str = None):
    """
    Merge per-pair timeline CSVs into a wide table keyed by t_center.
    Note: dt is excluded because it increases rows; only window-series metrics are kept.
    """
    metrics = ["z_S_t", "S_excess", "nA", "nB", "nM", "episode"]  # dt is excluded

    all_times = None
    parts = []
    for r in rows:
        p = r.get("timeline_csv")
        if not p or not os.path.exists(p):
            continue
        pid = str(r.get("pair_id"))
        sid = pid if pid else "id"
        df = pd.read_csv(p)
        core = df.drop_duplicates(subset=["t_center"]).copy()
        if "t_center" not in core.columns:
            continue
        kept = [m for m in metrics if m in core.columns]
        sub = core[["t_center"] + kept].copy()
        sub.columns = ["t_center"] + [f"{m}__{sid}" for m in kept]
        parts.append(sub)
        all_times = sub["t_center"].values if all_times is None else np.union1d(all_times, sub["t_center"].values)

    if not parts:
        print("[WARN] No timelines to merge for wide export.")
        return None

    merged = pd.DataFrame({"t_center": all_times})
    for sub in parts:
        merged = merged.merge(sub, on="t_center", how="left")

    cols = ["t_center"]
    for m in metrics:
        cols.extend([c for c in merged.columns if c.startswith(f"{m}__")])
    merged = merged[[c for c in cols if c in merged.columns]]

    tag_safe = regex.sub(r'\W+', '_', tag)
    out_name = filename or f"timelines_all_{tag_safe}__wide.csv"
    timelines_path = os.path.join(outdir, out_name)
    merged.to_csv(timelines_path, index=False)
    print(f"[DONE] Saved wide timelines to {timelines_path}")
    return timelines_path

def merge_all_under(outdir: str, pattern: str, master_name: str):
    """Concatenate CSV files matching outdir/pattern to build a master CSV."""
    paths = sorted(glob.glob(os.path.join(outdir, pattern)))
    out_master = os.path.join(outdir, master_name)
    paths = [p for p in paths if os.path.abspath(p) != os.path.abspath(out_master)]
    if len(paths) <= 1:
        print(f"[INFO] Skip creating {master_name}: {len(paths)} source file(s).")
        return None
    frames = [pd.read_csv(p) for p in paths if os.path.getsize(p) > 0]
    if not frames:
        print(f"[INFO] Skip creating {master_name}: no non-empty sources.")
        return None
    pd.concat(frames, ignore_index=True).to_csv(out_master, index=False)
    print(f"[DONE] Saved master table to {out_master}")
    return out_master

def infer_tag_from_indir(indir: str, fallback: str = "run") -> str:
    base = os.path.basename(os.path.abspath(indir))
    return regex.sub(r'\W+', '_', base) if base else fallback

def write_rows(rows, outdir, tag_from_cli, indir_for_tag):
    """Write per-pair summary rows to CSV, then update wide and master tables."""
    tag = tag_from_cli or infer_tag_from_indir(indir_for_tag or outdir)
    for r in rows:
        r.setdefault("condition", tag)

    tag_safe = regex.sub(r'\W+', '_', tag)
    out_csv = os.path.join(outdir, f"chew_sync_summary_{tag_safe}.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"[DONE] Saved summary to {out_csv}")

    _ = merge_timelines(rows, outdir, tag)  # Exclude dt in wide table
    merge_all_under(outdir, "chew_sync_summary_*.csv", "chew_sync_summary_master.csv")
    merge_all_under(outdir, "timelines_all_*.csv", "timelines_all_master.csv")

# ------------------------ CLI ------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=str, help="CSV for person A")
    ap.add_argument("--b", type=str, help="CSV for person B")
    ap.add_argument("--pairs", type=str, help="CSV mapping with columns pair_id,file_a,file_b")
    ap.add_argument("--indir", type=str, help="Directory to auto-pair CSVs by filename base with trailing _A/_B")
    ap.add_argument("--delta", type=float, default=0.5, help="tolerance window in seconds")
    ap.add_argument("--perms", type=int, default=1000, help="number of circular-shift permutations")
    ap.add_argument("--t", type=float, default=None, help="optional analysis duration in seconds (e.g., 180)")
    ap.add_argument("--outdir", type=str, default=".", help="output directory")
    ap.add_argument("--seed", type=int, default=0, help="random seed for permutation")
    ap.add_argument("--win", type=float, default=10.0, help="sliding window size in seconds")
    ap.add_argument("--step", type=float, default=1.0, help="sliding window step in seconds")
    ap.add_argument("--episode_z", type=float, default=1.64, help="episode z-threshold (≈ one-sided 95%)")
    ap.add_argument("--episode_minmatch", type=int, default=3, help="minimum matched pairs per window for episode")
    ap.add_argument("--episode_minsec", type=float, default=3.0, help="minimum episode duration in seconds")
    ap.add_argument("--tag", type=str, default=None,
                    help="condition tag to attach into outputs (e.g., vis / invis). "
                         "If omitted, will use the last folder name of --indir.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    def short_id_from_base(base: str) -> str:
        """Example: 20170706_01_02_01 → concatenate the last four 2-digit chunks to make a short ID."""
        digs = regex.findall(r'(\d{2})', base)
        if len(digs) >= 4:
            return ''.join(digs[-4:])
        return base

    if args.indir:
        pairs = auto_discover_pairs(args.indir)
        if not pairs:
            print(f"[WARN] No A/B pairs found in {args.indir}. Filenames must end with _A.csv/_B.csv or -A/-B.")
            return
        print(f"[INFO] Found {len(pairs)} pairs:")
        for pid, fa, fb in pairs:
            print(f"  - {pid}: A={os.path.basename(fa)}  B={os.path.basename(fb)}")
        rows = []
        for pid, fa, fb in pairs:
            out_prefix = os.path.join(args.outdir, f"pair_{pid}")
            res = analyze_pair(
                fa, fb,
                delta=args.delta, T_total=args.t, perms=args.perms, seed=args.seed,
                win=args.win, step=args.step,
                episode_z=args.episode_z, episode_minmatch=args.episode_minmatch,
                episode_minsec=args.episode_minsec,
                out_prefix=out_prefix
            )
            res["pair_id"] = pid
            res["file_a"] = fa
            res["file_b"] = fb
            rows.append(res)
        write_rows(rows, args.outdir, args.tag, args.indir)
        return

    if args.pairs:
        df = pd.read_csv(args.pairs)
        rows = []
        for _, r in df.iterrows():
            pid = str(r["pair_id"])
            fa = str(r["file_a"])
            fb = str(r["file_b"])
            delta = float(r["delta"]) if "delta" in df.columns and not pd.isna(r["delta"]) else args.delta
            T_total = float(r["T"]) if "T" in df.columns and not pd.isna(r["T"]) else args.t
            out_prefix = os.path.join(args.outdir, f"pair_{pid}")
            res = analyze_pair(
                fa, fb,
                delta=delta, T_total=T_total, perms=args.perms, seed=args.seed,
                win=args.win, step=args.step,
                episode_z=args.episode_z, episode_minmatch=args.episode_minmatch,
                episode_minsec=args.episode_minsec,
                out_prefix=out_prefix
            )
            res["pair_id"] = pid
            res["file_a"] = fa
            res["file_b"] = fb
            rows.append(res)
        write_rows(rows, args.outdir, args.tag, None)
        return

    if args.a and args.b:
        out_prefix = os.path.join(args.outdir, "pair_single")
        res = analyze_pair(
            args.a, args.b,
            delta=args.delta, T_total=args.t, perms=args.perms, seed=args.seed,
            win=args.win, step=args.step,
            episode_z=args.episode_z, episode_minmatch=args.episode_minmatch,
            episode_minsec=args.episode_minsec,
            out_prefix=out_prefix
        )
        print(res)
    else:
        raise SystemExit("Provide --a and --b, or --pairs, or --indir.")

if __name__ == "__main__":
    main()
