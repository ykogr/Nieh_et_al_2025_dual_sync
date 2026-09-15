#!/usr/bin/env python3
"""
02_verify_fix.py  (Trial-column fix, step 2: verification)

Runs three checks together and writes them to verify_report.txt.

 (1) Code check: among the .py files in the given folders, do all files that contain
     load_times_from_csv / NLOPM_match / nlopm_matches carry the patched markers?
 (2) Data check: does annotation_raw_data_chewonly/ still contain any Trial column,
     and does the file count match the original?
 (3) Numeric check (per pair video; no permutation test, so a few seconds):
       - difference in event counts between the old loader (all cells) and the new loader (Chew columns only)
       - how much the NLOPM match count changes between the old comparison (no tolerance) and the new one (FLOAT_TOL)
     If a patched script (e.g. ../gridnull_main/03_Micro_Analysis_NOLPM_gridnull_v2.py) is passed
     via --script, its load_times_from_csv / NLOPM_match are actually imported and
     checked to give the same results as this verification implementation.

Usage (run inside iscience_revision_analysis/00_trialfix/)
  python 02_verify_fix.py --raw ../annotation_raw_data --chewonly ../annotation_raw_data_chewonly \
      --folders ../gridnull_main ../R1_2 ../R2_2_2_5_2_10_3_2_3_3 ../R2_4_3_1 \
      --script ../gridnull_main/03_Micro_Analysis_NOLPM_gridnull_v2.py
"""
import argparse
import glob
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

FLOAT_TOL = 1e-6


def read_csv_smart(path):
    for enc in ("utf-8", "utf-8-sig", "cp932", "latin1"):
        try:
            return pd.read_csv(path, header=0, encoding=enc)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, header=0)


def load_legacy(path):
    df = read_csv_smart(path)
    s = pd.to_numeric(pd.Series(df.to_numpy().ravel()), errors="coerce")
    arr = s.dropna().astype(float).values
    return np.unique(np.sort(arr[arr >= 0]))


def load_fixed(path):
    df = read_csv_smart(path)
    # Value-bearing Unnamed columns (continuation of rows that ran out of Chew_k) are renamed to Chew_k and kept by 00, so include them here too
    cols = [c for c in df.columns if str(c).strip().lower().startswith("chew")
            or (str(c).startswith("Unnamed") and pd.to_numeric(df[c], errors="coerce").notna().any())]
    s = pd.to_numeric(pd.Series(df[cols].to_numpy().ravel()), errors="coerce")
    arr = s.dropna().astype(float).values
    return np.unique(np.sort(arr[arr >= 0]))


def nlopm(a, b, delta, tol):
    """Original NLOPM_match with an optional tolerance (tol=0 -> legacy)."""
    i = j = 0
    m, n = len(a), len(b)
    pairs = []
    while i < m and j < n:
        while j < n and b[j] < a[i] - delta - tol:
            j += 1
        if j >= n:
            break
        while i < m and a[i] < b[j] - delta - tol:
            i += 1
        if i >= m:
            break
        js = j
        while (js + 1) < n and (b[js + 1] <= a[i] + delta + tol) \
                and (abs(a[i] - b[js + 1]) <= abs(a[i] - b[js]) + tol):
            js += 1
        if abs(a[i] - b[js]) <= delta + tol:
            pairs.append((a[i], b[js]))
            i += 1
            j = js + 1
        else:
            if a[i] < b[j]:
                i += 1
            else:
                j += 1
    return pairs


def restrict(a, b):
    start, end = max(a[0], b[0]), min(a[-1], b[-1])
    return a[(a >= start) & (a <= end)], b[(b >= start) & (b <= end)], start, end


def pair_files(cond_dir):
    out = []
    for pa in sorted(glob.glob(os.path.join(cond_dir, "A", "*.csv"))):
        base = os.path.basename(pa)[:-6]  # strip _A.csv
        pb = os.path.join(cond_dir, "B", base + "_B.csv")
        if os.path.exists(pb):
            out.append((base, pa, pb))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="../annotation_raw_data")
    ap.add_argument("--chewonly", default="../annotation_raw_data_chewonly")
    ap.add_argument("--folders", nargs="*", default=[])
    ap.add_argument("--script", default=None)
    ap.add_argument("--delta", type=float, default=0.1)
    args = ap.parse_args()

    L = []
    ok_all = True

    # (1) code markers -----------------------------------------------------
    L.append("(1) code check")
    for folder in args.folders:
        for fn in sorted(os.listdir(folder)):
            if not fn.endswith(".py"):
                continue
            src = open(os.path.join(folder, fn), encoding="utf-8", errors="replace").read()
            need_loader = "def load_times_from_csv(" in src
            need_match = ("def NLOPM_match(" in src) or ("def nlopm_matches(" in src)
            need_time = "def sliding_window_with_dt(" in src
            need_period = "seg_start = start + t0" in src
            if not (need_loader or need_match or need_period):
                continue
            ok_loader = (not need_loader) or ("[EVENT-COLUMN PATCH" in src)
            ok_match = (not need_match) or ("[FLOAT-TOL PATCH" in src)
            ok_time = (not need_time) or ("[TIME-AXIS PATCH" in src)
            ok_period = (not need_period) or ("[PERIOD-ORIGIN PATCH" in src)
            ok = ok_loader and ok_match and ok_time and ok_period
            st = "OK " if ok else "NG "
            ok_all &= ok
            parts = [f"loader:{'patched' if ok_loader else 'UNPATCHED'}" if need_loader else None,
                     f"match:{'patched' if ok_match else 'UNPATCHED'}" if need_match else None,
                     f"time_axis:{'patched' if ok_time else 'UNPATCHED'}" if need_time else None,
                     f"period_origin:{'patched' if ok_period else 'UNPATCHED'}" if need_period else None]
            L.append(f"  {st} {os.path.join(folder, fn)}  " + "  ".join(x for x in parts if x))

    # (2) data ----------------------------------------------------------------
    L.append("(2) data check")
    raw_files = sorted(glob.glob(os.path.join(args.raw, "**", "*.csv"), recursive=True))
    co_files = sorted(glob.glob(os.path.join(args.chewonly, "**", "*.csv"), recursive=True))
    L.append(f"  raw csv files: {len(raw_files)}   chewonly csv files: {len(co_files)}")
    if len(raw_files) != len(co_files):
        L.append("  NG file counts differ")
        ok_all = False
    bad = []
    for p in co_files:
        cols = [str(c).strip().lower() for c in pd.read_csv(p, nrows=0).columns]
        if "trial" in cols or any(c.startswith("unnamed") for c in cols):
            bad.append(p)
    L.append(f"  chewonly files still containing Trial/Unnamed columns: {len(bad)}")
    # Do the Unnamed columns of the raw data hold values? (if so, 00 should have kept them as Chew_k)
    n_ovf = 0
    for f in raw_files:
        d = read_csv_smart(f)
        for c in d.columns:
            if str(c).startswith("Unnamed"):
                n_ovf += int(pd.to_numeric(d[c], errors="coerce").notna().sum())
    n_co = sum(int(pd.to_numeric(read_csv_smart(f).stack(), errors="coerce").notna().sum()) for f in co_files) if co_files else 0
    L.append(f"  chew times in unnamed (overflow) columns of raw data: {n_ovf}  -> kept in chewonly copy as extra Chew_k columns")
    L.append(f"  total numeric cells in chewonly copy: {n_co}")
    if bad:
        ok_all = False
        L.extend("    " + b for b in bad[:10])

    # (3) numbers -------------------------------------------------------------
    L.append(f"(3) numeric check (delta = {args.delta})")
    mod = None
    if args.script:
        spec = importlib.util.spec_from_file_location("patched_mod", args.script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        L.append(f"  imported {args.script}")
    rows = []
    for cond in ("visible_pair_Human", "invisible_pair_Human"):
        for base, pa, pb in pair_files(os.path.join(args.raw, cond)):
            a_leg, b_leg = load_legacy(pa), load_legacy(pb)
            a_fix, b_fix = load_fixed(pa), load_fixed(pb)
            a2, b2, s_leg, _ = restrict(a_leg, b_leg)
            m_leg = len(nlopm(a2, b2, args.delta, 0.0))
            a3, b3, s_fix, _ = restrict(a_fix, b_fix)
            m_fix_notol = len(nlopm(a3, b3, args.delta, 0.0))
            m_fix = len(nlopm(a3, b3, args.delta, FLOAT_TOL))
            n_spur = len(np.setdiff1d(a_leg, a_fix)) + len(np.setdiff1d(b_leg, b_fix))
            row = dict(condition=cond.split("_")[0], pair=base, n_spurious_removed=n_spur,
                       nA_legacy=len(a2), nA_fixed=len(a3), nB_legacy=len(b2), nB_fixed=len(b3),
                       start_legacy=s_leg, start_fixed=s_fix,
                       matched_legacy=m_leg, matched_chewonly_no_tol=m_fix_notol, matched_fixed=m_fix)
            if mod is not None:
                pa2 = os.path.join(args.chewonly, cond, "A", os.path.basename(pa))
                pb2 = os.path.join(args.chewonly, cond, "B", os.path.basename(pb))
                ma = mod.load_times_from_csv(pa2)
                mb = mod.load_times_from_csv(pb2)
                ma, mb, _, _ = restrict(ma, mb)
                M, _ = mod.NLOPM_match(ma, mb, args.delta)
                row["matched_patched_script"] = len(M)
                row["script_agrees"] = (len(M) == m_fix) and (len(ma) == len(a3)) and (len(mb) == len(b3))
            rows.append(row)
    df = pd.DataFrame(rows)
    here = os.path.dirname(os.path.abspath(__file__))
    df.to_csv(os.path.join(here, "verify_pairs.csv"), index=False)
    L.append(f"  pairs checked: {len(df)}")
    L.append(f"  spurious (Trial/Unnamed) events removed, pair files only (A+B, whole file): "
             f"{int(df.n_spurious_removed.sum())}   [the total in chewonly_report.txt also includes the solo files]")
    L.append(f"  events inside the overlap window, legacy -> fixed: "
             f"{int((df.nA_legacy + df.nB_legacy).sum())} -> {int((df.nA_fixed + df.nB_fixed).sum())} "
             f"(difference = spurious events + real chews lost because the overlap now starts later)")
    L.append(f"  overlap start: legacy always {df.start_legacy.min():g}-{df.start_legacy.max():g} s; "
             f"fixed {df.start_fixed.min():g}-{df.start_fixed.max():g} s (mean {df.start_fixed.mean():.1f})")
    L.append(f"  matched events, all pairs: legacy {int(df.matched_legacy.sum())} -> chew-only w/o tolerance "
             f"{int(df.matched_chewonly_no_tol.sum())} -> chew-only + FLOAT_TOL {int(df.matched_fixed.sum())}")
    L.append(f"  pairs whose match count changes with FLOAT_TOL alone: "
             f"{int((df.matched_fixed != df.matched_chewonly_no_tol).sum())} / {len(df)}")
    if mod is not None:
        agree = bool(df.script_agrees.all())
        ok_all &= agree
        L.append(f"  patched script reproduces this verification implementation on all pairs: {agree}")

    L.append("")
    L.append("RESULT: " + ("ALL OK" if ok_all else "PROBLEMS FOUND (see NG lines)"))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(here, "verify_report.txt"), "w") as f:
        f.write(txt + "\n")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
