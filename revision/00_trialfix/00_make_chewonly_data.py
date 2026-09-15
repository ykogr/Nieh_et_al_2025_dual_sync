#!/usr/bin/env python3
"""
00_make_chewonly_data.py  (Trial-column fix, step 0: data-level fix)

Purpose
  From each CSV in annotation_raw_data/ (columns: Trial, Chew_1, Chew_2, ...),
  remove the Trial column (trial number) and write a copy that keeps only the
  Chew_* columns to annotation_raw_data_chewonly/ with the same folder structure.

  The old pipeline's load_times_from_csv() numericised "every CSV cell" and
  treated the values as event times, so the Trial integers 1, 2, ..., k leaked
  in as "chewing events at t = 1 s, 2 s, ..., k s".
  With data in which the Trial column is physically absent, this contamination
  cannot occur regardless of the code version (a second safety net on top of the code fix 01).

  At the same time, the number of "spurious events that had leaked in" per file
  is written to a verification report (chewonly_report.csv / chewonly_report.txt).

Usage (run inside iscience_revision_analysis/00_trialfix/)
  python 00_make_chewonly_data.py \
      --src ../annotation_raw_data --dst ../annotation_raw_data_chewonly

  If dst already exists, add --overwrite to rebuild it.
"""
import argparse
import os
import shutil
import sys

import numpy as np
import pandas as pd


def read_csv_smart(path):
    for enc in ("utf-8", "utf-8-sig", "cp932", "latin1"):
        try:
            return pd.read_csv(path, header=0, encoding=enc)
        except UnicodeDecodeError:
            continue
    return pd.read_csv(path, header=0)


def legacy_times(df):
    """Same as the old load_times_from_csv: numericise all cells -> non-negative -> unique."""
    s = pd.to_numeric(pd.Series(df.to_numpy().ravel()), errors="coerce")
    arr = s.dropna().astype(float).values
    arr = arr[arr >= 0]
    return np.unique(np.sort(arr))


def is_chew_col(c):
    n = str(c).strip().lower()
    return n.startswith("chew")   # legacy non-English headers are normalised beforehand by 00_fix_legacy_header.py


def overflow_cols(df):
    """Columns without a header (Unnamed: k) that contain numbers = chew times continuing a trial that ran out of Chew_k columns.
    Do not discard them; keep them as Chew_<n+1>... (the old loader read every cell, so they were part of the old results too)."""
    out = []
    for c in df.columns:
        if str(c).startswith("Unnamed") and pd.to_numeric(df[c], errors="coerce").notna().any():
            out.append(c)
    return out


def chew_times(df):
    """After the fix: Chew* columns (+ value-bearing Unnamed columns) only -> numericise -> non-negative -> unique."""
    cols = [c for c in df.columns if is_chew_col(c)] + overflow_cols(df)
    if not cols:
        raise ValueError("No Chew columns found")
    s = pd.to_numeric(pd.Series(df[cols].to_numpy().ravel()), errors="coerce")
    arr = s.dropna().astype(float).values
    arr = arr[arr >= 0]
    return np.unique(np.sort(arr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="../annotation_raw_data")
    ap.add_argument("--dst", default="../annotation_raw_data_chewonly")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    src = os.path.abspath(args.src)
    dst = os.path.abspath(args.dst)
    if not os.path.isdir(src):
        sys.exit(f"[ERROR] src not found: {src}")
    if os.path.exists(dst):
        if not args.overwrite:
            sys.exit(f"[ERROR] dst already exists: {dst}  (use --overwrite to rebuild)")
        shutil.rmtree(dst)

    rows = []
    ovf_rows = []
    n_files = 0
    for root, dirs, files in os.walk(src):
        dirs.sort()
        rel = os.path.relpath(root, src)
        out_root = os.path.join(dst, rel) if rel != "." else dst
        os.makedirs(out_root, exist_ok=True)
        for fn in sorted(files):
            if not fn.lower().endswith(".csv"):
                continue
            p = os.path.join(root, fn)
            df = read_csv_smart(p)
            chew_cols = [c for c in df.columns if is_chew_col(c)]
            ovf = overflow_cols(df)
            other_cols = [c for c in df.columns if c not in chew_cols and c not in ovf]
            if not chew_cols:
                print(f"[WARN] no Chew columns, copied as is: {p}")
                shutil.copy2(p, os.path.join(out_root, fn))
                continue
            out = df[chew_cols + ovf].copy()
            if ovf:
                # Rename Unnamed: k -> Chew_<n+1>, Chew_<n+2>, ... and keep them
                ren = {c: f"Chew_{len(chew_cols) + i + 1}" for i, c in enumerate(ovf)}
                out = out.rename(columns=ren)
                n_ovf_cells = int(pd.to_numeric(df[ovf].stack(), errors="coerce").notna().sum())
                print(f"[INFO] {rel}/{fn}: {n_ovf_cells} chew times in {len(ovf)} header-less column(s) -> kept as {list(ren.values())[0]}..")
                ovf_rows.append({"rel_path": os.path.join(rel, fn), "n_overflow_cols": len(ovf), "n_overflow_cells": n_ovf_cells,
                                 "values": ";".join(f"{v:g}" for v in pd.to_numeric(df[ovf].stack(), errors="coerce").dropna())})
            out.to_csv(os.path.join(out_root, fn), index=False)
            n_files += 1

            leg = legacy_times(df)
            chw = chew_times(df)
            leaked = np.setdiff1d(leg, chw)          # spurious events that had leaked in
            # Chew times coinciding with a Trial number (absorbed by unique, so they vanished from the count)
            trial_vals = pd.to_numeric(df[other_cols].to_numpy().ravel(), errors="coerce") if other_cols else np.array([])
            trial_vals = np.unique(trial_vals[~np.isnan(trial_vals)]) if len(trial_vals) else trial_vals
            absorbed = np.intersect1d(trial_vals, chw)
            rows.append({
                "rel_path": os.path.relpath(p, src),
                "n_trial_rows": len(df),
                "dropped_columns": ";".join(map(str, other_cols)),
                "n_events_legacy": len(leg),
                "n_events_chew_only": len(chw),
                "n_spurious_events_removed": len(leaked),
                "spurious_values": ";".join(f"{v:g}" for v in leaked),
                "n_trial_values_coinciding_with_real_chews": len(absorbed),
                "first_event_legacy": leg[0] if len(leg) else np.nan,
                "first_event_chew_only": chw[0] if len(chw) else np.nan,
            })

    rep = pd.DataFrame(rows)
    rep_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)), "chewonly_report.csv")
    rep.to_csv(rep_csv, index=False)

    lines = []
    lines.append(f"src: {src}")
    lines.append(f"dst: {dst}")
    lines.append(f"CSV files processed: {n_files}")
    allcols = sorted({c for v in rep["dropped_columns"] for c in str(v).split(";") if c})
    named = [c for c in allcols if not c.startswith("Unnamed")]
    n_unnamed = sum(c.startswith("Unnamed") for c in allcols)
    lines.append(f"dropped named columns (unique): {named}"
                 + (f"   + {n_unnamed} kinds of empty 'Unnamed: k' columns (empty trailing header columns with no values)" if n_unnamed else ""))
    lines.append(f"files with empty trailing columns dropped: {int(rep['dropped_columns'].astype(str).str.contains('Unnamed').sum())}")
    lines.append(f"files with value-bearing unnamed columns kept as extra Chew_k columns: {len(ovf_rows)} "
                 f"({sum(r['n_overflow_cells'] for r in ovf_rows)} chew times recovered)")
    for r in ovf_rows:
        lines.append(f"    {r['rel_path']}: {r['n_overflow_cells']} values ({r['values']})")
    lines.append(f"total events (legacy, all cells)   : {int(rep['n_events_legacy'].sum())}")
    lines.append(f"total events (Chew columns only)   : {int(rep['n_events_chew_only'].sum())}")
    lines.append(f"total spurious events removed      : {int(rep['n_spurious_events_removed'].sum())}")
    lines.append(f"files with >=1 spurious event      : {int((rep['n_spurious_events_removed'] > 0).sum())} / {len(rep)}")
    lines.append(f"spurious events per file: min {rep['n_spurious_events_removed'].min()}, "
                 f"median {rep['n_spurious_events_removed'].median():.1f}, max {rep['n_spurious_events_removed'].max()}")
    lines.append(f"first event (legacy) : always {rep['first_event_legacy'].min():g}-{rep['first_event_legacy'].max():g} s "
                 f"(= Trial number 1)")
    lines.append(f"first event (chew)   : {rep['first_event_chew_only'].min():g}-{rep['first_event_chew_only'].max():g} s")
    txt = "\n".join(lines)
    with open(rep_csv.replace(".csv", ".txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)
    print(f"[DONE] report -> {rep_csv}")


if __name__ == "__main__":
    main()
