#!/usr/bin/env python3
"""
05_prevalence_table2.py

Dyad-level prevalence of significant chewing synchrony (Table 2).

Reads the per-session summary CSVs produced by 03_Micro_Analysis_NOLPM_gridnull.py
({delta}_chew_sync_summary_{visible|invisible}_pair_Human.csv) from --indir and,
for each condition x delta family, reports the number of dyads (out of 10) whose
circular-shift Monte Carlo p-values (p_S) are significant (p < .05) in
>=1, >=2, or all 3 sessions, under:
  - Holm correction across the 3 sessions within each dyad
    (uncorrected counts in parentheses)
  - Benjamini-Hochberg FDR correction across all 30 session-level
    tests within each condition x delta family
Additionally, for each condition at delta = 0.1, an exact one-sided binomial
test compares the uncorrected >=1-session count against the chance expectation
1 - 0.95^3.

Fully deterministic (no RNG). Output: prints the table and writes
prevalence_table2.csv into --indir.

Usage:
  python 05_prevalence_table2.py --indir Micro_Analysis_LMM
"""
import argparse
import glob
import math
import os
import re

import numpy as np
import pandas as pd

DELTAS = ["0.1", "0.2", "0.3", "0.5", "1.0"]
CONDITIONS = ["visible", "invisible"]
ALPHA = 0.05


def holm_adjust(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    running_max = 0.0
    for rank, idx in enumerate(order):
        val = (n - rank) * p[idx]
        running_max = max(running_max, val)
        adj[idx] = min(1.0, running_max)
    return adj


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n)
    running_min = 1.0
    for rank in range(n - 1, -1, -1):
        idx = order[rank]
        val = p[idx] * n / (rank + 1)
        running_min = min(running_min, val)
        adj[idx] = running_min
    return adj


def binom_sf_geq(k, n, p0):
    """Exact P(X >= k), X ~ Binomial(n, p0)."""
    return sum(math.comb(n, i) * p0**i * (1 - p0) ** (n - i) for i in range(k, n + 1))


def counts_by_dyad(df, col):
    sig = df[col] < ALPHA
    k = sig.groupby(df["dyad"]).sum()
    return [int((k >= m).sum()) for m in (1, 2, 3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", type=str, default="Micro_Analysis_LMM",
                    help="folder containing {delta}_chew_sync_summary_*.csv")
    args = ap.parse_args()

    rows = []
    binom_lines = []
    for d in DELTAS:
        for cond in CONDITIONS:
            path = os.path.join(args.indir, f"{d}_chew_sync_summary_{cond}_pair_Human.csv")
            if not os.path.exists(path):
                raise SystemExit(f"[ERROR] not found: {path}")
            df = pd.read_csv(path)
            if len(df) != 30:
                raise SystemExit(f"[ERROR] {path}: expected 30 rows, got {len(df)}")
            df["dyad"] = df["pair_id"].str.replace(r"_\d{2}$", "", regex=True)
            if df["dyad"].nunique() != 10:
                raise SystemExit(f"[ERROR] {path}: expected 10 dyads, got {df['dyad'].nunique()}")
            # Holm within each dyad (3 sessions)
            df["p_holm"] = np.nan
            for dy, idx in df.groupby("dyad").groups.items():
                df.loc[idx, "p_holm"] = holm_adjust(df.loc[idx, "p_S"].values)
            # BH across the 30 session-level tests in this condition x delta family
            df["p_bh"] = bh_adjust(df["p_S"].values)

            c_unc = counts_by_dyad(df, "p_S")
            c_holm = counts_by_dyad(df, "p_holm")
            c_bh = counts_by_dyad(df, "p_bh")
            rows.append({
                "delta": d, "condition": cond.capitalize(),
                "holm_ge1": f"{c_holm[0]} ({c_unc[0]})",
                "holm_ge2": f"{c_holm[1]} ({c_unc[1]})",
                "holm_all3": f"{c_holm[2]} ({c_unc[2]})",
                "bh_ge1": c_bh[0], "bh_ge2": c_bh[1], "bh_all3": c_bh[2],
            })
            if d == "0.1":
                p0 = 1.0 - 0.95**3
                pval = binom_sf_geq(c_unc[0], 10, p0)
                binom_lines.append(
                    f"delta=0.1 {cond}: uncorrected >=1-session count = {c_unc[0]}/10, "
                    f"chance p0 = {p0:.6f}, exact binomial (one-sided) p = {pval:.3g}")

    out = pd.DataFrame(rows)
    print("=" * 78)
    print("Table 2: dyads (out of 10) with significant synchrony (p_S < .05)")
    print("Holm: within-dyad across 3 sessions (uncorrected in parentheses)")
    print("BH-FDR: across all 30 session-level tests per condition x delta family")
    print("=" * 78)
    print(out.to_string(index=False))
    print()
    for line in binom_lines:
        print(line)
    dst = os.path.join(args.indir, "prevalence_table2.csv")
    out.to_csv(dst, index=False)
    print(f"\n[SAVED] {dst}")


if __name__ == "__main__":
    main()
