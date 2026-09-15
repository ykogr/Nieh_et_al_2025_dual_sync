#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R2-9: Multiple-comparison correction for dyad x session Monte Carlo p-values
(circular-shift test, p_S) underlying Table 2.

Reads the per-session summary CSVs used for Table 1/2
({delta}_chew_sync_summary_{visible|invisible}_pair_Human.csv, as in the
repository folder Micro_Analysis_LMM/) and recomputes the Table 2 prevalence
counts (>=1 / >=2 / >=3 significant sessions per dyad) under:

  (1) uncorrected p < .05            (as currently reported in Table 2)
  (2) Holm correction WITHIN each dyad across its 3 sessions
      (family = the 3 sessions supporting a dyad-level claim; directly
       addresses the reviewer's concern about the liberal >=1-session criterion)
  (3) Benjamini-Hochberg FDR across all 30 sessions WITHIN each
      condition x delta family (global descriptive correction)

Also reports, per condition x delta, an exact binomial test of the number of
dyads with >=1 uncorrected-significant session against the chance expectation
P(>=1 of 3 | H0) = 1 - 0.95^3 = 0.142625 (descriptive context).

Outputs (written to --outdir):
  session_pvalues_adjusted.csv  : all session-level p_S with p_holm_within_dyad
                                  and p_bh_within_family columns
  table2_corrected_counts.csv   : Table 2 counts under (1)-(3)
  multiplicity_summary.txt      : human-readable summary

Usage:
  python 02_R2-9_multiplicity_correction.py --indir <path/to/Micro_Analysis_LMM> --outdir multiplicity_out
"""
import argparse
import glob
import math
import os
import re

import numpy as np
import pandas as pd

ALPHA = 0.05


def holm_adjust(pvals):
    """Holm step-down adjusted p-values (monotone, capped at 1)."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m, dtype=float)
    running_max = 0.0
    for rank, idx in enumerate(order):
        val = (m - rank) * p[idx]
        running_max = max(running_max, val)
        adj[idx] = min(1.0, running_max)
    return adj


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values (monotone, capped at 1)."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m, dtype=float)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        val = min(prev, p[idx] * m / (rank + 1))
        adj[idx] = val
        prev = val
    return adj


def binom_sf_geq(k, n, prob):
    """Exact P(X >= k) for X ~ Binomial(n, prob)."""
    return sum(math.comb(n, i) * prob**i * (1 - prob) ** (n - i) for i in range(k, n + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True,
                    help="directory containing {delta}_chew_sync_summary_{cond}_pair_Human.csv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--alpha", type=float, default=ALPHA)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    alpha = args.alpha

    pat = re.compile(r"([\d.]+)_chew_sync_summary_(visible|invisible)_pair_Human\.csv$")
    rows = []
    for f in sorted(glob.glob(os.path.join(args.indir, "*_chew_sync_summary_*_pair_Human.csv"))):
        m = pat.search(os.path.basename(f))
        if not m:
            continue
        d = pd.read_csv(f)
        d["delta"] = float(m.group(1))
        d["condition"] = m.group(2)
        rows.append(d)
    if not rows:
        raise SystemExit(f"[ERROR] no summary CSVs found in {args.indir}")
    df = pd.concat(rows, ignore_index=True)

    df["session"] = df["pair_id"].str.extract(r"_(\d{2})$")[0].astype(int)
    df["dyad_id"] = df["pair_id"].str.replace(r"_\d{2}$", "", regex=True)

    # --- adjusted p-values ---
    df["p_holm_within_dyad"] = np.nan
    for _, idx in df.groupby(["delta", "condition", "dyad_id"]).groups.items():
        df.loc[idx, "p_holm_within_dyad"] = holm_adjust(df.loc[idx, "p_S"].values)
    df["p_bh_within_family"] = np.nan
    for _, idx in df.groupby(["delta", "condition"]).groups.items():
        df.loc[idx, "p_bh_within_family"] = bh_adjust(df.loc[idx, "p_S"].values)

    keep = ["delta", "condition", "dyad_id", "session", "pair_id",
            "z_S", "p_S", "p_holm_within_dyad", "p_bh_within_family"]
    df[keep].sort_values(["delta", "condition", "dyad_id", "session"]).to_csv(
        os.path.join(args.outdir, "session_pvalues_adjusted.csv"), index=False)

    # --- Table 2 counts under each scheme ---
    def dyad_counts(pcol):
        nsig = (df[pcol] < alpha).groupby(
            [df["delta"], df["condition"], df["dyad_id"]]).sum()
        out = nsig.reset_index(name="n_sig").groupby(["delta", "condition"])["n_sig"].agg(
            ge1=lambda s: int((s >= 1).sum()),
            ge2=lambda s: int((s >= 2).sum()),
            ge3=lambda s: int((s >= 3).sum()))
        return out

    raw = dyad_counts("p_S").rename(columns=lambda c: f"raw_{c}")
    holm = dyad_counts("p_holm_within_dyad").rename(columns=lambda c: f"holm_{c}")
    bh = dyad_counts("p_bh_within_family").rename(columns=lambda c: f"bh_{c}")
    tab = raw.join(holm).join(bh).reset_index()

    # binomial context: dyads with >=1 uncorrected-significant session vs chance
    p_chance = 1.0 - (1.0 - alpha) ** 3
    n_dyads = df.groupby(["delta", "condition"])["dyad_id"].nunique().reset_index(name="n_dyads")
    tab = tab.merge(n_dyads, on=["delta", "condition"])
    tab["binom_p_ge1_vs_chance"] = [
        binom_sf_geq(int(k), int(n), p_chance)
        for k, n in zip(tab["raw_ge1"], tab["n_dyads"])]
    tab.to_csv(os.path.join(args.outdir, "table2_corrected_counts.csv"), index=False)

    # --- text summary ---
    lines = []
    lines.append("R2-9 multiple-comparison correction of Table 2 session-level Monte Carlo p-values")
    lines.append(f"alpha = {alpha}; Holm family = 3 sessions within dyad; "
                 f"BH family = all sessions within condition x delta")
    lines.append(f"chance P(>=1 of 3 sessions significant | H0) = {p_chance:.6f}")
    lines.append("")
    lines.append(tab.to_string(index=False))
    txt = "\n".join(lines)
    with open(os.path.join(args.outdir, "multiplicity_summary.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
