#!/usr/bin/env python3
"""
03_R2-11_dip_test.py
R2-11: Hartigan's dip test on raw matched-event latency differences (dt),
with dequantization (jitter) sensitivity analysis.

Background
----------
Reviewer R2-11 asked that the dip test be repeated on the raw latency
differences rather than histogram counts, noting that the reported results
(p ~ 1.0 with Dip = 0.09-0.12) are unusual and that 0.1-s binning may be
too coarse. The raw dt values (dt.csv, one column per condition x delta,
pooled over all dyads and sessions) are quantized at the 0.1-s annotation
resolution, so each column contains only 3-7 distinct values. Applied
directly to such data, the dip test responds to the atoms (ties) created
by quantization rather than to the shape of the underlying latency
distribution.

This script therefore reports, for each column of dt.csv:
  (1) the dip test on the raw quantized values (as in the original
      analysis code), and
  (2) a dequantization sensitivity analysis: uniform jitter
      U(-r/2, +r/2) with r = 0.1 s (the annotation resolution) is added
      to each value, the dip test is recomputed, and this is repeated
      --n-jitter times; the distribution of dip statistics and p-values
      across replicates is summarized.

Requirements
------------
  conda activate nieh_iscience_2026
  pip install diptest        # if not yet installed (same package used in
                             # the original analysis, diptest >= 0.11)

Usage
-----
  python 03_R2-11_dip_test.py --dtcsv dt.csv --outdir dip_test_out \
      [--n-jitter 1000] [--seed 20260815]

Outputs (written to --outdir)
-----------------------------
  dip_test_raw_quantized.csv   : dip and table p on the raw quantized dt
                                 (reproduction of the original analysis)
  dip_test_jitter_summary.csv  : summary of dip/p over jitter replicates
  dip_test_summary.txt         : human-readable summary of both analyses
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

try:
    from diptest import diptest
except ImportError:
    sys.exit("The 'diptest' package is required. Install it with:\n"
             "  pip install diptest")

RESOLUTION = 0.1  # annotation grid (s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dtcsv", required=True,
                    help="dt.csv: one column per condition x delta, "
                         "pooled raw latency differences")
    ap.add_argument("--outdir", default="dip_test_out")
    ap.add_argument("--n-jitter", type=int, default=1000,
                    help="number of jitter replicates (default 1000)")
    ap.add_argument("--seed", type=int, default=20260815)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.dtcsv)
    cols = [c for c in df.columns
            if isinstance(c, str) and ("vis" in c or "invis" in c)]
    if not cols:
        sys.exit("No 'vis'/'invis' columns found in " + args.dtcsv)

    rng = np.random.default_rng(args.seed)
    half = RESOLUTION / 2.0

    raw_rows, jit_rows, lines = [], [], []
    lines.append("R2-11 Hartigan's dip test on raw latency differences (dt)")
    lines.append(f"input: {args.dtcsv}; jitter replicates: {args.n_jitter}; "
                 f"jitter: U(-{half}, +{half}) s; seed: {args.seed}")
    lines.append("")

    for col in cols:
        x = pd.to_numeric(df[col], errors="coerce").dropna().values
        n = len(x)
        uniq, cnts = np.unique(x, return_counts=True)

        # (1) raw quantized values (original analysis)
        dip_raw, p_raw = diptest(x)
        raw_rows.append({"column": col, "n": n, "n_unique": len(uniq),
                         "dip": dip_raw, "p_value": p_raw})

        # (2) jitter dequantization
        dips = np.empty(args.n_jitter)
        ps = np.empty(args.n_jitter)
        for i in range(args.n_jitter):
            xj = x + rng.uniform(-half, half, size=n)
            dips[i], ps[i] = diptest(xj)
        jit_rows.append({
            "column": col, "n": n, "n_unique_raw": len(uniq),
            "dip_median": np.median(dips),
            "dip_p2.5": np.percentile(dips, 2.5),
            "dip_p97.5": np.percentile(dips, 97.5),
            "p_median": np.median(ps),
            "p_p2.5": np.percentile(ps, 2.5),
            "p_p97.5": np.percentile(ps, 97.5),
            "prop_p_below_0.05": float(np.mean(ps < 0.05)),
        })

        lines.append(f"--- {col} (n = {n}) ---")
        lines.append("  raw value counts: " +
                     ", ".join(f"{v:+.1f}: {c}" for v, c in zip(uniq, cnts)))
        lines.append(f"  raw quantized:  dip = {dip_raw:.6f}, p = {p_raw:.6g} "
                     f"({len(uniq)} distinct values)")
        lines.append(f"  jittered:       dip median = {np.median(dips):.6f} "
                     f"[{np.percentile(dips, 2.5):.6f}, "
                     f"{np.percentile(dips, 97.5):.6f}],")
        lines.append(f"                  p   median = {np.median(ps):.4f} "
                     f"[{np.percentile(ps, 2.5):.4f}, "
                     f"{np.percentile(ps, 97.5):.4f}], "
                     f"prop(p < 0.05) = {np.mean(ps < 0.05):.4f}")
        lines.append("")

    pd.DataFrame(raw_rows).to_csv(
        os.path.join(args.outdir, "dip_test_raw_quantized.csv"), index=False)
    pd.DataFrame(jit_rows).to_csv(
        os.path.join(args.outdir, "dip_test_jitter_summary.csv"), index=False)
    txt = "\n".join(lines)
    with open(os.path.join(args.outdir, "dip_test_summary.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
