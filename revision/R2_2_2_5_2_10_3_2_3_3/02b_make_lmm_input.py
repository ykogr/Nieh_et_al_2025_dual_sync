#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2b of the R2-10/R2-5/R2-2 branch (R2-5): Convert the jitter-surrogate master
(from 02_jitter_surrogates.py) into the input format expected by
04_LMM_condition_by_period.R, so that the condition x period LMM and the
emmeans contrasts can be re-estimated on z_S values computed against a chosen
rate-preserving null.

It filters jitter_null_master.csv to ONE null configuration (e.g. null_type
jitter, w = 0.35) and writes <outdir>/chew_sync_summary_periods_master.csv
with the columns 04 expects (STR / Max_episode added as empty; episode metrics
are not defined for the surrogate reruns). All delta values present in the
master are kept (04 iterates over them).

Usage:
  python 02b_make_lmm_input.py --master out_jitter/jitter_null_master.csv \
      --null_type jitter --w 0.35 --outdir out_jitter/lmm_in_jitter035
  Rscript 04_LMM_condition_by_period.R out_jitter/lmm_in_jitter035

  # optional sensitivity configs:
  python 02b_make_lmm_input.py --master out_jitter/jitter_null_master.csv \
      --null_type jitter --w 0.2 --outdir out_jitter/lmm_in_jitter020
  python 02b_make_lmm_input.py --master out_jitter/jitter_null_master.csv \
      --null_type jitter --w 0.5 --outdir out_jitter/lmm_in_jitter050
  python 02b_make_lmm_input.py --master out_jitter/jitter_null_master.csv \
      --null_type isi --outdir out_jitter/lmm_in_isi
"""
import argparse
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=str, required=True,
                    help="path to jitter_null_master.csv")
    ap.add_argument("--null_type", type=str, required=True,
                    choices=["jitter", "isi", "circ"])
    ap.add_argument("--w", type=float, default=None,
                    help="jitter half-width to select (required for null_type=jitter)")
    ap.add_argument("--outdir", type=str, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.master)
    sel = df["null_type"] == args.null_type
    if args.null_type == "jitter":
        if args.w is None:
            raise SystemExit("[ERROR] --w is required for null_type=jitter")
        sel &= np.isclose(df["jitter_w"], args.w)
    out = df[sel].copy()
    if out.empty:
        raise SystemExit("[ERROR] no rows match the requested null configuration")

    # columns expected by 04_LMM_condition_by_period.R
    for col in ("STR", "Max_episode"):
        out[col] = np.nan
    cols = ["pair_id", "condition", "delta", "period", "t0_rel", "t1_rel", "T_seg",
            "n_A", "n_B", "sufficient", "n_matched", "S_obs",
            "null_mean", "null_sd", "z_S", "p_right", "p_left", "STR", "Max_episode"]
    out = out[[c for c in cols if c in out.columns]]

    os.makedirs(args.outdir, exist_ok=True)
    path = os.path.join(args.outdir, "chew_sync_summary_periods_master.csv")
    out.to_csv(path, index=False)
    print(f"[DONE] {path}")
    print(f"  null_type={args.null_type}"
          + (f", w={args.w}" if args.null_type == "jitter" else ""))
    print(f"  rows={len(out)}, deltas={sorted(out['delta'].unique().tolist())}, "
          f"conditions={sorted(out['condition'].unique().tolist())}, "
          f"periods={sorted(out['period'].unique().tolist())}, "
          f"sufficient=1: {int((out['sufficient'] == 1).sum())}")


if __name__ == "__main__":
    main()
