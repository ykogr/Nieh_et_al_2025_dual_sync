#!/usr/bin/env python3
"""
make_dt_csv.py  -- rebuild the latency table (dt.csv) for R2-11 (dip test) and Fig. 2

For each delta and condition, load the chew-only annotation CSVs with the
PATCHED loader / matcher of 03_Micro_Analysis_NOLPM_gridnull_v2.py (same
folder), restrict to the overlap interval (and the --t cap, default 180 s),
run NLOPM once on the whole interval and collect dt = t_A - t_B of every
matched pair. Each matched pair is counted once (whole-interval matching).

Output: one wide CSV, columns "<delta>_vis" and "<delta>_invis"
(e.g. 0.1_vis, 0.1_invis, 0.2_vis, ...), padded with NaN to equal length.
This is the column convention expected by
R2_8_2_9_2_11/03_R2-11_dip_test.py and 04_R2-11_dip_figure.py.

Usage (in gridnull_main/):
  python make_dt_csv.py --raw ../annotation_raw_data_chewonly \
      --deltas 0.1,0.2,0.3 --out dt_trialfix.csv
"""
import argparse
import glob
import importlib.util
import os

import numpy as np
import pandas as pd


def load_main_module(path):
    spec = importlib.util.spec_from_file_location("nlopm_main", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert "[EVENT-COLUMN PATCH" in open(path, encoding="utf-8").read(), \
        f"{path} is not patched; run 00_trialfix/01_patch_scripts.py first"
    return mod


def pair_files(cond_dir):
    out = []
    for pa in sorted(glob.glob(os.path.join(cond_dir, "A", "*.csv"))):
        base = os.path.basename(pa)[:-6]
        pb = os.path.join(cond_dir, "B", base + "_B.csv")
        if os.path.exists(pb):
            out.append((base, pa, pb))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default="../annotation_raw_data_chewonly")
    ap.add_argument("--script", default="03_Micro_Analysis_NOLPM_gridnull_v2.py")
    ap.add_argument("--deltas", default="0.1,0.2,0.3")
    ap.add_argument("--t", type=float, default=180.0)
    ap.add_argument("--out", default="dt_trialfix.csv")
    args = ap.parse_args()

    mod = load_main_module(args.script)
    cols = {}
    rows = []
    for delta in [float(x) for x in args.deltas.split(",")]:
        for cond, short in (("visible_pair_Human", "vis"), ("invisible_pair_Human", "invis")):
            dts = []
            for base, pa, pb in pair_files(os.path.join(args.raw, cond)):
                a = mod.load_times_from_csv(pa)
                b = mod.load_times_from_csv(pb)
                a2, b2, start, end = mod.restrict_overlap(a, b)
                if args.t and args.t > 0:
                    end = min(end, start + args.t)
                    a2 = a2[(a2 >= start) & (a2 <= end)]
                    b2 = b2[(b2 >= start) & (b2 <= end)]
                _, dt = mod.NLOPM_match(a2, b2, delta)
                dt = np.round(np.asarray(dt, dtype=float), 3)   # remove float noise (0.1-s grid)
                dts.extend(dt.tolist())
                rows.append(dict(delta=delta, condition=short, pair_id=base, n_matched=len(dt),
                                 median_abs_dt=float(np.median(np.abs(dt))) if len(dt) else np.nan))
            cols[f"{delta:g}_{short}"] = dts
            print(f"delta={delta:g} {short:5s}: {len(dts)} matched events; "
                  f"share |dt|=0: {np.mean(np.abs(dts) < 1e-6):.3f}" if dts else f"delta={delta:g} {short}: 0")
    n = max(len(v) for v in cols.values())
    wide = pd.DataFrame({k: v + [np.nan] * (n - len(v)) for k, v in cols.items()})
    wide.to_csv(args.out, index=False)
    pd.DataFrame(rows).to_csv(os.path.splitext(args.out)[0] + "_per_pair.csv", index=False)
    print(f"[DONE] {args.out}  (columns: {', '.join(wide.columns)})")


if __name__ == "__main__":
    main()
