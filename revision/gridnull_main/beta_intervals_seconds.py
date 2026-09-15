#!/usr/bin/env python3
"""
beta_intervals_seconds.py  (verification utility; not part of the analysis pipeline)

Convert the beta_state[i] intervals reported by check_beta_intervals_v2.py
into video seconds and total them per period (P1 0-20 s, P2 20-60 s,
P3 60-180 s), so that the numbers quoted in the response letter
("... s of the 180-s trial, of which ... s fall in P3", Fig. 4e caption)
can be read directly from the output of the user's own Stan fits.

Bin -> time mapping
  06_bayes_gridnull_fit.R builds  bin_id = as.integer(factor(t_center))
  over the long table made from Micro_Analysis_Bayesian/vis_<delta>.csv and
  invis_<delta>.csv.  Hence bin i is the i-th smallest distinct t_center that
  occurs in either file.  This script re-derives that mapping from the same
  two CSVs instead of assuming an offset, and stops with an error if the
  number of distinct t_center values differs from the number of beta_state
  bins in the Stan summary.

Definitions (identical to 04_Micro_Analysis_Bayesian.Rmd / check_beta_intervals_v2.py)
  positive interval : consecutive bins with q2.5  > 0   (visible > invisible)
  negative interval : consecutive bins with q97.5 < 0   (invisible > visible)
Each bin is 1 s wide, centred on t_center (edges t_center -/+ 0.5 s).
A bin is assigned to the period that contains its t_center.

Usage
  python beta_intervals_seconds.py --summary stan_fits/summary_z_S_t_delta01_human.csv \
      --vis Micro_Analysis_Bayesian/vis_0.1.csv --invis Micro_Analysis_Bayesian/invis_0.1.csv \
      [--periods P1:0-20,P2:20-60,P3:60-180] [--out beta_intervals_seconds_d0.1.txt]
"""
import argparse
import csv
import re
import sys

import pandas as pd


def runs(idx_sorted):
    out = []
    for i in idx_sorted:
        if out and i == out[-1][1] + 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return [(a, b) for a, b in out]


def read_beta(summary_csv):
    pat = re.compile(r"^beta_state\[(\d+)\]$")
    rows = {}
    with open(summary_csv, newline="") as f:
        for r in csv.DictReader(f):
            m = pat.match(r["variable"])
            if m:
                rows[int(m.group(1))] = (float(r["q2.5"]), float(r["q97.5"]), float(r["rhat"]))
    n = len(rows)
    if set(rows) != set(range(1, n + 1)):
        sys.exit(f"{summary_csv}: beta_state indices are not contiguous 1..{n}")
    return rows


def bin_centers(vis_csv, invis_csv):
    tc = set()
    for p in (vis_csv, invis_csv):
        tc.update(pd.read_csv(p, usecols=["t_center"])["t_center"].dropna().astype(float).tolist())
    return sorted(tc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--vis", required=True)
    ap.add_argument("--invis", required=True)
    ap.add_argument("--periods", default="P1:0-20,P2:20-60,P3:60-180")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    periods = []
    for tok in a.periods.split(","):
        name, rg = tok.split(":")
        lo, hi = rg.split("-")
        periods.append((name, float(lo), float(hi)))

    def period_of(t):
        for name, lo, hi in periods:
            if lo <= t < hi or (hi == periods[-1][2] and t == hi):
                return name
        return "none"

    beta = read_beta(a.summary)
    centers = bin_centers(a.vis, a.invis)
    if len(centers) != len(beta):
        sys.exit(f"ERROR: {len(centers)} distinct t_center values in {a.vis}/{a.invis} "
                 f"but {len(beta)} beta_state bins in {a.summary}. "
                 "Check that these CSVs are the ones the Stan fit was run on.")
    t_of = {i: centers[i - 1] for i in beta}

    lines = []
    w = lines.append
    w("=====================================================================")
    w("beta_state(t) 95% credible intervals excluding zero, in video seconds")
    w("=====================================================================")
    w(f"summary : {a.summary}")
    w(f"bins    : {len(beta)}  (bin 1 -> t_center {centers[0]:.1f} s, bin {len(beta)} -> {centers[-1]:.1f} s)")
    w(f"mapping : bin i -> i-th distinct t_center of {a.vis} + {a.invis} "
      f"(same as bin_id = as.integer(factor(t_center)) in 06_bayes_gridnull_fit.R)")
    w(f"periods : {a.periods}   (bin assigned by its t_center)")
    w(f"max Rhat over beta_state: {max(v[2] for v in beta.values()):.4f}")
    w("")
    for label, sel in (("positive (q2.5 > 0; visible > invisible)", lambda lo, hi: lo > 0),
                       ("negative (q97.5 < 0; invisible > visible)", lambda lo, hi: hi < 0)):
        idx = sorted(i for i, (lo, hi, _) in beta.items() if sel(lo, hi))
        per = {name: 0 for name, _, _ in periods}
        for i in idx:
            per[period_of(t_of[i])] = per.get(period_of(t_of[i]), 0) + 1
        w(f"--- {label} ---")
        w(f"  total: {len(idx)} bins = {len(idx)} s")
        w("  per period: " + ", ".join(f"{k} {v} s" for k, v in per.items()))
        for s, e in runs(idx):
            w(f"    bins {s:3d}-{e:3d}  t_center {t_of[s]:6.1f}-{t_of[e]:6.1f} s  "
              f"(edges {t_of[s]-0.5:6.1f}-{t_of[e]+0.5:6.1f} s, {e-s+1:3d} s, "
              f"{period_of(t_of[s])}{'' if period_of(t_of[s]) == period_of(t_of[e]) else '->' + period_of(t_of[e])})")
        w("")
    txt = "\n".join(lines)
    print(txt)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write(txt + "\n")
        print(f"Saved: {a.out}")


if __name__ == "__main__":
    main()
