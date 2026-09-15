#!/usr/bin/env python
"""check_beta_intervals_v2.py (verification utility; not part of the analysis pipeline)

List the contiguous bins in which the 95% credible interval of beta_state(t)
excludes zero, from Stan summary CSVs (summary_z_S_t_delta*_human.csv),
using exactly the same definition as the plotting code in
04_Micro_Analysis_Bayesian.Rmd:
  - positive: consecutive bins with q2.5  > 0 (Rmd: sigpos + rle)
  - negative: consecutive bins with q97.5 < 0 (Rmd: signeg + rle)
Bins are reported as the index i of beta_state[i] (1-based). No conversion
to seconds is performed here; with the confirmed timeline mapping,
bin i corresponds to t_center = i + 2.5 s (bin 1 = 3.5 s, bin 174 = 176.5 s,
1-s bins; t_center 0 and 180 rows of the timeline CSVs are all-missing and
are dropped by the model input preparation).

v2: all output and comments in English (v1 printed Japanese; superseded).

Example (run locally in the directory containing the CSVs):
  python check_beta_intervals_v2.py summary_z_S_t_delta01_human.csv summary_z_S_t_delta02_human.csv
"""
import csv
import re
import sys


def runs(idx_sorted):
    """Collapse a sorted list of 1-based bin indices into [(start, end), ...] (R rle equivalent)."""
    out = []
    for i in idx_sorted:
        if out and i == out[-1][1] + 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return [(a, b) for a, b in out]


def main(paths):
    pat = re.compile(r"^beta_state\[(\d+)\]$")
    for p in paths:
        rows = {}
        lp_rhat = None
        with open(p, newline="") as f:
            for r in csv.DictReader(f):
                v = r["variable"]
                if v == "lp__":
                    lp_rhat = float(r["rhat"])
                m = pat.match(v)
                if not m:
                    continue
                i = int(m.group(1))
                rows[i] = (float(r["q2.5"]), float(r["q97.5"]),
                           float(r["rhat"]), float(r["ess_bulk"]), float(r["ess_tail"]))
        n = len(rows)
        assert set(rows) == set(range(1, n + 1)), f"{p}: bin indices are not contiguous 1..{n}"
        pos = sorted(i for i, (lo, hi, rh, eb, et) in rows.items() if lo > 0)
        neg = sorted(i for i, (lo, hi, rh, eb, et) in rows.items() if hi < 0)
        max_rhat = max(rh for _, _, rh, _, _ in rows.values())
        min_ess_bulk = min(eb for _, _, _, eb, _ in rows.values())
        min_ess_tail = min(et for _, _, _, _, et in rows.values())
        print(f"== {p} ==")
        print(f"  beta_state bins: {n}")
        print(f"  beta_state max Rhat: {max_rhat:.4f}   "
              f"min ess_bulk: {min_ess_bulk:.0f}   min ess_tail: {min_ess_tail:.0f}   "
              f"lp__ Rhat: {lp_rhat}")
        print(f"  positive (q2.5 > 0), {len(pos)} bins: {runs(pos)}")
        print(f"  negative (q97.5 < 0), {len(neg)} bins: {runs(neg)}")


if __name__ == "__main__":
    main(sys.argv[1:])
