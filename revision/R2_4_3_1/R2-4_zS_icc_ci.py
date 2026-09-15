#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
R2-4_zS_icc_ci.py  (trialfix_kit v4.3)

Confidence intervals for the inter-/intra-coder reliability of session-level zS.

Reads the per-video tables written by R2-4_zS_reproduction.py
    R2-4_zS_<mode>_d<delta>_per_video.csv   (mode = inter | intra)
(columns: video, condition, coder, period, n_a, n_b, n_match, S, zS, p)
and, for every mode x delta x period (and for the three periods pooled), reports

  * ICC(A,1)  (two-way random effects, absolute agreement, single measures;
               the same estimator as R2-4_zS_reproduction.py)
  * a bootstrap 95% CI (videos resampled with replacement, percentile method)
  * the McGraw & Wong (1996) F-based 95% CI, for reference only (it assumes
    normality and is not meaningful when the point estimate is negative)
  * Pearson r for comparison

Nothing is re-simulated: the zS values are taken as written by the reliability
run, so this step takes a few seconds and can be re-run at any time.

Usage (run by run_R2_4_trialfix.sh; STAGE=icc_ci runs only this step):
  python R2-4_zS_icc_ci.py [--n_boot 10000] [--seed 20260827]
                           [--pattern 'R2-4_zS_*_d*_per_video.csv']
                           [--out_prefix R2-4_zS_icc_ci]

Outputs:
  <out_prefix>.txt   human-readable summary (the file the letter cites)
  <out_prefix>.csv   one row per mode x delta x period
"""

import argparse
import glob
import re
import sys

import numpy as np
import pandas as pd

try:
    from scipy import stats as _st
except Exception:  # scipy is optional (only for the F-based CI)
    _st = None


# ------------------------------------------------------------------ ICC(A,1)
def icc_a1(d):
    """d: (n, 2) array of paired ratings.  Returns (icc, msr, msc, mse, n, k)."""
    d = np.asarray(d, float)
    d = d[np.all(np.isfinite(d), axis=1)]
    n, k = d.shape
    if n < 3:
        return np.nan, np.nan, np.nan, np.nan, n, k
    mean_r = d.mean(axis=1)
    mean_c = d.mean(axis=0)
    grand = d.mean()
    msr = k * np.sum((mean_r - grand) ** 2) / (n - 1)
    msc = n * np.sum((mean_c - grand) ** 2) / (k - 1)
    sse = np.sum((d - mean_r[:, None] - mean_c[None, :] + grand) ** 2)
    mse = sse / ((n - 1) * (k - 1))
    denom = msr + (k - 1) * mse + k * (msc - mse) / n
    icc = (msr - mse) / denom if denom > 0 else np.nan
    return icc, msr, msc, mse, n, k


def icc_a1_f_ci(icc, msr, msc, mse, n, k, alpha=0.05):
    """McGraw & Wong (1996) confidence interval for ICC(A,1)."""
    if _st is None or not np.isfinite(icc) or icc >= 1:
        return np.nan, np.nan
    a = k * icc / (n * (1 - icc))
    b = 1 + k * icc * (n - 1) / (n * (1 - icc))
    num = (a * msc + b * mse) ** 2
    den = (a * msc) ** 2 / (k - 1) + (b * mse) ** 2 / ((n - 1) * (k - 1))
    if den <= 0:
        return np.nan, np.nan
    v = num / den
    f_l = _st.f.ppf(1 - alpha / 2, n - 1, v)
    f_u = _st.f.ppf(1 - alpha / 2, v, n - 1)
    lo = n * (msr - f_l * mse) / (f_l * (k * msc + (k * n - k - n) * mse) + n * msr)
    hi = n * (f_u * msr - mse) / (k * msc + (k * n - k - n) * mse + n * f_u * msr)
    return lo, hi


def boot_ci(d, n_boot, rng, alpha=0.05):
    d = np.asarray(d, float)
    n = d.shape[0]
    vals = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        vals[b] = icc_a1(d[idx])[0]
    vals = vals[np.isfinite(vals)]
    if len(vals) < 100:
        return np.nan, np.nan, len(vals)
    return (float(np.percentile(vals, 100 * alpha / 2)),
            float(np.percentile(vals, 100 * (1 - alpha / 2))), len(vals))


def pearson(d):
    d = np.asarray(d, float)
    if len(d) < 3 or d[:, 0].std() == 0 or d[:, 1].std() == 0:
        return np.nan
    return float(np.corrcoef(d[:, 0], d[:, 1])[0, 1])


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pattern", default="R2-4_zS_*_d*_per_video.csv")
    ap.add_argument("--n_boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--out_prefix", default="R2-4_zS_icc_ci")
    args = ap.parse_args()

    files = sorted(glob.glob(args.pattern))
    if not files:
        sys.exit(f"[ERROR] no file matches {args.pattern!r}. Run the reliability stage first "
                 "(bash run_R2_4_trialfix.sh, or STAGE=reliability bash run_R2_4_trialfix.sh).")

    rng = np.random.default_rng(args.seed)
    rows = []
    for f in files:
        m = re.search(r"R2-4_zS_(inter|intra)_d([\d.]+)_per_video\.csv$", f)
        if not m:
            print(f"[SKIP] {f}: name does not match R2-4_zS_<mode>_d<delta>_per_video.csv")
            continue
        mode, delta = m.group(1), float(m.group(2))
        df = pd.read_csv(f)
        need = {"video", "condition", "coder", "period", "zS"}
        if not need.issubset(df.columns):
            print(f"[SKIP] {f}: missing columns {need - set(df.columns)}")
            continue
        periods = list(dict.fromkeys(df["period"]))
        blocks = [(p, df[df["period"] == p]) for p in periods] + [("all periods pooled", df)]
        for pname, sub in blocks:
            w = sub.pivot_table(index=["video", "period", "condition"], columns="coder",
                                values="zS", aggfunc="first").reset_index()
            if not {"original", "second"}.issubset(w.columns):
                continue
            d = w[["original", "second"]].to_numpy(float)
            d = d[np.all(np.isfinite(d), axis=1)]
            icc, msr, msc, mse, n, k = icc_a1(d)
            blo, bhi, nb = boot_ci(d, args.n_boot, rng)
            flo, fhi = icc_a1_f_ci(icc, msr, msc, mse, n, k)
            rows.append(dict(mode=mode, delta=delta, period=pname, n=n,
                             n_visible=int((w["condition"] == "visible").sum()),
                             n_invisible=int((w["condition"] == "invisible").sum()),
                             icc_a1=icc, boot_lo=blo, boot_hi=bhi, n_boot_ok=nb,
                             f_lo=flo, f_hi=fhi, pearson_r=pearson(d), source=f))

    out = pd.DataFrame(rows).sort_values(["mode", "delta", "period"],
                                        key=lambda s: s if s.name != "period" else s.map(
                                            lambda x: (1, x) if x.startswith("all") else (0, x)))
    out.to_csv(f"{args.out_prefix}.csv", index=False)

    with open(f"{args.out_prefix}.txt", "w", encoding="utf-8") as fh:
        fh.write("=====================================================\n")
        fh.write("R2-4 reliability of session-level zS: ICC(A,1) with 95% CIs\n")
        fh.write("=====================================================\n\n")
        fh.write(f"Input files           : {len(files)} (pattern {args.pattern})\n")
        fh.write(f"Bootstrap replicates  : {args.n_boot} (videos resampled with replacement, percentile CI)\n")
        fh.write(f"Seed                  : {args.seed}\n")
        fh.write("F-based CI            : McGraw & Wong (1996) ICC(A,1); reference only, "
                 "not meaningful for negative point estimates\n\n")
        for (mode, delta), g in out.groupby(["mode", "delta"], sort=True):
            fh.write(f"--- {mode}  delta = {delta:g} s ---\n")
            fh.write("  period               n   ICC(A,1)   boot 95% CI          F 95% CI             Pearson r\n")
            for _, r in g.iterrows():
                fh.write(f"  {r['period']:<18s} {int(r['n']):3d}   {r['icc_a1']:7.3f}   "
                         f"[{r['boot_lo']:6.2f}, {r['boot_hi']:6.2f}]     "
                         f"[{r['f_lo']:6.2f}, {r['f_hi']:6.2f}]     {r['pearson_r']:6.3f}\n")
            fh.write("\n")
        # compact range lines for the letter (per-period CIs at each delta, both modes)
        fh.write("Range of the per-period bootstrap CIs (P1-P3 only, pooled row excluded):\n")
        for (mode, delta), g in out.groupby(["mode", "delta"], sort=True):
            gg = g[~g["period"].str.startswith("all")]
            fh.write(f"  {mode} d{delta:g}: ICC point estimates {gg['icc_a1'].min():.2f} to {gg['icc_a1'].max():.2f}; "
                     f"CI lower bounds {gg['boot_lo'].min():.2f} to {gg['boot_lo'].max():.2f}; "
                     f"CI upper bounds {gg['boot_hi'].min():.2f} to {gg['boot_hi'].max():.2f}\n")
        both = out[(out["delta"] == 0.1) & (~out["period"].str.startswith("all"))]
        if len(both):
            fh.write(f"  both modes, d0.1, P1-P3: CI lower bounds min {both['boot_lo'].min():.2f}, "
                     f"CI upper bounds max {both['boot_hi'].max():.2f}\n")
        fh.write("\nNote: ICC(A,1) = two-way random effects, absolute agreement, single measures,\n"
                 "computed on the zS values written by R2-4_zS_reproduction.py (nothing re-simulated).\n"
                 "With 10-12 videos the CIs are necessarily wide; the pooled row treats the\n"
                 "video x period cells as independent units and is given for orientation only.\n")

    print(open(f"{args.out_prefix}.txt", encoding="utf-8").read())
    print(f"Saved: {args.out_prefix}.txt, {args.out_prefix}.csv")


if __name__ == "__main__":
    main()
