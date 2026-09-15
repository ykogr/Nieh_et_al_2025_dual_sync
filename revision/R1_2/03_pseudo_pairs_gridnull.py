#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
[GRID-NULL PATCH 2026-08-28] Imports the corrected null from
02_Micro_Analysis_NLOPM_periods_gridnull.py; --null_shift {grid,continuous}.

Step 4 (R1-2): Pseudo-pair (surrogate pairing) control for the onset peak.

Rationale:
  The onset peak in micro-scale synchrony may reflect the shared task structure
  (both participants start eating at the same time) rather than genuine interpersonal
  coordination. To test this, we recombine participant A of one dyad with participant
  B of a DIFFERENT dyad recorded under the same condition. If pseudo-pairs show a
  comparable onset-period z_S, the onset synchrony is attributable to common task
  onset; if real pairs exceed pseudo-pairs, it indicates dyad-specific coordination.

Method:
  - Real pairs (i==j) and pseudo pairs (i!=j, optionally subsampled) are analyzed
    with the identical period-restricted pipeline as 02_* (NLOPM + within-segment
    circular-shift permutation).
  - Alignment: by default event times are used as annotated ('none'; valid when t=0
    corresponds to task onset in every recording). If annotation origins differ
    across recordings, use --align first_event to re-zero each series at its own
    first chew (approximation of individual eating onset).

Output:
  pseudo_pairs_{delta}_{tag}.csv : one row per (A-source, B-source) x period with
    columns a_id, b_id, is_real, period, n_A, n_B, n_matched, S, z_S, p_right, p_left.
  Statistical comparison real vs pseudo is done downstream in R
  (04_LMM_condition_by_period.R, crossed random effects for a_id and b_id).

Usage example:
  python 03_pseudo_pairs.py --indir ".../visible_pair_Human" --tag visible \
    --periods "P1=0:20,P2=20:60,P3=60:180" --delta 0.1 --perms 1000 \
    --max_pseudo 300 --align none --outdir out_pseudo
"""
import argparse
import itertools
import os
import re as regex

import numpy as np
import pandas as pd

# reuse the core routines from script 02 (must sit in the same directory)
import importlib.util
_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "nlopm_periods", os.path.join(_here, "02_Micro_Analysis_NLOPM_periods_gridnull.py"))  # [GRID-NULL PATCH]
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

load_times_from_csv = _mod.load_times_from_csv
restrict_overlap = _mod.restrict_overlap
segment_perm_stats = _mod.segment_perm_stats
parse_periods = _mod.parse_periods
auto_discover_pairs = _mod.auto_discover_pairs


# [PERIOD-ORIGIN PATCH 2026-09-08] "absolute" (default): periods t0-t1 are video
# time, intersected with the pair's overlap interval; "relative": seconds since
# the first shared chew (overlap start), which is what the code did before.
PERIOD_ORIGIN = os.environ.get("NLOPM_PERIOD_ORIGIN", "absolute")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", type=str, required=True)
    ap.add_argument("--tag", type=str, default=None)
    ap.add_argument("--periods", type=str, required=True,
                    help='e.g. "P1=0:20,P2=20:60,P3=60:180"')
    ap.add_argument("--delta", type=float, default=0.1)
    ap.add_argument("--perms", type=int, default=1000,
                    help="within-segment circular-shift permutations per combo")
    ap.add_argument("--max_pseudo", type=int, default=300,
                    help="max number of pseudo combinations (random subsample); 0 = all")
    ap.add_argument("--align", choices=["none", "first_event"], default="none",
                    help="'first_event': re-zero each series at its own first chew")
    ap.add_argument("--min_events", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--null_shift", type=str, default="grid", choices=["grid", "continuous"],
                    help="[GRID-NULL PATCH] 'grid' = 0.1-s integer-multiple shifts (corrected, default); "
                         "'continuous' = U(0,T) (previous behaviour)")
    ap.add_argument("--outdir", type=str, default="out_pseudo")
    args = ap.parse_args()
    _mod.DEFAULT_NULL_SHIFT = args.null_shift  # [GRID-NULL PATCH]

    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or regex.sub(r'\W+', '_', os.path.basename(os.path.abspath(args.indir)))
    periods = parse_periods(args.periods)
    pairs = auto_discover_pairs(args.indir)
    if len(pairs) < 2:
        raise SystemExit("[ERROR] need at least 2 dyads for pseudo pairing")
    print(f"[INFO] {len(pairs)} dyads found")

    # load all series once
    A_series, B_series = {}, {}
    for pid, fa, fb in pairs:
        a = load_times_from_csv(fa)
        b = load_times_from_csv(fb)
        if args.align == "first_event":
            a = a - a[0] if len(a) else a
            b = b - b[0] if len(b) else b
        A_series[pid] = a
        B_series[pid] = b

    ids = [pid for pid, _, _ in pairs]
    real_combos = [(i, i) for i in ids]
    pseudo_combos = [(i, j) for i, j in itertools.product(ids, ids) if i != j]
    rng = np.random.default_rng(args.seed)
    if args.max_pseudo and len(pseudo_combos) > args.max_pseudo:
        idx = rng.choice(len(pseudo_combos), size=args.max_pseudo, replace=False)
        pseudo_combos = [pseudo_combos[k] for k in sorted(idx)]
    combos = [(i, j, 1) for i, j in real_combos] + [(i, j, 0) for i, j in pseudo_combos]
    print(f"[INFO] analyzing {len(real_combos)} real + {len(pseudo_combos)} pseudo combos")

    rows = []
    for n, (ai, bj, is_real) in enumerate(combos, 1):
        a2, b2, start, end = restrict_overlap(A_series[ai], B_series[bj])
        for label, t0, t1 in periods:
            if PERIOD_ORIGIN == 'absolute':  # [PERIOD-ORIGIN PATCH]
                seg_start = max(float(t0), start)
                seg_end = min(float(t1), end)
            else:
                seg_start = start + t0
                seg_end = min(start + t1, end)
            T_seg = seg_end - seg_start
            if T_seg <= 0:
                rows.append({"a_id": ai, "b_id": bj, "is_real": is_real,
                             "condition": tag, "delta": args.delta, "period": label,
                             "t0_rel": t0, "t1_rel": t1, "T_seg": T_seg,
                             "n_A": 0, "n_B": 0, "n_matched": 0, "S": np.nan,
                             "z_S": np.nan, "p_right": np.nan, "p_left": np.nan,
                             "sufficient": 0})
                continue
            aa = a2[(a2 >= seg_start) & (a2 <= seg_end)]
            bb = b2[(b2 >= seg_start) & (b2 <= seg_end)]
            st = segment_perm_stats(aa, bb, seg_start, seg_end, args.delta,
                                    args.perms, seed=args.seed)
            rows.append({"a_id": ai, "b_id": bj, "is_real": is_real,
                         "condition": tag, "delta": args.delta, "period": label,
                         "t0_rel": t0, "t1_rel": t1, "T_seg": T_seg,
                         "n_A": len(aa), "n_B": len(bb),
                         "n_matched": st["n_matched"], "S": st["S_obs"],
                         "z_S": st["z_S"], "p_right": st["p_right"],
                         "p_left": st["p_left"],
                         "sufficient": int(min(len(aa), len(bb)) >= args.min_events)})
        if n % 20 == 0:
            print(f"[INFO] {n}/{len(combos)} combos done")

    out = os.path.join(args.outdir, f"pseudo_pairs_{args.delta}_{tag}.csv")
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"[DONE] {out}")

    # quick descriptive check (formal test is done in R with crossed random effects)
    print("\n[Descriptive summary: mean z_S by period x is_real]")
    print(df[df["sufficient"] == 1]
          .groupby(["period", "is_real"])["z_S"]
          .agg(["mean", "std", "count"]).round(3).to_string())


if __name__ == "__main__":
    main()
