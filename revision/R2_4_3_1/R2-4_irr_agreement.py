#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
R2-4_irr_agreement.py
Inter-rater (or intra-rater) reliability for chewing-event annotations
(Nieh et al. 2026, iScience revision; reviewer comment R2-4).

Inputs are per-video wide-format CSVs produced by 01_Video_Annotator.py
or ICC_annotate.py (header: Trial, Chew1/Chew_1, Chew2/Chew_2, ...; one
row per reach/trial; cells are chew timestamps in seconds). All chew
timestamps in a file are pooled into one event series per video pass.

Three ways to specify file pairs:
  (A) RECOMMENDED one-shot mode: pair the repository's original
      annotations with the second coder's output tree automatically:
        python R2-4_irr_agreement.py --annotated_data ../annotation_raw_data \
            --irr_root IRR_final --out_prefix R2-4_irr_results
      Expected second-coder layout (as produced by ICC_annotate_v2.py):
        IRR_final/<video>/A/*.csv and IRR_final/<video>/B/*.csv  (pair videos;
            A = person on the LEFT of the frame, B = person on the RIGHT,
            same convention as the repository's A/B folders)
        IRR_final/<video>/*.csv                                  (solo videos)
      The condition (visible/invisible/solo) is inferred from which
      repository folder contains the video. The generated pairing is
      written to <out_prefix>_pairs.csv for the record.
  (B) Two directories with identical filenames:
        python R2-4_irr_agreement.py --dir_a CODER1_DIR --dir_b CODER2_DIR \
            --out_prefix R2-4_irr
  (C) Explicit manifest CSV with columns:
        file_a,file_b,label,condition
      (label = free text, e.g. video/participant id; condition optional)
        python R2-4_irr_agreement.py --pairs manifest.csv --out_prefix R2-4_irr

Options:
  --grid 0.1        rounding grid in seconds (default 0.1, the annotation
                    resolution of the primary pipeline; set 0 to disable)
  --tolerances "0.05,0.1,0.15,0.2,0.3,0.5"
  --duration 180    session duration (s), used for window rates and kappa
  --window 20       window length (s) for windowed chewing-rate ICC
  --n_boot 10000    bootstrap iterations for across-video 95% CIs
  --seed 20260823

Outputs:
  <out_prefix>_per_video.csv : per video x tolerance agreement metrics
  <out_prefix>_summary.csv   : pooled + across-video summary per tolerance
  <out_prefix>_summary.txt   : human-readable report (counts ICC, windowed
                               rate ICC, kappa, timing-discrepancy stats)

Terminology note: |dt| here is the CODER timing discrepancy between two
annotation passes of the SAME participant. It is unrelated to the
inter-participant latency Delta-t analyzed in the manuscript.
"""

import argparse
import os
import sys
import glob
import math
import numpy as np
import pandas as pd


# ---------------------------------------------------------------- I/O
def read_event_series(path, grid):
    """Read a wide-format annotation CSV; return sorted 1-D array of chew
    times (s), pooled across trial rows."""
    df = pd.read_csv(path, dtype=str)
    chew_cols = [c for c in df.columns if c.lower().startswith("chew")]
    if not chew_cols:
        raise ValueError(f"No Chew columns found in {path}")
    vals = []
    for c in chew_cols:
        v = pd.to_numeric(df[c], errors="coerce")
        vals.append(v)
    t = pd.concat(vals).dropna().to_numpy(dtype=float)
    if grid and grid > 0:
        t = np.round(t / grid) * grid
    t.sort()
    return t


def strip_participant_suffix(fname):
    """'20170706_01_02_01_A.csv' -> '20170706_01_02_01' (also _B)."""
    base = fname[:-4] if fname.lower().endswith(".csv") else fname
    if base.endswith(("_A", "_B")):
        base = base[:-2]
    return base


def find_original(dirpath, video):
    """Find the original annotation CSV for <video> in dirpath, whatever
    participant suffix it carries (_A/_B/none)."""
    hits = [p for p in glob.glob(os.path.join(dirpath, "*.csv"))
            if strip_participant_suffix(os.path.basename(p)) == video]
    if len(hits) != 1:
        sys.exit(f"ERROR: expected exactly one original CSV for '{video}' "
                 f"in {dirpath}, found {len(hits)}: {hits}")
    return hits[0]


def find_second_coder(dirpath, video):
    """Find the second coder's CSV for <video> in dirpath (ICC_annotate
    output; excludes non-CSV backups)."""
    hits = sorted(glob.glob(os.path.join(dirpath, "*.csv")))
    if len(hits) != 1:
        sys.exit(f"ERROR: expected exactly one second-coder CSV in "
                 f"{dirpath}, found {len(hits)}: {hits}")
    return hits[0]


def build_pairs_auto(annotated_data, irr_root):
    """Pair repository originals (coder A = reference) with the second
    coder's IRR tree (coder B). Returns list of pair dicts."""
    layout = [("visible", "visible_pair_Human", True),
              ("invisible", "invisible_pair_Human", True),
              ("solo", "solo_Human", False)]
    # index all repository videos by condition
    video_condition = {}
    for cond, folder, is_pair in layout:
        base = os.path.join(annotated_data, folder)
        scan = os.path.join(base, "A") if is_pair else base
        if not os.path.isdir(scan):
            sys.exit(f"ERROR: folder not found: {scan}")
        for p in glob.glob(os.path.join(scan, "*.csv")):
            v = strip_participant_suffix(os.path.basename(p))
            video_condition[v] = (cond, folder, is_pair)
    pairs = []
    videos = sorted(d for d in os.listdir(irr_root)
                    if os.path.isdir(os.path.join(irr_root, d)))
    if not videos:
        sys.exit(f"ERROR: no video folders found under {irr_root}")
    for v in videos:
        if v not in video_condition:
            sys.exit(f"ERROR: video '{v}' in {irr_root} not found in the "
                     f"repository at {annotated_data}")
        cond, folder, is_pair = video_condition[v]
        vdir = os.path.join(irr_root, v)
        if is_pair:
            for part in ("A", "B"):
                pdir = os.path.join(vdir, part)
                if not os.path.isdir(pdir):
                    sys.exit(f"ERROR: missing participant folder: {pdir}")
                pairs.append(dict(
                    file_a=find_original(
                        os.path.join(annotated_data, folder, part), v),
                    file_b=find_second_coder(pdir, v),
                    label=f"{v}_{part}", condition=cond))
        else:
            pairs.append(dict(
                file_a=find_original(
                    os.path.join(annotated_data, folder), v),
                file_b=find_second_coder(vdir, v),
                label=v, condition=cond))
    return pairs


def build_pairs_repo(annotated_data, rerun_root):
    """Pair repository originals (coder A = reference, first pass) with a
    re-annotation tree that uses the SAME repository layout
    (visible_pair_Human/{A,B}, invisible_pair_Human/{A,B}, solo_Human).
    Only videos present in the re-annotation tree are analysed."""
    layout = [("visible", "visible_pair_Human", True),
              ("invisible", "invisible_pair_Human", True),
              ("solo", "solo_Human", False)]
    pairs = []
    for cond, folder, is_pair in layout:
        base = os.path.join(rerun_root, folder)
        if not os.path.isdir(base):
            continue
        if is_pair:
            for part in ("A", "B"):
                pdir = os.path.join(base, part)
                if not os.path.isdir(pdir):
                    sys.exit(f"ERROR: missing participant folder: {pdir}")
                for p in sorted(glob.glob(os.path.join(pdir, "*.csv"))):
                    v = strip_participant_suffix(os.path.basename(p))
                    pairs.append(dict(
                        file_a=find_original(
                            os.path.join(annotated_data, folder, part), v),
                        file_b=p, label=f"{v}_{part}", condition=cond))
        else:
            for p in sorted(glob.glob(os.path.join(base, "*.csv"))):
                v = strip_participant_suffix(os.path.basename(p))
                pairs.append(dict(
                    file_a=find_original(
                        os.path.join(annotated_data, folder), v),
                    file_b=p, label=v, condition=cond))
    if not pairs:
        sys.exit(f"ERROR: no annotation CSVs found under {rerun_root}")
    return pairs


def build_pairs(args):
    """Return list of dicts: file_a, file_b, label, condition."""
    pairs = []
    if args.annotated_data:
        if args.rerun_root and args.irr_root:
            sys.exit("Use either --irr_root or --rerun_root, not both")
        if args.rerun_root:
            return build_pairs_repo(args.annotated_data, args.rerun_root)
        if not args.irr_root:
            sys.exit("--annotated_data requires --irr_root or --rerun_root")
        return build_pairs_auto(args.annotated_data, args.irr_root)
    if args.pairs:
        m = pd.read_csv(args.pairs)
        need = {"file_a", "file_b"}
        if not need.issubset(m.columns):
            sys.exit("Manifest must contain columns: file_a, file_b "
                     "(optional: label, condition)")
        for _, r in m.iterrows():
            pairs.append(dict(
                file_a=r["file_a"], file_b=r["file_b"],
                label=r.get("label", os.path.basename(str(r["file_a"]))),
                condition=r.get("condition", "NA")))
    else:
        if not (args.dir_a and args.dir_b):
            sys.exit("Provide either --pairs MANIFEST or --dir_a and --dir_b")
        fa = {os.path.basename(p): p
              for p in glob.glob(os.path.join(args.dir_a, "*.csv"))}
        fb = {os.path.basename(p): p
              for p in glob.glob(os.path.join(args.dir_b, "*.csv"))}
        common = sorted(set(fa) & set(fb))
        only_a = sorted(set(fa) - set(fb))
        only_b = sorted(set(fb) - set(fa))
        if only_a:
            print(f"[warn] {len(only_a)} file(s) only in dir_a "
                  f"(skipped): {', '.join(only_a)}")
        if only_b:
            print(f"[warn] {len(only_b)} file(s) only in dir_b "
                  f"(skipped): {', '.join(only_b)}")
        if not common:
            sys.exit("No matching filenames between the two directories.")
        for name in common:
            pairs.append(dict(file_a=fa[name], file_b=fb[name],
                              label=name, condition="NA"))
    return pairs


# ------------------------------------------------------- event matching
def match_events_full(t_a, t_b, tol):
    """Greedy chronological 1-to-1 matching of two sorted event series
    within tolerance tol (same algorithm family as the NLOPM matcher:
    each event used at most once, earliest-first).
    Returns (index_pairs, dt_array) with dt = t_b - t_a of matched pairs."""
    i, j = 0, 0
    pairs_idx = []
    dts = []
    eps = 1e-9
    while i < len(t_a) and j < len(t_b):
        dt = t_b[j] - t_a[i]
        if abs(dt) <= tol + eps:
            pairs_idx.append((i, j))
            dts.append(dt)
            i += 1
            j += 1
        elif dt < 0:
            j += 1
        else:
            i += 1
    return pairs_idx, np.asarray(dts)


def match_events(t_a, t_b, tol):
    pairs_idx, dts = match_events_full(t_a, t_b, tol)
    return len(pairs_idx), dts


def local_gap(t, k):
    """Distance (s) from event k to its nearest neighbor within the SAME
    series."""
    gaps = []
    if k > 0:
        gaps.append(t[k] - t[k - 1])
    if k < len(t) - 1:
        gaps.append(t[k + 1] - t[k])
    return min(gaps) if gaps else np.inf


def prf(n_match, n_a, n_b):
    recall = n_match / n_a if n_a else np.nan       # A = reference coder
    precision = n_match / n_b if n_b else np.nan    # B = comparison coder
    if precision and recall and (precision + recall) > 0:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = np.nan
    return precision, recall, f1


# ------------------------------------------------------------------ ICC
def icc_a1(x, y):
    """ICC(A,1): two-way random effects, absolute agreement, single
    measures (McGraw & Wong case 2A). x, y = paired vectors."""
    d = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    n, k = d.shape
    if n < 2:
        return np.nan
    mean_r = d.mean(axis=1)
    mean_c = d.mean(axis=0)
    grand = d.mean()
    msr = k * np.sum((mean_r - grand) ** 2) / (n - 1)
    msc = n * np.sum((mean_c - grand) ** 2) / (k - 1)
    sse = np.sum((d - mean_r[:, None] - mean_c[None, :] + grand) ** 2)
    mse = sse / ((n - 1) * (k - 1))
    denom = msr + (k - 1) * mse + k * (msc - mse) / n
    return (msr - mse) / denom if denom > 0 else np.nan


def pearson(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def mean_nn_lag(t_left, t_right):
    """Mean absolute nearest-neighbor lag (s) from each left-participant
    chew to the closest right-participant chew (inter-participant
    chew-time lag; measure-level diagnostic for pair videos)."""
    if len(t_left) == 0 or len(t_right) == 0:
        return np.nan
    idx = np.searchsorted(t_right, t_left)
    idx_lo = np.clip(idx - 1, 0, len(t_right) - 1)
    idx_hi = np.clip(idx, 0, len(t_right) - 1)
    d = np.minimum(np.abs(t_left - t_right[idx_lo]),
                   np.abs(t_left - t_right[idx_hi]))
    return float(np.mean(d))


def mean_ici(t):
    """Mean inter-chew interval (s) within one pass."""
    return float(np.mean(np.diff(t))) if len(t) >= 2 else np.nan


def windowed_counts(t, duration, window):
    edges = np.arange(0, duration + 1e-9, window)
    if edges[-1] < duration:
        edges = np.append(edges, duration)
    cnt, _ = np.histogram(t, bins=edges)
    return cnt


# ---------------------------------------------------------------- kappa
def cohens_kappa_binned(t_a, t_b, duration, bin_s):
    edges = np.arange(0, duration + 1e-9, bin_s)
    a = (np.histogram(t_a, bins=edges)[0] > 0).astype(int)
    b = (np.histogram(t_b, bins=edges)[0] > 0).astype(int)
    po = np.mean(a == b)
    pa1, pb1 = a.mean(), b.mean()
    pe = pa1 * pb1 + (1 - pa1) * (1 - pb1)
    return (po - pe) / (1 - pe) if pe < 1 else np.nan


# ------------------------------------------------------------ bootstrap
def boot_ci(values, n_boot, rng, stat=np.mean):
    v = np.asarray([x for x in values if np.isfinite(x)], float)
    if len(v) == 0:
        return (np.nan, np.nan)
    if len(v) == 1:
        return (v[0], v[0])
    bs = np.array([stat(rng.choice(v, size=len(v), replace=True))
                   for _ in range(n_boot)])
    return tuple(np.percentile(bs, [2.5, 97.5]))


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description="Annotation reliability (R2-4)")
    ap.add_argument("--annotated_data",
                    help="path to the repository's Annotated_Data folder "
                         "(original coder = reference); use with --irr_root")
    ap.add_argument("--irr_root",
                    help="root of the second coder's output tree "
                         "(e.g. IRR_final)")
    ap.add_argument("--rerun_root",
                    help="root of a re-annotation tree in repository layout "
                         "(visible_pair_Human/{A,B}, invisible_pair_Human/"
                         "{A,B}, solo_Human), e.g. the same coder's second "
                         "pass (ICC_result); coder B = this tree")
    ap.add_argument("--dir_a", help="directory of coder A (reference) CSVs")
    ap.add_argument("--dir_b", help="directory of coder B CSVs")
    ap.add_argument("--pairs", help="manifest CSV: file_a,file_b,label,condition")
    ap.add_argument("--out_prefix", default="R2-4_irr")
    ap.add_argument("--grid", type=float, default=0.1)
    ap.add_argument("--tolerances", default="0.05,0.1,0.15,0.2,0.3,0.5")
    ap.add_argument("--duration", type=float, default=180.0)
    ap.add_argument("--window", type=float, default=20.0)
    ap.add_argument("--kappa_bin", type=float, default=0.1)
    ap.add_argument("--n_boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260823)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    tols = sorted(float(x) for x in args.tolerances.split(","))
    pairs = build_pairs(args)
    print(f"{len(pairs)} annotation pair(s) found.")
    pd.DataFrame(pairs).to_csv(f"{args.out_prefix}_pairs.csv", index=False)

    per_rows = []
    corr_rows = []
    gran_rows = []
    counts_a, counts_b = [], []
    ici_a, ici_b = [], []
    win_a, win_b = [], []
    kappas = []
    dt_ref_all = []          # signed dt at the widest tolerance, pooled
    dt_ref_tol = max(tols)
    offsets = []
    series = {}              # label -> (t_a, t_b) for measure-level stats

    for p in pairs:
        t_a = read_event_series(p["file_a"], args.grid)
        t_b = read_event_series(p["file_b"], args.grid)
        series[p["label"]] = (t_a, t_b)
        counts_a.append(len(t_a))
        counts_b.append(len(t_b))
        ici_a.append(mean_ici(t_a))
        ici_b.append(mean_ici(t_b))
        win_a.append(windowed_counts(t_a, args.duration, args.window))
        win_b.append(windowed_counts(t_b, args.duration, args.window))
        kappas.append(cohens_kappa_binned(t_a, t_b, args.duration,
                                          args.kappa_bin))
        # signed-offset estimate for this pass: median signed dt of pairs
        # matched at the widest tolerance (dt = t_b - t_a)
        pairs_ref, dts_ref = match_events_full(t_a, t_b, dt_ref_tol)
        matched_a = {i for i, _ in pairs_ref}
        matched_b = {j for _, j in pairs_ref}
        for k in range(len(t_a)):
            gran_rows.append(dict(
                label=p["label"], condition=p["condition"], coder="original",
                t=t_a[k], matched=int(k in matched_a),
                local_gap=local_gap(t_a, k)))
        for k in range(len(t_b)):
            gran_rows.append(dict(
                label=p["label"], condition=p["condition"],
                coder="second", t=t_b[k], matched=int(k in matched_b),
                local_gap=local_gap(t_b, k)))
        offset = float(np.median(dts_ref)) if len(dts_ref) else 0.0
        offsets.append(offset)
        dt_ref_all.append(dts_ref)
        t_b_corr = t_b - offset
        for tol in tols:
            n_m, dts = match_events(t_a, t_b, tol)
            precision, recall, f1 = prf(n_m, len(t_a), len(t_b))
            per_rows.append(dict(
                label=p["label"], condition=p["condition"], tolerance=tol,
                n_events_a=len(t_a), n_events_b=len(t_b), n_matched=n_m,
                precision=precision, recall=recall, f1=f1,
                median_abs_dt=(np.median(np.abs(dts)) if n_m else np.nan),
                p95_abs_dt=(np.percentile(np.abs(dts), 95) if n_m else np.nan),
                mean_dt_signed=(float(np.mean(dts)) if n_m else np.nan),
                median_dt_signed=(float(np.median(dts)) if n_m else np.nan),
                offset_applied=offset,
            ))
            # offset-corrected (per-pass constant shift), diagnostic only
            n_mc, dtsc = match_events(t_a, t_b_corr, tol)
            pc, rc, f1c = prf(n_mc, len(t_a), len(t_b_corr))
            corr_rows.append(dict(
                label=p["label"], condition=p["condition"], tolerance=tol,
                n_events_a=len(t_a), n_events_b=len(t_b), n_matched=n_mc,
                precision=pc, recall=rc, f1=f1c,
                median_abs_dt=(np.median(np.abs(dtsc)) if n_mc else np.nan),
                offset_applied=offset,
            ))

    per = pd.DataFrame(per_rows)
    per.to_csv(f"{args.out_prefix}_per_video.csv", index=False)
    perc = pd.DataFrame(corr_rows)
    perc.to_csv(f"{args.out_prefix}_per_video_offset_corrected.csv",
                index=False)
    pd.DataFrame(gran_rows).to_csv(f"{args.out_prefix}_events.csv",
                                   index=False)

    # ---- summary per tolerance (pooled over events + across-video CI)
    sum_rows = []
    for tol in tols:
        d = per[per["tolerance"] == tol]
        n_a, n_b = d["n_events_a"].sum(), d["n_events_b"].sum()
        n_m = d["n_matched"].sum()
        precision_p, recall_p, f1_p = prf(n_m, n_a, n_b)
        lo, hi = boot_ci(d["f1"], args.n_boot, rng)
        sum_rows.append(dict(
            tolerance=tol, n_videos=len(d),
            total_events_a=n_a, total_events_b=n_b, total_matched=n_m,
            precision_pooled=precision_p, recall_pooled=recall_p,
            f1_pooled=f1_p,
            f1_video_mean=d["f1"].mean(),
            f1_video_ci_lo=lo, f1_video_ci_hi=hi))
    summ = pd.DataFrame(sum_rows)
    summ.to_csv(f"{args.out_prefix}_summary.csv", index=False)

    sum_rows_c = []
    for tol in tols:
        d = perc[perc["tolerance"] == tol]
        n_a, n_b = d["n_events_a"].sum(), d["n_events_b"].sum()
        n_m = d["n_matched"].sum()
        precision_p, recall_p, f1_p = prf(n_m, n_a, n_b)
        lo, hi = boot_ci(d["f1"], args.n_boot, rng)
        sum_rows_c.append(dict(
            tolerance=tol, n_videos=len(d),
            total_events_a=n_a, total_events_b=n_b, total_matched=n_m,
            precision_pooled=precision_p, recall_pooled=recall_p,
            f1_pooled=f1_p,
            f1_video_mean=d["f1"].mean(),
            f1_video_ci_lo=lo, f1_video_ci_hi=hi))
    summc = pd.DataFrame(sum_rows_c)
    summc.to_csv(f"{args.out_prefix}_summary_offset_corrected.csv",
                 index=False)

    # ---- measure-level reliability (paper-relevant derived quantities)
    #  1) per-pass chew counts (ICC, existing)
    #  2) per-pass mean inter-chew interval
    #  3) pair videos: mean absolute nearest-neighbor inter-participant
    #     chew-time lag, computed within each coder, compared across coders
    pair_bases = {}
    for lab in series:
        if lab.endswith(("_A", "_B")):
            pair_bases.setdefault(lab[:-2], {})[lab[-1]] = series[lab]
    lag_a, lag_b, lag_labels = [], [], []
    for base in sorted(pair_bases):
        d = pair_bases[base]
        if "A" in d and "B" in d:
            ta_left, tb_left = d["A"]
            ta_right, tb_right = d["B"]
            lag_a.append(mean_nn_lag(ta_left, ta_right))   # coder A's data
            lag_b.append(mean_nn_lag(tb_left, tb_right))   # coder B's data
            lag_labels.append(base)

    # ---- global reliability stats
    icc_counts = icc_a1(counts_a, counts_b)
    wa = np.concatenate(win_a).astype(float)
    wb = np.concatenate(win_b).astype(float)
    icc_windows = icc_a1(wa, wb)
    kap = np.asarray(kappas, float)

    dt_pool = (np.concatenate(dt_ref_all) if dt_ref_all else np.array([]))
    abs_dt = np.abs(dt_pool)

    with open(f"{args.out_prefix}_summary.txt", "w") as fh:
        fh.write("=====================================================\n")
        fh.write("R2-4 annotation reliability report\n")
        fh.write("=====================================================\n\n")
        if args.rerun_root:
            fh.write("MODE: intra-rater  (coder B = re-annotation tree "
                     f"'{args.rerun_root}', same coder, repository "
                     "layout)\n")
        elif args.irr_root:
            fh.write("MODE: inter-rater  (coder B = independent second "
                     f"coder tree '{args.irr_root}')\n")
        fh.write(f"Videos compared      : {len(pairs)}\n")
        cond_counts = pd.Series([p['condition'] for p in pairs]) \
            .value_counts().to_dict()
        fh.write(f"By condition         : {cond_counts}\n")
        fh.write(f"Rounding grid        : {args.grid} s\n")
        fh.write(f"Events coder A / B   : {sum(counts_a)} / {sum(counts_b)}\n\n")

        fh.write("Tolerance-agreement curve (pooled over events; F1 CI is\n"
                 "a bootstrap 95% CI across videos):\n")
        fh.write(summ.to_string(index=False,
                                float_format=lambda x: f"{x:.3f}"))
        fh.write("\n\n")

        fh.write(f"Coder timing discrepancy (matched pairs at tol = "
                 f"{dt_ref_tol} s; n = {len(abs_dt)}):\n")
        if len(abs_dt):
            for q, v in [("median |dt|", np.median(abs_dt)),
                         ("mean |dt|", abs_dt.mean()),
                         ("P90 |dt|", np.percentile(abs_dt, 90)),
                         ("P95 |dt|", np.percentile(abs_dt, 95))]:
                fh.write(f"  {q:12s} = {v:.3f} s\n")
            for thr in (0.0, 0.1, 0.2):
                fh.write(f"  share |dt| <= {thr:.1f} s = "
                         f"{np.mean(abs_dt <= thr + 1e-9):.3f}\n")
            fh.write("  (dt = coder discrepancy; distinct from the\n"
                     "   inter-participant latency Delta-t in the paper)\n")
        fh.write("\n")

        fh.write("Signed timing discrepancy (dt = t_secondcoder - "
                 "t_original, same matches):\n")
        if len(dt_pool):
            fh.write(f"  mean dt      = {dt_pool.mean():+.3f} s\n")
            fh.write(f"  median dt    = {np.median(dt_pool):+.3f} s\n")
            fh.write(f"  SD dt        = {dt_pool.std(ddof=1):.3f} s\n")
            fh.write(f"  share dt > 0 = {np.mean(dt_pool > 0):.3f}\n")
            off = np.asarray(offsets, float)
            fh.write(f"  per-pass median offset: mean = {off.mean():+.3f} s, "
                     f"range {off.min():+.3f} to {off.max():+.3f} s, "
                     f"share > 0 = {np.mean(off > 0):.3f}\n")
        fh.write("\n")

        fh.write("Offset-corrected tolerance-agreement curve (DIAGNOSTIC:\n"
                 "each second-coder series shifted by its per-pass median\n"
                 f"signed dt estimated at tol = {dt_ref_tol} s; "
                 "raw curve above remains primary):\n")
        fh.write(summc.to_string(index=False,
                                 float_format=lambda x: f"{x:.3f}"))
        fh.write("\n\n")

        fh.write("Granularity diagnostic (unmatched events at tol = "
                 f"{dt_ref_tol} s vs. spacing to the nearest event of the\n"
                 "SAME coder; tests whether disagreements concentrate in "
                 "fast chew bursts):\n")
        gr = pd.DataFrame(gran_rows)
        for coder in ("original", "second"):
            for m, tag in ((0, "unmatched"), (1, "matched")):
                g = gr[(gr["coder"] == coder) & (gr["matched"] == m)]
                g = g[np.isfinite(g["local_gap"])]
                if len(g) == 0:
                    continue
                lg = g["local_gap"].to_numpy()
                fh.write(f"  {coder:8s} {tag:9s} (n = {len(g):5d}): "
                         f"median local gap = {np.median(lg):.3f} s; "
                         f"share with gap <= 0.3 s = {np.mean(lg <= 0.3):.3f},"
                         f" <= 0.4 s = {np.mean(lg <= 0.4):.3f},"
                         f" <= 0.5 s = {np.mean(lg <= 0.5):.3f}\n")
        fh.write("\n")

        fh.write("Measure-level reliability (paper-relevant derived "
                 "quantities):\n")
        fh.write(f"  ICC(A,1) per-pass mean inter-chew interval "
                 f"(n = {len(ici_a)} passes) = "
                 f"{icc_a1(ici_a, ici_b):.3f}  "
                 f"(Pearson r = {pearson(ici_a, ici_b):.3f})\n")
        if lag_labels:
            fh.write(f"  Pair videos (n = {len(lag_labels)}): mean absolute "
                     "nearest-neighbor inter-participant chew-time lag,\n"
                     "  computed within each coder's own annotations:\n")
            fh.write(f"    ICC(A,1) = {icc_a1(lag_a, lag_b):.3f}  "
                     f"(Pearson r = {pearson(lag_a, lag_b):.3f})\n")
            for base, va, vb in zip(lag_labels, lag_a, lag_b):
                fh.write(f"    {base:22s} original = {va:.3f} s   "
                         f"second coder = {vb:.3f} s\n")
        fh.write("\n")

        fh.write("Event-count reliability:\n")
        fh.write(f"  ICC(A,1) per-video event counts   = {icc_counts:.3f}\n")
        fh.write(f"  ICC(A,1) {args.window:.0f}-s window counts "
                 f"(n = {len(wa)} windows) = {icc_windows:.3f}\n\n")

        fh.write(f"Cohen's kappa on {args.kappa_bin} s bins "
                 f"(reference value only):\n")
        fh.write(f"  mean kappa = {np.nanmean(kap):.3f} "
                 f"(range {np.nanmin(kap):.3f} - {np.nanmax(kap):.3f})\n\n")

        fh.write("Per-video results saved to: "
                 f"{args.out_prefix}_per_video.csv\n")

    print(f"Saved: {args.out_prefix}_per_video.csv, "
          f"{args.out_prefix}_summary.csv, {args.out_prefix}_summary.txt")


if __name__ == "__main__":
    main()
