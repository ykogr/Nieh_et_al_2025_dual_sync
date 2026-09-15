#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
R2-4_irr_sample_selection.py

Stratified pseudo-random selection of videos for the inter-rater
reliability (IRR) re-annotation subset (Reviewer 2, comment 4).

Design
------
- Input : EITHER (recommended)
            --annotated_data <path to the repository's Annotated_Data
            directory>  -- the manifest of all candidate videos is then
            built automatically from the folder structure
            (visible_pair_Human/{A,B}, invisible_pair_Human/{A,B},
            solo_Human) and saved as <out_prefix>_manifest.csv;
          OR
            --manifest <CSV> listing ALL candidate videos with columns:
            file,condition   (condition in {visible, invisible, solo};
            case-insensitive; extra columns are ignored).
- Output: (1) <out_prefix>_selection.csv  -- the ordered annotation list
          (2) <out_prefix>_selection_log.txt -- seed, allocation, and a
              reporting sentence for STAR Methods.

Key properties
--------------
1. Stratified: a fixed number of videos is drawn per condition
   (default 4 visible / 3 invisible / 3 solo; override with flags).
2. Reproducible: a fixed random seed (default 20260823) is used and
   recorded in the log. Anyone can regenerate the identical list.
3. Truncation-robust ordering: the selected videos are interleaved by
   cycling through conditions (the cycle order itself is randomized),
   so if annotation stops early, any completed prefix remains a
   stratified random sample with near-balanced conditions.
4. Optional exclusions (e.g., videos unsuitable for IRR) via --exclude,
   a text file with one video name per line; exclusions are logged.

Usage
-----
  python R2-4_irr_sample_selection.py \
      --annotated_data /path/to/Nieh_et_al_2025_dual_sync/Annotated_Data \
      --exclude_videos 20170710_07_08_01,20170710_07_08_02,20170710_07_08_03 \
      --out_prefix R2-4_irr_sample
  python R2-4_irr_sample_selection.py --manifest all_videos.csv \
      --n_visible 4 --n_invisible 3 --n_solo 3 --seed 20260823 \
      --exclude excluded_videos.txt --out_prefix R2-4_irr_sample
"""

import argparse
import csv
import os
import random
import sys

VALID_CONDITIONS = ("visible", "invisible", "solo")


def strip_participant_suffix(fname):
    """'20170706_01_02_01_A.csv' -> '20170706_01_02_01' (also _B)."""
    base = fname[:-4] if fname.lower().endswith(".csv") else fname
    if base.endswith(("_A", "_B")):
        base = base[:-2]
    return base


def build_manifest_from_annotated_data(root):
    """Scan the repository's Annotated_Data folder and return
    [{'file':..., 'condition':...}, ...] at the session-video level."""
    layout = [("visible", "visible_pair_Human", ("A", "B")),
              ("invisible", "invisible_pair_Human", ("A", "B")),
              ("solo", "solo_Human", None)]
    rows = []
    for cond, folder, subs in layout:
        fdir = os.path.join(root, folder)
        if not os.path.isdir(fdir):
            sys.exit(f"ERROR: folder not found: {fdir}")
        if subs is None:
            names = {strip_participant_suffix(f)
                     for f in os.listdir(fdir) if f.lower().endswith(".csv")}
            rows += [{"file": v, "condition": cond} for v in sorted(names)]
        else:
            per_sub = {}
            for s in subs:
                sdir = os.path.join(fdir, s)
                if not os.path.isdir(sdir):
                    sys.exit(f"ERROR: folder not found: {sdir}")
                per_sub[s] = {strip_participant_suffix(f)
                              for f in os.listdir(sdir)
                              if f.lower().endswith(".csv")}
            if per_sub["A"] != per_sub["B"]:
                diff = sorted(per_sub["A"] ^ per_sub["B"])
                sys.exit(f"ERROR: A/B file mismatch in {folder}: {diff}")
            rows += [{"file": v, "condition": cond}
                     for v in sorted(per_sub["A"])]
    if not rows:
        sys.exit(f"ERROR: no CSV annotation files found under {root}")
    return rows


def read_manifest(path):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        cols = [c.strip().lower() for c in (reader.fieldnames or [])]
        if "file" not in cols or "condition" not in cols:
            sys.exit("ERROR: manifest must have columns: file,condition")
        for raw in reader:
            row = {k.strip().lower(): (v or "").strip() for k, v in raw.items()}
            if not row.get("file"):
                continue
            cond = row["condition"].lower()
            if cond not in VALID_CONDITIONS:
                sys.exit(
                    f"ERROR: unknown condition '{row['condition']}' for "
                    f"file '{row['file']}' (expected visible/invisible/solo)"
                )
            rows.append({"file": row["file"], "condition": cond})
    if not rows:
        sys.exit("ERROR: manifest is empty.")
    files = [r["file"] for r in rows]
    dup = {f for f in files if files.count(f) > 1}
    if dup:
        sys.exit(f"ERROR: duplicated file names in manifest: {sorted(dup)}")
    return rows


def read_exclude(path):
    if not path:
        return set()
    with open(path, encoding="utf-8-sig") as f:
        return {ln.strip() for ln in f if ln.strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--annotated_data",
                     help="path to the repository's Annotated_Data folder "
                          "(manifest is built automatically)")
    src.add_argument("--manifest",
                     help="CSV listing all candidate videos: file,condition")
    ap.add_argument("--n_visible", type=int, default=4)
    ap.add_argument("--n_invisible", type=int, default=3)
    ap.add_argument("--n_solo", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20260823,
                    help="random seed (recorded in the log; default 20260823)")
    ap.add_argument("--exclude", default=None,
                    help="optional text file, one video name per line")
    ap.add_argument("--exclude_videos", default=None,
                    help="optional comma-separated video names to exclude, "
                         "e.g. 20170710_07_08_01,20170710_07_08_02")
    ap.add_argument("--out_prefix", default="R2-4_irr_sample")
    args = ap.parse_args()

    if args.annotated_data:
        rows = build_manifest_from_annotated_data(args.annotated_data)
        man_path = f"{args.out_prefix}_manifest.csv"
        with open(man_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["file", "condition"])
            for r in rows:
                w.writerow([r["file"], r["condition"]])
        print(f"Manifest built from {args.annotated_data} "
              f"({len(rows)} videos) -> {man_path}")
        manifest_label = f"auto-built from {args.annotated_data}"
    else:
        rows = read_manifest(args.manifest)
        manifest_label = args.manifest

    excluded = read_exclude(args.exclude)
    if args.exclude_videos:
        excluded |= {v.strip() for v in args.exclude_videos.split(",")
                     if v.strip()}

    pool = {c: [] for c in VALID_CONDITIONS}
    n_excluded = 0
    for r in rows:
        if r["file"] in excluded:
            n_excluded += 1
            continue
        pool[r["condition"]].append(r["file"])

    want = {"visible": args.n_visible,
            "invisible": args.n_invisible,
            "solo": args.n_solo}
    for c in VALID_CONDITIONS:
        if want[c] > len(pool[c]):
            sys.exit(
                f"ERROR: requested {want[c]} {c} videos but only "
                f"{len(pool[c])} are available after exclusions."
            )
        if want[c] < 0:
            sys.exit("ERROR: negative sample sizes are not allowed.")

    rng = random.Random(args.seed)

    # 1) stratified draw (sample from a sorted pool so the result does not
    #    depend on the row order of the manifest)
    drawn = {}
    for c in VALID_CONDITIONS:
        drawn[c] = rng.sample(sorted(pool[c]), want[c])
        rng.shuffle(drawn[c])

    # 2) truncation-robust interleaving: randomize the condition cycle
    #    order, then take one video per condition per cycle until all
    #    strata are exhausted.
    cycle = [c for c in VALID_CONDITIONS if want[c] > 0]
    rng.shuffle(cycle)
    queues = {c: list(drawn[c]) for c in cycle}
    ordered = []
    while any(queues[c] for c in cycle):
        for c in cycle:
            if queues[c]:
                ordered.append((queues[c].pop(0), c))

    sel_path = f"{args.out_prefix}_selection.csv"
    with open(sel_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["order", "file", "condition"])
        for i, (fname, cond) in enumerate(ordered, 1):
            w.writerow([i, fname, cond])

    log_path = f"{args.out_prefix}_selection_log.txt"
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("R2-4 IRR subset selection log\n")
        f.write("=============================\n")
        f.write(f"Manifest              : {manifest_label}\n")
        f.write(f"Random seed           : {args.seed}\n")
        f.write(f"Candidates (after excl.): "
                f"{sum(len(pool[c]) for c in VALID_CONDITIONS)} "
                f"({ {c: len(pool[c]) for c in VALID_CONDITIONS} })\n")
        f.write(f"Excluded videos       : {n_excluded} "
                f"({sorted(excluded) if excluded else 'none'})\n")
        f.write(f"Requested allocation  : {want}\n")
        f.write(f"Condition cycle order : {cycle}\n")
        f.write("Selected (annotation order):\n")
        for i, (fname, cond) in enumerate(ordered, 1):
            f.write(f"  {i:2d}. {fname}  [{cond}]\n")
        f.write(
            "\nSTAR Methods reporting sentence (fill in as appropriate):\n"
            "  'A stratified pseudo-random subset of videos "
            f"({want['visible']} visible, {want['invisible']} invisible, "
            f"{want['solo']} solo) was selected from all session videos "
            f"using a fixed random seed ({args.seed}); videos were "
            "annotated in a condition-interleaved order so that an early "
            "stop would preserve a stratified random sample.'\n"
        )

    print(f"Saved: {sel_path}, {log_path}")
    for i, (fname, cond) in enumerate(ordered, 1):
        print(f"{i:2d}. {fname}  [{cond}]")


if __name__ == "__main__":
    main()
