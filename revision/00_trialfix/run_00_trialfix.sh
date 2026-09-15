#!/usr/bin/env bash
# =====================================================================
# run_00_trialfix.sh  -- apply the Trial-column contamination fix (data + code) in one go
#
# Location : iscience_revision_analysis/00_trialfix/   (place this whole folder there)
# Usage    : cd iscience_revision_analysis/00_trialfix
#            bash run_00_trialfix.sh              # run
#            bash run_00_trialfix.sh --dry_run    # only show what would change (nothing is rewritten)
#
# What it does
#   step -1 (v5, 2026-09-15) normalise the legacy header of ../annotation_raw_data/solo_Human/20170706_07_03_B.csv
#           (trial-number column -> Trial, k-th-bite columns -> Chew_k; data rows untouched; original saved to legacy_header_backup/)
#           report: legacy_header_report.txt.  Re-running does nothing once the header is standard.
#   step 0  ../annotation_raw_data  ->  ../annotation_raw_data_chewonly
#           (copy keeping only the Chew_* columns. Apart from step -1 the original data are not touched)
#           report: chewonly_report.csv / chewonly_report.txt
#   step 1  rewrite the .py files in the sibling folders (gridnull_main, R1_2, R2_2_..., R2_4_3_1, R2_8_...)
#           that contain load_times_from_csv / NLOPM_match / nlopm_matches /
#           sliding_window_with_dt
#           (the original file is kept as <name>.pre_trialfix.bak)
#   step 2  verification (verify_report.txt, verify_pairs.csv)
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"
ROOT=".."
DRY=""
[ "${1:-}" = "--dry_run" ] && DRY="--dry_run"

# Target folders: pass all sibling folders except data and backups (nothing happens if there are no matching files)
FOLDERS=()
for d in "$ROOT"/*/; do
  d="${d%/}"; b="$(basename "$d")"
  case "$b" in
    annotation_raw_data*|__pycache__|00_trialfix|Micro_Analysis_LMM|Micro_Analysis_Bayesian|k_ICC|*backup*) continue;;
  esac
  FOLDERS+=("$d")
done
echo "[INFO] patch targets: ${FOLDERS[*]}"

echo "================ step -1: legacy header normalisation ========"
python 00_fix_legacy_header.py --src "$ROOT/annotation_raw_data" $DRY

echo "================ step 0: chew-only data copy ================"
if [ -n "$DRY" ]; then
  echo "(dry run) python 00_make_chewonly_data.py --src $ROOT/annotation_raw_data --dst $ROOT/annotation_raw_data_chewonly"
else
  python 00_make_chewonly_data.py --src "$ROOT/annotation_raw_data" --dst "$ROOT/annotation_raw_data_chewonly" --overwrite
fi

echo "================ step 1: patch scripts ======================"
python 01_patch_scripts.py $DRY "${FOLDERS[@]}"

if [ -z "$DRY" ]; then
  echo "================ step 2: verify ============================="
  MAIN="$ROOT/gridnull_main/03_Micro_Analysis_NOLPM_gridnull_v2.py"
  [ -f "$MAIN" ] || MAIN="$ROOT/gridnull_main/03_Micro_Analysis_NOLPM_gridnull.py"
  if [ -f "$MAIN" ]; then
    python 02_verify_fix.py --raw "$ROOT/annotation_raw_data" --chewonly "$ROOT/annotation_raw_data_chewonly" \
      --folders "${FOLDERS[@]}" --script "$MAIN"
  else
    python 02_verify_fix.py --raw "$ROOT/annotation_raw_data" --chewonly "$ROOT/annotation_raw_data_chewonly" \
      --folders "${FOLDERS[@]}"
  fi
  echo
  echo "[NEXT] Run the run_*_trialfix.sh in each folder in the order given in the README"
fi
