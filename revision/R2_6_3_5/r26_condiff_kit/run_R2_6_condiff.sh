#!/usr/bin/env bash
# =====================================================================
# run_R2_6_condiff.sh  (r26_condiff_kit v2.2) -- R2_6_3_5: post-fit analyses on the saved macro fits
#
# Layout (v2.1): the kit may sit either directly in R2_6_3_5 (next to 03/04/05) or as the sub-folder
# R2_6_3_5/r26_condiff_kit/ (unzipped as is).  ANALYSIS_DIR (= the R2_6_3_5 folder) is found automatically:
# the first of <this folder>, <parent folder> that contains sensitivity_out/.  Override with ANALYSIS_DIR=.
# All default paths below are relative to ANALYSIS_DIR, so the outputs land in R2_6_3_5/sensitivity_out and
# R2_6_3_5/xscale_out in both layouts (the reported values are read from these locations).
#
#   STAGE=condiff  bash run_R2_6_condiff.sh   # 06: condition difference posterior (R2-6 / R3-5)
#   STAGE=rhat     bash run_R2_6_condiff.sh   # 07: observed max R-hat of the main (raw) fits  (co-author item 10)
#   STAGE=xscale   bash run_R2_6_condiff.sh   # 08: dyad-level micro P3 z_S vs macro coupling (co-author, R2-3)
#   STAGE=all      bash run_R2_6_condiff.sh   # 06 -> 07 -> 08   (default)
#
# Environment variables (all optional):
#   ANALYSIS_DIR  the R2_6_3_5 folder                          default: auto (see above)
#   FIT_DIR   folder with fit_<mode>_<condition>.rds          default: $ANALYSIS_DIR/sensitivity_out
#   MODE      raw | detrend | diff for STAGE=rhat            default: raw
#   DATA_DIR  Annotated_Data folder ({visible,invisible}_pair_Human/{A,B}/*.csv)
#             default: $ANALYSIS_DIR/../../Annotated_Data  (then $ANALYSIS_DIR/../Annotated_Data; $HOME/... is NOT guessed:
#             if none exists the script stops and asks for DATA_DIR=)
#   MICRO_DIR folder with out_periods_grid_visible/ and out_periods_grid_invisible/   default: $ANALYSIS_DIR/../R1_2
#   XSCALE_OUT output folder of stage xscale                 default: $ANALYSIS_DIR/xscale_out
#   DELTA     micro delta for stage xscale                   default: 0.1
#   N_PERM    permutations for stage xscale                  default: 10000
# =====================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

log() { printf '[run_R2_6_condiff] %s\n' "$*"; }

# ---- locate the R2_6_3_5 folder (v2.1) ----------------------------------------------------------
if [ -z "${ANALYSIS_DIR:-}" ]; then
  if   [ -d "$HERE/sensitivity_out" ];    then ANALYSIS_DIR="$HERE"
  elif [ -d "$HERE/../sensitivity_out" ]; then ANALYSIS_DIR="$(cd "$HERE/.." && pwd)"
  else
    log "STOP: sensitivity_out/ exists neither in this folder nor in the parent folder."
    log "      Place this kit inside R2_6_3_5 (directly, or as R2_6_3_5/r26_condiff_kit/)."
    log "      If it is elsewhere, add ANALYSIS_DIR=<R2_6_3_5 folder>."
    exit 1
  fi
fi
ANALYSIS_DIR="$(cd "$ANALYSIS_DIR" && pwd)"
log "ANALYSIS_DIR (R2_6_3_5) = $ANALYSIS_DIR"

STAGE="${STAGE:-all}"
FIT_DIR="${FIT_DIR:-$ANALYSIS_DIR/sensitivity_out}"
MODE="${MODE:-raw}"
MICRO_DIR="${MICRO_DIR:-$ANALYSIS_DIR/../R1_2}"
XSCALE_OUT="${XSCALE_OUT:-$ANALYSIS_DIR/xscale_out}"
DELTA="${DELTA:-0.1}"
N_PERM="${N_PERM:-10000}"

stage_condiff() {
  log "STAGE condiff : Rscript 06_condition_difference_posterior.R $FIT_DIR"
  Rscript 06_condition_difference_posterior.R "$FIT_DIR"
}

stage_rhat() {
  for c in visible_pair invisible_pair; do
    f="$FIT_DIR/fit_${MODE}_${c}.rds"
    if [ ! -f "$f" ]; then
      log "STOP: $f does not exist. A saved fit from running 03_sensitivity_drift_beta.R with mode=$MODE is required."
      log "      If it is elsewhere, add FIT_DIR=<folder>."
      exit 1
    fi
  done
  log "STAGE rhat : Rscript 07_rhat_main.R $FIT_DIR $MODE $FIT_DIR"
  Rscript 07_rhat_main.R "$FIT_DIR" "$MODE" "$FIT_DIR"
  log "-> $FIT_DIR/rhat_${MODE}.txt  (the observed R-hat values in the manuscript STAR Methods come from this file)"
}

resolve_data_dir() {
  if [ -n "${DATA_DIR:-}" ]; then
    [ -d "$DATA_DIR/visible_pair_Human" ] || { log "STOP: no visible_pair_Human in DATA_DIR=$DATA_DIR"; exit 1; }
    return
  fi
  for cand in "$ANALYSIS_DIR/../../Annotated_Data" "$ANALYSIS_DIR/../Annotated_Data" "$ANALYSIS_DIR/../../../Annotated_Data"; do
    if [ -d "$cand/visible_pair_Human" ]; then DATA_DIR="$cand"; return; fi
  done
  log "STOP: Annotated_Data not found. Run with DATA_DIR=<…/Annotated_Data>."
  log "      (it must contain visible_pair_Human/A, visible_pair_Human/B, invisible_pair_Human/A, invisible_pair_Human/B)"
  exit 1
}

stage_xscale() {
  resolve_data_dir
  for c in visible invisible; do
    f="$MICRO_DIR/out_periods_grid_${c}/chew_sync_summary_periods_master.csv"
    [ -f "$f" ] || { log "STOP: $f does not exist. Check MICRO_DIR=<R1_2 folder>."; exit 1; }
  done
  log "STAGE xscale : Rscript 08_xscale_micro_macro.R $DATA_DIR $MICRO_DIR $XSCALE_OUT $DELTA $N_PERM"
  Rscript 08_xscale_micro_macro.R "$DATA_DIR" "$MICRO_DIR" "$XSCALE_OUT" "$DELTA" "$N_PERM"
  log "-> $XSCALE_OUT/xscale_micro_macro.txt / .csv / .pdf  (the numbers in the manuscript and letter come from this file)"
}

case "$STAGE" in
  condiff) stage_condiff ;;
  rhat)    stage_rhat ;;
  xscale)  stage_xscale ;;
  all)     stage_condiff; stage_rhat; stage_xscale ;;
  *) log "unknown STAGE=$STAGE (condiff | rhat | xscale | all)"; exit 1 ;;
esac
log "done (STAGE=$STAGE)"
