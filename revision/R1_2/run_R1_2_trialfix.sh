#!/usr/bin/env bash
# =====================================================================
# run_R1_2_trialfix.sh
#   Redo all of R1-2 (period-wise analysis P1/P2/P3, pseudo pairs, latency, sensitivity analyses)
#   with the data cleaned of the Trial-column contamination (../annotation_raw_data_chewonly)
#   and the fixed scripts.
#
# Location     : iscience_revision_analysis/R1_2/
# Prerequisite : 00_trialfix/run_00_trialfix.sh has been run
#            (PATCH markers in 02_Micro_Analysis_NLOPM_periods_gridnull.py)
#            06_sensitivity_analyses_gridnull.sh has been replaced by the version bundled
#            in this kit (VIS_DIR/INVIS_DIR point to _chewonly)
# Usage    : bash run_R1_2_trialfix.sh          # python + R (if Rscript is available)
#            PERMS=1000 bash run_R1_2_trialfix.sh   # smoke test
#            SKIP_R=1 bash run_R1_2_trialfix.sh     # skip R
#            RELATIVE=0 bash run_R1_2_trialfix.sh   # skip sensitivity analysis [C] (relative-time origin)
#            STAGE=influence bash run_R1_2_trialfix.sh
#              -> (v4.4) run only 08_R1-2_p3_dyad_influence.R (leave-one-dyad-out of the P3 condition difference,
#                 sign distribution of per-dyad z_S, figure). Only reads the existing out_periods_grid_*;
#                 the permutation computation is not redone. Output dyad_influence_grid/. Requires R.
#
# Period origin: default is video time (P1 = video 0-20 s intersected with the overlap interval of each pair).
#   As sensitivity analysis [C], a relative-time version with the origin at the "first shared chew"
#   (NLOPM_PERIOD_ORIGIN=relative) is written to out_periods_grid_<cond>_relorigin /
#   out_pseudo_grid_<cond>_relorigin, and lmm_grid_out_relorigin is built with the same R scripts.
#   This shows that the P1 conclusion does not depend on the choice of origin.
#
# Output folder names are unchanged (out_periods_grid_*, out_pseudo_grid_*,
# lmm_grid_out, dt_grid_out, pseudo_interaction_grid, ..._excl_dyad, _minus5, _plus5),
# so the 04/05/07 R scripts work as is. Old outputs are moved aside.
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

RAW="${RAW:-../annotation_raw_data_chewonly}"
PERMS="${PERMS:-10000}"
PERIODS="${PERIODS:-P1=0:20,P2=20:60,P3=60:180}"
DELTA="${DELTA:-0.1}"
SKIP_R="${SKIP_R:-0}"
RELATIVE="${RELATIVE:-1}"
STAGE="${STAGE:-all}"          # all | influence
STAMP=$(date +%Y%m%d_%H%M)

# ---- STAGE=influence (v4.4): run only the dyad influence analysis and exit --------------------------
if [ "$STAGE" = "influence" ]; then
  [ -f 08_R1-2_p3_dyad_influence.R ] || { echo "[ERROR] 08_R1-2_p3_dyad_influence.R does not exist. Copy it from R1_2/ of trialfix_kit v4.4"; exit 1; }
  command -v Rscript >/dev/null 2>&1 || { echo "[ERROR] Rscript not found"; exit 1; }
  for c in visible invisible; do
    [ -f "out_periods_grid_${c}/chew_sync_summary_periods_master.csv" ] || { echo "[ERROR] out_periods_grid_${c}/chew_sync_summary_periods_master.csv does not exist. Run STAGE=all first"; exit 1; }
  done
  if [ -d dyad_influence_grid ]; then BK2="pre_redo_backup_${STAMP}"; mkdir -p "$BK2"; mv dyad_influence_grid "$BK2/"; echo "[INFO] old dyad_influence_grid/ moved aside to $BK2/"; fi
  echo "[RUN] Rscript 08_R1-2_p3_dyad_influence.R out_periods_grid_visible,out_periods_grid_invisible dyad_influence_grid P3"
  Rscript 08_R1-2_p3_dyad_influence.R out_periods_grid_visible,out_periods_grid_invisible dyad_influence_grid P3
  echo "[DONE] STAGE=influence -> dyad_influence_grid/{dyad_influence_results.txt, loo_estimates.csv, dyad_level_zS.csv, sign_tests.csv, fig_dyad_influence_delta01.pdf/.png}"
  echo "       Upload: dyad_influence_results.txt, loo_estimates.csv, sign_tests.csv, fig_dyad_influence_delta01.png"
  exit 0
fi
BK="pre_trialfix_backup_${STAMP}"

[ -d "$RAW" ] || { echo "[ERROR] $RAW does not exist (run 00_trialfix first)"; exit 1; }
grep -q "EVENT-COLUMN PATCH" 02_Micro_Analysis_NLOPM_periods_gridnull.py || { echo "[ERROR] 02_..._gridnull.py is not patched"; exit 1; }
grep -q "PERIOD-ORIGIN PATCH" 02_Micro_Analysis_NLOPM_periods_gridnull.py || { echo "[ERROR] 02_..._gridnull.py lacks the period-origin fix (re-run 00_trialfix)"; exit 1; }
grep -q "PERIOD-ORIGIN PATCH" 03_pseudo_pairs_gridnull.py || { echo "[ERROR] 03_pseudo_pairs_gridnull.py lacks the period-origin fix (re-run 00_trialfix)"; exit 1; }
grep -q "annotation_raw_data_chewonly" 06_sensitivity_analyses_gridnull.sh || { echo "[ERROR] 06_sensitivity_analyses_gridnull.sh is the old version (data path not fixed)"; exit 1; }

# ---- backup old outputs -------------------------------------------------
mkdir -p "$BK"
shopt -s nullglob
for d in out_periods_grid_* out_pseudo_grid_* lmm_grid_out* dt_grid_out* pseudo_interaction_grid* dyad_influence_grid*; do
  [ -d "$d" ] && mv "$d" "$BK/"
done
shopt -u nullglob
echo "[INFO] old outputs moved aside to $BK/"

# ---- 02: period-wise NLOPM (grid null) ------------------------------------
for cond in visible invisible; do
  echo "=============================================================="
  echo "[RUN] 02 periods  cond=$cond delta=$DELTA perms=$PERMS"
  python 02_Micro_Analysis_NLOPM_periods_gridnull.py \
    --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
    --periods "$PERIODS" --delta "$DELTA" --perms "$PERMS" \
    --episodes --winperms 300 --seed 0 \
    --outdir "out_periods_grid_${cond}"
done

# ---- 03: pseudo pairs ----------------------------------------------------
for cond in visible invisible; do
  echo "[RUN] 03 pseudo pairs  cond=$cond"
  python 03_pseudo_pairs_gridnull.py \
    --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
    --periods "$PERIODS" --delta "$DELTA" --perms "$PERMS" \
    --max_pseudo 300 --seed 0 \
    --outdir "out_pseudo_grid_${cond}"
done

# ---- [C] sensitivity: period origin = first shared chew (relative) ------------
if [ "$RELATIVE" = "1" ]; then
  for cond in visible invisible; do
    echo "[RUN] 02 periods (relative origin)  cond=$cond"
    NLOPM_PERIOD_ORIGIN=relative python 02_Micro_Analysis_NLOPM_periods_gridnull.py \
      --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
      --periods "$PERIODS" --delta "$DELTA" --perms "$PERMS" \
      --episodes --winperms 300 --seed 0 \
      --outdir "out_periods_grid_${cond}_relorigin"
    echo "[RUN] 03 pseudo pairs (relative origin)  cond=$cond"
    NLOPM_PERIOD_ORIGIN=relative python 03_pseudo_pairs_gridnull.py \
      --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
      --periods "$PERIODS" --delta "$DELTA" --perms "$PERMS" \
      --max_pseudo 300 --seed 0 \
      --outdir "out_pseudo_grid_${cond}_relorigin"
  done
fi

# ---- R: 04 LMM, 05 latency, 07 pseudo interaction, 06 sensitivity -------
if [ "$SKIP_R" = "1" ] || ! command -v Rscript >/dev/null 2>&1; then
  cat <<EOF
==============================================================
Python side finished. Run the following in an R environment (in this folder):
  Rscript 04_LMM_condition_by_period_gridnull.R out_periods_grid_visible,out_periods_grid_invisible out_pseudo_grid_visible,out_pseudo_grid_invisible lmm_grid_out
  Rscript 05_dt_latency_by_period_gridnull.R    out_periods_grid_visible,out_periods_grid_invisible dt_grid_out
  Rscript 07_pseudo_interaction_gridnull.R      out_pseudo_grid_visible,out_pseudo_grid_invisible pseudo_interaction_grid
  Rscript 08_R1-2_p3_dyad_influence.R           out_periods_grid_visible,out_periods_grid_invisible dyad_influence_grid P3   # v4.4
  bash    06_sensitivity_analyses_gridnull.sh   # [A] dyad exclusion, [B] boundaries +/-5 s (re-runs 02, so it takes time)
  # [C] version with the origin at the first shared chew (already written when RELATIVE=1):
  Rscript 04_LMM_condition_by_period_gridnull.R out_periods_grid_visible_relorigin,out_periods_grid_invisible_relorigin out_pseudo_grid_visible_relorigin,out_pseudo_grid_invisible_relorigin lmm_grid_out_relorigin
  Rscript 05_dt_latency_by_period_gridnull.R    out_periods_grid_visible_relorigin,out_periods_grid_invisible_relorigin dt_grid_out_relorigin
==============================================================
EOF
  exit 0
fi

echo "[RUN] 04 LMM by period"
Rscript 04_LMM_condition_by_period_gridnull.R \
  out_periods_grid_visible,out_periods_grid_invisible \
  out_pseudo_grid_visible,out_pseudo_grid_invisible lmm_grid_out
echo "[RUN] 05 latency by period"
Rscript 05_dt_latency_by_period_gridnull.R out_periods_grid_visible,out_periods_grid_invisible dt_grid_out
echo "[RUN] 07 pseudo interaction"
Rscript 07_pseudo_interaction_gridnull.R out_pseudo_grid_visible,out_pseudo_grid_invisible pseudo_interaction_grid
if [ -f 08_R1-2_p3_dyad_influence.R ]; then
  echo "[RUN] 08 dyad influence (leave-one-dyad-out, sign distribution; v4.4)"
  Rscript 08_R1-2_p3_dyad_influence.R out_periods_grid_visible,out_periods_grid_invisible dyad_influence_grid P3 || echo "[WARN] 08 failed"
fi
echo "[RUN] 06 sensitivity analyses [A][B]"
PERMS="$PERMS" bash 06_sensitivity_analyses_gridnull.sh
if [ "$RELATIVE" = "1" ]; then
  echo "[RUN] [C] LMM / latency with relative period origin"
  Rscript 04_LMM_condition_by_period_gridnull.R \
    out_periods_grid_visible_relorigin,out_periods_grid_invisible_relorigin \
    out_pseudo_grid_visible_relorigin,out_pseudo_grid_invisible_relorigin lmm_grid_out_relorigin
  Rscript 05_dt_latency_by_period_gridnull.R out_periods_grid_visible_relorigin,out_periods_grid_invisible_relorigin dt_grid_out_relorigin
fi
echo "[DONE] run_R1_2_trialfix.sh"
