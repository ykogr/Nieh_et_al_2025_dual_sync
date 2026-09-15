#!/usr/bin/env bash
# =====================================================================
# run_gridnull_main_trialfix.sh
#   Redo the whole main analysis of gridnull_main/ with the data cleaned of the
#   Trial-column contamination (../annotation_raw_data_chewonly) and the fixed 03.
#
# Location     : iscience_revision_analysis/gridnull_main/
# Prerequisite : 00_trialfix/run_00_trialfix.sh has already been run
#            (03_Micro_Analysis_NOLPM_gridnull_v2.py contains the PATCH markers)
# Usage    : bash run_gridnull_main_trialfix.sh            # all deltas
#            DELTAS="0.1" bash run_gridnull_main_trialfix.sh   # restrict deltas
#            PERMS=1000 bash ...                          # fewer permutations for a smoke test
#            STAGE=downstream bash run_gridnull_main_trialfix.sh
#              -> do not redo the permutation computation (03); from the existing out_grid_*_d*
#                 rebuild only step 4 (Table 2, changepoints, Fig S2, Fig 4a-d, coverage, dt.csv).
#                 Use e.g. after replacing a plotting script.
#                 Old figures are moved aside to pre_redo_backup_<date>/ (not part of the public export).
#            STAGE=fig4ad bash run_gridnull_main_trialfix.sh
#              -> (v4.3.3) redraw only Fig 4a-d (09_plot_figure4ad_gridnull_v2.py, vis/invis share the
#                 y axis per top/bottom row). The old fig4ad_gridnull/ is moved aside to pre_redo_backup_<date>/.
#                 With FIG4_YLIM_MIN_PAIRS=15 the y range is determined only from time points covered by >= 15 pairs
#                 (default 0 = all time points. The curves themselves are always drawn in full).
#            STAGE=fig2 bash run_gridnull_main_trialfix.sh
#              -> (v4.4.1) draw only manuscript Figure 2 (dt distribution, delta = 0.1/0.2/0.3, vis/invis)
#                 (10_plot_figure2_timelags.py, input dt_trialfix.csv). Output fig2_timelags/timelags.{pdf,png}
#                 -> figs/timelags.pdf on Overleaf. Neither the permutation computation nor step 4 is run.
#            STAGE=table1_ci bash run_gridnull_main_trialfix.sh
#              -> (v4.3.4) compute only the 95% CIs of Table 1 (02_table1_ci.R, Wald intervals using the
#                 Kenward-Roger df/SE = the same quantities as the t and p in the table). Output df_verification_out/table1_ci.{csv,txt}.
#                 Neither the permutation computation nor step 4 is run. Requires R (lme4, lmerTest, pbkrtest).
#            STAGE=bayes_check bash run_gridnull_main_trialfix.sh
#              -> (v4.2) only verify the output of 06_bayes_gridnull_fit.R (stan_fits/summary_*.csv).
#                 Neither the permutation computation nor step 4 is run.
#                 (i)  check_beta_intervals_v2.py  (list intervals by bin number)      -> beta_intervals_check.txt
#                 (ii) beta_intervals_seconds.py   (convert to seconds and sum per period) -> beta_intervals_seconds_d<delta>.txt
#                 (iii) if Rscript is available, 07_plot_figure4e_gridnull.R (Fig 4e, x axis in video seconds) -> fig4e_gridnull/
#                 The "... s" values in the letter and the Fig 4e caption are taken from this output.
#
# What it does
#   1. Move the old outputs (out_grid_*_d*, ../Micro_Analysis_LMM/*.csv,
#      Micro_Analysis_Bayesian/*.csv) aside to pre_trialfix_backup_<date>/
#   2. Run 03 for each delta x condition (same arguments as the original run)
#   3. Rebuild
#        ../Micro_Analysis_LMM/<delta>_chew_sync_summary_<cond>_pair_Human.csv
#      read by 04_Micro_Analysis_LMM.Rmd / 01_R2-8 / 02_R2-9, and
#        Micro_Analysis_Bayesian/{vis,invis}_<delta>.csv   (t_center + z_S_t__* columns only)
#      read by 06_bayes_gridnull_fit.R
#   4. Run 05 (Table 2), 08 (changepoints, K=2/3), 08_plot (Fig S2), 09 (Fig 4a-d),
#      make_dt_csv (dt.csv for R2-11)
#   5. Print the R steps (if Rscript is available, only 01_R2-8 is run automatically)
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

# ---- kit version check (v4.2.1): if an old 00_trialfix is left in place, the new patches are not applied ----
NEED_KIT="${NEED_KIT:-4.2}"
KIT_VER=$(cat ../00_trialfix/KIT_VERSION 2>/dev/null || echo "0")
if [ "$(printf '%s\n%s\n' "$NEED_KIT" "$KIT_VER" | sort -V | head -1)" != "$NEED_KIT" ]; then
  echo "[ERROR] ../00_trialfix is kit v$KIT_VER, but this kick script requires v$NEED_KIT or later."
  echo "        Re-extract trialfix_kit.zip (latest) and run  cp -R trialfix_kit/. .  in iscience_revision_analysis/"
  echo "        (check: cat 00_trialfix/KIT_VERSION must be $NEED_KIT or later)"
  exit 1
fi
echo "[INFO] trialfix_kit v$KIT_VER (00_trialfix/KIT_VERSION)"

RAW="${RAW:-../annotation_raw_data_chewonly}"
PERMS="${PERMS:-10000}"
DELTAS="${DELTAS:-0.1 0.2 0.3 0.5 1.0}"
LMM_DIR="${LMM_DIR:-../Micro_Analysis_LMM}"
BAYES_DIR="Micro_Analysis_Bayesian"
STAGE="${STAGE:-all}"          # all | downstream | fig4ad | fig2 | table1_ci | bayes_check
FIG4_SHAREY="${FIG4_SHAREY:-row}"              # v4.3.3: row (a/b and c/d share the y axis) | none (independent, as in the old version)
FIG4_YLIM_MIN_PAIRS="${FIG4_YLIM_MIN_PAIRS:-0}" # v4.3.3: minimum number of covering pairs used to compute the y range (0 = all time points)
# Script names work with or without the _v2 suffix (whichever exists is used)
pick() { for f in "$@"; do [ -f "$f" ] && { echo "$f"; return; }; done; echo "$1"; }
MAIN=$(pick 03_Micro_Analysis_NOLPM_gridnull_v2.py 03_Micro_Analysis_NOLPM_gridnull.py)
CP=$(pick 08_changepoint_gridnull_v2.py 08_changepoint_gridnull.py)
F4=$(pick 09_plot_figure4ad_gridnull_v2.py 09_plot_figure4ad_gridnull.py)
echo "[INFO] scripts: $MAIN / $CP / $F4"
STAMP=$(date +%Y%m%d_%H%M)
BK="pre_trialfix_backup_${STAMP}"

# ---- checks -----------------------------------------------------------
[ -f "$MAIN" ] || { echo "[ERROR] 03_Micro_Analysis_NOLPM_gridnull(_v2).py not found"; exit 1; }
[ -d "$RAW" ] || { echo "[ERROR] $RAW does not exist. Run 00_trialfix/run_00_trialfix.sh first"; exit 1; }
grep -q "EVENT-COLUMN PATCH" "$MAIN" || { echo "[ERROR] $MAIN is not patched. Run 00_trialfix/01_patch_scripts.py"; exit 1; }
grep -q "TIME-AXIS PATCH"    "$MAIN" || { echo "[ERROR] $MAIN lacks the time-axis patch"; exit 1; }
mkdir -p "$LMM_DIR" "$BAYES_DIR"

# ---- STAGE=fig4ad (v4.3.3): redraw only Fig 4a-d and exit ----------------------------
if [ "$STAGE" = "fig4ad" ]; then
  TL="out_grid_visible_d0.1/timelines_all_master.csv,out_grid_invisible_d0.1/timelines_all_master.csv"
  [ -f out_grid_visible_d0.1/timelines_all_master.csv ] || { echo "[ERROR] out_grid_visible_d0.1/timelines_all_master.csv does not exist. Run the main analysis with STAGE=all first"; exit 1; }
  grep -q -- "--sharey" "$F4" || { echo "[ERROR] $F4 is the old version (v2). Overwrite it with gridnull_main/09_plot_figure4ad_gridnull_v2.py from trialfix_kit v4.3.3"; exit 1; }
  BK="pre_redo_backup_${STAMP}"; mkdir -p "$BK"
  [ -d fig4ad_gridnull ] && mv fig4ad_gridnull "$BK/" && echo "[INFO] old fig4ad_gridnull/ moved aside to $BK/"
  echo "[RUN] $F4 (Fig 4a-d; --sharey $FIG4_SHAREY --ylim_min_pairs $FIG4_YLIM_MIN_PAIRS)"
  python "$F4" --timelines "$TL" --outdir fig4ad_gridnull --sharey "$FIG4_SHAREY" --ylim_min_pairs "$FIG4_YLIM_MIN_PAIRS"
  echo "[DONE] STAGE=fig4ad  -> fig4ad_gridnull/figure4ad_gridnull.{pdf,png}, figure4ad_ylim.txt (y range), figure4ad_peaks.txt"
  echo "       Upload: fig4ad_gridnull/figure4ad_gridnull.png, figure4ad_ylim.txt"
  exit 0
fi

# ---- STAGE=fig2 (v4.4.1): redraw only Fig 2 (dt distribution) and exit ----------------------------
if [ "$STAGE" = "fig2" ]; then
  [ -f dt_trialfix.csv ] || { echo "[ERROR] dt_trialfix.csv does not exist. Run make_dt_csv first with STAGE=all or STAGE=downstream"; exit 1; }
  [ -f 10_plot_figure2_timelags.py ] || { echo "[ERROR] 10_plot_figure2_timelags.py does not exist. Copy it from gridnull_main/ of trialfix_kit v4.4.1"; exit 1; }
  BK="pre_redo_backup_${STAMP}"; mkdir -p "$BK"
  [ -d fig2_timelags ] && mv fig2_timelags "$BK/" && echo "[INFO] old fig2_timelags/ moved aside to $BK/"
  echo "[RUN] 10_plot_figure2_timelags.py (Fig 2; dt_trialfix.csv)"
  python 10_plot_figure2_timelags.py --dtcsv dt_trialfix.csv --outdir fig2_timelags
  echo "[DONE] STAGE=fig2  -> fig2_timelags/timelags.{pdf,png}, timelags_source_data.csv"
  echo "       Overleaf: overwrite figs/timelags.pdf with fig2_timelags/timelags.pdf. Upload: timelags.png, timelags_source_data.csv"
  exit 0
fi

# ---- STAGE=table1_ci (v4.3.4): compute only the 95% CIs of Table 1 and exit --------------------
if [ "$STAGE" = "table1_ci" ]; then
  [ -f 02_table1_ci.R ] || { echo "[ERROR] 02_table1_ci.R does not exist. Copy gridnull_main/02_table1_ci.R from trialfix_kit v4.3.4 into this folder"; exit 1; }
  command -v Rscript >/dev/null 2>&1 || { echo "[ERROR] Rscript not found. Run in an R environment"; exit 1; }
  ls "$LMM_DIR"/0.1_chew_sync_summary_visible_pair_Human.csv >/dev/null 2>&1 || { echo "[ERROR] no summary CSV in $LMM_DIR. Run STAGE=all first"; exit 1; }
  echo "[RUN] Rscript 02_table1_ci.R $LMM_DIR df_verification_out"
  Rscript 02_table1_ci.R "$LMM_DIR" df_verification_out
  echo "[DONE] STAGE=table1_ci -> df_verification_out/table1_ci.csv, table1_ci.txt"
  echo "       Upload: df_verification_out/table1_ci.csv, table1_ci.txt (source of the CI column of manuscript Table 1)"
  exit 0
fi

# ---- STAGE=bayes_check: only verify the Stan output and exit ------------------------
if [ "$STAGE" = "bayes_check" ]; then
  [ -d stan_fits ] || { echo "[ERROR] stan_fits/ does not exist. Run Rscript 06_bayes_gridnull_fit.R fit ... first"; exit 1; }
  shopt -s nullglob
  SUMS=(stan_fits/summary_z_S_t_delta*_human.csv)
  shopt -u nullglob
  [ ${#SUMS[@]} -gt 0 ] || { echo "[ERROR] stan_fits/summary_z_S_t_delta*_human.csv does not exist"; exit 1; }
  echo "[RUN] check_beta_intervals_v2.py -> beta_intervals_check.txt  (${#SUMS[@]} summary files)"
  { echo "# generated $(date +%Y-%m-%d\ %H:%M) by run_gridnull_main_trialfix.sh STAGE=bayes_check"
    echo "# python check_beta_intervals_v2.py ${SUMS[*]}"
    python check_beta_intervals_v2.py "${SUMS[@]}"; } | tee beta_intervals_check.txt
  for f in "${SUMS[@]}"; do
    tag=$(basename "$f" | sed -E 's/summary_z_S_t_delta([0-9]+)_human\.csv/\1/')   # 01, 02, 03, 05, 1
    case "$tag" in 01) d=0.1;; 02) d=0.2;; 03) d=0.3;; 05) d=0.5;; 1) d=1.0;; *) echo "[WARN] unknown delta tag $tag ($f) -> skip"; continue;; esac
    if [ -f "$BAYES_DIR/vis_${d}.csv" ] && [ -f "$BAYES_DIR/invis_${d}.csv" ]; then
      echo "[RUN] beta_intervals_seconds.py delta=$d -> beta_intervals_seconds_d${d}.txt"
      python beta_intervals_seconds.py --summary "$f" --vis "$BAYES_DIR/vis_${d}.csv" \
        --invis "$BAYES_DIR/invis_${d}.csv" --out "beta_intervals_seconds_d${d}.txt" > /dev/null
      tail -n +1 "beta_intervals_seconds_d${d}.txt" | grep -E "total|per period|bins " || true
      if command -v Rscript >/dev/null 2>&1; then
        echo "[RUN] Rscript 07_plot_figure4e_gridnull.R $tag $d fig4e_gridnull"
        Rscript 07_plot_figure4e_gridnull.R "$tag" "$d" fig4e_gridnull || echo "[WARN] Fig 4e plot failed (delta $d)"
      fi
    else
      echo "[WARN] $BAYES_DIR/vis_${d}.csv / invis_${d}.csv missing, skipping the conversion to seconds (delta $d)"
    fi
  done
  command -v Rscript >/dev/null 2>&1 || echo "[INFO] Rscript not available, so Fig 4e was not drawn. In an R environment: Rscript 07_plot_figure4e_gridnull.R 01 0.1 fig4e_gridnull"
  echo "[DONE] STAGE=bayes_check  -> upload: beta_intervals_check.txt, beta_intervals_seconds_d*.txt (, fig4e_gridnull/*.png)"
  exit 0
fi

if [ "$STAGE" = "all" ]; then
# ---- 1. backup ----------------------------------------------------------
mkdir -p "$BK/Micro_Analysis_LMM" "$BK/Micro_Analysis_Bayesian"
shopt -s nullglob
for d in out_grid_*_d*; do [ -d "$d" ] && mv "$d" "$BK/"; done
for f in "$LMM_DIR"/*_chew_sync_summary_*_pair_Human.csv; do cp -p "$f" "$BK/Micro_Analysis_LMM/"; done
for f in "$BAYES_DIR"/vis_*.csv "$BAYES_DIR"/invis_*.csv; do cp -p "$f" "$BK/Micro_Analysis_Bayesian/"; done
for d in out_changepoint_grid* fig_S2_gridnull fig4ad_gridnull coverage_grid fig4e_gridnull fig2_timelags; do [ -d "$d" ] && mv "$d" "$BK/"; done
for f in beta_intervals_check.txt beta_intervals_seconds_d*.txt; do [ -f "$f" ] && mv "$f" "$BK/"; done
# Also copy the R-side outputs (overwritten by manual runs): Bayesian stan_fits/, LMM knit results
[ -d stan_fits ] && cp -rp stan_fits "$BK/stan_fits"
for f in "$LMM_DIR"/*.html "$LMM_DIR"/*.nb.html; do cp -p "$f" "$BK/Micro_Analysis_LMM/"; done
shopt -u nullglob
echo "[INFO] old outputs moved aside to $BK/"

# ---- 2-3. main runs -----------------------------------------------------
for d in $DELTAS; do
  for cond in visible invisible; do
    OUT="out_grid_${cond}_d${d}"
    echo "=============================================================="
    echo "[RUN] delta=$d cond=$cond perms=$PERMS -> $OUT"
    python "$MAIN" \
      --indir "$RAW/${cond}_pair_Human" --tag "${cond}_pair_Human" \
      --delta "$d" --t 180 --perms "$PERMS" \
      --win 5 --step 1 --episode_z 1.64 --episode_minmatch 3 --episode_minsec 3.0 \
      --seed 0 --outdir "$OUT"
    # LMM / R2-8 / R2-9 input
    cp "$OUT/${d}_chew_sync_summary_${cond}_pair_Human.csv" "$LMM_DIR/"
    # Bayesian input: t_center + z_S_t__<pair> only (same shape as before)
    short=vis; [ "$cond" = "invisible" ] && short=invis
    python - "$OUT/timelines_all_master.csv" "$BAYES_DIR/${short}_${d}.csv" <<'PY'
import sys, pandas as pd
src, dst = sys.argv[1:3]
df = pd.read_csv(src)
keep = ["t_center"] + [c for c in df.columns if c.startswith("z_S_t__")]
df[keep].sort_values("t_center").to_csv(dst, index=False)
print(f"  Bayesian input -> {dst}: {len(keep)-1} series x {len(df)} bins")
PY
  done
done

else
# ---- downstream only: move aside only the old figures/tables and go to step 4 ----------------------------
[ -d out_grid_visible_d0.1 ] || { echo "[ERROR] out_grid_visible_d0.1 does not exist. Run the main analysis with STAGE=all first"; exit 1; }
BK="pre_redo_backup_${STAMP}"; mkdir -p "$BK"
for d in out_changepoint_grid out_changepoint_grid_k3 fig_S2_gridnull fig4ad_gridnull coverage_grid fig2_timelags; do [ -d "$d" ] && mv "$d" "$BK/"; done
for f in dt_trialfix.csv dt_trialfix_per_pair.csv; do [ -f "$f" ] && mv "$f" "$BK/"; done
echo "[INFO] STAGE=downstream: the permutation computation is not rerun; old figures/tables are moved aside to $BK/ and only step 4 is rebuilt"
fi

# ---- 4. downstream (python) --------------------------------------------
TL="out_grid_visible_d0.1/timelines_all_master.csv,out_grid_invisible_d0.1/timelines_all_master.csv"
echo "=============================================================="
echo "[RUN] 05_prevalence_table2 (Table 2)"
python 05_prevalence_table2.py --indir "$LMM_DIR" || echo "[WARN] run 05 only after the summaries for all deltas (0.1,0.2,0.3,0.5,1.0) are available"
if [ -f out_grid_visible_d0.1/timelines_all_master.csv ]; then
  echo "[RUN] 08_changepoint (K=2, K=3)"
  python "$CP" --timelines "$TL" --n_cp 2 --outdir out_changepoint_grid
  python "$CP" --timelines "$TL" --n_cp 3 --outdir out_changepoint_grid_k3
  # Pick the newest version (v2 > v1). Duplicate copies with "(1)" are ignored
  S2=$(ls 08_plot_figure_S2_gridnull*.py 2>/dev/null | grep -v ' (' | sort -V | tail -1)
  echo "[RUN] $S2 (Fig S2)"
  python "$S2" --timelines "$TL" --labels visible,invisible \
    --group_curve out_changepoint_grid/changepoint_group_curve.csv \
    --boundaries 20 60 --delta 0.1 --outdir fig_S2_gridnull || echo "[WARN] Fig S2 plot failed (check --group_curve file name)"
  echo "[RUN] 09_plot_figure4ad (Fig 4a-d; --sharey $FIG4_SHAREY --ylim_min_pairs $FIG4_YLIM_MIN_PAIRS)"
  if grep -q -- "--sharey" "$F4"; then
    python "$F4" --timelines "$TL" --outdir fig4ad_gridnull --sharey "$FIG4_SHAREY" --ylim_min_pairs "$FIG4_YLIM_MIN_PAIRS"
  else
    echo "[WARN] $F4 is the old version (no shared y axis). Overwrite it with 09 from kit v4.3.3"
    python "$F4" --timelines "$TL" --outdir fig4ad_gridnull
  fi
  echo "[RUN] coverage_table (number of covering pairs per bin; basis for deciding where to start drawing)"
  python coverage_table.py --timelines "$TL" --labels visible,invisible --min_pairs 15 --outdir coverage_grid
fi
echo "[RUN] make_dt_csv (R2-11 / Fig 2 latency table)"
python make_dt_csv.py --raw "$RAW" --script "$MAIN" --deltas 0.1,0.2,0.3 --out dt_trialfix.csv
if [ -f 10_plot_figure2_timelags.py ]; then
  echo "[RUN] 10_plot_figure2_timelags.py (Fig 2, dt distribution -> fig2_timelags/timelags.pdf)"
  python 10_plot_figure2_timelags.py --dtcsv dt_trialfix.csv --outdir fig2_timelags || echo "[WARN] Fig 2 plot failed"
fi

# ---- 5. R steps -----------------------------------------------------------
cat <<EOF
==============================================================
Python side finished. Next, in R (cmdstanr environment):
  (a) LMM (Table 1, main-text z_S comparison):  knit 04_Micro_Analysis_LMM.Rmd
      (input $LMM_DIR/*_chew_sync_summary_*.csv has already been overwritten)
  (b) df verification (R2-8):   Rscript 01_R2-8_df_verification.R $LMM_DIR df_verification_out
      or ../R2_8_2_9_2_11/run_R2_8_trialfix.sh
  (b') 95% CIs of Table 1:  Rscript 02_table1_ci.R $LMM_DIR df_verification_out
      or STAGE=table1_ci bash run_gridnull_main_trialfix.sh
  (c) Bayesian state space (Fig 4e, β(t)):  Rscript 06_bayes_gridnull_fit.R fit $(echo $DELTAS | tr ' ' ',')
      then  STAGE=bayes_check bash run_gridnull_main_trialfix.sh
        -> beta_intervals_check.txt (bin numbers), beta_intervals_seconds_d*.txt (seconds, totals per period),
           fig4e_gridnull/fig4e_beta_delta01.{pdf,png} (Fig 4e, x axis = video seconds)
      Note: t_center is video time (same x.5 s grid as the old version). bin i is the i-th t_center in
            ascending order among those appearing in vis/invis_<delta>.csv (bin_id = as.integer(factor(t_center)) in 06).
==============================================================
EOF
if command -v Rscript >/dev/null 2>&1; then
  echo "[RUN] Rscript 01_R2-8_df_verification.R"
  Rscript 01_R2-8_df_verification.R "$LMM_DIR" df_verification_out || echo "[WARN] 01_R2-8 failed"
  if [ -f 02_table1_ci.R ]; then
    echo "[RUN] Rscript 02_table1_ci.R (95% CIs of Table 1; v4.3.4)"
    Rscript 02_table1_ci.R "$LMM_DIR" df_verification_out || echo "[WARN] 02_table1_ci failed"
  fi
fi
echo "[DONE] run_gridnull_main_trialfix.sh"
