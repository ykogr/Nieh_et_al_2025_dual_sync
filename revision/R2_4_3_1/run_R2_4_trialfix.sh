#!/usr/bin/env bash
# =====================================================================
# run_R2_4_trialfix.sh
#   Re-run R2-4 / 3-1 (grid null vs continuous null comparison, jitter simulation,
#   zS reproduction by an independent implementation, inter-rater agreement) after the fix.
#
# Location     : iscience_revision_analysis/R2_4_3_1/
# Prerequisite : 00_trialfix/run_00_trialfix.sh has been run
#            (PATCH markers in R2-4_zS_*.py, R2-4_jitter_simulation.py, 03_..._gridnull.py)
# Usage    : bash run_R2_4_trialfix.sh
#            NPERM=500 bash run_R2_4_trialfix.sh   # smoke test
#            STAGE=reliability bash run_R2_4_trialfix.sh
#              -> do not re-run the null comparison (nullcheck) and the jitter simulation;
#                 rebuild only rater reliability (zS reproduction inter/intra, agreement rate) and R3-1 power.
#                 For re-running after the zS-reproduction null was changed to grid-aligned in v4.2.
#            STAGE=icc_ci bash run_R2_4_trialfix.sh   (PER_VIDEO_DIR=... may specify the location of per_video.csv)
#              -> (v4.3) recompute nothing; from the existing R2-4_zS_{inter,intra}_d<delta>_per_video.csv
#                 write the ICC(A,1) of session-level zS and its 95% confidence interval (bootstrap, F method for reference)
#                 to R2-4_zS_icc_ci.txt / .csv. Finishes in seconds. Also runs automatically at the end of all / reliability.
#
# Notes
#  - R2-4_*.py originally read only the Chew column, so there was no Trial-row contamination, but
#    the comparison had no tolerance, so the handling of "events exactly delta apart" fluctuated
#    with the floating-point representation. The patch fixes that, so the numbers change slightly.
#  - Rater reliability comes in two parts:
#      inter (independent second coder)      IRR_ROOT   = IRR_final          (ICC_annotate layout)
#      intra (re-annotation by the original coder 8 months later) INTRA_ROOT = ../k_ICC/ICC_result (repository layout)
#    In both cases R2-4_zS_reproduction.py / R2-4_irr_agreement.py read only the Chew column, so
#    Trial-row contamination has no effect and the data are passed as is, without preprocessing. Outputs are separated by _inter / _intra.
#    If the location differs:  INTRA_ROOT=path/to/ICC_result bash run_R2_4_trialfix.sh
#  - The variance components (resid_sd, dyad_sd) and the observed P3 contrast (estimate and SE) for R3-1_power_sensitivity.R
#    are read automatically from ../R1_2/lmm_grid_out/lmm_results.txt (post-fix z_S, delta=0.1 model),
#    recorded in R3-1_power_params.txt and then passed to R. No manual entry.
#    (v4.2) The "observed P3 contrast" note at the end of R3-1 was hard-coded, so
#    01_patch_scripts.py rewrites it to be passed via the --observed/--observed_se arguments.
#  - (v4.2) The null distribution of the zS reproduction (R2-4_zS_reproduction.py) is changed to the same
#    grid-aligned shift (--null grid) as in the General note of the letter, and run for all of delta = 0.1/0.2/0.3/0.5/1.0.
#    Output names are R2-4_zS_{inter,intra}_d<delta>_summary.txt.
#    These patches are applied by calling 01_patch_scripts.py at the start of this kick script (safe to run repeatedly).
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
NPERM="${NPERM:-2000}"           # default of R2-4_zS_nullcheck
JPERM="${JPERM:-1000}"           # default of R2-4_jitter_simulation
IRR_ROOT="${IRR_ROOT:-IRR_final}"                 # inter-rater: second coder (ICC_annotate layout)
# intra-rater: re-annotation by the original coder (repository layout). The first candidate found is used
if [ -z "${INTRA_ROOT:-}" ]; then
  for c in ../k_ICC/ICC_result ICC_result ../ICC_result; do [ -d "$c" ] && { INTRA_ROOT="$c"; break; }; done
fi
INTRA_ROOT="${INTRA_ROOT:-}"
DELTAS="${DELTAS:-0.1 0.2 0.3 0.5 1.0}"   # deltas for which the zS reproduction (inter/intra) is run
REPRO_PERM="${REPRO_PERM:-10000}"          # default of R2-4_zS_reproduction
STAGE="${STAGE:-all}"                      # all | reliability | icc_ci
ICC_BOOT="${ICC_BOOT:-10000}"               # (v4.3) number of bootstrap replicates for R2-4_zS_icc_ci.py
STAMP=$(date +%Y%m%d_%H%M)
BK="pre_trialfix_backup_${STAMP}"

if [ "$STAGE" = "icc_ci" ]; then
  # (v4.3.1) confidence intervals only: existing outputs are not moved aside; only R2-4_zS_icc_ci.* is moved aside before overwriting.
  # Location of per_video.csv: specify with PER_VIDEO_DIR=... . If unset, search this folder, then pre_redo_backup_* / pre_trialfix_backup_* (newest first).
  PV_DIR="${PER_VIDEO_DIR:-}"
  if [ -z "$PV_DIR" ]; then
    if ls R2-4_zS_*_d*_per_video.csv >/dev/null 2>&1; then
      PV_DIR="."
    else
      for c in $(ls -d pre_redo_backup_* pre_trialfix_backup_* 2>/dev/null | sort -r); do
        if ls "$c"/R2-4_zS_*_d*_per_video.csv >/dev/null 2>&1; then PV_DIR="$c"; break; fi
      done
    fi
  fi
  if [ -z "$PV_DIR" ] || ! ls "$PV_DIR"/R2-4_zS_*_d*_per_video.csv >/dev/null 2>&1; then
    echo "[ERROR] R2-4_zS_*_d*_per_video.csv not found."
    echo "  working folder: $(pwd)"
    if [ ! -f R2-4_zS_reproduction.py ]; then
      echo "  -> R2-4_zS_reproduction.py is not in this folder. Are you running inside trialfix_kit/R2_4_3_1 (the kit side)?"
      echo "     After  cp -R trialfix_kit/. iscience_revision_analysis/ , run it in iscience_revision_analysis/R2_4_3_1."
      exit 1
    fi
    echo "  Please paste the following diagnostics:"
    echo "  --- R2-4_zS_* in this folder ---"; ls -la R2-4_zS_* 2>/dev/null || echo "  (none)"
    echo "  --- backup folders ---"; ls -d pre_redo_backup_* pre_trialfix_backup_* 2>/dev/null || echo "  (none)"
    for c in $(ls -d pre_redo_backup_* pre_trialfix_backup_* 2>/dev/null); do echo "  [$c]"; ls "$c" 2>/dev/null | grep 'R2-4_zS_' | sed 's/^/     /'; done
    echo "  --- does R2-4_zS_reproduction.py write per_video.csv ---"
    grep -c "per_video" R2-4_zS_reproduction.py 2>/dev/null | sed 's/^/     occurrences of per_video: /' || echo "     R2-4_zS_reproduction.py does not exist"
    echo "  --- candidates in subfolders ---"; find . -maxdepth 3 -name 'R2-4_zS_*per_video.csv' 2>/dev/null | head -20
    echo "  Fix: (a) if found in a folder listed above:  PER_VIDEO_DIR=<that folder> STAGE=icc_ci bash run_R2_4_trialfix.sh"
    echo "        (b) if found nowhere:  STAGE=reliability bash run_R2_4_trialfix.sh  (re-runs the zS reproduction for all deltas; rebuilds both per_video.csv and summary.txt)"
    exit 1
  fi
  [ "$PV_DIR" != "." ] && echo "[INFO] reading per_video.csv from $PV_DIR/ instead of this folder (make sure it is the output of the same run as summary.txt)"
  [ -f R2-4_zS_icc_ci.txt ] && { mkdir -p "pre_redo_backup_${STAMP}"; mv R2-4_zS_icc_ci.txt R2-4_zS_icc_ci.csv "pre_redo_backup_${STAMP}/" 2>/dev/null || true; }
  echo "[RUN] R2-4_zS_icc_ci (ICC(A,1) 95% CI, bootstrap ${ICC_BOOT}; input: $PV_DIR)"
  python R2-4_zS_icc_ci.py --pattern "$PV_DIR/R2-4_zS_*_d*_per_video.csv" --n_boot "$ICC_BOOT" --seed 20260827 --out_prefix R2-4_zS_icc_ci
  echo "[DONE] run_R2_4_trialfix.sh (STAGE=icc_ci)"
  exit 0
fi

[ -d "$RAW" ] || { echo "[ERROR] $RAW does not exist (run 00_trialfix first)"; exit 1; }
# v4.2: apply the additional patches (GRID-NULL / OBSERVED-ARG) to this folder (skipped if already done)
python ../00_trialfix/01_patch_scripts.py . || { echo "[ERROR] 01_patch_scripts.py returned MANUAL. Check the output above"; exit 1; }
for s in R2-4_zS_reproduction.py R2-4_zS_nullcheck.py R2-4_jitter_simulation.py; do
  grep -q "FLOAT-TOL PATCH" "$s" || { echo "[ERROR] $s is not patched"; exit 1; }
done
grep -q "GRID-NULL PATCH"    R2-4_zS_reproduction.py  || { echo "[ERROR] R2-4_zS_reproduction.py lacks the GRID-NULL patch"; exit 1; }
grep -q "OBSERVED-ARG PATCH" R3-1_power_sensitivity.R || { echo "[ERROR] R3-1_power_sensitivity.R lacks the OBSERVED-ARG patch"; exit 1; }

mkdir -p "$BK"
shopt -s nullglob
if [ "$STAGE" = "all" ]; then
  TO_BK="R2-4_nullcheck_* R2-4_jitter_sim_* R2-4_zS_* R2-4_irr_results* R2-4_nullcheck_lmm* R3-1_power*"
else
  rmdir "$BK" 2>/dev/null || true
  BK="pre_redo_backup_${STAMP}"; mkdir -p "$BK"
  TO_BK="R2-4_zS_* R2-4_irr_results* R3-1_power*"
fi
for f in $TO_BK; do
  [ -f "$f" ] && case "$f" in *.py|*.R|*.bak) ;; *) mv "$f" "$BK/";; esac
done
shopt -u nullglob
echo "[INFO] old outputs moved aside to $BK/ (STAGE=$STAGE)"

if [ "$STAGE" = "all" ]; then
echo "=============================================================="
echo "[RUN] R2-4_zS_nullcheck (grid vs continuous null, all deltas)"
python R2-4_zS_nullcheck.py --annotated_data "$RAW" --deltas 0.1,0.2,0.3,0.5,1.0 \
  --n_perm "$NPERM" --periods "P1:0-20,P2:20-60,P3:60-180" --seed 20260827 --out_prefix R2-4_nullcheck

echo "[RUN] R2-4_jitter_simulation (delta 0.1)"
python R2-4_jitter_simulation.py --annotated_data "$RAW" --jitter_sd 0.24 --delta 0.1 \
  --reps 10 --n_perm "$JPERM" --periods "P1:0-20,P2:20-60,P3:60-180" --seed 20260827 --out_prefix R2-4_jitter_sim
else
echo "[INFO] STAGE=reliability: nullcheck / jitter are not re-run (existing R2-4_nullcheck_*, R2-4_jitter_sim_* are used as is)"
fi

# ---- rater reliability (inter: second coder / intra: re-annotation by the same coder) ---------------------
# zS reproduction: null is grid-aligned (--null grid, same as the General note), outputs separated per delta
if [ -d "$IRR_ROOT" ]; then
  echo "[RUN] inter-rater: R2-4_zS_reproduction (null=grid, deltas: $DELTAS) / R2-4_irr_agreement  (coder B = $IRR_ROOT)"
  for d in $DELTAS; do
    python R2-4_zS_reproduction.py --annotated_data "$RAW" --irr_root "$IRR_ROOT" \
      --delta "$d" --null grid --n_perm "$REPRO_PERM" --seed 20260827 --out_prefix "R2-4_zS_inter_d${d}"
  done
  python R2-4_irr_agreement.py   --annotated_data "$RAW" --irr_root "$IRR_ROOT" --out_prefix R2-4_irr_results_inter
else
  echo "[SKIP] inter-rater: $IRR_ROOT does not exist (can be specified with IRR_ROOT=...)"
fi
if [ -n "$INTRA_ROOT" ] && [ -d "$INTRA_ROOT" ]; then
  echo "[RUN] intra-rater: R2-4_zS_reproduction (null=grid, deltas: $DELTAS) / R2-4_irr_agreement  (coder B = $INTRA_ROOT)"
  for d in $DELTAS; do
    python R2-4_zS_reproduction.py --annotated_data "$RAW" --rerun_root "$INTRA_ROOT" \
      --delta "$d" --null grid --n_perm "$REPRO_PERM" --seed 20260827 --out_prefix "R2-4_zS_intra_d${d}"
  done
  python R2-4_irr_agreement.py   --annotated_data "$RAW" --rerun_root "$INTRA_ROOT" --out_prefix R2-4_irr_results_intra
else
  echo "[SKIP] intra-rater: re-annotation folder (ICC_result) not found (specify with INTRA_ROOT=...)"
fi
# (v4.3) attach 95% confidence intervals to the ICC(A,1) of session-level zS (only reads per_video.csv)
if ls R2-4_zS_*_d*_per_video.csv >/dev/null 2>&1; then
  echo "[RUN] R2-4_zS_icc_ci (ICC(A,1) 95% CI, bootstrap ${ICC_BOOT})"
  python R2-4_zS_icc_ci.py --n_boot "$ICC_BOOT" --seed 20260827 --out_prefix R2-4_zS_icc_ci
else
  echo "[SKIP] R2-4_zS_icc_ci: per_video.csv does not exist"
fi

# ---- R3-1: read the variance components from the new LMM output of R1_2 ----------------------------
LMM_TXT="${LMM_TXT:-../R1_2/lmm_grid_out/lmm_results.txt}"
if [ -f "$LMM_TXT" ]; then
  python - "$LMM_TXT" R3-1_power_params.txt <<'PY'
import re, sys, datetime
src, dst = sys.argv[1:3]
txt = open(src, encoding="utf-8").read()
# look only at the z_S, delta = 0.1 block
i = txt.find("DV = z_S | delta = 0.1")
blk = txt[i:] if i >= 0 else txt
blk = blk[:blk.find("DV = ", 10)] if blk.find("DV = ", 10) > 0 else blk
resid = re.search(r"Residual\s+[\d.eE+-]+\s+([\d.eE+-]+)", blk).group(1)
dyad  = re.search(r"dyad_id\s+\(Intercept\)\s+[\d.eE+-]+\s+([\d.eE+-]+)", blk).group(1)
def fe(name):
    m = re.search(r"^" + re.escape(name) + r"\s+([\-\d.eE+]+)", blk, re.M)
    return float(m.group(1))
contrast = fe("conditionvisible") + fe("conditionvisible:periodP3")   # visible - invisible in P3
# SE of the P3 contrast: emmeans block "period = P3:" -> "invisible - visible  est  SE ..."
m_se = re.search(r"period = P3:.*?invisible - visible\s+[\-\d.]+\s+([\d.]+)", blk, re.S)
se = m_se.group(1) if m_se else "NA"
grid = sorted(set([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0] + [round(contrast,3)]))
with open(dst, "w", encoding="utf-8") as f:
    f.write(f"# generated {datetime.datetime.now():%Y-%m-%d %H:%M} from {src} (DV = z_S, delta = 0.1)\n")
    f.write(f"resid_sd={resid}\ndyad_sd={dyad}\nobserved_P3_contrast={contrast:.5f}\nobserved_P3_se={se}\n")
    f.write("grid=" + ",".join(f"{g:g}" for g in grid) + "\n")
print(open(dst).read())
PY
  RESID_SD=$(grep '^resid_sd=' R3-1_power_params.txt | cut -d= -f2)
  DYAD_SD=$(grep '^dyad_sd='  R3-1_power_params.txt | cut -d= -f2)
  GRID=$(grep '^grid='        R3-1_power_params.txt | cut -d= -f2)
  OBS=$(grep '^observed_P3_contrast=' R3-1_power_params.txt | cut -d= -f2)
  OBS_SE=$(grep '^observed_P3_se='    R3-1_power_params.txt | cut -d= -f2)
else
  echo "[WARN] $LMM_TXT does not exist. Finish step 2 (R1_2) first. Skipping R3-1"
fi

cat <<EOF
==============================================================
R side (in this folder; run automatically below if Rscript is available):
  Rscript R2-4_nullcheck_lmm.R --per_video R2-4_nullcheck_per_video.csv --out_prefix R2-4_nullcheck_lmm
  Rscript R3-1_power_sensitivity.R --out_prefix R3-1_power --nsim 1000 --seed 20260823 \
      --resid_sd ${RESID_SD:-<not obtained>} --dyad_sd ${DYAD_SD:-<not obtained>} --grid ${GRID:-<not obtained>} \
      --observed ${OBS:-<not obtained>} --observed_se ${OBS_SE:-NA}
    -> the values are recorded in R3-1_power_params.txt (auto-generated from ../R1_2/lmm_grid_out/lmm_results.txt).
       The old values (0.91338 / 0.08451) are based on the Trial-contaminated data and are not used.
==============================================================
EOF
if command -v Rscript >/dev/null 2>&1; then
  if [ -f R2-4_nullcheck_per_video.csv ]; then
    Rscript R2-4_nullcheck_lmm.R --per_video R2-4_nullcheck_per_video.csv --out_prefix R2-4_nullcheck_lmm || echo "[WARN] nullcheck_lmm failed"
  fi
  if [ -n "${RESID_SD:-}" ]; then
    Rscript R3-1_power_sensitivity.R --out_prefix R3-1_power --nsim 1000 --seed 20260823 \
      --resid_sd "$RESID_SD" --dyad_sd "$DYAD_SD" --grid "$GRID" \
      --observed "$OBS" --observed_se "$OBS_SE" || echo "[WARN] R3-1 power failed"
  fi
fi
echo "[DONE] run_R2_4_trialfix.sh"
