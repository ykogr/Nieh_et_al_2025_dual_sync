#!/usr/bin/env bash
# =====================================================================
# run_R2_2_trialfix.sh  (v3: integrates the δ=0.2 surrogates and the LMM step)
#   Re-run R2-2/2-5/2-10/3-2/3-3 (null-distribution symmetry check, jitter/ISI alternative nulls,
#   reactive shift test, band test, condition x period LMM of surrogate z_S) in one go
#   with the data cleaned of the Trial-column contamination.
#
# Location     : iscience_revision_analysis/R2_2_2_5_2_10_3_2_3_3/
# Prerequisite : 00_trialfix/run_00_trialfix.sh has been run
#            (the 01_/02_/03_ *.py in this folder contain the PATCH markers)
# Usage    : bash run_R2_2_trialfix.sh
#            PERMS=1000 bash run_R2_2_trialfix.sh   # smoke test
#            SKIP_R=1   bash run_R2_2_trialfix.sh   # Python part only (R manually, with the commands printed at the end)
#
# Script names are chosen automatically from what is in the folder:
#   *_gridnull.py / *_gridnull.R if present, otherwise the unsuffixed ones.
# Arguments are the same as in the Usage at the top of each script (perms=10000, seed=0, P1/P2/P3).
#
# Steps and outputs (all go to out_*/ directly under this folder):
#   01 null symmetry        δ=0.1                      -> out_symmetry/
#   02 jitter/ISI/circ      δ=0.1 (w=0.35,0.2,0.5)     -> out_jitter/jitter_null_0.1_<cond>.csv
#                           δ=0.2 (w=0.35,0.5)         -> out_jitter/jitter_null_0.2_<cond>.csv
#                                                         out_jitter/jitter_null_master.csv (both δ)
#                                                         out_jitter/jitter_null_group_summary.txt
#   03 shift / band test    δ=0.1, 0.2                 -> out_shift/
#   02b LMM input           jitter035 / jitter020 / jitter050 / isi -> out_jitter/lmm_in_<config>/
#   04  LMM (R)             same 4 configs             -> out_jitter/lmm_out_<config>/lmm_results.txt etc.
# Existing out_symmetry* / out_jitter* / out_shift* are moved aside to pre_trialfix_backup_<timestamp>/.
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

RAW="${RAW:-../annotation_raw_data_chewonly}"
PERMS="${PERMS:-10000}"
PERIODS="${PERIODS:-P1=0:20,P2=20:60,P3=60:180}"
SKIP_R="${SKIP_R:-0}"
STAMP=$(date +%Y%m%d_%H%M)
BK="pre_trialfix_backup_${STAMP}"

pick() {  # pick <stem>  -> prints "<stem>_gridnull.py" if present else "<stem>.py"
  if [ -f "$1_gridnull.py" ]; then echo "$1_gridnull.py"; else echo "$1.py"; fi
}
S01=$(pick 01_null_symmetry_check)
S02=$(pick 02_jitter_surrogates)
S03=$(pick 03_reactive_shift_test)
S02B="02b_make_lmm_input.py"
if   [ -f 04_LMM_condition_by_period_gridnull.R ]; then S04=04_LMM_condition_by_period_gridnull.R
elif [ -f 04_LMM_condition_by_period.R ];          then S04=04_LMM_condition_by_period.R
else S04=""; fi

[ -d "$RAW" ] || { echo "[ERROR] $RAW does not exist (run 00_trialfix first)"; exit 1; }
for s in "$S01" "$S02" "$S03"; do
  grep -q "EVENT-COLUMN PATCH" "$s" || { echo "[ERROR] $s is not patched (run 00_trialfix/01_patch_scripts.py on this folder)"; exit 1; }
done
echo "[INFO] using: $S01 / $S02 / $S03 / $S02B / ${S04:-(no R)}"

mkdir -p "$BK"
shopt -s nullglob
for d in out_symmetry* out_jitter* out_shift*; do [ -d "$d" ] && mv "$d" "$BK/"; done
shopt -u nullglob
echo "[INFO] old outputs moved aside to $BK/"

# ---- 01: equivalence of the A-shift / B-shift nulls (R2-10) ----
for cond in visible invisible; do
  echo "=============================================================="
  echo "[RUN] 01 null symmetry  cond=$cond"
  python "$S01" --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
    --periods "$PERIODS" --delta 0.1 --perms "$PERMS" --seed 0 --outdir out_symmetry
done

# ---- 02: rate-preserving surrogates (R2-2 / R2-5) ----
#   δ=0.1: w = 0.35 (main setting), 0.2, 0.5.  δ=0.2: w = 0.35, 0.5 (w ≤ δ is omitted as it carries no information)
for cond in visible invisible; do
  echo "[RUN] 02 jitter / ISI surrogates  cond=$cond delta=0.1"
  python "$S02" --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
    --periods "$PERIODS" --delta 0.1 --w "0.35,0.2,0.5" --nulls "jitter,isi,circ" \
    --perms "$PERMS" --seed 0 --outdir out_jitter
done
for cond in visible invisible; do
  echo "[RUN] 02 jitter / ISI surrogates  cond=$cond delta=0.2"
  python "$S02" --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
    --periods "$PERIODS" --delta 0.2 --w "0.35,0.5" --nulls "jitter,isi,circ" \
    --perms "$PERMS" --seed 0 --outdir out_jitter
done

# ---- 03: reactive realignment / band test (R2-2) ----
for d in 0.1 0.2; do
  for cond in visible invisible; do
    echo "[RUN] 03 reactive shift / band test  cond=$cond delta=$d"
    python "$S03" --indir "$RAW/${cond}_pair_Human" --tag "$cond" \
      --periods "$PERIODS" --delta "$d" --taus "0.15,0.2,0.25" \
      --perms "$PERMS" --seed 0 --outdir out_shift
  done
done

# ---- 02b: convert surrogate z_S to LMM input (all δ present in the master are included) ----
if [ -f "$S02B" ]; then
  echo "[RUN] 02b make LMM input (4 configs)"
  python "$S02B" --master out_jitter/jitter_null_master.csv --null_type jitter --w 0.35 --outdir out_jitter/lmm_in_jitter035
  python "$S02B" --master out_jitter/jitter_null_master.csv --null_type jitter --w 0.2  --outdir out_jitter/lmm_in_jitter020
  python "$S02B" --master out_jitter/jitter_null_master.csv --null_type jitter --w 0.5  --outdir out_jitter/lmm_in_jitter050
  python "$S02B" --master out_jitter/jitter_null_master.csv --null_type isi            --outdir out_jitter/lmm_in_isi
else
  echo "[WARN] $S02B not found, skipping 02b/04"
fi

# ---- 04: condition x period LMM (R) ----
if [ -f "$S02B" ]; then
  if [ "$SKIP_R" = "1" ] || [ -z "$S04" ]; then
    echo "[INFO] R step skipped. To run manually (in this folder):"
    for k in jitter035 jitter020 jitter050 isi; do
      echo "  Rscript ${S04:-04_LMM_condition_by_period_gridnull.R} out_jitter/lmm_in_$k \"\" out_jitter/lmm_out_$k"
    done
  else
    for k in jitter035 jitter020 jitter050 isi; do
      echo "[RUN] 04 LMM  $k"
      Rscript "$S04" "out_jitter/lmm_in_$k" "" "out_jitter/lmm_out_$k"
    done
  fi
fi

echo "[DONE] run_R2_2_trialfix.sh"
echo "  out_symmetry/null_symmetry_summary.txt"
echo "  out_jitter/jitter_null_group_summary.txt   (δ=0.1, 0.2)"
echo "  out_shift/shift_test_group_summary.txt"
echo "  out_jitter/lmm_out_<jitter035|jitter020|jitter050|isi>/lmm_results.txt"
