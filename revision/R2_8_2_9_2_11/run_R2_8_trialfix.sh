#!/usr/bin/env bash
# =====================================================================
# run_R2_8_trialfix.sh
#   Redo R2-8 (degrees-of-freedom verification), R2-9 (multiple-comparison correction),
#   R2-11 (dip test of the latency distribution) using the re-run results of gridnull_main.
#
# Location     : iscience_revision_analysis/R2_8_2_9_2_11/
# Prerequisite : gridnull_main/run_gridnull_main_trialfix.sh has been run
#            -> ../Micro_Analysis_LMM/*_chew_sync_summary_*.csv are new
#            -> ../gridnull_main/dt_trialfix.csv (output of make_dt_csv.py) exists
# Usage    : bash run_R2_8_trialfix.sh
#
# The scripts in this folder do not read the raw data themselves, so they need no fix.
# Only the inputs are swapped. dt.csv is the one rebuilt by make_dt_csv.py
# (dt = t_A - t_B of the NLOPM matches over the whole interval, columns 0.1_vis / 0.1_invis ...).
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

LMM_DIR="${LMM_DIR:-../Micro_Analysis_LMM}"
DTCSV="${DTCSV:-../gridnull_main/dt_trialfix.csv}"
STAMP=$(date +%Y%m%d_%H%M)
BK="pre_trialfix_backup_${STAMP}"

[ -d "$LMM_DIR" ] || { echo "[ERROR] $LMM_DIR does not exist"; exit 1; }
[ -f "$DTCSV" ]   || { echo "[ERROR] $DTCSV does not exist (run gridnull_main/make_dt_csv.py first)"; exit 1; }
# Check that the summary is the new one (the old summary has the Trial count in n_A and T_overlap of 178.x)
python - "$LMM_DIR" <<'PY'
import glob, os, sys, pandas as pd
d = sys.argv[1]
fs = sorted(glob.glob(os.path.join(d, "0.1_chew_sync_summary_visible_pair_Human.csv")))
if not fs: sys.exit("[ERROR] 0.1_chew_sync_summary_visible_pair_Human.csv does not exist")
s = pd.read_csv(fs[0])
if s["T_overlap"].max() > 178.0:
    sys.exit("[ERROR] the summary appears to be the old version (T_overlap>178: start point still at Trial=1). Re-run gridnull_main")
print(f"[OK] {fs[0]}: T_overlap {s.T_overlap.min():.1f}-{s.T_overlap.max():.1f} s, n_A {s.n_A.min()}-{s.n_A.max()}")
PY

mkdir -p "$BK"
shopt -s nullglob
for d in multiplicity_out dip_test_out dip_figure_out df_verification_out; do [ -d "$d" ] && mv "$d" "$BK/"; done
shopt -u nullglob

echo "[RUN] 02_R2-9 multiplicity correction"
python 02_R2-9_multiplicity_correction.py --indir "$LMM_DIR" --outdir multiplicity_out
echo "[RUN] 03_R2-11 dip test"
python 03_R2-11_dip_test.py --dtcsv "$DTCSV" --outdir dip_test_out
echo "[RUN] 04_R2-11 dip figure"
JS=dip_test_out/dip_test_jitter_summary.csv
if [ -f "$JS" ]; then
  python 04_R2-11_dip_figure.py --dtcsv "$DTCSV" --jitter-summary "$JS" --outdir dip_figure_out
else
  echo "[WARN] jitter summary csv not found in dip_test_out, so make the figure manually: python 04_R2-11_dip_figure.py --dtcsv $DTCSV --jitter-summary <csv>"
fi

if command -v Rscript >/dev/null 2>&1; then
  echo "[RUN] 01_R2-8 df verification"
  Rscript 01_R2-8_df_verification.R "$LMM_DIR" df_verification_out
else
  echo "In an R environment: Rscript 01_R2-8_df_verification.R $LMM_DIR df_verification_out"
fi
echo "[DONE] run_R2_8_trialfix.sh"
