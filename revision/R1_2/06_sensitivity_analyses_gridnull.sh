#!/usr/bin/env bash
# =====================================================================
# Step 6 (R1-2, grid-aligned null): Sensitivity analyses
#   [A] Exclude the visible dyad with sparse annotations in the final
#       minute (default: 20170710_07_08, all 3 sessions)
#       -> re-run 04 (LMM + pseudo), 07 (interaction), 05 (latency)
#   [B] Shift period boundaries by +/-5 s (delta = 0.1 only)
#       -> re-run 02 (grid null) and 04
#
# Layout: run from the directory that contains the *_gridnull scripts and
#   out_periods_grid_{visible|invisible} / out_pseudo_grid_{visible|invisible}
#
# Usage:
#   bash 06_sensitivity_analyses_gridnull.sh          # run [A] then [B]
#   bash 06_sensitivity_analyses_gridnull.sh A        # run [A] only
#   bash 06_sensitivity_analyses_gridnull.sh B        # run [B] only
# =====================================================================
set -euo pipefail

MODE="${1:-all}"

# ========= config =========
VIS_DIR="../annotation_raw_data_chewonly/visible_pair_Human"    # [TRIALFIX] data with the Trial columns removed
INVIS_DIR="../annotation_raw_data_chewonly/invisible_pair_Human"  # [TRIALFIX]
PERIODS_DIRS="out_periods_grid_visible,out_periods_grid_invisible"
PSEUDO_DIRS="out_pseudo_grid_visible,out_pseudo_grid_invisible"
EXCLUDE_DYAD="20170710_07_08"
PERMS="${PERMS:-10000}"   # [TRIALFIX] can be overridden via environment variable
# ==========================

if [ "$MODE" = "A" ] || [ "$MODE" = "all" ]; then
echo "===================================================="
echo "[A] Sensitivity: exclude dyad $EXCLUDE_DYAD"
echo "===================================================="
python3 - "$PERIODS_DIRS" "$PSEUDO_DIRS" "$EXCLUDE_DYAD" <<'PY'
import glob, os, sys
import pandas as pd

periods_dirs, pseudo_dirs, dyad = sys.argv[1:4]

for d in periods_dirs.split(","):
    out = d.rstrip("/") + "_excl_dyad"
    os.makedirs(out, exist_ok=True)
    for stem in ["chew_sync_summary_periods_master.csv", "dt_long_periods_master.csv"]:
        p = os.path.join(d, stem)
        if not os.path.exists(p):
            print(f"[WARN] not found, skipped: {p}")
            continue
        df = pd.read_csv(p)
        n0 = len(df)
        df = df[~df["pair_id"].astype(str).str.startswith(dyad)]
        df.to_csv(os.path.join(out, stem), index=False)
        print(f"{d}/{stem}: {n0} -> {len(df)} rows (excluded {n0 - len(df)})")

for d in pseudo_dirs.split(","):
    out = d.rstrip("/") + "_excl_dyad"
    os.makedirs(out, exist_ok=True)
    files = glob.glob(os.path.join(d, "pseudo_pairs_*.csv"))
    if not files:
        print(f"[WARN] no pseudo_pairs_*.csv in {d}")
    for p in files:
        df = pd.read_csv(p)
        n0 = len(df)
        mask = (df["a_id"].astype(str).str.startswith(dyad) |
                df["b_id"].astype(str).str.startswith(dyad))
        df = df[~mask]
        df.to_csv(os.path.join(out, os.path.basename(p)), index=False)
        print(f"{d}/{os.path.basename(p)}: {n0} -> {len(df)} rows")
PY

PV="${PERIODS_DIRS%%,*}_excl_dyad"; PI="${PERIODS_DIRS##*,}_excl_dyad"
SV="${PSEUDO_DIRS%%,*}_excl_dyad";  SI="${PSEUDO_DIRS##*,}_excl_dyad"
Rscript 04_LMM_condition_by_period_gridnull.R "$PV,$PI" "$SV,$SI" lmm_grid_out_excl_dyad
Rscript 07_pseudo_interaction_gridnull.R "$SV,$SI" pseudo_interaction_grid_excl_dyad
Rscript 05_dt_latency_by_period_gridnull.R "$PV,$PI" dt_grid_out_excl_dyad
fi

if [ "$MODE" = "B" ] || [ "$MODE" = "all" ]; then
echo "===================================================="
echo "[B] Sensitivity: period boundaries +/-5 s (delta=0.1, grid null)"
echo "===================================================="
for VAR in minus5 plus5; do
  if [ "$VAR" = "minus5" ]; then
    P="P1=0:15,P2=15:55,P3=55:180"
  else
    P="P1=0:25,P2=25:65,P3=65:180"
  fi
  echo "--- variant $VAR : $P ---"
  python 02_Micro_Analysis_NLOPM_periods_gridnull.py --indir "$VIS_DIR" --tag visible \
    --periods "$P" --delta 0.1 --perms "$PERMS" --episodes --winperms 300 \
    --outdir "out_periods_grid_visible_${VAR}"
  python 02_Micro_Analysis_NLOPM_periods_gridnull.py --indir "$INVIS_DIR" --tag invisible \
    --periods "$P" --delta 0.1 --perms "$PERMS" --episodes --winperms 300 \
    --outdir "out_periods_grid_invisible_${VAR}"
  Rscript 04_LMM_condition_by_period_gridnull.R \
    "out_periods_grid_visible_${VAR},out_periods_grid_invisible_${VAR}" "" "lmm_grid_out_${VAR}"
done
fi

echo "===================================================="
echo "Done. Outputs:"
echo "  [A] lmm_grid_out_excl_dyad/ , pseudo_interaction_grid_excl_dyad/ , dt_grid_out_excl_dyad/"
echo "  [B] lmm_grid_out_minus5/ , lmm_grid_out_plus5/"
echo "===================================================="
