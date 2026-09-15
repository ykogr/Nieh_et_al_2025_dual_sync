# Revision analyses (iScience, 2026) — scripts, data, and raw logs

This folder contains everything used for the revision of Nieh et al.: the corrected annotation data,
the analysis scripts, and the raw output logs that every number in the revised manuscript and the
response letter was taken from. Nothing here overwrites the original (pre-revision) pipeline in the
repository root.

## How the folders map to the reviewer comments

| Folder | Reviewer comments | Entry point | Main outputs (raw logs) |
|---|---|---|---|
| `00_trialfix/` | General note (data correction) | `run_00_trialfix.sh` | `legacy_header_report.txt`, `chewonly_report.txt`, `verify_report.txt` |
| `annotation_raw_data/` | — | — | original event annotations (input; one solo file had its header renamed to the standard `Trial, Chew_k`, see `00_trialfix/legacy_header_backup/`) |
| `annotation_raw_data_chewonly/` | General note | written by `00_trialfix` | corrected annotations (input to all analyses) |
| `gridnull_main/` | main pipeline; Table 1, Figures 2, 4, S2; R2-8 | `run_gridnull_main_trialfix.sh` | `out_grid_*_d*/`, `df_verification_out/table1_ci.csv`, `coverage_grid/`, `out_changepoint_grid*/`, `Micro_Analysis_Bayesian/`, `beta_intervals_seconds_d*.txt`, `fig*_gridnull/`, `fig2_timelags/` |
| `Micro_Analysis_LMM/` | Table S4 | written by `gridnull_main` | `*_chew_sync_summary_*.csv`, `prevalence_table2.csv` |
| `R1_2/` | R1-2 (period-resolved analyses, Table S1, Figure 5) | `run_R1_2_trialfix.sh` | `out_periods_grid_*/`, `lmm_grid_out*/`, `dt_grid_out*/`, `out_pseudo_grid_*/`, `pseudo_interaction_grid*/`, `dyad_influence_grid/` |
| `R2_2_2_5_2_10_3_2_3_3/` | R2-2, R2-5, R2-10, R3-2, R3-3 (Tables S2, S3) | `run_R2_2_trialfix.sh` | `out_symmetry/`, `out_jitter/`, `out_shift/` |
| `R2_4_3_1/` | R2-4, R3-1 (inter-/intra-rater reliability, power) | `run_R2_4_trialfix.sh` | `R2-4_irr_results_*`, `R2-4_zS_*`, `R2-4_nullcheck_*`, `R2-4_jitter_sim_*`, `R3-1_power_*`; inputs: `IRR_final/` (second coder), `w1/`–`w15/` |
| `k_ICC/` | R2-4 (intra-rater re-annotation) | `ICC_annotate.py` | `ICC_result/` |
| `R2_6_3_5/` | R2-6, R3-5 (macro-scale state-space model) | `r26_condiff_kit/` (see its README) | `sensitivity_out/`, `xscale_out/`, `div_check_*.csv` |
| `R2_8_2_9_2_11/` | R2-8, R2-9, R2-11 | `run_R2_8_trialfix.sh` | `df_verification_out/`, `multiplicity_out/`, `dip_test_out/`, `dip_figure_out/` |

`README_procedure.md` documents the order in which the folders were run and how the corrected data were
verified. The original annotation sheet handed to the second coder (`R2_4_3_1/R2-4_annotation_cheatsheet_final.txt`,
Japanese) is kept verbatim; an English translation is provided next to it (`*_en.txt`).

## Not included

- Stan fit objects (`stan_fits*/`) and other files above GitHub's 100 MB limit.
- Working copies and superseded outputs (`pre_*backup*/`, `upload_*/`, `old/`).

`TREE.txt` is the complete file tree of this folder.
