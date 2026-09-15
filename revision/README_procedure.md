# Procedure for the Trial-column contamination fix kit (trialfix_kit)

Created: 2026-09-07 (v2: 2026-09-08 — unified the time axis and the period origin to video time; added the relative-origin sensitivity analysis [C] and the coverage table. v2.1: keep the 22 chewing times found in header-less columns; support script names without `_v2`; state explicitly how the counts in the reports are defined. v3: merged all R2_2 steps into a single script; added 99_tidy, gitignore and 98_export. v4: added Micro_Analysis_LMM / R2_6_3_5 / k_ICC to 98_export; fixed the Fig S2 script selection in step 1 and added `STAGE=downstream`. v4.1: the R3-1 variance components in step 4 are now taken automatically from the R1_2 LMM output; coder reliability split into two runs, inter / intra. v4.2: the null distribution for the zS reproduction in step 4 changed to the grid-aligned `--null grid` and run for all δ; the hard-coded R3-1 note turned into a command-line argument (both are patches applied by 01_patch_scripts.py); added `STAGE=reliability`. Added `STAGE=bayes_check` to step 6 (converts β(t) intervals to seconds; 07_plot_figure4e_gridnull.R draws Fig 4e in video seconds). v4.2.1: added `00_trialfix/KIT_VERSION`; the kick scripts of steps 1 and 4 now stop if they detect an outdated 00_trialfix. v4.3: added `R2-4_zS_icc_ci.py` and `STAGE=icc_ci` to step 4 (attaches a bootstrap 95% confidence interval to the ICC(A,1) of session-level z_S; the source of the "ICC 95% CI" in letter R2-4). Added a note on the Figure S1/S2 numbering mismatch on the manuscript side (Section 4). v4.3.1: when `STAGE=icc_ci` cannot find per_video.csv in this folder, it searches `pre_redo_backup_*` / `pre_trialfix_backup_*` (newest first); if still not found it prints diagnostics (file list, script version) and stops. The location can be given explicitly with `PER_VIDEO_DIR=`. v4.3.2: detects when it is run inside the kit-side folder `trialfix_kit/R2_4_3_1` and gives guidance. v4.3.3: the Fig 4a–d plotting script `gridnull_main/09_plot_figure4ad_gridnull_v2.py` updated to v2.1 (top row a/b and bottom row c/d each share the y axis between vis and invis; numerical output unchanged). Added `STAGE=fig4ad` (redraw Fig 4a–d only) to step 1. v4.3.4: added `gridnull_main/02_table1_ci.R` and `STAGE=table1_ci` to step 1 (the source of the 95% CI column of manuscript Table 1; Wald intervals based on the Kenward-Roger df/SE, i.e. the same quantities as the t and p in the table). v4.4: added `R1_2/08_R1-2_p3_dyad_influence.R` and `STAGE=influence` to step 2 (leave-one-dyad-out for the P3 condition difference, sign distribution of dyad-level z_S with Wilcoxon tests, figure; the source of "the effect is a small shift spread across many pairs rather than driven by a few dyads"). For 02_table1_ci.R and 08 we only verified that they run without error in our R environment (only the PI-side run results are used for the numbers). v4.4.1: added `gridnull_main/10_plot_figure2_timelags.py` and `STAGE=fig2` to step 1 (manuscript Figure 2 = bar chart of the Δt distribution, δ = 0.1/0.2/0.3, vis/invis. Input is `dt_trialfix.csv` from make_dt_csv. Output `fig2_timelags/timelags.pdf` → Overleaf `figs/timelags.pdf`). `STAGE=all`/`downstream` also draw it automatically right after make_dt_csv)

**How to apply each update (always the same)**: copy the entire contents of the unpacked `trialfix_kit/` over `iscience_revision_analysis/`.
```bash
cp -R trialfix_kit/. iscience_revision_analysis/
```
The kit contains only the run scripts (`run_*.sh`, `00_trialfix/`, helper .py files, README, 98/99); it does not contain the analysis scripts themselves or any output, so copying everything over breaks nothing.
Scope: the whole micro-level analysis (NLOPM synchrony analysis) under `iscience_revision_analysis/`

---

## 1. What had happened (summary)

### 1-1. Trial-column contamination
The annotation CSV files hold one trial (Trial) per row, in wide format: `Trial, Chew_1, Chew_2, ...`.
The loader `load_times_from_csv`, inherited from the published `03_Micro_Analysis_NOLPM.py`,
**converts every cell of the table to a number and flattens them into a single column** (`df.stack()` or `to_numpy().ravel()`) when there is no column named `t`,
so the values 1, 2, ..., k in the `Trial` column were mixed in as "bite events at t = 1 s, 2 s, ..., k s".

This produced the following three distortions.

| Distortion | What it is | Effect |
|---|---|---|
| Spurious events | 2–41 spurious events per file (median 14) at t = 1..k s | n_A, n_B inflated |
| Spurious perfect synchrony | A and B have nearly the same number of Trials, so **simultaneous** spurious events for both appear at t = 1..k s and always match at Δt = 0 | S and z_S biased upward. In particular the "synchrony right after onset" in P1 (0–20 s) is almost entirely this artefact |
| Fixed start point | Start of the overlap interval = max(first of A, first of B) was always 1.0 s (Trial=1) | The actual first common bite is at 3.2–19.5 s. The time axis was shifted for every pair by about 1 s plus an individual offset |

**A second, incidental finding (v2.1)**: in two B files of the invisible condition (`20170707_01_02_02_B.csv`, `20170707_07_08_01_B.csv`),
the chewing in the last trial exceeds the number of `Chew_k` columns, and the remaining times (12 values at 165.8–179.5 s and 10 values at 173.8–179.5 s)
sit in header-less columns (`Unnamed: 15` etc. in pandas). The old loader read every cell, so it picked these up too.
The v2 "keep only Chew columns" step dropped these 22 values, so in v2.1 header-less columns that contain values are renamed to `Chew_15, …`
and kept. This can be checked in the line "files with value-bearing unnamed columns kept" in `chewonly_report.txt`
(no files other than these two are affected).

### 1-2. Floating-point comparison (secondary, but affects the numbers)
Comparisons such as `abs(a - b) <= delta` in `NLOPM_match` had no tolerance, so events **exactly δ apart** on the 0.1 s grid
matched or failed to match depending on the floating-point representation (e.g. 0.8 − 0.7 = 0.10000000000000009 > 0.1).
Counting the total number of matches at δ = 0.1 (60 pairs, whole interval) on the data fixed to Chew columns only gives
3,669 without tolerance → 5,127 with tolerance, and **the count changes in all 60 pairs**.
Because the null distribution goes through the same comparison, the effect on z_S is smaller than on S, but it still has to be fixed.

### 1-3. Time axis and period origin (policy agreed with the co-authors)
The window centres of the windowed analysis were computed as `start of each pair's overlap interval + 2.5 + k` seconds (the original design was relative to the start point).
Later it was decided to "align on video time", and the fact that the Trial values 1, 2, 3, … fall exactly on the integer-second grid was used
to effectively align on video time (3.5, 4.5, …) — that is the history. However, the Trial values did not merely act as grid markers;
they entered both A and B as spurious events at t = 1..k s, always matched at Δt = 0 and pushed z_S up.

The fixed version achieves the same alignment without the spurious events.

- **Time series (windowed analysis in 03)**: t_center stays in video time. The edge of each pair's first window is aligned to "the start point rounded up to an integer second"
  (for a pair starting at 4.3 s the windows start at 5–10 s, t_center = 7.5, 8.5, …). All pairs fall on the same x.5 s grid as in the old version,
  and the t_center values equal those of the old results = the same video time. Only the first < 1 s of each pair is lost.
- **Period-wise analyses (02/03 in R1_2, the period-based analyses in R2_2)**: P1 = 0–20 s etc. are defined in **video time** and intersected with each pair's overlap interval
  (`seg_start = max(t0, start)`, `seg_end = min(t1, end)`). The original code had `seg_start = start + t0` (relative to the start point), but
  since the start was always 1.0 s in the old data, this was effectively "video time 1–21 s". In video time, pairs that start late (up to 19.5 s)
  have almost no P1 data and get `sufficient = 0` (report that count).
- **Sensitivity analysis [C]**: the R1_2 kick also produces automatically a relative-time version whose origin is "the first common bite of both partners"
  (`out_*_relorigin`, `lmm_grid_out_relorigin`). Its purpose is to show in both coordinate systems that "whether measured in session time or in time since
  the start of eating, there is no above-chance synchrony right after onset (P1)".
- Switchable via environment variables: `NLOPM_TIME_AXIS=relative` (03), `NLOPM_PERIOD_ORIGIN=relative` (02/03/R2_2).
- The time of the first common bite does not differ between conditions (visible 5.9 s, invisible 6.6 s, Mann–Whitney p = 0.78),
  so the condition comparison is not distorted under either origin.

---

## 2. Kit contents and where they go

The kit mirrors the existing folder structure. Unpack the zip inside `iscience_revision_analysis/` (or
copy the contents of each folder into the corresponding folder). **The only existing files with the same name that get overwritten are
`R1_2/06_sensitivity_analyses_gridnull.sh` (data paths changed to `_chewonly`) and, from v4.3.3,
`gridnull_main/09_plot_figure4ad_gridnull_v2.py` (the shared-y-axis version of Fig 4a–d; the old version is identical to the one we sent on 2026-09-02)**.

```
iscience_revision_analysis/
├── 00_trialfix/                      ← new folder (place it as a whole)
│   ├── run_00_trialfix.sh            one-shot run (copy data → patch code → verify)
│   ├── 00_fix_legacy_header.py       (v5) renames the legacy header of one solo file to Trial/Chew_k (data rows untouched; original kept in legacy_header_backup/)
│   ├── 00_make_chewonly_data.py      annotation_raw_data → annotation_raw_data_chewonly (Chew columns only)
│   ├── 01_patch_scripts.py           patches the .py files in each folder automatically (originals saved as *.pre_trialfix.bak)
│   └── 02_verify_fix.py              checks for missed patches, data and numbers
├── annotation_raw_data/              (do not touch)
├── annotation_raw_data_chewonly/     ← generated automatically in step 0
├── gridnull_main/
│   ├── run_gridnull_main_trialfix.sh ← added
│   ├── make_dt_csv.py                ← added (rebuilds dt.csv for R2-11)
│   ├── coverage_table.py             ← added (number of covered pairs per bin; basis for deciding where to start drawing the figures)
│   ├── beta_intervals_seconds.py     ← added v4.2 (converts the β(t) intervals from bin numbers to video seconds and sums them per period)
│   ├── 07_plot_figure4e_gridnull.R   ← added v4.2 (Fig 4e. The plotting chunk of the public Rmd as is, only the x axis in video seconds)
│   ├── 02_table1_ci.R                ← added v4.3.4 (95% CI for Table 1. Same LMM as 01_R2-8, Kenward-Roger df/SE)
│   ├── 10_plot_figure2_timelags.py   ← added v4.4.1 (manuscript Figure 2. Bar chart of the Δt distribution from dt_trialfix.csv, δ = 0.1/0.2/0.3, vis/invis)
│   └── 09_plot_figure4ad_gridnull_v2.py ← replaced v4.3.3 (Fig 4a–d. v2.1: vis/invis share the y axis per row. Overwrites the existing file of the same name)
├── R1_2/
│   ├── run_R1_2_trialfix.sh          ← added
│   ├── 08_R1-2_p3_dyad_influence.R   ← added v4.4 (leave-one-dyad-out for the P3 condition difference, sign distribution of dyad-level z_S, figure)
│   └── 06_sensitivity_analyses_gridnull.sh ← replaced (VIS_DIR/INVIS_DIR changed to _chewonly)
├── R2_2_2_5_2_10_3_2_3_3/
│   └── run_R2_2_trialfix.sh          ← added
├── R2_4_3_1/
│   ├── run_R2_4_trialfix.sh          ← added
│   └── R2-4_zS_icc_ci.py             ← added v4.3 (95% CI of ICC(A,1) from R2-4_zS_*_per_video.csv. Called at the end of the kick / STAGE=icc_ci)
├── R2_8_2_9_2_11/
│   └── run_R2_8_trialfix.sh          ← added
├── Micro_Analysis_LMM/               (the kick replaces the summary csv files inside; old versions are moved aside into each backup)
├── k_ICC/, R2_6_3_5/                 (not affected, because they do not read bite times from the raw data. If k_ICC
│                                      contains a home-made script that reads Chew columns, add the folder to
│                                      01_patch_scripts.py and check)
```

### What 01_patch_scripts.py rewrites (all target files are detected automatically)

| Patch | Marker | Target function | Target files (as found) |
|---|---|---|---|
| Read Chew columns only | `[EVENT-COLUMN PATCH` | `load_times_from_csv` | gridnull_main/03_…_gridnull_v2.py, R1_2/02_…periods(_gridnull).py, R2_2…/01,02,03(_gridnull).py, R2_4_3_1/03_…_gridnull.py |
| Tolerance 1e-6 s in comparisons | `[FLOAT-TOL PATCH` | `NLOPM_match` (5 places), `nlopm_matches` (independent R2-4 implementation, 5 places) | the above + R2_4_3_1/R2-4_zS_reproduction.py, R2-4_zS_nullcheck.py, R2-4_jitter_simulation.py |
| Time axis on a common grid (video time) | `[TIME-AXIS PATCH` | `sliding_window_with_dt` | gridnull_main/03_…_v2.py, R2_4_3_1/03_…_gridnull.py |
| Period origin in video time | `[PERIOD-ORIGIN PATCH` | the `seg_start = start + t0` lines | R1_2/02_…periods(_gridnull).py, R1_2/03_pseudo_pairs(_gridnull).py, R2_2…/01,02,03(_gridnull).py |
| (v4.2) Grid-aligned null | `[GRID-NULL PATCH` | `zs_period`; adds the argument `--null {continuous,grid}` (default is continuous as before; the kick passes grid) | R2_4_3_1/R2-4_zS_reproduction.py only |
| (v4.2) Hard-coded note turned into an argument | `[OBSERVED-ARG PATCH` | the trailing `cat("Note: observed P3 contrast … = +0.xxx …")` → prints the values of `--observed/--observed_se` | R2_4_3_1/R3-1_power_sensitivity.R only (the only place where a .R file is handled) |

- Files already patched (marker present) are skipped, so it is safe to run any number of times.
- 01_patch_scripts.py prints its own version (`trialfix_kit v4.2.1`) on the first line when it runs. If this line does not appear, or
  the kick stops with "GRID-NULL patch missing" after `[DONE] files changed: 0`, then `00_trialfix/` is still the old version
  (check that `cat 00_trialfix/KIT_VERSION` is 4.2 or higher and redo `cp -R trialfix_kit/. .` from the latest zip).
- Files written in an unexpected way are reported as `MANUAL` and left untouched (exit code 2). Please contact us in that case.
- Scripts that do not read the raw data (05_prevalence, 08_changepoint, 08/09 plot, check_beta_intervals, all R scripts, all of R2_8) are not modified. Only their inputs change.
- The `R2_2_…` zip contained only the plain versions (01_null_symmetry_check.py etc.), but the `_gridnull` versions used for the letter have the same loader and the same `NLOPM_match`, so everything in the folder is patched automatically. The kick prefers the `_gridnull` version when present.

---

## 3. Run order

Assumed environment: conda env `nieh_iscience_2026` (pandas, numpy, scipy, ruptures, matplotlib), R (lme4/lmerTest, cmdstanr).
The number of permutations is 10,000 as in the original runs (it can be reduced for smoke tests with an environment variable such as `PERMS=1000`).

### step 0: fix data and code (tens of seconds)

Step -1 (added in kit v5, 2026-09-15): `annotation_raw_data/solo_Human/20170706_07_03_B.csv` was exported by an older version of the annotation tool
and was the only file whose header was not `Trial, Chew_1, Chew_2, ...`. `00_fix_legacy_header.py` rewrites that header line to the standard names;
the 11 data rows are copied byte for byte (MD5 recorded in `legacy_header_report.txt`) and the original file is kept in `00_trialfix/legacy_header_backup/`.
All counts in `chewonly_report.csv` for this file (174 events, 163 chews, 11 Trial values removed) are unchanged; only the header and the
"dropped named columns" line of `chewonly_report.txt` change. Running the step again does nothing once the header is standard.

How to read the counts in the two reports:
- "total spurious events removed" in `chewonly_report.txt` is the total number of Trial values removed from **all 177 files** (pair 120 + solo 57).
- Item (3) of `verify_report.txt` covers **only the 120 pair files**. "spurious … removed, pair files only" is the count with the same definition,
  and equals the total of `chewonly_report.txt` minus the solo share. "events inside the overlap window" is a different quantity
  (the count restricted to the overlap interval; the difference also includes real chews lost because the start moved later), so it is normal that they do not agree.
```bash
cd iscience_revision_analysis/00_trialfix
bash run_00_trialfix.sh --dry_run    # see what would change
bash run_00_trialfix.sh              # run → check that chewonly_report.txt and verify_report.txt say "ALL OK"
```

### step 1: main analysis (gridnull_main)

The kick automatically distinguishes script names with and without `_v2` (`03_Micro_Analysis_NOLPM_gridnull.py` vs `…_gridnull_v2.py` etc.).
At start it prints the file names it will use as `[INFO] scripts: …`; please check them.  Note: 10 runs (5 δ × 2 conditions) take about 30–60 minutes
```bash
cd ../gridnull_main
bash run_gridnull_main_trialfix.sh
```
- Old outputs are moved aside into `pre_trialfix_backup_<datetime>/`.
- Rebuilds `../Micro_Analysis_LMM/<δ>_chew_sync_summary_<cond>_pair_Human.csv` and
  `Micro_Analysis_Bayesian/{vis,invis}_<δ>.csv`.
- Then runs 05 (Table 2), 08 (change points K=2/3), 08_plot (Fig S2), 09 (Fig 4a–d), coverage_table (coverage table), make_dt_csv (dt_trialfix.csv).
- Check "the first bin with 15 or more pairs" in `coverage_grid/coverage_summary.txt` and decide whether to start Fig 4/S2 there or
  to show the coverage counts alongside the figure (at δ = 0.1, 15 or more pairs from 8.5 s, all 30 pairs from 22.5 s).
- Finally prints the R steps (knit 04_Micro_Analysis_LMM.Rmd, run 06_bayes_gridnull_fit.R).
- For Fig S2 the newest version of `08_plot_figure_S2_gridnull*.py` (v2 > v1) is used. Duplicate copies with "(1)" in the name are ignored.

**Rebuild figures and tables only (do not redo the permutations)**  Note: a few minutes
```bash
STAGE=downstream bash run_gridnull_main_trialfix.sh
```
Using the existing `out_grid_*_d*` as input, re-runs only 05 (Table 2), 08 (change points), 08_plot (Fig S2), 09 (Fig 4a–d), coverage_table and make_dt_csv.
Old figures and tables are moved aside into `pre_redo_backup_<datetime>/` (excluded from the public export). Use this when a plotting script has been replaced.

**Draw Fig 2 (Δt distribution) only (v4.4.1)**  Note: a few seconds
```bash
STAGE=fig2 bash run_gridnull_main_trialfix.sh
```
- Input is `dt_trialfix.csv` (created by make_dt_csv in `STAGE=all` / `downstream`; the same Δt as in the R2-11 dip test). Stops if it is missing.
- Output `fig2_timelags/timelags.pdf` (→ overwrite `figs/timelags.pdf` in Overleaf), `timelags.png`, `timelags_source_data.csv` (count and proportion per bar).
- The figure layout is the same as Figure 2 of the 2026-08-23 manuscript (panels (a)(b)(c) = δ 0.1/0.2/0.3, x axis from +δ to −δ, vis orange, invis blue, legend on the right).
- Caution: the figure in `R2_8_2_9_2_11/dip_figure_out/` is **Figure S1** (Δt after jitter and the dip test), not Figure 2.

**Redraw Fig 4a–d only (v4.3.3)**  Note: a few seconds
```bash
STAGE=fig4ad bash run_gridnull_main_trialfix.sh
```
- `09_plot_figure4ad_gridnull_v2.py` (v2.1) shares the z_S(t) y axis between the top-row panels (a) vis and (b) invis, and the S_excess(t) y axis between the bottom-row panels (c) and (d)
  (`--sharey row`, default). The shared range is "the union of the 95% CI bands of both conditions + 5% margin". To keep them independent as in the old version, use `FIG4_SHAREY=none`.
- Right after onset (before 8.5 s at δ = 0.1) few pairs are covered and the CI band is very wide, so the y range is pulled by it if left alone.
  With `FIG4_YLIM_MIN_PAIRS=15 STAGE=fig4ad bash run_gridnull_main_trialfix.sh` the y range is determined only from time points with 15 or more covered pairs
  (the same criterion as `--min_pairs 15` in coverage_table. The curves themselves are always drawn over all time points; only the initial band that falls outside the range is clipped).
  Decide which one to use by looking at the figure. The setting used is recorded in `fig4ad_gridnull/figure4ad_ylim.txt`.
- Means, CIs and peak values (`figure4ad_peaks.txt`, `figure4ad_source_data.csv`) are identical to v2 (we compared both settings and confirmed no difference).
- Old figures are moved aside into `pre_redo_backup_<datetime>/fig4ad_gridnull/`. `STAGE=all` / `STAGE=downstream` draw with the same settings (environment variables).

**Compute only the 95% CI for Table 1 (v4.3.4)**  Note: a few seconds (requires R)
```bash
STAGE=table1_ci bash run_gridnull_main_trialfix.sh
```
- `02_table1_ci.R` fits the same model as `01_R2-8_df_verification.R` (`dv ~ condition + (1 | pair_id)`) to each DV × δ and
  produces Wald-type 95% intervals (estimate ± t(0.975, df_KR) × SE_KR) using the Kenward-Roger degrees of freedom and standard error.
  Because these are based on the same quantities as the t and p columns of Table 1 (Kenward-Roger), the intervals and p values in the table cannot contradict each other.
  Profile-likelihood intervals are also listed for reference (may be NA in cells where a variance component is 0).
- Output: `df_verification_out/table1_ci.csv`, `table1_ci.txt`. The CI column of Table 1 is taken from this csv.
- Also run automatically at the end of `STAGE=all` (in environments with R).

**Verify the Stan output only (v4.2)**  Note: a few seconds
```bash
STAGE=bayes_check bash run_gridnull_main_trialfix.sh
```
Use after running `06_bayes_gridnull_fit.R`. Details in step 6.

### step 2: period-wise analyses (R1_2)  Note: 02 ×2 + 03 ×2 take 10–20 minutes; sensitivity analysis [B] runs 02 four more times, [C] runs 02/03 twice each
```bash
cd ../R1_2
bash run_R1_2_trialfix.sh        # if Rscript is available, runs 04/05/07 and 06 (sensitivity analyses [A][B]) and [C] (relative origin) automatically
```
- Main results (video-time origin): `out_periods_grid_<cond>`, `out_pseudo_grid_<cond>`, `lmm_grid_out`, `dt_grid_out`, `pseudo_interaction_grid`
- Sensitivity analysis [C] (first common bite as origin): `out_periods_grid_<cond>_relorigin`, `out_pseudo_grid_<cond>_relorigin`, `lmm_grid_out_relorigin`, `dt_grid_out_relorigin`
- `RELATIVE=0` skips [C].

**Run only the dyad influence analysis (v4.4)**  Note: 1–2 minutes (requires R; does not redo the permutations)
```bash
STAGE=influence bash run_R1_2_trialfix.sh
```
`08_R1-2_p3_dyad_influence.R` reads the existing `out_periods_grid_{visible,invisible}` and, with the same model as 04 (`z_S ~ condition * period + (1 | dyad_id)`), produces the following.
- (1) leave-one-dyad-out: re-estimates the P3 condition difference (visible − invisible, emmeans) leaving out each of the 20 dyads in turn. `loo_estimates.csv`. The summary shows the range of estimates, whether all are positive, the maximum p and the most influential dyad.
- (2) Sign distribution: number of dyads with z_S > 0 in P3 (mean over 3 sessions; 10 vs 10) and number of sessions (30 vs 30), exact binomial (sign) test, Wilcoxon signed-rank test (vs 0), Wilcoxon rank-sum test of dyad means between conditions. `sign_tests.csv`, `dyad_level_zS.csv`.
- (3) Figure `fig_dyad_influence_delta01.pdf/.png`: (a) dyad-mean z_S (large dots) and session values (small dots) per period × condition, (b) leave-one-dyad-out estimates with 95% CI. In environments without `patchwork`, (a) and (b) are written to separate files.
- The summary is `dyad_influence_results.txt`. The statement in the manuscript and letter that "the effect is a small shift spread across many pairs rather than strong synchrony in a few dyads" is taken from this output (do not write it if the results do not support it).
- Also run automatically after 07 in the R step of `STAGE=all`. Old outputs are moved aside into `pre_redo_backup_<datetime>/`.

### step 3: null-distribution checks (R2_2_2_5_2_10_3_2_3_3)  Note: tens of minutes with 10,000 permutations
```bash
cd ../R2_2_2_5_2_10_3_2_3_3
bash run_R2_2_trialfix.sh          # 01 → 02 (δ=0.1: w=0.35/0.2/0.5, δ=0.2: w=0.35/0.5) → 03 → 02b ×4 → 04 (R) ×4, all in one go
# SKIP_R=1 bash run_R2_2_trialfix.sh   # environments without R. The Rscript commands to type by hand are printed at the end
```
Since 02b and 04 (LMM) were also included in v3, there are no manual steps. The outputs fit in the 3 folders `out_symmetry/`, `out_jitter/` (including `lmm_in_<setting>/`, `lmm_out_<setting>/`) and `out_shift/`.
The former `run_R2_2_delta02_addon.sh` has been merged into v3 and is no longer needed (if you already ran it, the results are identical, because the same calls are made with the same seed).

### step 4: independent implementation, grid null, IRR (R2_4_3_1)
```bash
cd ../R2_4_3_1
bash run_R2_4_trialfix.sh        # if IRR_final/ exists, also runs zS_reproduction and irr_agreement
STAGE=reliability bash run_R2_4_trialfix.sh   # v4.2: does not rerun nullcheck / jitter; rebuilds only the zS reproduction, agreement rates and R3-1
STAGE=icc_ci bash run_R2_4_trialfix.sh        # v4.3: recomputes nothing; only writes the 95% CI of ICC(A,1) from the existing per_video.csv (seconds)
```
- Addition in v4.3 (`R2-4_zS_icc_ci.py`): letter R2-4 states that "the inter-coder ICC of session-level z_S is low (ICC ≤ 0.27)", and
  a reviewer could ask "with only point estimates and n = 10–12, how wide are the intervals?". So this script reads the zS reproduction output
  `R2-4_zS_{inter,intra}_d<δ>_per_video.csv` (z_S per video × period × coder) and attaches to the ICC(A,1) per period (P1/P2/P3, plus the 3 periods pooled for reference;
  the same estimator as in R2-4_zS_reproduction.py) a **bootstrap 95% confidence interval that resamples videos with replacement** (default 10,000 replicates,
  seed 20260827, changeable with `ICC_BOOT=`); for reference it also lists the McGraw & Wong F-method interval and Pearson r.
  The F method is meaningless when the point estimate is negative, so the bootstrap interval is used in the letter. No re-simulation is done, so it finishes in seconds,
  and it also runs automatically at the end of `bash run_R2_4_trialfix.sh` / `STAGE=reliability`. If step 4 has already been completed, `STAGE=icc_ci` alone is enough
  (`STAGE=icc_ci` does not move aside or delete any other output; only the existing `R2-4_zS_icc_ci.*` are moved aside into `pre_redo_backup_*/`).
  - (v4.3.1) The input `R2-4_zS_{inter,intra}_d<δ>_per_video.csv` is written by R2-4_zS_reproduction.py together with summary.txt.
    If it is not in this folder, `pre_redo_backup_*` / `pre_trialfix_backup_*` are searched newest first and used (this is reported with [INFO];
    please check that it is output from the same run as summary.txt); if it is nowhere, diagnostics are printed and the script stops. If you know the location, use
    `PER_VIDEO_DIR=path STAGE=icc_ci bash run_R2_4_trialfix.sh`. If it truly exists nowhere, rebuild the zS reproduction with `STAGE=reliability`
    (summary.txt is then regenerated as well, so check in the audit that the ICC point estimates in the letter still agree).
  - Output: `R2-4_zS_icc_ci.txt` (a table of mode × δ × period, and for the letter one line each with "the range of the lower and upper limits of the P1–P3 intervals"), `R2-4_zS_icc_ci.csv`.
  - Upload: `R2-4_zS_icc_ci.txt` (the ICC 95% CI values reported for R2-4 are taken from this file).
- Changes in v4.2 (the kick calls `../00_trialfix/01_patch_scripts.py .` at the start to apply them automatically; safe to repeat):
  - The null distribution of the zS reproduction (R2-4_zS_reproduction.py) is changed to the **grid-aligned shift** used in the General note of the letter (`--null grid`,
    shifts by integer multiples of 0.1 s only; the same design as R2-4_zS_nullcheck.py). The old continuous uniform shift remains available
    as `--null continuous` but is not used by the kick.
  - Run for **all of 0.1 / 0.2 / 0.3 / 0.5 / 1.0**, not only δ = 0.1. Output names are `R2-4_zS_{inter,intra}_d<δ>_summary.txt`
    (the old `R2-4_zS_inter_summary.txt` is the continuous-null, δ = 0.1 result and is not used; it is moved aside into `pre_redo_backup_*/`).
  - The hard-coded "observed P3 contrast = +0.xxx (SE 0.xxx)" at the end of R3-1_power_sensitivity.R is changed so that
    `observed_P3_contrast` / `observed_P3_se` from `R3-1_power_params.txt` (extracted from the emmeans table in lmm_results.txt)
    are passed via `--observed / --observed_se` and printed.
- Prerequisite: step 2 (R1_2) must be finished (reads `../R1_2/lmm_grid_out/lmm_results.txt`).
- What it runs: R2-4_zS_nullcheck (grid null vs continuous null, all δ) → R2-4_jitter_simulation (δ = 0.1) →
  coder reliability in two runs (inter: second coder `IRR_final/`, intra: re-annotation by the original coder `../k_ICC/ICC_result/`),
  each running R2-4_zS_reproduction and R2-4_irr_agreement, with `_inter` / `_intra` appended to the outputs to keep them apart
  (both read only Chew columns, so there is no Trial-row contamination and no preprocessing is needed. If the re-annotation folder is elsewhere, give it with `INTRA_ROOT=...`) →
  reads the R3-1 variance components (resid_sd, dyad_sd) and the observed P3 contrast automatically from `lmm_results.txt` and records them in `R3-1_power_params.txt` →
  if Rscript is available, also runs R2-4_nullcheck_lmm.R and R3-1_power_sensitivity.R.
- There is no place where values are entered by hand. If R lives in a different environment, run the 2 Rscript commands printed on screen as they are
  (the arguments are printed with the values from `R3-1_power_params.txt` already filled in).
- Output: `R2-4_nullcheck_*`, `R2-4_jitter_sim_*`, `R2-4_zS_{inter,intra}_d<δ>_*`, `R2-4_irr_results_{inter,intra}_*`, `R2-4_nullcheck_lmm*`, `R3-1_power*`, `R3-1_power_params.txt`
- Upload (after the v4.2 rerun): `R2-4_zS_inter_d{0.1,0.2,0.3,0.5,1.0}_summary.txt`, `R2-4_zS_intra_d{0.1,0.2,0.3,0.5,1.0}_summary.txt`,
  `R3-1_power_params.txt`, `R3-1_power_summary.txt` (the R3-1 output names are unchanged).

### step 5: degrees of freedom, multiple comparisons, dip (R2_8_2_9_2_11)  Note: after step 1
```bash
cd ../R2_8_2_9_2_11
bash run_R2_8_trialfix.sh        # input: ../Micro_Analysis_LMM, ../gridnull_main/dt_trialfix.csv
```

- 03_R2-11_dip_test.py requires the `diptest` package (`pip install diptest`).
- Note: at δ = 0.3 the Δt distribution is asymmetric at the ±0.3 s edges (more on the −0.3 s side). This comes from the NLOPM
  rule "choose the later B when equidistant" and remains even with tolerance 0 (a property of the original algorithm).
  At δ = 0.1 it is symmetric. For reference if this is mentioned in the R2-11 text.

### step 6: R-only steps
- Knit `04_Micro_Analysis_LMM.Rmd` (Table 1, z_S condition comparison in the main text).
- `gridnull_main/06_bayes_gridnull_fit.R fit 0.1,0.2,0.3,0.5,1.0` (the plain version; do not use v2/v3).
- Afterwards, verify the Stan output (v4.2):
```bash
cd gridnull_main
STAGE=bayes_check bash run_gridnull_main_trialfix.sh
```
  - Runs `check_beta_intervals_v2.py` on all `stan_fits/summary_z_S_t_delta*_human.csv` and saves the result to `beta_intervals_check.txt` (bin numbers).
  - `beta_intervals_seconds.py` converts the bin numbers to video seconds and gives the total seconds per P1/P2/P3 → `beta_intervals_seconds_d<δ>.txt`.
    Bin i is the i-th of the t_center values that appear in `Micro_Analysis_Bayesian/{vis,invis}_<δ>.csv`, sorted ascending
    (the same construction as `bin_id = as.integer(factor(t_center))` in 06). Stops with an error if the number of bins does not match.
  - If Rscript is available, draws Fig 4e with `07_plot_figure4e_gridnull.R` → `fig4e_gridnull/fig4e_beta_delta01.{pdf,png}`, `_bins.csv`.
    The plotting code is the corresponding chunk of `04_Micro_Analysis_Bayesian.Rmd` from the public repository as is, with only the x axis changed from bin number → video seconds.
  - The "… s of which … s in P3" and "reverse direction … s" in the letter and the Fig 4e caption are replaced with the values from this output.
  - Upload: `beta_intervals_check.txt`, `beta_intervals_seconds_d*.txt`, `fig4e_gridnull/fig4e_beta_delta01.png`.
- The R steps are **untested**, because there is no R in our environment. Input file names and column layouts are kept the same as in the old version
  (the Bayesian input has only the `t_center` + `z_S_t__<pair>` columns; the number of bins changes from 176 to about 171).

---

## 4. Results that change with the fix (checklist for revising the manuscript and letter)

| Item | Status |
|---|---|
| Table 1 / z_S condition comparison in the main text (all δ) | Redo (n_A, n_B, n_matched, z_S all change) |
| Table 2 (proportion with z_S > 1.64) | Redo (05) |
| Fig 2 (Δt distribution), R2-11 dip test | Rebuild from `dt_trialfix.csv`. The spurious Δt = 0 matches (from Trial) disappear |
| Fig 4a–d (z_S(t) time series), Fig 4e (β(t)), Fig S2 (change points) | Time axis stays in video time. Each pair's series starts at its actual start point (3.2–19.5 s), so few pairs are covered before 20 s. "Synchrony right after onset" must be re-examined |
| R1-2: LMM for P1/P2/P3, pseudo pairs, lag gradient, sensitivity analyses | Redo. In particular z_S ≈ 2.5 in P1 is very likely an artefact. Unify the period definition to "video time from session start" and report the number of pairs with `sufficient = 0` in P1 together with the results of sensitivity analysis [C] (relative origin) |
| R2-2/2-5/2-10/3-2/3-3: null symmetry, jitter/ISI null, reactive shift | Redo |
| R2-4/3-1: grid null vs continuous null, jitter simulation, independent-implementation reproduction, IRR of z_S | Tolerance fix only (no Trial-row contamination), but the numbers change. R3-1 power recomputed with the new variance components |
| R2-8/2-9: df verification, multiple-comparison correction | Rerun with replaced inputs. **The values in Table 1 (0.645 / t 2.400 etc.), its footnote and STAR Methods lines 1082–1092 are still the old version**: replace them with the values of `df_verification_table1.txt` (z_S δ = 0.1: 0.671 / SE 0.278 / t(18) = 2.42 / p = 0.027) and the Kenward–Roger footnote (df = 18 for all 15 tests) (same wording as the Author Action in letter R2-8) |
| Figure S1 / S2 numbering in the manuscript | The 20260823 manuscript calls the period-boundary figure "Figure S1" (around lines 241 and 373) and still uses the figure from the old `R1_2/01_segment_boundaries.py` + `08_plot_figure_S1.py` (with the data-driven text about change points at 20 s / 60 s). In the letter, S1 = the R2-11 dip figure (`R2_8_2_9_2_11/dip_figure_out`, same new dt as Fig 2) and S2 = the period figure (`gridnull_main/fig_S2_gridnull`, 08_plot_figure_S2_gridnull_v2). On the manuscript side, **replace the figure with the fig_S2_gridnull output, renumber S1 → S2, and align the text in lines 236–241 with letter R1-2 (paragraph 39)**. The old R1_2/01 and 08 scripts are not called from the run scripts, so they are not part of the 98_export public set either |
| STAR Methods | State explicitly: "only Chew columns are read as event times", "a tolerance of 1e-6 s in the comparison with δ", "windows start at the integer second obtained by rounding up the start of each pair's overlap interval, in 1 s steps (video time)", "periods are defined in video time from session start and intersected with the overlap interval" |
| k_ICC (κ, ICC), R2_6_3_5 | Not affected as long as bite times are not passed through NLOPM (to be confirmed) |

For reference (numbers verified on our side, δ = 0.1, 60 pairs):
- 2,593 spurious events removed (177 files; 1,909 in the 120 pair files); start of the overlap interval 1.0 s → 3.2–19.5 s;
  total matches over the whole interval old 4,560 → fixed 5,132 (breakdown: about 890 spurious matches from Trial disappear, about 1,460 are recovered by the tolerance).
  The 22 chews recovered from header-less columns add a few matches near the end for 2 invisible-condition pairs (5,127 in v2 → 5,132).
- Summary of the main analysis at δ = 0.1 (step 1, 10,000 permutations, as of v2, before the 22 values were recovered): visible condition mean z_S 0.35 (published version 2.17), proportion with z_S > 1.64 0.07 (published version 0.73);
  invisible condition mean z_S −0.34, proportion 0.03.
- Change points (08, video time): 62 s and 161 s for K=2; 62 / 137 / 146 s for K=3. The old 9 s / 20 s disappear and
  a boundary around 60 s re-emerges in a data-driven way (please confirm the interpretation with the co-authors).
- Fig 4a–d: peak of the visible condition at 79.5 s (z_S 0.43). The invisible condition shows a peak of 0.83 at 6.5 s, but
  this bin is covered by only 8 pairs and is not trustworthy (see below). The large early peak of the old version is gone.
- Start of the time axis: because each pair's series starts at its actual start point, few pairs are covered before 10 s
  (2–8 pairs at 6.5 s, 23–28 pairs at 10.5 s, almost all pairs from 15.5 s). The figure from 09 is not truncated by coverage, so
  please decide whether to drop the first few bins with a `--tmin`-like option or to show the number of covered pairs alongside the figure.
  08_changepoint has always excluded low-coverage bins with `--min_frac 0.5`. The end (up to 176.5 s) has 25–28 pairs and is not a problem.
- The x-axis label "Time in session (s)" refers to video time, so it can stay as it is.

---

## 5. Safeguards

**The most reliable method**: before running, duplicate the whole of `iscience_revision_analysis/` with a date suffix
(e.g. `cp -r iscience_revision_analysis iscience_backup_20260908`). This follows the same convention as the earlier `iscience_backup_20260829`.
If you do this, the old results are guaranteed to survive even if one of the individual backups below is missed.

What each kick moves aside before overwriting (moved or copied into `pre_trialfix_backup_<datetime>/`):

| Folder | Moved aside | Not moved aside, overwritten (manual R steps) |
|---|---|---|
| gridnull_main | `out_grid_*_d*`, `out_changepoint_grid*`, `fig_S2_gridnull`, `fig4ad_gridnull`, `fig2_timelags`, `coverage_grid`, `../Micro_Analysis_LMM/*_chew_sync_summary_*.csv` (copy), `Micro_Analysis_Bayesian/{vis,invis}_*.csv` (copy), `stan_fits/` (copy), knitted LMM html (copy) | output of `check_beta_intervals_v2.py` (rename when running, or save under another name first) |
| R1_2 | `out_periods_grid_*`, `out_pseudo_grid_*`, `lmm_grid_out*` (including `lmm_results.txt`), `dt_grid_out*`, `pseudo_interaction_grid*` | the sensitivity-analysis outputs of 06, `*_excl_dyad`, `*_minus5`, `*_plus5`, are moved aside because their names fall under `out_periods_grid_*` / `lmm_grid_out*` |
| R2_2 | `out_symmetry*`, `out_jitter*`, `out_shift*` | if 02b/04 were run manually into a different outdir, it depends on that name |
| R2_4_3_1 | `R2-4_nullcheck_*`, `R2-4_jitter_sim_*`, `R2-4_zS_{inter,intra}_*`, `R2-4_irr_results_{inter,intra}_*`, `R2-4_nullcheck_lmm*`, `R3-1_power*` | none |
| R2_8 | `multiplicity_out`, `dip_test_out`, `dip_figure_out`, `df_verification_out` | none |

- `Micro_Analysis_LMM/*_chew_sync_summary_*.csv` are read by the downstream R code under the same names, so a copy is moved aside and the files are **replaced under the same name**.
  To tell them apart from the old version, look in the backup location (`gridnull_main/pre_trialfix_backup_<datetime>/Micro_Analysis_LMM/`).

- The original data `annotation_raw_data/` is read-only. The fixed version goes into the separate folder `annotation_raw_data_chewonly/`.
- Each rewritten .py leaves `<name>.pre_trialfix.bak` in the same place.
- Each kick moves the old outputs into `pre_trialfix_backup_<datetime>/` before running.
- Each kick checks at start "is the data the _chewonly version" and "do the scripts carry the markers", and stops if not yet patched.
- `00_trialfix/02_verify_fix.py` actually imports the functions of the patched 03 and confirms that they agree with an independently written
  verification implementation on all 60 pairs (ALL OK means the fix is effective).

---

## 6. Checks on our side

- Python steps: step 0 (all folders), step 1 (δ = 0.1, 10,000 permutations, both conditions, video time), steps 2 (including [C]), 3, 4 and 5
  run with a reduced number of permutations; all completed without error.
- v4.2: the GRID-NULL / OBSERVED-ARG patches of 01_patch_scripts.py (confirmed that the second run skips), step 4 with `STAGE=reliability`
  (synthetic re-annotation folder, 2 values of δ, 200 permutations), `STAGE=bayes_check` (with the summary csv and vis/invis_0.1.csv you sent)
  were run on our side and completed. **The numbers obtained here are for smoke testing only; all values used in the letter are taken from your run logs** (as agreed).
- v4.3: `R2-4_zS_icc_ci.py` was tested with a per_video.csv reconstructed from the per-video z_S in the `R2-4_zS_{inter,intra}_d<δ>_summary.txt` you sent
  (confirmed that the point estimates agree with the ICC in the summary; 4 seconds for 2,000 bootstrap replicates). The `STAGE=icc_ci` kick
  (old output moved aside into pre_redo_backup on the second run) was also checked. **The interval values are taken from your run results** (as agreed).
- R steps: not run (no R environment). Input/output names and columns are kept the same as in the old version, but please inspect the output by eye the first time.
  07_plot_figure4e_gridnull.R was also not run (the chunk of the public Rmd moved over as is; `size =` is written so that it also works with older ggplot2).
- v4.3.3: 09 (v2.1) was run on the δ = 0.1 `timelines_all_master.csv` you sent, in the 3 combinations `--sharey row` / `none` and `--ylim_min_pairs 0` / `15`,
  confirming that peaks.txt and source_data.csv agree with the old version and that the y axes are aligned per row. The `STAGE=fig4ad` kick was also checked.
- When inspecting the products, the quickest way is to look at `chewonly_report.txt` (removal counts per file), `verify_report.txt` and `coverage_grid/coverage_summary.txt`.

## 7. Assembling the folder for GitHub release (98_export_for_github.py)

The working folder mixes intermediate products of the re-analysis, letter drafts and so on, so **instead of cleaning the working folder,
only the items to be released are copied next to it and assembled there** (the working folder is never modified).

```bash
cd iscience_revision_analysis
python 98_export_for_github.py            # only prints the list of what would be copied (export_manifest.txt)
python 98_export_for_github.py --apply    # actually copies into ../revision_export/
```

Only 4 kinds of items are copied.

| Kind | Selection rule |
|---|---|
| Scripts | `run_*_trialfix.sh` in each folder and the .py / .R / .Rmd / .stan / .sh files called from them (extracted automatically; if a `_gridnull` version exists, the plain one is excluded) |
| Outputs | files and folders in that folder newer than `pre_trialfix_backup_<datetime>` (= created after the fix) |
| Data | `annotation_raw_data/`, `annotation_raw_data_chewonly/` |
| Procedure | the whole of `00_trialfix/`, `README_procedure.md`, `.gitignore` |

`.docx`, `.zip`, `.log`, `.bak`, `upload_*`, backup folders, the discarded-scripts folder, prompt .md files and Stan fit objects are
not copied. Items that cannot be decided automatically (scripts updated after the fix but not called from `run_*.sh`, `stan_fits*`) are
listed at the end of the manifest as `[CHECK]`; copy by hand only those you need.

`99_tidy_for_github.sh` is a helper that moves the backup folders, `.bak` files and `upload_*` on the working-folder side into `../_local_only/`;
it is not required for building the release folder.

- The principle for `run_*_trialfix.sh` is "one script per folder, all steps top to bottom".
  Steps added later should be merged into this one script rather than becoming separate scripts (R2_2 in v3 is an example).
- `00_trialfix/` and `annotation_raw_data_chewonly/` are included in the release. They are the very data-correction procedure disclosed in the
  General note of the letter, so placing them next to the original data (`annotation_raw_data/`) is the most transparent option.
