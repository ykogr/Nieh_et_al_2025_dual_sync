# R2_6_3_5 additional scripts (r26_condiff_kit v2.2, 2026-09-14)

Either of the following locations is fine (v2.1).

- (A) `R2_6_3_5/r26_condiff_kit/` — the zip left unpacked inside R2_6_3_5 (the current folder layout).
- (B) directly under `R2_6_3_5/` — the 5 files placed at the same level as 03_sensitivity_drift_beta.R, 04 and 05.

`run_R2_6_condiff.sh` looks for `sensitivity_out/` first in its own folder and then in the parent folder to determine the location of R2_6_3_5 (`ANALYSIS_DIR`),
and prints `ANALYSIS_DIR (R2_6_3_5) = …` on the first line. The defaults of `FIT_DIR`, `MICRO_DIR`, `XSCALE_OUT` and `DATA_DIR` are all
relative to `ANALYSIS_DIR`, so with either (A) or (B) the outputs go to `R2_6_3_5/sensitivity_out/rhat_raw.*` and `R2_6_3_5/xscale_out/`
(letter_v4_kit and manuscript_kit read these locations). If `sensitivity_out/` exists in neither place, the script stops and asks for `ANALYSIS_DIR=<R2_6_3_5>`.

v2.1 → v2.2 (07_rhat_main.R only): because `fit$summary()` was applied to all variables, R̂ was also computed for the generated quantities `y_full_out`, `y_rep` and `log_lik`
(each of length N = number of windows, 3N variables in total), which took tens of minutes. v2.2 summarises only the variables obtained from `fit$metadata()$stan_variables`
after removing these 3 (all variables of the parameters block = y_mis, mu0, mu_pair, sigma_pair, mu_sess_init, sigma_sess_init, mu_state, sigma_mu,
beta0, beta_tan, sigma_beta, sigma_obs; beta_state from transformed parameters; and lp__). The generated quantities are
deterministic functions of the parameters and are not covered by the manuscript statement "all parameters showing R̂ < 1.1". The log records the summarised variable names, the excluded variable names and the elapsed seconds.
If an Rscript run had been left hanging, terminate it with Ctrl-C and rerun (no partial results are written).

The only change from v2 to v2.1 is this folder detection (v2 assumed that `sensitivity_out/` is in the same folder as the run script, so with layout (A)
it stopped with "sensitivity_out/fit_raw_visible_pair.rds not found"). The 3 R scripts are unchanged.

| File | Role | STAGE |
|---|---|---|
| `run_R2_6_condiff.sh` | Entry point (`STAGE=… bash run_R2_6_condiff.sh` in this folder) | — |
| `06_condition_difference_posterior.R` | Posterior distribution of the condition difference Δ (R2-6 / R3-5). Unchanged since v1 | `condiff` |
| `07_rhat_main.R` | Measured maximum R̂ of the main-analysis (raw mode) fit (co-author comment 10) | `rhat` |
| `08_xscale_micro_macro.R` | Correlation between dyad-level micro P3 z_S and macro coupling (co-author comment; measurement for R2-3) | `xscale` |

## Running (zsh, in the folder containing run_R2_6_condiff.sh)

    STAGE=rhat   bash run_R2_6_condiff.sh
    STAGE=xscale DATA_DIR=/path/to/Annotated_Data bash run_R2_6_condiff.sh
    STAGE=all    DATA_DIR=/path/to/Annotated_Data bash run_R2_6_condiff.sh   # 06 → 07 → 08

- `DATA_DIR` is the GitHub `Annotated_Data` (containing `visible_pair_Human/{A,B}`, `invisible_pair_Human/{A,B}`).
  If omitted, `$ANALYSIS_DIR/../../Annotated_Data` and then `$ANALYSIS_DIR/../Annotated_Data` are tried; if neither exists the script stops and asks (it does not guess).
- `MICRO_DIR` (default `$ANALYSIS_DIR/../R1_2`) must contain `out_periods_grid_visible/chew_sync_summary_periods_master.csv` and
  `out_periods_grid_invisible/…` (outputs of `run_R1_2_trialfix.sh` in R1_2).
- `FIT_DIR` (default `$ANALYSIS_DIR/sensitivity_out`) must contain `fit_raw_visible_pair.rds` and `fit_raw_invisible_pair.rds` (saved by 03; the same files 06 reads).

## STAGE=rhat (07_rhat_main.R)

Collects R̂ of all parameters from the main-analysis fit and writes it to `sensitivity_out/rhat_raw.txt`.

- "max R-hat, all sampled parameters" = maximum over all variables of the parameters block (including y_mis = imputed values of missing data) + beta_state + lp__ (v2.2: generated quantities excluded). The measured value to accompany the manuscript statement "all parameters showing R̂ < 1.1"
- "model parameters" = maximum excluding y_mis and lp__. "beta_state" = maximum over the 13 values of β(t) only.
- At the end it prints draft wording for the manuscript STAR Methods (3 decimal places). Flow of the value: letter_v4_kit (logs_v46.py) reads `R2_6_3_5/sensitivity_out/rhat_raw_summary.csv` to generate letter paragraph 715, and manuscript_kit builds E61 (the parenthetical in the MCMC paragraph) from that paragraph 715. Nothing is transcribed by hand.
- R̂ for the micro side (β(t) of Figure 4e) is already given as "max Rhat over beta_state" in `gridnull_main/beta_intervals_seconds_d0.1.txt`, so no additional run is needed.

## STAGE=xscale (08_xscale_micro_macro.R)

For each dyad (n = 20, 10 per condition), builds the following 2 quantities and correlates them.

- micro index: δ = 0.1 s, z_S in P3 averaged over the 3 sessions (from the R1_2 master table).
- macro index: with the same loader, 10 s bins, 60 s window (bins 6–18) and missing-data rules as the main analysis, the mean Fisher-z window correlation of the real pair
  minus the mean of time-preserving shuffles that pair the same dyad's A with the B of other pairs (a dyad-level version of β(t) = E_real − E_shuffled).
  For reference, the same tests are also reported for the raw real-pair mean (`macro_real`).

Output `xscale_out/`:

| File | Content |
|---|---|
| `xscale_dyad_table.csv` | The 2 indices per dyad (the points in the figure) |
| `xscale_micro_macro.txt` | Human-readable log. Numbers for the manuscript and letter come from here |
| `xscale_micro_macro.csv` | One row per test (pooled / within-condition / visible / invisible × macro_index / macro_real) |
| `xscale_micro_macro.pdf/.png` | Scatter plot |

The main result to report is **the within-condition Spearman ρ and permutation p (condition means subtracted, permutation within condition)**.
The pooled correlation would pick up the condition difference itself ("visible has higher micro and lower macro"),
so it is listed only as a reference value. 10,000 permutations, seed 20260914 (results are reproducible).

The `macro_index` in the condition-means rows (`condition means`) should nearly coincide with the 13-bin mean of β(t) in Figure 3 (mean_beta in reported_values.csv).
If they differ, suspect a wrong DATA_DIR.

On our side we only confirmed that 08 runs to completion on the GitHub public data (20 dyads matched). Only the PI-side run results are used for the numbers in the manuscript and letter.

---

# 06_condition_difference_posterior.R (R2-6 / R3-5: posterior distribution of the condition difference) — unchanged from v1

## Purpose
Replace the old manuscript's "visual cue deprivation significantly enhanced rhythmic coupling (p < 0.001)"
(main text around Figure 3, Figure 3 legend) with a quantity computed directly from the posterior of the Bayesian state-space model.

    Δ = mean_b β_invisible(b) − mean_b β_visible(b)   (b = 1..13 bins, window end 60–180 s)

Reported are the posterior mean of Δ, the 95% credible interval and the posterior probability P(Δ > 0).
The t test treating bins as independent samples (the presumed source of the old p < 0.001) is not used.

## Location
`R2_6_3_5/06_condition_difference_posterior.R`
(same level as 03_sensitivity_drift_beta.R, 04_check_divergence_impact.R, 05_reported_values.R)

## Prerequisites
`fit_<mode>_<condition>.rds` saved by `03_sensitivity_drift_beta.R`
(`sensitivity_out/fit_raw_visible_pair.rds` etc.). No refitting is done.

## Running (zsh, inside R2_6_3_5)
    STAGE=condiff bash run_R2_6_condiff.sh      # = Rscript 06_condition_difference_posterior.R sensitivity_out

- The 2nd argument restricts the mode (e.g. `raw` only → `Rscript 06_condition_difference_posterior.R sensitivity_out raw`)
- The 3rd argument changes the output location (if omitted, written to the folder containing the fits)

## Output (sensitivity_out/)
| File | Content |
|---|---|
| `condition_difference_<mode>.txt` | Human-readable log. Mean, median, SD, 95% CrI, P(Δ > 0) of Δ and the per-bin differences |
| `condition_difference_all.csv` | One row per mode (values for the manuscript and letter come from this row) |
| `condition_difference_bins_all.csv` | One row per mode × bin (CrI of Δ_b and P(Δ_b > 0)) |

## Use in the manuscript (values always from the raw row of condition_difference_all.csv)
- Main text: `(p < 0.001)` → `(posterior difference Δ = X.XX [X.XX, X.XX], P(Δ > 0) = 0.9xx)`
  - When `n_draws_delta_le0` is 0, write e.g. `P(Δ > 0) > 0.9999` (within the resolution of 20,000 draws)
- Figure 3 legend: add the same Δ and P(Δ > 0), delete "p < 0.001"
- STAR Methods: add one sentence to the effect that "the condition difference was evaluated as the posterior distribution of Δ by pairing the posterior draws of the two independently fitted conditions"

## Key points of the computation (for undergraduates)
- The 2 conditions are fitted separately, so the posterior draws are independent of each other. For each draw index d, forming
  Δ(d) = (13-bin mean of invisible) − (13-bin mean of visible) gives a sample from the posterior of Δ.
- The order of pairing should not affect the result, so P(Δ > 0) for a randomly permuted pairing (seed 20260910) is also
  output as a check. The two should nearly agree.
- The posterior mean of Δ should equal the difference of mean_beta in reported_values.csv (raw: 0.1556 − (−0.0347) = 0.190).
  If it does not, suspect that the wrong fits were used.
