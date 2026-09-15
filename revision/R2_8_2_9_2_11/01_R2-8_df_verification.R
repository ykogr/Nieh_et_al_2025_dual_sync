#!/usr/bin/env Rscript
# =====================================================================
# 09 (R2-8): Verification of denominator degrees of freedom in Table 1
#
# Purpose:
#   R2-8 questions t(18) for the zS LMM at delta = 0.1 s, suspecting that
#   pair-level random intercepts were not properly accounted for.
#   This script verifies, for every Table 1 cell (3 DVs x 5 deltas), that:
#     (1) the LMM  dv ~ condition + (1 | pair_id)  fitted to the 60
#         session-level values (20 dyads x 3 sessions) yields the reported
#         estimate/SE/t with SATTERTHWAITE denominator df;
#     (2) KENWARD-ROGER df (as used in the period-resolved analyses added
#         for R1-2) give numerically identical results in this balanced
#         design (10 dyads per condition, 3 sessions each);
#     (3) a dyad-level analysis -- averaging the 3 sessions per dyad and
#         running a two-sample t-test on the 20 dyad means -- reproduces
#         t with df = 18, confirming that the effective unit of inference
#         is the dyad (df = 20 - 2 = 18), not the 59 individuals.
#
# Input:
#   arg1 = directory containing
#          {delta}_chew_sync_summary_{visible|invisible}_pair_Human.csv
#          (the folder used by 04_Micro_Analysis_LMM.Rmd, e.g. Micro_Analysis_LMM)
#   arg2 = (optional) output directory (default: arg1/df_verification_out)
#
# Output:
#   df_verification_table1.csv : one row per DV x delta with Satterthwaite,
#                                Kenward-Roger, and dyad-mean t-test results
#   df_verification_table1.txt : human-readable summary incl. equality checks
#
# Requirements: lme4, lmerTest, pbkrtest (for Kenward-Roger),
#   dplyr/tidyr/readr/purrr/stringr/tibble (all part of tidyverse).
#   If pbkrtest is missing: install.packages("pbkrtest")
#
# Usage:
#   Rscript 01_R2-8_df_verification.R Micro_Analysis_LMM df_verification_out
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(lme4)
  library(lmerTest)
})

HAS_PBKR <- requireNamespace("pbkrtest", quietly = TRUE)
if (!HAS_PBKR) {
  message("[NOTE] 'pbkrtest' not installed; Kenward-Roger columns will be NA. ",
          "For KR df: install.packages('pbkrtest')")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 01_R2-8_df_verification.R <indir> [outdir]")
indir <- args[1]
outdir <- if (length(args) >= 2) args[2] else file.path(indir, "df_verification_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- load session-level summaries (same convention as 04_Micro_Analysis_LMM.Rmd) ----
files <- list.files(indir,
                    pattern = "^\\d+\\.\\d+_chew_sync_summary_(visible|invisible)_pair_Human\\.csv$",
                    full.names = TRUE)
if (length(files) == 0) stop(paste("no summary CSVs found in:", indir))

d <- map_dfr(files, function(file) {
  fname <- basename(file)
  delta <- str_extract(fname, "^\\d+\\.\\d+")
  condition <- str_remove(str_extract(fname, "(visible|invisible)_pair_Human"), "_Human")
  df <- read_csv(file, show_col_types = FALSE)
  df$delta <- as.numeric(delta)
  df$condition <- condition
  df
}) %>%
  mutate(
    session = as.integer(str_extract(pair_id, "\\d{2}$")),
    pair_id = str_replace(pair_id, "_\\d{2}$", "")
  )

dvs <- c("z_S", "STR", "Max_episode")

extract_cond_row <- function(coefs) {
  i <- grep("^condition", rownames(coefs))
  coefs[i, , drop = FALSE]
}

rows <- list()
for (dv in dvs) {
  for (delta_val in sort(unique(d$delta))) {
    dd <- d %>% filter(delta == delta_val)
    f <- reformulate("condition + (1 | pair_id)", response = dv)
    m <- lmer(f, data = dd)

    # (1) Satterthwaite (lmerTest default; what Table 1 reports)
    cs <- extract_cond_row(coef(summary(m, ddf = "Satterthwaite")))

    # (2) Kenward-Roger
    if (HAS_PBKR) {
      ck <- extract_cond_row(coef(summary(m, ddf = "Kenward-Roger")))
    } else {
      ck <- matrix(NA_real_, 1, 5)
    }

    # (3) dyad-mean two-sample t-test (Student, equal variances -> df = 18)
    dm <- dd %>%
      group_by(pair_id, condition) %>%
      summarise(dv_mean = mean(.data[[dv]]), .groups = "drop")
    tt <- t.test(dv_mean ~ condition, data = dm, var.equal = TRUE)
    # condition levels sort as invisible_pair < visible_pair; t.test reports
    # mean(invisible) - mean(visible), while the LMM coefficient is
    # visible - invisible, so flip signs for comparability.
    t_dyad <- -unname(tt$statistic)
    est_dyad <- unname(tt$estimate[2] - tt$estimate[1])  # mean(visible) - mean(invisible)

    rows[[paste(dv, delta_val)]] <- tibble(
      dv = dv, delta = delta_val,
      n_sessions = nrow(dd), n_dyads = n_distinct(dd$pair_id),
      estimate_satt = cs[1, "Estimate"], se_satt = cs[1, "Std. Error"],
      df_satt = cs[1, "df"], t_satt = cs[1, "t value"], p_satt = cs[1, "Pr(>|t|)"],
      df_kr = ck[1, 3], t_kr = ck[1, 4], p_kr = ck[1, 5],
      estimate_dyadmean = est_dyad, t_dyadmean = t_dyad,
      df_dyadmean = unname(tt$parameter), p_dyadmean = tt$p.value
    )
  }
}
res <- bind_rows(rows)

# equality checks
res <- res %>%
  mutate(
    # Satterthwaite df are obtained by numerical approximation, so allow
    # tiny rounding differences (df to 1e-3, t and p to 1e-6).
    kr_equals_satt = if (HAS_PBKR)
      abs(df_kr - df_satt) < 1e-3 & abs(t_kr - t_satt) < 1e-6 & abs(p_kr - p_satt) < 1e-6
      else NA,
    dyadmean_df_is_18 = abs(df_dyadmean - 18) < 1e-9
  )

write_csv(res, file.path(outdir, "df_verification_table1.csv"))

sink(file.path(outdir, "df_verification_table1.txt"), split = TRUE)
cat("=====================================================================\n")
cat("R2-8 df verification:", format(Sys.time()), "\n")
cat("Model: dv ~ condition + (1 | pair_id); coefficient = visible - invisible\n")
cat("=====================================================================\n\n")
res %>%
  select(dv, delta, n_sessions, n_dyads,
         estimate_satt, se_satt, df_satt, t_satt, p_satt,
         df_kr, t_kr, p_kr,
         t_dyadmean, df_dyadmean, p_dyadmean,
         kr_equals_satt, dyadmean_df_is_18) %>%
  as.data.frame() %>%
  print(digits = 4, row.names = FALSE)
cat("\nChecks:\n")
cat(sprintf("  All Kenward-Roger results identical to Satterthwaite: %s\n",
            ifelse(all(res$kr_equals_satt %in% TRUE), "YES",
                   ifelse(HAS_PBKR, "NO -- inspect CSV", "pbkrtest not installed"))))
cat(sprintf("  All dyad-mean t-tests have df = 18: %s\n",
            ifelse(all(res$dyadmean_df_is_18), "YES", "NO -- inspect CSV")))
cat(sprintf("  Max |t_LMM - t_dyadmean| across cells: %.4g\n",
            max(abs(res$t_satt - res$t_dyadmean))))
sink()
message("Results written to: ", outdir)
