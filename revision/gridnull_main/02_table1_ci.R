#!/usr/bin/env Rscript
# =====================================================================
# 02_table1_ci.R  (trialfix_kit v4.3.4)
#   Compute 95% confidence intervals for the condition difference (visible - invisible) in Table 1.
#
#   The same model as 01_R2-8_df_verification.R
#       dv ~ condition + (1 | pair_id)      (60 sessions = 20 dyads x 3)
#   is fitted for each DV (z_S, STR, Max_episode) x each delta, and a Wald-type interval
#   using the Kenward-Roger degrees of freedom df_KR and standard error SE_KR
#       estimate +/- qt(0.975, df_KR) * SE_KR
#   is reported. Because it is based on the same quantities as the t and p values
#   (Kenward-Roger) of Table 1, t, p and the interval are mutually consistent within the table.
#   For reference, the profile-likelihood interval (lme4::confint, method = "profile") is also
#   given (it may fail in cells where a variance component is 0; NA in that case).
#
# Input:  arg1 = folder containing {delta}_chew_sync_summary_{visible|invisible}_pair_Human.csv
#               (same as 04_Micro_Analysis_LMM.Rmd; e.g. ../Micro_Analysis_LMM)
#        arg2 = output folder (default df_verification_out)
# Output: table1_ci.csv  1 row = DV x delta. estimate, se_kr, df_kr, t_kr, p_kr,
#                       t_crit, ci_lo, ci_hi, ci_profile_lo, ci_profile_hi, ci_text
#        table1_ci.txt  human-readable summary (with "[lo, hi]" ready to paste into manuscript Table 1)
#
# Required packages: lme4, lmerTest, pbkrtest, dplyr, readr, purrr, stringr, tibble
# Usage: Rscript 02_table1_ci.R ../Micro_Analysis_LMM df_verification_out
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(stringr); library(tibble)
  library(lme4); library(lmerTest)
})
if (!requireNamespace("pbkrtest", quietly = TRUE))
  stop("pbkrtest is required: install.packages('pbkrtest')")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 02_table1_ci.R <indir> [outdir]")
indir  <- args[1]
outdir <- if (length(args) >= 2) args[2] else "df_verification_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- load data (same conventions as 01_R2-8) ----
files <- list.files(indir,
                    pattern = "^\\d+\\.\\d+_chew_sync_summary_(visible|invisible)_pair_Human\\.csv$",
                    full.names = TRUE)
if (length(files) == 0) stop(paste("summary CSV not found:", indir))
d <- map_dfr(files, function(file) {
  fname <- basename(file)
  df <- read_csv(file, show_col_types = FALSE)
  df$delta <- as.numeric(str_extract(fname, "^\\d+\\.\\d+"))
  df$condition <- str_remove(str_extract(fname, "(visible|invisible)_pair_Human"), "_Human")
  df
}) %>%
  mutate(session = as.integer(str_extract(pair_id, "\\d{2}$")),
         pair_id = str_replace(pair_id, "_\\d{2}$", ""))

dvs <- c("z_S", "STR", "Max_episode")
rows <- list()
for (dv in dvs) {
  for (delta_val in sort(unique(d$delta))) {
    dd <- d %>% filter(delta == delta_val)
    m  <- lmer(reformulate("condition + (1 | pair_id)", response = dv), data = dd)
    ck <- coef(summary(m, ddf = "Kenward-Roger"))
    i  <- grep("^condition", rownames(ck))
    est <- ck[i, "Estimate"]; se <- ck[i, "Std. Error"]; dfk <- ck[i, "df"]
    tk  <- ck[i, "t value"];  pk <- ck[i, "Pr(>|t|)"]
    tcrit <- qt(0.975, dfk)
    prof <- tryCatch({
      ci <- suppressMessages(confint(m, parm = rownames(ck)[i], method = "profile", quiet = TRUE))
      c(ci[1, 1], ci[1, 2])
    }, error = function(e) c(NA_real_, NA_real_), warning = function(w) c(NA_real_, NA_real_))
    rows[[paste(dv, delta_val)]] <- tibble(
      dv = dv, delta = delta_val, n_sessions = nrow(dd), n_dyads = n_distinct(dd$pair_id),
      coef = rownames(ck)[i], estimate = est, se_kr = se, df_kr = dfk, t_kr = tk, p_kr = pk,
      t_crit = tcrit, ci_lo = est - tcrit * se, ci_hi = est + tcrit * se,
      ci_profile_lo = prof[1], ci_profile_hi = prof[2]
    )
  }
}
res <- bind_rows(rows) %>%
  mutate(ci_text = sprintf("[%s, %s]", formatC(ci_lo, format = "f", digits = 3),
                           formatC(ci_hi, format = "f", digits = 3)))
write_csv(res, file.path(outdir, "table1_ci.csv"))

sink(file.path(outdir, "table1_ci.txt"), split = TRUE)
cat("=====================================================================\n")
cat("Table 1: 95% CI of the condition effect (visible - invisible):", format(Sys.time()), "\n")
cat("Model: dv ~ condition + (1 | pair_id); CI = estimate +/- qt(0.975, df_KR) * SE_KR\n")
cat("(Kenward-Roger df and SE, the same quantities as the t and p columns of Table 1)\n")
cat("=====================================================================\n\n")
res %>%
  transmute(dv, delta, estimate = round(estimate, 3), se_kr = round(se_kr, 3), df_kr = round(df_kr, 1),
            t_kr = round(t_kr, 3), p_kr = round(p_kr, 3), t_crit = round(t_crit, 3), ci_95 = ci_text,
            profile_lo = round(ci_profile_lo, 3), profile_hi = round(ci_profile_hi, 3)) %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nNotes:\n")
cat("  ci_95     : Wald-type interval with Kenward-Roger df/SE -> paste into Table 1 (95% CI column).\n")
cat("  profile_* : profile-likelihood interval (lme4::confint) for reference; NA where the profile could not be computed.\n")
sink()
message("written: ", file.path(outdir, "table1_ci.csv"), " and table1_ci.txt")
