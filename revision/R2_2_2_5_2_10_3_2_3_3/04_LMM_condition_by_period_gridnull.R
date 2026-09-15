#!/usr/bin/env Rscript
# =====================================================================
# Step 3 (R1-2): condition x period LMMs on period-restricted synchrony
#
# Inputs:
#   arg1 = directory (or comma-separated directories, one per condition)
#          containing chew_sync_summary_periods_master.csv
#          (output of 02_Micro_Analysis_NLOPM_periods_gridnull.py)
#   arg2 = (optional) directory (or comma-separated directories) containing
#          pseudo_pairs_*.csv (output of 03_pseudo_pairs_gridnull.py)
#   arg3 = (optional) output directory (default: first arg1 dir/lmm_out)
#
# Models:
#   z_S        ~ condition * period + (1 | pair_id)   [per delta]
#   STR        ~ condition * period + (1 | pair_id)   [periods long enough only]
#   Max_episode~ condition * period + (1 | pair_id)
#   emmeans: condition contrast within each period (Holm), and period contrasts
#            within each condition (Holm).
#   Pseudo-pair test (per condition x period):
#   z_S ~ is_real + (1 | a_id) + (1 | b_id)           [crossed random effects]
#
# Usage:
#   Rscript 04_LMM_condition_by_period_gridnull.R \
#     out_periods_grid_visible,out_periods_grid_invisible \
#     out_pseudo_grid_visible,out_pseudo_grid_invisible lmm_periods_grid_out
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 04_LMM_condition_by_period_gridnull.R <summary_dir[,summary_dir2]> [pseudo_dir[,pseudo_dir2]] [outdir]")
summary_dirs <- strsplit(args[1], ",")[[1]]
pseudo_dirs  <- if (length(args) >= 2 && nzchar(args[2])) strsplit(args[2], ",")[[1]] else NA
outdir       <- if (length(args) >= 3) args[3] else file.path(summary_dirs[1], "lmm_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------- load ----------------
master_paths <- file.path(summary_dirs, "chew_sync_summary_periods_master.csv")
for (mp in master_paths) if (!file.exists(mp)) stop(paste("not found:", mp))
d <- map_dfr(master_paths, read_csv, show_col_types = FALSE)

# session handling as in the original notebook: pair_id like 20170706_01_02_01
d <- d %>%
  mutate(
    session = str_extract(pair_id, "\\d{2}$"),
    dyad_id = str_replace(pair_id, "_\\d{2}$", ""),
    session = as.integer(session),
    condition = factor(condition),
    period = factor(period, levels = unique(period))  # keep P1, P2, ... order
  )

message("conditions: ", paste(levels(d$condition), collapse = ", "))
message("periods:    ", paste(levels(d$period), collapse = ", "))

n_excluded <- sum(d$sufficient == 0, na.rm = TRUE)
message(sprintf("excluding %d pair-x-period rows with sufficient == 0 (min event criterion)", n_excluded))
d_ok <- d %>% filter(sufficient == 1, !is.na(z_S))

sink(file.path(outdir, "lmm_results.txt"), split = TRUE)
cat("=====================================================\n")
cat("R1-2 period-resolved LMM analysis:", format(Sys.time()), "\n")
cat("Excluded rows (insufficient events):", n_excluded, "\n")
cat("=====================================================\n")

emm_tables <- list()

fit_one <- function(dat, dv, label, delta_val) {
  cat("\n-----------------------------------------------------\n")
  cat(sprintf("DV = %s | delta = %s | n rows = %d\n", dv, delta_val, nrow(dat)))
  cat("-----------------------------------------------------\n")
  f <- as.formula(paste(dv, "~ condition * period + (1 | dyad_id)"))
  m <- tryCatch(lmer(f, data = dat), error = function(e) { cat("MODEL FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) return(invisible(NULL))
  print(summary(m))
  cat("\n--- Type III ANOVA (Satterthwaite) ---\n")
  print(anova(m, type = 3))

  cat("\n--- condition contrast within each period (Holm) ---\n")
  emm_c <- emmeans(m, ~ condition | period)
  ct <- summary(contrast(emm_c, method = "pairwise"), adjust = "holm", infer = TRUE)
  print(ct)
  emm_tables[[paste(label, "cond_by_period", sep = "_")]] <<-
    as.data.frame(ct) %>% mutate(dv = dv, delta = delta_val)

  cat("\n--- period contrasts within each condition (Holm) ---\n")
  emm_p <- emmeans(m, ~ period | condition)
  pt <- summary(contrast(emm_p, method = "pairwise"), adjust = "holm", infer = TRUE)
  print(pt)
  emm_tables[[paste(label, "period_by_cond", sep = "_")]] <<-
    as.data.frame(pt) %>% mutate(dv = dv, delta = delta_val)
  invisible(m)
}

for (delta_val in sort(unique(d_ok$delta))) {
  dd <- d_ok %>% filter(delta == delta_val)

  ## z_S: primary DV, all periods
  fit_one(dd, "z_S", paste0("zS_d", delta_val), delta_val)

  ## STR / Max_episode: only periods long enough for episode detection (T_seg >= 2*win = 10 s
  ## is enforced upstream; additionally require non-missing values)
  dd_ep <- dd %>% filter(!is.na(STR))
  if (nrow(dd_ep) > 0 && n_distinct(dd_ep$period) >= 2) {
    fit_one(dd_ep, "STR", paste0("STR_d", delta_val), delta_val)
    fit_one(dd_ep, "Max_episode", paste0("MaxEp_d", delta_val), delta_val)
  } else {
    cat(sprintf("\n[skip] STR/Max_episode at delta=%s: not enough episode data\n", delta_val))
  }

  ## per-period one-sample check: is z_S > 0 at the group level within each period x condition?
  cat("\n--- group-level one-sample t: z_S vs 0, per condition x period (Holm within delta) ---\n")
  os <- dd %>%
    group_by(condition, period) %>%
    summarise(mean_zS = mean(z_S), sd = sd(z_S), n = n(),
              t = mean(z_S) / (sd(z_S) / sqrt(n())),
              p = 2 * pt(abs(mean(z_S) / (sd(z_S) / sqrt(n()))), df = n() - 1, lower.tail = FALSE),
              .groups = "drop") %>%
    mutate(p_holm = p.adjust(p, method = "holm"))
  print(as.data.frame(os), digits = 3)
  emm_tables[[paste0("onesample_d", delta_val)]] <- os %>% mutate(dv = "z_S_vs_0", delta = delta_val)
}

# ---------------- pseudo-pair comparison ----------------
if (!all(is.na(pseudo_dirs))) {
  pp_files <- list.files(pseudo_dirs, pattern = "^pseudo_pairs_.*\\.csv$", full.names = TRUE)
  if (length(pp_files) == 0) {
    cat("\n[skip] no pseudo_pairs_*.csv in", paste(pseudo_dirs, collapse = ", "), "\n")
  } else {
    pp <- map_dfr(pp_files, read_csv, show_col_types = FALSE) %>%
      filter(sufficient == 1, !is.na(z_S)) %>%
      mutate(is_real = factor(is_real, levels = c(0, 1), labels = c("pseudo", "real")),
             period = factor(period, levels = levels(d$period)))
    cat("\n=====================================================\n")
    cat("Pseudo-pair control: z_S ~ is_real + (1|a_id) + (1|b_id)\n")
    cat("=====================================================\n")
    for (delta_val in sort(unique(pp$delta))) {
      for (cond in unique(pp$condition)) {
        for (per in levels(droplevels(pp$period))) {
          sub <- pp %>% filter(delta == delta_val, condition == cond, period == per)
          if (nrow(sub) < 8 || n_distinct(sub$is_real) < 2) next
          cat(sprintf("\n--- delta=%s | condition=%s | period=%s (n=%d, real=%d) ---\n",
                      delta_val, cond, per, nrow(sub), sum(sub$is_real == "real")))
          m <- tryCatch(lmer(z_S ~ is_real + (1 | a_id) + (1 | b_id), data = sub),
                        error = function(e) { cat("MODEL FAILED:", conditionMessage(e), "\n"); NULL })
          if (!is.null(m)) {
            print(coef(summary(m)))
            emm_tables[[paste("pseudo", delta_val, cond, per, sep = "_")]] <-
              as.data.frame(coef(summary(m))) %>%
              rownames_to_column("term") %>%
              mutate(dv = "z_S_real_vs_pseudo", delta = delta_val,
                     condition = cond, period = per)
          }
          # descriptive: where does the real-pair mean sit in the pseudo distribution?
          real_mean <- mean(sub$z_S[sub$is_real == "real"])
          pseudo_z  <- sub$z_S[sub$is_real == "pseudo"]
          emp_p <- (sum(pseudo_z >= real_mean) + 1) / (length(pseudo_z) + 1)
          cat(sprintf("real mean z_S = %.3f | pseudo mean = %.3f | empirical P(pseudo >= real mean) = %.4f\n",
                      real_mean, mean(pseudo_z), emp_p))
        }
      }
    }
  }
}

sink()

# ---------------- export tables ----------------
if (length(emm_tables) > 0) {
  bind_rows(emm_tables, .id = "table_id") %>%
    write_csv(file.path(outdir, "emmeans_and_tests_all.csv"))
}
message("Results written to: ", outdir)
