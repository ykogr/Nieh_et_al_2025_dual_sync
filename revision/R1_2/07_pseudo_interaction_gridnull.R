#!/usr/bin/env Rscript
# =====================================================================
# Step 7 (R1-2): direct test of the condition difference in the
#                pseudo-pair effect (difference of differences)
#
# Motivation:
#   04 tests real vs pseudo SEPARATELY per condition x period. A reviewer
#   may ask: "did you test whether the real-pseudo difference differs
#   between visible and invisible?" This script answers that with the
#   interaction model, pooling both conditions:
#
#     z_S ~ is_real * condition + (1 | a_id) + (1 | b_id)
#
#   fitted per delta x period. The is_real:condition interaction (and the
#   emmeans interaction contrast) is the direct test of
#     (real - pseudo | visible) - (real - pseudo | invisible).
#
# Inputs:
#   arg1 = directory (or comma-separated directories, one per condition)
#          containing pseudo_pairs_*.csv (output of 03_pseudo_pairs_gridnull.py)
#   arg2 = (optional) output directory   (default: arg1/interaction_out)
#   arg3 = (optional) dyad prefix to exclude, e.g. 20170710_07_08
#          (rows are dropped if a_id OR b_id starts with this prefix)
#
# Outputs:
#   <outdir>/pseudo_interaction_results.txt
#   <outdir>/pseudo_interaction_tests.csv
#
# Usage:
#   Rscript 07_pseudo_interaction_gridnull.R \
#     out_pseudo_grid_visible,out_pseudo_grid_invisible pseudo_interaction_grid_out
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 07_pseudo_interaction_gridnull.R <pseudo_dir[,pseudo_dir2]> [outdir] [exclude_dyad]")
pseudo_dirs  <- strsplit(args[1], ",")[[1]]
outdir       <- if (length(args) >= 2 && nzchar(args[2])) args[2] else file.path(pseudo_dirs[1], "interaction_out")
exclude_dyad <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NA
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pp_files <- list.files(pseudo_dirs, pattern = "^pseudo_pairs_.*\\.csv$", full.names = TRUE)
if (length(pp_files) == 0) stop(paste("no pseudo_pairs_*.csv in", paste(pseudo_dirs, collapse = ", ")))

pp <- map_dfr(pp_files, read_csv, show_col_types = FALSE) %>%
  filter(sufficient == 1, !is.na(z_S))

if (!is.na(exclude_dyad)) {
  n0 <- nrow(pp)
  pp <- pp %>% filter(!startsWith(a_id, exclude_dyad), !startsWith(b_id, exclude_dyad))
  message(sprintf("excluded dyad prefix '%s': %d -> %d rows", exclude_dyad, n0, nrow(pp)))
}

if (n_distinct(pp$condition) < 2) stop("need BOTH conditions in pseudo_dir for the interaction test")

# order periods by their start time, keep is_real labelled
per_levels <- pp %>% distinct(period, t0_rel) %>% arrange(t0_rel) %>% pull(period)
pp <- pp %>%
  mutate(is_real   = factor(is_real, levels = c(0, 1), labels = c("pseudo", "real")),
         condition = factor(condition),
         period    = factor(period, levels = per_levels))

sink(file.path(outdir, "pseudo_interaction_results.txt"), split = TRUE)
cat("=====================================================\n")
cat("R1-2 pseudo-pair interaction test:", format(Sys.time()), "\n")
cat("model: z_S ~ is_real * condition + (1|a_id) + (1|b_id)\n")
cat("files:", paste(basename(pp_files), collapse = ", "), "\n")
if (!is.na(exclude_dyad)) cat("excluded dyad prefix:", exclude_dyad, "\n")
cat("conditions:", paste(levels(pp$condition), collapse = ", "), "\n")
cat("periods:   ", paste(levels(pp$period), collapse = ", "), "\n")
cat("=====================================================\n")

res_rows <- list()

for (delta_val in sort(unique(pp$delta))) {
  for (per in levels(droplevels(pp$period))) {
    sub <- pp %>% filter(delta == delta_val, period == per) %>% droplevels()
    if (nrow(sub) < 16 || n_distinct(sub$is_real) < 2 || n_distinct(sub$condition) < 2) next
    cat(sprintf("\n--- delta=%s | period=%s (n=%d, real=%d) ---\n",
                delta_val, per, nrow(sub), sum(sub$is_real == "real")))

    m <- tryCatch(lmer(z_S ~ is_real * condition + (1 | a_id) + (1 | b_id), data = sub),
                  error = function(e) { cat("MODEL FAILED:", conditionMessage(e), "\n"); NULL })
    if (is.null(m)) next

    cat("\nType III ANOVA (Satterthwaite):\n")
    print(anova(m, type = 3))

    # real - pseudo within each condition (KR df)
    emm <- emmeans(m, ~ is_real | condition, lmer.df = "kenward-roger")
    rp  <- contrast(emm, method = list(real_minus_pseudo = c(-1, 1)))
    cat("\nreal - pseudo within each condition:\n")
    print(summary(rp, infer = TRUE))

    # difference of differences: (real-pseudo | cond2) - (real-pseudo | cond1)
    dd <- contrast(emmeans(m, ~ is_real * condition, lmer.df = "kenward-roger"),
                   interaction = c("revpairwise", "revpairwise"))
    cat("\ninteraction contrast (difference of the real-pseudo effect between conditions):\n")
    dd_sum <- summary(dd, infer = TRUE)
    print(dd_sum)

    res_rows[[paste(delta_val, per, sep = "_")]] <-
      as.data.frame(dd_sum) %>%
      mutate(delta = delta_val, period = per,
             n = nrow(sub), n_real = sum(sub$is_real == "real"))
  }
}

# Holm across periods within each delta for the interaction contrast
if (length(res_rows) > 0) {
  res <- bind_rows(res_rows) %>%
    group_by(delta) %>%
    mutate(p_holm_within_delta = p.adjust(p.value, method = "holm")) %>%
    ungroup()
  cat("\n=====================================================\n")
  cat("Summary of interaction contrasts (Holm across periods within delta):\n")
  cat("positive estimate = (real-pseudo) larger in the alphabetically LATER condition\n")
  cat("(with conditions invisible/visible: positive = larger dyad-specific effect in VISIBLE)\n")
  cat("=====================================================\n")
  print(as.data.frame(res), digits = 4)
  write_csv(res, file.path(outdir, "pseudo_interaction_tests.csv"))
}

sink()
message("Results written to: ", outdir)
