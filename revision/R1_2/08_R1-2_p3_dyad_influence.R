#!/usr/bin/env Rscript
# =====================================================================
# 08_R1-2_p3_dyad_influence.R  (trialfix_kit v4.4)
#   Three analyses and one figure examining whether "the P3 condition difference (visible - invisible)
#    is driven by a few dyads or is a small shift spread across many pairs".
#
#   (1) leave-one-dyad-out: the same model as 04_LMM_condition_by_period_gridnull.R
#         z_S ~ condition * period + (1 | dyad_id)
#       is re-fitted with each dyad removed in turn, and the condition contrast in the target
#       period (default P3) (visible - invisible, emmeans) estimate, SE, df and p are tabulated for all dyads.
#   (2) Sign distribution: for z_S in the target period, at the dyad-mean level (mean of 3 sessions) and
#       at the session level, report the number/proportion of z_S > 0, the exact binomial (sign) test and
#       the Wilcoxon signed-rank test (vs 0) per condition. In addition, a between-condition
#       Wilcoxon rank-sum test on dyad means (non-parametric check of the P3 condition difference).
#   (3) Figure: (a) dot plot of dyad-mean z_S per period x condition (session values faint),
#           (b) leave-one-dyad-out estimates with 95% CI (full-sample estimate as a dashed line).
#
# Input: arg1 = folder containing chew_sync_summary_periods_master.csv (several allowed, comma-separated)
#            e.g. out_periods_grid_visible,out_periods_grid_invisible
#       arg2 = output folder (default dyad_influence_grid)
#       arg3 = target period (default P3)
# Output: dyad_influence_results.txt   human-readable summary (everything)
#       loo_estimates.csv            leave-one-dyad-out, 1 row = dropped dyad
#       dyad_level_zS.csv            dyad-mean z_S in the target period (by condition)
#       sign_tests.csv               table of sign tests / Wilcoxon
#       fig_dyad_influence.pdf/.png  figure (a)(b)
#
# Required packages: lme4, lmerTest, emmeans, dplyr, readr, purrr, stringr, tidyr, ggplot2
# Usage: Rscript 08_R1-2_p3_dyad_influence.R out_periods_grid_visible,out_periods_grid_invisible dyad_influence_grid P3
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(stringr); library(tidyr)
  library(lme4); library(lmerTest); library(emmeans); library(ggplot2)
})

options(width = 200)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 08_R1-2_p3_dyad_influence.R <summary_dir[,summary_dir2]> [outdir] [period]")
summary_dirs <- str_split(args[1], ",")[[1]]
outdir <- if (length(args) >= 2) args[2] else "dyad_influence_grid"
target <- if (length(args) >= 3) args[3] else "P3"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------- load (same conventions as 04) ----------------
master_paths <- file.path(summary_dirs, "chew_sync_summary_periods_master.csv")
for (mp in master_paths) if (!file.exists(mp)) stop(paste("not found:", mp))
d <- map_dfr(master_paths, read_csv, show_col_types = FALSE) %>%
  mutate(session = as.integer(str_extract(pair_id, "\\d{2}$")),
         dyad_id = str_replace(pair_id, "_\\d{2}$", ""),
         condition = factor(condition),
         period = factor(period, levels = unique(period)))
d_ok <- d %>% filter(sufficient == 1, !is.na(z_S))
if (!(target %in% levels(d_ok$period))) stop(paste("period", target, "not in data:", paste(levels(d_ok$period), collapse = ",")))
lv <- levels(d_ok$condition)
if (!all(c("visible", "invisible") %in% lv)) stop(paste("conditions must be visible/invisible, got:", paste(lv, collapse = ",")))
w <- ifelse(lv == "visible", 1, ifelse(lv == "invisible", -1, 0))   # contrast weights: visible - invisible

contrast_target <- function(m) {
  emm <- emmeans(m, ~ condition | period)
  ct <- summary(contrast(emm, method = list("visible - invisible" = w)), infer = TRUE, adjust = "none")
  ct <- as.data.frame(ct)
  ct[ct$period == target, ]
}

sink(file.path(outdir, "dyad_influence_results.txt"), split = TRUE)
cat("=====================================================================\n")
cat("R1-2 dyad influence / distribution analysis:", format(Sys.time()), "\n")
cat("model: z_S ~ condition * period + (1 | dyad_id); contrast = visible - invisible in", target, "(emmeans, unadjusted)\n")
cat("inputs:", paste(master_paths, collapse = " ; "), "\n")
cat("=====================================================================\n")

loo_all <- list(); dyad_all <- list(); sign_all <- list()
for (delta_val in sort(unique(d_ok$delta))) {
  dd <- d_ok %>% filter(delta == delta_val) %>% mutate(dyad_id = factor(dyad_id))
  cat(sprintf("\n############## delta = %s | n rows = %d | n dyads = %d ##############\n", delta_val, nrow(dd), n_distinct(dd$dyad_id)))

  ## ---- (1) full model and leave-one-dyad-out ----
  m_full <- lmer(z_S ~ condition * period + (1 | dyad_id), data = dd)
  full <- contrast_target(m_full)
  cat(sprintf("\n[full sample] %s contrast: est = %.3f, SE = %.3f, df = %.1f, t = %.2f, p = %.4f, 95%% CI [%.3f, %.3f]\n",
              target, full$estimate, full$SE, full$df, full$t.ratio, full$p.value, full$lower.CL, full$upper.CL))
  loo <- map_dfr(levels(dd$dyad_id), function(dy) {
    sub <- dd %>% filter(dyad_id != dy) %>% mutate(dyad_id = droplevels(dyad_id))
    m <- lmer(z_S ~ condition * period + (1 | dyad_id), data = sub)
    ct <- contrast_target(m)
    tibble(delta = delta_val, dropped_dyad = dy,
           dropped_condition = as.character(unique(dd$condition[dd$dyad_id == dy]))[1],
           n_rows = nrow(sub), estimate = ct$estimate, se = ct$SE, df = ct$df, t = ct$t.ratio, p = ct$p.value,
           ci_lo = ct$lower.CL, ci_hi = ct$upper.CL)
  })
  loo_all[[as.character(delta_val)]] <- loo
  cat(sprintf("\n[leave-one-dyad-out, %d refits] %s contrast (visible - invisible)\n", nrow(loo), target))
  print(as.data.frame(loo %>% mutate(across(c(estimate, se, ci_lo, ci_hi), ~ round(.x, 3)), df = round(df, 1), t = round(t, 2), p = round(p, 4))), row.names = FALSE)
  cat(sprintf("  range of estimates: %.3f to %.3f (full sample %.3f); all positive: %s; max p = %.4f; n with p < 0.05: %d / %d\n",
              min(loo$estimate), max(loo$estimate), full$estimate, all(loo$estimate > 0), max(loo$p), sum(loo$p < 0.05), nrow(loo)))
  infl <- loo %>% mutate(change = estimate - full$estimate) %>% arrange(change)
  cat(sprintf("  most influential dyad (largest drop in estimate when removed): %s (%s), est -> %.3f (change %.3f)\n",
              infl$dropped_dyad[1], infl$dropped_condition[1], infl$estimate[1], infl$change[1]))
  cat(sprintf("  largest increase when removed: %s (%s), est -> %.3f (change %+.3f)\n",
              infl$dropped_dyad[nrow(infl)], infl$dropped_condition[nrow(infl)], infl$estimate[nrow(infl)], infl$change[nrow(infl)]))

  ## ---- (2) sign distribution in the target period ----
  tp <- dd %>% filter(period == target)
  dyad_means <- tp %>% group_by(condition, dyad_id) %>%
    summarise(n_sessions = n(), zS_mean = mean(z_S), zS_min = min(z_S), zS_max = max(z_S), .groups = "drop") %>%
    mutate(delta = delta_val, period = target)
  dyad_all[[as.character(delta_val)]] <- dyad_means
  sign_row <- function(x, unit, cond) {
    n <- length(x); k <- sum(x > 0)
    bt <- binom.test(k, n, p = 0.5)
    wt <- suppressWarnings(wilcox.test(x, mu = 0, exact = TRUE))
    tibble(delta = delta_val, period = target, condition = cond, unit = unit, n = n, n_positive = k,
           prop_positive = k / n, mean_zS = mean(x), median_zS = median(x),
           sign_test_p = bt$p.value, wilcoxon_vs0_p = wt$p.value)
  }
  st <- bind_rows(
    map_dfr(lv, ~ sign_row(dyad_means$zS_mean[dyad_means$condition == .x], "dyad mean (3 sessions)", .x)),
    map_dfr(lv, ~ sign_row(tp$z_S[tp$condition == .x], "session", .x)))
  sign_all[[as.character(delta_val)]] <- st
  cat(sprintf("\n[sign distribution in %s]\n", target))
  print(as.data.frame(st %>% mutate(across(c(prop_positive, mean_zS, median_zS), ~ round(.x, 3)), across(c(sign_test_p, wilcoxon_vs0_p), ~ round(.x, 4)))), row.names = FALSE)
  rs <- wilcox.test(zS_mean ~ condition, data = dyad_means, exact = TRUE)
  cat(sprintf("  Wilcoxon rank-sum, dyad means visible vs invisible in %s: W = %.0f, p = %.4f (n = %d vs %d)\n",
              target, rs$statistic, rs$p.value, sum(dyad_means$condition == "visible"), sum(dyad_means$condition == "invisible")))
  cat(sprintf("  dyad means, visible:   %s\n", paste(sprintf("%+.2f", sort(dyad_means$zS_mean[dyad_means$condition == "visible"])), collapse = " ")))
  cat(sprintf("  dyad means, invisible: %s\n", paste(sprintf("%+.2f", sort(dyad_means$zS_mean[dyad_means$condition == "invisible"])), collapse = " ")))

  ## ---- (3) figure (one per delta; the delta used in the paper is 0.1) ----
  dm_all <- dd %>% group_by(delta, condition, period, dyad_id) %>% summarise(zS_mean = mean(z_S), .groups = "drop")
  cond_mean <- dm_all %>% group_by(condition, period) %>% summarise(m = mean(zS_mean), .groups = "drop")
  pa <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    geom_jitter(data = dd, aes(condition, z_S, colour = condition), width = 0.12, height = 0, alpha = 0.25, size = 1.2) +
    geom_jitter(data = dm_all, aes(condition, zS_mean, colour = condition), width = 0.12, height = 0, size = 2.6) +
    geom_crossbar(data = cond_mean, aes(condition, m, ymin = m, ymax = m), width = 0.5, colour = "black", linewidth = 0.4) +
    facet_wrap(~ period, nrow = 1) +
    scale_colour_manual(values = c(visible = "#E69F00", invisible = "#0072B2"), guide = "none") +
    labs(x = NULL, y = expression(italic(z)[S]), title = "(a) Dyad-level synchrony by period (large: dyad means, small: sessions, bar: condition mean)") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(size = 10))
  pb <- ggplot(loo, aes(x = reorder(dropped_dyad, estimate), y = estimate, colour = dropped_condition)) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
    geom_hline(yintercept = full$estimate, colour = "black", linetype = 3) +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.4) +
    scale_colour_manual(values = c(visible = "#E69F00", invisible = "#0072B2"), name = "Dropped dyad") +
    labs(x = "Dropped dyad", y = paste0(target, " contrast (visible - invisible)"),
         title = paste0("(b) Leave-one-dyad-out estimates of the ", target, " difference (dotted: full sample)")) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7), plot.title = element_text(size = 10))
  tag <- gsub("\\.", "", sprintf("%s", delta_val))
  fn <- file.path(outdir, paste0("fig_dyad_influence_delta", tag))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    g <- patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(1, 1.1))
    ggsave(paste0(fn, ".pdf"), g, width = 8, height = 8); ggsave(paste0(fn, ".png"), g, width = 8, height = 8, dpi = 200)
  } else {
    ggsave(paste0(fn, "_a.pdf"), pa, width = 8, height = 3.8); ggsave(paste0(fn, "_a.png"), pa, width = 8, height = 3.8, dpi = 200)
    ggsave(paste0(fn, "_b.pdf"), pb, width = 8, height = 4.2); ggsave(paste0(fn, "_b.png"), pb, width = 8, height = 4.2, dpi = 200)
    cat("  (patchwork not installed: panels (a) and (b) saved as separate files)\n")
  }
  cat("  figure ->", fn, "\n")
}
sink()
write_csv(bind_rows(loo_all), file.path(outdir, "loo_estimates.csv"))
write_csv(bind_rows(dyad_all), file.path(outdir, "dyad_level_zS.csv"))
write_csv(bind_rows(sign_all), file.path(outdir, "sign_tests.csv"))
message("written: ", outdir, "/dyad_influence_results.txt, loo_estimates.csv, dyad_level_zS.csv, sign_tests.csv, fig_dyad_influence_*")
