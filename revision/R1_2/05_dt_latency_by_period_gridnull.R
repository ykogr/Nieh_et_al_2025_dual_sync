#!/usr/bin/env Rscript
# =====================================================================
# Step 5 (R1-2): Period-resolved latency (Delta-t) distribution analysis
#
# Input:
#   arg1 = directory (or comma-separated directories, one per condition)
#          containing dt_long_periods_master.csv
#          (output of 02_Micro_Analysis_NLOPM_periods_gridnull.py)
#   arg2 = (optional) output directory (default: first arg1 dir/dt_out)
#
# Analyses (per delta):
#   1) Pair-level summaries per period: n matches, median dt, median |dt|,
#      proportion of dt > 0 (A lagging B), share of |dt| in the visual-response
#      band 150-300 ms (only meaningful for delta >= 0.2).
#   2) LMM: median_abs_dt ~ condition * period + (1 | dyad_id), emmeans (Holm).
#   3) Leader-follower asymmetry: Wilcoxon signed-rank of pair-level median dt
#      vs 0, per condition x period (Holm), + one-sample test of prop(dt>0) vs 0.5.
#   4) Distribution shift across periods: pairwise Kolmogorov-Smirnov tests on the
#      pooled dt distributions between periods within each condition (Holm),
#      reported as a descriptive complement to the pair-level LMM.
#   5) Plots: dt histograms/densities by period x condition.
#
# Usage:
#   Rscript 05_dt_latency_by_period_gridnull.R \
#     out_periods_grid_visible,out_periods_grid_invisible dt_grid_out
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# Exact Wilcoxon signed-rank p-value that handles ties and zeros.
# Priority: exactRankTests::wilcox.exact (shift algorithm, exact with ties/zeros)
#        -> coin::wilcoxsign_test (exact distribution, Pratt zero handling)
#        -> base wilcox.test (normal approximation with continuity correction).
HAS_EXACTRT <- requireNamespace("exactRankTests", quietly = TRUE)
HAS_COIN    <- requireNamespace("coin", quietly = TRUE)
if (!HAS_EXACTRT && !HAS_COIN) {
  message("[NOTE] Neither 'exactRankTests' nor 'coin' installed; ",
          "falling back to normal-approximation wilcox.test. ",
          "For exact p-values: install.packages('exactRankTests')")
}
wilcox_p_exact <- function(x, mu = 0) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  if (HAS_EXACTRT) {
    return(tryCatch(exactRankTests::wilcox.exact(x, mu = mu, exact = TRUE)$p.value,
                    error = function(e) NA_real_))
  }
  if (HAS_COIN) {
    return(tryCatch({
      d <- x - mu
      as.numeric(coin::pvalue(coin::wilcoxsign_test(
        d ~ rep(0, length(d)), distribution = "exact", zero.method = "Pratt")))
    }, error = function(e) NA_real_))
  }
  tryCatch(suppressWarnings(wilcox.test(x, mu = mu)$p.value),
           error = function(e) NA_real_)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 05_dt_latency_by_period_gridnull.R <dt_dir[,dt_dir2]> [outdir]")
dt_dirs <- strsplit(args[1], ",")[[1]]
outdir <- if (length(args) >= 2) args[2] else file.path(dt_dirs[1], "dt_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

paths <- file.path(dt_dirs, "dt_long_periods_master.csv")
for (p in paths) if (!file.exists(p)) stop(paste("not found:", p))
dt <- map_dfr(paths, read_csv, show_col_types = FALSE) %>%
  mutate(
    dyad_id = str_replace(pair_id, "_\\d{2}$", ""),
    condition = factor(condition),
    period = factor(period, levels = unique(period))
  )

sink(file.path(outdir, "dt_latency_results.txt"), split = TRUE)
cat("=====================================================\n")
cat("R1-2 period-resolved Delta-t analysis:", format(Sys.time()), "\n")
cat("Sign convention: dt = t_A - t_B > 0 means A lags behind B\n")
cat("=====================================================\n")

# ---------------- 1) pair-level summaries ----------------
pair_sum <- dt %>%
  group_by(delta, condition, period, dyad_id, pair_id) %>%
  summarise(
    n_matches = n(),
    median_dt = median(dt),
    median_abs_dt = median(abs(dt)),
    prop_pos = mean(dt > 0),
    prop_band_150_300 = mean(abs(dt) >= 0.15 & abs(dt) <= 0.30),
    .groups = "drop"
  )
write_csv(pair_sum, file.path(outdir, "dt_pair_level_summary.csv"))

cat("\n--- pair-level summary (means across pairs) ---\n")
pair_sum %>%
  group_by(delta, condition, period) %>%
  summarise(n_pairs = n(),
            mean_n_matches = mean(n_matches),
            mean_median_dt = mean(median_dt),
            mean_median_abs_dt = mean(median_abs_dt),
            mean_prop_pos = mean(prop_pos),
            mean_prop_band = mean(prop_band_150_300),
            .groups = "drop") %>%
  as.data.frame() %>% print(digits = 3)

min_matches <- 5
cat(sprintf("\nNote: pair x period cells with fewer than %d matches are excluded from tests below.\n", min_matches))
ps <- pair_sum %>% filter(n_matches >= min_matches)

results <- list()

for (delta_val in sort(unique(ps$delta))) {
  pd <- ps %>% filter(delta == delta_val)
  cat(sprintf("\n=====================================================\ndelta = %s\n=====================================================\n", delta_val))

  # ---------------- 2) LMM on median |dt| ----------------
  cat("\n--- LMM: median_abs_dt ~ condition * period + (1|dyad_id) ---\n")
  m <- tryCatch(lmer(median_abs_dt ~ condition * period + (1 | dyad_id), data = pd),
                error = function(e) { cat("MODEL FAILED:", conditionMessage(e), "\n"); NULL })
  if (!is.null(m)) {
    print(summary(m))
    print(anova(m, type = 3))
    ct <- summary(contrast(emmeans(m, ~ condition | period), "pairwise"),
                  adjust = "holm", infer = TRUE)
    cat("\ncondition contrast within period (Holm):\n"); print(ct)
    results[[paste0("lmm_absdt_d", delta_val)]] <- as.data.frame(ct) %>%
      mutate(dv = "median_abs_dt", delta = delta_val)
    pt <- summary(contrast(emmeans(m, ~ period | condition), "pairwise"),
                  adjust = "holm", infer = TRUE)
    cat("\nperiod contrast within condition (Holm):\n"); print(pt)
    results[[paste0("lmm_absdt_periods_d", delta_val)]] <- as.data.frame(pt) %>%
      mutate(dv = "median_abs_dt", delta = delta_val)
  }

  # ---------------- 3) leader-follower asymmetry ----------------
  cat("\n--- Wilcoxon signed-rank: pair-level median dt vs 0 (Holm per delta) ---\n")
  wt <- pd %>%
    group_by(condition, period) %>%
    summarise(n = n(),
              mean_median_dt = mean(median_dt),
              p = wilcox_p_exact(median_dt, mu = 0),
              .groups = "drop") %>%
    mutate(p_holm = p.adjust(p, method = "holm"))
  print(as.data.frame(wt), digits = 3)
  results[[paste0("wilcox_dt_d", delta_val)]] <- wt %>% mutate(dv = "median_dt_vs_0", delta = delta_val)

  cat("\n--- t-test: pair-level prop(dt>0) vs 0.5 (Holm per delta) ---\n")
  bt <- pd %>%
    group_by(condition, period) %>%
    summarise(n = n(),
              mean_prop_pos = mean(prop_pos),
              p = tryCatch(t.test(prop_pos, mu = 0.5)$p.value, error = function(e) NA_real_),
              .groups = "drop") %>%
    mutate(p_holm = p.adjust(p, method = "holm"))
  print(as.data.frame(bt), digits = 3)
  results[[paste0("proppos_d", delta_val)]] <- bt %>% mutate(dv = "prop_pos_vs_0.5", delta = delta_val)

  # ---------------- 4) KS tests between periods (descriptive) ----------------
  cat("\n--- pairwise KS tests on pooled dt between periods, within condition (Holm; descriptive) ---\n")
  dd <- dt %>% filter(delta == delta_val)
  ks_rows <- list()
  for (cond in levels(droplevels(dd$condition))) {
    pers <- levels(droplevels(dd$period))
    combs <- combn(pers, 2, simplify = FALSE)
    for (cb in combs) {
      x <- dd %>% filter(condition == cond, period == cb[1]) %>% pull(dt)
      y <- dd %>% filter(condition == cond, period == cb[2]) %>% pull(dt)
      if (length(x) < 10 || length(y) < 10) next
      kt <- suppressWarnings(ks.test(x, y))
      ks_rows[[paste(cond, cb[1], cb[2])]] <- tibble(
        condition = cond, period_1 = cb[1], period_2 = cb[2],
        n1 = length(x), n2 = length(y), D = unname(kt$statistic), p = kt$p.value)
    }
  }
  if (length(ks_rows) > 0) {
    ks_df <- bind_rows(ks_rows) %>%
      group_by(condition) %>% mutate(p_holm = p.adjust(p, "holm")) %>% ungroup()
    print(as.data.frame(ks_df), digits = 3)
    results[[paste0("ks_d", delta_val)]] <- ks_df %>% mutate(dv = "dt_KS", delta = delta_val)
  }

  # ---------------- 5) plots ----------------
  p1 <- dd %>%
    ggplot(aes(x = dt, fill = condition)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30,
                   position = "identity", alpha = 0.45) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_grid(condition ~ period) +
    labs(title = sprintf("Delta-t distribution by period (delta = %s s)", delta_val),
         x = expression(Delta * t == t[A] - t[B] ~ "(s)"), y = "density") +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(outdir, sprintf("dt_hist_by_period_delta%s.png", delta_val)),
         p1, width = 10, height = 5.5, dpi = 200)
}

sink()
bind_rows(results, .id = "table_id") %>%
  write_csv(file.path(outdir, "dt_tests_all.csv"))
message("Results written to: ", outdir)
