#!/usr/bin/env Rscript
# =====================================================================
# 06_condition_difference_posterior.R
# R2-6 / R3-5 (macro scale): posterior of the between-condition difference in
# mean rhythmic coupling, replacing the "p < 0.001" of the previous submission.
#
#   Delta = mean_b beta_invisible(b) - mean_b beta_visible(b),   b = 1..13 bins
#
# No refitting.  Works on the saved cmdstanr fits written by
# 03_sensitivity_drift_beta.R (fit_<mode>_<condition>.rds).  The two conditions
# are fitted independently, so a Monte-Carlo sample of Delta is obtained by
# pairing the posterior draws of the two fits draw-by-draw (same chains x
# iterations, same seed).  A random re-pairing (fixed seed) is reported as a
# check that the pairing order does not matter.
#
# Usage (zsh, from R2_6_3_5/):
#   Rscript 06_condition_difference_posterior.R <fit_dir> [modes] [out_dir]
#     fit_dir : folder containing fit_<mode>_visible_pair.rds and
#               fit_<mode>_invisible_pair.rds            (e.g. sensitivity_out)
#     modes   : comma-separated subset of raw,detrend,diff (default: raw,detrend,diff)
#     out_dir : output folder (default: <fit_dir>)
#
# Outputs (in out_dir):
#   condition_difference_<mode>.txt   : human-readable log (one per mode)
#   condition_difference_all.csv      : one row per mode (overall Delta)
#   condition_difference_bins_all.csv : one row per mode x bin (bin-wise Delta_b)
#
# Reported quantities per mode:
#   mean_beta_invisible, mean_beta_visible : posterior means of the across-bin means
#   delta_mean, delta_median, delta_sd, delta_q2.5, delta_q97.5
#   p_delta_gt0        : posterior probability P(Delta > 0)
#   n_draws_delta_le0  : number of draws with Delta <= 0 (0 -> report P > 1 - 1/n_draws)
#   p_delta_gt0_repair : same after random re-pairing of draws (seed 20260910)
#   n_bins_delta_excl_zero, bins_delta_excl_zero : bin-wise 95% CrI of Delta_b excluding 0
# =====================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript 06_condition_difference_posterior.R <fit_dir> [modes] [out_dir]")
}
fit_dir <- args[1]
modes   <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("raw", "detrend", "diff")
out_dir <- if (length(args) >= 3) args[3] else fit_dir
stopifnot(all(modes %in% c("raw", "detrend", "diff")))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

REPAIR_SEED <- 20260910
N_BIN_EXPECTED <- 13

get_beta_draws <- function(fit_rds) {
  stopifnot(file.exists(fit_rds))
  fit <- readRDS(fit_rds)
  m <- as_draws_matrix(fit$draws("beta_state"))   # draws x bins
  # order columns by bin index (beta_state[1] .. beta_state[13])
  idx <- as.integer(gsub("\\D", "", colnames(m)))
  m <- m[, order(idx), drop = FALSE]
  if (ncol(m) != N_BIN_EXPECTED) {
    warning(sprintf("%s: %d beta_state columns (expected %d)", fit_rds, ncol(m), N_BIN_EXPECTED))
  }
  m
}

fmt <- function(x, d = 6) formatC(x, digits = d, format = "f")

all_rows <- list()
bin_rows <- list()

for (mode in modes) {
  f_inv <- file.path(fit_dir, sprintf("fit_%s_invisible_pair.rds", mode))
  f_vis <- file.path(fit_dir, sprintf("fit_%s_visible_pair.rds", mode))
  message(sprintf("==== mode = %s ====", mode))
  B_inv <- get_beta_draws(f_inv)
  B_vis <- get_beta_draws(f_vis)
  n_inv <- nrow(B_inv); n_vis <- nrow(B_vis)
  n_use <- min(n_inv, n_vis)
  if (n_inv != n_vis) {
    warning(sprintf("[%s] draw counts differ (invisible %d, visible %d); first %d draws used",
                    mode, n_inv, n_vis, n_use))
    B_inv <- B_inv[seq_len(n_use), , drop = FALSE]
    B_vis <- B_vis[seq_len(n_use), , drop = FALSE]
  }
  n_bin <- ncol(B_inv)

  # ---- across-bin mean per draw, then the difference ---------------------------------------
  m_inv  <- rowMeans(B_inv)
  m_vis  <- rowMeans(B_vis)
  delta  <- m_inv - m_vis
  q      <- quantile(delta, c(0.025, 0.5, 0.975), names = FALSE)
  n_le0  <- sum(delta <= 0)
  p_gt0  <- mean(delta > 0)

  # ---- random re-pairing check (fixed seed) --------------------------------------------------
  set.seed(REPAIR_SEED)
  perm      <- sample.int(n_use)
  delta_rp  <- m_inv - m_vis[perm]
  q_rp      <- quantile(delta_rp, c(0.025, 0.975), names = FALSE)
  p_gt0_rp  <- mean(delta_rp > 0)

  # ---- bin-wise difference -------------------------------------------------------------------
  D_b   <- B_inv - B_vis
  b_mean <- colMeans(D_b)
  b_lo   <- apply(D_b, 2, quantile, probs = 0.025, names = FALSE)
  b_hi   <- apply(D_b, 2, quantile, probs = 0.975, names = FALSE)
  b_pgt0 <- colMeans(D_b > 0)
  excl   <- (b_lo > 0) | (b_hi < 0)
  bin_model <- seq_len(n_bin)
  bin_data  <- bin_model + 5
  window_end_s <- bin_data * 10

  bin_rows[[mode]] <- data.frame(
    mode = mode, bin_model = bin_model, bin_data = bin_data, window_end_s = window_end_s,
    mean_beta_invisible_bin = colMeans(B_inv), mean_beta_visible_bin = colMeans(B_vis),
    delta_mean = b_mean, delta_q2.5 = b_lo, delta_q97.5 = b_hi, p_delta_gt0 = b_pgt0,
    cri_excl_zero = as.integer(excl), stringsAsFactors = FALSE)

  all_rows[[mode]] <- data.frame(
    mode = mode, n_draws = n_use, n_bins = n_bin,
    fit_invisible = basename(f_inv), fit_visible = basename(f_vis),
    mean_beta_invisible = mean(m_inv), mean_beta_visible = mean(m_vis),
    delta_mean = mean(delta), delta_median = q[2], delta_sd = sd(delta),
    delta_q2.5 = q[1], delta_q97.5 = q[3],
    p_delta_gt0 = p_gt0, n_draws_delta_le0 = n_le0,
    p_delta_gt0_repair = p_gt0_rp, delta_q2.5_repair = q_rp[1], delta_q97.5_repair = q_rp[2],
    n_bins_delta_excl_zero = sum(excl),
    bins_delta_excl_zero = paste(bin_model[excl], collapse = ";"),
    windows_delta_excl_zero_s = paste(window_end_s[excl], collapse = ";"),
    reported_rounded = sprintf("%.2f [%.2f, %.2f]", mean(delta), q[1], q[3]),
    stringsAsFactors = FALSE)

  # ---- human-readable log ----------------------------------------------------------------------
  txt <- file.path(out_dir, sprintf("condition_difference_%s.txt", mode))
  con <- file(txt, open = "wt")
  w <- function(...) writeLines(sprintf(...), con)
  w("R2-6 / R3-5 macro-scale condition difference: Delta = mean_b beta_invisible(b) - mean_b beta_visible(b)")
  w("run: %s", format(Sys.time()))
  w("mode = %s", mode)
  w("fits: %s, %s", f_inv, f_vis)
  w("draws per fit used: %d (invisible %d, visible %d); bins: %d; fits independent, paired by draw index",
    n_use, n_inv, n_vis, n_bin)
  w("")
  w("mean_beta invisible = %s   mean_beta visible = %s", fmt(mean(m_inv)), fmt(mean(m_vis)))
  w("Delta mean = %s   median = %s   sd = %s", fmt(mean(delta)), fmt(q[2]), fmt(sd(delta)))
  w("Delta 95%% CrI = [%s, %s]", fmt(q[1]), fmt(q[3]))
  w("P(Delta > 0) = %s   (draws with Delta <= 0: %d of %d)", fmt(p_gt0), n_le0, n_use)
  if (n_le0 == 0) w("  -> report as P(Delta > 0) > %s (no draw <= 0 among %d)", fmt(1 - 1 / n_use, 5), n_use)
  w("check, random re-pairing (seed %d): P(Delta > 0) = %s, 95%% CrI = [%s, %s]",
    REPAIR_SEED, fmt(p_gt0_rp), fmt(q_rp[1]), fmt(q_rp[2]))
  w("reported (2 decimals): Delta = %s", all_rows[[mode]]$reported_rounded)
  w("")
  w("bin-wise difference Delta_b = beta_invisible(b) - beta_visible(b):")
  w("%9s %8s %12s %10s %10s %10s %10s %9s", "bin_model", "bin_data", "window_end_s",
    "mean", "q2.5", "q97.5", "P(>0)", "excl_zero")
  for (b in seq_len(n_bin)) {
    w("%9d %8d %12d %10.4f %10.4f %10.4f %10.4f %9d", bin_model[b], bin_data[b], window_end_s[b],
      b_mean[b], b_lo[b], b_hi[b], b_pgt0[b], as.integer(excl[b]))
  }
  w("bins with 95%% CrI of Delta_b excluding zero: %d of %d (bins: %s; windows ending: %s s)",
    sum(excl), n_bin, paste(bin_model[excl], collapse = ";"), paste(window_end_s[excl], collapse = ";"))
  close(con)
  cat(readLines(txt), sep = "\n"); cat("\n")
}

res  <- do.call(rbind, all_rows)
bins <- do.call(rbind, bin_rows)
write.csv(res,  file.path(out_dir, "condition_difference_all.csv"),      row.names = FALSE)
write.csv(bins, file.path(out_dir, "condition_difference_bins_all.csv"), row.names = FALSE)
cat("Wrote", file.path(out_dir, "condition_difference_all.csv"), "and condition_difference_bins_all.csv\n")
