#!/usr/bin/env Rscript
# 05_reported_values.R
# Aggregate bin-wise beta(t) posterior summaries into the values reported in the
# manuscript / response letter, and save them as a CSV so that every reported
# number has a direct one-to-one source in a log file.
#
# Reported values per (mode, condition):
#   - mean_beta        : mean over 13 bins of the bin-wise posterior means
#   - avg_cri_low/high : mean over 13 bins of the bin-wise 95% CrI bounds (q2.5, q97.5)
#   - n_bins_cri_excl_zero, bins_excl_zero (bin_model), windows_excl_zero_s (window end)
#   - max_rhat over the 13 beta_state parameters
#   - reported_rounded : the exact string as reported (2 decimals)
#
# Usage:
#   Rscript 05_reported_values.R <out_csv> <input1> [<input2> ...]
#
# Each <input> is either:
#   (a) a beta_summary_all.csv produced by 03_sensitivity_drift_beta.R
#       (columns: variable, mean, q2.5, q97.5, rhat, bin_model, window_end_s, mode, condition)
#       -> may contain several mode x condition blocks; all are aggregated.
#   (b) a div_check CSV produced by 04_check_divergence_impact.R, which has no
#       mode/condition columns. For these, append labels to the path:
#           <path>:<mode>:<condition>
#       e.g.  div_check_visible.csv:raw:visible_pair
#       (columns used: variable, mean_all, q2.5_all, q97.5_all, bin_model)
#
# Example:
#   Rscript 05_reported_values.R sensitivity_out/reported_values.csv \
#       sensitivity_out/beta_summary_all.csv \
#       div_check_visible.csv:raw:visible_pair \
#       div_check_invisible.csv:raw:invisible_pair

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 05_reported_values.R <out_csv> <input1> [<input2> ...]")
}
out_csv <- args[1]
inputs  <- args[-1]

read_one <- function(spec) {
  parts <- strsplit(spec, ":", fixed = TRUE)[[1]]
  if (length(parts) == 3) {
    path <- parts[1]; mode <- parts[2]; cond <- parts[3]
  } else {
    path <- spec; mode <- NA; cond <- NA
  }
  stopifnot(file.exists(path))
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (all(c("mean", "q2.5", "q97.5", "mode", "condition") %in% names(d))) {
    # beta_summary_all.csv format
    out <- data.frame(
      mode = d$mode, condition = d$condition,
      bin_model = d$bin_model,
      window_end_s = if ("window_end_s" %in% names(d)) d$window_end_s else d$bin_model * 10 + 50,
      mean = d$mean, q_low = d$`q2.5`, q_high = d$`q97.5`,
      rhat = if ("rhat" %in% names(d)) d$rhat else NA,
      source_file = basename(path),
      stringsAsFactors = FALSE
    )
  } else if (all(c("mean_all", "q2.5_all", "q97.5_all") %in% names(d))) {
    # div_check format: needs explicit labels
    if (is.na(mode) || is.na(cond)) {
      stop(sprintf("div_check file '%s' needs labels: %s:<mode>:<condition>", path, path))
    }
    out <- data.frame(
      mode = mode, condition = cond,
      bin_model = d$bin_model,
      window_end_s = d$bin_model * 10 + 50,
      mean = d$mean_all, q_low = d$`q2.5_all`, q_high = d$`q97.5_all`,
      rhat = NA,
      source_file = basename(path),
      stringsAsFactors = FALSE
    )
  } else {
    stop(sprintf("Unrecognized column layout in '%s'", path))
  }
  out
}

d <- do.call(rbind, lapply(inputs, read_one))

groups <- unique(d[, c("mode", "condition")])
res <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
  g <- d[d$mode == groups$mode[i] & d$condition == groups$condition[i], ]
  g <- g[order(g$bin_model), ]
  excl <- g$q_low > 0 | g$q_high < 0
  data.frame(
    mode = groups$mode[i],
    condition = groups$condition[i],
    n_bins = nrow(g),
    mean_beta = mean(g$mean),
    avg_cri_low = mean(g$q_low),
    avg_cri_high = mean(g$q_high),
    n_bins_cri_excl_zero = sum(excl),
    bins_excl_zero = paste(g$bin_model[excl], collapse = ";"),
    windows_excl_zero_s = paste(g$window_end_s[excl], collapse = ";"),
    max_rhat_beta = if (all(is.na(g$rhat))) NA else max(g$rhat, na.rm = TRUE),
    reported_rounded = sprintf("%.2f [%.2f, %.2f]",
                               mean(g$mean), mean(g$q_low), mean(g$q_high)),
    source_file = paste(unique(g$source_file), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n\n")
print(res, row.names = FALSE)
