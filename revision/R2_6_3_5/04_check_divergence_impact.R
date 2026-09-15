#!/usr/bin/env Rscript
# =====================================================================
# 04_check_divergence_impact.R
# Assess whether divergent transitions materially affect beta_state.
# No refitting: works on a saved cmdstanr fit (.rds).
#
# Usage (zsh):
#   Rscript 04_check_divergence_impact.R <fit_rds> <out_csv>
#
# Output:
#   - Console: divergence counts per chain, max |mean shift| across bins,
#     max across-chain spread of bin means
#   - CSV: per-bin comparison of beta_state with all draws vs
#     non-divergent draws only
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04_check_divergence_impact.R <fit_rds> <out_csv>")
}
fit_rds <- args[1]
out_csv <- args[2]

fit <- readRDS(fit_rds)

diag_df <- as_draws_df(fit$sampler_diagnostics())
draws_df <- as_draws_df(fit$draws("beta_state"))
stopifnot(nrow(diag_df) == nrow(draws_df))

div <- diag_df$divergent__ == 1
cat(sprintf("Fit: %s\n", fit_rds))
cat(sprintf("Total draws: %d | divergent: %d (%.2f%%)\n",
            length(div), sum(div), 100 * mean(div)))
cat("Divergences per chain:\n")
print(tapply(div, diag_df$.chain, sum))

# attach divergence flag by .draw index
div_tbl <- tibble(.draw = draws_df$.draw, divergent = div)
long <- draws_df %>%
  pivot_longer(starts_with("beta_state"),
               names_to = "variable", values_to = "value") %>%
  left_join(div_tbl, by = ".draw")

sum_all <- long %>%
  group_by(variable) %>%
  summarise(
    mean_all = mean(value),
    q2.5_all = quantile(value, 0.025),
    q97.5_all = quantile(value, 0.975),
    .groups = "drop"
  )

sum_nodiv <- long %>%
  filter(!divergent) %>%
  group_by(variable) %>%
  summarise(
    mean_nodiv = mean(value),
    q2.5_nodiv = quantile(value, 0.025),
    q97.5_nodiv = quantile(value, 0.975),
    .groups = "drop"
  )

chain_spread <- long %>%
  group_by(variable, .chain) %>%
  summarise(chain_mean = mean(value), .groups = "drop") %>%
  group_by(variable) %>%
  summarise(chain_mean_spread = max(chain_mean) - min(chain_mean),
            .groups = "drop")

res <- sum_all %>%
  left_join(sum_nodiv, by = "variable") %>%
  left_join(chain_spread, by = "variable") %>%
  mutate(
    bin_model = as.integer(gsub("\\D", "", variable)),
    mean_shift = mean_all - mean_nodiv
  ) %>%
  arrange(bin_model)

write_csv(res, out_csv)

cat(sprintf("\nMax |mean shift| (all vs non-divergent): %.4f\n",
            max(abs(res$mean_shift))))
cat(sprintf("Max across-chain spread of bin means:    %.4f\n",
            max(res$chain_mean_spread)))
cat(sprintf("Output written to: %s\n", out_csv))
