#!/usr/bin/env Rscript
# =====================================================================
# 07_rhat_main.R  (co-author item 10)
# Report the observed convergence diagnostics of the MAIN macro-scale fit
# (mode "raw" = the pipeline whose beta(t) is shown in Figure 3), so that the
# STAR Methods sentence "all parameters showing R-hat < 1.1" can carry the
# actually observed maximum.  No refitting: works on the saved cmdstanr fits.
#
# Usage (zsh):
#   Rscript 07_rhat_main.R <fit_dir> [mode] [out_dir]
#     fit_dir : folder containing fit_<mode>_visible_pair.rds and
#               fit_<mode>_invisible_pair.rds          (e.g. sensitivity_out)
#     mode    : raw (default) | detrend | diff
#     out_dir : where to write (default: fit_dir)
#
# Outputs (out_dir):
#   rhat_<mode>.txt          : human-readable log (numbers for the manuscript come from here)
#   rhat_<mode>_params.csv   : one row per parameter and condition (rhat, ess_bulk, ess_tail)
#   rhat_<mode>_summary.csv  : one row per condition (max rhat over all parameters, over
#                              beta_state only, over the non-imputation parameters, divergences)
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 07_rhat_main.R <fit_dir> [mode] [out_dir]")
fit_dir <- args[1]
mode    <- if (length(args) >= 2) args[2] else "raw"
out_dir <- if (length(args) >= 3) args[3] else fit_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(out_dir, sprintf("rhat_%s.txt", mode))
cat("", file = log_file)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = log_file, append = TRUE) }

logf("=====================================================================")
logf("07_rhat_main.R : convergence of the %s-mode macro fits : %s", mode, format(Sys.time()))
logf("fit_dir = %s", fit_dir)
logf("=====================================================================")

all_params <- list(); summ <- list()
for (cond in c("visible_pair", "invisible_pair")) {
  f <- file.path(fit_dir, sprintf("fit_%s_%s.rds", mode, cond))
  if (!file.exists(f)) stop("fit not found: ", f, "\n  -> run 03_sensitivity_drift_beta.R with mode '", mode, "' first (it saves fit_<mode>_<condition>.rds)")
  fit <- readRDS(f)
  # v2.2: summarise only the sampled parameters, the transformed parameters (beta_state) and lp__.
  # The generated quantities y_full_out / y_rep / log_lik have 3 x N elements (N = number of windows) and
  # made fit$summary(NULL, ...) run for tens of minutes; they are deterministic functions of the parameters
  # and not part of the convergence statement ("all parameters").
  GQ_EXCLUDE <- c("y_full_out", "y_rep", "log_lik")
  meta0 <- tryCatch(fit$metadata(), error = function(e) NULL)
  vars  <- if (!is.null(meta0)) setdiff(meta0$stan_variables, GQ_EXCLUDE) else NULL
  if (is.null(vars)) stop("cannot read fit$metadata()$stan_variables from ", f)
  logf("[%s] variables summarised: %s", cond, paste(vars, collapse = ", "))
  logf("[%s] excluded generated quantities: %s", cond, paste(GQ_EXCLUDE, collapse = ", "))
  t0 <- Sys.time()
  s <- fit$summary(variables = vars, rhat, ess_bulk, ess_tail) %>%
    mutate(condition = cond,
           block = str_replace(variable, "\\[.*$", ""),
           is_imputed = block == "y_mis",
           is_lp = variable == "lp__")
  all_params[[cond]] <- s

  ds <- fit$diagnostic_summary(quiet = TRUE)
  meta <- tryCatch(fit$metadata(), error = function(e) NULL)
  n_chain <- if (!is.null(meta)) length(meta$id) else NA_integer_
  n_iter  <- if (!is.null(meta)) meta$iter_sampling else NA_integer_
  n_warm  <- if (!is.null(meta)) meta$iter_warmup else NA_integer_

  s_model <- s %>% filter(!is_imputed, !is_lp)
  s_beta  <- s %>% filter(block == "beta_state")
  i_all <- which.max(s$rhat); i_model <- which.max(s_model$rhat); i_beta <- which.max(s_beta$rhat)

  logf("")
  logf("[%s]  file = %s", cond, basename(f))
  logf("  chains = %s | iter_sampling = %s | iter_warmup = %s", n_chain, n_iter, n_warm)
  logf("  summary time: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  logf("  parameters summarised: %d (of which y_mis imputations: %d)", nrow(s), sum(s$is_imputed))
  logf("  max R-hat, all sampled parameters incl. y_mis and lp__ : %.4f  (%s)", s$rhat[i_all], s$variable[i_all])
  logf("  max R-hat, model parameters (no y_mis, no lp__): %.4f  (%s)", s_model$rhat[i_model], s_model$variable[i_model])
  logf("  max R-hat, beta_state[1..%d]                    : %.4f  (%s)", nrow(s_beta), s_beta$rhat[i_beta], s_beta$variable[i_beta])
  logf("  min ESS bulk / tail (model parameters)         : %.0f / %.0f", min(s_model$ess_bulk, na.rm = TRUE), min(s_model$ess_tail, na.rm = TRUE))
  logf("  divergent transitions per chain                : %s (total %d)", paste(ds$num_divergent, collapse = "/"), sum(ds$num_divergent))
  logf("  max-treedepth hits per chain                   : %s", paste(ds$num_max_treedepth, collapse = "/"))
  logf("  all R-hat < 1.1 : %s | all R-hat < 1.05 : %s | all R-hat <= 1.01 : %s",
       all(s$rhat < 1.1, na.rm = TRUE), all(s$rhat < 1.05, na.rm = TRUE), all(s$rhat <= 1.01, na.rm = TRUE))

  summ[[cond]] <- tibble(
    mode = mode, condition = cond, n_params = nrow(s), n_y_mis = sum(s$is_imputed),
    max_rhat_all = s$rhat[i_all], max_rhat_all_param = s$variable[i_all],
    max_rhat_model = s_model$rhat[i_model], max_rhat_model_param = s_model$variable[i_model],
    max_rhat_beta = s_beta$rhat[i_beta], max_rhat_beta_param = s_beta$variable[i_beta],
    min_ess_bulk_model = min(s_model$ess_bulk, na.rm = TRUE), min_ess_tail_model = min(s_model$ess_tail, na.rm = TRUE),
    n_divergent = sum(ds$num_divergent), n_chain = n_chain, iter_sampling = n_iter, iter_warmup = n_warm
  )
}

params <- bind_rows(all_params)
write_csv(params, file.path(out_dir, sprintf("rhat_%s_params.csv", mode)))
S <- bind_rows(summ)
write_csv(S, file.path(out_dir, sprintf("rhat_%s_summary.csv", mode)))

i <- which.max(S$max_rhat_all)
logf("")
logf("--- across both conditions ---")
logf("  overall max R-hat (all variables) : %.4f  (%s, %s)", S$max_rhat_all[i], S$max_rhat_all_param[i], S$condition[i])
j <- which.max(S$max_rhat_model)
logf("  overall max R-hat (model params)  : %.4f  (%s, %s)", S$max_rhat_model[j], S$max_rhat_model_param[j], S$condition[j])
logf("  total divergent transitions       : %d", sum(S$n_divergent))
logf("")
logf("Suggested STAR Methods wording (numbers from this log; 3 decimals):")
logf("  \"... with all parameters showing R-hat < 1.1 (observed maximum %.3f, for %s in the %s condition) ...\"",
     S$max_rhat_all[i], S$max_rhat_all_param[i], str_remove(S$condition[i], "_pair"))
logf("Done.")
