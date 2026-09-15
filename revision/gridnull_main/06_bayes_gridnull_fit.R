# 06_bayes_gridnull_fit.R
# Purpose : STEP 6 (Figure 4e, Bayesian beta(t)) -- fit and/or recover.
#           The fitting code below replicates the fitting chunks of
#           04_Micro_Analysis_Bayesian.Rmd VERBATIM (same data prep, same
#           model, same sampler settings, same seed = 123 per fit), with
#           three operational fixes only:
#             (1) dir.create("stan_fits") so save_object() cannot fail;
#             (2) output_dir passed to $sample() so CmdStan chain CSVs
#                 persist on disk (a later crash can never lose sampling);
#             (3) the set of deltas to fit is taken from the command line,
#                 so already-completed fits are not redone.
#           Because seed = 123 is set independently inside every $sample()
#           call, fitting deltas one at a time yields results identical to
#           the original all-in-one loop.
#
# Usage (run in ROOT/gridnull_main, where 04_ssm_chewing_missing_Bayesian.stan
#        and Micro_Analysis_Bayesian/ live):
#
#   Fit mode    : Rscript 06_bayes_gridnull_fit.R fit 0.1,0.2,0.3,0.5,1.0
#   Recover mode: Rscript 06_bayes_gridnull_fit.R recover <dir_with_chain_csvs> 0.1
#
# Recover mode rebuilds the fit object from surviving CmdStan chain CSVs of a
# finished run (no re-sampling) and writes the same two outputs as fit mode:
#   stan_fits/fit_<nm>_human.rds  and  stan_fits/summary_<nm>_human.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 06_bayes_gridnull_fit.R fit <deltas-comma-separated>\n",
       "   or: Rscript 06_bayes_gridnull_fit.R recover <csv_dir> <delta>")
}
mode <- args[1]

output_dir <- "stan_fits"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

nm_of_delta <- function(delta) {
  paste0("z_S_t_delta", gsub("\\.", "", as.character(delta)))
}

write_outputs <- function(fit, nm) {
  fit$save_object(file = file.path(output_dir, paste0("fit_", nm, "_human.rds")))
  summary_df <- fit$summary(
    variables = NULL,
    posterior::default_summary_measures(),
    posterior::default_convergence_measures(),
    extra_quantiles = ~posterior::quantile2(., probs = c(.025, .975))
  )
  summary_df$condition <- nm
  write.csv(summary_df,
            file = file.path(output_dir, paste0("summary_", nm, "_human.csv")),
            row.names = FALSE)
  message(sprintf("Wrote %s/fit_%s_human.rds and %s/summary_%s_human.csv",
                  output_dir, nm, output_dir, nm))
}

# ---------------------------------------------------------------------------
# Recover mode: rebuild fit from surviving CmdStan chain CSVs (no sampling)
# ---------------------------------------------------------------------------
if (mode == "recover") {
  csv_dir <- args[2]
  delta   <- as.numeric(args[3])
  nm      <- nm_of_delta(delta)
  files <- list.files(csv_dir, pattern = "^model-.*\\.csv$", full.names = TRUE)
  message(sprintf("Found %d chain CSV file(s) in %s:", length(files), csv_dir))
  for (f in files) message("  ", f)
  if (length(files) != 4) {
    stop("Expected exactly 4 chain CSVs (one per chain). ",
         "If other model CSVs are mixed in, move the 4 target files ",
         "into an empty folder and pass that folder.")
  }
  fit <- cmdstanr::as_cmdstan_fit(files)
  write_outputs(fit, nm)
  quit(save = "no", status = 0)
}

if (mode != "fit") stop("First argument must be 'fit' or 'recover'.")
delta_vals <- as.numeric(strsplit(args[2], ",")[[1]])
if (any(is.na(delta_vals))) stop("Could not parse deltas: ", args[2])

# ---------------------------------------------------------------------------
# Data loading -- verbatim from 04_Micro_Analysis_Bayesian.Rmd
# ---------------------------------------------------------------------------
files <- list.files("Micro_Analysis_Bayesian", pattern = "vis_|invis_", full.names = TRUE)

d_timeline_20251002 <- map_dfr(files, function(f) {
  fname <- basename(f)
  condition <- ifelse(str_detect(fname, "^vis"), "visible_pair", "invisible_pair")
  delta <- as.numeric(str_extract(fname, "\\d+\\.\\d+"))
  dat <- read_csv(f, show_col_types = FALSE) %>%
    pivot_longer(
      cols = contains("z_S_t"),
      names_to = "key",
      values_to = "z_S_t"
    ) %>%
    mutate(
      pair_session = str_remove(key, "z_S_t__"),
      pair_id = str_extract(pair_session, "^\\d{8}_\\d{2}_\\d{2}"),
      session = str_extract(pair_session, "\\d{2}$"),
      condition = condition,
      delta = delta
    ) %>%
    dplyr::select(condition, pair_id, session, t_center, z_S_t, delta)
  dat
})

d_timeline_20251002 <- d_timeline_20251002 %>% mutate(session = as.integer(session))

make_stan_data <- function(df, resp_var) {
  df <- df %>%
    mutate(
      pair_id_int = as.integer(factor(pair_id)),
      session_id_int = as.integer(factor(session)),
      bin_id = as.integer(factor(t_center)),
      is_visible = if_else(condition == "visible_pair", 1L, 0L)
    )
  y_vec <- df[[resp_var]]
  obs_idx <- which(!is.na(y_vec))
  mis_idx <- which(is.na(y_vec))
  stan_data <- list(
    N         = nrow(df),
    N_obs     = length(obs_idx),
    N_mis     = length(mis_idx),
    N_pair    = max(df$pair_id_int),
    N_session = max(df$session_id_int),
    N_bin     = max(df$bin_id),
    obs_idx   = obs_idx,
    mis_idx   = mis_idx,
    pair_id   = df$pair_id_int,
    session_id= df$session_id_int,
    bin       = df$bin_id,
    is_visible= df$is_visible,
    y_obs     = y_vec[obs_idx]
  )
  return(stan_data)
}

make_init_fun <- function(stan_data) {
  y_init <- if (length(stan_data$y_obs) > 0) mean(stan_data$y_obs, na.rm = TRUE) else 0
  list(
    y_mis = rep(y_init, stan_data$N_mis),
    mu0 = y_init,
    mu_pair = rep(y_init, stan_data$N_pair),
    mu_sess_init = matrix(y_init, stan_data$N_pair, stan_data$N_session),
    sigma_pair = 1,
    sigma_sess_init = 1,
    mu_state = rep(y_init, stan_data$N_bin),
    sigma_mu = 1,
    beta0 = 0,
    beta_tan = rep(0, stan_data$N_bin - 1),
    sigma_beta = 1,
    sigma_obs = 1
  )
}

resp_vars <- c("z_S_t")

stan_data_list <- list()
for (delta in delta_vals) {
  for (resp in resp_vars) {
    df_sub <- d_timeline_20251002 %>% filter(delta == !!delta)
    if (nrow(df_sub) == 0) stop("No rows for delta = ", delta,
                                " -- check Micro_Analysis_Bayesian/ inputs.")
    nm <- paste0(resp, "_delta", gsub("\\.", "", as.character(delta)))
    stan_data_list[[nm]] <- make_stan_data(df_sub, resp)
  }
}

mod <- cmdstan_model("04_ssm_chewing_missing_Bayesian.stan")

options(mc.cores = parallel::detectCores())
set.seed(1234)

fit_list <- list()

for (nm in names(stan_data_list)) {
  message(sprintf("Fitting: %s", nm))
  csv_dir_nm <- file.path(output_dir, paste0("csv_", nm))
  dir.create(csv_dir_nm, showWarnings = FALSE, recursive = TRUE)
  fit <- mod$sample(
    data = stan_data_list[[nm]],
    chains = 4, parallel_chains = 4,
    iter_sampling = 5000, iter_warmup = 1000,
    seed = 123,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 15,
    init = function() make_init_fun(stan_data_list[[nm]]),
    output_dir = csv_dir_nm
  )
  fit_list[[nm]] <- fit
  write_outputs(fit, nm)
}

message("Done.")
