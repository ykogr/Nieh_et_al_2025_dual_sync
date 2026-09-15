#!/usr/bin/env Rscript
# =====================================================================
# 03_sensitivity_drift_beta.R
# Sensitivity analysis for R2-6: shared slow-drift control for beta(t)
#
# Modes:
#   raw     : replication of the original pipeline (sanity check)
#   detrend : linear trend removed per person x session before windowing
#   diff    : first-differenced bin counts (stronger drift control)
#
# Usage (zsh):
#   Rscript 03_sensitivity_drift_beta.R <data_dir> <stan_file> <out_dir> [modes] [legacy_complete]
#     data_dir        : folder containing {visible_pair,invisible_pair}_Human/{A,B}/*.csv
#     stan_file       : path to 02_ssm_chewing_sync_missing.stan
#     out_dir         : output folder (created if absent)
#     modes           : comma-separated subset of raw,detrend,diff (default: raw,detrend,diff)
#     legacy_complete : TRUE (default) reproduces the original complete() call, which
#                       crosses condition x pair_id fully; FALSE restricts zero-filling
#                       to observed condition/pair/session combinations
#
# Outputs (in out_dir):
#   beta_summary_<mode>_<condition>.csv : posterior summary of beta_state per bin
#   condition_means_<mode>.csv          : across-bin mean of posterior means + CIs
#   beta_summary_all.csv                : combined summaries for all modes run
#   fit_<mode>_<condition>.rds          : saved cmdstanr fit objects
#   diagnostics.txt                     : data dimensions, Rhat, divergences
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 03_sensitivity_drift_beta.R <data_dir> <stan_file> <out_dir> [modes] [legacy_complete]")
}
base_dir  <- args[1]
stan_file <- args[2]
out_dir   <- args[3]
modes     <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else c("raw", "detrend", "diff")
legacy_complete <- if (length(args) >= 5) as.logical(args[5]) else TRUE

stopifnot(all(modes %in% c("raw", "detrend", "diff")))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
diag_file <- file.path(out_dir, "diagnostics.txt")
log_diag <- function(...) cat(sprintf(...), "\n", file = diag_file, append = TRUE)
cat(sprintf("Run started: %s\n", Sys.time()), file = diag_file)
log_diag("modes = %s | legacy_complete = %s", paste(modes, collapse = ","), legacy_complete)

# ==== 1. Load annotated chewing events (identical to original pipeline) ====
conds   <- c("visible_pair", "invisible_pair")
persons <- c("A", "B")

file_tbl <- expand_grid(condition = conds, person = persons) %>%
  mutate(
    folder = file.path(base_dir, paste0(condition, "_Human"), person),
    files  = map(folder, ~list.files(.x, pattern = "\\.csv$", full.names = TRUE))
  ) %>%
  unnest(files)

if (nrow(file_tbl) == 0) stop("No CSV files found under ", base_dir)

header    <- names(read_csv(file_tbl$files[1], n_max = 0, show_col_types = FALSE))
chew_cols <- header[str_detect(header, "^Chew_")]

time_to_sec <- function(t) {
  if (is.na(t) || t == "") return(NA_real_)
  m <- str_match(t, "^(\\d+):(\\d+\\.\\d+)$")
  if (!is.na(m[1, 1])) {
    return(as.numeric(m[1, 2]) * 60 + as.numeric(m[1, 3]))
  }
  as.numeric(t)
}

chew_long <- map2_dfr(
  file_tbl$files, file_tbl$condition, function(f, condition) {
    fname <- basename(f)
    id_match <- str_match(fname, "^(\\d{8}_\\d+_\\d+)_([0-9]+)_([AB])\\.csv$")
    if (is.na(id_match[1, 1])) return(tibble())
    id_main <- id_match[1, 2]
    session <- as.integer(id_match[1, 3])
    person  <- id_match[1, 4]

    df <- read_csv(
      f,
      col_types = cols(.default = col_character(), Trial = col_integer()),
      na = c("", "NA"),
      show_col_types = FALSE
    )
    chew_cols_in_file <- intersect(names(df), chew_cols)
    df %>%
      mutate(trial = row_number()) %>%
      pivot_longer(
        cols = all_of(chew_cols_in_file),
        names_to = "chew_num", values_to = "time_str"
      ) %>%
      mutate(
        time_sec = map_dbl(time_str, time_to_sec),
        pair_id = id_main, session = session,
        person = person, condition = condition
      ) %>%
      filter(!is.na(time_sec)) %>%
      dplyr::select(condition, pair_id, session, person, time_sec)
  }
)

# ==== 2. 10-s bin counts, zero-filled ====
bin_width <- 10
n_bin_all <- 180 / bin_width  # 18

bin_counts_base <- chew_long %>%
  mutate(
    bin = cut(time_sec, breaks = seq(0, bin_width * n_bin_all, by = bin_width),
              labels = FALSE, include.lowest = TRUE, right = FALSE)
  ) %>%
  filter(!is.na(bin)) %>%
  group_by(condition, pair_id, session, bin, person) %>%
  summarise(n_chew = n(), .groups = "drop")

if (legacy_complete) {
  bin_counts_long <- bin_counts_base %>%
    complete(
      condition, pair_id, session, bin = 1:n_bin_all, person = c("A", "B"),
      fill = list(n_chew = 0)
    )
} else {
  bin_counts_long <- bin_counts_base %>%
    complete(
      nesting(condition, pair_id, session), bin = 1:n_bin_all, person = c("A", "B"),
      fill = list(n_chew = 0)
    )
}
bin_counts_long <- bin_counts_long %>% arrange(condition, pair_id, session, person, bin)

# ==== 3. Per-mode value series ====
add_value <- function(df, mode) {
  df <- df %>% group_by(condition, pair_id, session, person) %>% arrange(bin, .by_group = TRUE)
  if (mode == "raw") {
    df <- df %>% mutate(value = as.numeric(n_chew))
  } else if (mode == "detrend") {
    df <- df %>% mutate(value = as.numeric(resid(lm(n_chew ~ bin))))
  } else if (mode == "diff") {
    df <- df %>% mutate(value = n_chew - lag(n_chew))
  }
  df %>% ungroup()
}

# Windowed correlation with missing rules matching the original pipeline:
#  - NA if any value in the window is NA
#  - NA if BOTH persons' RAW counts are zero across the whole window
#  - NA if correlation is undefined (zero variance)
win_cor <- function(m) {
  if (any(is.na(m[, 1:2]))) return(NA_real_)
  if (all(m[, 3] + m[, 4] == 0)) return(NA_real_)
  r <- suppressWarnings(cor(m[, 1], m[, 2]))
  if (is.na(r)) return(NA_real_)
  r
}

roll_series <- function(df_wide, width) {
  df_wide %>%
    group_by(condition, pair_id, session) %>%
    arrange(bin, .by_group = TRUE) %>%
    mutate(
      rolling_cor = rollapply(
        data = cbind(val_A, val_B, raw_A, raw_B),
        width = width, FUN = win_cor,
        by.column = FALSE, align = "right", fill = NA_real_
      ),
      bin = as.integer(bin)
    ) %>%
    ungroup() %>%
    mutate(rolling_zcor = atanh(rolling_cor))
}

build_cor_long <- function(mode) {
  # For diff mode a 60-s window of 6 counts yields 5 first differences,
  # so the rolling width is 5; valid time points remain bins 6-18.
  width_use <- if (mode == "diff") 5 else 6

  vals <- add_value(bin_counts_long, mode) %>%
    dplyr::select(condition, pair_id, session, bin, person, n_chew, value)

  wide <- vals %>%
    pivot_wider(names_from = person, values_from = c(value, n_chew),
                names_glue = "{.value}_{person}") %>%
    rename(val_A = value_A, val_B = value_B, raw_A = n_chew_A, raw_B = n_chew_B)

  real <- roll_series(wide, width_use) %>%
    mutate(pair_type = "real_pair") %>%
    dplyr::select(condition, session, pair_id, bin, rolling_cor, rolling_zcor, pair_type)

  # Shuffled pairs: A of pair j with B of pair k (k != j),
  # matched on condition, session, and bin (timeline-preserving)
  b_side <- wide %>%
    dplyr::select(condition, session, bin, pair_id, val_B, raw_B) %>%
    rename(other_pair_id = pair_id, val_B_shuf = val_B, raw_B_shuf = raw_B)

  shuf_wide <- wide %>%
    dplyr::select(condition, session, bin, pair_id, val_A, raw_A) %>%
    expand_grid(other_pair_id = unique(wide$pair_id)) %>%
    filter(pair_id != other_pair_id) %>%
    left_join(b_side, by = c("condition", "session", "bin", "other_pair_id")) %>%
    filter(!is.na(raw_B_shuf))

  shuf <- shuf_wide %>%
    group_by(condition, session, pair_id, other_pair_id) %>%
    arrange(bin, .by_group = TRUE) %>%
    mutate(
      rolling_cor = rollapply(
        data = cbind(val_A, val_B_shuf, raw_A, raw_B_shuf),
        width = width_use, FUN = win_cor,
        by.column = FALSE, align = "right", fill = NA_real_
      ),
      bin = as.integer(bin)
    ) %>%
    ungroup() %>%
    mutate(rolling_zcor = atanh(rolling_cor), pair_type = "shuffled_pair") %>%
    dplyr::select(condition, session, pair_id, bin, rolling_cor, rolling_zcor, pair_type)

  bind_rows(real, shuf)
}

# ==== 4. Stan data construction (identical to original) ====
make_stan_data_miss <- function(dat) {
  dat <- dat %>% filter(bin >= 6)
  dat <- dat %>% mutate(rolling_zcor = ifelse(is.infinite(rolling_zcor), NA, rolling_zcor))
  dat <- dat %>% arrange(pair_id, session, bin, pair_type)

  y_vec  <- dat$rolling_zcor
  obs_idx <- which(!is.na(y_vec))
  mis_idx <- which(is.na(y_vec))

  list(
    N = nrow(dat),
    N_obs = length(obs_idx),
    N_mis = length(mis_idx),
    N_pair = n_distinct(dat$pair_id),
    N_session = n_distinct(dat$session),
    N_bin = 13,
    obs_idx = obs_idx,
    mis_idx = mis_idx,
    pair_id = as.integer(as.factor(dat$pair_id)),
    session_id = as.integer(as.factor(dat$session)),
    bin = dat$bin - 5,
    is_real = as.integer(dat$pair_type == "real_pair"),
    y_obs = y_vec[obs_idx]
  )
}

make_init_fun <- function(stan_data) {
  y_init <- if (length(stan_data$y_obs) > 0) mean(stan_data$y_obs, na.rm = TRUE) else 0
  list(
    y_mis = rep(y_init, stan_data$N_mis),
    mu0 = y_init,
    mu_pair = rep(y_init, stan_data$N_pair),
    mu_sess_init = matrix(y_init, stan_data$N_pair, stan_data$N_session),
    sigma_pair = 1, sigma_sess_init = 1,
    mu_state = rep(y_init, stan_data$N_bin),
    sigma_mu = 1,
    beta0 = 0,
    beta_tan = rep(0, stan_data$N_bin - 1),
    sigma_beta = 1,
    sigma_obs = 1
  )
}

# ==== 5. Fit per mode x condition ====
mod <- cmdstan_model(stan_file)
options(mc.cores = parallel::detectCores())
set.seed(1234)

all_beta <- list()

for (mode in modes) {
  message(sprintf("==== Mode: %s ====", mode))
  cor_long <- build_cor_long(mode)
  cond_means <- list()

  for (cond in conds) {
    dat_sub <- cor_long %>% filter(condition == cond)
    sd_ <- make_stan_data_miss(dat_sub)

    log_diag("[%s | %s] N=%d N_obs=%d N_mis=%d N_pair=%d N_session=%d",
             mode, cond, sd_$N, sd_$N_obs, sd_$N_mis, sd_$N_pair, sd_$N_session)
    message(sprintf("  %s: N=%d, N_obs=%d, N_mis=%d, N_pair=%d",
                    cond, sd_$N, sd_$N_obs, sd_$N_mis, sd_$N_pair))

    fit <- mod$sample(
      data = sd_,
      chains = 4, parallel_chains = 4,
      iter_sampling = 5000, iter_warmup = 1000,
      seed = 123, refresh = 500,
      adapt_delta = 0.99, max_treedepth = 15,
      init = function() make_init_fun(sd_)
    )
    fit$save_object(file = file.path(out_dir, paste0("fit_", mode, "_", cond, ".rds")))

    beta_sum <- fit$summary(
      variables = "beta_state",
      mean, median, sd,
      ~quantile2(.x, probs = c(0.025, 0.975)),
      rhat, ess_bulk, ess_tail
    ) %>%
      mutate(
        bin_model = as.integer(str_extract(variable, "\\d+")),
        bin_data = bin_model + 5,
        window_end_s = bin_data * 10,
        mode = mode, condition = cond
      )
    write_csv(beta_sum, file.path(out_dir, paste0("beta_summary_", mode, "_", cond, ".csv")))
    all_beta[[paste(mode, cond)]] <- beta_sum

    ds <- fit$diagnostic_summary()
    log_diag("[%s | %s] divergences=%s max_treedepth_hits=%s max_rhat(beta_state)=%.4f",
             mode, cond, paste(ds$num_divergent, collapse = "/"),
             paste(ds$num_max_treedepth, collapse = "/"),
             max(beta_sum$rhat, na.rm = TRUE))

    m  <- beta_sum$mean
    nb <- length(m)
    se <- sd(m) / sqrt(nb)
    cond_means[[cond]] <- tibble(
      mode = mode, condition = cond, n_bins = nb,
      mean_beta = mean(m), sd_across_bins = sd(m), se_across_bins = se,
      ci95_normal_low = mean(m) - 1.96 * se,
      ci95_normal_high = mean(m) + 1.96 * se,
      ci95_t_low = mean(m) + qt(0.025, nb - 1) * se,
      ci95_t_high = mean(m) + qt(0.975, nb - 1) * se,
      n_bins_cri_excl_zero = sum(beta_sum$q2.5 > 0 | beta_sum$q97.5 < 0),
      max_rhat = max(beta_sum$rhat, na.rm = TRUE)
    )
  }
  write_csv(bind_rows(cond_means), file.path(out_dir, paste0("condition_means_", mode, ".csv")))
}

write_csv(bind_rows(all_beta), file.path(out_dir, "beta_summary_all.csv"))
log_diag("Run finished: %s", format(Sys.time()))
message("Done. Outputs in: ", out_dir)
