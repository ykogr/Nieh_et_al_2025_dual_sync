#!/usr/bin/env Rscript
# =====================================================================
# 08_xscale_micro_macro.R  (co-author request, R2-3 follow-up)
# Dyad-level association between micro-scale synchrony and macro-scale
# rate coupling.  n = 20 dyads (10 visible, 10 invisible).
#
#   micro index (per dyad) : mean over the 3 sessions of the P3 z_S
#                            (delta = 0.1 s, grid-aligned null) taken from
#                            R1_2/out_periods_grid_<cond>/chew_sync_summary_periods_master.csv
#   macro index (per dyad) : mean over sessions and windows (bins 6-18, window
#                            end 60-180 s) of the Fisher-z windowed correlation
#                            of 10-s chew counts of the REAL pair, minus the
#                            same mean over all timeline-preserving SHUFFLED
#                            pairings that use this dyad's person A
#                            (dyad-level analogue of beta(t) = E_real - E_shuffled).
#                            Loader, binning, window rule and shuffle scheme are
#                            copied verbatim from 03_sensitivity_drift_beta.R (mode raw).
#
# Usage (zsh):
#   Rscript 08_xscale_micro_macro.R <data_dir> <micro_dir> <out_dir> [delta] [n_perm]
#     data_dir  : folder containing {visible_pair,invisible_pair}_Human/{A,B}/*.csv
#     micro_dir : folder containing out_periods_grid_visible/ and out_periods_grid_invisible/
#                 (normally ../R1_2)
#     out_dir   : output folder (created if absent), e.g. xscale_out
#     delta     : micro delta to use (default 0.1)
#     n_perm    : permutations for the permutation p-values (default 10000)
#
# Outputs (in out_dir):
#   xscale_dyad_table.csv    : one row per dyad (condition, micro_P3_zS, macro_real,
#                              macro_shuffled, macro_index, n_windows)
#   xscale_micro_macro.txt   : human-readable log (all numbers reported in the
#                              manuscript / letter come from this file)
#   xscale_micro_macro.csv   : one row per test (pooled, within-condition, per condition)
#   xscale_micro_macro.pdf/.png : scatter plot
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 08_xscale_micro_macro.R <data_dir> <micro_dir> <out_dir> [delta] [n_perm]")
}
base_dir  <- args[1]
micro_dir <- args[2]
out_dir   <- args[3]
delta_use <- if (length(args) >= 4) as.numeric(args[4]) else 0.1
n_perm    <- if (length(args) >= 5) as.integer(args[5]) else 10000L
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(out_dir, "xscale_micro_macro.txt")
cat("", file = log_file)
logf <- function(...) { s <- sprintf(...); cat(s, "\n", sep = ""); cat(s, "\n", sep = "", file = log_file, append = TRUE) }

logf("=====================================================================")
logf("08_xscale_micro_macro.R : dyad-level micro (P3 z_S) vs macro coupling : %s", format(Sys.time()))
logf("data_dir = %s | micro_dir = %s | delta = %.1f | n_perm = %d", base_dir, micro_dir, delta_use, n_perm)
logf("=====================================================================")

# ==== 1. Macro: load chew events (verbatim from 03_sensitivity_drift_beta.R) ====
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
  if (!is.na(m[1, 1])) return(as.numeric(m[1, 2]) * 60 + as.numeric(m[1, 3]))
  as.numeric(t)
}

chew_long <- map2_dfr(
  file_tbl$files, file_tbl$condition, function(f, condition) {
    fname <- basename(f)
    id_match <- str_match(fname, "^(\\d{8}_\\d+_\\d+)_([0-9]+)_([AB])\\.csv$")
    if (is.na(id_match[1, 1])) return(tibble())
    id_main <- id_match[1, 2]; session <- as.integer(id_match[1, 3]); person <- id_match[1, 4]
    df <- read_csv(f, col_types = cols(.default = col_character(), Trial = col_integer()),
                   na = c("", "NA"), show_col_types = FALSE)
    chew_cols_in_file <- intersect(names(df), chew_cols)
    df %>%
      mutate(trial = row_number()) %>%
      pivot_longer(cols = all_of(chew_cols_in_file), names_to = "chew_num", values_to = "time_str") %>%
      mutate(time_sec = map_dbl(time_str, time_to_sec),
             pair_id = id_main, session = session, person = person, condition = condition) %>%
      filter(!is.na(time_sec)) %>%
      dplyr::select(condition, pair_id, session, person, time_sec)
  }
)

# ==== 2. 10-s bin counts (legacy complete, as in the main pipeline) ====
bin_width <- 10; n_bin_all <- 180 / bin_width
bin_counts_long <- chew_long %>%
  mutate(bin = cut(time_sec, breaks = seq(0, bin_width * n_bin_all, by = bin_width),
                   labels = FALSE, include.lowest = TRUE, right = FALSE)) %>%
  filter(!is.na(bin)) %>%
  group_by(condition, pair_id, session, bin, person) %>%
  summarise(n_chew = n(), .groups = "drop") %>%
  complete(condition, pair_id, session, bin = 1:n_bin_all, person = c("A", "B"),
           fill = list(n_chew = 0)) %>%
  arrange(condition, pair_id, session, person, bin)

win_cor <- function(m) {
  if (any(is.na(m[, 1:2]))) return(NA_real_)
  if (all(m[, 3] + m[, 4] == 0)) return(NA_real_)
  r <- suppressWarnings(cor(m[, 1], m[, 2]))
  if (is.na(r)) return(NA_real_)
  r
}

wide <- bin_counts_long %>%
  mutate(value = as.numeric(n_chew)) %>%
  pivot_wider(names_from = person, values_from = c(value, n_chew), names_glue = "{.value}_{person}") %>%
  rename(val_A = value_A, val_B = value_B, raw_A = n_chew_A, raw_B = n_chew_B)

real <- wide %>%
  group_by(condition, pair_id, session) %>% arrange(bin, .by_group = TRUE) %>%
  mutate(rolling_cor = rollapply(cbind(val_A, val_B, raw_A, raw_B), width = 6, FUN = win_cor,
                                 by.column = FALSE, align = "right", fill = NA_real_),
         bin = as.integer(bin)) %>%
  ungroup() %>% mutate(rolling_zcor = atanh(rolling_cor), pair_type = "real_pair")

b_side <- wide %>% dplyr::select(condition, session, bin, pair_id, val_B, raw_B) %>%
  rename(other_pair_id = pair_id, val_B_shuf = val_B, raw_B_shuf = raw_B)
shuf <- wide %>% dplyr::select(condition, session, bin, pair_id, val_A, raw_A) %>%
  expand_grid(other_pair_id = unique(wide$pair_id)) %>%
  filter(pair_id != other_pair_id) %>%
  left_join(b_side, by = c("condition", "session", "bin", "other_pair_id")) %>%
  filter(!is.na(raw_B_shuf)) %>%
  group_by(condition, session, pair_id, other_pair_id) %>% arrange(bin, .by_group = TRUE) %>%
  mutate(rolling_cor = rollapply(cbind(val_A, val_B_shuf, raw_A, raw_B_shuf), width = 6, FUN = win_cor,
                                 by.column = FALSE, align = "right", fill = NA_real_),
         bin = as.integer(bin)) %>%
  ungroup() %>% mutate(rolling_zcor = atanh(rolling_cor), pair_type = "shuffled_pair")

keep_ok <- function(d) d %>% filter(bin >= 6, is.finite(rolling_zcor))

# The legacy complete() also creates rows for pair_ids that do not belong to a
# condition (all-zero counts); those windows are NA by the window rule and drop out here.
macro_real <- keep_ok(real) %>% group_by(condition, pair_id) %>%
  summarise(macro_real = mean(rolling_zcor), n_windows_real = n(), .groups = "drop")
macro_shuf <- keep_ok(shuf) %>% group_by(condition, pair_id) %>%
  summarise(macro_shuffled = mean(rolling_zcor), n_windows_shuf = n(), .groups = "drop")
macro <- inner_join(macro_real, macro_shuf, by = c("condition", "pair_id")) %>%
  filter(n_windows_real > 0) %>%
  mutate(macro_index = macro_real - macro_shuffled,
         condition = str_remove(condition, "_pair"))

# ==== 3. Micro: dyad mean of P3 z_S ====
read_master <- function(cond) {
  f <- file.path(micro_dir, sprintf("out_periods_grid_%s", cond), "chew_sync_summary_periods_master.csv")
  if (!file.exists(f)) stop("micro master table not found: ", f)
  read_csv(f, show_col_types = FALSE) %>% mutate(src = f)
}
micro_raw <- bind_rows(read_master("visible"), read_master("invisible"))
logf("micro inputs: %s", paste(unique(micro_raw$src), collapse = " ; "))
micro <- micro_raw %>%
  filter(abs(delta - delta_use) < 1e-9, period == "P3", sufficient == 1) %>%
  mutate(dyad_id = str_replace(pair_id, "_\\d+$", "")) %>%
  group_by(condition, dyad_id) %>%
  summarise(micro_P3_zS = mean(z_S), n_sessions = n(), .groups = "drop")

tab <- inner_join(micro, macro, by = c("condition", "dyad_id" = "pair_id")) %>%
  arrange(condition, dyad_id)
write_csv(tab, file.path(out_dir, "xscale_dyad_table.csv"))

logf("")
logf("dyads matched: %d (visible %d, invisible %d)", nrow(tab), sum(tab$condition == "visible"), sum(tab$condition == "invisible"))
logf("micro-only dyads: %d | macro-only dyads: %d", nrow(anti_join(micro, macro, by = c("condition", "dyad_id" = "pair_id"))),
     nrow(anti_join(macro, micro, by = c("condition", "pair_id" = "dyad_id"))))
logf("")
logf("--- per-dyad table (also xscale_dyad_table.csv) ---")
for (i in seq_len(nrow(tab))) {
  logf("  %-10s %-16s micro_P3_zS = %+.3f (n_sess %d) | macro_real = %+.3f  macro_shuffled = %+.3f  macro_index = %+.3f (windows %d)",
       tab$condition[i], tab$dyad_id[i], tab$micro_P3_zS[i], tab$n_sessions[i], tab$macro_real[i],
       tab$macro_shuffled[i], tab$macro_index[i], tab$n_windows_real[i])
}
logf("")
logf("--- condition means (sanity check against Figure 3 / Figure 5a) ---")
cm <- tab %>% group_by(condition) %>% summarise(across(c(micro_P3_zS, macro_real, macro_shuffled, macro_index), mean), .groups = "drop")
for (i in seq_len(nrow(cm))) logf("  %-10s micro_P3_zS %+.3f | macro_real %+.3f | macro_shuffled %+.3f | macro_index %+.3f",
                                  cm$condition[i], cm$micro_P3_zS[i], cm$macro_real[i], cm$macro_shuffled[i], cm$macro_index[i])

# ==== 4. Correlations with permutation p-values ====
set.seed(20260914)
perm_p <- function(x, y, method, strata = NULL) {
  obs <- suppressWarnings(cor(x, y, method = method))
  cnt <- 0L
  for (b in seq_len(n_perm)) {
    if (is.null(strata)) {
      yp <- sample(y)
    } else {
      yp <- y
      for (s in unique(strata)) { idx <- which(strata == s); yp[idx] <- y[sample(idx)] }
    }
    if (abs(suppressWarnings(cor(x, yp, method = method))) >= abs(obs) - 1e-12) cnt <- cnt + 1L
  }
  list(r = obs, p = (cnt + 1) / (n_perm + 1))
}
pearson_ci <- function(x, y) {
  n <- length(x); r <- cor(x, y)
  if (n < 4) return(c(NA, NA))
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  tanh(z + c(-1, 1) * qnorm(0.975) * se)
}

run_block <- function(label, d, strata = NULL, macro_col = "macro_index") {
  x <- d$micro_P3_zS; y <- d[[macro_col]]
  sp <- perm_p(x, y, "spearman", strata)
  pe <- perm_p(x, y, "pearson", strata)
  ci <- pearson_ci(x, y)
  logf("  [%s | macro = %s] n = %d : Spearman rho = %+.3f (perm p = %.4f) | Pearson r = %+.3f [%+.3f, %+.3f] (perm p = %.4f)",
       label, macro_col, nrow(d), sp$r, sp$p, pe$r, ci[1], ci[2], pe$p)
  tibble(analysis = label, macro_measure = macro_col, n = nrow(d),
         spearman_rho = sp$r, spearman_perm_p = sp$p,
         pearson_r = pe$r, pearson_ci_low = ci[1], pearson_ci_high = ci[2], pearson_perm_p = pe$p)
}

logf("")
logf("--- association between dyad-level micro P3 z_S and dyad-level macro coupling ---")
logf("  permutation p: two-sided, %d permutations of the macro values (within condition for the 'within-condition' rows); seed 20260914", n_perm)
centered <- tab %>% group_by(condition) %>%
  mutate(micro_P3_zS = micro_P3_zS - mean(micro_P3_zS),
         macro_index = macro_index - mean(macro_index),
         macro_real  = macro_real  - mean(macro_real)) %>% ungroup()
res <- bind_rows(
  run_block("pooled (20 dyads)", tab),
  run_block("within-condition (condition-centred, strata-permuted)", centered, strata = centered$condition),
  run_block("visible only", filter(tab, condition == "visible")),
  run_block("invisible only", filter(tab, condition == "invisible")),
  run_block("pooled (20 dyads)", tab, macro_col = "macro_real"),
  run_block("within-condition (condition-centred, strata-permuted)", centered, strata = centered$condition, macro_col = "macro_real")
)
write_csv(res, file.path(out_dir, "xscale_micro_macro.csv"))

prim <- res[2, ]  # within-condition, macro_index = primary (avoids the condition-mean confound)
logf("")
logf("PRIMARY (within-condition, macro_index): Spearman rho = %+.3f, perm p = %.4f, n = %d", prim$spearman_rho, prim$spearman_perm_p, prim$n)
logf("reading: the pooled correlation mixes the between-condition difference (visible: higher micro, lower macro) with the")
logf("         dyad-level association; the within-condition row removes the condition means and is the one to report.")

# ==== 5. Figure ====
p <- ggplot(tab, aes(x = micro_P3_zS, y = macro_index, colour = condition, shape = condition)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c(visible = "#20808D", invisible = "#A84B2F")) +
  labs(x = expression(paste("Dyad-level micro-scale synchrony, P3 ", italic(z)[S], " (", delta, " = 0.1 s)")),
       y = "Dyad-level macro-scale coupling\n(real - shuffled Fisher-z windowed correlation)",
       colour = NULL, shape = NULL,
       caption = sprintf("within-condition Spearman rho = %+.2f, permutation p = %.3f (n = %d dyads)",
                         prim$spearman_rho, prim$spearman_perm_p, prim$n)) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out_dir, "xscale_micro_macro.pdf"), p, width = 4.6, height = 4.2)
ggsave(file.path(out_dir, "xscale_micro_macro.png"), p, width = 4.6, height = 4.2, dpi = 200)
logf("")
logf("figure: %s (.pdf/.png)", file.path(out_dir, "xscale_micro_macro"))
logf("Done.")
