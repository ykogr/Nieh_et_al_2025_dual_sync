#!/usr/bin/env Rscript
## =====================================================================
## 07_plot_figure4e_gridnull.R
##   Fig. 4e: posterior mean and 95% credible interval of beta_state(t)
##   (visible - invisible contrast of z_S(t)) from the Stan fit produced by
##   06_bayes_gridnull_fit.R, with shaded bands where the interval excludes 0.
##
##   The plotting code is the chunk "Plot / delta=0.1" of
##   04_Micro_Analysis_Bayesian.Rmd in the public repository
##   (https://github.com/ykogr/Nieh_et_al_2025_dual_sync), kept verbatim
##   except for two additions:
##     (1) the x axis is video time (s) instead of the bin index: bin i is
##         mapped to the i-th distinct t_center of
##         Micro_Analysis_Bayesian/{vis,invis}_<delta>.csv, exactly as
##         bin_id = as.integer(factor(t_center)) does inside the fit script;
##     (2) the figure is saved to files instead of being printed in a notebook.
##
## Usage (in gridnull_main/):
##   Rscript 07_plot_figure4e_gridnull.R [delta_tag] [delta] [outdir]
##     delta_tag : suffix of the summary file (default "01"  -> summary_z_S_t_delta01_human.csv)
##     delta     : suffix of the input CSVs   (default "0.1" -> vis_0.1.csv / invis_0.1.csv)
##     outdir    : output folder              (default "fig4e_gridnull")
##   -> fig4e_gridnull/fig4e_beta_delta01.pdf / .png
##      fig4e_gridnull/fig4e_beta_delta01_bins.csv  (bin, t_center, mean, q2.5, q97.5, sigpos, signeg)
## =====================================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

args      <- commandArgs(trailingOnly = TRUE)
delta_tag <- if (length(args) >= 1) args[1] else "01"
delta     <- if (length(args) >= 2) args[2] else "0.1"
outdir    <- if (length(args) >= 3) args[3] else "fig4e_gridnull"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

summary_file <- file.path("stan_fits", sprintf("summary_z_S_t_delta%s_human.csv", delta_tag))
vis_file     <- file.path("Micro_Analysis_Bayesian", sprintf("vis_%s.csv", delta))
invis_file   <- file.path("Micro_Analysis_Bayesian", sprintf("invis_%s.csv", delta))
for (f in c(summary_file, vis_file, invis_file))
  if (!file.exists(f)) stop("missing input: ", f)

## ---- bin -> t_center mapping (same construction as make_stan_data) -------
t_centers <- sort(unique(c(
  read_csv(vis_file,   show_col_types = FALSE)$t_center,
  read_csv(invis_file, show_col_types = FALSE)$t_center)))

## ---- 1. Loading summary files  (verbatim from the Rmd) --------------------
df_sum <- read_csv(summary_file, show_col_types = FALSE)

## ---- 2. Extract beta_state, bin number, and determine the confidence interval
df_beta <- df_sum %>%
  filter(str_detect(variable, "^beta_state\\[")) %>%
  mutate(
    bin = as.integer(str_extract(variable, "\\d+")),
    mean = mean,
    q025 = q2.5,
    q975 = q97.5,
    sigpos = if_else(q025 > 0, 1, 0),    # Positive direction significant
    signeg = if_else(q975 < 0, 1, 0)     # Negative direction significant
  ) %>%
  arrange(bin)

if (nrow(df_beta) != length(t_centers))
  stop(sprintf("%d beta_state bins but %d distinct t_center values in %s/%s",
               nrow(df_beta), length(t_centers), vis_file, invis_file))
df_beta$t_center <- t_centers[df_beta$bin]          # addition (1)

## ---- 3. Function to obtain consecutive intervals using run-length encoding
get_sig_ranges <- function(sigvec) {
  rle_sig <- rle(sigvec)
  lens <- rle_sig$lengths
  vals <- rle_sig$values
  ends <- cumsum(lens)
  starts <- ends - lens + 1
  tibble(start = starts[vals == 1], end = ends[vals == 1])
}
sig_ranges_pos <- get_sig_ranges(df_beta$sigpos)
sig_ranges_neg <- get_sig_ranges(df_beta$signeg)
## bin index -> seconds for the shaded bands (bin edges = t_center -/+ 0.5 s)
to_sec <- function(r) r %>% mutate(start = t_centers[start] - 0.5, end = t_centers[end] + 0.5)
sig_ranges_pos <- to_sec(sig_ranges_pos)
sig_ranges_neg <- to_sec(sig_ranges_neg)

## ---- 4. Plot  (verbatim except x = t_center) ------------------------------
p <- ggplot(df_beta, aes(x = t_center, y = mean)) +
  # Top (red) band
  geom_rect(
    data = sig_ranges_pos,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#D55E00", alpha = 0.10
  ) +
  # Bottom (blue) band
  geom_rect(
    data = sig_ranges_neg,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#0072B2", alpha = 0.10
  ) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2, fill = "grey40") +
  labs(x = "Time (s)", y = expression(beta(t) ~ "mean ± 95%CI")) +
  theme_classic(base_size = 16) +
  geom_hline(yintercept = 0)

base <- file.path(outdir, sprintf("fig4e_beta_delta%s", delta_tag))
ggsave(paste0(base, ".pdf"), p, width = 9, height = 4.5)
ggsave(paste0(base, ".png"), p, width = 9, height = 4.5, dpi = 300)
write_csv(df_beta %>% select(bin, t_center, mean, q025, q975, sigpos, signeg),
          paste0(base, "_bins.csv"))
cat(sprintf("Saved: %s.pdf / .png / _bins.csv  (positive %d s, negative %d s)\n",
            base, sum(df_beta$sigpos), sum(df_beta$signeg)))
