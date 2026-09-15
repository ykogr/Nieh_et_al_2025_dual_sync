#!/usr/bin/env Rscript
# =====================================================================
# R3-1_power_sensitivity.R
# Simulation-based sensitivity power analysis for the period-resolved
# micro-scale LMM (Nieh et al. 2026, iScience revision; reviewer R3-1).
#
# Design (fixed by the original study, Ogura et al. 2020):
#   2 conditions (visible / invisible) x 10 dyads per condition
#   x 3 sessions per dyad x 3 periods (P1, P2, P3)
#   DV = session-level standardized synchrony z_S
#   Model: z_S ~ condition * period + (1 | dyad_id)
#
# Question answered: what is the minimal detectable effect size (MDES)
# for the primary estimand -- the visible - invisible contrast in P3 --
# at 80% power, given the observed variance components?
#
# Parameter defaults are taken from the fitted R1-2 period-resolved LMM
# with the grid-aligned null (DV = z_S, delta = 0.1, n = 180 rows, 20 dyads):
#   dyad random-intercept SD = 0.1152
#   residual SD              = 0.9090
#   fixed effects: intercept 2.2928, conditionvisible 0.1797,
#     periodP2 -2.3961, periodP3 -2.6832, cond:P2 -0.1209
#   observed P3 contrast (visible - invisible) = +0.829 (SE 0.240)
#
# Usage (zsh, conda env not required; base R with lme4/lmerTest):
#   Rscript R3-1_power_sensitivity.R --out_prefix R3-1_power \
#     --nsim 1000 --seed 20260823
# Optional overrides:
#   --grid "0.2,0.3,0.4,0.5,0.6,0.83,0.7,0.8,0.9,1.0"
#   --resid_sd 0.9090 --dyad_sd 0.1152 --alpha 0.05
#   --n_dyads_per_cond 10 --n_sessions 3
#
# Outputs:
#   <out_prefix>_curve.csv   : power per effect-size grid point
#   <out_prefix>_summary.txt : parameters, MDES (80% power), notes
# =====================================================================

suppressMessages({
  library(lme4)
  library(lmerTest)
})

## ---------------- argument parsing ----------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}

out_prefix       <- get_arg("--out_prefix", "R3-1_power")
nsim             <- as.integer(get_arg("--nsim", "1000"))
seed             <- as.integer(get_arg("--seed", "20260823"))
alpha            <- as.numeric(get_arg("--alpha", "0.05"))
resid_sd         <- as.numeric(get_arg("--resid_sd", "0.9090"))
dyad_sd          <- as.numeric(get_arg("--dyad_sd", "0.1152"))
n_dyads_per_cond <- as.integer(get_arg("--n_dyads_per_cond", "10"))
n_sessions       <- as.integer(get_arg("--n_sessions", "3"))
grid_str         <- get_arg("--grid",
                            "0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.829,0.9,1.0")
effect_grid <- sort(unique(as.numeric(strsplit(grid_str, ",")[[1]])))

## Fixed effects held at their observed values (delta = 0.1 grid-null fit).
## Nuisance terms only: in this balanced design they shift cell means but
## do not affect the sampling distribution of the P3 contrast test.
b_intercept <- 2.2928    # invisible, P1
d_P1        <- 0.1797    # visible - invisible contrast in P1 (observed)
b_P2        <- -2.3961
b_P3        <- -2.6832
d_P2        <- 0.1797 + (-0.1209)  # observed P2 contrast (= 0.0589)
## The P3 contrast (visible - invisible in P3) is the varied quantity.

set.seed(seed)

## ---------------- design skeleton ----------------
make_design <- function() {
  dyads <- data.frame(
    dyad_id   = factor(sprintf("d%02d", 1:(2 * n_dyads_per_cond))),
    condition = factor(rep(c("invisible", "visible"),
                           each = n_dyads_per_cond),
                       levels = c("invisible", "visible"))
  )
  dat <- expand.grid(dyad_id = dyads$dyad_id,
                     session = 1:n_sessions,
                     period  = factor(c("P1", "P2", "P3"),
                                      levels = c("P1", "P2", "P3")))
  merge(dat, dyads, by = "dyad_id")
}
design <- make_design()
stopifnot(nrow(design) == 2 * n_dyads_per_cond * n_sessions * 3)

## Contrast vectors on the fixed-effect scale
## fixef order: (Intercept), conditionvisible, periodP2, periodP3,
##              conditionvisible:periodP2, conditionvisible:periodP3
L_P1 <- c(0, 1, 0, 0, 0, 0)
L_P2 <- c(0, 1, 0, 0, 1, 0)
L_P3 <- c(0, 1, 0, 0, 0, 1)

simulate_once <- function(delta_P3) {
  d <- design
  is_vis <- d$condition == "visible"
  mu <- b_intercept +
    ifelse(d$period == "P2", b_P2, 0) +
    ifelse(d$period == "P3", b_P3, 0) +
    ifelse(is_vis & d$period == "P1", d_P1, 0) +
    ifelse(is_vis & d$period == "P2", d_P2, 0) +
    ifelse(is_vis & d$period == "P3", delta_P3, 0)
  b_d <- rnorm(nlevels(d$dyad_id), 0, dyad_sd)
  d$y <- mu + b_d[as.integer(d$dyad_id)] + rnorm(nrow(d), 0, resid_sd)

  fit <- suppressMessages(suppressWarnings(
    lmer(y ~ condition * period + (1 | dyad_id), data = d,
         control = lmerControl(check.conv.singular = "ignore"))
  ))
  p1 <- contest1D(fit, L_P1)[["Pr(>|t|)"]]
  p2 <- contest1D(fit, L_P2)[["Pr(>|t|)"]]
  p3 <- contest1D(fit, L_P3)[["Pr(>|t|)"]]
  p_holm3 <- p.adjust(c(p1, p2, p3), method = "holm")[3]
  c(p3 = p3, p3_holm = p_holm3)
}

## ---------------- run grid ----------------
cat(sprintf("Sensitivity power analysis: nsim = %d per grid point, %d grid points\n",
            nsim, length(effect_grid)))
res <- lapply(effect_grid, function(es) {
  pv <- t(vapply(seq_len(nsim), function(i) simulate_once(es),
                 c(p3 = 0, p3_holm = 0)))
  pw  <- mean(pv[, "p3"] < alpha)
  pwh <- mean(pv[, "p3_holm"] < alpha)
  cat(sprintf("  effect = %.3f : power (uncorrected) = %.3f | power (Holm) = %.3f\n",
              es, pw, pwh))
  data.frame(effect_size = es,
             power_uncorrected = pw,
             power_holm = pwh,
             mc_se = sqrt(pw * (1 - pw) / nsim),
             nsim = nsim)
})
res <- do.call(rbind, res)

## ---------------- MDES by linear interpolation ----------------
interp_mdes <- function(x, y, target = 0.80) {
  if (all(y < target)) return(NA_real_)
  if (y[1] >= target) return(x[1])
  i <- max(which(y < target))
  x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
}
mdes_unc  <- interp_mdes(res$effect_size, res$power_uncorrected)
mdes_holm <- interp_mdes(res$effect_size, res$power_holm)

## ---------------- outputs ----------------
csv_path <- paste0(out_prefix, "_curve.csv")
txt_path <- paste0(out_prefix, "_summary.txt")
write.csv(res, csv_path, row.names = FALSE)

sink(txt_path)
cat("=====================================================\n")
cat("R3-1 sensitivity power analysis:", format(Sys.time()), "\n")
cat("=====================================================\n\n")
cat("Design: 2 conditions x", n_dyads_per_cond, "dyads/condition x",
    n_sessions, "sessions x 3 periods\n")
cat("Model : z_S ~ condition * period + (1 | dyad_id)\n")
cat("Estimand: visible - invisible contrast in P3 (Satterthwaite t)\n\n")
cat("Parameters (from fitted delta = 0.1 z_S LMM):\n")
cat(sprintf("  residual SD = %.5f | dyad intercept SD = %.5f\n",
            resid_sd, dyad_sd))
cat(sprintf("  fixed effects: intercept %.5f, P1 contrast %.5f, P2 %.5f, P3 %.5f, P2 contrast %.5f\n",
            b_intercept, d_P1, b_P2, b_P3, d_P2))
cat(sprintf("  alpha = %.3f | nsim = %d | seed = %d\n\n", alpha, nsim, seed))
cat("Power curve:\n")
print(res, row.names = FALSE)
cat("\n")
cat(sprintf("MDES at 80%% power (uncorrected P3 contrast) : %.3f SD units of z_S\n",
            mdes_unc))
cat(sprintf("MDES at 80%% power (Holm across 3 periods)   : %.3f SD units of z_S\n",
            mdes_holm))
# [OBSERVED-ARG PATCH] value comes from the command line, not from a hard-coded string
observed    <- suppressWarnings(as.numeric(get_arg("--observed", NA)))
observed_se <- suppressWarnings(as.numeric(get_arg("--observed_se", NA)))
if (is.finite(observed)) {
  cat(sprintf("\nNote: observed P3 contrast (visible - invisible) in the re-run LMM = %+.3f%s.\n",
              observed, if (is.finite(observed_se)) sprintf(" (SE %.3f)", observed_se) else ""))
  cat("      (from ../R1_2/lmm_grid_out/lmm_results.txt via R3-1_power_params.txt)\n")
} else {
  cat("\nNote: observed P3 contrast not supplied (--observed); see R3-1_power_params.txt.\n")
}
cat("Post hoc observed power is intentionally NOT reported (circular).\n")
sink()

cat("\nSaved:", csv_path, "and", txt_path, "\n")
