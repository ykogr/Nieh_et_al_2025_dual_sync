#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# R2-4 / discretization check, step 2:
# Fit the manuscript's LMM (zS ~ condition * period + (1 | dyad)) to the
# per-video zS values produced by R2-4_zS_nullcheck.py, separately for each
# tolerance delta and each null convention (continuous = published pipeline,
# grid = corrected), and report the visible - invisible contrast per period
# (Kenward-Roger df, as in the manuscript; emmeans).
#
# This answers: does the Table-1 condition contrast survive the corrected
# (grid-aligned) null convention?
#
# NOTE: this uses the independent reimplementation's zS values (n_perm as in
# the nullcheck run), not the original pipeline's outputs.  It is a fast
# directional check; final citable corrected numbers should come from
# re-running the original pipeline with a grid-aligned shift.
#
# Run (Mac / zsh, conda env nieh_iscience_2026, from iscience_revision_analysis/R2_4_3_1):
#
#   Rscript R2-4_nullcheck_lmm.R --per_video R2-4_nullcheck_per_video.csv \
#       --out_prefix R2-4_nullcheck_lmm
#
# Outputs:
#   <out_prefix>_contrasts.csv : delta x null x period: V-I estimate, SE,
#                                df, t, p (unadjusted) + Holm across periods
#   <out_prefix>_summary.txt   : the same as a readable table
# ---------------------------------------------------------------------------

suppressMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}
per_video_csv <- get_arg("--per_video", "R2-4_nullcheck_per_video.csv")
out_prefix    <- get_arg("--out_prefix", "R2-4_nullcheck_lmm")

stopifnot(file.exists(per_video_csv))
d <- read.csv(per_video_csv, stringsAsFactors = FALSE)

# dyad id = video name minus the trailing session token
# (pair videos are <date>_<pA>_<pB>_<session>)
d$dyad <- sub("_[^_]+$", "", d$video)
d$condition <- factor(d$condition, levels = c("invisible", "visible"))
d$period <- factor(d$period)
d <- d[is.finite(d$zS), ]

cat(sprintf("Loaded %d rows | %d videos | %d dyads\n",
            nrow(d), length(unique(d$video)), length(unique(d$dyad))))
cat("Dyads:", paste(sort(unique(d$dyad)), collapse = " "), "\n\n")

res <- NULL
for (delta_i in sort(unique(d$delta))) {
  for (nul in sort(unique(d$null))) {
    dd <- droplevels(d[d$delta == delta_i & d$null == nul, ])
    if (nrow(dd) < 10) next
    fit <- tryCatch(
      lmer(zS ~ condition * period + (1 | dyad), data = dd),
      error = function(e) { message("lmer failed: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    emm <- emmeans(fit, ~ condition | period, lmer.df = "kenward-roger")
    ct  <- as.data.frame(pairs(emm, reverse = TRUE))  # visible - invisible
    ct$p_holm <- p.adjust(ct$p.value, method = "holm")
    ct$delta <- delta_i
    ct$null  <- nul
    res <- rbind(res, ct)
  }
}

res <- res[, c("delta", "null", "period", "estimate", "SE", "df",
               "t.ratio", "p.value", "p_holm")]
names(res) <- c("delta", "null", "period", "V_minus_I", "SE", "df",
                "t", "p", "p_holm")

write.csv(res, paste0(out_prefix, "_contrasts.csv"), row.names = FALSE)

sink(paste0(out_prefix, "_summary.txt"))
cat("R2-4 nullcheck LMM: zS ~ condition * period + (1 | dyad)\n")
cat("Contrast: visible - invisible per period (Kenward-Roger df)\n")
cat("Input:", per_video_csv, "|", nrow(d), "rows |",
    length(unique(d$dyad)), "dyads\n")
cat("NOTE: reimplementation zS, directional check only; final corrected\n")
cat("numbers require re-running the original pipeline with grid shifts.\n\n")
for (nul in unique(res$null)) {
  cat("=====", nul, "null =====\n")
  sub <- res[res$null == nul, ]
  sub <- sub[order(sub$delta, sub$period), ]
  cat(sprintf("%6s %4s %10s %7s %6s %7s %8s %8s\n",
              "delta", "per", "V-I", "SE", "df", "t", "p", "p_holm"))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("%6.1f %4s %+10.3f %7.3f %6.1f %+7.2f %8.4f %8.4f\n",
                sub$delta[i], as.character(sub$period[i]),
                sub$V_minus_I[i], sub$SE[i], sub$df[i],
                sub$t[i], sub$p[i], sub$p_holm[i]))
  }
  cat("\n")
}
sink()

cat(readLines(paste0(out_prefix, "_summary.txt")), sep = "\n")
cat("\nWrote:", paste0(out_prefix, "_contrasts.csv"), "\n")
cat("Wrote:", paste0(out_prefix, "_summary.txt"), "\n")
