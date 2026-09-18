# Pooled null calibration per arm: ECDF at alpha with binomial CI for SKAT as shipped, modified Liu,
# the saddlepoint (direct form) and the adaptive-Davies reference; screen misses; per-cell agreement in the tail.
W <- commandArgs(trailingOnly = TRUE)[1]
arms <- list(quantitative = "y_null_r", binary_K0.01 = "y_nb01_r", binary_K0.002 = "y_nb002_r")
ci <- function(k, n) { b <- binom.test(k, n)$conf.int; sprintf("%d (%.1e, %.1e-%.1e)", k, k / n, b[1], b[2]) }
for (arm in names(arms)) {
  fs <- Sys.glob(file.path(W, "runs", paste0(arms[[arm]], "*.variants.rds"))); if (length(fs) == 0) next
  d <- do.call(rbind, lapply(fs, function(f) { x <- readRDS(f); x$rep <- sub(".*_r([0-9]+)\\.variants\\.rds$", "\\1", f); x }))
  d <- d[!is.na(d$p_skat), ]; n <- nrow(d)
  cat(sprintf("\n===== %s: %d replicates, %d cells (multi-rho SKAT-O cells; single-rho cells excluded)\n", arm, length(fs), n))
  spa <- ifelse(d$deep == 1 & !is.na(d$p_spa_direct), d$p_spa_direct, d$p_skat)   # above 1e-2 the saddlepoint cannot reach the alphas below
  ref <- ifelse(d$deep == 1 & !is.na(d$p_ref1), d$p_ref1, d$p_skat)
  cat("alpha      expected   SKAT shipped                Liu                          saddlepoint-direct           adaptive-Davies ref\n")
  for (a in c(0.05, 0.01, 1e-3, 1e-4, 1e-5)) cat(sprintf("%-8g %8.1f   %-28s %-28s %-28s %s\n", a, n * a, ci(sum(d$p_skat < a), n), ci(sum(d$p_liu < a, na.rm = TRUE), n), ci(sum(spa < a), n), ci(sum(ref < a), n)))
  cat("reference available on", sum(!is.na(d$p_ref1)), "of", sum(d$deep == 1), "cells below 1e-2; saddlepoint-direct on", sum(!is.na(d$p_spa_direct)), "\n")
  cat("screen misses: Liu p > 0.01 with SKAT p < 1e-3:", sum(d$p_liu > 0.01 & d$p_skat < 1e-3, na.rm = TRUE), "; with reference p < 1e-3:", sum(d$p_liu > 0.01 & !is.na(d$p_ref1) & d$p_ref1 < 1e-3, na.rm = TRUE), "; cells passing the screen:", sum(d$p_liu <= 0.01, na.rm = TRUE), sprintf("(%.1f%%)", 100 * mean(d$p_liu <= 0.01, na.rm = TRUE)), "\n")
  E <- d[!is.na(d$p_ref1), ]
  if (nrow(E)) { r <- function(p) round(quantile(p / E$p_ref1, c(0, .05, .5, .95, 1), na.rm = TRUE), 3)
    cat("ratio to reference on the referenced null cells (min, 5%, median, 95%, max): SKAT", r(E$p_skat), "; Liu", r(E$p_liu), "; saddlepoint-direct", r(E$p_spa_direct), "\n")
    cat("smallest reference p in this arm:", signif(min(E$p_ref1), 3), "; smallest SKAT p:", signif(min(d$p_skat), 3), "; smallest Liu p:", signif(min(d$p_liu, na.rm = TRUE), 3), "\n") }
  cat("per-replicate SKAT count below 1e-3 (expected", round(n / length(fs) * 1e-3, 1), "each):", paste(tapply(d$p_skat < 1e-3, d$rep, sum), collapse = " "), "\n")
}
