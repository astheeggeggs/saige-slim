#!/usr/bin/env Rscript
# Negative-control check: p-values from a permuted-phenotype run must be
# consistent with Uniform(0,1). Halts (exit 1) if the null misbehaves.
# Usage: Rscript check_null.R <step2_output.txt>
#
# MAC STRATIFICATION IS NOT OPTIONAL FOR A BINARY TRAIT. SAIGE routes markers
# with MAC <= max_MAC_for_ER (default 4) to the exact ER test, whose p-value
# distribution is ATOMIC -- a singleton against K*N cases can only take a
# handful of values. A KS test against a CONTINUOUS uniform is not a valid test
# of those p-values and will reject a perfectly calibrated run once there are
# enough markers for it to have power. The bundled fixture has 70 markers, so
# the KS test has almost none and this passed for the wrong reason; at 200,000
# markers it would fail a correct run outright.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
d <- read.table(args[1], header = TRUE, stringsAsFactors = FALSE)
stopifnot("p.value" %in% names(d))

# MAC is not an output column; derive it from the allele count and the sample
# size, which are (under either trait type).
mac <- NULL
if ("MAC" %in% names(d)) {
  mac <- as.numeric(d$MAC)
} else if ("AC_Allele2" %in% names(d)) {
  n <- if (all(c("N_case", "N_ctrl") %in% names(d))) d$N_case + d$N_ctrl else
       if ("N" %in% names(d)) d$N else NA
  if (!all(is.na(n))) mac <- pmin(as.numeric(d$AC_Allele2), 2 * n - as.numeric(d$AC_Allele2))
}

p_all <- as.numeric(d$p.value)
keep <- is.finite(p_all)
p_all <- p_all[keep]
if (!is.null(mac)) mac <- mac[keep]
n <- length(p_all)
lam <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (!length(p)) return(NA_real_)
  median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

cat(sprintf("n markers        : %d\n", n))

# The min-p floor is valid on every marker regardless of discreteness.
minp <- min(p_all)
minp_crit <- 0.001 / n
cat(sprintf("min p-value      : %.2e         (flag if < %.1e)\n", minp, minp_crit))

ks_p <- NA_real_; ks_d <- NA_real_; gated_on <- "all markers"
if (!is.null(mac)) {
  br <- list(c(0, 4), c(5, 19), c(20, Inf))
  nm <- c("MAC<=4 (ER, atomic)", "MAC 5-19", "MAC>=20 (near-continuous)")
  cat("stratum                     n   lambda_GC    min p\n")
  for (i in seq_along(br)) {
    s <- p_all[mac >= br[[i]][1] & mac <= br[[i]][2]]
    cat(sprintf("%-24s %5d  %9s %8.2e\n", nm[i], length(s),
                if (length(s)) sprintf("%.3f", lam(s)) else "-",
                if (length(s)) min(s) else NA))
  }
  cont <- p_all[mac >= 20]
  if (length(cont) >= 30) {
    gated_on <- sprintf("MAC>=20 (%d markers)", length(cont))
    k <- suppressWarnings(ks.test(cont, "punif")); ks_p <- k$p.value; ks_d <- k$statistic
  } else {
    cat(sprintf("NOTE: only %d markers with MAC>=20; too few to gate on, falling back to all.\n",
                length(cont)))
  }
}
if (is.na(ks_p)) {
  k <- suppressWarnings(ks.test(p_all, "punif")); ks_p <- k$p.value; ks_d <- k$statistic
}
cat(sprintf("lambda_GC (all)  : %.3f            (informational; n may be small)\n", lam(p_all)))
cat(sprintf("KS vs uniform    : D=%.3f p=%.3f  on %s  (need p > 0.01)\n",
            ks_d, ks_p, gated_on))

fail <- (ks_p <= 0.01) || (minp < minp_crit)
if (fail) { cat("== NEGATIVE CONTROL: FAIL - null run shows non-uniform p-values\n"); quit(status = 1) }
cat("== NEGATIVE CONTROL: PASS\n")
