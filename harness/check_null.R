#!/usr/bin/env Rscript
# Negative-control check: p-values from a permuted-phenotype run must be
# consistent with Uniform(0,1). Halts (exit 1) if the null misbehaves.
# Usage: Rscript check_null.R <step2_output.txt>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
d <- read.table(args[1], header = TRUE, stringsAsFactors = FALSE)
stopifnot("p.value" %in% names(d))
p <- as.numeric(d$p.value)
p <- p[is.finite(p)]
n <- length(p)

ks <- suppressWarnings(ks.test(p, "punif"))
chisq <- qchisq(1 - p, df = 1)
lambda <- median(chisq) / qchisq(0.5, df = 1)
minp <- min(p)
# Bonferroni-style floor: under the null, min p ~ Beta(1, n)
minp_crit <- 0.001 / n

cat(sprintf("n markers        : %d\n", n))
cat(sprintf("KS vs uniform    : D=%.3f p=%.3f  (need p > 0.01)\n", ks$statistic, ks$p.value))
cat(sprintf("lambda_GC        : %.3f            (informational; n is small)\n", lambda))
cat(sprintf("min p-value      : %.2e         (flag if < %.1e)\n", minp, minp_crit))

fail <- (ks$p.value <= 0.01) || (minp < minp_crit)
if (fail) { cat("== NEGATIVE CONTROL: FAIL — null run shows non-uniform p-values\n"); quit(status = 1) }
cat("== NEGATIVE CONTROL: PASS\n")
