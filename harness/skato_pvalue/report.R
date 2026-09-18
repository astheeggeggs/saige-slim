A <- commandArgs(trailingOnly = TRUE); S <- A[1]; CD <- if (length(A) > 1) A[2] else "cells"
d <- readRDS(file.path(S, CD, "variants.rds")); d <- d[!is.na(d$p_def), ]   # multi-rho cells only
cat("multi-rho cells:", nrow(d), " deep:", sum(d$deep == 1), "\n")
cat("CONTROL reimplementation == SKAT exactly:", sum(d$p_def == d$p_skat), "of", nrow(d), "; max |rel diff|", signif(max(abs(d$p_def / d$p_skat - 1)), 3), "\n")
cat("SKAT branches: outer Liu fallback", sum(d$def_note == "liu-fallback"), "; pmin*8 cap", sum(d$def_cap == 1), "; p<=0 -> 3*min", sum(d$def_bonf == 1), "\n")
cat("per-rho layer: cells keeping >=1 NON-CONVERGED Davies value", sum(d$def_noconv_kept > 0), "; cells with >=1 Liu substitution", sum(d$def_liu_each > 0), "\n")
d$p_ref <- d$p_ref1; d$ref_src <- ifelse(!is.na(d$p_ref1), "davies-adaptive", NA)
att <- !is.na(d$ref1_note) | !is.na(d$p_ref1)
cat("\nreference (adaptive Davies, direct form, both layers strict): attempted on", sum(att), "cells (SKAT p < 1e-2); available on", sum(!is.na(d$p_ref1)), "; failed:", paste(names(table(d$ref1_note[att & is.na(d$p_ref1)])), table(d$ref1_note[att & is.na(d$p_ref1)]), collapse = ", "), "\n")
cat("SKAT p of the cells WITHOUT a reference:", paste(signif(sort(unique(d$p_skat[att & is.na(d$p_ref1)])), 2)[1:min(12, sum(att & is.na(d$p_ref1)))], collapse = " "), "\n")
E <- d[!is.na(d$p_ref), ]; ref <- E$p_ref
cat("reference p range:", signif(min(ref), 3), "to", signif(max(ref), 3), "\n")
band <- cut(ref, c(0, 1e-12, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2), right = TRUE)
rat <- function(p) { r <- p / ref; out <- do.call(rbind, lapply(split(seq_along(r), band), function(i) if (length(i) == 0) c(n = 0, fails = 0, median = NA, min = NA, max = NA) else c(n = length(i), fails = sum(is.na(r[i])), median = round(median(r[i], na.rm = TRUE), 3), min = round(suppressWarnings(min(r[i], na.rm = TRUE)), 3), max = round(suppressWarnings(max(r[i], na.rm = TRUE)), 3)))); out[out[, "n"] > 0, , drop = FALSE] }
reg <- ifelse(E$p_skat >= 1e-5, E$p_skat, ifelse(!is.na(E$p_spa), E$p_spa, E$p_tight))
regd <- ifelse(E$p_skat >= 1e-5, E$p_skat, ifelse(!is.na(E$p_spa_direct), E$p_spa_direct, E$p_tight))
wu <- ifelse(!is.na(E$p_tight), E$p_tight, E$p_spa_direct)
V <- list("SKAT as shipped" = E$p_skat, "modified Liu" = E$p_liu, "SKAT inner Davies, direct form" = E$p_def_direct, "Davies acc1e-9 lim1e5, 1-int" = E$p_tight, "saddlepoint, 1-int" = E$p_spa, "saddlepoint, direct" = E$p_spa_direct, "REGENIE tier (skat if >=1e-5, else spa, else tight)" = reg, "REGENIE tier with direct forms" = regd, "Wu-like (tight Davies, spa fallback)" = wu)
cat("\n== ratio to the reference, by reference p band (n, fails = NA, median, min, max) ==\n")
for (nm in names(V)) { cat("--", nm, "\n"); print(rat(V[[nm]])) }
cat("\n== per-rho layer alone: SKAT's min-over-rho p / reference min-over-rho p (n, median, min, max) ==\n"); rr <- E$pmin_def / E$pmin_ref1; print(do.call(rbind, lapply(split(seq_along(rr), band), function(i) if (length(i) == 0) NULL else c(n = length(i), median = round(median(rr[i]), 3), min = round(min(rr[i]), 3), max = round(max(rr[i]), 3)))))
cat("saddlepoint per-rho min / reference min:"); rr2 <- E$pmin_spa / E$pmin_ref1; cat(" median", round(median(rr2, na.rm=TRUE), 3), " range", round(range(rr2, na.rm=TRUE), 3), "\n")
cat("\n== timing, ms per cell, mean over deep cells ==\n"); D <- d[d$deep == 1, ]
print(round(1000 * colMeans(D[, c("t_skat", "t_def", "t_liu", "t_spa", "t_spa_direct", "t_tight")], na.rm = TRUE), 1))
cat("reference cost, s per cell:", round(mean(E$t_ref1, na.rm = TRUE), 2), "\n")
cat("\n== whole run: SKAT kernel", round(sum(d$t_skat), 1), "s; Liu", round(sum(d$t_liu), 1), "s; cells with SKAT p<1e-5:", sum(d$p_skat < 1e-5), " (spa cost there", round(sum(d$t_spa[d$p_skat < 1e-5], na.rm = TRUE), 1), "s)\n")
cat("Liu screen at 0.01: skip", sum(d$p_liu > 0.01), "cells; misses with reference p < 1e-3:", sum(d$p_liu > 0.01 & !is.na(d$p_ref) & d$p_ref < 1e-3), "; with SKAT p < 1e-3:", sum(d$p_liu > 0.01 & d$p_skat < 1e-3), "\n")
cat("\n== cells with reference < 1e-5, unique ==\n")
x <- unique(E[E$p_ref < 1e-5, c("q", "cond", "p_skat", "def_cap", "def_bonf", "def_noconv_kept", "p_liu", "p_spa_direct", "p_tight", "p_ref", "ref_src")]); x <- x[order(x$p_ref), ]
print(data.frame(q = x$q, cond = signif(x$cond, 2), skat = signif(x$p_skat, 3), br = ifelse(x$def_bonf > 0, "3min", ifelse(x$def_cap == 1, "cap8", "int")), kept = x$def_noconv_kept, liu = signif(x$p_liu, 3), spa_d = signif(x$p_spa_direct, 3), tight = signif(x$p_tight, 3), ref = signif(x$p_ref, 3), src = x$ref_src, r_skat = round(x$p_skat / x$p_ref, 2), r_spa = round(x$p_spa_direct / x$p_ref, 3), r_liu = round(x$p_liu / x$p_ref, 2))[1:min(40, nrow(x)), ], row.names = FALSE)
