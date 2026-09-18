# Refined screen: per-rho layer = SKAT's Davies (unchanged, so Pvalue_SKAT / Pvalue_Burden / SE_Burden are byte-identical),
# outer integral = Liu form first; Davies outer only when the Liu omnibus p <= threshold.
S <- commandArgs(trailingOnly = TRUE)[1]; CD <- commandArgs(trailingOnly = TRUE)[2]
src <- readLines("/Users/dpalmer/Repositories/SAIGE/harness/skato_pvalue/variants.R"); cut <- grep("^set.seed", src); A <- c(S, CD); eval(parse(text = src[4:(cut - 1)]))
layers <- function(Score, Phi, r.corr) {
  p.m <- nrow(Phi); if (length(Phi) <= 1 || (ncol(Phi) <= 10 && qr(Phi)$rank <= 1)) return(NULL)
  r.all <- r.corr; r.all[r.all >= 0.999] <- 0.999
  t0 <- proc.time()[[3]]
  Q.r <- SKAT_META_Optimal_Get_Q(Score, r.all)$Q.r; Phi2 <- Phi / 2
  lambda.all <- lapply(r.all, function(rc) { R.M <- diag(rep(1 - rc, p.m)) + matrix(rc, p.m, p.m); L <- chol(R.M, pivot = TRUE); Get_Lambda(L %*% (Phi2 %*% t(L))) })
  param.m <- SKAT_META_Optimal_Param(Phi2, r.all); n.r <- length(r.all); pval <- numeric(n.r); pmat <- matrix(0, n.r, 3); pliu_each <- numeric(n.r)
  for (i in 1:n.r) { lam <- lambda.all[[i]]; c1 <- sapply(1:4, function(k) sum(lam^k)); pt <- Get_Liu_Params_Mod(c1); pl <- Get_Liu_PVal.MOD.Lambda(Q.r[i], lam); pliu_each[i] <- pl
    if (length(lam) == 1) p <- pl else { s <- S_dav_skat(lam, Q.r[i]); p <- s$p; if (is.na(p) || p > 1 || p <= 0) p <- pl }
    pval[i] <- p; pmat[i, ] <- c(pt$muQ, pt$sigmaQ^2, pt$l) }
  pmin <- min(pval); pmin.q <- sapply(1:n.r, function(i) { q.org <- qchisq(1 - pmin, df = pmat[i, 3]); (q.org - pmat[i, 3]) / sqrt(2 * pmat[i, 3]) * sqrt(pmat[i, 2]) + pmat[i, 1] })
  t1 <- proc.time()[[3]]
  p_liu_outer <- SKAT_Optimal_PValue_Liu(pmin.q, param.m, r.all, pmin)      # includes the pmin*8 cap
  t2 <- proc.time()[[3]]
  p_dav_outer <- get("SKAT_Optimal_PValue_Davies", ns)(pmin.q, param.m, r.all, pmin)
  t3 <- proc.time()[[3]]
  list(t_each = t1 - t0, t_liu_outer = t2 - t1, t_dav_outer = t3 - t2, p_liu_outer = p_liu_outer, p_dav_outer = p_dav_outer, p_each_rho0 = pval[1], p_each_rho1 = pval[n.r], pliu_rho0 = pliu_each[1], pliu_rho1 = pliu_each[n.r])
}
rows <- list()
for (ce in cells) { if (length(ce$Score) == 1) next; l <- layers(ce$Score, ce$Phi, ce$r.corr); if (is.null(l)) next; rows[[length(rows) + 1]] <- c(l, p_skat = ce$p, q = length(ce$Score)) }
d <- as.data.frame(do.call(rbind, lapply(rows, unlist)))
cat("multi-rho cells:", nrow(d), "\n")
cat(sprintf("time per cell, ms (mean): per-rho Davies layer %.1f | Liu outer %.2f | Davies outer %.1f  -> outer is %.0f%% of the kernel\n", 1000*mean(d$t_each), 1000*mean(d$t_liu_outer), 1000*mean(d$t_dav_outer), 100*sum(d$t_dav_outer)/sum(d$t_each + d$t_dav_outer)))
cat("control: Davies-outer reproduces SKAT p exactly on", sum(d$p_dav_outer == d$p_skat), "of", nrow(d), "(max rel diff", signif(max(abs(d$p_dav_outer/d$p_skat - 1)), 2), ")\n")
for (th in c(0.05, 0.01, 0.001)) cat(sprintf("screen on Liu-OUTER p > %g: skip %.1f%% of cells, %.1f%% of Davies-outer time; misses (Davies p < %g among skipped): %d; worst skipped Davies p %.3g\n", th, 100*mean(d$p_liu_outer > th), 100*sum(d$t_dav_outer[d$p_liu_outer > th])/sum(d$t_dav_outer), th/10, sum(d$p_liu_outer > th & d$p_skat < th/10), min(d$p_skat[d$p_liu_outer > th])))
r <- d$p_liu_outer / d$p_skat; b <- cut(d$p_skat, c(0, 1e-3, 1e-2, 0.05, 1))
cat("Liu-outer / Davies omnibus by Davies p band (median, min, max):\n"); print(do.call(rbind, tapply(r, b, function(x) round(c(n = length(x), median = median(x), min = min(x), max = max(x)), 3))))
cat("per-rho: Liu vs Davies at rho=0 (SKAT component) ratio range", round(range(d$pliu_rho0 / d$p_each_rho0), 3), "; at rho=1 (burden) range", round(range(d$pliu_rho1 / d$p_each_rho1), 4), "\n")
