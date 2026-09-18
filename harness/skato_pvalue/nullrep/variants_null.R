# Offline SKAT-O p-value variants on dumped cells (Score, Phi, r.corr per cell).
# References for the tail: (R1) Davies with an ADAPTIVE absolute tolerance,
# acc = 1e-3 * pilot (pilot = saddlepoint), lim 1e8, both layers, direct form --
# valid where every inner call converges (ifault 0), i.e. down to ~1e-14;
# (R2) Ruben's positive series where the spectrum's condition allows it.
A <- commandArgs(trailingOnly = TRUE); S <- A[1]; CD <- A[2]   # cells at file.path(S, CD, "cells.rds")
suppressPackageStartupMessages(library(SKAT)); ns <- asNamespace("SKAT")
for (f in c("Met_SKAT_Get_Pvalue","Get_Lambda","SKAT_davies","SKAT_META_Optimal_Param","SKAT_META_Optimal_Get_Q","Get_Liu_Params_Mod","Get_Liu_PVal.MOD.Lambda","SKAT_Optimal_PValue_Liu")) assign(f, get(f, ns))
cells <- readRDS(file.path(S, CD, "cells.rds")); cat("cells:", length(cells), "\n")
# ---- inner survival functions S(lambda, q) -> P(sum lambda_j chi2_1 > q). Return list(p, ok).
S_dav_skat <- function(lambda, q) {  # exactly what SKAT does inside Get_PValue.Lambda / the integrand: acc 1e-6, lim 1e4, value kept on ifault
  o <- SKAT_davies(q, lambda, lim = 1e4, acc = 1e-6); list(p = o$Qq, ok = o$ifault == 0) }
S_dav <- function(acc, lim) function(lambda, q) { if (q <= 0) return(list(p = 1, ok = TRUE)); o <- SKAT_davies(q, lambda, lim = lim, acc = acc); list(p = if (o$ifault == 0) o$Qq else NA_real_, ok = o$ifault == 0) }
spa_raw <- function(lambda, q) {
  if (q <= 0) return(1); mu <- sum(lambda); lmax <- max(lambda)
  if (abs(q - mu) < 1e-7 * mu) return(NA_real_)
  Kp <- function(z) sum(lambda / (1 - 2 * z * lambda)) - q
  if (q > mu) br <- c(0, (1 / (2 * lmax)) * (1 - 1e-12)) else { lo <- -1 / lmax; while (Kp(lo) > 0) lo <- lo * 2; br <- c(lo, 0) }
  zh <- uniroot(Kp, br, tol = 1e-12, maxiter = 5000)$root
  K <- -0.5 * sum(log(1 - 2 * zh * lambda)); Kpp <- sum(2 * lambda^2 / (1 - 2 * zh * lambda)^2)
  w <- sign(zh) * sqrt(2 * (zh * q - K)); v <- zh * sqrt(Kpp)
  if (!is.finite(w) || abs(w) < 1e-6) return(NA_real_)
  pnorm(w + log(v / w) / w, lower.tail = FALSE)
}
S_spa <- function(lambda, q) { p <- spa_raw(lambda, q); list(p = p, ok = !is.na(p)) }
# R1: Davies with acc from a saddlepoint pilot. Per-rho layer: strict (NA when Davies cannot deliver).
# Outer integrand: a node whose survival value is below 1e-3 of the cell's own pilot p contributes
# < 1e-3 relative to the integral, and takes the saddlepoint value instead of aborting the cell.
P_CELL <- 1
S_adapt_each <- function(lambda, q) {
  if (q <= 0) return(list(p = 1, ok = TRUE)); pil <- spa_raw(lambda, q); if (is.na(pil)) pil <- 0.5
  if (pil < 1e-17) return(list(p = NA_real_, ok = FALSE))
  o <- SKAT_davies(q, lambda, lim = 1e8, acc = max(min(pil * 1e-4, 1e-6), 1e-15)); list(p = if (o$ifault == 0 && o$Qq > 0) o$Qq else NA_real_, ok = o$ifault == 0 && o$Qq > 0) }
S_adapt_outer <- function(lambda, q) {
  if (q <= 0) return(list(p = 1, ok = TRUE)); pil <- spa_raw(lambda, q); if (is.na(pil)) pil <- 0.5
  if (pil < max(1e-2 * P_CELL, 1e-18)) return(list(p = pil, ok = TRUE))   # <= 1e-2 of the cell p at ~10% SPA error: < 0.1% on the cell
  o <- SKAT_davies(q, lambda, lim = 1e8, acc = max(min(pil * 1e-4, 1e-6), 1e-15)); list(p = if (o$ifault == 0 && o$Qq > 0) o$Qq else NA_real_, ok = o$ifault == 0 && o$Qq > 0) }
ruben_cache <- new.env(parent = emptyenv())
ruben_coef <- function(lambda, Kmax = 15000) {
  key <- paste(signif(lambda, 12), collapse = ","); if (!is.null(ruben_cache[[key]])) return(ruben_cache[[key]])
  beta <- 2 * min(lambda) * max(lambda) / (min(lambda) + max(lambda)); gam <- 1 - beta / lambda; gmax <- max(abs(gam))
  K <- ceiling(95 / (-log(gmax))) + 100; if (K > Kmax) return(NULL)
  g <- vapply(1:K, function(r) sum(gam^r), 0); a <- numeric(K + 1); a[1] <- prod(sqrt(beta / lambda))
  for (k in 1:K) a[k + 1] <- sum(g[1:k] * a[k:1]) / (2 * k)
  rc <- list(a = a, beta = beta, n = length(lambda), K = K); ruben_cache[[key]] <- rc; rc }
S_ruben <- function(lambda, q) { if (q <= 0) return(list(p = 1, ok = TRUE)); rc <- ruben_coef(lambda); if (is.null(rc)) return(list(p = NA_real_, ok = FALSE))
  list(p = sum(rc$a * pchisq(q / rc$beta, df = rc$n + 2 * (0:rc$K), lower.tail = FALSE)), ok = TRUE) }
# ---- the SKAT-O path (mirrors SKAT_META_Optimal -> _Get_Pvalue, method optimal.adj, no resampling), inner S pluggable.
# skat_semantics = TRUE reproduces SKAT: per-rho value kept whatever ifault says, Liu only if outside (0,1]; outer integral
# aborts to the Liu-form p if any inner call has ifault != 0. Otherwise: a failed inner call makes the cell NA (no silent fallback).
skato <- function(Score, Phi, r.corr, S_each, S_outer, direct = FALSE, skat_semantics = FALSE, rel.tol = 1e-8) {
  p.m <- nrow(Phi)
  if (sum(abs(Phi)) == 0) return(list(p = 1, note = "nopoly"))
  if (length(Phi) <= 1) r.corr <- 0 else if (ncol(Phi) <= 10 && qr(Phi)$rank <= 1) r.corr <- 0
  if (length(r.corr) == 1) return(list(p = NA_real_, note = "single-rho"))
  r.all <- r.corr; r.all[r.all >= 0.999] <- 0.999
  Q.r <- SKAT_META_Optimal_Get_Q(Score, r.all)$Q.r; Phi2 <- Phi / 2
  lambda.all <- lapply(r.all, function(rc) { R.M <- diag(rep(1 - rc, p.m)) + matrix(rc, p.m, p.m); L <- chol(R.M, pivot = TRUE); Get_Lambda(L %*% (Phi2 %*% t(L))) })
  param.m <- SKAT_META_Optimal_Param(Phi2, r.all)
  n.r <- length(r.all); pval <- numeric(n.r); pmat <- matrix(0, n.r, 3); n_noconv_kept <- 0L; n_liu <- 0L
  for (i in 1:n.r) {
    lam <- lambda.all[[i]]; c1 <- sapply(1:4, function(k) sum(lam^k)); pt <- Get_Liu_Params_Mod(c1)
    pl <- Get_Liu_PVal.MOD.Lambda(Q.r[i], lam)
    if (length(lam) == 1) p <- pl else {
      s <- S_each(lam, Q.r[i]); p <- s$p
      if (skat_semantics) { if (!s$ok && !is.na(p) && p > 0 && p <= 1) n_noconv_kept <- n_noconv_kept + 1L
        if (is.na(p) || p > 1 || p <= 0) { p <- pl; n_liu <- n_liu + 1L } }
      else if (is.na(p) || !s$ok || p > 1 || p <= 0) return(list(p = NA_real_, note = "each-fail"))
    }
    pval[i] <- p; pmat[i, ] <- c(pt$muQ, pt$sigmaQ^2, pt$l)
  }
  pmin <- min(pval)
  pmin.q <- sapply(1:n.r, function(i) { q.org <- qchisq(1 - pmin, df = pmat[i, 3]); (q.org - pmat[i, 3]) / sqrt(2 * pmat[i, 3]) * sqrt(pmat[i, 2]) + pmat[i, 1] })
  sd1 <- sqrt(param.m$VarQ - param.m$VarRemain) / sqrt(param.m$VarQ)
  Sfun <- function(x) { temp <- (pmin.q - param.m$tau %x% t(x)) / (1 - r.all); temp.min <- apply(temp, 2, min)
    vapply(temp.min, function(min1) { if (min1 > sum(param.m$lambda) * 1e4) return(0)
      s <- S_outer(param.m$lambda, (min1 - param.m$MuQ) * sd1 + param.m$MuQ)
      if (skat_semantics) { if (!s$ok) stop("dav.re$ifault is not 0") } else if (!s$ok || is.na(s$p)) stop("inner failure")
      min(s$p, 1) }, 0) }
  note <- ""
  if (!direct) {
    re <- try(integrate(function(x) (1 - Sfun(x)) * dchisq(x, df = 1), lower = 0, upper = 40, subdivisions = 1000, abs.tol = 1e-25), silent = TRUE)
    if (inherits(re, "try-error")) { if (skat_semantics) { pvalue <- SKAT_Optimal_PValue_Liu(pmin.q, param.m, r.all, pmin); note <- "liu-fallback" } else return(list(p = NA_real_, note = "outer-fail", pmin = pmin)) } else pvalue <- 1 - re$value
  } else {
    xstar <- max(min(pmin.q / param.m$tau), 0)
    # x = t^2 removes the 1/sqrt(x) singularity of dchisq(x, 1) at 0: dchisq(t^2, 1) * 2t = 2 * dnorm(t)
    re <- try(integrate(function(t) Sfun(t^2) * 2 * dnorm(t), lower = 0, upper = sqrt(xstar), subdivisions = 1000, rel.tol = rel.tol, abs.tol = 0), silent = TRUE)
    if (inherits(re, "try-error")) return(list(p = NA_real_, note = "outer-fail", pmin = pmin))
    pvalue <- pchisq(xstar, df = 1, lower.tail = FALSE) + re$value
  }
  p_unc <- pvalue; cap <- 0L
  if (pmin * n.r < pvalue) { pvalue <- pmin * n.r; cap <- 1L }
  IDX <- which(pval > 0); bonf <- 0L
  if (pvalue <= 0 || length(IDX) < n.r) { pvalue <- min(pval) * 3; bonf <- 1L }
  if (pvalue == 0 && length(IDX) > 0) { pvalue <- min(pval[IDX]); bonf <- 2L }
  list(p = pvalue, p_unc = p_unc, pmin = pmin, cap = cap, bonf = bonf, n_noconv_kept = n_noconv_kept, n_liu = n_liu, note = note, cond = max(param.m$lambda) / min(param.m$lambda), cond_each = max(sapply(lambda.all, function(l) max(l) / min(l))))
}
tm <- function(expr) { t0 <- proc.time()[[3]]; v <- try(expr, silent = TRUE); t <- proc.time()[[3]] - t0; if (inherits(v, "try-error")) list(p = NA_real_, note = paste("ERR", conditionMessage(attr(v, "condition"))), t = t) else c(v, t = t) }
g <- function(v, k) { x <- v[[k]]; if (is.null(x) || length(x) == 0) NA else x }
S_tight <- S_dav(1e-9, 1e5)
rows <- list()
for (i in seq_along(cells)) {
  ce <- cells[[i]]; if (length(ce$r.corr) == 1) next
  li <- tm(Met_SKAT_Get_Pvalue(Score = ce$Score, Phi = ce$Phi, r.corr = ce$r.corr, method = "liu.mod", Score.Resampling = NULL))
  row <- list(i = i, q = length(ce$Score), p_skat = ce$p, t_skat = ce$t, p_liu = g(li, "p.value"), t_liu = g(li, "t"), deep = 0L)
  if (!is.na(row$p_skat) && row$p_skat < 1e-2) {
    row$deep <- 1L
    sp <- tm(skato(ce$Score, ce$Phi, ce$r.corr, S_spa, S_spa)); sd <- tm(skato(ce$Score, ce$Phi, ce$r.corr, S_spa, S_spa, direct = TRUE, rel.tol = 1e-6))
    P_CELL <<- if (!is.na(g(sd, "p"))) g(sd, "p") else row$p_skat
    r1 <- tm(skato(ce$Score, ce$Phi, ce$r.corr, S_adapt_each, S_adapt_outer, direct = TRUE, rel.tol = 1e-4))
    d <- tm(skato(ce$Score, ce$Phi, ce$r.corr, S_dav_skat, S_dav_skat, skat_semantics = TRUE))
    row <- c(row, list(p_spa = g(sp, "p"), p_spa_direct = g(sd, "p"), sd_note = g(sd, "note"), t_spa_direct = g(sd, "t"), p_ref1 = g(r1, "p"), ref1_note = g(r1, "note"), p_def = g(d, "p"), def_cap = g(d, "cap"), def_liu_each = g(d, "n_liu")))
  }
  rows[[length(rows) + 1]] <- row
}
allk <- unique(unlist(lapply(rows, names)))
df <- as.data.frame(do.call(rbind, lapply(rows, function(r) { r <- r[allk]; names(r) <- allk; lapply(r, function(x) if (is.null(x)) NA else x) })), stringsAsFactors = FALSE)
for (k in names(df)) { v <- unlist(lapply(df[[k]], function(x) if (is.null(x)) NA else x)); df[[k]] <- if (all(is.na(v)) || is.character(v)) v else as.numeric(v) }
saveRDS(df, sub("cells.rds$", "variants.rds", file.path(S, CD, "cells.rds"))); cat("wrote", nrow(df), "rows\n")
