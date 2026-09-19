# Validate the OUTER rho integral of the saddlepoint tail (the half
# handoff-2026-09-19-saddlepoint-tail-validation.md left open). The inner
# survival function is validated; this asks how its 3.6-12% conservative,
# depth-driven bias COMPOSES through the direct piecewise outer integral.
#
# Method: run the shipped skatoSaddlepointTail against an identical copy whose
# ONLY difference is that every survival call goes to the cancellation-free
# Ruben reference instead of Kuonen's saddlepoint. Same envelope, same knots,
# same quadrature, same post-rules -- so any difference is the survival function
# propagating, which is exactly the question.
suppressPackageStartupMessages({library(SAIGE); library(SKAT)})
ns <- asNamespace("SKAT")
GetLambda <- get("Get_Lambda", ns); MetaParam <- get("SKAT_META_Optimal_Param", ns)
MetaGetQ <- get("SKAT_META_Optimal_Get_Q", ns); GetPV <- get("Get_PValue.Lambda", ns)
LiuMod <- get("Get_Liu_Params_Mod", ns)

# Ruben with the coefficient recursion CACHED per lambda: a_k depends only on
# lambda, not on q, so the outer integrand's many nodes share one recursion.
rubenFactory <- function(lambda, kmax = 3000L) {
  n <- length(lambda); beta <- min(lambda) * (1 - 1e-13)
  gam <- 1 - beta / lambda
  a <- numeric(kmax + 1L); g <- numeric(kmax + 1L)
  a[1] <- exp(0.5 * sum(log(beta / lambda))); asum <- a[1]; pw <- rep(1, n); K <- 0L
  for (k in seq_len(kmax)) {
    pw <- pw * gam; g[k + 1L] <- sum(pw)
    a[k + 1L] <- sum(g[2:(k + 1L)] * a[k:1]) / (2 * k)
    asum <- asum + a[k + 1L]; K <- k
    if (k %% 32L == 0L && (1 - asum) < 1e-15) break
  }
  ak <- a[seq_len(K + 1L)]; dfs <- n + 2 * (0:K)
  function(q) {
    if (q <= 0) return(1)
    sum(ak * pchisq(q / beta, df = dfs, lower.tail = FALSE))
  }
}

# A byte-for-byte copy of skatoSaddlepointTail with the survival swapped.
mkRubenTail <- function() {
  src <- deparse(SAIGE:::skatoSaddlepointTail)
  src <- gsub("skatoSaddlepointSurvival(param.m$lambda, q)", "outerS(q)", src, fixed = TRUE)
  f <- eval(parse(text = paste(src, collapse = "\n")))
  environment(f) <- new.env(parent = environment(SAIGE:::skatoSaddlepointTail))
  f
}
rubenTailTemplate <- mkRubenTail()

cellTail <- function(Score, Phi, r.corr) {
  p.m <- nrow(Phi)
  if (length(r.corr) <= 1 || sum(abs(Phi)) == 0 || length(Phi) <= 1 ||
        (ncol(Phi) <= 10 && qr(Phi)$rank <= 1)) return(NULL)
  r.all <- r.corr; r.all[r.all >= 0.999] <- 0.999; n.r <- length(r.all)
  Qr <- MetaGetQ(Score, r.all)$Q.r; Phi2 <- Phi / 2
  lambda.all <- vector("list", n.r); param.mat <- matrix(0, n.r, 3); pvalDavies <- numeric(n.r)
  for (i in seq_len(n.r)) {
    RM <- diag(rep(1 - r.all[i], p.m)) + matrix(rep(r.all[i], p.m * p.m), ncol = p.m)
    L <- chol(RM, pivot = TRUE)
    lam <- GetLambda(L %*% (Phi2 %*% t(L))); lambda.all[[i]] <- lam
    c1 <- c(sum(lam), sum(lam^2), sum(lam^3), sum(lam^4)); pt <- LiuMod(c1)
    param.mat[i, ] <- c(pt$muQ, pt$sigmaQ^2, pt$l)
    pvalDavies[i] <- GetPV(lam, Qr[i])$p.value
  }
  param.m <- MetaParam(Phi2, r.all)
  shipped <- SAIGE:::skatoSaddlepointTail(Qr, lambda.all, param.m, r.all, pvalDavies, param.mat)
  env <- environment(rubenTailTemplate)
  env$outerS <- rubenFactory(param.m$lambda)
  ref <- rubenTailTemplate(Qr, lambda.all, param.m, r.all, pvalDavies, param.mat)
  if (is.null(shipped) || is.null(ref)) return(NULL)
  list(shipped = shipped$p.value, ref = ref$p.value, q = length(Score),
       davies = min(pvalDavies))
}
