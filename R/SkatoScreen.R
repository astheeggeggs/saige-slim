# LEDGER #158: the SKAT-O kernel's cost is almost entirely the outer rho
# integral (95% of the kernel on sim50k: 4.4 ms for the eight per-rho Davies
# calls against 80 ms for the integral, with Davies' inversion at every
# quadrature node). SKAT's own modified-Liu form of that integral
# (SKAT_Optimal_PValue_Liu, a closed form per node) costs 5 ms and, above a
# threshold, sits within 10% of the Davies value (0.95-1.10 on 3,151 cells of
# two traits; median 1.00 above 0.05, 1.02 in (0.01, 0.05]). Below the
# threshold it is anticonservative by 3-7x, which is why it is a SCREEN and
# never the answer there.
#
# The design (docs/state/handoff-2026-09-17-skato-pvalue-methods.md, section
# 10) is therefore: the per-rho layer exactly as SKAT computes it -- Davies at
# acc 1e-6, Liu only when the value leaves (0, 1] -- so `Pvalue_SKAT`,
# `Pvalue_Burden`, `BETA_Burden` and `SE_Burden` (the meta-analysis's signed
# burden Z is back-derived from that burden p) are byte-identical for every
# cell of a run; then the Liu-form outer integral first, reported as `Pvalue`
# when it is above the threshold, and SKAT's Davies path, verbatim, otherwise.
# The screen was measured with ZERO misses at 0.01 (no screened cell with a
# Davies p below 0.001) over ~250k null cells of three trait/prevalence arms
# and ~4.5k cells with signal; the per-rho burden value it does NOT touch
# differs between Liu and Davies by 0.13x-4.6x, which is why the screen is
# outer-only.
#
# `SKATO_screen_threshold = 0` and `SKATO_tail_method = "davies"` are the
# upstream kernel, byte for byte (test-step2-skato-screen.R pins it). The
# default threshold is 0.05: in a multi-biobank meta-analysis of SKAT-O
# p-values a per-cohort bias adds coherently, and 0.05 confines Liu to the
# band where its median bias is exactly 1.00.
#
# Provenance is a SIDECAR, `<SAIGEOutputFile>.skatoMethod.txt`, in the shape
# of `.pooledTests.txt` (LEDGER #127): one row per tested cell naming which
# path produced `Pvalue` and the screen's Liu value. No result-file schema
# moves.

skato_screen_env <- new.env(parent = emptyenv())
skato_screen_env$threshold <- 0
skato_screen_env$tail <- "davies"
skato_screen_env$tailThreshold <- 1e-5
skato_screen_env$rows <- list()
skato_screen_env$counts <- c(screened = 0L, davies = 0L, single = 0L, saddlepoint = 0L)

#' Set the run-wide SKAT-O kernel settings. Called once per SPAGMMATtest
#' region run, BEFORE the SKAT-O pool forks, so every worker inherits them.
#' @noRd
skatoScreenSet <- function(threshold, tailMethod, tailThreshold = 1e-5) {
  if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) ||
        threshold < 0 || threshold > 1) {
    stop("SKATO_screen_threshold must be a single number in [0, 1]; 0 disables the screen")
  }
  if (!is.numeric(tailThreshold) || length(tailThreshold) != 1 || is.na(tailThreshold) ||
        tailThreshold < 0 || tailThreshold > 1) {
    stop("SKATO_tail_threshold must be a single number in [0, 1]")
  }
  tailMethod <- match.arg(tailMethod, c("saddlepoint", "davies"))
  skato_screen_env$threshold <- as.numeric(threshold)
  skato_screen_env$tail <- tailMethod
  skato_screen_env$tailThreshold <- as.numeric(tailThreshold)
  skatoScreenReset()
  invisible(NULL)
}

#' Forget the previous run's records (one dataset per process, but the tests
#' drive the region scan repeatedly).
#' @noRd
skatoScreenReset <- function() {
  skato_screen_env$rows <- list()
  skato_screen_env$counts <- c(screened = 0L, davies = 0L, single = 0L, saddlepoint = 0L)
  invisible(NULL)
}

skatoScreenOn <- function() {
  skato_screen_env$threshold > 0 || skato_screen_env$tail != "davies"
}

#' The SKAT-O kernel with the outer integral screened. Returns what
#' SKAT:::Met_SKAT_Get_Pvalue returns (p.value, param, p.value.resampling),
#' plus attributes `skato_method` ("liu-screen", "davies" or "single-rho")
#' and `skato_screen_p` (the Liu-form omnibus, NA when not computed).
#'
#' The preamble and the per-rho layer are transcriptions of SKAT 2.2.5's
#' Met_SKAT_Get_Pvalue -> SKAT_META_Optimal -> SKAT_META_Optimal_Get_Pvalue ->
#' SKAT_Optimal_Each_Q for method "optimal.adj" without resampling; the same
#' SKAT internals are called on the same inputs, so the per-rho values are
#' the ones SKAT would compute. Anything SKAT routes away from SKAT-O (an
#' all-zero Phi, a 1 x 1 Phi, a rank <= 1 Phi at q <= 10, a single rho) goes
#' to SKAT itself untouched.
#' @noRd
skatoScreenKernel <- function(Score, Phi, r.corr) {
  skatNs <- asNamespace("SKAT")
  original <- function() {
    get("Met_SKAT_Get_Pvalue", skatNs)(Score = Score, Phi = Phi, r.corr = r.corr,
                                       method = "optimal.adj", Score.Resampling = NULL)
  }
  p.m <- nrow(Phi)
  singleRho <- length(r.corr) <= 1 || sum(abs(Phi)) == 0 || length(Phi) <= 1 ||
    (ncol(Phi) <= 10 && qr(Phi)$rank <= 1)
  if (singleRho) {
    res <- original()
    attr(res, "skato_method") <- "single-rho"
    attr(res, "skato_screen_p") <- NA_real_
    return(res)
  }
  # ---- SKAT_META_Optimal: rho >= 0.999 -> 0.999, Q per rho, Phi / 2 -------
  r.all <- r.corr
  r.all[r.all >= 0.999] <- 0.999
  n.r <- length(r.all)
  Qr <- get("SKAT_META_Optimal_Get_Q", skatNs)(Score, r.all)$Q.r
  Phi2 <- Phi / 2
  # ---- SKAT_META_Optimal_Get_Pvalue: the per-rho spectra and param.m --------
  lambda.all <- vector("list", n.r)
  for (i in seq_len(n.r)) {
    rc <- r.all[i]
    RM <- diag(rep(1 - rc, p.m)) + matrix(rep(rc, p.m * p.m), ncol = p.m)
    L <- chol(RM, pivot = TRUE)
    lambda.all[[i]] <- get("Get_Lambda", skatNs)(L %*% (Phi2 %*% t(L)))
  }
  param.m <- get("SKAT_META_Optimal_Param", skatNs)(Phi2, r.all)
  # ---- SKAT_Optimal_Each_Q, method "optimal.adj": Davies per rho ------------
  pval <- numeric(n.r)
  param.mat <- matrix(0, n.r, 3)
  for (i in seq_len(n.r)) {
    lam <- lambda.all[[i]]
    c1 <- c(sum(lam), sum(lam^2), sum(lam^3), sum(lam^4))
    pt <- get("Get_Liu_Params_Mod", skatNs)(c1)
    pval[i] <- get("Get_PValue.Lambda", skatNs)(lam, Qr[i])$p.value
    param.mat[i, ] <- c(pt$muQ, pt$sigmaQ^2, pt$l)
  }
  pmin <- min(pval)
  pmin.q <- numeric(n.r)
  for (i in seq_len(n.r)) {
    df <- param.mat[i, 3]
    q.org <- qchisq(1 - pmin, df = df)
    pmin.q[i] <- (q.org - df) / sqrt(2 * df) * sqrt(param.mat[i, 2]) + param.mat[i, 1]
  }
  # ---- the screen: SKAT's Liu-form outer integral, with SKAT's post-rules ---
  pLiu <- get("SKAT_Optimal_PValue_Liu", skatNs)(pmin.q, param.m, r.all, pmin)
  multi <- if (n.r < 3) 2 else 3
  IDX <- which(pval > 0)
  if (pLiu <= 0 || length(IDX) < n.r) pLiu <- min(pval) * multi
  if (pLiu == 0 && length(IDX) > 0) pLiu <- min(pval[IDX])
  # threshold 0 = no screen: every cell takes the Davies path (and, below the
  # tail threshold, the tail)
  if (skato_screen_env$threshold == 0 || pLiu <= skato_screen_env$threshold) {
    res <- original()
    method <- "davies"
    if (skato_screen_env$tail == "saddlepoint" && !is.null(res$p.value) &&
          is.finite(res$p.value) && res$p.value < skato_screen_env$tailThreshold) {
      # the tail: both layers again with the saddlepoint, direct integral form
      tail <- skatoSaddlepointTail(Qr, lambda.all, param.m, r.all, pval, param.mat)
      if (!is.null(tail)) {
        res <- tail
        method <- "saddlepoint"
      }
    }
    attr(res, "skato_method") <- method
    attr(res, "skato_screen_p") <- pLiu
    return(res)
  }
  # ---- screened: SKAT_META_Optimal's result, with the Liu omnibus as p.value
  param <- list(p.val.each = pval, q.val.each = Qr)
  param$rho <- r.all
  param$minp <- pmin
  id_temp <- which(param$p.val.each == pmin)
  id_temp1 <- which(param$rho >= 0.999)
  if (length(id_temp1) > 0) param$rho[id_temp1] <- 1
  param$rho_est <- param$rho[id_temp]
  res <- list(p.value = pLiu, param = param, p.value.resampling = NULL)
  attr(res, "skato_method") <- "liu-screen"
  attr(res, "skato_screen_p") <- pLiu
  res
}

#' Record which path produced one cell's `Pvalue`. Called from
#' region_postprocess where the labels are in scope; the collect pass of the
#' parallel path hands a stub answer with no method and records nothing.
#' @noRd
skatoMethodRecord <- function(region, group, maxMAF, weight, test, method, screenP) {
  if (is.null(method)) {
    return(invisible(NULL))
  }
  k <- switch(method, "liu-screen" = "screened", "davies" = "davies",
              "saddlepoint" = "saddlepoint", "single")
  skato_screen_env$counts[k] <- skato_screen_env$counts[k] + 1L
  skato_screen_env$rows[[length(skato_screen_env$rows) + 1L]] <- data.frame(
    Region = as.character(region),
    Group = as.character(group),
    max_MAF = as.character(maxMAF),
    Weight = if (is.null(weight) || length(weight) == 0) NA_character_ else as.character(weight),
    Test = test,
    Method = method,
    Liu_p = if (is.null(screenP) || length(screenP) == 0) NA_real_ else as.numeric(screenP),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

#' Print the end-of-run line (always, so the mode a run took is on record --
#' a knob that failed to arrive would otherwise pass as green, LEDGER #30's
#' shape) and, when the screen is on, write the sidecar.
#' @noRd
skatoMethodSummary <- function(OutputFile) {
  if (!skatoScreenOn()) {
    cat("SKAT-O screen: off (SKATO_screen_threshold = 0, SKATO_tail_method = davies); ",
        "the SKAT kernel ran unchanged.\n", sep = "")
    return(invisible(NULL))
  }
  n <- skato_screen_env$counts
  path <- paste0(OutputFile, ".skatoMethod.txt")
  cat("SKAT-O screen: threshold ", skato_screen_env$threshold, ", tail ",
      skato_screen_env$tail, " below ", skato_screen_env$tailThreshold, "; ",
      n[["screened"]], " cell(s) took the Liu-form outer integral, ", n[["davies"]],
      " the Davies integral, ", n[["saddlepoint"]], " the saddlepoint tail, ",
      n[["single"]], " were single-rho. Per-cell provenance in '", path, "'.\n", sep = "")
  if (length(skato_screen_env$rows) > 0) {
    tab <- do.call(rbind, skato_screen_env$rows)
    utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# LEDGER #158, step 2: the tail. Below SKATO_tail_threshold, SKAT's own p-value
# is unreliable in BOTH directions: its per-rho Davies at acc 1e-6 (an
# ABSOLUTE tolerance) returns a value that is noise, or non-positive -- in
# which case modified Liu, 3-7x anticonservative there, is substituted -- and
# the outer integral inherits it. Measured against an adaptive-Davies
# reference on the binary trait, SKAT as shipped was 0.05x-7.5x off in
# (1e-8, 1e-6]. Kuonen's saddlepoint (Lugannani-Rice form) for a weighted sum
# of chi-squares has bounded relative error in the far tail: 5-12%
# conservative and stable in every band on both traits, 8-13% against Ruben's
# exact series on synthetic spectra to 1e-23. REGENIE uses the same switch at
# the same threshold (Davies, then Kuonen below 1e-5, after Wu et al. 2016).
#
# Both layers are recomputed with the saddlepoint, and the outer integral is
# taken in the direct form P = P(chi2_1 > x*) + int_0^{x*} S(delta_min(x)) f(x)
# dx (Wu's formula (3) shape) rather than SKAT's 1 - int, which underflows
# below ~1e-16. x = t^2 removes the density's singularity at 0. A rho or a
# node where the saddlepoint is undefined (the statistic at its mean) takes
# the Davies value already computed. SKAT's post-rules apply as on every path.

#' P(sum_j lambda_j chi2_1 > q) by Kuonen's saddlepoint; NA where undefined.
#' @noRd
skatoSaddlepointSurvival <- function(lambda, q) {
  if (q <= 0) {
    return(1)
  }
  mu <- sum(lambda)
  lmax <- max(lambda)
  if (abs(q - mu) < 1e-7 * mu) {
    return(NA_real_)
  }
  kPrime <- function(z) sum(lambda / (1 - 2 * z * lambda)) - q
  if (q > mu) {
    bracket <- c(0, (1 / (2 * lmax)) * (1 - 1e-12))
  } else {
    lo <- -1 / lmax
    while (kPrime(lo) > 0) lo <- lo * 2
    bracket <- c(lo, 0)
  }
  zh <- tryCatch(uniroot(kPrime, bracket, tol = 1e-12, maxiter = 5000)$root,
                 error = function(e) NA_real_)
  if (is.na(zh)) {
    return(NA_real_)
  }
  k0 <- -0.5 * sum(log(1 - 2 * zh * lambda))
  k2 <- sum(2 * lambda^2 / (1 - 2 * zh * lambda)^2)
  w <- sign(zh) * sqrt(2 * (zh * q - k0))
  v <- zh * sqrt(k2)
  if (!is.finite(w) || abs(w) < 1e-6 || !is.finite(v)) {
    return(NA_real_)
  }
  pnorm(w + log(v / w) / w, lower.tail = FALSE)
}

#' The saddlepoint tier over the layers the screen already built. `pvalDavies`
#' are SKAT's per-rho values (the fallback per rho); `param.mat` their Liu
#' parameters (muQ, varQ, df) per rho, which the pmin -> quantile map needs.
#' Returns a result in Met_SKAT_Get_Pvalue's shape, or NULL if the outer
#' integral could not be evaluated (the caller then keeps SKAT's answer).
#' @noRd
skatoSaddlepointTail <- function(Qr, lambda.all, param.m, r.all, pvalDavies, param.mat) {
  n.r <- length(r.all)
  pval <- numeric(n.r)
  for (i in seq_len(n.r)) {
    lam <- lambda.all[[i]]
    p <- if (length(lam) == 1) pvalDavies[i] else skatoSaddlepointSurvival(lam, Qr[i])
    if (is.na(p) || p <= 0 || p > 1) p <- pvalDavies[i]
    pval[i] <- p
  }
  pmin <- min(pval)
  pmin.q <- numeric(n.r)
  for (i in seq_len(n.r)) {
    df <- param.mat[i, 3]
    q.org <- qchisq(1 - pmin, df = df)
    pmin.q[i] <- (q.org - df) / sqrt(2 * df) * sqrt(param.mat[i, 2]) + param.mat[i, 1]
  }
  sd1 <- sqrt(param.m$VarQ - param.m$VarRemain) / sqrt(param.m$VarQ)
  daviesNode <- function(q) {
    o <- get("SKAT_davies", asNamespace("SKAT"))(q, param.m$lambda, acc = 1e-6)
    if (o$ifault != 0 || o$Qq <= 0 || o$Qq > 1) NA_real_ else o$Qq
  }
  sFun <- function(x) {
    temp <- (pmin.q - param.m$tau %x% t(x)) / (1 - r.all)
    tempMin <- apply(temp, 2, min)
    vapply(tempMin, function(min1) {
      if (min1 > sum(param.m$lambda) * 1e4) {
        return(0)
      }
      q <- (min1 - param.m$MuQ) * sd1 + param.m$MuQ
      s <- skatoSaddlepointSurvival(param.m$lambda, q)
      if (is.na(s)) s <- daviesNode(q)
      if (is.na(s)) stop("SKAT-O saddlepoint tail: no survival value at a node")
      min(s, 1)
    }, 0)
  }
  xstar <- max(min(pmin.q / param.m$tau), 0)
  # delta_min(x) is the lower envelope of the eight lines
  # (pmin.q_i - tau_i x) / (1 - r_i), so the integrand has a kink wherever the
  # minimising rho changes. Integrate piecewise between the pairwise crossings
  # inside (0, x*), where it is smooth; QUADPACK otherwise reports "extremely
  # bad integrand behaviour" on some cells. In t = sqrt(x) throughout.
  slope <- -param.m$tau / (1 - r.all)
  icpt <- pmin.q / (1 - r.all)
  knots <- c()
  for (i in seq_len(n.r - 1)) {
    for (j in (i + 1):n.r) {
      d <- slope[i] - slope[j]
      if (d != 0) {
        x <- (icpt[j] - icpt[i]) / d
        if (is.finite(x) && x > 0 && x < xstar) knots <- c(knots, x)
      }
    }
  }
  edges <- sqrt(unique(c(0, sort(knots), xstar)))
  integrand <- function(t) sFun(t^2) * 2 * dnorm(t)
  piecewise <- function(relTol) {
    total <- 0
    for (k in seq_len(length(edges) - 1)) {
      if (edges[k + 1] <= edges[k]) next
      re <- integrate(integrand, lower = edges[k], upper = edges[k + 1],
                      subdivisions = 1000, rel.tol = relTol, abs.tol = 0)
      total <- total + re$value
    }
    total
  }
  area <- tryCatch(piecewise(1e-6), error = function(e) NULL)
  if (is.null(area)) area <- tryCatch(piecewise(1e-4), error = function(e) NULL)
  if (is.null(area)) {
    return(NULL)
  }
  pvalue <- pchisq(xstar, df = 1, lower.tail = FALSE) + area
  if (pmin * n.r < pvalue) pvalue <- pmin * n.r
  multi <- if (n.r < 3) 2 else 3
  IDX <- which(pval > 0)
  if (pvalue <= 0 || length(IDX) < n.r) pvalue <- min(pval) * multi
  if (pvalue == 0 && length(IDX) > 0) pvalue <- min(pval[IDX])
  param <- list(p.val.each = pval, q.val.each = Qr)
  param$rho <- r.all
  param$minp <- pmin
  id_temp <- which(param$p.val.each == pmin)
  id_temp1 <- which(param$rho >= 0.999)
  if (length(id_temp1) > 0) param$rho[id_temp1] <- 1
  param$rho_est <- param$rho[id_temp]
  list(p.value = pvalue, param = param, p.value.resampling = NULL)
}
