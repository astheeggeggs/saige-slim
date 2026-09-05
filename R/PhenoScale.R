# LEDGER #65 FIXED (2026-09-05, owner decision: standardise).
#
# The quantitative null model has several rules that are absolute in phenotype
# units: in the AI-REML loop a floor of 0.1 on var(Y) that was a stop(), the
# `+ tol` in the convergence rule's denominator and the divergence cap
# max(tau) > tol^-2 = 2500, plus a fixed tau start of (1, 0); and on the
# variance-ratio path two 1e-4 floors on the diagonal of Sigma = tau0/W +
# tau1*GRM in src/SAIGE_fitGLMM_fast.cpp (gen_sp_Sigma, getDiagOfSigma) and the
# absolute PCG stop rule (LEDGER #6). So a trait in millimetres and the same
# trait in metres took different paths: on the bundled fixture y * 100 hit the
# cap and y / 100 hit the stop(), and with only the fit standardised, y / 100
# still produced a variance ratio 33% off because tau ~ 2e-5 sat under the
# 1e-4 floor.
#
# fitNULLGLMM therefore divides a quantitative phenotype by the sd of its
# covariate-adjusted values right after the inverse normalisation, runs the
# glm, the mixed-model fit AND the variance-ratio estimation on that
# standardised trait, and scales the saved model(s) and the variance-ratio
# file back to phenotype units at the end. Every quantity a user or step 2
# reads is in phenotype units; the arithmetic in between is scale-free.
#
# Scaling rules, for a factor f applied to y: variances (theta, cov) by f^2;
# everything linear in y (coefficients, linear predictors, fitted values, the
# working vector Y, residuals, y, offsets, alpha0) by f; obj.noK is recomputed
# from the rescaled quantities (its pieces mix 1/theta[1] and residuals). The
# null and null_noXadj variance ratios, whose numerator variance is in 1/y^2
# units and whose denominator g'g is in genotype units only, scale by 1/f^2 --
# measured 2026-09-05 on the fixture: y * 100 gives ratios 1e-4 x the unscaled
# ones to 5e-8 relative, and step-2 p-values identical to 1.5e-7. The sparse
# ratio is scale-free (both sides go through Sigma^-1) and is left alone
# (LEDGER #97: scaling it too was an 11-12% error on the N = 50,000 fixture).

# sd of the covariate-adjusted phenotype: the residual of the plain linear model
# `formula` on `data`. Stops on a constant phenotype, the one case nothing can
# fit; everything else is a valid scale.
quantPhenoScale <- function(formula, data, phenoCol) {
  res <- stats::residuals(stats::lm(formula, data = data))
  s <- stats::sd(res)
  if (!is.finite(s) || s <= 0) {
    stop(
      "The quantitative phenotype ", phenoCol, " has zero variance after ",
      "covariate adjustment (sd of the linear-model residual is ", s,
      "); nothing can be fitted.\n"
    )
  }
  s
}

# Multiply every phenotype-unit quantity of a fitted quantitative null model by
# `factor` (f = yScale to return to phenotype units, 1 / yScale to standardise
# a model read from disk). Fields absent from the model are left absent.
scaleQuantModel <- function(modglmm, factor) {
  if (!identical(modglmm$traitType, "quantitative")) {
    stop("scaleQuantModel: the model is not a quantitative null model\n")
  }
  modglmm <- scaleQuantFit(modglmm, factor, X = modglmm$X)
  if (isTRUE(modglmm$LOCO) && !is.null(modglmm$LOCOResult)) {
    for (j in seq_along(modglmm$LOCOResult)) {
      r <- modglmm$LOCOResult[[j]]
      if (is.list(r) && isTRUE(r$isLOCO)) {
        # a LOCO entry has no y or theta of its own: its obj.noK is rebuilt
        # from the model's, which the call above has already put in the
        # target units (LEDGER #99: these used to be copied INTO the entry and
        # then scaled a second time with the entry's own fields)
        modglmm$LOCOResult[[j]] <- scaleQuantFit(r, factor, X = modglmm$X,
                                                 y = modglmm$y, theta = modglmm$theta)
      }
    }
  }
  modglmm
}

# One fit's fields (the top-level model or one LOCOResult entry): the variance
# terms by factor^2, everything linear in y by factor, and obj.noK rebuilt from
# the rescaled pieces. `X` is the design the score-test pieces were built on
# (modglmm$X: X or Xorig, whichever the fitter used). `y` and `theta` are what
# obj.noK is rebuilt from, ALREADY in the target units; left NULL they are the
# fit's own after scaling (the top-level model), and a LOCO entry, which
# carries neither, is handed the model's.
scaleQuantFit <- function(m, factor, X, y = NULL, theta = NULL) {
  if (!is.null(m$theta)) m$theta <- m$theta * factor^2
  if (!is.null(m$cov)) m$cov <- m$cov * factor^2
  for (k in c("coefficients", "linear.predictors", "fitted.values", "Y",
              "residuals", "y", "offset", "alpha0")) {
    if (!is.null(m[[k]])) m[[k]] <- m[[k]] * factor
  }
  if (is.null(y)) y <- m$y
  if (is.null(theta)) theta <- m$theta
  if (!is.null(m$obj.noK) && !is.null(m$fitted.values) && !is.null(y) &&
        !is.null(theta) && !is.null(X)) {
    mu2 <- rep(1 / theta[1], length(m$fitted.values))
    m$obj.noK <- ScoreTest_NULL_Model(m$fitted.values, mu2, y, X)
  }
  m
}

# Rescale the saved quantitative model files of a step-1 run in place: exactly
# the files it wrote (the main model file and, under --isLowMemLOCO, the
# _noLOCO file and each per-chromosome file). A file the run says it wrote and
# cannot find is an error, not a skip (LEDGER #98).
rescaleSavedQuantModels <- function(files, factor) {
  for (f in files) {
    if (!file.exists(f)) {
      stop("rescaleSavedQuantModels: ", f, " was written by this run but is missing\n")
    }
    e <- new.env()
    load(f, envir = e)
    modglmm <- scaleQuantModel(e$modglmm, factor)
    save(modglmm, file = f)
  }
  invisible(files)
}
