# Gates that REFUSE a null fit which cannot be validly tested, instead of
# fitting it and leaving the problem to be discovered in step 2 -- or not
# discovered at all.
#
# Why they exist, with the measurement behind each:
#
# * Case count. This branch's 32-replicate null study
#   (`harness/work/nullrep-20260918/`, 100k cells per arm, N = 50,000) measured
#   the BINARY SKAT-O null tail at 2.7x its nominal rate at p = 1e-4 with 100
#   cases, and 1.5x at 1e-3 with 500; `Pvalue_SKAT` is 2x at 1e-3 at 100 cases.
#   The quantitative arm is calibrated. The cause is the kernel's SPA/ER Phi
#   adjustment at few cases -- present at tau2 = 0 and tau2 > 0 alike, so it is
#   not variance-component collapse, and it is identical for every p-value
#   option (LEDGER #159). So a small binary fit yields step-2 p-values that are
#   liberal exactly where they would be reported.
#
# * Covariate cell counts. A categorical covariate level with no cases is
#   COMPLETE SEPARATION: the logistic MLE for that level is infinite. Under
#   `isCovariateOffset = TRUE` (the default) that divergent coefficient is
#   frozen into `covoffset` and never re-estimated, so `mu` pins at 0 or 1,
#   `mu * (1 - mu)` collapses, the score variance goes to zero, and step-2
#   p-values underflow. A level with a handful of cases is the quasi-separated
#   version of the same failure. `minCovariateCellCount` defaults to 10, the
#   events-per-variable convention.
#
# These compose with the pre-existing `minCovariateCount` (default -1, off),
# which EXCLUDES a low-count binary covariate from the model rather than
# refusing to fit. The exclusion runs first; these gates then judge what
# survives, so setting both means "drop what you can, refuse what you cannot".

# Whether the eligibility gates are active at all. SAIGE_FIT_GATES=0 turns
# every gate in this file off at once, for step 1 and step 2 alike.
#
# It exists for one reason: the bundled fixture is deliberately BELOW both
# thresholds (98 cases against minCaseCount = 100; a10 level 2 has 7 cases
# against minCovariateCellCount = 10), so the test suite and the golden
# harness would otherwise have to carry an opt-out flag on ~150 command lines,
# and every command added later would have to remember it. harness/env.sh sets
# it, and run_tests.sh and every harness script already source that file. This
# is the same kill-switch shape the branch already uses for SAIGE_PGEN_SPARSE,
# SAIGE_PLINK_NO_MMAP and SAIGE_REGION_CHUNK_MEMORY_MB.
#
# It is NOT a quiet switch. A run with the gates disabled says so on stdout,
# because a safety gate that can be off without anyone noticing is not a
# safety gate. `announce` is TRUE at the one call per fit that should print.
fitGatesEnabled <- function(announce = FALSE) {
  raw <- Sys.getenv("SAIGE_FIT_GATES", unset = "1")
  enabled <- !(tolower(trimws(raw)) %in% c("0", "false", "off", "no"))
  if (!enabled && announce) {
    cat("WARNING: eligibility gates DISABLED by SAIGE_FIT_GATES=", raw,
        " -- case-count, covariate-cell, separation and convergence checks ",
        "are NOT being applied to this fit.\n", sep = "")
  }
  enabled
}

# Case count below which a binary fit is warned about but still run. Not an
# option: it is the level the null study above actually measured, and it is
# quoted in the warning so the user can see what it refers to.
saigeCaseCountAdvisory <- 500L

# One (covariate, level) tabulation against the phenotype, as offending rows.
cellsOfCovariate <- function(name, lev, y, threshold, isBinary) {
  out <- list()
  if (isBinary) {
    tab <- table(lev, as.character(y))
    for (l in rownames(tab)) {
      for (g in colnames(tab)) {
        if (tab[l, g] < threshold) {
          out[[length(out) + 1L]] <- data.frame(
            covariate = name, level = l,
            group = if (g %in% c("1", "TRUE")) "cases" else "controls",
            count = as.integer(tab[l, g]), stringsAsFactors = FALSE
          )
        }
      }
    }
  } else {
    tab <- table(lev)
    for (l in names(tab)) {
      if (tab[[l]] < threshold) {
        out[[length(out) + 1L]] <- data.frame(
          covariate = name, level = l, group = "samples",
          count = as.integer(tab[[l]]), stringsAsFactors = FALSE
        )
      }
    }
  }
  out
}

# The (level x phenotype) cells of every covariate for which a cell count means
# anything. Returns one row per cell strictly below `threshold`; a zero-row
# frame when there is nothing to report.
#
# Two things this has to get right, both measured:
#
# 1. By the time fitNULLGLMM calls it, preparePhenotypeFrame has expanded the
#    design with model.matrix, so a declared categorical is present as its
#    DUMMIES and its REFERENCE LEVEL HAS NO COLUMN AT ALL. A per-dummy scan
#    therefore cannot see a reference level with no cases: on a 3-level factor
#    whose reference level had zero cases, a dummy scan at threshold 10 fired
#    FALSE and a factor scan fired TRUE. So when the original column is gone,
#    the factor is RECONSTRUCTED -- a sample is in the reference level exactly
#    when every dummy of that covariate is 0. (The dummies are found by name
#    prefix, which is how model.matrix builds them; a covariate whose name is a
#    prefix of another's, e.g. `a1` beside `a10`, would over-match. Such a pair
#    is caught either way, only the label would be wrong.)
# 2. A covariate that is effectively binary is checked even when it was NOT
#    declared in qCovarCol -- that is the "someone forgot to list sex as
#    categorical" case, which is exactly the one that matters.
smallCovariateCells <- function(data, phenoCol, covarColList, qCovarCol,
                                threshold, isBinary) {
  empty <- data.frame(
    covariate = character(0), level = character(0),
    group = character(0), count = integer(0),
    stringsAsFactors = FALSE
  )
  if (!is.finite(threshold) || threshold <= 0) {
    return(empty)
  }
  y <- data[[phenoCol]]
  found <- list()
  claimed <- character(0)

  for (q in qCovarCol) {
    if (q %in% names(data)) {
      lev <- as.character(data[[q]])
    } else {
      # Search covarColList, NOT names(data): the frame also holds the
      # phenotype, IID and IndexGeno, and nothing outside the model's own
      # covariates may ever be counted here.
      dum <- grep(paste0("^", q), covarColList, value = TRUE)
      dum <- dum[dum %in% names(data)]
      dum <- dum[vapply(dum, function(z) all(data[[z]] %in% c(0, 1, NA)),
                        logical(1))]
      if (length(dum) == 0L) {
        next
      }
      claimed <- c(claimed, dum)
      m <- as.matrix(data[, dum, drop = FALSE])
      labels <- sub(paste0("^", q), "", dum)
      lev <- ifelse(rowSums(m) == 0, "<reference>",
                    labels[max.col(m, ties.method = "first")])
    }
    if (length(unique(lev[!is.na(lev)])) > 20L) {
      next
    }
    found <- c(found, cellsOfCovariate(q, lev, y, threshold, isBinary))
  }

  for (cc in setdiff(covarColList, c(qCovarCol, claimed))) {
    if (!cc %in% names(data)) {
      next
    }
    v <- data[[cc]]
    if (length(unique(v[!is.na(v)])) != 2L) {
      next
    }
    found <- c(found, cellsOfCovariate(cc, as.character(v), y, threshold,
                                       isBinary))
  }

  if (length(found) == 0L) {
    return(empty)
  }
  do.call(rbind, found)
}

# One printable line per offending cell.
formatSmallCells <- function(bad) {
  paste0(
    "  ", bad$covariate, " = ", bad$level, ": ", bad$count, " ", bad$group,
    collapse = "\n"
  )
}

# Refuses a fit that cannot be validly tested. `data` must be the FINAL fitted
# sample set (dataMerge_sort -- already intersected with the genotypes, the
# sparse GRM and SampleIDIncludeFile), so every count here is a count the fit
# will actually see.
checkFitEligibility <- function(data, phenoCol, traitType, covarColList, qCovarCol,
                                minCaseCount, minSampleCount, minCovariateCellCount) {
  if (!fitGatesEnabled(announce = TRUE)) {
    return(invisible(TRUE))
  }
  isBinary <- identical(traitType, "binary")
  y <- data[[phenoCol]]
  nObs <- sum(!is.na(y))
  if (is.finite(minSampleCount) && minSampleCount > 0 && nObs < minSampleCount) {
    stop(
      "ERROR! ", phenoCol, " has ", nObs, " samples in the fitted set, below ",
      "--minSampleCount=", minSampleCount, ". Not fitting. Lower or disable the ",
      "gate with --minSampleCount=0 to fit anyway.\n"
    )
  }
  if (isBinary) {
    nCase <- sum(y == 1, na.rm = TRUE)
    cat("cases in the fitted set: ", nCase, " of ", nObs, "\n", sep = "")
    if (is.finite(minCaseCount) && minCaseCount > 0 && nCase < minCaseCount) {
      stop(
        "ERROR! ", phenoCol, " has ", nCase, " cases in the fitted set, below ",
        "--minCaseCount=", minCaseCount, ". A binary null this small gives ",
        "step-2 p-values that are liberal in the tail (measured: 2.7x the ",
        "nominal rate at 1e-4 with 100 cases). Not fitting. Use ",
        "--minCaseCount=0 to fit anyway.\n"
      )
    }
    if (nCase < saigeCaseCountAdvisory) {
      warning(
        phenoCol, " has ", nCase, " cases. This branch's null study measured the ",
        "SKAT-O tail at 1.5x nominal at 1e-3 with ", saigeCaseCountAdvisory,
        " cases and 2.7x at 1e-4 with 100, so treat step-2 tail p-values for ",
        "this phenotype with caution.",
        call. = FALSE
      )
    }
  }
  bad <- smallCovariateCells(
    data, phenoCol, covarColList, qCovarCol,
    minCovariateCellCount, isBinary
  )
  if (nrow(bad) > 0) {
    zero <- any(bad$count == 0L)
    stop(
      "ERROR! a categorical covariate has a level with fewer than ",
      "--minCovariateCellCount=", minCovariateCellCount,
      if (isBinary) " in cases or controls" else " samples", ":\n",
      formatSmallCells(bad), "\n",
      if (zero) {
        paste0(
          "A level with 0 is COMPLETE SEPARATION: the coefficient is infinite, ",
          "and with --isCovariateOffset=TRUE it is frozen into the offset, ",
          "pinning fitted values at 0 or 1 and collapsing the score variance.\n"
        )
      } else {
        ""
      },
      "Drop the covariate, collapse its levels, restrict the analysis (see ",
      "--MaleOnly/--FemaleOnly for a single-sex phenotype), or set ",
      "--minCovariateCellCount=0 to fit anyway.\n"
    )
  }
  invisible(TRUE)
}

# R's own separation diagnostic, on a GLM that is being fitted anyway
# (`modwitcov` under isCovariateOffset = TRUE, `fit0` otherwise): the IRLS
# failing to converge, or a fitted probability at the 0/1 boundary. Catches
# separation on a CONTINUOUS covariate -- an interaction such as age_sex --
# which no cell count can see. `action` is "error", "warn" or "none".
checkSeparation <- function(fit, phenoCol, action, maxFraction = 0.02) {
  if (!fitGatesEnabled()) {
    return(invisible(TRUE))
  }
  # A value that is none of the three is an error, not a silent downgrade: an
  # unrecognised mode that quietly disables a check is how a gate stops gating
  # without anyone noticing (cf. SAIGE_COMPARE_MODE=bytes).
  if (!identical(length(action), 1L) || !action %in% c("error", "warn", "none")) {
    stop(
      "ERROR! --separationCheck must be one of \"error\", \"warn\" or ",
      "\"none\"; got ", paste(deparse(action), collapse = " "), "\n"
    )
  }
  if (identical(action, "none")) {
    return(invisible(TRUE))
  }
  eps <- 1e-8
  mu <- fit$fitted.values
  frac <- sum(mu < eps | mu > 1 - eps, na.rm = TRUE) / length(mu)
  # Same fraction rule as the step-2 check: one pinned sample in thousands is
  # a rounding error, not a collapse. Non-convergence of the IRLS is not a
  # fraction and is always reported.
  if (isTRUE(fit$converged) && frac <= maxFraction) {
    return(invisible(TRUE))
  }
  cf <- coef(fit)
  cf <- cf[!is.na(cf) & names(cf) != "(Intercept)"]
  worst <- if (length(cf) > 0) {
    ord <- order(abs(cf), decreasing = TRUE)[seq_len(min(3L, length(cf)))]
    paste0("  ", names(cf)[ord], ": ", signif(cf[ord], 4), collapse = "\n")
  } else {
    "  (no covariates)"
  }
  msg <- paste0(
    "the covariate GLM for ", phenoCol, " shows separation (",
    if (!isTRUE(fit$converged)) "did not converge" else "converged",
    "; ", signif(100 * frac, 3), "% of fitted probabilities at the 0/1 ",
    "boundary, range ", signif(min(mu, na.rm = TRUE), 3), " .. ",
    signif(max(mu, na.rm = TRUE), 3), "). Largest coefficients:\n", worst, "\n",
    "With --isCovariateOffset=TRUE these are frozen into the offset and never ",
    "re-estimated, which pins fitted values at 0 or 1 and collapses the step-2 ",
    "score variance. Drop or collapse the covariate, or set ",
    "--separationCheck=warn or none to fit anyway."
  )
  if (identical(action, "warn")) {
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  stop("ERROR! ", msg, "\n")
}

# mu at the 0/1 boundary, on a model that has already been fitted. Shared by
# the step-2 gate; `label` names the model in the message.
checkMuSeparation <- function(mu, label, action, maxFraction = 0.02) {
  if (!identical(length(action), 1L) || !action %in% c("error", "warn", "none")) {
    stop(
      "ERROR! --separationCheck must be one of \"error\", \"warn\" or ",
      "\"none\"; got ", paste(deparse(action), collapse = " "), "\n"
    )
  }
  if (identical(action, "none")) {
    return(invisible(TRUE))
  }
  eps <- 1e-8
  atBound <- mu < eps | mu > 1 - eps
  nBound <- sum(atBound, na.rm = TRUE)
  if (nBound == 0L) {
    return(invisible(TRUE))
  }
  frac <- nBound / length(mu)
  msg <- paste0(
    "the null model in ", label, " has ", nBound, " of ", length(mu),
    " fitted values (", signif(100 * frac, 3), "%) at the 0/1 boundary ",
    "(range ", signif(min(mu, na.rm = TRUE), 3), " .. ",
    signif(max(mu, na.rm = TRUE), 3), "). That is what separation in the ",
    "covariate model does: a pinned sample contributes ~0 to the score ",
    "variance instead of mu*(1 - mu), so it carries no information and the ",
    "variance these tests divide by shrinks. Refit without the separating ",
    "covariate (see --MaleOnly/--FemaleOnly for a single-sex phenotype), ",
    "raise --separationMaxFraction, or set --separationCheck=warn or none."
  )
  # Below the fraction the loss is real but small -- a pinned sample costs one
  # sample's worth of information, so 0.5% pinned is a 0.5% loss -- and saying
  # nothing would hide it, so it warns rather than passing silently. Refusing
  # on ANY pinned sample (the rule this replaces) refused a model over fewer
  # than 20 samples in 2,737, which is not variance collapse by any reading.
  if (frac <= maxFraction) {
    warning(msg, " (below --separationMaxFraction=", maxFraction,
            ", so this is a warning and the run continues)", call. = FALSE)
    return(invisible(FALSE))
  }
  if (identical(action, "warn")) {
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  stop("ERROR! ", msg, "\n")
}

# The step-2 half of the gates, applied to a null model as LOADED.
#
# This is NOT redundant with the step-1 gate. A model reaching step 2 may have
# been fitted by a SAIGE that had no gates at all -- everything fitted before
# they existed, and everything from upstream -- so for such a model this is the
# only check it ever meets. It also catches a .rda produced elsewhere.
#
# Everything it needs is already in the file, with nothing inferred:
#   modglmm$y              the phenotype, so the case count is exact;
#   modglmm$fitted.values  the POST-FIT mu, which is the quantity the tests
#                          divide by. mu pinned at the 0/1 boundary IS the
#                          damage separation does, which makes this strictly
#                          more complete than step 1's covariate-cell scan --
#                          that scan sees only dummy columns, so it cannot see
#                          a zero-case reference level, while this can.
checkModelEligibility <- function(modglmm, modelFile, minCaseCount,
                                  separationCheck, convergenceCheck = "error",
                                  separationMaxFraction = 0.02) {
  if (!fitGatesEnabled(announce = TRUE)) {
    return(invisible(TRUE))
  }
  label <- basename(modelFile)
  # Convergence is not a binary-trait question, so it is checked first.
  checkConvergence(modglmm$converged, label, convergenceCheck)
  if (!identical(modglmm$traitType, "binary")) {
    return(invisible(TRUE))
  }
  y <- modglmm$y
  nCase <- sum(y == 1, na.rm = TRUE)
  cat("null model ", label, ": ", nCase, " cases of ", length(y), "\n", sep = "")
  if (is.finite(minCaseCount) && minCaseCount > 0 && nCase < minCaseCount) {
    stop(
      "ERROR! the null model in ", label, " was fitted on ", nCase,
      " cases, below --minCaseCount=", minCaseCount, ". Region and marker ",
      "p-values from a binary null this small are liberal in the tail ",
      "(measured: 2.7x the nominal rate at 1e-4 with 100 cases). Not testing. ",
      "Use --minCaseCount=0 to test anyway.\n"
    )
  }
  if (nCase < saigeCaseCountAdvisory) {
    warning(
      "the null model in ", label, " was fitted on ", nCase, " cases. This ",
      "branch's null study measured the SKAT-O tail at 1.5x nominal at 1e-3 ",
      "with ", saigeCaseCountAdvisory, " cases and 2.7x at 1e-4 with 100, so ",
      "treat tail p-values from it with caution.",
      call. = FALSE
    )
  }
  checkMuSeparation(as.vector(modglmm$fitted.values), label, separationCheck,
                    separationMaxFraction)
}

# Did the null fit actually converge? `modglmm$converged` is set TRUE only at
# the two genuine exits of the AI-REML loop -- the tau2 = 0 optimum, and the
# criterion being met (LEDGER #123 made the flag record the exit rather than
# the iteration count). FALSE therefore covers TWO different things, and the
# flag alone cannot tell them apart:
#   (a) the loop ran out, completing maxiter iterations without meeting its
#       criterion;
#   (b) the loop BAILED OUT, because max(tau) exceeded tol^(-2), and printed
#       "Large variance estimate observed in the iterations".
# The step-1 log distinguishes them. Note that tol does double duty there, so
# a loose tol makes convergence easier while making the bail-out far more
# aggressive: tol = 0.02 aborts above tau 2500, tol = 0.2 above tau 25.
#
# Nothing in step 2 read this field before 2026-09-23, so an unconverged model
# was tested exactly like a converged one. A model whose variance components
# never settled has a variance ratio, and therefore test statistics, that mean
# nothing in particular. `action` is "error", "warn" or "none".
#
# A model file that predates the field entirely reports NA rather than failing:
# absence of the flag is not evidence of non-convergence.
checkConvergence <- function(converged, label, action) {
  if (!identical(length(action), 1L) || !action %in% c("error", "warn", "none")) {
    stop(
      "ERROR! --convergenceCheck must be one of \"error\", \"warn\" or ",
      "\"none\"; got ", paste(deparse(action), collapse = " "), "\n"
    )
  }
  if (identical(action, "none")) {
    return(invisible(TRUE))
  }
  if (is.null(converged) || length(converged) != 1L || is.na(converged)) {
    cat("note: ", label, " carries no usable convergence flag; not checked\n",
        sep = "")
    return(invisible(TRUE))
  }
  if (isTRUE(converged)) {
    return(invisible(TRUE))
  }
  msg <- paste0(
    "the null model ", label, " did NOT converge, so the variance components ",
    "-- and the variance ratio derived from them -- did not settle, and test ",
    "statistics scaled by them are not interpretable. The AI-REML loop either ",
    "ran out at --maxiter without meeting its criterion, or bailed out with ",
    "\"Large variance estimate observed in the iterations\" because max(tau) ",
    "exceeded tol^(-2); the step-1 log says which. Note that --tol sets both ",
    "the criterion and that bail-out, so a loose tol converges more easily ",
    "while aborting far sooner (0.02 -> tau 2500, 0.2 -> tau 25). Refit, or ",
    "set --convergenceCheck=warn or none to proceed anyway."
  )
  if (identical(action, "warn")) {
    warning(msg, call. = FALSE)
    return(invisible(FALSE))
  }
  stop("ERROR! ", msg, "\n")
}
