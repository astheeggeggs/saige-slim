# The step-1 eligibility gates (R/FitGates.R), added 2026-09-23.
#
# What they are for, and why the numbers are what they are, is on the
# functions. What matters here is that each gate is drilled BOTH ways: a
# configuration it must refuse and one it must pass. A gate that only ever
# sees the refusing case cannot tell you it is not refusing everything.
#
# The bundled fixture is BELOW both shipped defaults -- 98 cases against
# minCaseCount = 100, and a10 level 2 has 7 cases against
# minCovariateCellCount = 10 -- which is why every other step-1 test and every
# golden command turns them off. That is not a workaround: this fixture is a
# mechanical fixture, not a valid analysis, and the gates are saying so.
# harness/env.sh sets SAIGE_FIT_GATES=0 for the whole suite, so EVERY test
# here turns the gates back on for itself -- in-process with
# withr::local_envvar, and for a CLI run by passing the variable through.
# Tier: fast for the pure functions, standard for the three CLI runs.

fg_pheno <- function() {
  utils::read.table(
    file.path(saige_extdata, "input",
              "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"),
    header = TRUE, stringsAsFactors = FALSE
  )
}

# ---- smallCovariateCells ---------------------------------------------------

test_that("smallCovariateCells finds the fixture's thin cell and nothing else", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  # Established by direct tabulation: a9's smallest (level x phenotype) cell is
  # 14 and a10's is 7, so a threshold of 10 must find a10 and not a9.
  bad <- SAIGE:::smallCovariateCells(d, "y_binary", c("x1", "x2", "a9", "a10"),
                                     c("a9", "a10"), 10, isBinary = TRUE)
  expect_identical(sort(unique(bad$covariate)), "a10")
  expect_true(all(bad$count < 10))
  expect_true(all(bad$group %in% c("cases", "controls")))

  # A threshold below every cell finds nothing: the positive control's twin.
  expect_identical(
    nrow(SAIGE:::smallCovariateCells(d, "y_binary", c("x1", "x2", "a9", "a10"),
                                     c("a9", "a10"), 5, isBinary = TRUE)),
    0L
  )
  # 0 disables the gate outright.
  expect_identical(
    nrow(SAIGE:::smallCovariateCells(d, "y_binary", c("x1", "x2", "a9", "a10"),
                                     c("a9", "a10"), 0, isBinary = TRUE)),
    0L
  )
})

test_that("a continuous covariate has no cells to count", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  # x1 is continuous (1,000 distinct values) and must never be reported,
  # however small the threshold -- otherwise the gate would refuse every fit.
  cont <- SAIGE:::smallCovariateCells(d, "y_binary", "x1", character(0), 50,
                                      isBinary = TRUE)
  expect_identical(nrow(cont), 0L)
  # x2 IS binary in this fixture (461/29 | 441/69) and is NOT declared
  # categorical, yet must still be checked: that is the "someone forgot to put
  # sex in qCovarColList" case, which is the one that matters.
  bad <- SAIGE:::smallCovariateCells(d, "y_binary", "x2", character(0),
                                     50, isBinary = TRUE)
  expect_identical(unique(bad$covariate), "x2")
})

test_that("the quantitative arm counts level sizes, not cases", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  # a10 level 2 has 7 CASES but 50 SAMPLES, so it trips the binary gate at 10
  # and not the quantitative one -- the two arms measure different things.
  asBinary <- SAIGE:::smallCovariateCells(d, "y_binary", "a10", "a10", 10,
                                          isBinary = TRUE)
  asQuant <- SAIGE:::smallCovariateCells(d, "y_quantitative", "a10", "a10", 10,
                                         isBinary = FALSE)
  expect_true(nrow(asBinary) > 0L)
  expect_identical(nrow(asQuant), 0L)
})

# ---- checkFitEligibility ---------------------------------------------------

test_that("the case-count gate refuses below its threshold and passes above", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  args <- list(data = d, phenoCol = "y_binary", traitType = "binary",
               covarColList = c("x1", "x2"), qCovarCol = character(0),
               minSampleCount = 0, minCovariateCellCount = 0)
  expect_error(
    do.call(SAIGE:::checkFitEligibility, c(args, list(minCaseCount = 100))),
    "98 cases in the fitted set, below --minCaseCount=100"
  )
  # 98 cases is above 50 and the gate must then be silent about the count.
  expect_warning(
    do.call(SAIGE:::checkFitEligibility, c(args, list(minCaseCount = 50))),
    "98 cases"          # the sub-500 advisory, which is a warning not a stop
  )
  # 0 disables it entirely.
  expect_warning(
    do.call(SAIGE:::checkFitEligibility, c(args, list(minCaseCount = 0))),
    "98 cases"
  )
})

test_that("the sample-count gate is off by default and fires when asked", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  args <- list(data = d, phenoCol = "y_quantitative", traitType = "quantitative",
               covarColList = c("x1", "x2"), qCovarCol = character(0),
               minCaseCount = 0, minCovariateCellCount = 0)
  expect_silent(do.call(SAIGE:::checkFitEligibility,
                        c(args, list(minSampleCount = 0))))
  expect_error(
    do.call(SAIGE:::checkFitEligibility, c(args, list(minSampleCount = 5000))),
    "1000 samples in the fitted set, below --minSampleCount=5000"
  )
})

test_that("a zero-case covariate level is named as complete separation", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  # Plant the IFHern shape: a covariate level with cases and one with none.
  d$grp <- ifelse(seq_len(nrow(d)) %% 2 == 0, 1, 0)
  d$y_binary[d$grp == 1] <- 0                # level 1 now has zero cases
  err <- tryCatch(
    SAIGE:::checkFitEligibility(d, "y_binary", "binary", "grp", "grp",
                                minCaseCount = 0, minSampleCount = 0,
                                minCovariateCellCount = 10),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("COMPLETE SEPARATION", err, fixed = TRUE))
  expect_true(grepl("grp = 1: 0 cases", err, fixed = TRUE))
  expect_true(grepl("MaleOnly/--FemaleOnly", err, fixed = TRUE))
})

# ---- checkSeparation -------------------------------------------------------

test_that("checkSeparation is silent on a clean fit and fires on a separated one", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  clean <- stats::glm(y_binary ~ x1 + x2, data = d, family = stats::binomial)
  expect_silent(SAIGE:::checkSeparation(clean, "y_binary", "error"))

  d$sep <- as.numeric(d$y_binary == 1)        # perfectly predicts the outcome
  sepfit <- suppressWarnings(stats::glm(y_binary ~ x1 + sep, data = d,
                                        family = stats::binomial))
  expect_error(SAIGE:::checkSeparation(sepfit, "y_binary", "error"),
               "shows separation")
  expect_warning(SAIGE:::checkSeparation(sepfit, "y_binary", "warn"),
                 "shows separation")
  expect_silent(SAIGE:::checkSeparation(sepfit, "y_binary", "none"))
})

test_that("an unrecognised separationCheck is an error, not a silent no-op", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  d <- fg_pheno()
  clean <- stats::glm(y_binary ~ x1, data = d, family = stats::binomial)
  # The failure this guards against: a typo silently disabling the check.
  expect_error(SAIGE:::checkSeparation(clean, "y_binary", "None"),
               "must be one of")
  expect_error(SAIGE:::checkSeparation(clean, "y_binary", "off"),
               "must be one of")
})

# ---- the gates as the CLI actually applies them ----------------------------

fg_cli <- function(extra, phenoFile = NULL, phenoCol = "y_binary",
                   covar = "x1,x2,a9,a10", qcovar = "a9,a10") {
  prefix <- tempfile("s1fg_")
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  if (is.null(phenoFile)) {
    phenoFile <- "./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"
  }
  args <- c("--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
            paste0("--sparseGRMFile=", sg),
            paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            "--useSparseGRMtoFitNULL=TRUE", paste0("--phenoFile=", phenoFile),
            paste0("--phenoCol=", phenoCol), paste0("--covarColList=", covar),
            paste0("--qCovarColList=", qcovar), "--sampleIDColinphenoFile=IID",
            "--traitType=binary", "--IsOverwriteVarianceRatioFile=TRUE",
            paste0("--outputPrefix=", prefix), extra)
  st <- run_saige_cli("step1_fitNULLGLMM.R", args,
                      env = c(SAIGE_FIT_GATES = "1"))
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE),
       prefix = prefix)
}
# LEDGER #171 FIXED: a refused fit leaves NO file behind -- not the 0-byte .rda
# and .varianceRatio.txt it used to -- and says what it removed.
expect_no_leftovers <- function(r) {
  expect_false(file.exists(paste0(r$prefix, ".rda")))
  expect_false(file.exists(paste0(r$prefix, ".varianceRatio.txt")))
  expect_true(any(grepl("empty output file(s) left by a fit that did not complete", r$log,
                        fixed = TRUE)))
}

test_that("the shipped defaults refuse the bundled fixture, before fitting", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  skip_unless_level("standard")
  t0 <- proc.time()[["elapsed"]]
  r <- fg_cli(character(0))
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("below --minCaseCount=100", r$log, fixed = TRUE)))
  # The refusal must come BEFORE the fit, or it saves nothing: a full fit on
  # this fixture takes seconds, the refusal should be startup only.
  expect_false(any(grepl("Start fitting the NULL GLMM", r$log, fixed = TRUE)))
  expect_lt(proc.time()[["elapsed"]] - t0, 60)
  expect_no_leftovers(r)
})

test_that("the covariate-cell gate refuses on its own, with the case gate off", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  skip_unless_level("standard")
  r <- fg_cli("--minCaseCount=0")
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("--minCovariateCellCount=10", r$log, fixed = TRUE)))
  expect_true(any(grepl("a10 = 2", r$log, fixed = TRUE)))
  expect_no_leftovers(r)
})

test_that("with both gates off the same fit runs to completion", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  skip_unless_level("standard")
  # The positive control. Without it the two refusals above are consistent
  # with a gate that refuses everything.
  r <- fg_cli(c("--minCaseCount=0", "--minCovariateCellCount=0"))
  expect_identical(r$status, 0L)
  expect_true(any(grepl("Start fitting the NULL GLMM", r$log, fixed = TRUE)))
  # LEDGER #171: and a fit that completes keeps everything it wrote
  expect_gt(file.size(paste0(r$prefix, ".rda")), 0)
  expect_gt(file.size(paste0(r$prefix, ".varianceRatio.txt")), 0)
  expect_false(any(grepl("empty output file(s)", r$log, fixed = TRUE)))
})

# ---- the kill switch ------------------------------------------------------

test_that("SAIGE_FIT_GATES=0 disables the gates and says so out loud", {
  d <- fg_pheno()
  args <- list(data = d, phenoCol = "y_binary", traitType = "binary",
               covarColList = c("x1", "x2"), qCovarCol = character(0),
               minCaseCount = 1e6, minSampleCount = 0,
               minCovariateCellCount = 0)
  # ON: thresholds nothing can meet, so it must refuse.
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  expect_error(do.call(SAIGE:::checkFitEligibility, args), "below --minCaseCount")

  # OFF: the same impossible thresholds pass, AND the run announces it. A
  # switch that can turn a safety gate off silently is not acceptable, so the
  # announcement is asserted, not assumed.
  withr::local_envvar(c(SAIGE_FIT_GATES = "0"))
  expect_output(do.call(SAIGE:::checkFitEligibility, args),
                "eligibility gates DISABLED by SAIGE_FIT_GATES=0")
})

test_that("only the documented off-values disable the gates", {
  d <- fg_pheno()
  args <- list(data = d, phenoCol = "y_binary", traitType = "binary",
               covarColList = c("x1", "x2"), qCovarCol = character(0),
               minCaseCount = 1e6, minSampleCount = 0,
               minCovariateCellCount = 0)
  for (off in c("0", "false", "OFF", "No")) {
    withr::local_envvar(c(SAIGE_FIT_GATES = off))
    expect_output(do.call(SAIGE:::checkFitEligibility, args), "DISABLED")
  }
  # anything else leaves the gates ON -- an unrecognised value must not be
  # read as "off", which is the failure mode this asserts against.
  for (on in c("1", "yes", "", "true")) {
    withr::local_envvar(c(SAIGE_FIT_GATES = on))
    expect_error(do.call(SAIGE:::checkFitEligibility, args), "below --minCaseCount")
  }
})

test_that("a zero-case REFERENCE level survives model.matrix expansion", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  # preparePhenotypeFrame replaces a declared categorical with its DUMMIES, and
  # the reference level gets no column at all -- so a per-dummy scan is blind
  # to a reference level with no cases (measured: dummy scan FALSE, factor scan
  # TRUE). Build the frame exactly as the fitter sees it and assert the
  # reconstruction finds it.
  set.seed(1)
  n <- 900
  lv <- rep(0:2, each = n / 3)
  y <- c(rep(0, n / 3), rbinom(n / 3, 1, 0.3), rbinom(n / 3, 1, 0.3))
  x <- model.matrix(~f, data.frame(f = factor(lv)))[, -1, drop = FALSE]
  expanded <- data.frame(y = y, x)
  bad <- SAIGE:::smallCovariateCells(expanded, "y", colnames(x), "f",
                                     10, isBinary = TRUE)
  expect_true(any(bad$level == "<reference>" & bad$count == 0L))
  expect_identical(unique(bad$covariate), "f")
})

# ---- the step-2 half --------------------------------------------------------
# Not redundant with step 1: a model reaching step 2 may have been fitted by a
# SAIGE that had no gates at all, and for such a model this is the only check
# it ever meets. The committed golden is exactly that case -- 98 cases.

test_that("step 2 refuses a null model fitted on too few cases", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  f <- file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")
  expect_error(SAIGE:::ReadModel(f, chrom = "", LOCO = FALSE),
               "was fitted on 98 cases, below --minCaseCount=100")
  # The positive control: with the gate lowered the SAME file loads, so the
  # refusal above is about the count and not about the file.
  expect_warning(m <- SAIGE:::ReadModel(f, chrom = "", LOCO = FALSE,
                                        minCaseCount = 0))
  expect_identical(m$traitType, "binary")
  expect_identical(sum(m$y == 1), 98L)
})

test_that("the step-2 separation check reads the post-fit mu", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  expect_silent(SAIGE:::checkMuSeparation(c(0.01, 0.2, 0.5), "clean.rda", "error"))
  # mu at the boundary is the quantity step 2 divides by, which is why this is
  # checked on the LOADED model and not only on step 1's covariate GLM.
  expect_error(SAIGE:::checkMuSeparation(c(1e-16, 0.2, 0.5), "sep.rda", "error"),
               "at the 0/1 boundary")
  expect_warning(SAIGE:::checkMuSeparation(c(1e-16, 0.2), "sep.rda", "warn"),
                 "at the 0/1 boundary")
  expect_silent(SAIGE:::checkMuSeparation(c(1e-16, 0.2), "sep.rda", "none"))
  expect_error(SAIGE:::checkMuSeparation(c(0.1, 0.2), "x.rda", "nope"),
               "must be one of")
})

test_that("a quantitative model is not gated on cases", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  f <- file.path(saige_goldens, "step1", "quant_sparseGRM.rda")
  skip_if_not(file.exists(f))
  # converged = TRUE so the (trait-agnostic) convergence check is satisfied and
  # what this asserts is only that the CASE gate does not apply.
  quantModel <- list(traitType = "quantitative", converged = TRUE)
  expect_silent(invisible(
    SAIGE:::checkModelEligibility(quantModel, f, 1e6, "error")
  ))
  # ... while the convergence gate DOES apply to a quantitative model.
  quantBad <- list(traitType = "quantitative", converged = FALSE)
  expect_error(SAIGE:::checkModelEligibility(quantBad, f, 1e6, "error"),
               "did NOT converge")
})

# ---- convergence -----------------------------------------------------------
# `converged` is set TRUE only at the two genuine AI-REML exits (LEDGER #123),
# so a loop that ran out at maxiter leaves it FALSE. Nothing in step 2 read
# this field before these gates existed.

test_that("the convergence check fires only on a FALSE flag", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  expect_silent(SAIGE:::checkConvergence(TRUE, "ok.rda", "error"))
  expect_error(SAIGE:::checkConvergence(FALSE, "bad.rda", "error"),
               "did NOT converge")
  expect_warning(SAIGE:::checkConvergence(FALSE, "bad.rda", "warn"),
                 "did NOT converge")
  expect_silent(SAIGE:::checkConvergence(FALSE, "bad.rda", "none"))
  expect_error(SAIGE:::checkConvergence(TRUE, "ok.rda", "maybe"),
               "must be one of")
})

test_that("a model file with no convergence flag is not treated as failed", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  # Absence of the field is not evidence of non-convergence: a model written by
  # a SAIGE that predates it must still be testable.
  expect_output(SAIGE:::checkConvergence(NULL, "old.rda", "error"),
                "no usable convergence flag")
  expect_output(SAIGE:::checkConvergence(NA, "old.rda", "error"),
                "no usable convergence flag")
})

test_that("step 2 refuses a planted non-converged model, and accepts its twin", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  src <- file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")
  e <- new.env()
  load(src, envir = e)

  # the twin: untouched except for the flag, so the ONLY difference between
  # the two arms is the thing under test.
  good <- tempfile(fileext = ".rda")
  modglmm <- e$modglmm
  save(modglmm, file = good)
  expect_warning(SAIGE:::ReadModel(good, chrom = "", LOCO = FALSE,
                                   minCaseCount = 0))

  bad <- tempfile(fileext = ".rda")
  modglmm$converged <- FALSE
  save(modglmm, file = bad)
  expect_error(SAIGE:::ReadModel(bad, chrom = "", LOCO = FALSE,
                                 minCaseCount = 0),
               "did NOT converge")
  # and it is escapable, like every other gate
  expect_warning(SAIGE:::ReadModel(bad, chrom = "", LOCO = FALSE,
                                   minCaseCount = 0, convergenceCheck = "none"))
})

test_that("every committed step-1 golden converged, so this gate is inert", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  # If this ever fails, a golden was regenerated from a fit that did not
  # converge -- which is a finding about the golden, not about the gate.
  for (f in list.files(file.path(saige_goldens, "step1"), pattern = "rda$",
                       full.names = TRUE)) {
    e <- new.env()
    load(f, envir = e)
    expect_true(isTRUE(e$modglmm$converged), info = basename(f))
  }
})

# ---- the separation FRACTION -----------------------------------------------
# The rule this replaced refused on ANY pinned sample, which on real AoU data
# refused a model over fewer than 20 samples in 2,737 -- under 1% of the cohort,
# i.e. under 1% of the effective sample size, which is not variance collapse by
# any reading. The vector below is synthetic.

test_that("a tiny pinned fraction warns; a large one refuses", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  n <- 2737
  mu <- rep(0.05, n)

  # one planted sample of the synthetic 2,737: below the default fraction -> warning, run continues.
  one <- mu
  one[1] <- 1e-16
  expect_warning(SAIGE:::checkMuSeparation(one, "tiny.rda", "error"),
                 "below --separationMaxFraction")

  # 15% pinned (the shape of a real broken model) -> refused.
  many <- mu
  many[seq_len(round(0.15 * n))] <- 1e-16
  expect_error(SAIGE:::checkMuSeparation(many, "broken.rda", "error"),
               "at the 0/1 boundary")

  # and the whole cohort pinned, which is what a collapsed fit looks like
  expect_error(SAIGE:::checkMuSeparation(rep(2.22e-16, n), "dead.rda", "error"),
               "100% of fitted probabilities|at the 0/1 boundary")
})

test_that("the fraction is a knob, and both directions of it are drilled", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  n <- 1000
  mu <- rep(0.05, n)
  mu[seq_len(30)] <- 1e-16          # 3%

  # default 2%: 3% is above it -> refused
  expect_error(SAIGE:::checkMuSeparation(mu, "m.rda", "error"),
               "at the 0/1 boundary")
  # raised to 5%: the SAME vector now only warns
  expect_warning(SAIGE:::checkMuSeparation(mu, "m.rda", "error",
                                           maxFraction = 0.05),
                 "below --separationMaxFraction")
  # lowered to 0: reverts to the old refuse-on-any behaviour
  one <- rep(0.05, n)
  one[1] <- 1e-16
  expect_error(SAIGE:::checkMuSeparation(one, "m.rda", "error",
                                         maxFraction = 0),
               "at the 0/1 boundary")
})

# ---- LEDGER #175 (characterization): the check sees only what is below 1e-8 --
# A separated group is fitted wherever the covariate glm stopped iterating --
# near 1e-8 x deviance / (1.26 x its size) -- so the check counts the part of
# it that happens to be below 1e-8. The two vectors are the shapes the pinned
# fixture's two fits produced (harness/tools/make_pinned_fixture.py): 8
# separated samples fitted just above the threshold, and 60 straddling it.
# Invert when #175 is decided.
test_that("#175: a separated group fitted just above 1e-8 passes silently", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  mu <- rep(0.1, 1000)
  mu[1:8] <- 10^seq(log10(1.7e-8), log10(4.8e-7), length.out = 8)
  expect_silent(SAIGE:::checkMuSeparation(mu, "eight.rda", "error"))
})

test_that("#175: a separated group straddling 1e-8 is reported as the part below it", {
  withr::local_envvar(c(SAIGE_FIT_GATES = "1"))
  mu <- rep(0.1, 1000)
  mu[1:13] <- 10^seq(log10(1.7e-9), log10(9.9e-9), length.out = 13)
  mu[14:60] <- 10^seq(log10(1.1e-8), log10(1.2e-7), length.out = 47)
  expect_warning(SAIGE:::checkMuSeparation(mu, "sixty.rda", "error"),
                 "has 13 of 1000 fitted values")
})

# LEDGER #175, repair (b) (2026-10-01): the effective sample size implied by a
# binary model's fitted values. sum(mu(1-mu)) IS the Fisher information the
# score tests divide by, so the ratio to what a cohort of the same size at the
# same prevalence would carry is the effective-N fraction. It replaces counting
# fitted values below an arbitrary floor: dimensionless, no threshold to tune,
# and stable under prevalence, group size and N.
#
# It is a DIAGNOSTIC -- it refuses nothing and changes no p-value -- so its
# validation is a hand computation plus the invariants, not a null sweep.

test_that("LEDGER #175(b): effective N matches a hand computation", {
  f <- SAIGE:::effectiveSampleSize
  # a balanced, perfectly flat model: every mu equal, so no information is lost
  mu <- rep(0.3, 1000)
  e <- f(mu)
  expect_equal(e$n, 1000L)
  expect_equal(e$neff, 1000, tolerance = 1e-12)
  expect_equal(e$frac, 1, tolerance = 1e-12)

  # a hand case: half at 0.5, half at 0.1
  mu <- c(rep(0.5, 100), rep(0.1, 100))
  e <- f(mu)
  info <- 100 * 0.25 + 100 * 0.09          # 34
  pbar <- 0.3
  expect_equal(e$info, info, tolerance = 1e-12)
  expect_equal(e$pbar, pbar, tolerance = 1e-12)
  expect_equal(e$neff, info / (pbar * (1 - pbar)), tolerance = 1e-12)
  expect_lt(e$frac, 1)                      # dispersion always costs information
})

test_that("LEDGER #175(b): a pinned group costs exactly its share of the information", {
  f <- SAIGE:::effectiveSampleSize
  base <- rep(0.4, 900)
  # 100 samples pinned at ~0 contribute ~0 to sum(mu(1-mu))
  e_pin  <- f(c(base, rep(1e-12, 100)))
  e_none <- f(base)
  # the pinned samples add to n but almost nothing to the information
  expect_equal(e_pin$info, e_none$info, tolerance = 1e-8)
  expect_equal(e_pin$n, 1000L)
  expect_lt(e_pin$frac, e_none$frac)
  # and the loss is visible as a fraction without any threshold being chosen
  expect_lt(e_pin$frac, 0.95)
})

test_that("LEDGER #175(b): it is scale-free in N and defined only for binary", {
  f <- SAIGE:::effectiveSampleSize
  # the same mu spectrum at 10x the size gives the same FRACTION
  spec <- c(rep(0.5, 10), rep(0.05, 10))
  expect_equal(f(spec)$frac, f(rep(spec, 10))$frac, tolerance = 1e-12)
  # degenerate inputs return NULL rather than NaN
  expect_null(f(numeric(0)))
  expect_null(f(rep(0, 100)))               # pbar = 0, denominator vanishes
  expect_null(f(rep(1, 100)))
  # quantitative traits report nothing at all
  expect_null(SAIGE:::reportEffectiveN(spec, "quantitative", "lbl"))
})
