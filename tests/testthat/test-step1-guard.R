# LEDGER row 122, the guard around the binary AI-REML step. The rules live in
# the TauStepGuard source file and are pinned here on real fits of the bundled
# N = 1,000 fixture (100 families of 10, the sparse GRM in extdata/output).
#
# 1. A healthy fit is untouched. The bundled binary fit converges in four
#    undamped iterations with one direction reversal; the guard arms on the
#    second, so the printed trace must equal the committed golden log's to
#    the last digit and the iteration count must be 4 (not 5: the twice rule
#    applies only once armed).
# 2. A fit that cycled: planted from the fixture's own GRM, 3% prevalence
#    with a family effect of variance 1 on the logit scale, whose plain
#    iteration swings 1.47, 1.70, 1.13, 0.67, ... and reported 1.14 / 1.04 /
#    1.99 from starts 0 / 0.3 / 2 before the guard (the 1.99 is a cold start
#    accepted on its first iteration). The dense PQL fixed point is 1.056.
#    Under the guard and the WRONG AI curvature (before LEDGER #124) the three
#    starts agreed near 1.05; with the curvature fixed the same guard never
#    damped here and the starts reported 1.280 / 1.236 / 0.891 (LEDGER #125).
#    The mechanism was the inner (alpha, b) pass converging only to the shared
#    tol = 0.02, which left the outer REML score four times noisier than the
#    exact score; with the inner pass at tolCoef = 1e-5 (the #125 fix) the
#    three starts agree at 1.085-1.094 and --tolCoef=0.02 reproduces the old
#    start-dependent stop, which is how the flag is pinned below.
# A planted trait whose fixed point is the boundary (score at 0 negative)
# must still stop at 0, and now says so.
#
# Tier: standard (~19 step-1 fits of ~3 s, half of them in-process controls).

# The helpers below call run_saige_cli / expect_cli_ok / saige_extdata from
# helper-saige.R, which lintr does not see; only that linter is silenced.
# nolint start: object_usage_linter.
sg_grm <- function() {
  f <- file.path(saige_extdata, "output",
                 "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  as.matrix(methods::as(Matrix::readMM(f), "generalMatrix"))
}

# The planted phenotypes, seeds fixed; the recipe is the one the ledger row's
# measurements used (harness/work/flexrv_study, 2026-09-10).
sg_pheno_file <- local({
  path <- NULL
  function() {
    if (!is.null(path)) return(path)
    g <- sg_grm()
    n <- nrow(g)
    idFile <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
    ids <- readLines(file.path(saige_extdata, "output", idFile))
    l <- t(chol(g + diag(1e-8, n)))
    set.seed(31)
    x1 <- stats::rnorm(n)
    x2 <- stats::rbinom(n, 1, 0.5)
    ph <- data.frame(IID = ids, x1 = x1, x2 = x2)
    set.seed(32) # binary traits with an interior fixed point: b_tau0.5
    for (t2 in c(0, 0.5, 1, 2)) {
      b <- as.vector(l %*% stats::rnorm(n)) * sqrt(t2)
      eta <- -1.5 + 0.5 * x1 + 0.3 * x2 + b
      ph[[paste0("b_tau", t2)]] <- stats::rbinom(n, 1, stats::plogis(eta))
    }
    set.seed(33) # low prevalence: c_K03_t1 cycles, c_K03_t3's fixed point is the boundary
    for (k in c(0.03, 0.05)) {
      for (t2 in c(1, 3, 6)) {
        b <- as.vector(l %*% stats::rnorm(n)) * sqrt(t2)
        lp <- 0.5 * x1 + 0.3 * x2 + b
        y <- as.integer(lp + stats::rlogis(n) > stats::quantile(lp + stats::rlogis(n), 1 - k))
        ph[[sprintf("c_K%02d_t%d", round(100 * k), t2)]] <- y
      }
    }
    path <<- tempfile("planted_", fileext = ".txt")
    utils::write.table(ph, path, quote = FALSE, row.names = FALSE)
    path
  }
})

sg_fit <- function(pheno, col, covar, qcovar, tau_init = NULL, extra = character()) {
  prefix <- file.path(tempfile(paste0("guard_", col, "_")), "fit")
  dir.create(dirname(prefix), recursive = TRUE)
  st <- run_saige_cli("step1_fitNULLGLMM.R", c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed",
           ".sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE", "--useSparseGRMforVarRatio=TRUE",
    "--skipVarianceRatioEstimation=TRUE", "--traitType=binary",
    paste0("--phenoFile=", pheno), paste0("--phenoCol=", col),
    paste0("--covarColList=", covar), paste0("--qCovarColList=", qcovar),
    "--sampleIDColinphenoFile=IID", "--IsOverwriteVarianceRatioFile=TRUE",
    "--maxiter=40",
    if (!is.null(tau_init)) paste0("--tauInit=1,", tau_init),
    extra, # last value wins in optparse (measured 2026-09-10), so a caller's --maxiter overrides
    paste0("--outputPrefix=", prefix)
  ))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  e <- new.env()
  load(paste0(prefix, ".rda"), e)
  trace <- as.numeric(sub("^tau: +1 +", "", grep("^tau: ", log, value = TRUE)))
  list(theta = e$modglmm$theta, converged = e$modglmm$converged,
       iters = sum(grepl("^Iteration", log)), trace = trace, log = log)
}

# The guard-disabled CONTROL. The same fit in a child process with the CLI's
# conventions (covariates as an offset, relatednessCutoff 0, which is what
# makes it reproduce sg_fit's theta exactly), once as installed and once with
# tauGuardApply replaced in the namespace by a pass-through. "The guard is
# inert on this fit" then means the two .rda files are byte-identical, at any
# curvature; a positive control (a fit where the guard does act) shows the
# replacement takes effect. Returns the .rda path and theta.
sg_fit_inproc <- function(pheno, col, covar, qcovar, guard, offset = TRUE, tol_coef = NULL) {
  prefix <- file.path(tempfile(paste0("guardctl_", col, "_")), "fit")
  dir.create(dirname(prefix), recursive = TRUE)
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    if (!guard) c('ns <- asNamespace("SAIGE"); unlockBinding("tauGuardApply", ns)',
                  'assign("tauGuardApply",',
                  "       function(state, tau0, raw, step, ...) list(p = step, state = state), ns)",
                  'lockBinding("tauGuardApply", ns)'),
    paste0('sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed',
           '.sparseGRM.mtx"'),
    paste0('fitNULLGLMM(plinkFile = "./input/nfam_100_nindep_0_step1_',
           'includeMoreRareVariants_poly_22chr",'),
    '  sparseGRMFile = sg, sparseGRMSampleIDFile = paste0(sg, ".sampleIDs.txt"),',
    "  useSparseGRMtoFitNULL = TRUE, useSparseGRMforVarRatio = TRUE,",
    sprintf('  skipVarianceRatioEstimation = TRUE, phenoFile = "%s", phenoCol = "%s",', pheno, col),
    sprintf("  covarColList = c(%s), qCovarCol = c(%s),",
            paste0('"', strsplit(covar, ",")[[1]], '"', collapse = ", "),
            paste0('"', strsplit(qcovar, ",")[[1]], '"', collapse = ", ")),
    '  sampleIDColinphenoFile = "IID", traitType = "binary", IsOverwriteVarianceRatioFile = TRUE,',
    sprintf("  maxiter = 40, isCovariateOffset = %s, relatednessCutoff = 0,",
            if (offset) "TRUE" else "FALSE"),
    if (!is.null(tol_coef)) sprintf("  tolCoef = %s,", format(tol_coef)),
    sprintf('  outputPrefix = "%s")', prefix)
  ), script)
  owd <- setwd(saige_extdata)
  on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = FALSE, stderr = FALSE))
  expect_identical(as.integer(st), 0L)
  e <- new.env()
  load(paste0(prefix, ".rda"), e)
  list(rda = paste0(prefix, ".rda"), theta = e$modglmm$theta)
}

# Guard live vs guard neutralised, and the CLI fit they must agree with.
expect_guard_inert <- function(pheno, col, covar, qcovar, cli_theta2) {
  g <- sg_fit_inproc(pheno, col, covar, qcovar, guard = TRUE)
  u <- sg_fit_inproc(pheno, col, covar, qcovar, guard = FALSE)
  expect_files_byte_identical(g$rda, u$rda)
  expect_equal(g$theta[2], cli_theta2, tolerance = 1e-8)
}
# nolint end

test_that("a healthy fit is untouched: the bundled binary fit's trace equals the golden log's", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  gl <- readLines(file.path(saige_goldens, "logs", "step1_binary_sparseGRM_a.log"), warn = FALSE)
  golden <- as.numeric(sub("^tau: +1 +", "", grep("^tau: ", gl, value = TRUE)))
  f <- sg_fit("./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt", "y_binary",
              "x1,x2,a9,a10", "a9,a10")
  expect_equal(f$iters, length(golden))
  expect_equal(f$trace, golden, tolerance = 1e-6)
  expect_equal(f$theta[2], golden[length(golden)], tolerance = 1e-6)
  # the control: with the guard neutralised the fit is the same file
  expect_guard_inert("./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
                     "y_binary", "x1,x2,a9,a10", "a9,a10", f$theta[2])
})

test_that("the guard-disabled control is not vacuous: a fit the guard acts on differs", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # the planted 3% trait under the FULL covariate model (offset = FALSE) with
  # the OLD inner tolerance (tolCoef = 0.02) is a fit where the guard damps:
  # 1.089 with it, 1.208 without (2026-09-10)
  ph <- sg_pheno_file()
  g <- sg_fit_inproc(ph, "c_K03_t1", "x1,x2", "x2", guard = TRUE, offset = FALSE, tol_coef = 0.02)
  u <- sg_fit_inproc(ph, "c_K03_t1", "x1,x2", "x2", guard = FALSE, offset = FALSE, tol_coef = 0.02)
  expect_false(tools::md5sum(g$rda) == tools::md5sum(u$rda))
  expect_gt(abs(g$theta[2] - u$theta[2]), 0.05)
})

test_that("#125: with the inner pass converged, the planted 3% trait agrees from three starts", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  ph <- sg_pheno_file()
  fits <- lapply(list(NULL, 0.3, 2), function(t0) sg_fit(ph, "c_K03_t1", "x1,x2", "x2", t0))
  tau <- vapply(fits, function(f) f$theta[2], numeric(1))
  # no fit may accept its first proposal (the cold start at 2 did, before #122)
  expect_true(all(vapply(fits, function(f) f$iters, numeric(1)) > 1))
  expect_true(all(vapply(fits, function(f) isTRUE(f$converged), logical(1))))
  # 1.091 / 1.094 / 1.085 in 5 / 3 / 3 iterations (2026-09-10): the dense PQL
  # fixed point is 1.056 and the trace estimator's se here is ~0.025
  expect_lt(max(tau) / min(tau) - 1, 0.03)
  expect_true(all(abs(tau - 1.056) < 0.1))
  # (the guard still acts mildly on the default start here, 1.0914 with it and
  # 1.0878 without, both inside the trace noise; inertness is the HEALTHY fits'
  # claim, tested above and below, not this trait's)
  # the knob is wired: the old inner tolerance reproduces the old stop on a
  # single pair inside 2% while the trace still swings (7 iterations, 1.2797)
  old <- sg_fit(ph, "c_K03_t1", "x1,x2", "x2", extra = "--tolCoef=0.02")
  expect_equal(old$iters, 7)
  expect_equal(old$theta[2], 1.2797, tolerance = 1e-3)
})

test_that("a planted trait whose fixed point is the boundary stops at 0 and says why", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  f <- sg_fit(sg_pheno_file(), "c_K03_t3", "x1,x2", "x2")
  expect_equal(f$theta[2], 0)
  expect_true(any(grepl("tau2 = 0 is the optimum", f$log)))
})

test_that("a healthy planted trait keeps its pre-guard fit exactly", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  f <- sg_fit(sg_pheno_file(), "b_tau0.5", "x1,x2", "x2")
  # 3 iterations and 0.1475275 before LEDGER #124 corrected the AI curvature,
  # 4 and 0.14645 after it, 2 and 0.14676 with the inner pass at tolCoef =
  # 1e-5 (#125; trace 0.1411 0.1468, no reversal, no damping)
  expect_equal(f$iters, 2)
  expect_equal(f$theta[2], 0.14676, tolerance = 1e-4)
  expect_guard_inert(sg_pheno_file(), "b_tau0.5", "x1,x2", "x2", f$theta[2])
})

test_that("#123: a fit converging on its last allowed iteration is flagged converged", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # the planted 3% trait meets the criterion at iteration 5 from the default
  # start (12 before LEDGER #124, 7 before #125)
  f <- sg_fit(sg_pheno_file(), "c_K03_t1", "x1,x2", "x2", extra = "--maxiter=5")
  expect_equal(f$iters, 5)
  expect_true(f$theta[2] > 0.5) # it did converge: the iterate is at the fixed point
  expect_true(f$converged)
  # the contrast: one iteration fewer and the flag says so
  g <- sg_fit(sg_pheno_file(), "c_K03_t1", "x1,x2", "x2", extra = "--maxiter=4")
  expect_equal(g$iters, 4)
  expect_false(g$converged)
})
