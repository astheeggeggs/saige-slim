# LEDGER row 122, the guard around the binary AI-REML step. The rules live in
# the TauStepGuard source file and are pinned here on real fits of the bundled
# N = 1,000 fixture (100 families of 10, the sparse GRM in extdata/output).
#
# 1. A healthy fit is untouched. The bundled binary fit converges in four
#    undamped iterations with one direction reversal; the guard arms on the
#    second, so the printed trace must equal the committed golden log's to
#    the last digit and the iteration count must be 4 (not 5: the twice rule
#    applies only once armed).
# 2. A fit that cycled reports the same tau2 from different starts. Planted
#    from the fixture's own GRM: 3% prevalence with a family effect of
#    variance 1 on the logit scale, whose plain iteration swings 1.47, 1.70,
#    1.13, 0.67, ... and reported 1.14 / 1.04 / 1.99 from starts 0 / 0.3 / 2
#    before the guard (the 1.99 is a cold start accepted on its first
#    iteration). The dense PQL fixed point is 1.056; the guarded fits agree
#    with each other and sit within the trace estimator's own noise of it.
# A planted trait whose fixed point is the boundary (score at 0 negative)
# must still stop at 0, and now says so.
#
# Tier: standard (~8 step-1 fits of ~3 s).

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
    extra,
    paste0("--outputPrefix=", prefix)
  ))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  e <- new.env()
  load(paste0(prefix, ".rda"), e)
  trace <- as.numeric(sub("^tau: +1 +", "", grep("^tau: ", log, value = TRUE)))
  list(theta = e$modglmm$theta, iters = sum(grepl("^Iteration", log)), trace = trace, log = log)
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
})

test_that("#122: the planted 3% trait reports the same tau2 from three starts", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  ph <- sg_pheno_file()
  fits <- lapply(list(NULL, 0.3, 2), function(t0) sg_fit(ph, "c_K03_t1", "x1,x2", "x2", t0))
  tau <- vapply(fits, function(f) f$theta[2], numeric(1))
  # no fit may accept its first proposal (the cold start at 2 did, before)
  expect_true(all(vapply(fits, function(f) f$iters, numeric(1)) > 1))
  expect_lt(max(tau) / min(tau) - 1, 0.03)
  # the dense PQL fixed point is 1.056; the trace estimator's se here is ~0.025
  expect_true(all(abs(tau - 1.056) < 0.1))
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
  expect_equal(f$iters, 3)
  expect_equal(f$theta[2], 0.1475275, tolerance = 1e-6) # measured before the guard, 2026-09-10
})
