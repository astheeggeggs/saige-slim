# LEDGER #124: the covariate projection inside the AI term.
#
# P v = Sigma^-1 v - Sigma^-1 X (X' Sigma^-1 X)^-1 X' Sigma^-1 v, and the
# average-information term is APY' P APY with APY = G P Y. The kernels
# getAIScore (binary) and getAIScore_q (quantitative, the 2 x 2 matrix) solve
# Sigma^-1 APY by PCG and then applied the projection to THAT solved vector,
# X' Sigma^-2 APY, instead of to APY; PY1 two lines above and GetTrace project
# correctly. Only the curvature was wrong, so every Newton step length was off
# by under 1% on the bundled fixture and no fixed point moved.
#
# The kernels read the GRM from C++ globals, so each trait is fitted in a child
# process (the sparse-GRM fit of the bundled N = 1,000 fixture, ~3 s) and the
# kernel is then called at the fitted theta with the solved covariate vectors
# computed DENSELY in R from the same W, tau and GRM. Two dense evaluations of
# the AI term at those inputs, the projection applied to APY and applied to
# Sigma^-1 APY, bracket the kernel. Measured 2026-09-10: at unpatched HEAD the
# kernel equalled the second to 1e-8 and missed the first by 6.4e-3 (binary,
# tau2 = 0.39) and by 1.3e-3 / 2.6e-4 / 4.7e-4 (the quantitative entries);
# fixed, it equals the first to 1e-7. The control is the kernel's returned PY
# against the dense PY (1e-7), which fails if the sample order or Sigma were
# not the kernel's; the two dense values must themselves differ by more than
# the "differs" margin, or the fixture could not tell the formulas apart.
#
# Tier: standard (two child fits).

# nolint start: object_usage_linter.
ap_kernel_state <- function(trait) {
  out <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  pre <- tempfile("ap124_")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('trait <- "%s"; out <- "%s"; pre <- "%s"', trait, out, pre),
    'ns <- asNamespace("SAIGE"); rec <- new.env()',
    # the design matrix and W the kernel saw: recorded from its last call
    'for (fn in c("getAIScore", "getAIScore_q")) local({',
    "  f <- fn; orig <- get(f, ns)",
    "  wrap <- function(...) { re <- orig(...); rec$args <- list(...); re }",
    "  unlockBinding(f, ns); assign(f, wrap, ns); lockBinding(f, ns) })",
    paste0('sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed',
           '.sparseGRM.mtx"'),
    paste0('common <- list(plinkFile = "./input/nfam_100_nindep_0_step1_',
           'includeMoreRareVariants_poly_22chr",'),
    '  sparseGRMFile = sg, sparseGRMSampleIDFile = paste0(sg, ".sampleIDs.txt"),',
    "  useSparseGRMtoFitNULL = TRUE,",
    '  phenoFile = "./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",',
    '  covarColList = c("x1", "x2", "a9", "a10"), qCovarCol = c("a9", "a10"),',
    # isCovariateOffset = FALSE is the configuration the margins below were
    # measured on (2026-09-10). LEDGER #151 (2026-09-15) moved the API default
    # to TRUE, the CLI's value; with the covariates as an offset X loses those
    # columns and the quantitative fixture's right-vs-wrong margins fall to
    # 4.7e-5 / 5.6e-5 / 9.97e-5, under the 1e-4 "differs" bar (the kernel-vs-
    # right checks at 1e-6 still passed). Named, so a default cannot move it.
    '  sampleIDColinphenoFile = "IID", IsOverwriteVarianceRatioFile = TRUE, outputPrefix = pre,',
    "  isCovariateOffset = FALSE)",
    'spec <- if (trait == "binary") list(phenoCol = "y_binary", traitType = "binary") else',
    '  list(phenoCol = "y_quantitative", traitType = "quantitative", invNormalize = TRUE)',
    'do.call(fitNULLGLMM, c(common, spec)); load(paste0(pre, ".rda"))',
    "X <- rec$args[[2]]; W <- as.numeric(rec$args[[3]])",
    "Y <- as.numeric(modglmm$Y); tau <- as.numeric(modglmm$theta)",
    'if (trait == "binary") { mu <- as.numeric(modglmm$fitted.values); W <- mu * (1 - mu) }',
    'G <- as.matrix(methods::as(Matrix::readMM(sg), "generalMatrix"))',
    'o <- match(modglmm$sampleID, readLines(paste0(sg, ".sampleIDs.txt"))); G <- G[o, o]',
    "Sigma <- tau[1] * diag(1 / W) + tau[2] * G",
    "SiX <- solve(Sigma, X); cov <- solve(crossprod(X, SiX)); SiY <- as.numeric(solve(Sigma, Y))",
    'kern <- get(if (trait == "binary") "getAIScore" else "getAIScore_q", ns)',
    "re <- kern(Y, X, W, tau, SiY, SiX, cov, 30L, 500L, 1e-5, 0.0025)",
    "saveRDS(list(X = X, W = W, Y = Y, tau = tau, G = G, SiX = SiX, cov = cov, SiY = SiY,",
    "             re = re), out)"
  ), script)
  log <- tempfile(fileext = ".log")
  owd <- setwd(saige_extdata)
  on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = log, stderr = log))
  attr(st, "log") <- log
  expect_cli_ok(st)
  readRDS(out)
}

# The two dense evaluations of v' P v-type terms at the kernel's own inputs.
ap_dense <- function(s) {
  Sigma <- s$tau[1] * diag(1 / s$W) + s$tau[2] * s$G
  proj <- function(v, to_apy) {
    Sv <- solve(Sigma, v)
    Sv - s$SiX %*% (s$cov %*% (t(s$SiX) %*% (if (to_apy) v else Sv)))
  }
  PY <- as.numeric(s$SiY - s$SiX %*% (s$cov %*% crossprod(s$X, s$SiY)))
  APY <- s$G %*% PY
  A0PY <- PY / s$W
  list(PY = PY,
       right = c(a00 = sum(A0PY * proj(A0PY, TRUE)), a11 = sum(APY * proj(APY, TRUE)),
                 a01 = sum(A0PY * proj(APY, TRUE))),
       wrong = c(a00 = sum(A0PY * proj(A0PY, FALSE)), a11 = sum(APY * proj(APY, FALSE)),
                 a01 = sum(A0PY * proj(APY, FALSE))))
}

ap_rel <- function(a, b) abs(a - b) / abs(b)

ap_check <- function(code, d, entries) {
  for (k in entries) {
    expect_gt(ap_rel(d$right[[k]], d$wrong[[k]]), 1e-4, label = paste("fixture discriminates", k))
    expect_lt(ap_rel(code[[k]], d$right[[k]]), 1e-6,
              label = paste("kernel vs projection of APY", k))
    expect_gt(ap_rel(code[[k]], d$wrong[[k]]), 1e-4,
              label = paste("kernel vs projection of Sigma^-1 APY", k))
  }
}

test_that("the binary AI term projects APY, not Sigma^-1 APY", {
  skip_unless_level("standard")
  s <- ap_kernel_state("binary")
  d <- ap_dense(s)
  expect_lt(max(abs(as.numeric(s$re$PY) - d$PY)) / max(abs(d$PY)), 1e-5)
  ap_check(c(a11 = s$re$AI), d, "a11")
})

test_that("the quantitative AI matrix projects A0PY and APY, not their solves", {
  skip_unless_level("standard")
  s <- ap_kernel_state("quant")
  d <- ap_dense(s)
  expect_lt(max(abs(as.numeric(s$re$PY) - d$PY)) / max(abs(d$PY)), 1e-5)
  expect_equal(s$re$AI[1, 2], s$re$AI[2, 1])
  ap_check(c(a00 = s$re$AI[1, 1], a11 = s$re$AI[2, 2], a01 = s$re$AI[1, 2]), d,
           c("a00", "a11", "a01"))
})
# nolint end
