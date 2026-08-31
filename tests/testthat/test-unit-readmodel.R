# ReadModel + friends (R/readInGLMM.R): the step1 -> step2 model-loading
# contract (implicit and name-based; see CLAUDE.md).

golden_rda <- function(name) file.path(saige_goldens, "step1", paste0(name, ".rda"))

load_golden <- function(name) {
  e <- new.env(); load(golden_rda(name), envir = e); get("modglmm", envir = e)
}

test_that("ReadModel on the quantitative golden", {
  skip_if_no_goldens()
  raw <- load_golden("quant_sparseGRM")
  m <- ReadModel(golden_rda("quant_sparseGRM"), chrom = "", LOCO = FALSE)
  expect_null(m$Y); expect_null(m$linear.predictors)
  expect_null(m$coefficients); expect_null(m$cov)
  expect_equal(m$mu, as.vector(raw$fitted.values))
  expect_equal(m$mu2, rep(1 / raw$theta[1], length(m$mu)))
  expect_identical(m$Sigma_iXXSigma_iX, matrix(0, nrow = 1, ncol = 1))
  expect_identical(length(m$offset), nrow(m$X))
})

test_that("ReadModel on the binary golden (needs SKAT)", {
  skip_if_no_goldens()
  skip_if_not(requireNamespace("SKAT", quietly = TRUE), "SKAT not installed")
  m <- ReadModel(golden_rda("binary_sparseGRM_a"), chrom = "", LOCO = FALSE)
  expect_equal(m$mu2, m$mu * (1 - m$mu))
  expect_equal(m$obj_cc$pi_1, m$mu2)
  expect_equal(m$obj_cc$mu, m$mu)
})

test_that("ReadModel argument and LOCO validation", {
  skip_if_no_goldens()
  expect_error(ReadModel("/nope.rda"), "does not exsit")  # typo verbatim
  expect_error(ReadModel(golden_rda("quant_sparseGRM"), LOCO = "yes"),
               "LOCO should be TRUE or FALSE")
  # all goldens were fit without LOCO
  expect_error(ReadModel(golden_rda("quant_sparseGRM"), chrom = "1", LOCO = TRUE),
               "LOCO is TRUE but the null model file .rda does not contain LOCO results")
})

test_that("ReadModel LOCO=TRUE picks the per-chromosome result (synthetic fixture)", {
  n <- 8
  f0 <- rnorm(n); r0 <- rnorm(n); f3 <- rnorm(n); r3 <- rnorm(n)
  LOCOResult <- vector("list", 22)
  LOCOResult[[3]] <- list(fitted.values = f3, residuals = r3,
                          obj.noK = list(marker = "chr3"))
  modglmm <- list(LOCO = TRUE, traitType = "quantitative", theta = c(2, 1),
                  y = rnorm(n), X = matrix(1, n, 1),
                  fitted.values = f0, residuals = r0, LOCOResult = LOCOResult)
  f <- tempfile(fileext = ".rda"); save(modglmm, file = f)

  m <- ReadModel(f, chrom = "chr3", LOCO = TRUE)
  expect_equal(m$mu, as.vector(f3))
  expect_equal(m$residuals, r3)
  expect_identical(m$obj.noK$marker, "chr3")
  expect_true(all(vapply(m$LOCOResult, is.null, logical(1))))

  expect_error(ReadModel(f, chrom = "", LOCO = TRUE),
               "chrom needs to be specified in order to apply Leave-one-chromosome-out")
})

test_that("getChromNumber", {
  expect_error(getChromNumber(""), "chrom is not specified")
  expect_identical(getChromNumber("chr5"), 5)
  expect_identical(getChromNumber("CHR5"), 5)
  expect_identical(getChromNumber("5"), 5)
  suppressWarnings({
    expect_identical(getChromNumber("X"), "X")     # out of 1..22 -> original string
    expect_identical(getChromNumber("23"), "23")
  })
})

test_that("removeLOCOResult", {
  obj <- list(LOCOResult = as.list(1:22))
  out <- removeLOCOResult(setdiff(1:22, 3), obj)
  expect_identical(out$LOCOResult[[3]], 3L)
  expect_identical(length(out$LOCOResult), 22L)
  expect_true(all(vapply(out$LOCOResult[-3], is.null, logical(1))))
})

test_that("ReadModel_subsample repeats ReadModel's LOCO validation", {
  # ReadModel_subsample (reached from SPAGMMATtest when --subSampleFile is
  # given) carries its own copy of the LOCO checks. The copies in ReadModel
  # are covered above; these are the duplicates, which no golden run reaches
  # because no golden passes --subSampleFile. Keep the two in step.
  set.seed(101)
  n <- 8
  ids <- paste0("S", seq_len(n))
  fixture <- function(loco) {
    modglmm <- list(LOCO = loco, traitType = "quantitative", theta = c(2, 1),
                    y = rnorm(n), X = matrix(1, n, 1),
                    fitted.values = rnorm(n), residuals = rnorm(n),
                    sampleID = ids, LOCOResult = vector("list", 22))
    f <- tempfile(fileext = ".rda"); save(modglmm, file = f); f
  }
  ss <- tempfile(); writeLines(ids[1:5], ss)

  expect_error(ReadModel_subsample(fixture(FALSE), chrom = "", LOCO = "yes",
                                   subSampleFile = ss),
               "LOCO should be TRUE or FALSE.", fixed = TRUE)
  expect_error(ReadModel_subsample(fixture(FALSE), chrom = "1", LOCO = TRUE,
                                   subSampleFile = ss),
               "does not contain LOCO results")
  expect_error(ReadModel_subsample(fixture(TRUE), chrom = "", LOCO = TRUE,
                                   subSampleFile = ss),
               "chrom needs to be specified in order to apply Leave-one-chromosome-out")
})

test_that("ReadModel_subsample validation and prefix-subset happy path", {
  skip_if_no_goldens()
  expect_error(ReadModel_subsample(golden_rda("quant_sparseGRM"), chrom = "", LOCO = FALSE,
                                   subSampleFile = "/nope"), "does not exist")

  raw <- load_golden("quant_sparseGRM")
  bad <- tempfile(); writeLines(paste0("ZZZ", raw$sampleID[1:5]), bad)
  expect_error(ReadModel_subsample(golden_rda("quant_sparseGRM"), chrom = "", LOCO = FALSE,
                                   subSampleFile = bad),
               "No samples analyzed in Step 1 are included")

  k <- 500  # 100 leaves the covariate matrix singular in ScoreTest_NULL_Model
  keep <- tempfile(); writeLines(raw$sampleID[seq_len(k)], keep)
  m <- ReadModel_subsample(golden_rda("quant_sparseGRM"), chrom = "", LOCO = FALSE,
                           subSampleFile = keep)
  expect_identical(length(m$sampleID), as.integer(k))
  expect_identical(nrow(m$X), as.integer(k))
  expect_identical(m$Sigma_iXXSigma_iX, matrix(0, nrow = 1, ncol = 1))
})
