# Get_Variance_Ratio (R/readInGLMM.R): the parser of the step1->step2
# variance-ratio contract. All 6 args are positionally required; isGroupTest
# and useSparseGRMtoFitNULL are accepted but unused.

vr_file <- function(name) file.path(saige_goldens, "step1", paste0(name, ".varianceRatio.txt"))
CATE8 <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)
CATE7 <- c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)

test_that("modern 3-column file: categorical null + null_noXadj, no sparse", {
  skip_if_no_goldens()
  out <- Get_Variance_Ratio(vr_file("binary_sparseGRM_a"), c(10, 20), c(20), FALSE, FALSE, FALSE)
  expect_identical(out$ratioVec_sparse, c(-1))
  expect_equal(out$ratioVec_null, c(0.947747869047814, 0.942377359586025), tolerance = 1e-14)
  expect_equal(out$ratioVec_null_noXadj, c(0.939781574541214, 0.932242806034484), tolerance = 1e-14)
})

test_that("isSparseGRM=TRUE with no sparse rows is a hard error", {
  skip_if_no_goldens()
  # typos "estimatedin"/"Pleae" are verbatim in the source
  expect_error(
    Get_Variance_Ratio(vr_file("binary_sparseGRM_a"), c(10, 20), c(20), FALSE, TRUE, FALSE),
    "sparse GRM is specified but the variance ratio for sparse GRM was not estimatedin Step 1")
})

test_that("single-ratio file skips categorical validation", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("0.871784619442486 null 1", "0.868978655998129 null_noXadj 1"), f)
  out <- Get_Variance_Ratio(f, c(10), c(), FALSE, FALSE, FALSE)
  expect_equal(out$ratioVec_null, 0.871784619442486, tolerance = 1e-14)
  expect_identical(out$ratioVec_sparse, c(-1))
  # wrong-length cate vecs must NOT error when there is only one ratio
  expect_no_error(Get_Variance_Ratio(f, c(1, 2, 3), c(1, 2), FALSE, FALSE, FALSE))
})

test_that("categorical count validation", {
  skip_if_no_goldens()
  f <- vr_file("binary_sparseGRM_b")  # 8 null rows
  expect_no_error(Get_Variance_Ratio(f, CATE8, CATE7, FALSE, FALSE, FALSE))
  expect_error(Get_Variance_Ratio(f, CATE8[-1], CATE7, FALSE, FALSE, FALSE),
               "The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude")
  # "lenght" typo is verbatim in the source
  expect_error(Get_Variance_Ratio(f, CATE8, CATE7[-1], FALSE, FALSE, FALSE),
               "does not match with the lenght of cateVarRatioMinMACVecExclude")
})

test_that("empty path returns GRM-mode sentinels; missing path errors", {
  expect_identical(Get_Variance_Ratio("", c(10), c(), FALSE, TRUE, FALSE),
                   list(ratioVec_sparse = c(1), ratioVec_null = c(-1), ratioVec_null_noXadj = c(-1)))
  expect_identical(Get_Variance_Ratio("", c(10), c(), FALSE, FALSE, FALSE),
                   list(ratioVec_sparse = c(-1), ratioVec_null = c(1), ratioVec_null_noXadj = c(1)))
  expect_error(Get_Variance_Ratio("/nope/x.txt", c(10), c(), FALSE, FALSE, FALSE),
               "varianceRatioFile is specified but the file /nope/x.txt does not exist")
})

test_that("LEDGER #13 fixed: legacy (<1.0.6) files are rejected cleanly", {
  f <- tempfile(); writeLines("0.95", f)
  expect_error(Get_Variance_Ratio(f, c(10), c(), FALSE, FALSE, FALSE),
               "pre-1.0.6 legacy format")
  expect_error(Get_Variance_Ratio(f, c(10), c(), FALSE, TRUE, FALSE),
               "pre-1.0.6 legacy format")
})

test_that("missing null_noXadj rows yield numeric(0) from the parser (gated at run level)", {
  # Get_Variance_Ratio itself still returns numeric(0) (group tests need
  # that); SPAGMMATtest now rejects it for single-variant runs (LEDGER #10).
  f <- tempfile(); writeLines(c("0.95 null 1", "0.92 null 2"), f)
  out <- Get_Variance_Ratio(f, c(10, 20), c(20), FALSE, FALSE, FALSE)
  expect_identical(out$ratioVec_null_noXadj, numeric(0))
  expect_equal(out$ratioVec_null, c(0.95, 0.92))
})
