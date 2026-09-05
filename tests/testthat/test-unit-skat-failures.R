# LEDGER #57, the half left open by the chunk-file fix: a SKAT kernel error
# inside get_SKAT_pvalue used to become NA p-values and NOTHING else -- an
# error.code set and never read. Policy chosen 2026-09-05: keep the NA rows,
# carry the kernel's message out, count, and summarise at the end of the run
# (R/SkatFailures.R). These are the seams, fast tier; the end-to-end run with a
# planted kernel failure is test-step2-skat-failures.R.

test_that("get_SKAT_pvalue carries the kernel's error message out; success carries none", {
  # a NaN Phi makes SKAT:::Met_SKAT_Get_Pvalue error ("missing value where
  # TRUE/FALSE needed"). A ZERO Phi does not -- SKAT warns "No polymorphic
  # SNPs!" and returns p = 1 -- which is why the end-to-end test has to plant
  # its failure rather than construct one from a degenerate weight.
  bad <- suppressWarnings(get_SKAT_pvalue(Score = c(1, 2, 3), Phi = matrix(NaN, 3, 3),
                                          r.corr = 0, regionTestType = "SKAT-O"))
  expect_true(is.na(bad$Pvalue_SKATO))
  expect_true(is.na(bad$Pvalue_Burden))
  expect_true(is.na(bad$Pvalue_SKAT))
  expect_type(bad$error, "character")
  expect_length(bad$error, 1)
  expect_gt(nchar(bad$error), 0)
  expect_false(grepl("[\r\n]", bad$error)) # one line, so the sidecar stays one row per cell

  set.seed(1)
  A <- matrix(rnorm(9), 3)
  good <- get_SKAT_pvalue(Score = c(1, 2, 3), Phi = crossprod(A) + diag(3),
                          r.corr = 0, regionTestType = "SKAT-O")
  expect_null(good$error)
  expect_false("error" %in% names(good))   # the success shape is unchanged
  expect_true(is.finite(good$Pvalue_SKATO))
})

test_that("the accumulator counts, summarises, and writes the sidecar only when non-empty", {
  skatFailureReset()
  expect_identical(nrow(skatFailureTable()), 0L)
  out <- tempfile(fileext = ".txt")
  msg <- capture.output(tab <- skatFailureSummary(out))
  expect_identical(msg, "SKAT failures: 0 cells returned NA from the SKAT kernel.")
  expect_false(file.exists(paste0(out, ".skatFailures.txt")))

  skatFailureRecord("GENE1", "lof", 0.001, NULL, "marginal", "boom")
  skatFailureRecord("GENE1", "lof", 0.01, "CADD", "conditional", "boom")
  skatFailureRecord("GENE7", "missense;lof", 0.01, NA, "marginal", "other")
  tab <- skatFailureTable()
  expect_identical(dim(tab), c(3L, 6L))
  expect_identical(names(tab), c("Region", "Group", "max_MAF", "Weight", "Test", "Error"))
  expect_identical(tab$Weight, c(NA, "CADD", NA))
  expect_identical(tab$max_MAF, c("0.001", "0.01", "0.01"))

  msg <- capture.output(skatFailureSummary(out))
  expect_match(msg[1], "^SKAT failures: 3 cell\\(s\\) in 2 region\\(s\\)")
  expect_identical(msg[-1], c("  2 x boom", "  1 x other"))   # most frequent first
  side <- paste0(out, ".skatFailures.txt")
  expect_true(file.exists(side))
  back <- utils::read.table(side, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                            na.strings = "NA")
  expect_identical(back$Region, tab$Region)
  expect_identical(back$Error, tab$Error)
  expect_identical(back$Test, tab$Test)

  skatFailureReset()
  expect_identical(nrow(skatFailureTable()), 0L)
})
