# region_postprocess (R/SAIGE_SPATest_Region.R): the degenerate corner where
# the region has rare/ultra-rare carriers (so the core runs) but NO
# (anno x max_MAF x weight set) cell emits a row -- every flagged marker's
# p-value came back "NA", or a weighted ultra-rare collapse produced an
# all-zero pseudo-genotype.
#
# LEDGER #128 (FIXED 2026-09-13) CHANGED THE CONTRACT HERE, and these tests are
# inverted to match. The branch used to reproduce HEAD's legacy construction --
# get_CCT_pvalue(NULL) is NA for every field and rbind(NULL, cctVec) yields a
# one-row CHARACTER MATRIX with NO column names. That was bad two ways, both
# MEASURED on the bundled fixture before the fix landed:
#   * with two or more regions it KILLED the run one region later, in
#     regionFinish's cross-region rbind ("names do not match previous names");
#   * with ONE region it exited 0 and wrote an 89-byte results file whose header
#     was `V1 V2 V3 ... V13` -- R's placeholder names for an unnamed matrix --
#     over a single all-NA row. A file that looks parseable and has the wrong
#     name on every column is worse than no file, which is why "keep the legacy
#     shape" was not a contract worth holding.
# LEDGER #116 (FIXED later on 2026-09-13) then DECIDED the reporting question
# and these tests were inverted again: the region contributes exactly ONE row,
# "GENE1 Cauchy NA ...", built from the schema stated once in R/RegionSchema.R
# (the same row the BURDEN writer always emitted for such a region), with the
# CCT values still computed, since cctpval is the fastTest escalation signal, and
# the #127 .pooledTests.txt sidecar carrying the region with Rows = 0.
#
# Nothing in the golden fixtures reaches this path, which is why it is driven
# directly here.

# C++ writes p-values as strings and missing ones as the literal "NA";
# convert_str_to_log() is what turns them into log(p), and it requires that
# string form (R's NA_character_ makes its `a != "NA"` test error).
degenerate_outList <- function(q = 4L, nIndCol = 2L) {
  list(
    pvalVec = rep("NA", q),
    annoMAFIndicatorMat = matrix(1L, nrow = q, ncol = nIndCol),
    q_maf_for_annoVec = c(1, 1),
    NumRare_GroupVec = c(2, 2),
    NumUltraRare_GroupVec = c(1, 1),
    TstatVec_flip = rep(0.5, q),
    # every marker dropped by the NA filter, so the variance block is empty
    VarMat = matrix(numeric(0), 0, 0),
    MAFVec = rep(0.001, q),
    MAC_GroupVec = c(5, 5),
    MACCase_GroupVec = c(3, 3),
    MACCtrl_GroupVec = c(2, 2),
    numofUR = 0,
    gyVec = rep(0.1, q)
  )
}

# the default is the unweighted legacy run: one implicit beta set, beta(1, 25)
degenerate_run <- function(traitType = "binary",
                           weight_sets = weightSetDescriptor("beta", "Beta", matrix(c(1, 25), 1)),
                           WEIGHT = c(0), nIndCol = 2L) {
  region_postprocess(
    outList = degenerate_outList(nIndCol = nIndCol),
    mu = rep(0.3, 10), WEIGHT = WEIGHT, traitType = traitType, r.corr = 0,
    regionTestType = "SKAT-O", isCondition = FALSE,
    annolistsub = c("lof", "missense"), maxMAFlist = c(0.01),
    regionName = "GENE1", is_single_in_groupTest = TRUE, is_rerun = FALSE,
    weight_sets = weight_sets)
}

expect_empty_row <- function(df, traitType, isMultiWeight = FALSE) {
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), 1L)
  expect_identical(colnames(df), regionResultColnames(traitType, FALSE, isMultiWeight))
  expect_identical(df$Region, "GENE1")
  expect_identical(df$Group, "Cauchy")
  expect_true(all(is.na(unlist(df[, setdiff(colnames(df), c("Region", "Group"))]))))
}

test_that("a region with no emitting cell builds ONE NA Cauchy row, named, and does not error", {
  r <- degenerate_run("binary")
  # LEDGER #116: one named row from the schema, so the cross-region rbind in
  # regionFinish matches names and the writer sees the real header, never V1..V13
  expect_empty_row(r$pval.Region, "binary")
  # the CCT value is still computed and returned -- it is the fastTest
  # escalation signal, and get_CCT_pvalue(NULL) is NA
  expect_identical(r$cctpval, NA)
})

test_that("the trait type selects the NA row's schema: no MAC_case/MAC_control on quantitative", {
  # 13 columns on a binary trait, 11 on a quantitative one -- the legacy
  # vector's accident is now the schema's rule, stated in regionResultColnames
  r <- degenerate_run("quantitative")
  expect_empty_row(r$pval.Region, "quantitative")
  expect_false(any(c("MAC_case", "MAC_control") %in% colnames(r$pval.Region)))
  expect_identical(ncol(r$pval.Region), 11L)
  expect_identical(ncol(degenerate_run("binary")$pval.Region), 13L)
  expect_identical(r$cctpval, NA)
})

test_that("multi-weight builds exactly one NA row, whatever K is", {
  # It used to build exactly ONE shared legacy row however large K was; K must
  # not multiply the NA row into existence either.
  ws <- list(kind = c("custom", "custom"), customCol = c(1, 2),
             betaA = c(NA_real_, NA_real_), betaB = c(NA_real_, NA_real_),
             force = FALSE, name = c("CADD", "REVEL"))
  r <- degenerate_run("binary", weight_sets = ws, WEIGHT = matrix(1, 4, 2), nIndCol = 4L)
  # exactly ONE row under the multi-weight schema (Weight, ACAT-V, ACAT-O
  # columns present), not one per set
  expect_empty_row(r$pval.Region, "binary", isMultiWeight = TRUE)
  expect_true(all(c("Weight", "Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(r$pval.Region)))
  expect_identical(r$cctpval, NA)
})

test_that("convert_str_to_log parses the C++ 'mEe' form", {
  # Pinning the coupling the fixture above depends on. The helper splits on a
  # literal "E" and recombines mantissa and exponent in log space, which is
  # what lets C++ hand over p-values far below double precision as text.
  expect_equal(convert_str_to_log("5.0E-01"), log(0.5))
  expect_equal(convert_str_to_log("1.0E-300"), log(1e-300))
  # below the smallest normal double: the whole point of the string hand-off
  expect_equal(convert_str_to_log("1.0E-400"), -400 * log(10))
  expect_identical(convert_str_to_log("NA"), NA)
})

test_that("convert_str_to_log degrades silently on input it was not written for", {
  # A plain decimal has no "E", so the exponent is NA and the result is NA
  # with no warning -- a p-value that would then read as missing. Reachable
  # only if a caller builds pvalVec itself; C++ always writes the E form.
  # LEDGER #23 FIXED: a plain decimal is log()'d, R's NA is NA, the E form is
  # unchanged (kept on the log scale, which is the point of the function).
  expect_identical(convert_str_to_log("0.5"), log(0.5))
  expect_identical(convert_str_to_log(NA_character_), NA)
  expect_identical(convert_str_to_log("NA"), NA)
  expect_equal(convert_str_to_log("2.5E-3"), log(2.5e-3), tolerance = 1e-12)
  expect_equal(convert_str_to_log("1.0E-400"), -400 * log(10), tolerance = 1e-9)
})
