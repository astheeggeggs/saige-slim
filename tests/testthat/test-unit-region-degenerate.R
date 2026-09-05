# region_postprocess (R/SAIGE_SPATest_Region.R): the degenerate corner where
# the region has rare/ultra-rare carriers (so the core runs) but NO
# (anno x max_MAF x weight set) cell emits a row -- every flagged marker's
# p-value came back "NA", or a weighted ultra-rare collapse produced an
# all-zero pseudo-genotype.
#
# This is the branch added as the C6 review fix: `rep(NA, ncol(NULL))` used to
# error out mid-run, so the guard reproduces HEAD's legacy construction
# instead -- get_CCT_pvalue(NULL) is NA for every field and rbind(NULL, cctVec)
# yields a one-row CHARACTER MATRIX (not a data.frame). Nothing in the golden
# fixtures reaches this path, which is why it is driven directly here.
#
# The contract being pinned is "continue, do not crash, and keep the legacy
# shape". If a future change makes this return a data.frame, the writer
# downstream would silently emit different bytes for the degenerate region.

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

degenerate_run <- function(traitType = "binary", weight_sets = NULL,
                           WEIGHT = c(0), nIndCol = 2L) {
  region_postprocess(
    outList = degenerate_outList(nIndCol = nIndCol),
    mu = rep(0.3, 10), WEIGHT = WEIGHT, traitType = traitType, r.corr = 0,
    regionTestType = "SKAT-O", isCondition = FALSE,
    annolistsub = c("lof", "missense"), maxMAFlist = c(0.01),
    regionName = "GENE1", is_single_in_groupTest = TRUE, is_rerun = FALSE,
    weight_sets = weight_sets)
}

test_that("a region with no emitting cell returns a legacy Cauchy row instead of erroring", {
  r <- degenerate_run("binary")
  expect_identical(dim(r$pval.Region), c(1L, 13L))
  expect_true(is.matrix(r$pval.Region))
  expect_type(r$pval.Region, "character")
  expect_identical(unname(r$pval.Region[1, 1]), "GENE1")
  expect_identical(unname(r$pval.Region[1, 2]), "Cauchy")
  # every statistic field is NA: get_CCT_pvalue(NULL) is NA
  expect_true(all(is.na(r$pval.Region[1, 3:13])))
  expect_identical(r$cctpval, NA)
})

test_that("the degenerate row drops the two binary-only columns for a quantitative trait", {
  r <- degenerate_run("quantitative")
  expect_identical(dim(r$pval.Region), c(1L, 11L))
  expect_identical(unname(r$pval.Region[1, 2]), "Cauchy")
})

test_that("multi-weight reuses the same degenerate row, not one row per set", {
  # Exactly ONE pooled Cauchy row per region, whatever K is -- the degenerate
  # branch is shared, so K = 2 must not multiply it.
  ws <- list(kind = c("custom", "custom"), customCol = c(1, 2),
             betaA = c(NA_real_, NA_real_), betaB = c(NA_real_, NA_real_),
             force = FALSE, name = c("CADD", "REVEL"))
  r <- degenerate_run("binary", weight_sets = ws, WEIGHT = matrix(1, 4, 2), nIndCol = 4L)
  expect_identical(dim(r$pval.Region), c(1L, 13L))
  expect_identical(unname(r$pval.Region[1, 1:2]), c("GENE1", "Cauchy"))
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
