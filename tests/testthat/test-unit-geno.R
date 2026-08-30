# Geno.R: PLINK format detection, markerInfo construction, and conditioning-
# marker resolution. setGenoInput calls setPLINKobjInCPP, which replaces a
# file-static C++ pointer (leaking the previous object) — safe to call
# repeatedly in one session, results unaffected.

bed100 <- function() file.path(saige_extdata, "input", "genotype_100markers.bed")
fam_ids <- function() {
  fam <- utils::read.table(file.path(saige_extdata, "input", "genotype_100markers.fam"),
                           colClasses = "character")
  fam$V2
}

test_that("checkGenoInput validation", {
  expect_error(checkGenoInput(bedFile = bed100(), sampleInModel = NULL),
               "sampleInModel is not specified")
  expect_error(checkGenoInput(bedFile = "/nope.bed", sampleInModel = "a"),
               "does not exsit")  # typo verbatim
  # characterization: no bedFile at all fails messily, not with a clean stop
  expect_error(checkGenoInput(bedFile = "", sampleInModel = "a"),
               "'dosageFileType' not found")
})

test_that("setGenoInput builds the plink markerInfo contract", {
  obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids())
  expect_identical(obj$genoType, "plink")
  expect_false(obj$anyInclude)
  mi <- obj$markerInfo
  expect_identical(nrow(mi), 100L)
  expect_identical(colnames(mi), c("CHROM", "ID", "POS", "genoIndex"))
  expect_identical(data.table::key(mi), "ID")
  # C-locale sort by ID; genoIndex is the 0-based bim row
  expect_identical(mi$ID[1:3], c("rs1", "rs10", "rs100"))
  expect_identical(mi$genoIndex[1:3], c(0, 9, 99))

  # second call replaces the C++ reader and returns identical results
  obj2 <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids())
  expect_identical(obj$markerInfo, obj2$markerInfo)
})

test_that("chrom filter and include-file plumbing", {
  obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(), chrom = "2")
  expect_identical(nrow(obj$markerInfo), 0L)  # all markers are on chrom 1

  ids <- tempfile(); writeLines(c("rs13", "rs79"), ids)
  obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                      idstoIncludeFile = ids)
  expect_true(obj$anyInclude)
  expect_identical(sort(obj$markerInfo$ID), c("rs13", "rs79"))
})

test_that("extract_genoIndex_condition resolves rsIDs in condition order", {
  mi <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids())$markerInfo
  out <- extract_genoIndex_condition("rs13,rs79", mi, "plink")
  expect_equal(out$cond_genoIndex, c(12, 78))
  # order follows the condition string, not genoIndex
  expect_equal(extract_genoIndex_condition("rs79,rs13", mi, "plink")$cond_genoIndex,
               c(78, 12))
  expect_error(extract_genoIndex_condition("rs13,rs999", mi, "plink"),
               "1 conditioning markers are not found in the geno file")
  expect_error(extract_genoIndex_condition("", mi, "plink"), "condition is empty!")
})
