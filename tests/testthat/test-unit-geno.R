# Geno.R: PLINK format detection, markerInfo construction, and conditioning-
# marker resolution. setGenoInput calls setPLINKobjInCPP, which replaces a
# file-static C++ pointer (leaking the previous object) — safe to call
# repeatedly in one session, results unaffected.

bed100 <- function() file.path(saige_extdata, "input", "genotype_100markers.bed")
pgen100 <- function() file.path(saige_extdata, "input", "genotype_100markers.pgen")
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
  # LEDGER #15c fixed: missing genotype input now stops cleanly. Phase 7
  # widened this from "bedFile is required" to accept pgenFile too.
  expect_error(checkGenoInput(bedFile = "", sampleInModel = "a"),
               "is required as genotype input")
  # Phase 7: bedFile and pgenFile are mutually exclusive.
  expect_error(checkGenoInput(bedFile = bed100(), pgenFile = pgen100(),
                              sampleInModel = "a"),
               "mutually exclusive")
})

test_that("checkGenoInput accepts pgen input and defaults pvar/psam from pgenFile", {
  obj <- checkGenoInput(pgenFile = pgen100(), sampleInModel = "a")
  expect_identical(obj$genoType, "pgen")
  # E6: the resolved companions come back with the type
  expect_identical(obj$pvarFile, sub("pgen$", "pvar", pgen100()))
  expect_identical(obj$psamFile, sub("pgen$", "psam", pgen100()))
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

test_that("setGenoInput builds the pgen markerInfo contract, matching plink", {
  # PgenReader seam (Phase 7): a real, lightweight pgenlib load (not a full
  # CLI run), over the SAME fixture as the plink test above, so the two
  # markerInfo tables should agree on everything but genoType.
  plink_obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids())
  obj <- setGenoInput(pgenFile = pgen100(), sampleInModel = fam_ids())
  expect_identical(obj$genoType, "pgen")
  expect_false(obj$anyInclude)
  mi <- obj$markerInfo
  expect_identical(nrow(mi), 100L)
  expect_identical(colnames(mi), c("CHROM", "ID", "POS", "genoIndex"))
  expect_identical(data.table::key(mi), "ID")
  expect_identical(mi, plink_obj$markerInfo)

  # second call replaces the C++ reader and returns identical results
  obj2 <- setGenoInput(pgenFile = pgen100(), sampleInModel = fam_ids())
  expect_identical(obj$markerInfo, obj2$markerInfo)
})

test_that("chrom filter and include-file plumbing", {
  # LEDGER #76 FIXED: an empty chromosome selection stops in the reader now
  # (it used to hand an empty markerInfo on, and the region path wrote a
  # zero-byte results file from it).
  expect_error(setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(), chrom = "2"),
               "No markers on chrom 2 are found")

  ids <- tempfile(); writeLines(c("rs13", "rs79"), ids)
  obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                      idstoIncludeFile = ids)
  expect_true(obj$anyInclude)
  expect_identical(sort(obj$markerInfo$ID), c("rs13", "rs79"))

  # LEDGER #11 fixed: ranges now actually select markers on the PLINK path
  rng <- tempfile(); writeLines("1	10	20", rng)
  obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                      rangestoIncludeFile = rng)
  expect_true(obj$anyInclude)
  expect_identical(nrow(obj$markerInfo), 11L)  # POS 10..20 inclusive
  expect_true(all(obj$markerInfo$POS >= 10 & obj$markerInfo$POS <= 20))
})

test_that("a multi-column idstoIncludeFile is rejected and names the file", {
  ids <- tempfile(); writeLines(c("rs13\tjunk", "rs79\tjunk"), ids)
  expect_error(setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                            idstoIncludeFile = ids),
               "should only include one column", fixed = TRUE)
  expect_error(setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                            idstoIncludeFile = ids),
               ids, fixed = TRUE)
})

test_that("a malformed rangestoIncludeFile fails in fread, before its own guard", {
  # Geno.R:127 checks ncol(RangesToInclude) != 3, but fread is called with a
  # length-3 colClasses, so any other column count errors inside fread first
  # and that stop() cannot be reached. Pinned as the message the user actually
  # sees. LEDGER #24: FIXED, the pin below is inverted.
  # LEDGER #24 FIXED: the intended message is reachable now.
  rng <- tempfile(); writeLines("1\t10", rng)
  expect_error(setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                            rangestoIncludeFile = rng),
               "should only include three columns")
})

test_that("include files that match nothing leave an empty markerInfo", {
  ids <- tempfile(); writeLines(c("rsNOPE1", "rsNOPE2"), ids)
  obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                      idstoIncludeFile = ids)
  expect_true(obj$anyInclude)
  expect_identical(nrow(obj$markerInfo), 0L)

  rng <- tempfile(); writeLines("1\t100000\t200000", rng)
  obj <- setGenoInput(bedFile = bed100(), sampleInModel = fam_ids(),
                      rangestoIncludeFile = rng)
  expect_identical(nrow(obj$markerInfo), 0L)
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
