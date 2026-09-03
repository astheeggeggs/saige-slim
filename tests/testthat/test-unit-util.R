# R/Util.R: output-file/index protocol and small pure helpers.

# checkOutputFile gained a fourth element, `chunkBytes`, with the LEDGER #47 fix
# (S9): the size the output file had at the end of each recorded chunk, so a
# restart can truncate to it. These assertions were whole-list `expect_identical`
# comparisons and are now field-by-field -- pinning the three fields that existed
# plus the shape of the new one, rather than pinning "the list has exactly three
# elements", which is a property nothing depends on and which made an additive
# change look like a regression.
expect_output_state <- function(actual, Start, End, indexChunk) {
  expect_identical(actual$Start, Start)
  expect_identical(actual$End, End)
  expect_identical(actual$indexChunk, indexChunk)
}

test_that("checkOutputFile fresh/overwrite semantics", {
  f <- tempfile(); fi <- tempfile()
  out <- checkOutputFile(f, fi, "Region", 5, TRUE)
  expect_output_state(out, TRUE, FALSE, 0)
  # a fresh run has no recorded chunks at all
  expect_identical(out$chunkBytes, numeric(0))
  f2 <- tempfile(); fi2 <- tempfile()
  expect_output_state(checkOutputFile(f2, fi2, "Marker", 5, TRUE), TRUE, FALSE, 1)
})

test_that("index-file round trip: resume and completion states", {
  f <- tempfile(); fi <- tempfile(); file.create(f)
  writeOutputFileIndex(OutputFileIndex = fi, AnalysisType = "Region",
                       nEachChunk = "5", indexChunk = 3, Start = TRUE, End = FALSE)
  out <- checkOutputFile(f, fi, "Region", "5", FALSE)
  expect_output_state(out, FALSE, FALSE, 3)
  # The region path writes its index without a byte count, so a resume there
  # cannot be verified -- NA, named for the chunk it came from. Only the marker
  # path records one (LEDGER #47).
  expect_identical(out$chunkBytes, c("3" = NA_real_))

  writeOutputFileIndex(OutputFileIndex = fi, AnalysisType = "Region",
                       nEachChunk = "5", indexChunk = 3, Start = FALSE, End = TRUE)
  out <- checkOutputFile(f, fi, "Region", "5", FALSE)
  expect_true(out$End)
  expect_identical(out$indexChunk, 3)
})

test_that("checkOutputFile rejects mismatched or missing index files", {
  f <- tempfile(); fi <- tempfile(); file.create(f)
  expect_error(checkOutputFile(f, fi, "Region", "5", FALSE))  # output exists, no index

  writeOutputFileIndex(OutputFileIndex = fi, AnalysisType = "Marker",
                       nEachChunk = "5", indexChunk = 1, Start = TRUE, End = FALSE)
  expect_error(checkOutputFile(f, fi, "Region", "5", FALSE))  # wrong AnalysisType
})

test_that("updateChrStartEndIndexVec computes 0-based per-chromosome bounds", {
  out <- updateChrStartEndIndexVec(c(1, 1, 2, 2, 2, 3))
  expect_true(out$LOCO)
  expect_equal(out$chromosomeStartIndexVec[1:3], c(0, 2, 5))
  expect_equal(out$chromosomeEndIndexVec[1:3], c(1, 4, 5))
  expect_true(all(is.na(out$chromosomeStartIndexVec[4:22])))

  expect_error(updateChrStartEndIndexVec(c(2, 2, 1)), "ordered from 1 to 22")
  expect_false(updateChrStartEndIndexVec(rep(1, 5))$LOCO)  # single autosome
})

test_that("convert_str_to_log handles sub-underflow p-value strings", {
  expect_equal(convert_str_to_log("1.5E-10"), log(1.5e-10))
  expect_equal(convert_str_to_log("1E-400"), -400 * log(10))  # below double underflow
  expect_true(is.na(convert_str_to_log("NA")))
  # characterization: strings without a capital E silently return NA
  expect_true(is.na(suppressWarnings(convert_str_to_log("0.5"))))
})

test_that("small pure helpers", {
  expect_true(isContainValue(c(1, 2, 3), 2))
  expect_false(isContainValue(c(1, 2, 3), 9))
  expect_error(checkArgsBool("x", "flagname"), "flagname should be TRUE or FALSE")
  expect_error(checkArgsNumeric(0.7, "x", 0, 0.5), "0 to 0.5")
  # LEDGER #15a fixed: the message now reflects the actual bounds
  expect_error(checkArgsNumeric(5, "y", 1, 3), "ranging from 1 to 3")
})
