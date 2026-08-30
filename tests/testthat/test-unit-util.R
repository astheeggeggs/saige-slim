# R/Util.R: output-file/index protocol and small pure helpers.

test_that("checkOutputFile fresh/overwrite semantics", {
  f <- tempfile(); fi <- tempfile()
  expect_identical(checkOutputFile(f, fi, "Region", 5, TRUE),
                   list(Start = TRUE, End = FALSE, indexChunk = 0))
  f2 <- tempfile(); fi2 <- tempfile()
  expect_identical(checkOutputFile(f2, fi2, "Marker", 5, TRUE),
                   list(Start = TRUE, End = FALSE, indexChunk = 1))
})

test_that("index-file round trip: resume and completion states", {
  f <- tempfile(); fi <- tempfile(); file.create(f)
  writeOutputFileIndex(OutputFileIndex = fi, AnalysisType = "Region",
                       nEachChunk = "5", indexChunk = 3, Start = TRUE, End = FALSE)
  out <- checkOutputFile(f, fi, "Region", "5", FALSE)
  expect_identical(out, list(Start = FALSE, End = FALSE, indexChunk = 3))

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
  # characterization: message hardcodes "0 to 0.5" regardless of the bounds given
  expect_error(checkArgsNumeric(0.7, "x", 0, 0.5), "0 to 0.5")
})
