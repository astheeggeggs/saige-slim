# checkArgBool / checkArgNumeric (R/checkArgs.R): the generic CLI argument
# validators, called ~22 times across SPAGMMATtest and the region setup. Their
# error text reaches the user directly, so the messages are part of the
# contract and are asserted here alongside the boundary semantics.
#
# NOTE on the near-duplicate family: R/Util.R also defines checkArgsBool /
# checkArgsNumeric (PLURAL "Args"). Those have no callers anywhere and are not
# what the CLI uses -- do not "consolidate" the two without checking which is
# live.

test_that("checkArgBool accepts logicals and rejects non-logicals", {
  expect_silent(checkArgBool(TRUE, "isLOCO"))
  expect_silent(checkArgBool(FALSE, "isLOCO"))
  expect_error(checkArgBool(NA, "isLOCO"), "isLOCO should be TRUE or FALSE.", fixed = TRUE)
  expect_error(checkArgBool("yes", "isLOCO"), "isLOCO should be TRUE or FALSE.", fixed = TRUE)
})

test_that("checkArgBool rejects what %in% used to coerce (LEDGER #21 FIXED, a)", {
  expect_error(checkArgBool(1, "isLOCO"), "isLOCO should be TRUE or FALSE.", fixed = TRUE)
  expect_error(checkArgBool(0, "isLOCO"), "isLOCO should be TRUE or FALSE.", fixed = TRUE)
  expect_error(checkArgBool("TRUE", "isLOCO"), "isLOCO should be TRUE or FALSE.", fixed = TRUE)
  expect_error(checkArgBool(NA, "isLOCO"), "isLOCO should be TRUE or FALSE.", fixed = TRUE)
  expect_silent(checkArgBool(TRUE, "isLOCO"))
})

test_that("checkArgBool on a vector gives the intended message (LEDGER #21 FIXED, b)", {
  expect_error(checkArgBool(c(TRUE, TRUE), "isLOCO"), "isLOCO should be TRUE or FALSE.", fixed = TRUE)
})

test_that("checkArgNumeric rejects non-numerics by name", {
  expect_error(checkArgNumeric("x", "minMAF"), "minMAF should be a numeric value.", fixed = TRUE)
  expect_error(checkArgNumeric(TRUE, "minMAF"), "minMAF should be a numeric value.", fixed = TRUE)
  expect_silent(checkArgNumeric(0.5, "minMAF"))
})

test_that("checkArgNumeric minVal boundary honours incMin", {
  # inclusive (default): the endpoint passes; exclusive: it is rejected
  expect_silent(checkArgNumeric(0, "minMAF", minVal = 0))
  expect_silent(checkArgNumeric(0, "minMAF", minVal = 0, incMin = TRUE))
  expect_error(checkArgNumeric(0, "minMAF", minVal = 0, incMin = FALSE),
               "minMAF should be a numeric value greater than")
  expect_error(checkArgNumeric(-1e-12, "minMAF", minVal = 0),
               "minMAF should be a numeric value greater than or equal to 0")
  expect_silent(checkArgNumeric(1e-12, "minMAF", minVal = 0, incMin = FALSE))
})

test_that("checkArgNumeric maxVal boundary honours incMax", {
  expect_silent(checkArgNumeric(1, "maxMAF", maxVal = 1))
  expect_error(checkArgNumeric(1, "maxMAF", maxVal = 1, incMax = FALSE),
               "maxMAF should be a numeric value less than")
  expect_error(checkArgNumeric(1 + 1e-12, "maxMAF", maxVal = 1),
               "maxMAF should be a numeric value less than or equal to 1")
  expect_silent(checkArgNumeric(1 - 1e-12, "maxMAF", maxVal = 1, incMax = FALSE))
})

test_that("checkArgNumeric exclusive-bound messages are spaced (LEDGER #21 FIXED, c)", {
  expect_error(checkArgNumeric(0, "p", minVal = 0, incMin = FALSE), "greater than 0", fixed = TRUE)
  expect_error(checkArgNumeric(1, "p", maxVal = 1, incMax = FALSE), "less than 1", fixed = TRUE)
})

test_that("checkArgNumeric on NA or a vector names the argument (LEDGER #21 FIXED, d)", {
  expect_error(checkArgNumeric(NA_real_, "minMAF", minVal = 0),
               "minMAF should be a single non-missing numeric value.", fixed = TRUE)
  expect_error(checkArgNumeric(NA_real_, "minMAF"),
               "minMAF should be a single non-missing numeric value.", fixed = TRUE)
  expect_error(checkArgNumeric(c(1, 2), "minMAF"),
               "minMAF should be a single non-missing numeric value.", fixed = TRUE)
})

test_that("checkArgNumeric applies both bounds when given", {
  expect_silent(checkArgNumeric(0.5, "maxMAF", minVal = 0, maxVal = 1))
  expect_error(checkArgNumeric(2, "maxMAF", minVal = 0, maxVal = 1), "less than or equal to 1")
  expect_error(checkArgNumeric(-1, "maxMAF", minVal = 0, maxVal = 1), "greater than or equal to 0")
})
