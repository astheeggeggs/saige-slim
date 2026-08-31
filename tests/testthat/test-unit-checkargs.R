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

test_that("checkArgBool accepts anything %in% c(TRUE, FALSE) coerces to", {
  # `arg0 %in% c(TRUE, FALSE)` coerces both sides, so 1/0 and the STRINGS
  # "TRUE"/"FALSE" pass the gate that is meant to reject non-logicals.
  # LEDGER #21: recorded, not fixed -- the CLI's optparse layer already
  # supplies real logicals, so this only bites direct R API callers.
  expect_silent(checkArgBool(1, "isLOCO"))
  expect_silent(checkArgBool(0, "isLOCO"))
  expect_silent(checkArgBool("TRUE", "isLOCO"))
})

test_that("checkArgBool on a vector raises R's condition-length error", {
  # `if(!arg0 %in% ...)` gets a length-2 condition; since R 4.2 that is an
  # error from the language, not the intended "should be TRUE or FALSE".
  # LEDGER #21: recorded, not fixed.
  expect_error(checkArgBool(c(TRUE, TRUE), "isLOCO"), "condition has length > 1")
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

test_that("checkArgNumeric exclusive-bound messages lose the space before the bound", {
  # "greater than" / "less than" have no trailing space, so the exclusive
  # messages render as "greater than0" / "less than1" (the inclusive forms are
  # spaced correctly). LEDGER #21: recorded, not fixed -- cosmetic, but pinned
  # so a future edit does not change user-visible text unnoticed.
  expect_error(checkArgNumeric(0, "p", minVal = 0, incMin = FALSE), "greater than0", fixed = TRUE)
  expect_error(checkArgNumeric(1, "p", maxVal = 1, incMax = FALSE), "less than1", fixed = TRUE)
})

test_that("checkArgNumeric on NA fails with R's condition error, not its own message", {
  # is.numeric(NA_real_) is TRUE, so the guard passes and `NA < minVal` yields
  # NA, which `if` rejects. The user sees "missing value where TRUE/FALSE
  # needed" with no argument name. LEDGER #21: recorded, not fixed.
  expect_error(checkArgNumeric(NA_real_, "minMAF", minVal = 0),
               "missing value where TRUE/FALSE needed")
  # with no bounds to compare against, NA passes silently
  expect_silent(checkArgNumeric(NA_real_, "minMAF"))
})

test_that("checkArgNumeric applies both bounds when given", {
  expect_silent(checkArgNumeric(0.5, "maxMAF", minVal = 0, maxVal = 1))
  expect_error(checkArgNumeric(2, "maxMAF", minVal = 0, maxVal = 1), "less than or equal to 1")
  expect_error(checkArgNumeric(-1, "maxMAF", minVal = 0, maxVal = 1), "greater than or equal to 0")
})
