# get_ACATV_pvalue (R/SAIGE_SPATest_Region_Func.R): the multi-weight ACAT-V
# cell statistic (K > 1 only; the single-weight schema has no ACAT-V column).
# Zero-weight markers are dropped (an all-zero weight vector would NaN CCT's
# normalisation and abort the run); p-values that underflowed double
# precision are clamped to the smallest normal double instead of collapsing
# the statistic to exactly 0.

test_that("get_ACATV_pvalue drops zero-weight markers", {
  expect_identical(get_ACATV_pvalue(c(0.5, 0.01), c(0, 2)),
                   CCT(0.01, 2))
  # weight-0 markers never influence the value
  expect_identical(get_ACATV_pvalue(c(0.5, 0.01, 0.9), c(0, 2, 0)),
                   get_ACATV_pvalue(0.01, 2))
})

test_that("get_ACATV_pvalue returns NA when nothing survives", {
  expect_identical(get_ACATV_pvalue(c(0.5, 0.9), c(0, 0)), NA)
  expect_identical(get_ACATV_pvalue(c(NA_real_, NA_real_), c(1, 1)), NA)
  expect_identical(get_ACATV_pvalue(numeric(0), numeric(0)), NA)
})

test_that("get_ACATV_pvalue survives p-values that underflow double precision", {
  # exp(-900) == 0 exactly; without the clamp CCT returns exactly 0
  p <- get_ACATV_pvalue(c(exp(-900), 0.5), c(1, 1))
  expect_true(is.finite(p) && p > 0 && p < 1e-300)
  # denormal p: CCT's small-p 1/p branch overflows to Inf on these unclamped
  p2 <- get_ACATV_pvalue(c(1e-315, 0.5), c(1, 1))
  expect_true(is.finite(p2) && p2 > 0)
  # NA p-values are still dropped, as in get_CCT_pvalue
  expect_identical(get_ACATV_pvalue(c(NA, 0.3), c(1, 1)),
                   get_ACATV_pvalue(0.3, 1))
})
