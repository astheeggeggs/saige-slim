# LEDGER #85 FIXED: Covariate_Transform no longer declares a `traitType` formal
# it never read. extractNglmm called it without one and lived on the promise
# never being forced; a future read of the argument would have re-armed
# 'argument "traitType" is missing'. The signature is (formula, data) now, and
# the two callers pass exactly that. Fast tier.

test_that("Covariate_Transform takes (formula, data) and transforms a small design", {
  expect_identical(names(formals(SAIGE:::Covariate_Transform)), c("formula", "data"))
  set.seed(85L)
  d <- data.frame(y = rnorm(40), x1 = rnorm(40), x2 = rnorm(40))
  out <- SAIGE:::Covariate_Transform(y ~ x1 + x2, data = d)
  expect_true(is.list(out))
  expect_identical(nrow(out$X1), 40L)
  expect_identical(ncol(out$X1), 3L)
})
