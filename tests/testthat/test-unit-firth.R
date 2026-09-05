# fast_logistf_fit (src/Main.cpp, Rcpp-exported): the penalised-likelihood
# logistic fit. This is the R-callable free function, NOT the one the marker
# cascade uses -- getMarkerPval calls SAIGEClass::fast_logistf_fit_simple
# (src/SAIGE_test.cpp), a separate implementation reachable only through a
# step-2 run with --is_Firth_beta=TRUE (see test-step2-firth.R). Nothing in
# the package calls this export; it is public API via exportPattern, and
# these tests pin what it does today.
#
# The reconciliation that makes the rest meaningful: with firth = FALSE the
# routine is plain Newton-Raphson and must agree with stats::glm's MLE. With
# firth = TRUE the Jeffreys penalty pulls the estimate toward zero.

# Named for readability at the call sites; the last argument (isfirthconverge)
# is bool& in C++, so its mutation is NOT observable from R -- Rcpp materialises
# a temporary from the SEXP and discards it. Do not assert on it here.
flf <- function(x, y, firth, init = c(0, 0), maxit = 50L, maxstep = 15L,
                gconv = 1e-5, xconv = 1e-5,
                weight = rep(1, nrow(x)), offset = rep(0, nrow(x))) {
  as.vector(fast_logistf_fit(x, y, weight, offset, firth, init,
                             maxit, maxstep, gconv, xconv, FALSE))
}

# Deterministic well-separated-but-not-separable design.
fixture <- function(n = 300, seed = 42) {
  set.seed(seed)
  x2 <- rnorm(n)
  list(x = cbind(1, x2), y = rbinom(n, 1, plogis(-0.5 + 0.8 * x2)))
}

test_that("firth = FALSE reproduces the glm MLE", {
  d <- fixture()
  mle <- unname(coef(glm(d$y ~ d$x[, 2], family = binomial)))
  expect_equal(flf(d$x, d$y, firth = FALSE), mle, tolerance = 1e-10)
})

test_that("firth = TRUE shrinks the estimate toward zero", {
  d <- fixture()
  mle <- flf(d$x, d$y, firth = FALSE)
  fth <- flf(d$x, d$y, firth = TRUE)
  expect_true(all(is.finite(fth)))
  expect_lt(abs(fth[2]), abs(mle[2]))
  # the penalty perturbs, it does not overwhelm: same sign, same magnitude
  expect_identical(sign(fth[2]), sign(mle[2]))
  expect_lt(abs(fth[2] - mle[2]), 0.1 * abs(mle[2]))
})

test_that("firth = TRUE stays finite under complete separation", {
  # y is perfectly predicted by x2, so the true MLE is infinite; the penalised
  # fit is the reason Firth exists. (glm reports a large finite slope only
  # because it stops at its own iteration limit.)
  set.seed(19)
  n <- 100
  x2 <- c(rnorm(n / 2, -2), rnorm(n / 2, 2))
  y <- c(rep(0, n / 2), rep(1, n / 2))
  x <- cbind(1, x2)
  fth <- flf(x, y, firth = TRUE)
  unpen <- flf(x, y, firth = FALSE)
  expect_true(all(is.finite(fth)))
  expect_lt(abs(fth[2]), abs(unpen[2]))
})

test_that("LEDGER #20 FIXED: the dead logistf parameters are gone from the signature", {
  # The upstream logistf step-halving (maxhs/lconv) and the col_fit column
  # subsetting were dropped from this port, but the parameters stayed in the
  # signature and were silently ignored -- a caller copying the logistf
  # signature got no column restriction and no step-halving. Owner decision
  # 2026-09-05: drop them rather than implement them.
  f <- names(formals(fast_logistf_fit))
  expect_false(any(c("col_fit", "maxhs", "lconv") %in% f))
  expect_identical(f, c("x", "y", "weight", "offset", "firth", "init", "maxit",
                        "maxstep", "gconv", "xconv", "isfirthconverge"))
})

test_that("LEDGER #20 FIXED: maxit = 0 returns init untouched", {
  # `while(iter <= maxit)` used to test before the increment, so maxit = 0
  # executed one Newton step. The bound is now `iter < maxit`; for maxit >= 1
  # the iteration count is unchanged (the post-increment `iter == maxit` break
  # fires at the same step), which the glm reconciliation above keeps honest.
  d <- fixture(n = 80, seed = 3)
  init <- c(0.3, 0.7)
  expect_identical(flf(d$x, d$y, firth = TRUE, init = init, maxit = 0L), init)
  # one iteration really is one Newton step away from init
  one <- flf(d$x, d$y, firth = TRUE, init = init, maxit = 1L)
  expect_false(isTRUE(all.equal(one, init)))
  expect_true(all(is.finite(one)))
})

test_that("offset and weight enter the fit", {
  d <- fixture(n = 150, seed = 11)
  base <- flf(d$x, d$y, firth = TRUE)
  expect_false(isTRUE(all.equal(
    flf(d$x, d$y, firth = TRUE, offset = rep(1, nrow(d$x))), base)))
  expect_false(isTRUE(all.equal(
    flf(d$x, d$y, firth = TRUE, weight = rep(2, nrow(d$x))), base)))
  # a zero-weight sample contributes nothing: equals the fit without that row
  w <- rep(1, nrow(d$x)); w[1] <- 0
  expect_equal(flf(d$x, d$y, firth = TRUE, weight = w),
               flf(d$x[-1, ], d$y[-1], firth = TRUE), tolerance = 1e-8)
})
