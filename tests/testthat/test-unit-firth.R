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
                maxhs = 15L, lconv = 1e-5, gconv = 1e-5, xconv = 1e-5,
                col_fit = as.integer(seq_len(ncol(x)) - 1L),
                weight = rep(1, nrow(x)), offset = rep(0, nrow(x))) {
  as.vector(fast_logistf_fit(x, y, weight, offset, firth, col_fit, init,
                             maxit, maxstep, maxhs, lconv, gconv, xconv, FALSE))
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

test_that("col_fit, maxhs and lconv are accepted but ignored", {
  # The upstream logistf step-halving (maxhs/lconv) and the col_fit column
  # subsetting were dropped from this port: the parameters remain in the
  # signature but never appear in the loop body (src/Main.cpp:2377-2437).
  # LEDGER #20: recorded, not fixed -- callers may reasonably expect col_fit
  # to restrict which columns are updated, and it does not.
  set.seed(7)
  x <- cbind(1, rnorm(80)); y <- rbinom(80, 1, 0.4)
  base <- flf(x, y, firth = TRUE)
  expect_identical(flf(x, y, firth = TRUE, col_fit = 0L, maxhs = 99L, lconv = 1e-1),
                   base)
})

test_that("maxit = 0 still performs one update", {
  # `while(iter <= maxit)` tests before the increment, so maxit = 0 executes
  # exactly one Newton step rather than returning `init` untouched.
  # LEDGER #20: recorded, not fixed.
  d <- fixture(n = 80, seed = 3)
  init <- c(0.3, 0.7)
  out <- flf(d$x, d$y, firth = TRUE, init = init, maxit = 0L)
  expect_false(isTRUE(all.equal(out, init)))
  expect_true(all(is.finite(out)))
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
