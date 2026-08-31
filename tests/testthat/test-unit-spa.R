# Saddlepoint kernels for the binary trait (src/SPA_binary.cpp, all
# Rcpp-exported): the binomial cumulant generating function K, its adjusted
# first derivative K1 - q, its second derivative K2, the Newton root finder
# over K1, the saddlepoint tail probability, and the log-space adder.
#
# These are pure numeric functions, so they are pinned against their closed
# forms rather than against stored output: at t = 0 the binomial CGF has exact
# analytic values, and K1/K2 must agree with numerical derivatives of K. That
# is the reconciliation that makes the SPA cascade trustworthy; the golden
# tests only ever see the composed result.

spa_fixture <- function(n = 200, seed = 5) {
  set.seed(seed)
  list(mu = runif(n, 0.05, 0.95), g = rnorm(n))
}

test_that("the CGF and its derivatives match their closed forms at t = 0", {
  d <- spa_fixture()
  mu <- d$mu; g <- d$g
  # K(0) = sum(log(1 - mu + mu)) = 0
  expect_equal(Korg_Binom(0, mu, g), 0)
  # K'(0) = sum(mu * g); the exported form subtracts q
  expect_equal(K1_adj_Binom(0, mu, g, 0), sum(mu * g))
  expect_equal(K1_adj_Binom(0, mu, g, 2.5), sum(mu * g) - 2.5)
  # K''(0) = sum(mu (1 - mu) g^2), the score variance
  expect_equal(K2_Binom(0, mu, g), sum(mu * (1 - mu) * g^2))
})

test_that("K1 and K2 are the numerical derivatives of K", {
  d <- spa_fixture()
  mu <- d$mu; g <- d$g
  h <- 1e-6
  for (t0 in c(-0.2, 0.05, 0.3)) {
    expect_equal((Korg_Binom(t0 + h, mu, g) - Korg_Binom(t0 - h, mu, g)) / (2 * h),
                 K1_adj_Binom(t0, mu, g, 0), tolerance = 1e-6)
    expect_equal((K1_adj_Binom(t0 + h, mu, g, 0) - K1_adj_Binom(t0 - h, mu, g, 0)) / (2 * h),
                 K2_Binom(t0, mu, g), tolerance = 1e-6)
  }
})

test_that("K2 is non-negative (it is a variance) and zero for a null genotype", {
  d <- spa_fixture()
  expect_true(all(vapply(c(-1, -0.1, 0, 0.1, 1),
                         function(t) K2_Binom(t, d$mu, d$g), 0) >= 0))
  expect_identical(K2_Binom(0.5, d$mu, rep(0, length(d$g))), 0)
})

test_that("getroot_K1_Binom solves K1(root) = q", {
  d <- spa_fixture()
  m <- sum(d$mu * d$g)
  for (q in c(m + 3, m + 8, m - 4)) {
    r <- getroot_K1_Binom(0, d$mu, d$g, q, 1e-10, 100L)
    expect_true(r$Isconverge)
    expect_equal(K1_adj_Binom(r$root, d$mu, d$g, q), 0, tolerance = 1e-8)
  }
})

test_that("getroot_K1_Binom reports non-convergence when the iteration budget runs out", {
  d <- spa_fixture()
  r <- getroot_K1_Binom(0, d$mu, d$g, sum(d$mu * d$g) + 8, 1e-12, 1L)
  expect_false(r$Isconverge)
  expect_identical(r$niter, 1L)
})

test_that("getroot_K1_Binom reports an infinite root as converged for an all-zero genotype", {
  # K1 is identically -q, so Newton never moves and the routine walks the root
  # out to Inf while still flagging convergence. Unreachable through the CLI
  # (a monomorphic marker is dropped by minMAC before it gets here), but pinned
  # so the degenerate branch is not silently changed. LEDGER #22: recorded,
  # not fixed.
  d <- spa_fixture()
  r <- getroot_K1_Binom(0, d$mu, rep(0, length(d$g)), 0, 1e-8, 100L)
  expect_true(is.infinite(r$root))
  expect_true(r$Isconverge)
})

test_that("Get_Saddle_Prob_Binom returns a SIGNED half-tail, not a probability", {
  # The lower-tail branch returns -pnorm(...) deliberately; callers take
  # abs() (SPA_binary's linear path) or feed add_logp, which applies -abs()
  # itself. Do not "fix" a negative return here into a probability.
  d <- spa_fixture()
  m <- sum(d$mu * d$g)
  q_lo <- m - 6
  r <- getroot_K1_Binom(0, d$mu, d$g, q_lo, 1e-10, 100L)
  s <- Get_Saddle_Prob_Binom(r$root, d$mu, d$g, q_lo, FALSE)
  expect_true(s$isSaddle)
  expect_lt(s$pval, 0)
  expect_lt(abs(s$pval), 1)
})

test_that("Get_Saddle_Prob_Binom degenerates at zeta = 0", {
  # w == 0 fails the flagrun guard, so no saddlepoint is computed and the
  # caller falls back to pval_noadj/2.
  d <- spa_fixture()
  m <- sum(d$mu * d$g)
  lin <- Get_Saddle_Prob_Binom(0, d$mu, d$g, m, FALSE)
  expect_false(lin$isSaddle)
  expect_identical(lin$pval, 0)
  lg <- Get_Saddle_Prob_Binom(0, d$mu, d$g, m, TRUE)
  expect_false(lg$isSaddle)
  expect_identical(lg$pval, -Inf)
})

test_that("SPA_binary log and linear paths reconcile", {
  # Two ways to compute one quantity: abs(p1)+abs(p2) versus add_logp(p1,p2).
  # They agree only because add_logp applies -abs() to both arguments, which
  # is what makes the signed half-tails above safe.
  d <- spa_fixture()
  m <- sum(d$mu * d$g)
  for (q in c(m + 3, m + 6, m + 10)) {
    qinv <- 2 * m - q
    lin <- SPA_binary(d$mu, d$g, q, qinv, 0.05, 1e-8, FALSE)
    lg  <- SPA_binary(d$mu, d$g, q, qinv, log(0.05), 1e-8, TRUE)
    expect_true(lin$Isconverge && lg$Isconverge)
    expect_equal(exp(lg$pval), lin$pval, tolerance = 1e-10)
    expect_gt(lin$pval, 0)
    expect_lte(lin$pval, 1)
  }
})

test_that("SPA_binary returns pval_noadj unchanged for an unreachable q", {
  # A q beyond anything K1 can attain sends the root finder straight to Inf,
  # which it reports as CONVERGED in 0 iterations (as in the all-zero case
  # above). Get_Saddle_Prob_Binom then sees non-finite k1/k2, declines the
  # saddlepoint, and each tail falls back to pval_noadj/2 -- so the two halves
  # reassemble the unadjusted p-value exactly, with Isconverge still TRUE.
  # The Isconverge = FALSE branch needs the root finder to exhaust its own
  # (defaulted, not caller-supplied) iteration budget and is not reachable
  # through this export.
  d <- spa_fixture()
  m <- sum(d$mu * d$g)
  huge <- sum(abs(d$g)) + 1e3
  expect_true(is.infinite(getroot_K1_Binom(0, d$mu, d$g, huge, 1e-8, 100L)$root))
  r <- SPA_binary(d$mu, d$g, huge, 2 * m - huge, 0.0123, 1e-8, FALSE)
  expect_true(r$Isconverge)
  expect_identical(r$pval, 0.0123)
})

test_that("add_logp is log-sum-exp, symmetric, and insensitive to argument sign", {
  expect_equal(add_logp(log(0.1), log(0.2)), log(0.3))
  expect_equal(add_logp(log(0.2), log(0.1)), log(0.3))
  # -abs() on both arguments: a positive (sign-flipped) input is treated as
  # its own negation, which is what makes the SPA log path agree with linear
  expect_identical(add_logp(log(0.1), -log(0.2)), add_logp(log(0.1), log(0.2)))
  # far below the underflow point of exp(), where naive log(p1 + p2) gives -Inf
  expect_equal(add_logp(-700, -700), log(2) - 700)
  expect_true(is.finite(add_logp(-5000, -5000)))
  expect_equal(add_logp(-5000, -Inf), -5000)
})

test_that("sum_arma1 sums a vector", {
  expect_identical(sum_arma1(as.numeric(1:5)), 15)
  expect_identical(sum_arma1(numeric(0)), 0)
})
