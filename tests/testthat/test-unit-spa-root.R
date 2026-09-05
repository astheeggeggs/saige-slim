# LEDGER #45 FIXED: the NaN guard in both saddlepoint root finders was written
# `tnew == NA_REAL`, which is false for every operand because NA_REAL is a NaN;
# it is `std::isnan(tnew)` now. The provenance shows the branch is unreachable
# on the live inputs (K2 vanishes only for an all-zero g, which the bracket test
# above it catches first), so this pins the behaviour around it rather than the
# branch: an ordinary root is found and reported converged, and the all-zero
# corner returns its infinite root and, since LEDGER #22 was fixed, reports
# non-convergence. Fast tier.

test_that("getroot_K1_Binom finds an ordinary root and reports convergence", {
  set.seed(45L)
  n <- 200
  mu <- stats::runif(n, 0.05, 0.4)
  g <- sample(0:2, n, replace = TRUE, prob = c(0.7, 0.25, 0.05))
  q <- sum(g * mu) + 3            # a modest positive deviation from the mean
  r <- getroot_K1_Binom(0, mu, g, q, 1e-8, 100L)
  expect_true(isTRUE(r$Isconverge))
  expect_true(is.finite(r$root))
  # K1(root) == q, the root's defining identity, recomputed here
  K1 <- function(t) sum(mu * g * exp(t * g) / (1 - mu + mu * exp(t * g)))
  expect_lt(abs(K1(r$root) - q), 1e-6)
})

test_that("an all-zero genotype vector takes the bracket exit and reports non-convergence (LEDGER #22 FIXED)", {
  mu <- rep(0.2, 50)
  r <- getroot_K1_Binom(0, mu, rep(0, 50), 1, 1e-8, 100L)
  expect_true(is.infinite(r$root))
  expect_false(isTRUE(r$Isconverge))
})
