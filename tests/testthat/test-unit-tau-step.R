# tauStepBinary (src/SAIGE_fitGLMM_fast.cpp, Rcpp-exported): the binary
# AI-REML's variance-component update for its one free component, applied on
# every iteration by fitglmmaiRPCG. LEDGER #121 (FIXED): the tol floor used to
# run before the step-halving loop, so the halving was unreachable and an
# overshoot from any tau0 ended the fit at tau2 = 0 (R's `if (tau[2] == 0)
# break`). The kernel now follows the reference rule, transcribed below, and
# these tests hold it there; the expectations marked #121 were inverted by the
# fix and are red at the commit before it.

tol <- 0.02

# GMMAT glmmkin.ai's rule, which this file's quantitative copy (fitglmmaiRPCG_q,
# `zeroVec`) also follows: floor before halving only when tau0 was itself below
# tol, halve while negative, floor unconditionally only after.
referenceStep <- function(tau0, Dtau, tol) {
  tau <- tau0 + Dtau
  if (tau < tol && tau0 < tol) tau <- 0
  step <- 1
  while (tau < 0) {
    step <- step / 2
    tau <- tau0 + step * Dtau
    if (tau < tol && tau0 < tol) tau <- 0
  }
  if (tau < tol) tau <- 0
  tau
}

# The kernel computes in float; 1e-6 is far above float rounding on these
# magnitudes and far below any difference the rule makes.
expect_step <- function(tau0, Dtau, expected) {
  testthat::expect_equal(tauStepBinary(tau0, Dtau, tol), expected, tolerance = 1e-6)
}

test_that("a step that stays above tol is applied as is, under both rules", {
  expect_step(0.5, -0.45, 0.05)
  expect_step(0.26, 0.3, 0.56)
  expect_equal(referenceStep(0.5, -0.45, tol), 0.05)
  expect_equal(referenceStep(0.26, 0.3, tol), 0.56)
})

test_that("#121: an overshoot from above tol is halved, not floored to 0", {
  expect_step(0.26, -0.4, 0.06) # one halving; 0 before the fix
  expect_equal(referenceStep(0.26, -0.4, tol), 0.06)
})

test_that("an overshoot whose halved step still lands below tol is 0 under both rules", {
  expect_step(0.26, -0.5, 0) # halved: 0.01, under the final floor
  expect_equal(referenceStep(0.26, -0.5, tol), 0)
})

test_that("from below tol, a negative step is 0 under both rules", {
  expect_step(0.01, -0.5, 0)
  expect_equal(referenceStep(0.01, -0.5, tol), 0)
})

test_that("#121: the kernel agrees with the reference rule on a (tau0, Dtau) grid", {
  grid <- expand.grid(tau0 = c(0.01, 0.05, 0.26, 1, 3),
                      Dtau = c(-4, -1.5, -0.4, -0.27, -0.2, 0, 0.5))
  kernel <- mapply(function(a, b) tauStepBinary(a, b, tol), grid$tau0, grid$Dtau)
  reference <- mapply(function(a, b) referenceStep(a, b, tol), grid$tau0, grid$Dtau)
  expect_equal(kernel, reference, tolerance = 1e-6)
  # the grid does exercise the branch the fix opened: overshoots from tau0 >= tol
  # whose halved step clears tol, which were 0 before the fix
  opened <- grid$tau0 >= tol & grid$tau0 + grid$Dtau < tol & reference > 0
  expect_true(any(opened))
})
