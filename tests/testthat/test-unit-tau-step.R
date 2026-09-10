# tauStepBinary (src/SAIGE_fitGLMM_fast.cpp, Rcpp-exported): the binary
# AI-REML's variance-component update for its one free component, applied on
# every iteration by fitglmmaiRPCG. LEDGER #121: the tol floor runs before the
# step-halving loop, so the halving is unreachable and an overshoot from any
# tau0 ends the fit at tau2 = 0 (R's `if (tau[2] == 0) break`). These tests
# pin what the kernel does TODAY beside an R transcription of the reference
# rule; the expectations marked #121 invert when the row is fixed.

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

test_that("#121: an overshoot from above tol is floored to 0 before the halving can run", {
  expect_step(0.26, -0.4, 0) # one halving would give 0.06
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

test_that("#121: the two rules disagree exactly on overshoots whose halved step clears tol", {
  grid <- expand.grid(tau0 = c(0.01, 0.05, 0.26, 1, 3),
                      Dtau = c(-4, -1.5, -0.4, -0.27, -0.2, 0, 0.5))
  kernel <- mapply(function(a, b) tauStepBinary(a, b, tol), grid$tau0, grid$Dtau)
  reference <- mapply(function(a, b) referenceStep(a, b, tol), grid$tau0, grid$Dtau)
  differs <- abs(kernel - reference) > 1e-6
  expect_true(any(differs))
  expect_true(all(kernel[differs] == 0))
  expect_true(all(reference[differs] > 0))
  expect_true(all(grid$tau0[differs] >= tol))
  expect_true(all(grid$tau0[differs] + grid$Dtau[differs] < tol))
})
