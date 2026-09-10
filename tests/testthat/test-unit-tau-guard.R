# tauGuardApply (R/TauStepGuard.R) and tauProposalBinary (src/SAIGE_fitGLMM_fast.cpp,
# Rcpp-exported): LEDGER #122's guard around the binary AI-REML step. Pure rules,
# so they are pinned on synthetic score functions here; the fits themselves are
# pinned in test-step1-guard.R.

tol <- 0.02

test_that("tauProposalBinary reproduces the float Newton step and the #121 kernel", {
  p <- tauProposalBinary(0.26, 154.5, 154.1, 47, tol)
  expect_equal(p[["raw"]], 0.26 + 0.4 / 47, tolerance = 1e-6)
  expect_equal(p[["step"]], tauStepBinary(0.26, 0.4 / 47, tol), tolerance = 1e-6)
  q <- tauProposalBinary(0.26, 100, 120, 50, tol) # an overshoot to -0.14: halved by #121
  expect_equal(q[["raw"]], 0.26 - 0.4, tolerance = 1e-6)
  expect_equal(q[["step"]], 0.06, tolerance = 1e-6)
})

# Drive the guard along a sequence of (tau0, raw) pairs, the step being the raw
# proposal itself (no floor or halving needed on these values).
drive <- function(taus, raws) {
  st <- SAIGE::tauGuardInit()
  out <- numeric(0)
  for (k in seq_along(taus)) {
    sg <- SAIGE::tauGuardApply(st, taus[k], raws[k], raws[k])
    st <- sg$state
    out <- c(out, sg$p)
  }
  list(p = out, state = st)
}

test_that("the bundled binary fit's trace passes untouched: one reversal does not arm", {
  # golden log step1_binary_sparseGRM_a: 0 -> 0.336 -> 0.373 -> 0.336 -> 0.327
  taus <- c(0, 0.3357163, 0.372536, 0.3362476)
  raws <- c(0.3357163, 0.372536, 0.3362476, 0.3274886)
  d <- drive(taus, raws)
  expect_identical(d$p, raws)
  expect_equal(d$state$reversals, 1L)
  expect_false(d$state$armed)
  expect_false(d$state$fired)
})

test_that("the second reversal arms the guard, later reversals halve the step", {
  # up, down, up (armed, full step still), down (halved), up (quartered)
  taus <- c(0.5, 1.0, 0.6, 0.9, 0.7)
  raws <- c(1.0, 0.6, 0.9, 0.7, 0.8)
  d <- drive(taus, raws)
  expect_equal(d$p[1:3], raws[1:3])          # nothing applied before the second reversal
  expect_equal(d$p[4], 0.9 + 0.5 * (0.7 - 0.9)) # third reversal: half the step
  expect_equal(d$p[5], 0.7 + 0.25 * (0.8 - 0.7))   # fourth reversal: a quarter
  expect_true(d$state$armed)
  expect_true(d$state$fired)
})

test_that("the growth cap applies only when armed", {
  st <- SAIGE::tauGuardInit()
  sg <- SAIGE::tauGuardApply(st, 0.1, 1.88, 1.88) # first step, unarmed: taken as is
  expect_equal(sg$p, 1.88)
  st <- sg$state
  sg <- SAIGE::tauGuardApply(st, 1.88, 0.3, 0.3)  # reversal 1
  st <- sg$state
  sg <- SAIGE::tauGuardApply(st, 0.3, 2.0, 2.0)   # reversal 2: armed, step capped at 4 x 0.3
  expect_equal(sg$p, 1.2)
  expect_true(sg$state$fired)
})

# A synthetic AI-REML: the score is exactly -(tau - root) and the reported
# curvature is `curv` times the true one. With curv = 1 Newton lands on the
# root; with curv = 0.4 the step is 2.5 times the distance, so the plain
# iteration reflects through the root and oscillates, which is LEDGER #122's
# failure in its simplest form.
iterateSynthetic <- function(curv, guard, root = 1, start = 0.2, maxit = 40) {
  tau <- start
  st <- SAIGE::tauGuardInit()
  path <- tau
  for (i in seq_len(maxit)) {
    score <- -(tau - root)
    raw <- tau + score / curv
    p <- SAIGE::tauStepBinary(tau, score / curv, tol)
    if (guard) {
      sg <- SAIGE::tauGuardApply(st, tau, raw, p)
      p <- sg$p
      st <- sg$state
    }
    path <- c(path, p)
    if (abs(raw - tau) / (abs(raw) + abs(tau) + tol) < tol) {
      tau <- p
      break
    }
    tau <- p
  }
  list(tau = tau, path = path, its = i, fired = st$fired)
}

test_that("with the true curvature the guarded iteration is the plain one", {
  plain <- iterateSynthetic(1, guard = FALSE)
  guarded <- iterateSynthetic(1, guard = TRUE)
  expect_identical(guarded$path, plain$path)
  expect_false(guarded$fired)
  expect_equal(guarded$tau, 1, tolerance = 1e-6)
})

test_that("#122: under-estimated curvature: the plain step oscillates, the guarded one converges", {
  plain <- iterateSynthetic(0.4, guard = FALSE)
  guarded <- iterateSynthetic(0.4, guard = TRUE)
  expect_true(guarded$fired)
  expect_lt(guarded$its, plain$its)
  expect_equal(guarded$tau, 1, tolerance = 1e-2)
  expect_gt(abs(plain$tau - 1), abs(guarded$tau - 1))
})
