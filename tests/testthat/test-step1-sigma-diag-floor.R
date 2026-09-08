# Characterization of LEDGER #91: Sigma's diagonal is floored at an ABSOLUTE
# 1e-4 in phenotype units.
#
# Sigma = tau0/W + tau1*GRM is the matrix every step-1 solve inverts. Three
# sites in src/SAIGE_fitGLMM_fast.cpp clamp its diagonal entries up to 1e-4:
# gen_sp_Sigma (which builds the matrix the DIRECT sparse solve factorises, so
# the clamp changes the answer) and getDiagOfSigma / getDiagOfSigma_LOCO (which
# only build the Jacobi preconditioner 1/diag, so there the clamp reaches the
# answer only through the absolute PCG stopping rule of LEDGER #6).
#
# 1e-4 is a constant in phenotype units, so whether it fires depends on the
# scale of the trait and not on anything about the model. A quantitative trait
# whose covariate-adjusted sd is below ~0.01 has tau under the floor and gets a
# DIFFERENT Sigma, silently. A binary trait cannot reach it: tau0 = 1 and
# W = mu(1-mu) <= 1/4, so tau0/W >= 4.
#
# Since LEDGER #65, fitNULLGLMM standardises every quantitative trait to unit
# adjusted variance before any solve, so the floor was unreachable from the CLI.
# It stayed live for a direct caller of the fitters. These tests are that direct
# caller: gen_sp_Sigma is exported and reads the global sparse GRM, which
# setupSparseGRM installs, so no fit is needed to reach the clamp.
#
# LEDGER #91 FIXED: gen_sp_Sigma no longer clamps at all (a singular Sigma now
# reaches the "spsolve(): solution not found" its only caller already threw),
# and the two preconditioner sites clamp only a non-positive or non-finite
# entry, to 1. The assertions marked INVERTED pinned the defect; they now pin
# that Sigma is the model's own matrix at every scale.

# A tiny symmetric sparse GRM: identity plus one related pair.
tiny_grm <- function(n = 6) {
  i <- c(seq_len(n), 1L, 2L)
  j <- c(seq_len(n), 2L, 1L)
  x <- c(rep(1, n), 0.5, 0.5)
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
}

install_tiny_grm <- function(m) {
  a <- Matrix::summary(methods::as(m, "TsparseMatrix"))
  # 0-based, as R/SAIGE_fitGLMM_fast.R:1027 passes them
  SAIGE:::setupSparseGRM(nrow(m), rbind(a$i - 1L, a$j - 1L), a$x)
  invisible(NULL)
}

sigma_diag <- function(tau, w) {
  m <- tiny_grm(length(w))
  install_tiny_grm(m)
  s <- SAIGE:::gen_sp_Sigma(as.numeric(w), as.numeric(tau))
  Matrix::diag(s)
}

test_that("the control holds: above the floor, Sigma's diagonal is tau0/W + tau1*GRM", {
  n <- 6
  w <- rep(1, n)
  tau <- c(0.5, 0.25) # tau0/w = 0.5, well above 1e-4
  got <- sigma_diag(tau, w)
  expect_equal(got, rep(tau[1] / 1 + tau[2] * 1, n), tolerance = 1e-6)
})

test_that("the control holds on a binary-shaped W, which can never reach the floor", {
  n <- 6
  mu <- c(0.5, 0.3, 0.1, 0.01, 0.001, 1e-4)
  w <- mu * (1 - mu)
  tau <- c(1, 0.3)
  got <- sigma_diag(tau, w)
  expect_equal(got, tau[1] / w + tau[2] * 1, tolerance = 1e-6)
  expect_true(all(got >= 4)) # tau0 = 1, W <= 1/4
})

test_that("LEDGER #91 FIXED: a small-variance quantitative trait keeps its Sigma", {
  n <- 6
  w <- rep(1, n)
  # the bundled quantitative phenotype divided by 100: tau ~ (2.2e-5, 4.6e-5),
  # so every diagonal entry is under the absolute floor
  tau <- c(2.2e-5, 4.6e-5)
  exact <- tau[1] / w + tau[2] * 1
  expect_true(all(exact < 1e-4)) # the fixture really is under the floor
  got <- sigma_diag(tau, w)
  # INVERTED (#91 fixed): the model's own Sigma, not a clamped stand-in.
  expect_equal(got, exact, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(got, rep(1e-4, n), tolerance = 1e-6)))
  # the old clamp put this diagonal ~1.5x too high, and the ratio grew
  # without bound as the trait was scaled down; it is now exactly 1
  expect_equal(mean(got / exact), 1, tolerance = 1e-6)
})

test_that("LEDGER #91 FIXED: Sigma's diagonal is exactly scale-equivariant", {
  # The same model on a trait scaled by s has tau scaled by s^2, so every
  # diagonal entry scales by s^2 and Sigma^-1 by 1/s^2. That is exact
  # arithmetic -- until the floor fires on one of the two scales.
  n <- 6
  w <- rep(1, n)
  big <- c(0.5, 0.25)
  s2 <- 1e-4
  small <- big * s2
  d_big <- sigma_diag(big, w)
  d_small <- sigma_diag(small, w)
  # INVERTED (#91 fixed): the ratio is s2, exactly, at any scale.
  expect_equal(d_small / d_big, rep(s2, n), tolerance = 1e-6)
  expect_false(isTRUE(all.equal(d_small, rep(1e-4, n), tolerance = 1e-6)))
})
