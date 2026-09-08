# Characterization of LEDGER #104 (and #105, the same mechanism): on binary
# traits the region path's Phi adjustment computes its burden p-value on the
# NAIVE burden variance var1 = sum(mu (1 - mu) g.sum^2) -- no GRM, no variance
# ratio -- and then rescales the whole of Phi by 1 / min(VarQ / VarQ_2, 1),
# where VarQ = 1' D^1/2 Phi D^1/2 1 (the mixed-model burden variance as C++
# delivered it) and VarQ_2 = S^2 / qchisq(p_burden). Below the |z| < 2 cutoff
# p_burden is the plain chi-square p on var1, so VarQ_2 == var1 exactly and
# the adjusted burden variance is max(VarQ, var1): wherever the GRM and the
# ratio REDUCE the variance, the region path reports the naive burden test
# and inflates all of Phi with it (every cell of golden P6).
#
# Synthetic inputs, self-consistent the way region_postprocess builds them:
# q.sum - m1 == sum(Score) and p.new the normal p on Phi's own diagonal (so the
# per-marker calibration D is exactly 1). Phi is scaled so that 1'Phi 1 is a
# chosen multiple of var1.
#
# LEDGER #104/#105 FIXED (fix A: the burden statistic is standardised by
# 1'Phi 1 before the cutoff test and the saddlepoint, the C++ marker SPA's own
# device). The assertions marked INVERTED pinned the defect; they now pin that
# the adjusted burden variance is the mixed-model one wherever no saddlepoint
# fires, and that a scalar on Phi is no longer cancelled.

phi_case <- function(ratio, z_target) {
  set.seed(104)
  n <- 400; q <- 3
  mu <- runif(n, 0.05, 0.30)
  G <- matrix(rbinom(n * q, 2, 0.02), n, q)
  w <- c(1.0, 0.8, 1.3)
  Gw <- sweep(G, 2, w, "*")
  g.sum <- rowSums(Gw)
  var1 <- sum(mu * (1 - mu) * g.sum^2)
  # a PSD Phi with off-diagonal mass, scaled to 1'Phi 1 = ratio * var1
  corMat <- matrix(0.3, q, q); diag(corMat) <- 1
  sdVec <- sqrt(colSums(mu * (1 - mu) * Gw^2))
  Phi <- (sdVec %o% sdVec) * corMat
  Phi <- Phi * (ratio * var1 / sum(Phi))
  # scores with the requested naive z, split across markers
  S <- z_target * sqrt(var1)
  Score <- S * c(0.5, 0.3, 0.2)
  m1 <- sum(mu * g.sum)
  q.sum <- m1 + S
  p.new <- pchisq(Score^2 / diag(Phi), df = 1, lower.tail = FALSE, log.p = TRUE)
  list(q.sum = q.sum, mu = mu, g.sum = g.sum, p.new = p.new, Score = Score, Phi = Phi, var1 = var1)
}

adjusted_burden_var <- function(tc) {
  out <- SAIGE:::get_newPhi_scaleFactor_traitType(tc$q.sum, tc$mu, tc$g.sum, tc$p.new, tc$Score,
                                                    tc$Phi, "SKAT-O", "binary")
  sum(out$val)
}

test_that("the fixture is self-consistent: D = 1 and |z| below the cutoff", {
  tc <- phi_case(0.9, 1.0)
  expect_equal(sum(tc$Phi) / tc$var1, 0.9)
  expect_lt(abs(tc$q.sum - sum(tc$mu * tc$g.sum)) / sqrt(tc$var1), 2)
  D <- (tc$Score^2 / qchisq(tc$p.new, df = 1, lower.tail = FALSE, log.p = TRUE)) / diag(tc$Phi)
  expect_equal(D, rep(1, 3), tolerance = 1e-12)
})

test_that("where the mixed model REDUCES the burden variance, the region path keeps it", {
  tc <- phi_case(0.9, 1.0)
  # INVERTED (#104 fixed): sum(Phi), the mixed-model variance C++ delivered,
  # and NOT the naive var1 any more.
  expect_equal(adjusted_burden_var(tc), sum(tc$Phi), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(adjusted_burden_var(tc), tc$var1, tolerance = 1e-6)))
})

test_that("where the mixed model INFLATES it, the clamp leaves Phi alone (holds both ways)", {
  tc <- phi_case(1.1, 1.0)
  expect_equal(adjusted_burden_var(tc), sum(tc$Phi), tolerance = 1e-12)
})

test_that("negative control: with Phi the naive covariance nothing is rescaled (holds both ways)", {
  tc <- phi_case(1.0, 1.0)
  expect_equal(adjusted_burden_var(tc), sum(tc$Phi), tolerance = 1e-12)
  expect_equal(adjusted_burden_var(tc), tc$var1, tolerance = 1e-12)
})

test_that("a scalar on Phi is no longer cancelled by the rescale", {
  # Evidence 2 of the diagnosis (pairs R3/R6): before the fix, halving every
  # entry of Phi left the adjusted burden variance unchanged, because r halved
  # with it and the clamp then restored the naive value. INVERTED (#104
  # fixed): the adjusted variance follows Phi.
  tc <- phi_case(0.9, 1.0)
  half <- tc; half$Phi <- tc$Phi / 2
  half$p.new <- pchisq(half$Score^2 / diag(half$Phi), df = 1, lower.tail = FALSE, log.p = TRUE)
  expect_equal(adjusted_burden_var(half), adjusted_burden_var(tc) / 2, tolerance = 1e-12)
})
