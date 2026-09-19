# The SKAT-O saddlepoint tail (LEDGER #158 step 2, --SKATO_tail_method=saddlepoint)
# at SMALL numbers of terms: Kuonen's relative error is o(n^-3/2), i.e. it IMPROVES
# with n, and the published validation (Chen & Lumley 2019) spans only
# n = 305..11,259 -- while this pipeline's Davies time is concentrated at n = 4-10.
# docs/state/handoff-2026-09-19-saddlepoint-tail-validation.md measured it and the
# worry does not materialise. This file is the guard on that result.
#
# Davies CANNOT be the reference here: SKAT hardcodes acc = 1e-6, and over 840
# points it returned a NON-POSITIVE probability on 277 (33%) while reporting
# ifault = 0 every time. The reference is Ruben's series summed in the UPPER tail
# with beta <= min(lambda), which makes every coefficient nonnegative -- a sum of
# positive terms, with no 1-F cancellation anywhere.
#
# What is asserted is the SHAPE of the result (conservative, driven by depth, flat
# in n), never a pinned ratio: the measured bias is 1.036 at p = 1e-2 rising to
# 1.113 at 1e-12, and a test demanding more precision than that would measure the
# platform's libm, not the estimator.

# Ruben/Farebrother: F(x) = sum_k a_k P(chi2_{n+2k} <= x/beta), sum_k a_k = 1, so
# the survival is sum_k a_k P(chi2_{n+2k} > x/beta). beta <= min(lambda) puts every
# gamma_j in [0,1), hence every a_k >= 0. Truncating at K omits R_K = 1 - sum a_k,
# which bounds the omitted contribution because every chi-square factor is <= 1;
# the function REFUSES rather than guessing when that bound is not small enough.
COEF_MASS_FLOOR <- 1e-14

rubenSurvival <- function(lambda, q, kmax = 8000L, rtol = 1e-6) {
  n <- length(lambda)
  beta <- min(lambda) * (1 - 1e-13)
  gam <- 1 - beta / lambda
  x <- q / beta
  a <- numeric(kmax + 1L)
  g <- numeric(kmax + 1L)
  a[1] <- exp(0.5 * sum(log(beta / lambda)))
  surv <- a[1] * pchisq(x, df = n, lower.tail = FALSE)
  asum <- a[1]
  pw <- rep(1, n)
  for (k in seq_len(kmax)) {
    pw <- pw * gam
    g[k + 1L] <- sum(pw)
    a[k + 1L] <- sum(g[2:(k + 1L)] * a[k:1]) / (2 * k)
    surv <- surv + a[k + 1L] * pchisq(x, df = n + 2 * k, lower.tail = FALSE)
    asum <- asum + a[k + 1L]
    if (k %% 32L == 0L && (1 - asum) < max(rtol * surv, COEF_MASS_FLOOR)) break
  }
  bound <- max(1 - asum, 0)
  # The bound is `1 - sum a_k`, an ABSOLUTE bound on the omitted contribution
  # (every chi-square factor is <= 1). Relative to a survival of ~1e-13 that asks
  # for 1e-19, while `1 - sum a_k` bottoms out at double-precision noise near
  # 1e-16 -- so below COEF_MASS_FLOOR the bound is arithmetic noise, not error,
  # and demanding more would make the reference refuse every deep point. The
  # accuracy that actually matters is demonstrated by the two exact closed forms
  # and by the reconciliation with the independent Python implementation.
  list(p = surv, bound = bound,
       ok = surv > 0 && bound < max(rtol * surv, COEF_MASS_FLOOR))
}

# lambda_j == c  =>  Q = c * chi2_n, exact.
exactEqual <- function(cval, n, q) pchisq(q / cval, df = n, lower.tail = FALSE)

# lambda = (a,a,b,b) => Q = 2a E1 + 2b E2 with E ~ Exp(1), a != b.
exactTwoPairs <- function(a, b, q) {
  al <- 1 / (2 * a)
  be <- 1 / (2 * b)
  (al * exp(-be * q) - be * exp(-al * q)) / (al - be)
}

# Real per-rho spectra from the finished-build guard sweep's own cells, rebuilt as
# skatoScreenKernel builds them; one per term count, lowest condition number so the
# series converges in a few hundred terms. `ref` is the Python Ruben implementation
# (harness/skato_pvalue/tailref/) -- pinned so the two implementations must reconcile.
pinnedSpectra <- list(
  list(
    lambda = c(20.181256553630458, 18.482453812476169, 10.34956696608028,
               8.9915005866287707),
    level = c(1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12),
    q = c(205.114841892518, 294.75657920632398, 384.37354310662101,
          564.18363420578203, 744.93475557399597, 926.50957779779401,
          1108.7335619988601),
    ref = c(9.984883e-03, 9.950465e-04, 9.949188e-05, 9.937814e-07,
            9.867415e-09, 9.753778e-11, 9.624057e-13)
  ),
  list(
    lambda = c(41.69147531339636, 26.978369726287724, 24.126103172983871,
               15.032133945355913, 13.950108937986252, 10.173411242727655),
    level = c(1e-2, 1e-3, 1e-4, 1e-6, 1e-8),
    q = c(411.09422533297101, 587.00231939011303, 766.953375332527,
          1133.77069136968, 1505.60581074692),
    ref = c(1.031048e-02, 9.932149e-04, 9.610722e-05, 9.273206e-07, 9.128470e-09)
  ),
  list(
    lambda = c(75.851373269692985, 32.915971283705673, 28.056846237560798,
               19.165510263136831, 16.597607460479185, 11.642431614072333,
               11.566132125346591, 10.946918036480822),
    level = c(1e-2, 1e-3, 1e-4, 1e-6, 1e-8),
    q = c(669.86686562171894, 989.09063859023399, 1316.9217025677499,
          1985.2519053875201, 2662.3902220367499),
    ref = c(9.732425e-03, 9.442529e-04, 9.307936e-05, 9.164261e-07, 9.079750e-09)
  ),
  list(
    lambda = c(92.124960141268019, 58.322788696336488, 52.863570163765274,
               34.932654417838869, 24.373763125157922, 23.290006163939726,
               17.759233860536007, 14.929481543102561, 13.23397289368525,
               10.794317666257429),
    level = c(1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12),
    q = c(960.51471305366999, 1348.89274572656, 1746.4237417272,
          2556.9269737588002, 3378.5639309457001, 4206.3292016026398,
          5037.9513365250696),
    ref = c(1.030185e-02, 9.902815e-04, 9.582448e-05, 9.259298e-07,
            9.121967e-09, 9.044768e-11, 8.992707e-13)
  )
)

# ---- CONTROLS: the reference must reproduce two exact closed forms -------------

test_that("the Ruben reference reproduces the equal-eigenvalue closed form", {
  grid <- expand.grid(n = c(4L, 7L, 10L), target = c(1e-2, 1e-6, 1e-12))
  got <- exact <- ok <- numeric(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    qq <- 0.7 * qchisq(grid$target[i], df = grid$n[i], lower.tail = FALSE)
    r <- rubenSurvival(rep(0.7, grid$n[i]), qq)
    got[i] <- r$p
    ok[i] <- r$ok
    exact[i] <- exactEqual(0.7, grid$n[i], qq)
  }
  expect_true(all(ok == 1))
  expect_equal(got, exact, tolerance = 1e-10)
})

test_that("the Ruben reference reproduces the two-pair exponential closed form", {
  # distinct eigenvalues, so the gamma recursion is actually exercised (the
  # equal-eigenvalue control terminates the series at gamma = 0 and does not).
  pairs <- list(c(1, 0.5), c(1, 0.2), c(1, 0.1))
  targets <- c(1e-2, 1e-6, 1e-10)
  got <- exact <- ok <- numeric(length(pairs) * length(targets))
  i <- 0L
  for (ab in pairs) {
    for (target in targets) {
      i <- i + 1L
      qq <- uniroot(function(x) exactTwoPairs(ab[1], ab[2], x) - target,
                    c(1e-6, 1e4), tol = 1e-12)$root
      r <- rubenSurvival(c(ab[1], ab[1], ab[2], ab[2]), qq)
      got[i] <- r$p
      ok[i] <- r$ok
      exact[i] <- exactTwoPairs(ab[1], ab[2], qq)
    }
  }
  expect_true(all(ok == 1))
  expect_equal(got, exact, tolerance = 1e-6)
})

test_that("the R and Python Ruben implementations reconcile on the pinned spectra", {
  # the same quantity computed two ways, in two languages: they must agree.
  got <- ref <- ok <- c()
  for (s in pinnedSpectra) {
    for (i in seq_along(s$q)) {
      r <- rubenSurvival(s$lambda, s$q[i])
      got <- c(got, r$p)
      ok <- c(ok, r$ok)
      ref <- c(ref, s$ref[i])
    }
  }
  expect_true(all(ok == 1))
  expect_equal(got, ref, tolerance = 1e-4)
})

# ---- THE GUARD: the shipped saddlepoint against that reference ----------------

ratiosByDepth <- function() {
  out <- list()
  for (si in seq_along(pinnedSpectra)) {
    s <- pinnedSpectra[[si]]
    for (i in seq_along(s$q)) {
      got <- SAIGE:::skatoSaddlepointSurvival(s$lambda, s$q[i])
      out[[length(out) + 1L]] <- data.frame(
        n_lam = length(s$lambda), spectrum = si, level = s$level[i], ref = s$ref[i],
        spa = got, ratio = got / s$ref[i]
      )
    }
  }
  do.call(rbind, out)
}

test_that("the shipped saddlepoint tail is conservative and bounded at small q", {
  r <- ratiosByDepth()
  expect_gt(nrow(r), 20)
  # Measured range over 468 real points was 0.9665 .. 1.1203. The band is wider on
  # purpose: what is being pinned is that the method stays close and on the safe
  # side, not a number that would drift with libm.
  expect_true(all(r$ratio > 0.95), info = paste("min ratio", min(r$ratio)))
  expect_true(all(r$ratio < 1.30), info = paste("max ratio", max(r$ratio)))
  # Conservative where it matters: at and below 1e-6 the reported p is never
  # materially smaller than the truth (measured min there was 1.0063).
  deep <- r[r$ref <= 1e-6, ]
  expect_gt(nrow(deep), 8)
  expect_true(all(deep$ratio > 0.99), info = paste("min deep ratio", min(deep$ratio)))
})

test_that("the saddlepoint bias grows with tail depth", {
  r <- ratiosByDepth()
  shallow <- r[r$ref > 1e-3, ]
  deep <- r[r$ref < 1e-8, ]
  expect_gt(nrow(shallow), 3)
  expect_gt(nrow(deep), 3)
  # depth drives it: the bias grows as the tail deepens (1.036 -> 1.113 measured)
  expect_gt(median(deep$ratio), median(shallow$ratio))
  # At a MATCHED depth (1e-8, the deepest level all four spectra share) every
  # spectrum sits inside the measured within-level spread. NOTE what this is NOT:
  # it is NOT a test that the bias is independent of n. There is one spectrum per
  # term count here, so spectrum-to-spectrum variation (measured max-min 0.088
  # within a level) and any n effect are perfectly confounded. The n-independence
  # evidence is 468 points with spearman(ratio, n_lam) = 0.017, and it lives in
  # docs/state/handoff-2026-09-19-saddlepoint-tail-validation.md, not in a
  # four-spectrum fixture. Asserting it here would be asserting a property this
  # data cannot distinguish.
  at8 <- r[r$level == 1e-8, ]
  expect_equal(nrow(at8), length(pinnedSpectra))
  expect_lt(diff(range(at8$ratio)), 0.15)
})
