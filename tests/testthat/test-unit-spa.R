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
  # so the degenerate branch is not silently changed. LEDGER #22 FIXED: it
  # reports non-convergence now.
  d <- spa_fixture()
  r <- getroot_K1_Binom(0, d$mu, rep(0, length(d$g)), 0, 1e-8, 100L)
  expect_true(is.infinite(r$root))
  expect_false(r$Isconverge)   # LEDGER #22 FIXED: no root, so not "converged"
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
  # which since LEDGER #22 was FIXED it reports as NOT converged (it used to
  # say converged in 0 iterations). SPA_binary then takes its no-saddlepoint
  # fallback: the unadjusted p-value exactly, Isconverge FALSE. The composed
  # p is the same value the old route reassembled from two half-tails.
  d <- spa_fixture()
  m <- sum(d$mu * d$g)
  huge <- sum(abs(d$g)) + 1e3
  expect_true(is.infinite(getroot_K1_Binom(0, d$mu, d$g, huge, 1e-8, 100L)$root))
  r <- SPA_binary(d$mu, d$g, huge, 2 * m - huge, 0.0123, 1e-8, FALSE)
  expect_false(r$Isconverge)
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

# LEDGER #163 (2026-09-26, FIXED the same day): the Newton root finder
# starts at t = 0, where K2(0) = sum(mu (1 - mu) g^2) is the score variance.
# When the carriers' fitted probabilities are tiny, that first step
# (q - K1(0)) / K2(0) exceeds ~730, exp(-g t) underflows, K2 evaluates to
# exactly 0, the next iterate is -Inf, and 0 * Inf in the following K1 makes
# NaN: non-convergence for a root that exists (q is strictly inside
# (gneg, gpos), K1 is monotone). Two carriers at mu = 6e-4, one of them a
# case, is a SYNTHETIC shape of that failure; the same two at mu = 1e-3 never
# failed. The fix is a bracketed, safeguarded Newton that runs ONLY after the
# plain iteration has failed, so every root the plain Newton finds keeps its
# bytes. For two unit carriers the root is closed-form: log((1 - mu) / mu).
spa_tiny_mu_case <- function(m, n = 2000, centred = FALSE) {
  mu <- rep(0.07, n); g <- rep(0, n)
  mu[1:2] <- m; g[1:2] <- 1
  # centred = the shape of a covariate-ADJUSTED genotype (gtilde): the
  # projection on an intercept makes the non-carriers slightly negative so the
  # entries sum to ~0, which is what makes the MIRRORED q attainable
  # (qinv = 2 m1 - q > gneg). The uncentred vector has gneg = 0, so its
  # mirrored tail is unattainable by construction and lands in #164.
  if (centred) g[3:n] <- -sum(g[1:2]) / (n - 2)
  # one case carrier, one control carrier: q = sum(g * y) = 1 uncentred;
  # centred, the non-carrier cases contribute their (tiny) g's
  set.seed(163); y <- rbinom(n, 1, mu); y[1] <- 1; y[2] <- 0
  list(mu = mu, g = g, q = sum(g * y), m1 = sum(mu * g))
}

test_that("#163 control: the root finder converges for carriers at mu = 1e-3", {
  d <- spa_tiny_mu_case(1e-3)
  r <- getroot_K1_Binom(0, d$mu, d$g, d$q, .Machine$double.eps^0.25, 1000L)
  expect_true(r$Isconverge)
  expect_equal(r$root, log((1 - 1e-3) / 1e-3), tolerance = 1e-4)
})

test_that("#163 FIXED: the root finder converges for carriers at mu = 6e-4 and 1e-4", {
  for (m in c(6e-4, 1e-4)) {
    d <- spa_tiny_mu_case(m)
    # the first plain Newton step is past the exp underflow point ...
    expect_gt((d$q - d$m1) / sum(d$mu * (1 - d$mu) * d$g^2), 730)
    # ... and the rescue still lands on the closed-form root
    r <- getroot_K1_Binom(0, d$mu, d$g, d$q, .Machine$double.eps^0.25, 1000L)
    expect_true(r$Isconverge)
    expect_equal(r$root, log((1 - m) / m), tolerance = 1e-4)
    expect_equal(K1_adj_Binom(r$root, d$mu, d$g, d$q), 0, tolerance = 1e-6)
  }
})

test_that("#163 FIXED: the composed saddlepoint p is near the exact tail, not the normal p", {
  d <- spa_tiny_mu_case(6e-4, centred = TRUE)
  # both roots exist: q and its mirror are strictly inside (gneg, gpos)
  gpos <- sum(d$g[d$g > 0]); gneg <- sum(d$g[d$g < 0])
  expect_true(d$q < gpos && d$q > gneg)
  expect_true(2 * d$m1 - d$q < gpos && 2 * d$m1 - d$q > gneg)
  S <- d$q - d$m1
  v <- sum(d$mu * (1 - d$mu) * d$g^2)
  p_na <- pchisq(S^2 / v, 1, lower.tail = FALSE)
  qinv <- -abs(d$q - d$m1) + d$m1
  conv <- FALSE
  p <- SPA_pval(d$mu, d$g, d$q, qinv, p_na, .Machine$double.eps^0.25, FALSE, "binary", conv)
  exact <- 1 - prod(1 - d$mu[1:2])   # one-sided: at least one carrier is a case
  expect_lt(p_na, 1e-100)
  expect_gt(p, exact / 10)
  expect_lt(p, exact * 10)
})

# LEDGER #164 (OPEN): when the saddlepoint does not converge on a binary
# trait, SPA_pval returns pval_noadj -- the normal approximation -- as THE
# p-value, and getMarkerPval prints it in p.value with Is.SPA = FALSE. With
# #163 fixed the tiny-mu shape no longer reaches this; the unattainable-q
# rule (LEDGER #22: q outside (gneg, gpos), no root exists) still does.
# Characterization; invert when the row is decided and fixed.
test_that("#164: SPA non-convergence hands back the normal-approximation p", {
  d <- spa_tiny_mu_case(0.07)
  gpos <- sum(d$g[d$g > 0])
  q <- gpos + 0.5                       # unattainable: K1 never reaches it
  qinv <- -abs(q - d$m1) + d$m1
  p_na <- 1e-120
  conv <- TRUE
  p <- SPA_pval(d$mu, d$g, q, qinv, p_na, .Machine$double.eps^0.25, FALSE, "binary", conv)
  expect_equal(p, p_na)
})

# LEDGER #168 (2026-09-26, from the AoU fallback sidecar: 172 of 245 failed
# saddlepoints were "root(mirrored q) not found"). A REAL adjusted genotype is
# centred with the working weights W = mu (1 - mu): the non-carriers' negative
# entries sum to the carriers' W-weighted mass, and for carriers with mu ~ 1e-3
# that mass is ~1e-2 while the mirrored q is ~-1 -- outside the attainable
# range. That tail has probability exactly 0 (S cannot fall below gneg), which
# is how SPAtest treats it; LEDGER #22's rule made it a non-convergence and
# the whole test fell back. The earlier synthetic vector centred WITHOUT the
# weights (negative mass -2) never reached this, which is why the rescue
# looked complete on it.
spa_w_centred_case <- function(m, n = 2000) {
  mu <- rep(0.07, n); g <- rep(0, n)
  mu[1:2] <- m; g[1:2] <- 1
  W <- mu * (1 - mu)
  g[3:n] <- -sum(W[1:2] * g[1:2]) / sum(W[3:n])        # W-weighted centring: sum(W g) = 0
  set.seed(168); y <- rbinom(n, 1, mu); y[1] <- 1; y[2] <- 0
  list(mu = mu, g = g, q = sum(g * y), m1 = sum(mu * g))
}

test_that("#168 FIXED: an unattainable mirrored q yields the one-sided saddlepoint p, not a fallback", {
  d <- spa_w_centred_case(6e-4)
  gneg <- sum(d$g[d$g < 0]); qinv <- 2 * d$m1 - d$q
  expect_lt(qinv, gneg)                      # the mirror IS unattainable on this shape
  S <- d$q - d$m1; v <- sum(d$mu * (1 - d$mu) * d$g^2)
  p_na <- pchisq(S^2 / v, 1, lower.tail = FALSE)
  conv <- FALSE
  p <- SPA_pval(d$mu, d$g, d$q, qinv, p_na, .Machine$double.eps^0.25, FALSE, "binary", conv)
  expect_lt(p_na, 1e-100)
  expect_false(isTRUE(all.equal(p, p_na)))
  exact <- 1 - prod(1 - d$mu[1:2])            # one-sided: at least one carrier is a case
  expect_gt(p, exact / 10); expect_lt(p, exact * 10)
})

# LEDGER #163 amendment (from the same sidecar: 73 "root(q) not found"): the
# root scales as 1/(vector scale), so the rescue's bracket cap is in units of
# t * max|g|. Under weights of 1e-6 the root is near 1e7 and used to be past
# the 2^20 cap.
test_that("#163 amendment: the rescued root finder converges under weights of 1e-6", {
  d <- spa_tiny_mu_case(6e-4, centred = TRUE)
  for (w in c(1, 1e-3, 1e-6)) {
    r <- getroot_K1_Binom(0, d$mu, w * d$g, w * d$q, .Machine$double.eps^0.25, 1000L)
    expect_true(r$Isconverge)
    # converged means the remaining Newton step |K1/K2| is below the finder's
    # absolute tolerance on t (K1's residual itself scales with the weight)
    step <- abs(K1_adj_Binom(r$root, d$mu, w * d$g, w * d$q) / K2_Binom(r$root, d$mu, w * d$g))
    expect_lt(step, .Machine$double.eps^0.25)
    expect_equal(r$root * w, 7.4158, tolerance = 1e-4)   # the root scales as 1/w
  }
})

# LEDGER #168 (b): the mirrored q attainable by a HAIR. The vector is the
# W-centred shape with its negative mass scaled to just beyond the mirrored
# q, so the mirror's root exists but lies past any finder (K1 is within a
# hair of its limit only at |t| ~ 1e6 and beyond). The Chernoff bound on
# that tail is astronomically small, so the composed p is the q-side tail.
# From the second SAS rerun: 24 fallbacks, all "root(mirrored q) not found"
# by the finder, none by the attainability rule.
test_that("#168 (b) FIXED: a mirrored root beyond the finder is bounded away, not a fallback", {
  # The mirrored root's distance grows linearly with N when the mirror is
  # attainable by a fixed margin (measured: -1.6e4 at N = 2,000, -6.5e4 at
  # N = 8,000), so at N = 1.2e6 it sits near -1e7, beyond the bracket's cap
  # of 2^20 / max|g|. The finder fails; the Chernoff bound on that tail is
  # astronomically small; the composed p is the q-side tail.
  n <- 1200000
  mu <- rep(0.07, n); g <- rep(0, n); mu[1:2] <- 6e-4; g[1:2] <- 1
  W <- mu * (1 - mu); g[3:n] <- -sum(W[1:2] * g[1:2]) / sum(W[3:n])
  set.seed(168); y <- rbinom(n, 1, mu); y[1] <- 1; y[2] <- 0
  q <- sum(g * y)
  neg <- g < 0; marg <- 1.001
  A <- sum(mu[!neg] * g[!neg]); B <- sum(mu[neg] * g[neg]); g0 <- sum(g[neg])
  g[neg] <- g[neg] * (marg * (2 * A - q) / (g0 - 2 * marg * B))
  m1 <- sum(mu * g); qinv <- 2 * m1 - q
  expect_gt(qinv, sum(g[neg]))                                 # attainable ...
  r2 <- getroot_K1_Binom(0, mu, g, qinv, .Machine$double.eps^0.25, 1000L)
  expect_false(r2$Isconverge)                                  # ... beyond the finder
  r1 <- getroot_K1_Binom(0, mu, g, q, .Machine$double.eps^0.25, 1000L)
  expect_true(r1$Isconverge)
  S <- q - m1; v <- sum(mu * (1 - mu) * g^2)
  p_na <- pchisq(S^2 / v, 1, lower.tail = FALSE)
  conv <- FALSE
  p <- SPA_pval(mu, g, q, qinv, p_na, .Machine$double.eps^0.25, FALSE, "binary", conv)
  expect_false(isTRUE(all.equal(p, p_na)))
  exact <- 1 - prod(1 - mu[1:2])
  expect_gt(p, exact / 20); expect_lt(p, exact * 20)
})

# LEDGER #169 (2026-09-27, from the 20 BenIntNeo saddlepoint dumps analysed
# where the data lives): the q-side root converges, the mirror is
# unattainable (#168 (a)), and the saddle at that root is refused because
# K(t) = sum log(1 - mu + mu e^{g t}) overflows to Inf, while K1 and K2 --
# written in e^{-g t} -- stay finite, so the Newton finder is happy where the
# CGF cannot be evaluated. The root sits there whenever the cell's positive
# mass is dominated by carriers with tiny fitted probabilities and the
# observed S is within a hair of gpos: root * max|g| is of order 1/mu of the
# dominant carrier (the dumps: 730 to 6,000; exp overflows at ~709). The
# fallback was labelled "root(mirrored q) not found" by the final else, which
# is what hid it. The fix rewrites the overflowing terms as
# x + log((1 - mu) e^{-x} + mu), ONLY when the plain sum is non-finite.
# Reference: exact enumeration of every configuration of a 20-sample vector.
spa_boundary_case <- function(mu1, ncase_rest) {
  n <- 20
  mu <- c(mu1, seq(0.05, 0.15, length.out = n - 1)); W <- mu * (1 - mu)
  g <- c(1, seq(5e-4, 1.5e-3, length.out = n - 1))
  g[-1] <- g[-1] - 1.2 * sum(W * g) / sum(W[-1])       # over-centred: the rest is negative
  y <- rep(0, n); y[1] <- 1                             # the dominant carrier is a case
  if (ncase_rest > 0) y[n - seq_len(ncase_rest) + 1] <- 1
  q <- sum(g * y); m1 <- sum(mu * g)
  sums <- 0; pr <- 1
  for (i in seq_len(n)) { sums <- c(sums, sums + g[i]); pr <- c(pr * (1 - mu[i]), pr * mu[i]) }
  list(mu = mu, g = g, q = q, m1 = m1, qinv = 2 * m1 - q, gneg = sum(g[g < 0]),
       exact_up = sum(pr[sums >= q - 1e-12]), exact_low = sum(pr[sums <= 2 * m1 - q + 1e-12]))
}

test_that("#169 FIXED: the CGF is finite at a root past exp's range and the saddle is taken", {
  d <- spa_boundary_case(1e-3, 1)
  tol <- .Machine$double.eps^0.25
  expect_lt(d$qinv, d$gneg); expect_equal(d$exact_low, 0)       # #168 (a) shape: no lower tail
  r1 <- getroot_K1_Binom(0, d$mu, d$g, d$q, tol, 1000L)
  expect_true(r1$Isconverge)
  expect_gt(r1$root * max(abs(d$g)), 709)                        # past exp's range
  expect_true(is.infinite(sum(log(1 - d$mu + d$mu * exp(d$g * r1$root)))))   # the plain sum overflows
  K <- Korg_Binom(r1$root, d$mu, d$g)
  expect_true(is.finite(K))                                      # the kernel no longer does
  x <- d$g * r1$root; m <- pmax(x, 0)
  expect_equal(K, sum(m + log((1 - d$mu) * exp(-m) + d$mu * exp(x - m))), tolerance = 1e-12)
  expect_true(Get_Saddle_Prob_Binom(r1$root, d$mu, d$g, d$q, FALSE)$isSaddle)
  p_na <- 1e-100; conv <- FALSE
  p <- SPA_pval(d$mu, d$g, d$q, d$qinv, p_na, tol, FALSE, "binary", conv)
  expect_false(isTRUE(all.equal(p, p_na)))                       # not a fallback
  # measured: 2.16e-4 against an exact 1.57e-4 (1.37x); the same shape one notch
  # inside exp's range is 1.4-1.8x on the same side, so this is the saddlepoint's
  # own accuracy at the boundary of the support, not something the fix moved
  expect_gt(p, d$exact_up / 2); expect_lt(p, d$exact_up * 2)
})

test_that("#169: the fast kernel's CGF is finite at the same root", {
  d <- spa_boundary_case(1e-3, 1)
  tol <- .Machine$double.eps^0.25
  # the dominant carrier is the only raw-nonzero sample; the rest folds to a normal
  NAmu <- sum(d$mu[-1] * d$g[-1]); NAsigma <- sum(d$mu[-1] * (1 - d$mu[-1]) * d$g[-1]^2)
  rf <- getroot_K1_fast_Binom(0, d$mu, d$g, d$q, d$g[-1], d$g[1], d$mu[-1], d$mu[1], NAmu, NAsigma, tol, 1000L)
  expect_true(rf$Isconverge); expect_gt(rf$root * max(abs(d$g)), 709)
  Kf <- Korg_fast_Binom(rf$root, d$mu, d$g, d$g[-1], d$g[1], d$mu[-1], d$mu[1], NAmu, NAsigma)
  expect_true(is.finite(Kf))
  expect_equal(Kf, rf$root + log((1 - d$mu[1]) * exp(-rf$root) + d$mu[1]) + NAmu * rf$root + 0.5 * NAsigma * rf$root^2,
               tolerance = 1e-12)
  s <- Get_Saddle_Prob_fast_Binom(rf$root, d$mu, d$g, d$q, d$g[-1], d$g[1], d$mu[-1], d$mu[1], NAmu, NAsigma, FALSE)
  expect_true(s$isSaddle)
  expect_gt(abs(s$pval), d$exact_up / 3); expect_lt(abs(s$pval), d$exact_up * 3)
})

test_that("#169 control: inside exp's range the CGF is the plain sum and the composed p is unchanged", {
  d <- spa_boundary_case(1e-3, 2)                                 # same shape, root ~530
  tol <- .Machine$double.eps^0.25
  r1 <- getroot_K1_Binom(0, d$mu, d$g, d$q, tol, 1000L)
  expect_true(r1$Isconverge); expect_lt(r1$root * max(abs(d$g)), 709)
  plain <- sum(log(1 - d$mu + d$mu * exp(d$g * r1$root)))
  expect_true(is.finite(plain))
  expect_equal(Korg_Binom(r1$root, d$mu, d$g), plain, tolerance = 1e-12)
  conv <- FALSE
  p <- SPA_pval(d$mu, d$g, d$q, d$qinv, 1e-100, tol, FALSE, "binary", conv)
  expect_gt(p, d$exact_up / 2); expect_lt(p, d$exact_up * 2)      # measured 1.77x, as before the fix
})

test_that("#169: S at its maximum (every carrier a case) is still the exact test's, not a saddle", {
  d <- spa_boundary_case(1e-3, 0)                                 # q == gpos: no root
  tol <- .Machine$double.eps^0.25
  r1 <- getroot_K1_Binom(0, d$mu, d$g, d$q, tol, 1000L)
  expect_false(r1$Isconverge); expect_equal(r1$niter, 0L)
  conv <- TRUE
  p <- SPA_pval(d$mu, d$g, d$q, d$qinv, 1e-100, tol, FALSE, "binary", conv)
  expect_equal(p, 1e-100)
})
