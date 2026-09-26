# The exact resampling test (ER; src/ER_binary_func.cpp + the Binary_* machinery)
# on WEIGHTED vectors -- LEDGER #166 (every weight set is now ER-eligible at
# MAC <= --max_MAC_for_ER, as beta sets always were) and #167 (its absolute
# tie epsilon made every configuration a tie under tiny weights; the vector is
# now scaled so its largest entry is 1 before the test).
#
# The kernel ranks the 2^k case/control configurations of the k carriers by
# Q = (sum_i Z_i (y_i - mu_i))^2 and reports a mid-p: the probability of
# those with Q > Q_obs plus half of those tied with it; a configuration's probability is the product of the cases' odds
# mu/(1-mu), normalised within its case count j, times P(j cases among the
# carriers), which the kernel gets from a grouped conditional hypergeometric
# (weight-INDEPENDENT, and exercised by every golden exact-test row). The R
# reference below takes that P(j) from the kernel's own export and does the
# weight-dependent part -- enumeration, odds, ranking -- independently. So a
# match on weighted Z means the kernel treats weights exactly as the
# definition says; a match on unit Z is the control that the reference is
# right where the goldens already are.

er_reference <- function(Z, y, mu, probK, eps = 1e-6) {
  idx <- which(Z != 0); k <- length(idx)
  z <- Z[idx]; m <- mu[idx]; odds <- m / (1 - m)
  confs <- as.matrix(expand.grid(rep(list(0:1), k)))
  j <- rowSums(confs)
  w <- apply(confs, 1, function(cc) prod(odds[cc == 1]))
  prob <- numeric(nrow(confs))
  for (jj in 0:k) {
    sel <- j == jj
    prob[sel] <- w[sel] / sum(w[sel]) * probK[jj + 1]
  }
  Q <- as.numeric((confs %*% z - sum(z * m))^2)
  Qobs <- (sum(z * (y[idx] - m)))^2
  # the kernel returns a MID-p: P(Q > Qobs) + P(Q == Qobs) / 2, ties by the
  # absolute epsilon (SKATExactBin_Work: `pval[0] - pval1[0]/2`)
  sum(prob[Q > Qobs + eps]) + 0.5 * sum(prob[abs(Q - Qobs) <= eps])
}

er_case <- function(seed, k, n = 300, weighted = TRUE) {
  set.seed(seed)
  mu <- runif(n, 0.01, 0.3)
  y <- rbinom(n, 1, mu)
  Z <- numeric(n); idx <- sample(n, k)
  Z[idx] <- if (weighted) runif(k, 0.05, 1) * sample(1:2, k, replace = TRUE) else sample(1:2, k, replace = TRUE)
  list(Z = Z, y = y, mu = mu)
}

test_that("ER reference == kernel on UNIT genotype vectors (the control)", {
  for (seed in 1:20) for (k in c(1, 2, 4, 6)) {
    d <- er_case(seed, k, weighted = FALSE)
    r <- erExactTestInCPP(d$Z, d$y, d$mu, 2000000L, TRUE)
    expect_true(r$exact)
    expect_equal(r$pval, er_reference(d$Z, d$y, d$mu, r$probK), tolerance = 1e-10)
  }
})

test_that("#166: ER reference == kernel on WEIGHTED vectors", {
  for (seed in 1:30) for (k in c(1, 2, 3, 5, 8)) {
    d <- er_case(seed, k, weighted = TRUE)
    r <- erExactTestInCPP(d$Z, d$y, d$mu, 2000000L, TRUE)
    expect_true(r$exact)
    expect_equal(r$pval, er_reference(d$Z, d$y, d$mu, r$probK), tolerance = 1e-10)
  }
})

test_that("#167 FIXED: the kernel's p is invariant to the vector's scale over eleven orders", {
  d <- er_case(7, 5)
  base <- erExactTestInCPP(d$Z, d$y, d$mu, 2000000L, TRUE)$pval
  for (cc in c(1e-8, 1e-6, 1e-4, 1e-2, 0.5, 2, 7.8, 1e3)) {
    expect_identical(erExactTestInCPP(cc * d$Z, d$y, d$mu, 2000000L, TRUE)$pval, base)
  }
})

test_that("#167 characterization: without the normalisation, tiny weights make every configuration a tie", {
  base <- tiny <- numeric(6)
  for (seed in 1:6) {
    d <- er_case(seed, 5)
    base[seed] <- erExactTestInCPP(d$Z, d$y, d$mu, 2000000L, FALSE)$pval
    tiny[seed] <- erExactTestInCPP(1e-6 * d$Z, d$y, d$mu, 2000000L, FALSE)$pval
  }
  # the data-driven p's spread; the scaled ones are all the same number:
  # every configuration ties with the observed one, so the mid-p is exactly
  # half of the whole mass, 0.5 whatever the data
  expect_gt(max(abs(base - 0.5)), 0.1)
  expect_equal(tiny, rep(0.5, 6), tolerance = 1e-9)
})

test_that("#167: weights spanning six orders inside one vector agree with the reference", {
  set.seed(11); n <- 300; mu <- runif(n, 0.01, 0.3); y <- rbinom(n, 1, mu)
  Z <- numeric(n); idx <- sample(n, 6); Z[idx] <- c(1, 0.3, 1e-2, 1e-4, 1e-6, 0.7)
  r <- erExactTestInCPP(Z, y, mu, 2000000L, TRUE)
  # the reference applies the kernel's own tie rule after the same scaling
  expect_equal(r$pval, er_reference(Z, y, mu, r$probK), tolerance = 1e-10)
})
