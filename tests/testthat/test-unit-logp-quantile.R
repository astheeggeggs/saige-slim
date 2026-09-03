# The SE back-calculation's two-sided quantile on the LOG p-value scale
# (LEDGER #42), and the evidence that the regime reaching it is real.
#
# WHY THIS FILE EXISTS RATHER THAN A GOLDEN. The log format is chosen only when
# the score p-value UNDERFLOWS to exactly 0.0 in double -- chi-sq beyond ~1490,
# |z| beyond ~38.6. The bundled fixtures have 1,000 samples and their strongest
# marker is p ~ 1e-2 (|z| ~ 2.6); a perfectly separating marker at n = 1000
# caps out near chi-sq 500. So NO golden can reach these lines, and one that
# stayed byte-identical after the fix would say nothing at all. As with the S2
# seed counter and the S5 reader hooks, the change carries its own measurement.
#
# Tier: fast (pure kernels, no CLI).

test_that("the log-scale quantile takes log(p) - log(2), not log(p)/2", {
  # The defect: halving a LOG is a square root. qnorm(log(p)/2, log.p = TRUE)
  # evaluates the quantile at sqrt(p) rather than p/2, which is a different
  # number by ~sqrt(2) in this regime -- and the SE is |BETA| / qval, so SE came
  # out ~1.42x too large exactly on the strongest hits in a scan.
  for (lp in c(-710, -827, -1500, -3000)) {
    got <- testQnormTwoSidedLogPInCPP(lp)
    expect_equal(got, qnorm(lp - log(2), lower.tail = FALSE, log.p = TRUE),
                 tolerance = 1e-12, label = paste("log p =", lp))
    # ... and it is NOT what the code used to compute. Without this the test
    # would pass on a build that still halved the log at, say, lp = -1e-8.
    asCodedBefore <- qnorm(lp / 2, lower.tail = FALSE, log.p = TRUE)
    expect_gt(got / asCodedBefore, 1.4)
    expect_lt(got / asCodedBefore, 1.42)
  }
})

test_that("it agrees with the linear-scale sibling wherever both are defined", {
  # The two arms compute the same statistic and must not disagree. Only the log
  # arm was wrong, so this is the reconciliation that the underflow regime broke:
  # anything computable two ways must reconcile.
  for (p in c(0.5, 1e-3, 1e-30, 1e-200, 1e-300)) {
    expect_equal(testQnormTwoSidedLogPInCPP(log(p)),
                 abs(qnorm(p / 2, lower.tail = FALSE)),
                 tolerance = 1e-9, label = paste("p =", p))
  }
})

test_that("the saddlepoint really does converge in the regime that reaches it", {
  # REACHABILITY, which LEDGER #42 originally carried as an [ASSUMPTION]: the
  # defective line sits inside `if(t_isSPAConverge)`, so it matters whether the
  # SPA converges at all once the score p has underflowed. It does.
  #
  # Note this is also why a doctored variance-ratio file is NOT a way to reach
  # the path from the CLI: q is scaled by 1/sqrt(varRatio), so shrinking the
  # ratio pushes q outside the CGF's range and the root finder declines --
  # measured, 23 of 70 markers reached the log format that way and every one of
  # them came back Is.SPA = false. Self-consistent data is what reaches it, and
  # that needs a sample size the bundled fixtures do not have.
  n <- 20000
  mu <- rep(0.5, n)
  g <- rep(c(1, 0), each = n / 2)
  m1 <- sum(mu * g)
  q <- 7000                                  # inside (K1(0) = 5000, sum(g) = 10000)
  z <- (q - m1) / sqrt(sum(mu * (1 - mu) * g^2))
  lp <- pchisq(z^2, 1, lower.tail = FALSE, log.p = TRUE)
  expect_gt(abs(z), 38.6)                    # past the underflow threshold
  expect_lt(lp, log(.Machine$double.xmin))   # the score p really is unrepresentable
  conv <- TRUE
  spa <- SPA_pval(mu, g, q, -1 * abs(q - m1) + m1, lp, .Machine$double.eps^0.25,
                  TRUE, "binary", conv)
  expect_true(conv)                          # ... and the saddlepoint converges
  expect_lt(spa, log(.Machine$double.xmin))  # ... on a log p of about -827
  # so both conditions the defective line needs are satisfiable at once.
  expect_gt(testQnormTwoSidedLogPInCPP(spa), 40)
})
