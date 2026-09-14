# LEDGER #143 (units half FIXED 2026-09-14): the binary conditional region
# path fed the burden scale-factor kernel an UNWEIGHTED conditioning score
# against a WEIGHTED Phi, q.sum and g.sum (RegionSetUpConditional_binary_InCPP
# exported m_Tstat_cond raw; it now exports m_Tstat_cond % w0).
#
# The observable consequence was an invariance that failed. Conditioning on a
# set of markers is a linear reparameterisation, so the conditional columns
# cannot legitimately depend on the weight attached to a conditioning marker:
# for ONE conditioning marker any non-zero weight must give the same answer,
# and for several, a UNIFORM rescaling of all their weights must. Before the
# fix they did not: on P3's configuration conditioning on rs30 alone, weight 1
# against weight 5 moved every conditional p-value on GENE1 by 0.3-6 %
# (measured 2026-09-14 at e508a8e4, numbers on the ledger row), and
# SE_Burden_cond by ~2e-5 on the rows with real rare markers (r cancels out
# of VarMatAdjCond, the per-marker unit factor did not). After the fix the
# two runs agree to 1e-15 on every conditional column.
#
# Two layers, so a failure is LOCATED rather than observed:
#   1. fast tier: BOTH kernels are invariant under a CONSISTENT rescaling of
#      their four inputs (q.sum, g.sum, Score by w; Phi by w^2). So the kernel
#      is not where a weight dependence could come from. Exactly (1e-15) below
#      the |z| < 2 gate; above it the C++ saddlepoint root-finder's ABSOLUTE
#      tolerance (eps^0.25 on the root, which scales as 1/w) puts a x7 rescale
#      2.4e-4 off, so that regime is held at 1e-3. That is a property of
#      SPA_pval, not of these kernels, and is not #143.
#   2. standard tier: the CLI IS invariant (the inverted characterization).
#      A regression here with 1 still green means the C++ export or the R
#      consumer has put one input back into different units.
#
# Aliased at top level so the linter can see the bindings testthat supplies
# from helper-saige.R / helper-multiweight.R at run time.
l143_run <- function(...) get("mw_run")("group_new_snpid.txt", c(...))
l143_pnum <- function(x) suppressWarnings(as.numeric(x))

L143_UNCOND <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden",
                 "MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")
L143_COND <- c("Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond",
               "BETA_Burden_cond", "SE_Burden_cond")

# max relative difference over the non-Cauchy rows of a column
l143_maxrel <- function(a, b, col) {
  x <- l143_pnum(a[[col]][a$Group != "Cauchy"])
  y <- l143_pnum(b[[col]][b$Group != "Cauchy"])
  ok <- !is.na(x) & !is.na(y)
  max(abs(x[ok] - y[ok]) / pmax(abs(x[ok]), abs(y[ok]), .Machine$double.xmin))
}

# ---- 1. the kernels are invariant under consistent rescaling (fast) ---------

# A small binary null: n samples, q conditioning markers, g.sum residualised
# so that (q.sum - m1)^2 == sum(Score)^2 as on the unconditional path.
l143_inputs <- function(seed, z_target) {
  set.seed(seed)
  n <- 400L
  q <- 2L
  mu <- stats::plogis(stats::rnorm(n, -1.5, 0.7))
  g <- matrix(stats::rbinom(n * q, 2, 0.05), n, q)
  w <- c(3, 11)
  gw <- g %*% w
  g.sum <- as.vector(gw - mean(gw))
  var1 <- sum(mu * (1 - mu) * g.sum^2)
  m1 <- sum(mu * g.sum)
  # a score vector whose sum lands the standardised burden at z_target
  s_raw <- stats::rnorm(q)
  Score <- s_raw * (z_target * sqrt(var1) / sum(s_raw))
  q.sum <- m1 + sum(Score)
  # a PD Phi in the same (weighted) units as Score, and per-marker log p
  A <- matrix(stats::rnorm(q * q), q, q)
  Phi <- crossprod(A) + diag(q)
  Phi <- Phi * (var1 / sum(Phi)) * 0.9
  p.new <- stats::pchisq(Score^2 / diag(Phi), df = 1, lower.tail = FALSE, log.p = TRUE)
  list(q.sum = q.sum, mu.a = mu, g.sum = g.sum, p.new = p.new, Score = Score, Phi = Phi)
}

l143_scale <- function(inp, w) {
  m1 <- sum(inp$mu.a * inp$g.sum)
  inp$g.sum <- inp$g.sum * w
  inp$q.sum <- m1 * w + (inp$q.sum - m1) * w
  inp$Score <- inp$Score * w
  inp$Phi <- inp$Phi * w^2
  inp
}

l143_kernel <- function(fn, inp, ...) {
  do.call(fn, c(unname(inp[c("q.sum", "mu.a", "g.sum", "p.new", "Score", "Phi")]),
                list("SKAT-O", ...)))
}

test_that("both scale-factor kernels are invariant under a consistent rescaling of their inputs", {
  ns <- asNamespace("SAIGE")
  twin <- get("get_newPhi_scaleFactor", ns)
  fixed <- get("get_newPhi_scaleFactor_traitType", ns)
  for (z in c(0.8, 3.5)) {          # below and above the |z| < 2 saddlepoint gate
    inp <- l143_inputs(143L, z)
    for (w in c(0.2, 7)) {
      sc <- l143_scale(inp, w)
      lab <- sprintf("z=%g w=%g", z, w)
      tol <- if (z < 2) 1e-12 else 1e-3
      a <- l143_kernel(twin, inp)
      b <- l143_kernel(twin, sc)
      expect_equal(b$scaleFactor, a$scaleFactor, tolerance = tol, label = paste("twin", lab))
      expect_equal(b$val / w^2, a$val, tolerance = tol, label = paste("twin val", lab))
      a <- l143_kernel(fixed, inp, "binary")
      b <- l143_kernel(fixed, sc, "binary")
      expect_equal(b$scaleFactor, a$scaleFactor, tolerance = tol, label = paste("fixed", lab))
      expect_equal(b$val / w^2, a$val, tolerance = tol, label = paste("fixed val", lab))
    }
  }
  # and the inputs really straddle the gate: the two z regimes give different r
  inp_lo <- l143_inputs(143L, 0.8)
  inp_hi <- l143_inputs(143L, 3.5)
  expect_false(identical(l143_kernel(fixed, inp_lo, "binary")$scaleFactor,
                         l143_kernel(fixed, inp_hi, "binary")$scaleFactor))
})

# ---- 2. the CLI is not (standard) --------------------------------------------

test_that("conditioning on ONE marker: _cond does not depend on its weight (LEDGER #143)", {
  a <- l143_run("--condition=rs30", "--weights_for_condition=1")$region
  b <- l143_run("--condition=rs30", "--weights_for_condition=5")$region
  expect_identical(a$Region, b$Region)
  expect_identical(a$Group, b$Group)
  # the unconditional columns do not see the conditioning weight at all
  for (cc in L143_UNCOND) expect_identical(a[[cc]], b[[cc]], label = cc)
  # the inverted characterization: before the fix these moved by 1e-3..6e-2
  # (p-values, BETA) and ~2e-5 (SE); measured 1e-15 after it
  for (cc in L143_COND) expect_lt(l143_maxrel(a, b, cc), 1e-10, label = cc)
})

test_that("conditioning on TWO markers: rescaling both weights leaves _cond alone (LEDGER #143)", {
  a <- l143_run("--condition=rs30,rs79", "--weights_for_condition=1,1")$region
  b <- l143_run("--condition=rs30,rs79", "--weights_for_condition=5,5")$region
  expect_identical(a$Group, b$Group)
  for (cc in L143_UNCOND) expect_identical(a[[cc]], b[[cc]], label = cc)
  for (cc in L143_COND) expect_lt(l143_maxrel(a, b, cc), 1e-10, label = cc)
})
