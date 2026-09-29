# LEDGER #172 (2026-09-29): the atom-dominance measure that routes a weighted
# set to the exact test. Kish's effective carrier count (sum v)^2 / sum v^2 over
# the carriers' values, and the number of MATERIAL carriers -- those at or above
# eps x the largest value, eps being the exact test's tie epsilon (1e-6 absolute
# after its scaling to a largest entry of 1). Closed forms pinned here; the
# routing itself is drilled end to end in test-step2-exact-by-weight.R.
eps <- 1e-6

test_that("equal weights: the effective count is the carrier count and every carrier is material", {
  for (n in c(1L, 4L, 23L)) {
    r <- atomDominanceInCPP(rep(1, n), eps)
    expect_equal(r$n_eff, n)
    expect_equal(r$n_material, n)
  }
  r <- atomDominanceInCPP(rep(0.37, 8), eps)          # scale-free
  expect_equal(r$n_eff, 8)
})

test_that("one atom among negligible carriers: effective count ~1, only the atom material", {
  v <- c(1, 10^seq(-19, -7, length.out = 22))         # the FlexRV sbeta20 shape (synthetic)
  r <- atomDominanceInCPP(v, eps)
  expect_lt(r$n_eff, 1.001)
  expect_equal(r$n_material, 1)
})

test_that("the material count follows the tie epsilon exactly", {
  v <- c(1, 2e-6, 1e-6, 9.9e-7, 1e-9)
  r <- atomDominanceInCPP(v, eps)
  expect_equal(r$n_material, 3)                       # 1, 2e-6 and 1e-6 (at the epsilon) count
  expect_equal(atomDominanceInCPP(v * 1e-3, eps)$n_material, 3)   # relative to the largest
})

test_that("a dominant atom with real companions: Kish counts the companions fractionally", {
  v <- c(1, 0.135, 0.014, 1e-3, 1e-6, 1e-9)   # sbeta20 on scores 1, .9, .8, .7 ... (synthetic)
  r <- atomDominanceInCPP(v, eps)
  expect_equal(r$n_eff, sum(v)^2 / sum(v^2))
  expect_lt(r$n_eff, 1.5)
  expect_equal(r$n_material, 5)
})

test_that("homozygote dosages: two doses of one weight are one atom of twice the value", {
  r <- atomDominanceInCPP(c(2, 1, 1), eps)            # one homozygote and two heterozygotes
  expect_equal(r$n_eff, 16 / 6)
  expect_equal(r$n_material, 3)
})

test_that("an all-zero vector has no atoms", {
  r <- atomDominanceInCPP(c(0, 0), eps)
  expect_equal(r$n_eff, 0)
  expect_equal(r$n_material, 0)
})
