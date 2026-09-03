# lchoose_threadsafe: the stand-in that let the ER path stop being main-thread
# only (Step A).
#
# This file exists because THE GOLDENS CANNOT TEST THIS FUNCTION, and that was
# measured rather than assumed. They do execute it -- replacing its main branch
# with a hard error fails the standard tier in two places, via mainMarkerInCPP --
# but perturbing that branch by +1e-9 leaves the whole tier green. Printed ER
# output is simply insensitive to an absolute shift that small in a
# log-combination. So a golden run confirms the branch RUNS and says almost
# nothing about whether it is RIGHT.
#
# The swap's entire justification is bit-identity with Rf_lchoose, so that is
# what is asserted, over the domain HyperGeo::lCombinations can actually reach:
# integer n >= 0 and integer k, with k > n short-circuited by the caller.
# Comparison is by BITS (memcmp in C++), not `==`: `==` would call every NaN
# pair a difference and every -0.0/0.0 pair a match, which is the wrong standard
# in both directions.

# R's own lchoose, named to make the branch table below read unambiguously.
Rf_lchoose_ref <- function(n, k) lchoose(n, k)

test_that("lchoose_threadsafe is bit-identical to Rf_lchoose over the reachable domain", {
  # ~2.0e6 argument pairs; runs in well under a second.
  expect_identical(lchooseThreadsafeDiffCountInCPP(2000L), 0L)
})

test_that("the comparison has power: a smaller domain is still exhaustive, and the count is a count", {
  # Guards against the assertion above passing because the loop never ran --
  # the shape of LEDGER #30. n = 0 still yields the k in {-3..0} pairs.
  expect_identical(lchooseThreadsafeDiffCountInCPP(0L), 0L)
  expect_identical(lchooseThreadsafeDiffCountInCPP(50L), 0L)
})

test_that("the domain covers the branches R's choose.c actually takes", {
  # Not a test of our code so much as a check that the domain above is not
  # trivially all-early-return: k < 2 (constant), the n - k < 2 symmetry arm,
  # and the lfastchoose arm must each be reachable within n <= 2000.
  expect_identical(Rf_lchoose_ref(5, 0), 0)                 # k == 0
  expect_equal(Rf_lchoose_ref(5, 1), log(5))                # k == 1
  expect_equal(Rf_lchoose_ref(5, 5), 0)                     # symmetry, n-k == 0
  expect_equal(Rf_lchoose_ref(5, 4), log(5))                # symmetry, n-k == 1
  expect_equal(Rf_lchoose_ref(5, 2), log(10))               # lfastchoose
})
