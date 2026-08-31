# CCT (Cauchy combination test): reconcile the R implementation
# (R/CCT_modified.R, used for the per-region "Cauchy" rows) against the C++
# one (src/CCT.cpp CCT_cpp, used inside mainRegionInCPP) — the same quantity
# computed two ways must agree.
#
# WARNING: never pass NA or out-of-[0,1] values to CCT_cpp — it calls
# exit(EXIT_FAILURE) and kills the R process (untestable in-process).

test_that("R and C++ CCT agree on the equal-weights path", {
  cases <- list(
    c(2e-02, 4e-04, 0.2, 0.1, 0.8),  # reference: 0.0019534044057701871
    c(0.1, 0.1),
    0.5,
    1e-5,
    c(0.3, 1, 0.7),                  # p==1 present -> Bonferroni-on-min path
    c(1, 1),
    c(1e-20, 0.5),                   # tiny-p branch
    c(0, 0.5),                       # p==0 -> 0
    c(0, 1, 0.5)                     # 0 beats 1
  )
  # Tolerance 1e-10, not 1e-12: CCT round-trips through
  # 0.5 - atan(tan((0.5 - p) * pi)) / pi, and for small p the argument sits
  # next to pi/2 where tan is ill-conditioned, so the round trip loses
  # precision. Measured worst case across these cases is 1.55e-12, all of it
  # from the scalar p = 1e-5 (every vector case is <= 5.3e-15). This is
  # inherent conditioning, not a divergence between the two implementations;
  # 1e-10 keeps ~15x margin over the observed error. Edition 2 accepted the
  # old 1e-12 because its numeric comparison was looser than the number
  # written here; waldo enforces it literally.
  for (p in cases) {
    expect_equal(CCT(p), CCT_cpp(as.matrix(p)), tolerance = 1e-10,
                 label = paste0("CCT(", paste(p, collapse = ","), ")"))
  }
})

test_that("CCT reference values are stable", {
  expect_equal(CCT(c(2e-02, 4e-04, 0.2, 0.1, 0.8)), 0.0019534044057701871, tolerance = 1e-14)
  expect_equal(CCT(c(0.3, 1, 0.7)), 0.9, tolerance = 1e-14)          # min(1, 0.3*3)
  expect_equal(CCT(c(1e-20, 0.5)), 2e-20, tolerance = 1e-12)
  expect_identical(CCT(c(0, 0.5)), 0)
  expect_equal(CCT(c(0.01, 0.5), weights = c(3, 1)), 0.013329923861388204, tolerance = 1e-14)
})

test_that("scalar CCT is the identity on (1e-16, 1)", {
  # identity only to the tan/atan round-trip precision -- 6.55e-12 at p = 1e-5
  # (exact at 0.2, 0.5, 0.99). See the tolerance note above.
  for (p in c(1e-5, 0.2, 0.5, 0.99)) expect_equal(CCT(p), p, tolerance = 1e-10)
})

test_that("CCT weight handling", {
  p <- c(0.01, 0.5, 0.2)
  w <- c(3, 1, 2)
  expect_equal(CCT(p, w), CCT(p, 10 * w), tolerance = 1e-14)  # normalization invariance
  expect_error(CCT(p, c(1, 2)), "The length of weights should be the same")
  expect_error(CCT(p, c(-1, 1, 1)), "All the weights must be positive!")
})

test_that("CCT input validation (R side)", {
  expect_error(CCT(c(0.1, NA)), "Cannot have NAs in the p-values!")
  expect_error(CCT(c(0.1, 1.2)), "All p-values must be between 0 and 1!")
  expect_error(CCT(c(-0.1, 0.5)), "All p-values must be between 0 and 1!")
})
