# FlexRV weight kernels (src/FlexRV.cpp, exported as flexrvWeightsInCPP) and
# the R-side transform table (R/FlexRV.R), pinned to the closed forms of
# docs/design/flexrv-design.md. Two routes to every weight -- the C++ kernels
# on one side, the paper's dbeta ratios and cubic-root formula written out in
# R on the other -- must agree at 1e-12 on a grid that includes every boundary
# the design names: s = 0, s = 1, s = t, the MAF = maxMAF boundary from both
# sides, MAC = t and MAC = t + 1, and the 16 x 12 product ordering
# k = c*12 + d. Fast tier: pure kernels, like test-unit-spa.R.

flexrv_grid <- function() {
  s <- c(0, 0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
  q <- length(s)
  list(
    score = s,
    isLof = rep(c(0, 1), length.out = q),
    MAF = c(0, 5e-4, 0.001 - 1e-12, 0.001, 0.001 + 1e-12, 0.002, 0.004, 0.0001,
            0.00025, 0.0009, 0.00099, 0.0005),
    MAC = c(1, 3, 5, 10, 25, 50, 100, 200, 500, 2, 11, 501)
  )
}

test_that("the default specs give the paper's 16 x 12 transforms in score-major order", {
  tt <- flexrvTransformTable(flexrvDefaultScoreSpec, flexrvDefaultMAFSpec, "PAI3D")
  expect_identical(tt$nScore, 16L)
  expect_identical(tt$nMAF, 12L)
  expect_identical(dim(tt$table), c(192L, 4L))
  expect_identical(tt$scoreSuffix,
                   c("unw", "lof", paste0("cbrt", c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)),
                     paste0("sbeta", c(1.5, 2, 3, 5, 10, 20))))
  expect_identical(tt$mafSuffix,
                   c("unw", "mafbeta2", "mafbeta5",
                     paste0("MAC", c(1, 3, 5, 10, 25, 50, 100, 200, 500))))
  # k = c*12 + d (0-based): the MAF transform is the fastest-varying index
  for (k in c(0, 1, 11, 12, 13, 100, 191)) {
    cc <- k %/% 12
    d <- k %% 12
    expect_identical(tt$names[k + 1],
                     paste0("PAI3D_", tt$scoreSuffix[cc + 1], "_", tt$mafSuffix[d + 1]))
    kindName <- c("unw", "lof", "cbrt", "sbeta")[1 + (cc >= 1) + (cc >= 2) + (cc >= 10)]
    expect_identical(tt$table[k + 1, 1], as.numeric(flexrvScoreKinds[[kindName]]))
  }
  expect_identical(tt$names[1], "PAI3D_unw_unw")
  expect_identical(tt$names[192], "PAI3D_sbeta20_MAC500")
  expect_false(anyDuplicated(tt$names) > 0)
  expect_true(tt$usesLof)
})

test_that("spec parsing rejects what the design says it must", {
  expect_error(parseFlexRVTransformSpec("unw;lof", "maf"), "unknown transform 'lof'")
  expect_error(parseFlexRVTransformSpec("cbrt:0,0.5", "score"), "inside \\(0, 1\\)")
  expect_error(parseFlexRVTransformSpec("sbeta:1", "score"), "b > 1")
  expect_error(parseFlexRVTransformSpec("mafbeta:0.5", "maf"), "b > 1")
  expect_error(parseFlexRVTransformSpec("mac:0", "maf"), ">= 1")
  expect_error(parseFlexRVTransformSpec("mac:1,1", "maf"), "duplicated")
  expect_error(parseFlexRVTransformSpec("unw:1", "score"), "takes no parameters")
  expect_error(parseFlexRVTransformSpec("cbrt", "score"), "parameter list")
  expect_error(parseFlexRVTransformSpec("cbrt:a", "score"), "non-numeric")
  expect_error(parseFlexRVTransformSpec("", "score"), "non-empty")
  expect_error(parseFlexRVTransformSpec("unw;unw", "score"), "listed twice")
  expect_error(parseFlexRVTransformSpec("foo", "maf"), "unknown transform 'foo'")
  # a reduced spec is fine and keeps spec order
  p <- parseFlexRVTransformSpec("sbeta:5;unw", "score")
  expect_identical(p$suffix, c("sbeta5", "unw"))
  expect_identical(p$kind, c(3L, 0L))
  expect_error(flexrvTransformTable("unw", "unw", "S", maxMAF = 0), "flexRV_maxMAF")
  expect_error(flexrvTransformTable("unw", "unw", "S", maxMAF = 0.6), "flexRV_maxMAF")
})

test_that("C++ kernels and the R closed forms agree on the boundary grid (default table)", {
  g <- flexrv_grid()
  tt <- flexrvTransformTable(flexrvDefaultScoreSpec, flexrvDefaultMAFSpec, "S")
  wC <- flexrvWeightsInCPP(tt$table, g$score, g$isLof, g$MAF, g$MAC, 0.001)
  wR <- flexrvWeightsR(tt$table, g$score, g$isLof, g$MAF, g$MAC, 0.001)
  expect_identical(dim(wC), c(12L, 192L))
  expect_true(all(is.finite(wC)))
  expect_equal(wC, wR, tolerance = 1e-12)
  # positive control on the comparison itself: a perturbed table disagrees
  bad <- tt$table
  bad[25, 2] <- 0.25   # row k = 24 is cbrt0.2 x unw: t = 0.25 in the C++ table only
  expect_false(isTRUE(all.equal(flexrvWeightsInCPP(bad, g$score, g$isLof, g$MAF, g$MAC, 0.001),
                                wR, tolerance = 1e-12)))
})

test_that("score transforms: pinned values at s = 0, s = 1 and s = t", {
  tt <- flexrvTransformTable(flexrvDefaultScoreSpec, "unw", "S")
  s <- c(0, 1, 0.2, 0.5, 0.9)
  w <- flexrvWeightsInCPP(tt$table, s, rep(0, 5), rep(0, 5), rep(1, 5), 0.001)
  colnames(w) <- tt$scoreSuffix
  # every nonlinear transform maps 0 -> 0 and 1 -> 1
  nonlin <- grep("^cbrt|^sbeta", tt$scoreSuffix, value = TRUE)
  expect_identical(unname(w[1, nonlin]), rep(0, length(nonlin)))
  expect_equal(unname(w[2, nonlin]), rep(1, length(nonlin)), tolerance = 1e-15)
  # w(t) = t^(1/3) / ((1-t)^(1/3) + t^(1/3)): the design's 0.3865 / 0.5 / 0.6753
  expect_equal(w[[3, "cbrt0.2"]], 0.2^(1 / 3) / (0.8^(1 / 3) + 0.2^(1 / 3)), tolerance = 1e-12)
  expect_equal(round(w[[3, "cbrt0.2"]], 4), 0.3865)
  expect_equal(w[[4, "cbrt0.5"]], 0.5, tolerance = 1e-15)
  expect_equal(round(w[[5, "cbrt0.9"]], 4), 0.6753)
  # sbeta2 is the identity, which is why the paper's b = 2 panel is a line
  expect_equal(unname(w[, "sbeta2"]), s, tolerance = 1e-15)
  expect_equal(w[[4, "sbeta20"]], 0.5^19, tolerance = 1e-15)
  # unw is 1 everywhere; lof reads the indicator, not the score
  expect_identical(unname(w[, "unw"]), rep(1, 5))
  expect_identical(unname(w[, "lof"]), rep(0, 5))
  wl <- flexrvWeightsInCPP(tt$table, s, rep(1, 5), rep(0, 5), rep(1, 5), 0.001)
  expect_identical(unname(wl[, 2]), rep(1, 5))
  expect_identical(unname(wl[, 1]), rep(1, 5))
  # monotone increasing in s for every nonlinear transform
  fine <- seq(0, 1, by = 0.01)
  wf <- flexrvWeightsInCPP(tt$table, fine, rep(0, 101), rep(0, 101), rep(1, 101), 0.001)
  for (k in which(tt$scoreSuffix %in% nonlin)) {
    expect_true(all(diff(wf[, k]) >= 0), label = tt$scoreSuffix[k])
    expect_true(all(wf[, k] >= 0 & wf[, k] <= 1), label = tt$scoreSuffix[k])
  }
})

test_that("MAF transforms: the beta direction, the boundary from both sides, the MAC ladder", {
  tt <- flexrvTransformTable("unw", flexrvDefaultMAFSpec, "S")
  eps <- 1e-12
  MAF <- c(0, 1e-4, 5e-4, 0.001 - eps, 0.001, 0.001 + eps, 0.002)
  q <- length(MAF)
  w <- flexrvWeightsInCPP(tt$table, rep(0.5, q), rep(0, q), MAF, rep(1, q), 0.001)
  colnames(w) <- tt$mafSuffix
  # decreasing in MAF (A1): 1 at MAF 0, 0 at and above the boundary
  expect_identical(w[[1, "mafbeta2"]], 1)
  expect_identical(w[[1, "mafbeta5"]], 1)
  expect_true(all(diff(w[1:4, "mafbeta2"]) < 0))
  expect_true(all(diff(w[1:4, "mafbeta5"]) < 0))
  expect_identical(unname(w[5:7, "mafbeta2"]), c(0, 0, 0))
  expect_identical(unname(w[5:7, "mafbeta5"]), c(0, 0, 0))
  # the design's figure check: MAF/0.001 = 0.1 gives 0.9 (b=2) and 0.6561 (b=5)
  expect_equal(w[[2, "mafbeta2"]], 0.9, tolerance = 1e-12)
  expect_equal(w[[2, "mafbeta5"]], 0.9^4, tolerance = 1e-12)
  # one ulp below the boundary is still positive, the boundary itself is 0
  expect_true(w[[4, "mafbeta2"]] > 0)
  expect_identical(w[[5, "mafbeta2"]], 0)
  # the MAC ladder: 1 at MAC = t, 0 at MAC = t + 1 (and t + 0.5: MAC is a double)
  thr <- c(1, 3, 5, 10, 25, 50, 100, 200, 500)
  MAC <- c(thr, thr + 1, thr + 0.5)
  qm <- length(MAC)
  wm <- flexrvWeightsInCPP(tt$table, rep(0.5, qm), rep(0, qm), rep(1e-4, qm), MAC, 0.001)
  colnames(wm) <- tt$mafSuffix
  for (j in seq_along(thr)) {
    col <- paste0("MAC", thr[j])
    expect_identical(wm[[j, col]], 1)
    expect_identical(wm[[length(thr) + j, col]], 0)
    expect_identical(wm[[2 * length(thr) + j, col]], 0)
    # a MAC below every threshold up to and including this one qualifies
    expect_identical(unname(wm[[1, col]]), 1)
  }
  expect_identical(unname(wm[, "unw"]), rep(1, qm))
})

test_that("maxMAF != 0.001 rescales and renames the beta transforms; the MAC ladder is untouched", {
  tt <- flexrvTransformTable("unw", flexrvDefaultMAFSpec, "S", maxMAF = 0.005)
  expect_identical(tt$mafSuffix[2:3], c("mafbeta2@0.005", "mafbeta5@0.005"))
  expect_identical(tt$mafSuffix[4], "MAC1")
  expect_identical(tt$names[2], "S_unw_mafbeta2@0.005")
  MAF <- c(0.001, 0.0025, 0.005)
  MAC <- c(3, 10, 26)
  w <- flexrvWeightsInCPP(tt$table, rep(0.5, 3), rep(0, 3), MAF, MAC, 0.005)
  colnames(w) <- tt$mafSuffix
  expect_equal(unname(w[, "mafbeta2@0.005"]), c(0.8, 0.5, 0), tolerance = 1e-12)
  expect_equal(unname(w[, "mafbeta5@0.005"]), c(0.8^4, 0.5^4, 0), tolerance = 1e-12)
  expect_identical(unname(w[, "MAC3"]), c(1, 0, 0))
  expect_identical(unname(w[, "MAC25"]), c(1, 1, 0))
  # the R route rescales the same way
  expect_equal(unname(w), flexrvWeightsR(tt$table, rep(0.5, 3), rep(0, 3), MAF, MAC, 0.005),
               tolerance = 1e-12)
  # the default maxMAF keeps the paper's bare names
  expect_identical(flexrvTransformTable("unw", "mafbeta:2", "S", maxMAF = 0.001)$names,
                   "S_unw_mafbeta2")
})

test_that("the weight is the product of the two transforms; only unw_unw is all ones", {
  g <- flexrv_grid()
  tt <- flexrvTransformTable(flexrvDefaultScoreSpec, flexrvDefaultMAFSpec, "S")
  ts <- flexrvTransformTable(flexrvDefaultScoreSpec, "unw", "S")
  tm <- flexrvTransformTable("unw", flexrvDefaultMAFSpec, "S")
  w <- flexrvWeightsInCPP(tt$table, g$score, g$isLof, g$MAF, g$MAC, 0.001)
  ws <- flexrvWeightsInCPP(ts$table, g$score, g$isLof, g$MAF, g$MAC, 0.001)
  wm <- flexrvWeightsInCPP(tm$table, g$score, g$isLof, g$MAF, g$MAC, 0.001)
  for (k in seq_len(192)) {
    cc <- (k - 1) %/% 12 + 1
    d <- (k - 1) %% 12 + 1
    expect_identical(w[, k], ws[, cc] * wm[, d])
  }
  allOnes <- which(apply(w, 2, function(col) all(col == 1)))
  expect_identical(tt$names[allOnes], "S_unw_unw")
})

test_that("the C++ export validates its inputs", {
  bad <- matrix(c(9, 0, 0, 0), nrow = 1)
  expect_error(flexrvWeightsInCPP(bad, 0.5, 0, 1e-4, 1, 0.001), "unknown score transform kind 9")
  bad2 <- matrix(c(0, 0, 7, 0), nrow = 1)
  expect_error(flexrvWeightsInCPP(bad2, 0.5, 0, 1e-4, 1, 0.001), "unknown MAF transform kind 7")
  expect_error(flexrvWeightsInCPP(matrix(0, 1, 3), 0.5, 0, 1e-4, 1, 0.001), "4 columns")
  ok <- matrix(c(0, 0, 0, 0), nrow = 1)
  expect_error(flexrvWeightsInCPP(ok, c(0.5, 0.5), 0, 1e-4, 1, 0.001), "one entry per marker")
  expect_error(flexrvWeightsInCPP(ok, 0.5, 0, 1e-4, 1, 0), "maxMAF must be > 0")
  # a 0-row table is legal (no FlexRV sets) and gives a q x 0 matrix
  expect_identical(dim(flexrvWeightsInCPP(matrix(0, 0, 4), 0.5, 0, 1e-4, 1, 0.001)), c(1L, 0L))
})

test_that("parseFlexRVLofAnno treats ',', ':' and ';' alike", {
  expect_identical(parseFlexRVLofAnno("lof"), "lof")
  expect_identical(parseFlexRVLofAnno("pLoF,LC:ptv;lof"), c("pLoF", "LC", "ptv", "lof"))
  expect_identical(parseFlexRVLofAnno("lof,lof"), "lof")
  expect_identical(parseFlexRVLofAnno(NULL), character(0))
  expect_error(parseFlexRVLofAnno(",;"), "at least one")
})
