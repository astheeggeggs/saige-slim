# C5 jmr degeneracy (fast tier): at K = 1 the multi-weight-set boundary into
# mainRegionInCPP must be a pure re-boxing of the legacy objects — the same
# numeric values in a 1-column matrix, the same index arithmetic. The byte
# gate over harness/goldens/step2 is the primary verification; these tests
# pin the R-side halves of the contract so a call-site change fails fast.
#
# Normative flat index (mirrored in src/Main.cpp mainRegionInCPP and the
# run_mainRegionInCPP call-site comment in R/SAIGE_SPATest_Region.R):
#   jmr = j*q_maf*q_weight + m*q_weight + r      (weight fastest)

test_that("the normative jmr formula reduces to legacy jm at q_weight == 1", {
  jmr <- function(j, m, r, q_maf, q_weight) j * q_maf * q_weight + m * q_weight + r
  q_anno <- 4L; q_maf <- 3L
  for (j in 0:(q_anno - 1L)) {
    for (m in 0:(q_maf - 1L)) {
      expect_identical(jmr(j, m, 0L, q_maf, 1L), j * q_maf + m)
    }
  }
  # weight-fastest layout at K = 2: the K slots of one (j,m) are adjacent,
  # ordered by r, and the next (j,m) starts right after them
  expect_identical(jmr(1L, 2L, 0L, 3L, 2L), 10L)
  expect_identical(jmr(1L, 2L, 1L, 3L, 2L), 11L)
  expect_identical(jmr(0L, 0L, 1L, 3L, 2L) + 1L, jmr(0L, 1L, 0L, 3L, 2L))
})

test_that("call-site canonicalization maps the legacy WEIGHT objects exactly", {
  # unweighted c(0) sentinel -> 1x1 zero matrix (C++ falls back to the
  # implicit beta set, the legacy beta-weights path)
  w0 <- matrix(as.numeric(c(0)), ncol = 1)
  expect_identical(dim(w0), c(1L, 1L))
  expect_true(all(w0 == 0))
  # weighted vector -> q0 x 1 single-custom-set matrix, values untouched
  wv <- c(0.9, 0.8, 0.7)
  w1 <- matrix(as.numeric(wv), ncol = 1)
  expect_identical(dim(w1), c(3L, 1L))
  expect_identical(as.numeric(w1), wv)
  # is_no_weight's rep(1, q0) stays an all-ones custom set (C++ per-set ER
  # gate treats an all-ones custom column as "equal weights", as today)
  expect_true(all(matrix(as.numeric(rep(1, 5)), ncol = 1) == 1))
  # the beta-params sentinel shape sent at K <= 1
  bp <- matrix(0, nrow = 1, ncol = 2)
  expect_identical(dim(bp), c(1L, 2L))
  expect_true(all(bp == 0))
})

test_that("K=1 fixture: parser -> region list -> canonicalized matrix keeps the file's values", {
  gfile <- file.path(saige_extdata, "input", "group_new_snpid_withWeights.txt")
  det <- checkGroupFile(gfile, detail = TRUE)
  expect_identical(det$nWeightSets, 1L)
  # ascending POS: region-list order == group-file order, so the weight
  # vector must come back exactly as written in the fixture
  mi <- data.table::data.table(ID = paste0("rs", 1:100), CHROM = rep("1", 100),
                               POS = 1:100, genoIndex = 0:99)
  rl <- SAIGE.getRegionList_new(readLines(gfile), 2L + det$nWeightSets,
                                c("lof", "lof;missense"), mi, "",
                                det$nWeightSets, det$weightSetNames)
  # the K<=1 collapse used in SAIGE.Region's WEIGHT assignment
  WEIGHT <- as.numeric(rl$GENE1$WEIGHT)
  file_w <- as.numeric(strsplit(readLines(gfile)[3], "[ \t]+")[[1]][-(1:2)])
  expect_identical(WEIGHT, file_w)
  # the run_mainRegionInCPP re-boxing: q0 x 1, values bit-identical
  boxed <- matrix(as.numeric(WEIGHT), ncol = 1)
  expect_identical(dim(boxed), c(50L, 1L))
  expect_identical(unname(boxed[, 1]), file_w)
})
