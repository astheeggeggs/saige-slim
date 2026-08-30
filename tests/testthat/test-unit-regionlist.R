# SAIGE.getRegionList_new K-weight-set generalization
# (R/SAIGE_SPATest_Region.R). Weight columns ride inside RegionData through
# merge/filter/sort under sanitized internal names (SAIGEWT_k), so user-facing
# set names can never collide with markerInfo columns. Per-region return:
# K=0 keeps the c(0) sentinel; K>=1 is a numeric matrix (nrow = surviving
# markers) with the user-facing set names as colnames.

# Synthetic fixtures: POS deliberately out of group-file order so the
# setorderv(CHROM, POS) realignment is exercised; genoIndex distinct from POS
# so any weight/markerInfo column collision would be visible. No ID2 column.
mi_fake <- function() {
  data.table::data.table(
    ID = c("rsA", "rsB", "rsC"),
    CHROM = c("1", "1", "1"),
    POS = c(300, 100, 200),
    genoIndex = c(7, 3, 5)
  )
}

var_line  <- "G1 var rsA rsB rsC"
anno_line <- "G1 anno lof missense lof"

test_that("K=0 (2-line) matches the legacy signature; WEIGHT is exactly c(0)", {
  lines <- c(var_line, anno_line)
  legacy <- SAIGE.getRegionList_new(lines, 2, c("lof", "lof;missense"),
                                    mi_fake(), "")
  explicit <- SAIGE.getRegionList_new(lines, 2, c("lof", "lof;missense"),
                                      mi_fake(), "",
                                      nWeightSets = 0L,
                                      weightSetNames = character(0))
  expect_identical(explicit, legacy)
  expect_identical(legacy$G1$SNP, c("rsB", "rsC", "rsA"))  # POS order
  expect_identical(legacy$G1$WEIGHT, c(0))
})

test_that("K=1 legacy bare weight line: 1-col matrix, colname WEIGHT, legacy values", {
  lines <- c(var_line, anno_line, "G1 weight 0.9 0.8 0.7")

  # no anno filtering: all three markers survive, POS-sorted rsB, rsC, rsA
  rl <- SAIGE.getRegionList_new(lines, 3, c("lof", "lof;missense"),
                                mi_fake(), "")
  w <- rl$G1$WEIGHT
  expect_true(is.matrix(w) && is.numeric(w))
  expect_identical(dim(w), c(3L, 1L))
  expect_identical(colnames(w), "WEIGHT")
  # legacy path result: as.numeric of the same strings, CHROM/POS-sorted
  expect_identical(unname(w[, 1]), as.numeric(c("0.8", "0.7", "0.9")))

  # explicit new-style call is identical to the inferred legacy call
  explicit <- SAIGE.getRegionList_new(lines, 3, c("lof", "lof;missense"),
                                      mi_fake(), "",
                                      nWeightSets = 1L,
                                      weightSetNames = "WEIGHT")
  expect_identical(explicit, rl)

  # anno filtering drops rsB (missense); its weight goes with it
  rlf <- SAIGE.getRegionList_new(lines, 3, c("lof"), mi_fake(), "")
  expect_identical(rlf$G1$SNP, c("rsC", "rsA"))
  expect_identical(dim(rlf$G1$WEIGHT), c(2L, 1L))
  expect_identical(unname(rlf$G1$WEIGHT[, 1]), as.numeric(c("0.7", "0.9")))
})

test_that("K=2 with adversarial set names survives merge, marker drop and anno filter", {
  lines <- c("G1 var rsA rsB rsC rsX",
             "G1 anno lof missense lof lof",
             "G1 weight:POS 10 20 30 40",
             "G1 weight,genoIndex 0.1 0.2 0.3 0.4")
  # rsX is absent from markerInfo (dropped in the merge, with a message);
  # rsB is dropped by the anno filter (annoVec = "lof" only).
  # Survivors in POS order: rsC (200), rsA (300).
  out <- capture.output(
    rl <- SAIGE.getRegionList_new(lines, 4, c("lof"), mi_fake(), "",
                                  nWeightSets = 2L,
                                  weightSetNames = c("POS", "genoIndex"))
  )
  expect_true(any(grepl("not in 'GenoFile'", out)))
  g1 <- rl$G1
  expect_identical(g1$SNP, c("rsC", "rsA"))
  w <- g1$WEIGHT
  expect_identical(dim(w), c(2L, 2L))
  expect_identical(colnames(w), c("POS", "genoIndex"))
  # per-marker values from the input lines: rsC -> (30, 0.3), rsA -> (10, 0.1)
  expect_identical(unname(w[, "POS"]), as.numeric(c("30", "10")))
  expect_identical(unname(w[, "genoIndex"]), as.numeric(c("0.3", "0.1")))
  # markerInfo-derived fields are untouched by the would-be name collision
  expect_identical(g1$genoIndex, c("5", "7"))
})

test_that("markerInfo=NULL + chrom filter path works with K=2", {
  lines <- c("G1 var 1:100:A:C 1:200:G:T 2:300:A:G",
             "G1 anno lof lof missense",
             "G1 weight:W1 0.1 0.2 0.3",
             "G1 weight:W2 1 2 3",
             "G2 var 2:400:A:C 2:500:G:T",
             "G2 anno lof lof",
             "G2 weight:W1 0.4 0.5",
             "G2 weight:W2 4 5")
  rl <- SAIGE.getRegionList_new(lines, 4, c("lof", "lof;missense"), NULL, "1",
                                nWeightSets = 2L,
                                weightSetNames = c("W1", "W2"))
  # Pre-existing markerInfo=NULL contract (preserved): a region whose markers
  # all fall on another chromosome returns list(SNP = NULL); a region that
  # keeps markers gets NO RegionList entry (the list(...) assignment is
  # guarded by !is.null(markerInfo)). The weight machinery on this path is
  # therefore observable only as "no error" plus structural identity with
  # the weightless call on the same var/anno lines.
  expect_identical(rl$G2, list(SNP = NULL))
  expect_null(rl$G1)
  legacy <- SAIGE.getRegionList_new(lines[c(1, 2, 5, 6)], 2,
                                    c("lof", "lof;missense"), NULL, "1")
  expect_identical(rl, legacy)
})

test_that("layout-contract violations are rejected", {
  lines <- c(var_line, anno_line, "G1 weight 0.9 0.8 0.7")
  expect_error(
    SAIGE.getRegionList_new(lines, 3, c("lof"), mi_fake(), "",
                            nWeightSets = 0L, weightSetNames = character(0)),
    "must equal 2 \\+ nWeightSets")
  expect_error(
    SAIGE.getRegionList_new(c(lines, "G1 weight:B 1 2 3"), 4, c("lof"),
                            mi_fake(), ""),
    "nWeightSets must be provided")
  expect_error(
    SAIGE.getRegionList_new(c(lines, "G1 weight:B 1 2 3"), 4, c("lof"),
                            mi_fake(), "",
                            nWeightSets = 2L, weightSetNames = "onlyone"),
    "one name per weight set")
})
