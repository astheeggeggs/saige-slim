# Committed weighted group-file fixtures (C4):
#   extdata/input/group_new_snpid_withWeights.txt  — legacy bare 3-line format
#     (one "weight" line per region -> single set named WEIGHT)
#   extdata/input/group_new_snpid_multiWeights.txt — two named sets per region,
#     "weight:CADD" then "weight,REVEL" (both delimiter forms)
# plus the unweighted control group_new_snpid.txt. Characterization tests:
# the expected values below are literal transcriptions of the committed
# fixture content, and the getRegionList checks verify that weight/marker
# alignment survives merge + POS sort + anno filtering under a deliberately
# scrambled markerInfo. Fast tier.

gf_unweighted   <- function() file.path(saige_extdata, "input", "group_new_snpid.txt")
gf_withweights  <- function() file.path(saige_extdata, "input", "group_new_snpid_withWeights.txt")
gf_multiweights <- function() file.path(saige_extdata, "input", "group_new_snpid_multiWeights.txt")

# Literal transcriptions of the committed weight lines (fixture position k
# belongs to rs<k> for GENE1 and rs<50+k> for GENE2). The multiWeights CADD
# lines carry the same values as the withWeights bare-weight lines, and the
# fixture reuses one 50-value REVEL sequence for both regions.
w_gene1 <- c(7.8, 3.9, 8.9, 3.8, 4.4, 1.3, 4.9, 6.6, 4.1, 2.3,
             0.8, 4, 1.4, 7.4, 7.3, 1.9, 3.2, 6.1, 1.1, 9.9,
             5.5, 5.2, 1.7, 4.3, 0.5, 9.7, 7.9, 0.7, 3.3, 4.8,
             8.3, 3.7, 9.6, 8, 6.3, 1.8, 7.6, 4.2, 3.1, 2.7,
             7.5, 8.5, 2.1, 8.4, 10, 5.6, 6.7, 6, 5.9, 5.3)
w_gene2 <- c(4.1, 7.7, 9.5, 8.4, 2.3, 3.5, 5.1, 6.6, 7.4, 3.1,
             10, 9.4, 8.2, 7.2, 5.5, 7.8, 7.1, 5.4, 1, 0.8,
             5.8, 4.8, 9.1, 5.7, 6.7, 2.1, 1.5, 8.7, 8.9, 8.8,
             7, 1.9, 1.8, 3.2, 7.3, 9.8, 8.6, 2.6, 1.1, 9,
             6.2, 6.4, 8.1, 4, 4.7, 5.3, 6.5, 6.8, 3.9, 1.7)
revel_vals <- c(0.01, 0.38, 0.75, 0.15, 0.52, 0.89, 0.29, 0.66, 0.06, 0.43,
                0.8, 0.2, 0.57, 0.94, 0.34, 0.71, 0.11, 0.48, 0.85, 0.25,
                0.62, 0.02, 0.39, 0.76, 0.16, 0.53, 0.9, 0.3, 0.67, 0.07,
                0.44, 0.81, 0.21, 0.58, 0.95, 0.35, 0.72, 0.12, 0.49, 0.86,
                0.26, 0.63, 0.03, 0.4, 0.77, 0.17, 0.54, 0.91, 0.31, 0.68)
anno_gene1 <- c(rep("lof", 10), rep("missense", 20), rep("lof", 20))
anno_gene2 <- c(rep("missense", 30), rep("lof", 20))

# Fake markerInfo covering all 100 fixture rsIDs. POS is reversed relative to
# group-file order (rs_k -> 1010 - 10*k) so the setorderv(CHROM, POS)
# realignment is exercised: per-region POS order is rs50..rs1 / rs100..rs51.
# genoIndex (1001..1100) is disjoint from POS so a weight/markerInfo column
# mix-up would be visible. No ID2 column.
mi_scrambled <- function() {
  data.table(
    ID = paste0("rs", 1:100),
    CHROM = rep("1", 100),
    POS = seq(1000, 10, by = -10),
    genoIndex = 1000 + 1:100
  )
}

test_that("checkGroupFile(detail=TRUE) characterizes all three committed fixtures", {
  expect_identical(checkGroupFile(gf_unweighted(), detail = TRUE),
                   list(nRegions = 2, is_weight_included = FALSE,
                        nWeightSets = 0L, weightSetNames = character(0)))
  expect_identical(checkGroupFile(gf_withweights(), detail = TRUE),
                   list(nRegions = 2, is_weight_included = TRUE,
                        nWeightSets = 1L, weightSetNames = "WEIGHT"))
  expect_identical(checkGroupFile(gf_multiweights(), detail = TRUE),
                   list(nRegions = 2, is_weight_included = TRUE,
                        nWeightSets = 2L, weightSetNames = c("CADD", "REVEL")))
})

test_that("fixture invariant: every weight line has one numeric value per var rsID", {
  for (path in c(gf_withweights(), gf_multiweights())) {
    toks <- strsplit(readLines(path), "[ \t]+")
    genes <- unique(vapply(toks, `[`, character(1), 1))
    expect_identical(genes, c("GENE1", "GENE2"), info = basename(path))
    for (g in genes) {
      gl <- Filter(function(x) x[1] == g, toks)
      types <- vapply(gl, `[`, character(1), 2)
      nvar <- length(gl[[which(types == "var")]]) - 2L
      expect_identical(length(gl[[which(types == "anno")]]) - 2L, nvar,
                       info = paste(basename(path), g, "anno count"))
      widx <- grep("^weight", types)
      expect_gte(length(widx), 1L)
      for (j in widx) {
        vals <- suppressWarnings(as.numeric(gl[[j]][-(1:2)]))
        expect_identical(length(vals), nvar,
                         info = paste(basename(path), g, types[j], "value count"))
        expect_false(anyNA(vals),
                     info = paste(basename(path), g, types[j], "non-numeric value"))
      }
    }
  }
})

test_that("withWeights fixture content matches its desk-derived transcription", {
  toks <- strsplit(readLines(gf_withweights()), "[ \t]+")
  expect_identical(length(toks), 6L)
  expect_identical(vapply(toks, `[`, character(1), 1),
                   rep(c("GENE1", "GENE2"), each = 3))
  expect_identical(vapply(toks, `[`, character(1), 2),
                   rep(c("var", "anno", "weight"), 2))
  expect_identical(toks[[1]][-(1:2)], paste0("rs", 1:50))
  expect_identical(toks[[2]][-(1:2)], anno_gene1)
  expect_identical(as.numeric(toks[[3]][-(1:2)]), w_gene1)
  expect_identical(toks[[4]][-(1:2)], paste0("rs", 51:100))
  expect_identical(toks[[5]][-(1:2)], anno_gene2)
  expect_identical(as.numeric(toks[[6]][-(1:2)]), w_gene2)
})

test_that("multiWeights fixture content matches its desk-derived transcription", {
  toks <- strsplit(readLines(gf_multiweights()), "[ \t]+")
  expect_identical(length(toks), 8L)
  expect_identical(vapply(toks, `[`, character(1), 2),
                   rep(c("var", "anno", "weight:CADD", "weight,REVEL"), 2))
  # var/anno lines are token-identical to the unweighted control file
  ctrl <- strsplit(readLines(gf_unweighted()), "[ \t]+")
  expect_identical(toks[c(1, 2, 5, 6)], ctrl)
  # CADD carries the same values as the withWeights bare-weight lines
  expect_identical(as.numeric(toks[[3]][-(1:2)]), w_gene1)
  expect_identical(as.numeric(toks[[7]][-(1:2)]), w_gene2)
  # the same REVEL sequence appears in both regions (committed fixture fact)
  expect_identical(as.numeric(toks[[4]][-(1:2)]), revel_vals)
  expect_identical(as.numeric(toks[[8]][-(1:2)]), revel_vals)
})

test_that("getRegionList_new on withWeights: 1-col WEIGHT stays aligned through the POS scramble", {
  lines <- readLines(gf_withweights())
  rl <- SAIGE.getRegionList_new(lines, 3, c("lof", "lof;missense"),
                                mi_scrambled(), "")
  expect_identical(names(rl), c("GENE1", "GENE2"))

  g1 <- rl$GENE1
  expect_identical(g1$SNP, paste0("rs", 50:1))          # reversed-POS order
  expect_true(is.matrix(g1$WEIGHT) && is.numeric(g1$WEIGHT))
  expect_identical(dim(g1$WEIGHT), c(50L, 1L))
  expect_identical(colnames(g1$WEIGHT), "WEIGHT")
  expect_identical(unname(g1$WEIGHT[, 1]), rev(w_gene1))
  # spot alignment: the value at fixture position k follows rs<k>
  # single-element matrix extraction inherits the colname as a name in base R
  expect_identical(unname(g1$WEIGHT[match("rs1", g1$SNP), 1]), 7.8)
  expect_identical(unname(g1$WEIGHT[match("rs20", g1$SNP), 1]), 9.9)
  expect_identical(unname(g1$WEIGHT[match("rs45", g1$SNP), 1]), 10)
  expect_identical(g1$genoIndex[match("rs1", g1$SNP)], "1001")
  expect_identical(unname(colSums(g1$annoIndicatorMat)), c(30, 50))

  g2 <- rl$GENE2
  expect_identical(g2$SNP, paste0("rs", 100:51))
  expect_identical(dim(g2$WEIGHT), c(50L, 1L))
  expect_identical(colnames(g2$WEIGHT), "WEIGHT")
  expect_identical(unname(g2$WEIGHT[, 1]), rev(w_gene2))
  expect_identical(unname(g2$WEIGHT[match("rs61", g2$SNP), 1]), 10)
  expect_identical(unname(g2$WEIGHT[match("rs100", g2$SNP), 1]), 1.7)
  expect_identical(unname(colSums(g2$annoIndicatorMat)), c(20, 50))

  # explicit new-style call is identical to the inferred legacy call
  explicit <- SAIGE.getRegionList_new(lines, 3, c("lof", "lof;missense"),
                                      mi_scrambled(), "",
                                      nWeightSets = 1L, weightSetNames = "WEIGHT")
  expect_identical(explicit, rl)
})

test_that("getRegionList_new on multiWeights: 2-col CADD/REVEL matrix, all markers kept", {
  lines <- readLines(gf_multiweights())
  rl <- SAIGE.getRegionList_new(lines, 4, c("lof", "lof;missense"),
                                mi_scrambled(), "",
                                nWeightSets = 2L,
                                weightSetNames = c("CADD", "REVEL"))
  expect_identical(names(rl), c("GENE1", "GENE2"))
  for (g in c("GENE1", "GENE2")) {
    w <- rl[[g]]$WEIGHT
    expect_true(is.matrix(w) && is.numeric(w))
    expect_identical(dim(w), c(50L, 2L))
    expect_identical(colnames(w), c("CADD", "REVEL"))
  }
  expect_identical(unname(rl$GENE1$WEIGHT[, "CADD"]), rev(w_gene1))
  expect_identical(unname(rl$GENE1$WEIGHT[, "REVEL"]), rev(revel_vals))
  expect_identical(unname(rl$GENE2$WEIGHT[, "CADD"]), rev(w_gene2))
  expect_identical(unname(rl$GENE2$WEIGHT[, "REVEL"]), rev(revel_vals))
  # per-marker rows survive as (CADD, REVEL) pairs from the fixture position
  expect_identical(unname(rl$GENE1$WEIGHT[match("rs14", rl$GENE1$SNP), ]),
                   c(7.4, 0.94))
  expect_identical(unname(rl$GENE2$WEIGHT[match("rs74", rl$GENE2$SNP), ]),
                   c(5.7, 0.76))
})

test_that("anno filtering on multiWeights drops missense rows from both weight columns", {
  lines <- readLines(gf_multiweights())
  rl <- SAIGE.getRegionList_new(lines, 4, c("lof"), mi_scrambled(), "",
                                nWeightSets = 2L,
                                weightSetNames = c("CADD", "REVEL"))
  g1 <- rl$GENE1
  # lof markers in GENE1 are rs1-rs10 and rs31-rs50, POS-reversed
  expect_identical(g1$SNP, paste0("rs", c(50:31, 10:1)))
  expect_identical(dim(g1$WEIGHT), c(30L, 2L))
  expect_identical(unname(g1$WEIGHT[, "CADD"]), w_gene1[c(50:31, 10:1)])
  expect_identical(unname(g1$WEIGHT[, "REVEL"]), revel_vals[c(50:31, 10:1)])
  expect_identical(unname(g1$WEIGHT[match("rs31", g1$SNP), ]), c(8.3, 0.44))

  g2 <- rl$GENE2
  # lof markers in GENE2 are rs81-rs100 (fixture positions 31-50)
  expect_identical(g2$SNP, paste0("rs", 100:81))
  expect_identical(dim(g2$WEIGHT), c(20L, 2L))
  expect_identical(unname(g2$WEIGHT[, "CADD"]), w_gene2[50:31])
  expect_identical(unname(g2$WEIGHT[, "REVEL"]), revel_vals[50:31])
  expect_identical(unname(g2$WEIGHT[match("rs81", g2$SNP), ]), c(7, 0.44))
})

test_that("unweighted control keeps the c(0) WEIGHT sentinel under the same markerInfo", {
  rl <- SAIGE.getRegionList_new(readLines(gf_unweighted()), 2,
                                c("lof", "lof;missense"), mi_scrambled(), "")
  expect_identical(rl$GENE1$WEIGHT, c(0))
  expect_identical(rl$GENE2$WEIGHT, c(0))
  expect_identical(rl$GENE1$SNP, paste0("rs", 50:1))
})

test_that("multiWeights round-trips checkGroupFile(detail=TRUE) into getRegionList_new", {
  # the same handoff SPAGMMATtest performs: layout comes from the detail return
  d <- checkGroupFile(gf_multiweights(), detail = TRUE)
  expect_identical(d$nRegions, 2)
  lines <- readLines(gf_multiweights())
  rl <- SAIGE.getRegionList_new(lines, 2 + d$nWeightSets,
                                c("lof", "lof;missense"), mi_scrambled(), "",
                                nWeightSets = d$nWeightSets,
                                weightSetNames = d$weightSetNames)
  expect_identical(length(rl), 2L)
  expect_identical(colnames(rl$GENE1$WEIGHT), d$weightSetNames)
  expect_identical(colnames(rl$GENE2$WEIGHT), d$weightSetNames)
  # identical to the same call with hand-written layout literals
  literal <- SAIGE.getRegionList_new(lines, 4, c("lof", "lof;missense"),
                                     mi_scrambled(), "",
                                     nWeightSets = 2L,
                                     weightSetNames = c("CADD", "REVEL"))
  expect_identical(rl, literal)
})
