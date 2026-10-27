# Group-file parsing: checkGroupFile (R/SAIGE_SPATest_Region_Func.R) and
# SAIGE.getRegionList_new (R/SAIGE_SPATest_Region.R) — the region-test input
# contract. PLINK group files key on bim column 2 (rs IDs).

groupfile <- function() file.path(saige_extdata, "input", "group_new_snpid.txt")
mi100 <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      fam <- utils::read.table(file.path(saige_extdata, "input", "genotype_100markers.fam"),
                               colClasses = "character")
      cache <<- setGenoInput(bedFile = file.path(saige_extdata, "input", "genotype_100markers.bed"),
                             sampleInModel = fam$V2)$markerInfo
    }
    cache
  }
})

test_that("checkGroupFile on the bundled rsID group file", {
  expect_identical(checkGroupFile(groupfile()),
                   list(nRegions = 2, is_weight_included = FALSE))
})

test_that("checkGroupFile detects the 3-line weighted format", {
  lines <- readLines(groupfile())
  n1 <- length(strsplit(lines[1], "[ \t]+")[[1]]) - 2
  n2 <- length(strsplit(lines[3], "[ \t]+")[[1]]) - 2
  weighted <- c(lines[1:2], paste(c("GENE1", "weight", rep("1", n1)), collapse = " "),
                lines[3:4], paste(c("GENE2", "weight", rep("1", n2)), collapse = " "))
  f <- tempfile(); writeLines(weighted, f)
  expect_identical(checkGroupFile(f), list(nRegions = 2, is_weight_included = TRUE))
})

test_that("checkGroupFile error paths", {
  lines <- readLines(groupfile())
  f1 <- tempfile(); writeLines(lines[1:3], f1)   # GENE2 lost its anno line
  expect_error(checkGroupFile(f1), "Group file is incomplete for GENE2")

  bad <- lines; bad[2] <- sub(" anno ", " var ", bad[2])
  f2 <- tempfile(); writeLines(bad, f2)
  expect_error(checkGroupFile(f2), "No anno is specified")

  short <- lines
  short[2] <- paste(strsplit(short[2], "[ \t]+")[[1]][1:10], collapse = " ")
  f3 <- tempfile(); writeLines(short, f3)
  expect_error(checkGroupFile(f3),
               "not equal to the length of marker IDs")
})

test_that("SAIGE.getRegionList_new builds regions from the bundled fixtures", {
  rl <- SAIGE.getRegionList_new(readLines(groupfile()), 2,
                                c("lof", "lof;missense"), mi100(), "")
  expect_identical(names(rl), c("GENE1", "GENE2"))
  g1 <- rl$GENE1
  expect_identical(g1$SNP, paste0("rs", 1:50))         # restored to POS order
  expect_identical(g1$WEIGHT, c(0))                    # scalar 0 when unweighted
  expect_identical(dim(g1$annoIndicatorMat), c(50L, 2L))
  expect_identical(unname(colSums(g1$annoIndicatorMat)), c(30, 50))
  # genoIndex is character, format()-padded to common width (leading spaces!)
  expect_identical(g1$genoIndex[1], " 0")
  expect_identical(g1$genoIndex[50], "49")

  g2 <- rl$GENE2
  expect_identical(unname(colSums(g2$annoIndicatorMat)), c(20, 50))
  expect_identical(g2$genoIndex[1], "50")
})

test_that("annotation filtering drops unmatched markers and regions", {
  rl <- SAIGE.getRegionList_new(readLines(groupfile()), 2, c("lof"), mi100(), "")
  expect_identical(nrow(rl$GENE2$annoIndicatorMat), 20L)
  expect_identical(rl$GENE2$SNP[1], "rs81")
})

test_that("missing markers are dropped with a message; duplicates are errors", {
  mi_minus <- mi100()[!mi100()$ID == "rs1", ]
  rl <- SAIGE.getRegionList_new(readLines(groupfile()), 2,
                                c("lof", "lof;missense"), mi_minus, "")
  expect_identical(length(rl$GENE1$SNP), 49L)

  dup <- readLines(groupfile())[c(1, 2, 1, 2)]
  expect_error(SAIGE.getRegionList_new(dup, 2, c("lof"), mi100(), ""),
               "GENE1 is duplicated in the group File")
})

# --- multi-weight-set group files: checkGroupFile(detail = TRUE) ---
# Region 1 defines K weight lines per region; the 2nd field of a weight line
# is "weight" (set name WEIGHT), "weight:NAME" or "weight,NAME".

multiweight_groupfile <- function(tokens1, tokens2 = tokens1, value = "0.5") {
  lines <- readLines(groupfile())
  n1 <- length(strsplit(lines[1], "[ \t]+")[[1]]) - 2
  n2 <- length(strsplit(lines[3], "[ \t]+")[[1]]) - 2
  w1 <- vapply(tokens1, function(tk)
    paste(c("GENE1", tk, rep(value, n1)), collapse = " "), character(1))
  w2 <- vapply(tokens2, function(tk)
    paste(c("GENE2", tk, rep(value, n2)), collapse = " "), character(1))
  f <- tempfile()
  writeLines(c(lines[1:2], w1, lines[3:4], w2), f)
  f
}

test_that("detail=TRUE on the unweighted bundled file appends K=0 fields", {
  d <- checkGroupFile(groupfile(), detail = TRUE)
  expect_identical(d, list(nRegions = 2, is_weight_included = FALSE,
                           nWeightSets = 0L, weightSetNames = character(0)))
  expect_identical(d[1:2], checkGroupFile(groupfile()))
})

test_that("detail=TRUE on the legacy bare 3-line format names the set WEIGHT", {
  f <- multiweight_groupfile("weight", value = "1")
  expect_identical(checkGroupFile(f, detail = TRUE),
                   list(nRegions = 2, is_weight_included = TRUE,
                        nWeightSets = 1L, weightSetNames = "WEIGHT"))
})

test_that("two weight sets are parsed with both ':' and ',' delimiters", {
  f <- multiweight_groupfile(c("weight:CADD", "weight,REVEL"))
  expect_identical(checkGroupFile(f, detail = TRUE),
                   list(nRegions = 2, is_weight_included = TRUE,
                        nWeightSets = 2L, weightSetNames = c("CADD", "REVEL")))
})

test_that("multi-weight-set validation errors", {
  # weight lines permuted in region 2 (the fork silently mislabels; we error)
  f <- multiweight_groupfile(c("weight:CADD", "weight:REVEL"),
                             c("weight:REVEL", "weight:CADD"))
  expect_error(checkGroupFile(f),
               "must appear in the same order in every region")

  # duplicate set name within region 1 (via both delimiters)
  f <- multiweight_groupfile(c("weight:CADD", "weight,CADD"))
  expect_error(checkGroupFile(f), "Duplicate weight set name 'CADD'")

  # "Cauchy" is reserved
  f <- multiweight_groupfile("weight:Cauchy")
  expect_error(checkGroupFile(f), "reserved weight set name")

  # weight-count != marker-count on one line
  lines <- readLines(groupfile())
  n1 <- length(strsplit(lines[1], "[ \t]+")[[1]]) - 2
  n2 <- length(strsplit(lines[3], "[ \t]+")[[1]]) - 2
  bad <- c(lines[1:2],
           paste(c("GENE1", "weight:CADD", rep("1", n1 - 1)), collapse = " "),
           lines[3:4],
           paste(c("GENE2", "weight:CADD", rep("1", n2)), collapse = " "))
  f <- tempfile(); writeLines(bad, f)
  expect_error(checkGroupFile(f),
               "The length of weights for markers in region GENE1")

  # non-numeric weight values
  f <- multiweight_groupfile("weight:CADD", value = "abc")
  expect_error(checkGroupFile(f),
               "Non-numeric weight values in weight set 'CADD' for region GENE1")
})

test_that("detail=FALSE keeps the legacy two-element return for every format", {
  expect_identical(checkGroupFile(groupfile()),
                   list(nRegions = 2, is_weight_included = FALSE))
  f <- multiweight_groupfile("weight", value = "1")
  expect_identical(checkGroupFile(f),
                   list(nRegions = 2, is_weight_included = TRUE))
  f2 <- multiweight_groupfile(c("weight:CADD", "weight,REVEL"))
  expect_identical(checkGroupFile(f2),
                   list(nRegions = 2, is_weight_included = TRUE))
})

test_that("negative weights are rejected only in multi-weight mode", {
  # K = 1: historical behavior preserved -- weights were never sign-checked
  # upstream (they enter only multiplicatively; no ACAT-V at K = 1)
  f1 <- multiweight_groupfile("weight", value = "-1")
  expect_identical(checkGroupFile(f1),
                   list(nRegions = 2, is_weight_included = TRUE))

  # K > 1: ACAT-V hands these to CCT as Cauchy combination weights, which
  # stop()s mid-run on a negative -- so the parser fails fast with location
  f2 <- multiweight_groupfile(c("weight:CADD", "weight:REVEL"), value = "-1")
  expect_error(checkGroupFile(f2),
               "Negative weight values in weight set 'CADD' for region GENE1")
})

# ---------------------------------------------------------------------------
# FlexRV F1: the optional `score:NAME` line. One per region, values in [0, 1],
# fixed position (immediately after anno, before any weight line). Parsed and
# carried into RegionData as SAIGESC_1; NOTHING consumes it yet, so a file
# without one must behave exactly as before.

# build a group file from the bundled one with a score line spliced in, and
# optionally K weight lines after it
gf_with_score <- function(scoreType = "score:PAI3D", vals = NULL,
                          weights = 0L, vals2 = NULL) {
  lines <- readLines(groupfile())
  n1 <- length(strsplit(lines[1], "[ \t]+")[[1]]) - 2
  n2 <- length(strsplit(lines[3], "[ \t]+")[[1]]) - 2
  s1 <- if (is.null(vals)) rep("0.5", n1) else vals
  s2 <- if (is.null(vals2)) rep("0.25", n2) else vals2
  wl <- function(gene, n) {
    if (weights == 0L) {
      return(character(0))
    }
    mk <- function(k) {
      paste(c(gene, paste0("weight:W", k), rep("1", n)), collapse = " ")
    }
    vapply(seq_len(weights), mk, character(1))
  }
  out <- c(lines[1:2], paste(c("GENE1", scoreType, s1), collapse = " "), wl("GENE1", n1),
           lines[3:4], paste(c("GENE2", scoreType, s2), collapse = " "), wl("GENE2", n2))
  f <- tempfile()
  writeLines(out, f)
  f
}

test_that("a score line is parsed and named, and the detail=FALSE contract is unchanged", {
  f <- gf_with_score()
  # the historical return must stay EXACTLY list(nRegions, is_weight_included):
  # a score line is not a weight line, so is_weight_included stays FALSE
  expect_identical(checkGroupFile(f), list(nRegions = 2, is_weight_included = FALSE))
  d <- checkGroupFile(f, detail = TRUE)
  expect_identical(d$scoreName, "PAI3D")
  expect_identical(d$nWeightSets, 0L)
  # a bare `score` type is named SCORE, like a bare `weight` is named WEIGHT
  expect_identical(checkGroupFile(gf_with_score("score"), detail = TRUE)$scoreName, "SCORE")
  # comma form, as for weight lines
  expect_identical(checkGroupFile(gf_with_score("score,PAI3D"), detail = TRUE)$scoreName,
                   "PAI3D")
})

test_that("a file with no score line keeps the historical detail=TRUE return", {
  d <- checkGroupFile(groupfile(), detail = TRUE)
  # NOT character(0): the field is absent entirely, so the four-element
  # contract that three older tests pin is preserved. length(NULL) == 0 is
  # what the callers need.
  expect_null(d$scoreName)
  expect_identical(names(d), c("nRegions", "is_weight_included",
                               "nWeightSets", "weightSetNames"))
  expect_identical(d$nRegions, 2)
})

test_that("score lines coexist with weight lines, which keep their own names", {
  f <- gf_with_score(weights = 2L)
  d <- checkGroupFile(f, detail = TRUE)
  expect_identical(d$scoreName, "PAI3D")
  expect_identical(d$nWeightSets, 2L)
  expect_identical(d$weightSetNames, c("W1", "W2"))
  expect_true(d$is_weight_included)
})

test_that("score values are validated: numeric, non-NA, and inside [0, 1]", {
  n1 <- length(strsplit(readLines(groupfile())[1], "[ \t]+")[[1]]) - 2
  bad <- function(v) {
    vals <- rep("0.5", n1)
    vals[1] <- v
    gf_with_score(vals = vals)
  }
  expect_error(checkGroupFile(bad("1.5")), "outside \\[0, 1\\]")
  expect_error(checkGroupFile(bad("-0.1")), "outside \\[0, 1\\]")
  expect_error(checkGroupFile(bad("NA")), "Non-numeric or NA score")
  expect_error(checkGroupFile(bad("abc")), "Non-numeric or NA score")
  # the boundaries themselves are allowed
  vals <- rep("0.5", n1)
  vals[1] <- "0"
  vals[2] <- "1"
  expect_identical(checkGroupFile(gf_with_score(vals = vals), detail = TRUE)$scoreName,
                   "PAI3D")
})

test_that("every region must carry the same score, of the same length", {
  lines <- readLines(groupfile())
  n1 <- length(strsplit(lines[1], "[ \t]+")[[1]]) - 2
  n2 <- length(strsplit(lines[3], "[ \t]+")[[1]]) - 2
  # region 2 names a different score
  mixed <- c(lines[1:2], paste(c("GENE1", "score:A", rep("0.5", n1)), collapse = " "),
             lines[3:4], paste(c("GENE2", "score:B", rep("0.5", n2)), collapse = " "))
  f <- tempfile()
  writeLines(mixed, f)
  expect_error(checkGroupFile(f), "does not match score name")
  # region 2 has no score line at all
  missing <- c(lines[1:2], paste(c("GENE1", "score:A", rep("0.5", n1)), collapse = " "),
               lines[3:4])
  f2 <- tempfile()
  writeLines(missing, f2)
  expect_error(checkGroupFile(f2))
  # wrong number of scores
  short <- c(lines[1:2], paste(c("GENE1", "score:A", rep("0.5", n1 - 1)), collapse = " "),
             lines[3:4], paste(c("GENE2", "score:A", rep("0.5", n2)), collapse = " "))
  f3 <- tempfile()
  writeLines(short, f3)
  expect_error(checkGroupFile(f3), "length of scores")
})

test_that("several score lines are allowed, together after anno, with distinct names", {
  lines <- readLines(groupfile())
  n1 <- length(strsplit(lines[1], "[ \t]+")[[1]]) - 2
  n2 <- length(strsplit(lines[3], "[ \t]+")[[1]]) - 2
  two <- c(lines[1:2],
           paste(c("GENE1", "score:A", rep("0.5", n1)), collapse = " "),
           paste(c("GENE1", "score:B", rep("0.5", n1)), collapse = " "),
           lines[3:4],
           paste(c("GENE2", "score:A", rep("0.5", n2)), collapse = " "),
           paste(c("GENE2", "score:B", rep("0.5", n2)), collapse = " "))
  f <- tempfile()
  writeLines(two, f)
  d <- checkGroupFile(f, detail = TRUE)
  expect_identical(d$scoreName, c("A", "B"))
  expect_identical(d$nWeightSets, 0L)
  # the same name twice is refused
  dup <- sub("score:B", "score:A", two, fixed = TRUE)
  f1 <- tempfile()
  writeLines(dup, f1)
  expect_error(checkGroupFile(f1), "Duplicate score name 'A'")
  # a score line AFTER a weight line is refused
  after <- c(lines[1:2],
             paste(c("GENE1", "weight", rep("1", n1)), collapse = " "),
             paste(c("GENE1", "score:A", rep("0.5", n1)), collapse = " "),
             lines[3:4],
             paste(c("GENE2", "weight", rep("1", n2)), collapse = " "),
             paste(c("GENE2", "score:A", rep("0.5", n2)), collapse = " "))
  f2 <- tempfile()
  writeLines(after, f2)
  expect_error(checkGroupFile(f2), "score line after a weight line")
})

test_that("getRegionList_new carries the score into RegionData as SAIGESC_1", {
  f <- gf_with_score()
  d <- checkGroupFile(f, detail = TRUE)
  lines <- readLines(f)
  rl <- SAIGE.getRegionList_new(lines, 2 + length(d$scoreName) + d$nWeightSets,
                                c("lof", "missense", "synonymous"), mi100(),
                                nWeightSets = d$nWeightSets,
                                weightSetNames = d$weightSetNames,
                                scoreName = d$scoreName)
  expect_identical(names(rl), c("GENE1", "GENE2"))
  # the score reached the per-region structure and kept its values
  g1 <- rl[["GENE1"]]
  expect_identical(g1$SCORE, rep(0.5, length(g1$SNP)))
  # and a scoreless file must not gain the field at all
  d0 <- checkGroupFile(groupfile(), detail = TRUE)
  rl0 <- SAIGE.getRegionList_new(readLines(groupfile()), 2,
                                 c("lof", "missense", "synonymous"), mi100(),
                                 nWeightSets = d0$nWeightSets,
                                 weightSetNames = d0$weightSetNames,
                                 scoreName = d0$scoreName)
  expect_null(rl0[["GENE1"]]$SCORE)
})
