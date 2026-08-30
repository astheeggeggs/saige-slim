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
