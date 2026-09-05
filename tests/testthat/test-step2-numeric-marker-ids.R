# All-numeric marker IDs (LEDGER #74, FIXED).
#
# R/Geno.R used to read the bim with `fread(bimFile, header = FALSE,
# select = c(1, 2, 4, 5, 6))` and no `colClasses`, so data.table INFERRED the
# type of the ID column. Every bundled fixture uses rs-prefixed IDs, which infer
# character; a bim whose IDs are all numeric inferred INTEGER, and every join
# against it then failed, because group-file and --condition IDs come from
# strsplit and are always character. The fix pins bim column 2 (and pvar column
# 3) to character, so the ID's SPELLING no longer changes the analysis.
#
# The fixtures are derived at run time and reuse the committed .bed/.pgen and
# .fam/.psam UNCHANGED -- a genotype file is indexed by row, so rewriting the
# variant file's ID column changes no genotype. That is what makes the
# comparison exact: the same genotypes, the same model, only the ID spelling
# differs, so every non-ID column must agree to the byte.

# helper-saige.R supplies saige_extdata/saige_goldens; bound at top level so
# this file is self-contained for static analysis as well as at run time.
nid_src <- file.path(saige_extdata, "input", "genotype_100markers")
nid_group <- file.path(saige_extdata, "input", "group_new_snpid.txt")

num_ids_dir <- function() {
  d <- tempfile()
  dir.create(d)
  src <- nid_src
  bim <- utils::read.table(paste0(src, ".bim"), stringsAsFactors = FALSE)
  bim[[2]] <- sub("^rs", "", bim[[2]])
  stopifnot(!anyNA(suppressWarnings(as.numeric(bim[[2]]))))   # really all-numeric
  utils::write.table(bim, file.path(d, "num.bim"), sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  file.copy(paste0(src, ".bed"), file.path(d, "num.bed"))
  file.copy(paste0(src, ".fam"), file.path(d, "num.fam"))
  # the pgen mirror: same markers, same order, ID column rewritten the same way
  pvar <- readLines(paste0(src, ".pvar"))
  hdr <- startsWith(pvar, "#")
  body <- do.call(rbind, strsplit(pvar[!hdr], "\t", fixed = TRUE))
  body[, 3] <- sub("^rs", "", body[, 3])
  writeLines(c(pvar[hdr], apply(body, 1, paste, collapse = "\t")),
             file.path(d, "num.pvar"))
  file.copy(paste0(src, ".pgen"), file.path(d, "num.pgen"))
  file.copy(paste0(src, ".psam"), file.path(d, "num.psam"))
  grp <- readLines(nid_group)
  writeLines(gsub("rs([0-9]+)", "\\1", grp), file.path(d, "group_numeric.txt"))
  d
}

nid_model <- local({
  m <- file.path(saige_goldens, "step1", "quant_sparseGRM")
  c(paste0("--GMMATmodelFile=", m, ".rda"),
    paste0("--varianceRatioFile=", m, ".varianceRatio.txt"))
})
nid_base <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE")
nid_cli <- run_saige_cli
nid_read <- read_saige_table

# `numeric` picks the derived numeric-ID variant file; otherwise the committed
# one. `genoType` selects the reader, so the same assertions cover both halves
# of the fix (bim column 2 and pvar column 3).
nid_run <- function(dir, numeric_ids, extra, genoType = "plink") {
  pre <- if (numeric_ids) file.path(dir, "num") else nid_src
  geno <- if (genoType == "plink") {
    c(paste0("--bedFile=", pre, ".bed"), paste0("--bimFile=", pre, ".bim"),
      paste0("--famFile=", pre, ".fam"), "--AlleleOrder=alt-first")
  } else {
    c(paste0("--pgenFile=", pre, ".pgen"), paste0("--pvarFile=", pre, ".pvar"),
      paste0("--psamFile=", pre, ".psam"))
  }
  out <- tempfile(fileext = ".txt")
  st <- nid_cli("step2_SPAtests.R",
                c(geno, paste0("--SAIGEOutputFile=", out),
                  nid_base, nid_model, extra))
  list(status = as.integer(st), out = out,
       log = readLines(attr(st, "log"), warn = FALSE))
}

joined_ids <- function(r) any(grepl("Incompatible join types", r$log, fixed = TRUE))

# The two runs differ only in how the IDs are spelled, so every column except
# MarkerID must be identical -- and MarkerID itself must match after stripping
# the "rs" prefix, which proves the same markers were selected in the same order.
expect_same_but_ids <- function(num, orig) {
  dn <- nid_read(num$out)
  do <- nid_read(orig$out)
  testthat::expect_identical(nrow(dn), nrow(do))
  testthat::expect_identical(dn$MarkerID, sub("^rs", "", do$MarkerID))
  testthat::expect_identical(dn[setdiff(names(dn), "MarkerID")],
                             do[setdiff(names(do), "MarkerID")])
}

test_that("NEGATIVE CONTROL: a plain marker scan is unaffected by numeric IDs", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- num_ids_dir()
  num <- nid_run(d, TRUE, character())
  orig <- nid_run(d, FALSE, character())
  expect_identical(num$status, 0L)
  expect_identical(orig$status, 0L)
  # the ID column is not joined on this path, so it was never affected either
  expect_same_but_ids(num, orig)
})

test_that("LEDGER #74 FIXED: --condition against an all-numeric bim joins", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- num_ids_dir()
  num <- nid_run(d, TRUE, "--condition=13")
  expect_identical(num$status, 0L)
  expect_false(joined_ids(num))
  # and the conditional analysis is the SAME analysis: same conditioning
  # marker, same .bed, rs-spelled.
  orig <- nid_run(d, FALSE, "--condition=rs13")
  expect_identical(orig$status, 0L)
  expect_same_but_ids(num, orig)
  # the conditioning really happened -- otherwise the comparison above would
  # be satisfied by two unconditional runs
  expect_true("p.value_c" %in% names(nid_read(num$out)))
})

test_that("LEDGER #74 FIXED: a numeric-ID group file joins the same way", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- num_ids_dir()
  grp <- c("--annotation_in_groupTest=lof,missense:lof", "--maxMAF_in_groupTest=0.01")
  num <- nid_run(d, TRUE, c(paste0("--groupFile=", file.path(d, "group_numeric.txt")), grp))
  expect_identical(num$status, 0L)
  expect_false(joined_ids(num))
  orig <- nid_run(d, FALSE, c(paste0("--groupFile=", nid_group), grp))
  expect_identical(orig$status, 0L)
  # region results carry no marker IDs at all, so they must agree exactly
  expect_identical(readLines(num$out), readLines(orig$out))
  # the single-variant-in-group rows do, and differ only in the spelling
  expect_same_but_ids(list(out = paste0(num$out, ".singleAssoc.txt")),
                      list(out = paste0(orig$out, ".singleAssoc.txt")))
})

test_that("LEDGER #74 FIXED: the pgen half -- an all-numeric pvar ID column", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- num_ids_dir()
  num <- nid_run(d, TRUE, "--condition=13", genoType = "pgen")
  expect_identical(num$status, 0L)
  expect_false(joined_ids(num))
  orig <- nid_run(d, FALSE, "--condition=rs13", genoType = "pgen")
  expect_identical(orig$status, 0L)
  expect_same_but_ids(num, orig)
})
