# All-numeric marker IDs (LEDGER #74).
#
# R/Geno.R reads the bim with `fread(bimFile, header = FALSE,
# select = c(1, 2, 4, 5, 6))` and no `colClasses`, so data.table INFERS the type
# of the ID column. Every bundled fixture uses rs-prefixed IDs, which infer
# character; a bim whose IDs are all numeric infers INTEGER, and every join
# against it then fails, because group-file and --condition IDs are always
# character.
#
# The fixtures are derived at run time and reuse the committed .bed and .fam
# UNCHANGED -- a .bed is indexed by row, so rewriting bim column 2 changes no
# genotype. That is also what makes the negative control exact: the same
# genotypes, the same model, only the ID spelling differs.

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
  grp <- readLines(nid_group)
  writeLines(gsub("rs([0-9]+)", "\\1", grp), file.path(d, "group_numeric.txt"))
  d
}

nid_model <- local({
  m <- file.path(saige_goldens, "step1", "quant_sparseGRM")
  c(paste0("--GMMATmodelFile=", m, ".rda"),
    paste0("--varianceRatioFile=", m, ".varianceRatio.txt"))
})
nid_base <- c("--AlleleOrder=alt-first", "--chrom=1", "--minMAF=0",
              "--minMAC=0.5", "--LOCO=FALSE")
nid_cli <- run_saige_cli

# `numeric` picks the derived numeric-ID bim; otherwise the committed one.
nid_run <- function(dir, numeric_ids, extra) {
  pre <- if (numeric_ids) file.path(dir, "num") else nid_src
  out <- tempfile(fileext = ".txt")
  st <- nid_cli("step2_SPAtests.R",
                c(paste0("--bedFile=", pre, ".bed"),
                  paste0("--bimFile=", pre, ".bim"),
                  paste0("--famFile=", pre, ".fam"),
                  paste0("--SAIGEOutputFile=", out),
                  nid_base, nid_model, extra))
  list(status = as.integer(st), out = out,
       log = readLines(attr(st, "log"), warn = FALSE))
}

joined_ids <- function(r) any(grepl("Incompatible join types", r$log, fixed = TRUE))

test_that("NEGATIVE CONTROL: a plain marker scan is unaffected by numeric IDs", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- num_ids_dir()
  num <- nid_run(d, TRUE, character())
  orig <- nid_run(d, FALSE, character())
  expect_identical(num$status, 0L)
  expect_identical(orig$status, 0L)
  # the ID column is not joined on this path, so the two runs analyse the same
  # markers and differ only in how the IDs are spelled
  dn <- read_saige_table(num$out)
  do <- read_saige_table(orig$out)
  expect_identical(nrow(dn), nrow(do))
  expect_identical(dn$MarkerID, sub("^rs", "", do$MarkerID))
  expect_identical(dn$p.value, do$p.value)
})

test_that("LEDGER #74: --condition against an all-numeric bim fails the join", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- num_ids_dir()
  num <- nid_run(d, TRUE, "--condition=13")
  expect_false(identical(num$status, 0L))
  expect_true(joined_ids(num))
  # the message names both sides and their inferred types
  expect_true(any(grepl("x.ID (integer) and i.SNP (character)", num$log, fixed = TRUE)))
  # POSITIVE CONTROL: the same conditioning marker, same .bed, rs-spelled
  orig <- nid_run(d, FALSE, "--condition=rs13")
  expect_identical(orig$status, 0L)
  expect_false(joined_ids(orig))
})

test_that("LEDGER #74: a numeric-ID group file fails the same join", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- num_ids_dir()
  grp <- c("--annotation_in_groupTest=lof,missense:lof", "--maxMAF_in_groupTest=0.01")
  num <- nid_run(d, TRUE, c(paste0("--groupFile=", file.path(d, "group_numeric.txt")), grp))
  expect_false(identical(num$status, 0L))
  expect_true(joined_ids(num))
  orig <- nid_run(d, FALSE, c(paste0("--groupFile=", nid_group), grp))
  expect_identical(orig$status, 0L)
  expect_false(joined_ids(orig))
})
