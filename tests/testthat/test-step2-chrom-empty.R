# LEDGER #76 and #113(a) FIXED: two routes to an empty region run.
#
# #76 is the --chrom filter emptying the marker set; #113(a) is the ID JOIN
# between the group file and the variant file coming up empty while the
# filter is fine. Both used to reach the writer with nothing and produce a
# ZERO-BYTE results file at exit 0 with "Analysis done!".
#
# The original #76 note follows.
#
# Before the fix the REGION path wrote a ZERO-BYTE results file at exit 0 and
# printed "Analysis done!" -- markerInfo was filtered to nothing and nothing
# checked. The MARKER path happened to catch it downstream. Both readers now
# stop in R/Geno.R, where the filter is, naming the chromosome and the file.
#
# Every bundled fixture is chromosome 1 only, so --chrom=2 is the empty
# selection on all of them. Tier: standard (five step-2 runs, three of which
# fail fast).

chrom_run <- function(chrom, geno = c("plink", "pgen"), region = TRUE) {
  geno <- match.arg(geno)
  out <- tempfile(fileext = ".txt")
  # The region path needs a model whose variance-ratio file has `sparse` rows
  # (the frozen legacy model); the marker path needs `null_noXadj` rows (the
  # committed step-1 golden). Neither model serves both.
  m <- if (region) file.path(saige_extdata, "output", "example_binary_sparseGRM")
       else file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  gargs <- if (geno == "plink") {
    c("--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
      "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first")
  } else {
    c("--pgenFile=./input/genotype_100markers.pgen", "--pvarFile=./input/genotype_100markers.pvar",
      "--psamFile=./input/genotype_100markers.psam")
  }
  args <- c(gargs, paste0("--SAIGEOutputFile=", out), paste0("--chrom=", chrom),
            "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"))
  if (region) {
    args <- c(args, paste0("--sparseGRMFile=", sg),
              paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
              "--groupFile=./input/group_new_snpid.txt",
              "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
              "--maxMAF_in_groupTest=0.0001,0.001,0.01")
  }
  st <- run_saige_cli("step2_SPAtests.R", args)
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}

expect_empty_chrom_refused <- function(r, file_re) {
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl(paste0("No markers on chrom 2 are found in .*", file_re), r$log)))
  expect_false(any(grepl("Analysis done", r$log, fixed = TRUE)))
  # and no zero-byte results file is left looking like an answer
  expect_true(!file.exists(r$out) || file.size(r$out) > 0)
}

test_that("a region run at a --chrom with no markers stops instead of writing an empty file (PLINK)", {
  skip_unless_level("standard")
  expect_empty_chrom_refused(chrom_run("2", "plink"), "genotype_100markers\\.bim")
})

test_that("... and through the PGEN reader", {
  skip_unless_level("standard")
  expect_empty_chrom_refused(chrom_run("2", "pgen"), "genotype_100markers\\.pvar")
})

test_that("the marker path refuses the same way (it already did; the message is now the reader's)", {
  skip_unless_level("standard")
  expect_empty_chrom_refused(chrom_run("2", "plink", region = FALSE), "genotype_100markers\\.bim")
})

test_that("POSITIVE CONTROL: --chrom=1 still runs both paths to a non-empty file", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  for (region in c(TRUE, FALSE)) {
    r <- chrom_run("1", "plink", region = region)
    expect_identical(r$status, 0L)
    expect_true(file.exists(r$out) && file.size(r$out) > 0)
    expect_true(any(grepl("Analysis done", r$log, fixed = TRUE)))
  }
})


# ---------------------------------------------------------------------------
# LEDGER #113(a): the group file's marker IDs match NOTHING in the variant file.
# The --chrom filter is fine here (90 markers on chromosome 1); it is the join
# that empties. PLINK keys on bim column 2 and PGEN on the pvar ID column, so a
# group file written for another format matches nothing. The six dead
# chr:pos:ref:alt group files that made this easy to hit were retired by
# #113(b), so the fixture is written here instead -- a hand-written group file
# with a typo reaches the same path, which is why (a) stayed open.

group_file_with_ids <- function(ids) {
  f <- tempfile(fileext = ".txt")
  writeLines(c(paste("GENE1 var", paste(ids, collapse = " ")),
               paste("GENE1 anno", paste(rep("lof", length(ids)), collapse = " "))), f)
  f
}

nomatch_run <- function(gf, geno = c("plink", "pgen")) {
  geno <- match.arg(geno)
  out <- tempfile(fileext = ".txt")
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  gargs <- if (geno == "plink") {
    c("--bedFile=./input/genotype_100markers.bed",
      "--bimFile=./input/genotype_100markers.bim",
      "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first")
  } else {
    c("--pgenFile=./input/genotype_100markers.pgen",
      "--pvarFile=./input/genotype_100markers.pvar",
      "--psamFile=./input/genotype_100markers.psam")
  }
  st <- run_saige_cli("step2_SPAtests.R",
    c(gargs, paste0("--SAIGEOutputFile=", out), "--chrom=1",
      "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      paste0("--GMMATmodelFile=", m, ".rda"),
      paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
      paste0("--sparseGRMFile=", sg),
      paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
      paste0("--groupFile=", gf),
      "--annotation_in_groupTest=lof",
      "--maxMAF_in_groupTest=0.01"))
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}

test_that("LEDGER #113(a) FIXED: a group file matching no marker stops (PLINK)", {
  skip_unless_level("standard")
  gf <- group_file_with_ids(c("1:100:A:C", "1:200:G:T", "1:300:C:A"))
  r <- nomatch_run(gf, "plink")
  expect_false(identical(r$status, 0L))
  txt <- paste(r$log, collapse = "\n")
  expect_match(txt, "None of the 3 markers in the group file were found", fixed = TRUE)
  expect_match(txt, "PLINK keys on column 2 of the .bim", fixed = TRUE)
  expect_false(any(grepl("Analysis done", r$log, fixed = TRUE)))
  # The results file may still EXIST and be empty: openOutfile() runs at
  # R/SAIGE_SPATest_Region.R:952, long before the group file is ingested at
  # :1364, so the guard cannot prevent its creation -- and deleting it would be
  # unsafe, since under --is_overwrite_output=FALSE that same file can hold
  # earlier chunks' real results. What the guard removes is the SUCCESS: exit is
  # non-zero and "Analysis done!" is absent, so nothing downstream can read the
  # empty file as an answer. That is the difference from #76, whose guard runs
  # in R/Geno.R before any output file is opened.
})

test_that("... and through the PGEN reader", {
  skip_unless_level("standard")
  gf <- group_file_with_ids(c("1:100:A:C", "1:200:G:T"))
  r <- nomatch_run(gf, "pgen")
  expect_false(identical(r$status, 0L))
  expect_match(paste(r$log, collapse = "\n"),
               "None of the 2 markers in the group file were found", fixed = TRUE)
  expect_false(any(grepl("Analysis done", r$log, fixed = TRUE)))
})

test_that("POSITIVE CONTROL: a PARTIAL match still runs, and empty regions stay legal", {
  skip_unless_level("standard")
  # The guard must fire only when NOTHING matched. Twenty real markers among two
  # fakes must still produce an answer, and the pre-existing "N markers in
  # 'RegionFile' are not in 'GenoFile'" notice must still be the only complaint.
  # Twenty rather than two: a gene small enough to yield no qualifying cell
  # writes an empty file at exit 0 by a DIFFERENT route (the cell filters, not
  # the join), which this control must not accidentally exercise.
  gf <- group_file_with_ids(c(paste0("rs", 1:20), "1:100:A:C", "1:200:G:T"))
  r <- nomatch_run(gf, "plink")
  expect_identical(r$status, 0L)
  txt <- paste(r$log, collapse = "\n")
  expect_match(txt, "markers in 'RegionFile' are not in 'GenoFile'", fixed = TRUE)
  expect_false(grepl("None of the", txt, fixed = TRUE))
  expect_true(any(grepl("Analysis done", r$log, fixed = TRUE)))
  expect_true(file.exists(r$out) && file.size(r$out) > 0)
})
