# LEDGER #76 FIXED: a --chrom that matches no marker.
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
