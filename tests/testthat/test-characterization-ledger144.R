# LEDGER #144 (FIXED 2026-09-14): `.markerList.txt` survives a region restart.
#
# LEDGER #59 seeded the results-file flag from `isappend` so a resumed run
# appends; the adjacent marker-list flag was seeded FALSE regardless, so on a
# restart the first resumed chunk with marker-list rows took the
# `append = FALSE, col.names = TRUE` branch and TRUNCATED the file to the
# resumed regions' rows. #140's restart truncation covers the results file
# and `.singleAssoc.txt` and says nothing about this one. The fix mirrors #59:
# `iswriteMarkerList <- rep(isappend, nModels)`.
#
# The test is test-step2-region-restart.R's reconstruction with the flag on:
# the complete run's three files are cut back to GENE1, the index says chunk
# 1 is done, the same command re-runs with --is_overwrite_output=FALSE, and
# the restarted `.markerList.txt` must equal the complete run's. Before the
# fix it held GENE2's rows only. Tier: standard (two group runs).

l144_cli <- function(...) get("run_saige_cli")(...)
l144_ok <- function(...) get("expect_cli_ok")(...)
l144_extdata <- function() get("saige_extdata")

l144_run <- function(out, extra = character()) {
  m <- file.path(l144_extdata(), "output", "example_binary_sparseGRM")
  sg <- file.path(l144_extdata(), "output",
                  "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  args <- c(
    "--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
    paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0", "--minMAC=0.5",
    "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", m, ".rda"), paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
    paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=FALSE",
    "--is_output_markerList_in_groupTest=TRUE", extra
  )
  st <- l144_cli("step2_SPAtests.R", args)
  l144_ok(st)
  readLines(attr(st, "log"), warn = FALSE)
}

test_that("a restarted region run keeps the pre-restart .markerList.txt rows (LEDGER #144)", {
  skip_unless_level("standard")
  complete <- file.path(tempfile("step2_l144_full_"), "region.txt")
  dir.create(dirname(complete), recursive = TRUE)
  l144_run(complete)
  full <- readLines(complete)
  full_ml <- readLines(paste0(complete, ".markerList.txt"))
  idx <- readLines(paste0(complete, ".index"), warn = FALSE)
  expect_identical(tail(idx, 1), "Have completed the analyses of all chunks.")
  # the fixture really has marker-list rows for BOTH regions
  expect_gt(length(grep("^GENE1\t", full_ml)), 0)
  expect_gt(length(grep("^GENE2\t", full_ml)), 0)

  restarted <- file.path(tempfile("step2_l144_part_"), "region.txt")
  dir.create(dirname(restarted), recursive = TRUE)
  writeLines(full[c(1, grep("^GENE1\t", full))], restarted)
  writeLines(c(idx[1:3], "Have completed the analysis of chunk 1"), paste0(restarted, ".index"))
  for (sfx in c(".singleAssoc.txt", ".markerList.txt")) {
    s <- readLines(paste0(complete, sfx))
    writeLines(s[c(1, grep("^GENE1\t", s))], paste0(restarted, sfx))
  }
  lg <- l144_run(restarted, "--is_overwrite_output=FALSE")
  expect_true(any(grepl("we restart the analysis from the 2 chunk", lg, fixed = TRUE)))
  expect_identical(readLines(restarted), full)
  # the inverted characterization: before the fix this held the header and
  # GENE2's rows only
  expect_identical(readLines(paste0(restarted, ".markerList.txt")), full_ml)
})
