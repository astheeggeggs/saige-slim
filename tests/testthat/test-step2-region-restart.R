# LEDGER #25/#81 FIXED, the restart itself: a region run interrupted after its
# first region resumes at the second and reproduces the complete run.
#
# Reconstructed rather than killed (same reasoning as the marker restart test):
# the complete run's results file is cut back to region 1's rows and the index
# to "chunk 1", then the same command is re-run with --is_overwrite_output=
# FALSE. Before the fix no region restart could get past index validation.
# Tier: standard (two group runs).

test_that("a region run restarted after region 1 reproduces the complete run", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  run <- function(out, extra = character()) {
    st <- run_saige_cli("step2_SPAtests.R", c(
      "--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
      "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
      paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      paste0("--GMMATmodelFile=", m, ".rda"), paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
      paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=FALSE", extra))
    expect_cli_ok(st)
    readLines(attr(st, "log"), warn = FALSE)
  }
  complete <- file.path(tempfile("step2_rr_full_"), "region.txt"); dir.create(dirname(complete), recursive = TRUE)
  run(complete)
  full <- readLines(complete)
  idx <- readLines(paste0(complete, ".index"), warn = FALSE)
  expect_identical(tail(idx, 1), "Have completed the analyses of all chunks.")

  restarted <- file.path(tempfile("step2_rr_part_"), "region.txt"); dir.create(dirname(restarted), recursive = TRUE)
  # keep the header and every GENE1 row; the index says chunk 1 is done
  writeLines(full[c(1, grep("^GENE1\t", full))], restarted)
  writeLines(c(idx[1:3], "Have completed the analysis of chunk 1"), paste0(restarted, ".index"))
  for (sfx in c(".singleAssoc.txt", ".markerList.txt")) {
    if (file.exists(paste0(complete, sfx))) {
      s <- readLines(paste0(complete, sfx)); writeLines(s[c(1, grep("^GENE1\t", s))], paste0(restarted, sfx))
    }
  }
  lg <- run(restarted, "--is_overwrite_output=FALSE")
  expect_true(any(grepl("we restart the analysis from the 2 chunk", lg, fixed = TRUE)))
  expect_identical(readLines(restarted), full)
})
