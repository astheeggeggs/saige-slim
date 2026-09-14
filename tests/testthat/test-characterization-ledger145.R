# LEDGER #145 (FIXED 2026-09-14): --is_overwrite_output=TRUE removes a previous
# run's sidecars along with the results file.
#
# `.pooledTests.txt` (#127) and `.skatFailures.txt` (#57) are written only
# when a run has something to disclose, and `.markerList.txt` only under its
# flag, so an overwriting run that produces none of them used to leave the
# PREVIOUS run's sidecars on disk beside a log that says `SKAT failures: 0`
# and a pooled-row line saying every region was a plain Cauchy combination.
# Only the results file was removed (`Check_OutputFile_Create`). Now the
# region overwrite branch unlinks the three sidecars too. The inverted
# characterization: plant all three, run an overwriting configuration that
# writes none, assert they are gone. Tier: standard (one region run).

l145_cli <- function(...) get("run_saige_cli")(...)
l145_ok <- function(...) get("expect_cli_ok")(...)

L145_SIDECARS <- c(".pooledTests.txt", ".skatFailures.txt", ".markerList.txt")

test_that("an overwriting region run removes stale sidecars (LEDGER #145)", {
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_l145_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  sentinel <- "planted by test-characterization-ledger145.R"
  # a previous run's results file AND its three sidecars
  for (f in c("", L145_SIDECARS)) writeLines(sentinel, paste0(out, f))
  args <- c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--is_fastTest=TRUE", "--is_overwrite_output=TRUE",
    paste0("--SAIGEOutputFile=", out)
  )
  st <- l145_cli("step2_SPAtests.R", args)
  l145_ok(st)
  lg <- readLines(attr(st, "log"), warn = FALSE)
  # this configuration discloses nothing and lists no markers ...
  expect_true(any(grepl("SKAT failures: 0", lg, fixed = TRUE)))
  expect_true(any(grepl("every one a plain Cauchy combination", lg, fixed = TRUE)))
  # ... so the results file is fresh (the positive control the old code passed)
  expect_false(any(grepl(sentinel, readLines(out, warn = FALSE), fixed = TRUE)))
  # ... and the planted sidecars are gone (the inverted characterization: before
  # the fix all three survived with the sentinel in them)
  for (f in L145_SIDECARS) expect_false(file.exists(paste0(out, f)), label = f)
})
