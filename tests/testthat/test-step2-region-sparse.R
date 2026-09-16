# The carrier-space PGEN read on the REGION path (2026-09-16, second site of
# the marker path's prepareSparseMarker; test-step2-pgen-sparse.R is the first).
# A record the pgen writer stored as a difference list is prepared over its
# listed samples -- counts, impute/flip, QC, weights, pr.nz / pr.gnz, and the
# case/control tallies by rank -- and the dense vector is scattered only for a
# marker tested singly. The claim: sparse and dense are the SAME arithmetic,
# so the sparse run must reproduce the committed pgen region goldens byte for
# byte AND the dense run of the same command (SAIGE_PGEN_SPARSE=0), and the
# region profile line must REPORT how many markers took the sparse path.
#
# Cross-run comparisons: gated in BOTH compare modes. Tier: standard.

PGEN_G <- c("--pgenFile=./input/genotype_100markers.pgen",
            "--pvarFile=./input/genotype_100markers.pvar",
            "--psamFile=./input/genotype_100markers.psam")
MISS_PGEN_G <- c("--pgenFile=./input/genotype_100markers_missing.pgen",
                 "--pvarFile=./input/genotype_100markers_missing.pvar",
                 "--psamFile=./input/genotype_100markers_missing.psam")
GRP_MODEL <- c(
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
)
SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
QUANT_MODEL <- function() c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.varianceRatio.txt"))
)
GRP_FLAGS <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
               "--groupFile=./input/group_new_snpid.txt",
               "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
               "--maxMAF_in_groupTest=0.0001,0.001,0.01")

region_run <- function(args, tag, sparse = TRUE) {
  out <- file.path(tempfile(paste0("step2_rsparse_", tag, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  env <- c(SAIGE_PROFILE_REGION = "1")
  if (!sparse) env <- c(env, SAIGE_PGEN_SPARSE = "0")
  status <- run_saige_cli("step2_SPAtests.R",
                          c(args, paste0("--SAIGEOutputFile=", out)), env = env)
  list(status = status, path = out, log = attr(status, "log"))
}

# The cumulative profile line's sparse tally; the last line is the run's total.
region_sparse_reads <- function(log) {
  ln <- grep("^REGION PROFILE cum", readLines(log, warn = FALSE), value = TRUE)
  if (!length(ln)) return(NA_integer_)
  as.integer(sub("^.* sparse=([0-9]+) .*$", "\\1", ln[length(ln)]))
}

expect_region_sparse_equals_dense <- function(args, tag, golden = NULL) {
  s <- region_run(args, paste0(tag, "_sparse"), sparse = TRUE)
  d <- region_run(args, paste0(tag, "_dense"), sparse = FALSE)
  expect_cli_ok(s$status)
  expect_cli_ok(d$status)
  expect_gt(region_sparse_reads(s$log), 0L)
  expect_identical(region_sparse_reads(d$log), 0L)
  expect_files_byte_identical(d$path, s$path)
  if (!is.null(golden)) {
    expect_output_matches_golden(golden, s$path)
    expect_output_matches_golden(golden, d$path)
  }
  invisible(list(sparse = s, dense = d))
}

test_that("region sparse read reproduces golden P14 (pgen SKAT-O fastTest) and reports itself", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_region_sparse_equals_dense(
    c(PGEN_G, GRP_FLAGS, GRP_MODEL, SG, "--is_fastTest=TRUE"), "p14",
    golden = file.path(saige_goldens, "step2", "pgen_group_fastTest.txt"))
})

test_that("region sparse read reproduces golden P18 (pgen BURDEN, no single tests: never scatters)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_region_sparse_equals_dense(
    c(PGEN_G, GRP_FLAGS, GRP_MODEL, SG, "--r.corr=1", "--is_fastTest=FALSE"), "p18",
    golden = file.path(saige_goldens, "step2", "pgen_group_burden.txt"))
})

test_that("single-in-group with moreDetails: the case/control tallies by rank are the same bytes", {
  # AF_case / AF_ctrl (mean of the gathered vector), MAC_case / MAC_control
  # (its accu) and the hom/het counts (find over the gathered vector) are
  # computed over the carriers by rank on the sparse path.
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_region_sparse_equals_dense(
    c(PGEN_G, GRP_FLAGS, GRP_MODEL, SG, "--is_fastTest=TRUE",
      "--is_single_in_groupTest=TRUE", "--is_output_moreDetails=TRUE",
      "--is_output_markerList_in_groupTest=TRUE"), "single_more")
})

test_that("missing genotypes with mean imputation on the region path are the same bytes", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_region_sparse_equals_dense(
    c(MISS_PGEN_G, GRP_FLAGS, GRP_MODEL, SG, "--is_fastTest=TRUE",
      "--is_single_in_groupTest=TRUE", "--impute_method=mean"), "miss_mean")
})

test_that("quantitative model on the region path is the same bytes", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_region_sparse_equals_dense(
    c(PGEN_G, GRP_FLAGS, QUANT_MODEL(), "--is_fastTest=TRUE",
      "--is_single_in_groupTest=TRUE"), "quant")
})

test_that("four region threads over the sparse read are byte-identical to one and to the golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  one <- region_run(c(PGEN_G, GRP_FLAGS, GRP_MODEL, SG, "--is_fastTest=TRUE", "--nThreads=1"), "t1")
  four <- region_run(c(PGEN_G, GRP_FLAGS, GRP_MODEL, SG, "--is_fastTest=TRUE", "--nThreads=4"), "t4")
  expect_cli_ok(one$status)
  expect_cli_ok(four$status)
  # The profile line (and its sparse tally) is the SERIAL arm's instrument
  # and is not printed on a threaded run, so the count comes from the
  # one-thread run; the threaded run is gated by its bytes.
  expect_gt(region_sparse_reads(one$log), 0L)
  expect_files_byte_identical(one$path, four$path)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_fastTest.txt"), four$path)
})
