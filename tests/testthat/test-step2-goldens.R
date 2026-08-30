# End-to-end step-2 regression tests: run the four PLINK step-2 golden commands
# through the real CLI and compare outputs against harness/goldens/step2/.
# Byte-exact in default mode (step 2 is deterministic; see plan/LEDGER), or
# numeric-tolerance in SAIGE_COMPARE_MODE=numeric for cross-platform runs.
#
# Frozen model inputs: our committed step-1 goldens, plus
# extdata/output/example_binary_sparseGRM.{rda,varianceRatio.txt} for group
# tests (the only legacy committed artifact still compatible; step-1 sparse-fit
# models have no `sparse`-type variance ratio rows, so group tests need it).

PLINK_ARGS <- function() c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first"
)
SG_ARGS <- c(
  "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
  "--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
)
GRP_MODEL <- c(
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
)
bin_model <- function() c(
  paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt"))
)
quant_model <- function() c(
  paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "quant_sparseGRM.rda")),
  paste0("--varianceRatioFile=", file.path(saige_goldens, "step1", "quant_sparseGRM.varianceRatio.txt"))
)

run_step2_case <- function(golden_base, extra_args, extra_outputs = character()) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step2_"), golden_base)
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(extra_args, paste0("--SAIGEOutputFile=", out)))
  expect_cli_ok(status)
  for (suffix in c("", extra_outputs)) {
    expect_output_matches_golden(
      file.path(saige_goldens, "step2", paste0(golden_base, suffix)),
      paste0(out, suffix)
    )
  }
}

test_that("step2 single-variant binary (ER path, minMAC=1) matches golden", {
  run_step2_case("plink_marker_binary.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
      bin_model(), "--is_output_moreDetails=TRUE"))
})

test_that("step2 group fastTest binary matches golden", {
  run_step2_case("plink_group_fastTest.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--is_fastTest=TRUE"),
    extra_outputs = c(".singleAssoc.txt", ".singleAssoc.txt_temp"))
})

test_that("step2 group conditional binary (rs30,rs79) matches golden", {
  run_step2_case("plink_group_cond_binary.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--condition=rs30,rs79", "--is_fastTest=TRUE"),
    extra_outputs = ".singleAssoc.txt")
})

test_that("step2 group BURDEN (r.corr=1) matches golden", {
  run_step2_case("plink_group_burden.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--r.corr=1", "--is_fastTest=FALSE"))
  # BURDEN-only runs emit no .singleAssoc.txt (only SKAT-O forces singles on)
})

test_that("step2 group non-fastTest (sparse-GRM org path) matches golden", {
  run_step2_case("plink_group_nofast.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--is_fastTest=FALSE"),
    extra_outputs = ".singleAssoc.txt")
})

test_that("step2 chunked single-variant (multi-chunk append) matches golden", {
  run_step2_case("plink_marker_chunked.txt",
    c("--bedFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bed",
      "--bimFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim",
      "--famFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.fam",
      "--AlleleOrder=alt-first",
      "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
      bin_model(), "--markers_per_chunk=1000"))
})

test_that("step2 single-variant conditional quantitative (rs13,rs79) matches golden", {
  run_step2_case("plink_marker_cond_quantitative.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=20", "--LOCO=FALSE",
      quant_model(), "--is_output_moreDetails=TRUE", "--condition=rs13,rs79"))
})
