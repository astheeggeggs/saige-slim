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

test_that("step2 group BURDEN conditional (r.corr=1, rs30,rs79) matches golden", {
  # P25, added with the LEDGER #41 fix. The combination r.corr=1 + --condition
  # is the one through which a halved conditional p-value reached a GENE-level
  # statistic: the C++ BURDEN writer stores the printed string and
  # src/Main.cpp:2019 parses it back with std::stod into CCT_cpp, so the region
  # Cauchy row inherited the halving. Nothing frozen covered it before -- P3/P10
  # are conditional but SKAT-O, P5/P9 are BURDEN but unconditional -- which is
  # why #41 lived in every conditional golden without any of them failing.
  run_step2_case("plink_group_cond_burden.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--r.corr=1", "--is_fastTest=FALSE", "--condition=rs30,rs79"))
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

test_that("step2 single-variant conditional BINARY (rs13,rs79) matches golden", {
  # The regression pin for LEDGER #40: the conditional SPA block is gated on
  # traitType != "quantitative", so the quantitative sibling above cannot
  # reach it. minMAC=1 retains the ultra-rare ER-path markers (rs44, rs56)
  # whose conditional statistics exceed the cutoff while their marginal block
  # is skipped -- the exact combination that read uninitialised state.
  run_step2_case("plink_marker_cond_binary.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
      bin_model(), "--is_output_moreDetails=TRUE", "--condition=rs13,rs79"))
})

# ---- multi-weight goldens (Phase 5 C10) ------------------------------------
# P2's configuration throughout, so a diff against plink_group_fastTest.txt
# isolates the weight dimension.

test_that("step2 group multi-weight SKAT-O (K=2) matches golden", {
  run_step2_case("plink_group_multiweight.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid_multiWeights.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--is_fastTest=TRUE"),
    extra_outputs = c(".singleAssoc.txt", ".singleAssoc.txt_temp"))
})

test_that("step2 group multi-weight BURDEN (K=2, r.corr=1) matches golden", {
  # different emitter from the SKAT-O case: the Weight column here is written
  # by the C++ BURDEN writer (C7), not the R post-processing.
  # is_fastTest=FALSE mirrors the single-weight BURDEN golden; the TRUE
  # combination loses regions (LEDGER #27) and is characterized separately.
  run_step2_case("plink_group_multiweight_burden.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid_multiWeights.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--r.corr=1", "--is_fastTest=FALSE"))
})

test_that("step2 group conditional multi-weight (K=2, rs30,rs79) matches golden", {
  # P10 (C11, design D7): P3's conditional configuration on the K=2 custom
  # file. The conditioning-weight vector is shared across sets; the per-set
  # G1-side adjustment is formed in R from the raw exported matrix.
  run_step2_case("plink_group_cond_multiweight.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid_multiWeights.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      "--condition=rs30,rs79", "--is_fastTest=TRUE"),
    extra_outputs = ".singleAssoc.txt")
})

test_that("step2 group composed sets (2 beta pairs + No_Weights) matches golden", {
  run_step2_case("plink_group_composed.txt",
    c(PLINK_ARGS(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      GRP_MODEL, SG_ARGS,
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
      # system2 goes through a shell, so the semicolon needs quoting
      paste0("--weights.beta=", shQuote("1,25;1,1")),
      "--is_no_weight_in_groupTest=TRUE",
      "--is_fastTest=TRUE"),
    extra_outputs = c(".singleAssoc.txt", ".singleAssoc.txt_temp"))
})
