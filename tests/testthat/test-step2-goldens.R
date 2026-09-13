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

# ---- FlexRV goldens (F4, 2026-09-09; P58-P61 in harness/make_goldens.sh) ----
# One annotation group, one max MAF (0.01: at the model's 1,000 samples that
# is MAC <= 20, the only cell of this fixture with real markers at the default
# collapse cutoff -- helper-flexrv.R has the arithmetic).
FLEXRV_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                   "--groupFile=./input/group_new_snpid_flexrv.txt",
                   "--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
                   "--annotation_in_groupTest=lof:missense")

test_that("step2 FlexRV 192-set BURDEN (binary, collapse on) matches golden", {
  run_step2_case("plink_group_flexrv_burden.txt",
    c(PLINK_ARGS(), GRP_MODEL, SG_ARGS, FLEXRV_COMMON, "--r.corr=1", "--is_fastTest=FALSE"))
})

test_that("step2 FlexRV 192-set BURDEN (quantitative, no collapse) matches golden", {
  run_step2_case("plink_group_flexrv_burden_quant.txt",
    c(PLINK_ARGS(), quant_model(), FLEXRV_COMMON, "--MACCutoff_to_CollapseUltraRare=0",
      "--r.corr=1", "--is_fastTest=FALSE"))
})

test_that("step2 FlexRV reduced-spec SKAT-O (2 x 2 sets) matches golden", {
  # max MAF 0.02 and mafbeta rather than a MAC threshold, so no in-cell weight
  # is exactly 0. That was chosen when LEDGER #117 made such a cell NA out;
  # #117 is FIXED (2026-09-11, a zero weight now drops the marker), so the
  # avoidance is no longer load-bearing -- the golden is kept because it gates
  # the reduced spec, the transformWeightMat export and the per-set ACAT
  # columns, and test-step2-flexrv.R covers the zero-weight case directly.
  run_step2_case("plink_group_flexrv_skato.txt",
    c(PLINK_ARGS(), GRP_MODEL, SG_ARGS,
      "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      "--groupFile=./input/group_new_snpid_flexrv.txt",
      "--flexRV_score=REVEL", "--flexRV_maxMAF=0.02",
      "--annotation_in_groupTest=lof:missense",
      # system2 goes through a shell, so the semicolons need quoting
      paste0("--flexRV_scoreTransforms=", shQuote("unw;sbeta:2")),
      paste0("--flexRV_MAFTransforms=", shQuote("unw;mafbeta:2")),
      "--is_fastTest=FALSE"),
    extra_outputs = ".singleAssoc.txt")
})

test_that("step2 FlexRV 192-set BURDEN through the pgen reader matches golden", {
  run_step2_case("pgen_group_flexrv_burden.txt",
    c("--pgenFile=./input/genotype_100markers.pgen",
      "--pvarFile=./input/genotype_100markers.pvar",
      "--psamFile=./input/genotype_100markers.psam",
      GRP_MODEL, SG_ARGS, FLEXRV_COMMON, "--r.corr=1", "--is_fastTest=FALSE"))
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

# ---------------------------------------------------------------------------
# P42-P48: missing genotype calls. The step-2 tier re-runs goldens from THIS
# EXPLICIT LIST, not by globbing harness/goldens/step2, so a golden that is not
# named here gates nothing however carefully it was generated.
MISSING_ARGS <- function() c(
  "--bedFile=./input/genotype_100markers_missing.bed",
  "--bimFile=./input/genotype_100markers_missing.bim",
  "--famFile=./input/genotype_100markers_missing.fam",
  "--AlleleOrder=alt-first"
)
MISSING_PGEN_ARGS <- function() c(
  "--pgenFile=./input/genotype_100markers_missing.pgen",
  "--pvarFile=./input/genotype_100markers_missing.pvar",
  "--psamFile=./input/genotype_100markers_missing.psam",
  "--AlleleOrder=alt-first"
)
MISSING_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
                    "--is_output_moreDetails=TRUE")

test_that("step2 marker path with missing calls, CLI defaults, matches golden", {
  run_step2_case("plink_marker_missing.txt",
    c(MISSING_ARGS(), MISSING_COMMON, bin_model()))
})

test_that("step2 --impute_method=mean matches golden", {
  run_step2_case("plink_marker_missing_mean.txt",
    c(MISSING_ARGS(), MISSING_COMMON, bin_model(), "--impute_method=mean"))
})

test_that("step2 --impute_method=minor matches golden", {
  run_step2_case("plink_marker_missing_minor.txt",
    c(MISSING_ARGS(), MISSING_COMMON, bin_model(), "--impute_method=minor"))
})

test_that("step2 --maxMissing=1 (gate off) matches golden", {
  run_step2_case("plink_marker_missing_maxmiss1.txt",
    c(MISSING_ARGS(), MISSING_COMMON, bin_model(), "--maxMissing=1"))
})

test_that("step2 mean imputation with --dosage_zerod_cutoff=0 matches golden", {
  run_step2_case("plink_marker_missing_meannozerod.txt",
    c(MISSING_ARGS(), MISSING_COMMON, bin_model(),
      "--impute_method=mean", "--dosage_zerod_cutoff=0"))
})

test_that("step2 PGEN mirror under missing calls matches golden", {
  run_step2_case("pgen_marker_missing.txt",
    c(MISSING_PGEN_ARGS(), MISSING_COMMON, bin_model()))
})

# P48 retired 2026-09-09 (LEDGER #108): --is_imputed_data=TRUE is now refused,
# so there is no output to compare. The refusal itself is pinned by
# test-step2-missing.R; nothing here can run that flag any more.

# ---------------------------------------------------------------------------
# P49-P51: dosages. Only the pgen reader can carry them, and before these no
# golden ran the dosage decode path at all.
DOSAGE_ARGS <- function() {
  c("--pgenFile=./input/genotype_100markers_dosage.pgen",
    "--pvarFile=./input/genotype_100markers_dosage.pvar",
    "--psamFile=./input/genotype_100markers_dosage.psam",
    "--AlleleOrder=alt-first")
}
DOSAGE_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
                   "--is_output_moreDetails=TRUE")

test_that("step2 dosage pgen, CLI defaults, matches golden", {
  run_step2_case("pgen_marker_dosage.txt",
    c(DOSAGE_ARGS(), DOSAGE_COMMON, bin_model()))
})

test_that("step2 dosage pgen with --dosage_zerod_cutoff=0 matches golden", {
  run_step2_case("pgen_marker_dosage_nozerod.txt",
    c(DOSAGE_ARGS(), DOSAGE_COMMON, bin_model(), "--dosage_zerod_cutoff=0"))
})

test_that("step2 dosage pgen with --dosage_zerod_MAC_cutoff=100 matches golden", {
  run_step2_case("pgen_marker_dosage_mac100.txt",
    c(DOSAGE_ARGS(), DOSAGE_COMMON, bin_model(), "--dosage_zerod_MAC_cutoff=100"))
})

# ---------------------------------------------------------------------------
# P52-P53: LOCO. The same model and command with only --chrom changed, so the
# pair isolates ReadModel's LOCOResult[[chrom]] swap and nothing else.
LOCO_GOLD_G <- function() {
  c("--bedFile=./input/genotype_100markers_2chr.bed",
    "--bimFile=./input/genotype_100markers_2chr.bim",
    "--famFile=./input/genotype_100markers_2chr.fam",
    "--AlleleOrder=alt-first")
}
loco_gold_model <- function() {
  c(paste0("--GMMATmodelFile=",
           file.path(saige_goldens, "step1", "binary_loco_2chr.rda")),
    paste0("--varianceRatioFile=",
           file.path(saige_goldens, "step1", "binary_loco_2chr.varianceRatio.txt")))
}
LOCO_GOLD_COMMON <- c("--minMAF=0", "--minMAC=1", "--LOCO=TRUE",
                      "--is_output_moreDetails=TRUE")

test_that("step2 LOCO on chromosome 1 matches golden", {
  run_step2_case("plink_marker_loco_chr1.txt",
    c(LOCO_GOLD_G(), "--chrom=1", LOCO_GOLD_COMMON, loco_gold_model()))
})

test_that("step2 LOCO on chromosome 2 matches golden", {
  run_step2_case("plink_marker_loco_chr2.txt",
    c(LOCO_GOLD_G(), "--chrom=2", LOCO_GOLD_COMMON, loco_gold_model()))
})

# ---------------------------------------------------------------------------
# P54-P55: the REGION path under LOCO. Every other group golden runs
# --LOCO=FALSE, so ReadModel's per-chromosome swap was gated only on the
# single-variant path. GENE2 spans the chromosome boundary in this fixture, so
# it is tested against a different refit in each of the two runs.
LOCO_R_MODEL <- function() {
  c(paste0("--GMMATmodelFile=",
           file.path(saige_goldens, "step1", "binary_loco_vr_2chr.rda")),
    paste0("--varianceRatioFile=",
           file.path(saige_goldens, "step1", "binary_loco_vr_2chr.varianceRatio.txt")))
}
LOCO_R_COMMON <- c("--minMAF=0", "--minMAC=0.5", "--LOCO=TRUE", "--is_fastTest=FALSE",
                   "--groupFile=./input/group_new_snpid.txt",
                   "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                   "--maxMAF_in_groupTest=0.0001,0.001,0.01")

test_that("step2 region test under LOCO, chromosome 1, matches golden", {
  args <- c(LOCO_GOLD_G(), "--chrom=1", LOCO_R_COMMON, LOCO_R_MODEL(), SG_ARGS)
  run_step2_case("plink_group_loco_chr1.txt", args,
                 extra_outputs = ".singleAssoc.txt")
})

test_that("step2 region test under LOCO, chromosome 2, matches golden", {
  args <- c(LOCO_GOLD_G(), "--chrom=2", LOCO_R_COMMON, LOCO_R_MODEL(), SG_ARGS)
  run_step2_case("plink_group_loco_chr2.txt", args,
                 extra_outputs = ".singleAssoc.txt")
})

# ---------------------------------------------------------------------------
# P56-P57: the quantitative region path's multi-weight and conditional shapes.
# Binary had 14 fastTest, 9 multi-weight and 9 conditional region goldens and
# quantitative had none of the three; these close two of them. Both use the only
# quantitative model with `sparse` variance-ratio rows, so they run the same
# sparse-GRM region path LEDGER #106 was first-order on.
QUANT_VR_MODEL <- function() {
  c(paste0("--GMMATmodelFile=",
           file.path(saige_goldens, "step1", "quant_sparseGRM_vr.rda")),
    paste0("--varianceRatioFile=",
           file.path(saige_goldens, "step1", "quant_sparseGRM_vr.varianceRatio.txt")))
}
QUANT_GRP_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                      "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                      "--is_fastTest=TRUE")

test_that("step2 quantitative multi-weight region test (K=2) matches golden", {
  args <- c(PLINK_ARGS(), QUANT_VR_MODEL(), SG_ARGS, QUANT_GRP_COMMON,
            "--groupFile=./input/group_new_snpid_multiWeights.txt")
  run_step2_case("plink_group_multiweight_quant.txt", args,
                 extra_outputs = c(".singleAssoc.txt", ".singleAssoc.txt_temp"))
})

test_that("step2 quantitative conditional region test matches golden", {
  args <- c(PLINK_ARGS(), QUANT_VR_MODEL(), SG_ARGS, QUANT_GRP_COMMON,
            "--groupFile=./input/group_new_snpid.txt", "--condition=rs30,rs79")
  run_step2_case("plink_group_cond_quant.txt", args,
                 extra_outputs = ".singleAssoc.txt")
})

# ---- P62-P64: the ultra-rare pseudo-marker's two defects, both FIXED -------
# These three goldens were added at unpatched HEAD so the fixes for LEDGER #129
# and #130 would MOVE something. Two claims this block first made were wrong
# and are corrected in full at harness/make_goldens.sh (read it there): the
# pre-existing UR-row count is 359, not 269 -- none of them above AF 0.5, so
# #130 was indeed unreached -- and "no group golden carries a custom weight
# small enough" is FALSE, falsified by the #129 fix itself, which moved 17
# pre-existing goldens through group_new_snpid_multiWeights.txt's REVEL weights
# of 0.01-0.07. P64 is the minimal reproduction, not the only golden.
# The commands and the fixture arithmetic are recorded beside the runs in
# harness/make_goldens.sh; test-step2-region-flip.R and
# test-step2-region-ur-zerod.R are the characterization tests that say what the
# numbers in them MEAN.
URFLIP_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                   "--groupFile=./input/group_new_snpid.txt",
                   "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                   "--maxMAF_in_groupTest=0.5", "--MACCutoff_to_CollapseUltraRare=2000",
                   "--is_single_in_groupTest=TRUE", "--is_fastTest=FALSE")

test_that("step2 high-AF ultra-rare pseudo-marker, BURDEN path, matches golden", {
  run_step2_case("plink_group_urflip_burden.txt",
    c(PLINK_ARGS(), GRP_MODEL, SG_ARGS, URFLIP_COMMON, "--r.corr=1"),
    extra_outputs = ".singleAssoc.txt")
})

test_that("step2 high-AF ultra-rare pseudo-marker, SKAT-O path, matches golden", {
  run_step2_case("plink_group_urflip_skato.txt",
    c(PLINK_ARGS(), GRP_MODEL, SG_ARGS, URFLIP_COMMON),
    extra_outputs = ".singleAssoc.txt")
})

test_that("step2 MIXED cell (rare + high-AF pseudo-marker) matches golden", {
  # P65: the only bundled configuration where a cell holds BOTH surviving rare
  # markers and a collapsed pseudo-marker above AF 0.5 -- the regime in which
  # the pseudo-marker's coding changes the statistic rather than only its sign.
  # The cutoff of 700 is the fixture arithmetic (harness/make_goldens.sh); the
  # mixed-cell properties are asserted in test-step2-region-flip.R.
  run_step2_case("plink_group_urmixed_burden.txt",
    c(PLINK_ARGS(), GRP_MODEL, SG_ARGS,
      "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.5", "--MACCutoff_to_CollapseUltraRare=700",
      "--is_single_in_groupTest=TRUE", "--r.corr=1", "--is_fastTest=FALSE"),
    extra_outputs = ".singleAssoc.txt")
})

test_that("step2 MIXED cell at --minGroupMAC_in_BurdenTest=1500 matches golden", {
  # P66 (LEDGER #135 FIXED): P65 with the gate raised past N. The emitter used
  # to fold a cell's AGGREGATE count as a frequency, so cells of 3064, 1915 and
  # 2997 alleles gated as 1000, 85 and 1000 and this file was two `Cauchy NA`
  # rows. It now refuses only GENE1 lof (748); the mechanism is characterized in
  # test-step2-burden-aggregate-mac.R.
  run_step2_case("plink_group_urmixed_burden_gate1500.txt",
    c(PLINK_ARGS(), GRP_MODEL, SG_ARGS,
      "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.5", "--MACCutoff_to_CollapseUltraRare=700",
      "--is_single_in_groupTest=TRUE", "--r.corr=1", "--is_fastTest=FALSE",
      "--minGroupMAC_in_BurdenTest=1500"),
    extra_outputs = c(".singleAssoc.txt", ".pooledTests.txt"))
})

test_that("step2 all-ultra-rare cell with small custom weights matches golden", {
  # LEDGER #129 FIXED: the cell used to vanish to the dosage filter, and this
  # golden was a bare Cauchy NA row plus the #127 sidecar disclosing the zero.
  # It is now one tested cell at MAC 2 with Pvalue_Burden 2.522646E-01 -- the
  # value the row had already recorded for the SAME data at weights 0.5 and at
  # collapse 0 -- and no sidecar at all, since the region has nothing left to
  # disclose. That absence is asserted in test-step2-region-ur-zerod.R; here it
  # is simply not in extra_outputs, and compare_goldens.sh would flag an EXTRA
  # file if one came back.
  run_step2_case("plink_group_ur_smallweights.txt",
    c(PLINK_ARGS(), GRP_MODEL, SG_ARGS,
      "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      "--groupFile=./input/group_ur_smallWeights.txt",
      "--annotation_in_groupTest=lof", "--maxMAF_in_groupTest=0.01",
      "--MACCutoff_to_CollapseUltraRare=10", "--r.corr=1", "--is_fastTest=FALSE"))
})
