# Shared scaffolding for the multi-weight (K > 1) step-2 tests, split across
# test-step2-multiweight{,-cond,-behaviour}.R. It lives in a helper so all
# three files (which testthat may run in different parallel workers) share one
# definition; mw_cache is per-process, so a run needed by two of the files is
# executed once per worker rather than once overall -- a deliberate trade of a
# little duplicated work for a much shorter critical path.
#
# Nothing here runs at source time: mw_run() calls skip_unless_level("standard")
# before it spawns anything, so the fast tier pays only the definitions.

# Multi-weight (K > 1) region tests, end-to-end through the step-2 CLI.
#
# Two runs of the P2 group-fastTest configuration (same model, sparse GRM,
# annotations and max-MAF cutoffs as test-step2-goldens.R's "group fastTest"
# case), differing only in the group file:
#   - group_new_snpid_multiWeights.txt  (K = 2: weight:CADD, weight,REVEL)
#   - group_new_snpid_withWeights.txt   (K = 1: one bare weight line)
#
# What is pinned here:
#   K = 2 -> rows are (anno x max_MAF x weight set) with the weight innermost,
#            a Weight column right after max_MAF, additive Pvalue_ACATV /
#            Pvalue_ACATO columns (Pvalue stays SKAT-O), and exactly ONE
#            pooled Cauchy row per region.
#   K = 1 -> the legacy schema, byte-for-byte the same columns as the
#            committed unweighted golden: no Weight/ACATV/ACATO.
#   Either way a real marker is score-tested ONCE, no matter how many weight
#   sets are in play; only the collapsed ultra-rare pseudo-markers are per set.
#
# Tier: standard (runs the CLI), like the golden tests.

MW_ARGS <- function(groupFile, out) c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first",
  "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
  "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
  "--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
  paste0("--groupFile=./input/", groupFile),
  "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
  "--maxMAF_in_groupTest=0.0001,0.001,0.01",
  "--is_fastTest=TRUE",
  paste0("--SAIGEOutputFile=", out)
)

# One CLI run per (group file, extra args) pair, cached across the test_that
# blocks below.
mw_cache <- new.env(parent = emptyenv())
mw_run <- function(groupFile, extra = character(0)) {
  skip_unless_level("standard")
  key <- paste(c(groupFile, extra), collapse = " ")
  if (!is.null(mw_cache[[key]])) return(mw_cache[[key]])
  out <- file.path(tempfile("step2_mw_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R", c(MW_ARGS(groupFile, out), extra))
  expect_cli_ok(status)
  # BURDEN-only runs do not emit the single-variant file
  single_path <- paste0(out, ".singleAssoc.txt")
  res <- list(region = read_saige_table(out),
              single = if (file.exists(single_path)) read_saige_table(single_path) else NULL,
              # the file itself, for the runs that also have a committed golden
              path = out)
  mw_cache[[key]] <- res
  res
}

multi_run <- function() mw_run("group_new_snpid_multiWeights.txt")
single_run <- function() mw_run("group_new_snpid_withWeights.txt")
# BURDEN-only (--r.corr=1): the region file is written from C++
# (openOutfile/writeOutfile_BURDEN), not the R SKAT-O post-processing
burden_multi_run <- function() mw_run("group_new_snpid_multiWeights.txt", "--r.corr=1")
burden_single_run <- function() mw_run("group_new_snpid_withWeights.txt", "--r.corr=1")
# C8 set composition (D10): explicit beta sets from a weights.beta pair list
# on an UNWEIGHTED group file (the semicolon needs shell quoting), the
# appended No_Weights set on a weighted file, and the forced Weight column at
# a single set (SKAT-O and BURDEN writers)
beta_run <- function() mw_run("group_new_snpid.txt",
                              paste0("--weights.beta=", shQuote("1,25;1,1")))
noweights_run <- function() mw_run("group_new_snpid_multiWeights.txt",
                                   "--is_no_weight_in_groupTest=TRUE")
forced_run <- function() mw_run("group_new_snpid_withWeights.txt",
                                "--is_output_weight_column=TRUE")
burden_forced_run <- function() mw_run("group_new_snpid_withWeights.txt",
                                       c("--r.corr=1", "--is_output_weight_column=TRUE"))

WSETS <- c("CADD", "REVEL")
pnum <- function(x) suppressWarnings(as.numeric(x))
acat_p <- function(p) 0.5 - atan(mean(tan((0.5 - p) * pi))) / pi

# Legacy BURDEN header written by openOutfile for a binary trait, no
# conditioning (Region.txt is C++-written when --r.corr=1)
BURDEN_HDR_K1 <- c("Region", "Group", "max_MAF", "Pvalue_Burden", "BETA_Burden",
                   "SE_Burden", "MAC", "MAC_case", "MAC_control",
                   "Number_rare", "Number_ultra_rare")

# ---- C8: weights.beta pair list -> explicit beta sets ----------------------

BETA_SETS <- c("Beta_1_25", "Beta_1_1")

# ---- conditional multi-weight (C11, design D7) ------------------------------
# The C9 guard is lifted: --condition now composes with K > 1 and with the
# forced Weight column. The conditioning-marker weight vector stays a SINGLE
# vector across sets (weights_for_condition, else pdf(beta(1,25))); what
# varies per set is the G1-side weight, applied in region_postprocess from the
# raw exported G1tilde_P_G2tilde matrix (C++ bakes it only at q_weight == 1).
#
# Reconciliation strategy (anything computable two ways must reconcile): the
# per-set R adjustment at K = 2 against the K = 1 C++-baked path. There are now
# TWO routes for that comparison, and the history matters.
#
# The BETA route (Beta_1_25 of a composed run vs the plain implicit-beta run)
# was the only one available until the LEDGER #17 fix: a custom-weighted group
# file with --condition read an UNINITIALIZED w0 on every collapsed ultra-rare
# row, and on this platform the garbage zeroed the adjustment outright, so the
# custom fixture could not be a baseline for anything.
#
# **#17 is fixed** (w0 = 1.0 for a custom set -- its UR collapse is weighted, so
# the weight is already in the pseudo-marker), and the CUSTOM route is now the
# sharper of the two: group_new_snpid_withWeights.txt's weight line is
# byte-identical to the CADD line of group_new_snpid_multiWeights.txt, so the
# K = 1 run and the K = 2 run's CADD rows are the same statistic computed by two
# different paths, with no beta-distribution algebra in between.
#
# Both routes are kept. Statistics are compared at 1e-12 relative (algebraically
# identical, different float grouping -- never byte, per the cross-path rule in
# CLAUDE.md).

cond_multi_run <- function() mw_run("group_new_snpid_multiWeights.txt",
                                    "--condition=rs30,rs79")
cond_plain_run <- function() mw_run("group_new_snpid.txt",
                                    "--condition=rs30,rs79")
cond_beta2_run <- function() mw_run("group_new_snpid.txt",
                                    c(paste0("--weights.beta=", shQuote("1,25;1,1")),
                                      "--condition=rs30,rs79"))
# K = 1, CUSTOM weights (CADD only) + --condition: the path LEDGER #17 was in.
cond_single_run <- function() mw_run("group_new_snpid_withWeights.txt",
                                     "--condition=rs30,rs79")
cond_forced_run <- function() mw_run("group_new_snpid.txt",
                                     c("--is_output_weight_column=TRUE",
                                       "--condition=rs30,rs79"))

COND_STATS <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden",
                "SE_Burden", "Pvalue_cond", "Pvalue_Burden_cond",
                "Pvalue_SKAT_cond", "BETA_Burden_cond", "SE_Burden_cond",
                "MAC", "MAC_case", "MAC_control",
                "Number_rare", "Number_ultra_rare")
expect_stats_close <- function(a, b, cols = COND_STATS, tol = 1e-12) {
  for (cc in cols) {
    x <- pnum(a[[cc]])
    y <- pnum(b[[cc]])
    expect_equal(x, y, tolerance = tol, label = cc)
  }
}

# ---- the implicit_beta set (LEDGER #19) -------------------------------------
# An UNWEIGHTED group file with --is_output_weight_column=TRUE composes no set
# at all, so region_postprocess falls to kind "implicit_beta": the literal
# legacy dbeta(MAFVec, 1, 25). This is the only way to observe that arm, and
# it is where the R/C++ weight asymmetry becomes visible in the output.
implicit_beta_run <- function(beta) {
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_implbeta_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- c(MW_ARGS("group_new_snpid.txt", out),
            "--is_output_weight_column=TRUE",
            paste0("--weights.beta=", beta))
  # single anno/max_MAF keeps the run small; the asymmetry is per row
  args[grep("^--annotation_in_groupTest=", args)] <- "--annotation_in_groupTest=lof"
  args[grep("^--maxMAF_in_groupTest=", args)] <- "--maxMAF_in_groupTest=0.01"
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
  read_saige_table(out)
}
