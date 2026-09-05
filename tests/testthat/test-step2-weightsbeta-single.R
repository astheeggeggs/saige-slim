# LEDGER #19 FIXED: a single non-default --weights.beta pair on an unweighted
# group file is honoured on BOTH sides of the region test. C++ always used the
# pair (genoSumMat accumulation, ultra-rare collapse scaling); the R-side
# SKAT/Burden weights in region_postprocess were the literal dbeta(MAF, 1, 25),
# so `--weights.beta=0.5,0.5` produced a hybrid test that was neither.
#
# The gate is one statistic computed two ways. The multi-pair path
# (`--weights.beta="0.5,0.5;1,25"`) composes EXPLICIT beta sets whose R-side
# weights were always dbeta(MAF, a_r, b_r) (design D2), so its Beta_0.5_0.5
# rows are the correct answer for that pair. The single-pair run must agree
# with them cell for cell. Two algebraically identical code paths differ by
# ULPs, so statistics are compared at 1e-10 relative, counts exactly (CLAUDE.md,
# "Never assert bit-identity between two different code paths").
#
# Configuration: golden P6 (plink_group_nofast.txt), is_fastTest = FALSE, so no
# region can escalate differently in the two runs -- escalation depends on the
# pooled Cauchy p, which the K = 2 run pools over both sets.
#
# Tier: standard (runs the CLI).

WB_ARGS <- function(out, weights_beta) {
  c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed",
           ".sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--is_fastTest=FALSE",
    paste0("--weights.beta=", shQuote(weights_beta)),
    paste0("--SAIGEOutputFile=", out)
  )
}

wb_cache <- new.env(parent = emptyenv())
wb_run <- function(weights_beta) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  if (!is.null(wb_cache[[weights_beta]])) return(wb_cache[[weights_beta]])
  out <- file.path(tempfile("step2_wb_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", WB_ARGS(out, weights_beta)))
  res <- list(region = read_saige_table(out), path = out)
  wb_cache[[weights_beta]] <- res
  res
}

wb_num <- function(x) suppressWarnings(as.numeric(x))
STAT_COLS <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")
COUNT_COLS <- c("MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")

# The K = 2 run's rows for ONE named weight set, in (Region, Group, max_MAF)
# order, Cauchy rows dropped (the K = 2 Cauchy row pools both sets).
set_rows <- function(region, set) {
  r <- region[region$Group != "Cauchy" & region$Weight == set, , drop = FALSE]
  r[order(r$Region, r$Group, r$max_MAF), , drop = FALSE]
}
cell_rows <- function(region) {
  r <- region[region$Group != "Cauchy", , drop = FALSE]
  r[order(r$Region, r$Group, r$max_MAF), , drop = FALSE]
}

expect_same_cells <- function(single, multi_set) {
  expect_identical(nrow(single), nrow(multi_set))
  expect_gt(nrow(single), 10)
  for (cn in c("Region", "Group", "max_MAF", COUNT_COLS)) {
    expect_identical(single[[cn]], multi_set[[cn]], label = paste("same", cn))
  }
  for (cn in STAT_COLS) {
    expect_lt(rel_err(wb_num(single[[cn]]), wb_num(multi_set[[cn]])), 1e-10,
              label = paste("rel_err", cn))
  }
}

test_that("LEDGER #19 FIXED: a single non-default pair matches its explicit beta set", {
  single <- cell_rows(wb_run("0.5,0.5")$region)
  multi <- wb_run("0.5,0.5;1,25")$region
  expect_true("Weight" %in% colnames(multi))
  expect_same_cells(single, set_rows(multi, "Beta_0.5_0.5"))
})

test_that("LEDGER #19: the default pair through the explicit path reproduces golden P6 (control)", {
  # the comparison method has power only if the explicit path is itself right:
  # its Beta_1_25 rows must be the committed default run
  multi <- wb_run("0.5,0.5;1,25")$region
  p6 <- cell_rows(read_saige_table(file.path(saige_goldens, "step2", "plink_group_nofast.txt")))
  expect_same_cells(p6, set_rows(multi, "Beta_1_25"))
  # and the two pairs are different tests, so the assertion above is not vacuous
  b1 <- set_rows(multi, "Beta_1_25")
  b5 <- set_rows(multi, "Beta_0.5_0.5")
  expect_gt(sum(b1$Pvalue_SKAT != b5$Pvalue_SKAT), 5)
})

test_that("LEDGER #19: the default single pair is byte-identical to golden P6 (no golden moves)", {
  res <- wb_run("1,25")
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_group_nofast.txt"),
                               res$path)
})
