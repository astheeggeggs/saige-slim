# LEDGER #183 (2026-10-28): the sparse-GRM variance guard.
#
# On the sparse-GRM band the variance of a test vector is g' Sigma^-1 g from a
# conjugate-gradient solve against tau1 I + tau2 K_thinned, and a thinned GRM
# is not positive definite in general: on amr Height a pseudo-marker's variance
# walked from 10.7 through 0.39 to below zero across the FlexRV weight ladder
# with the solve converged, printing p 1.4e-73 one step before the #44 collapse.
# The guard computes the diagonal (no-GRM) variance of the same projected
# vector beside the solve's and flags a test whose solve is not positive or is
# below the threshold times the diagonal value; by default it only observes.
#
# The bundled Sigma is positive definite (harness/tools/aou_sigma_pd_check.R
# measures it), so the defect is reached by the plant SAIGE_SIGMA_PLANT_RATIO,
# which scales the solve's value before the guard sees it, exactly as an
# indefinite Sigma would have:
#   1. guard on (the default): every golden byte stays, the run REPORTS the
#      count of evaluated tests and 0 flagged;
#   2. the plant at 0.1: every test on the band is flagged and listed, the
#      ratio column reads the plant, nothing negative, p-values moved;
#   3. the plant at -1: every test is negative and the results carry the #44
#      collapse (p 1, BETA NA) on exactly those cells;
#   4. the plant at 0.1 with SAIGE_SIGMA_GUARD_APPLY=1: the flagged tests take
#      the diagonal variance, the sidecar says so, and their reported p is the
#      sidecar's;
#   5. the kill switch: golden bytes, nothing evaluated, the line says off.
skip_if_no_goldens()

sg_common <- c("--bedFile=./input/genotype_100markers.bed",
               "--bimFile=./input/genotype_100markers.bim",
               "--famFile=./input/genotype_100markers.fam",
               "--AlleleOrder=alt-first",
               "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
               "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
               "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
               paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
                      "randomMarkersUsed.sparseGRM.mtx"),
               paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
                      "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
               "--is_fastTest=FALSE",
               "--groupFile=./input/group_new_snpid_withWeights_abovecap.txt",
               "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
               "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--r.corr=1")
sg_gold <- file.path(saige_goldens, "step2", "plink_group_customweight_abovecap_burden.txt")
sg_run <- function(env = character()) {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", c(sg_common, paste0("--SAIGEOutputFile=", out)),
                      env = env)
  expect_cli_ok(st)
  list(out = out, log = readLines(attr(st, "log"), warn = FALSE))
}
sg_line <- function(log) grep("^Sigma variance: ", log, value = TRUE)
sg_counts <- function(log) {
  # the four counts of the line: evaluated, flagged, not positive, took the diagonal
  l <- sg_line(log)
  pat <- "[0-9]+(?= (test\\(s\\)|flagged|not positive|took))"
  as.numeric(regmatches(l, gregexpr(pat, l, perl = TRUE))[[1]])
}
sg_side <- function(out) {
  utils::read.table(paste0(out, ".sigmaGuard.txt"), header = TRUE, sep = "\t",
                    colClasses = "character", stringsAsFactors = FALSE)
}
# the results row a sidecar's burden row names (K = 1 custom weight: no Weight column)
sg_row <- function(tab, side_row) {
  tab[tab$Region == side_row$Region & tab$Group == side_row$Group &
        as.numeric(tab$max_MAF) == as.numeric(side_row$max_MAF), ]
}
sg_results <- function(out) read_saige_table(out)

test_that("guard on by default: golden bytes, the run reports the band's tests and 0 flagged", {
  skip_unless_level("standard")
  r <- sg_run()
  expect_output_matches_golden(sg_gold, r$out)
  expect_length(sg_line(r$log), 1L)
  cnt <- sg_counts(r$log)
  expect_length(cnt, 4L)
  expect_gt(cnt[1], 0)                      # tests on the sparse-GRM band were evaluated
  expect_identical(cnt[2:4], c(0, 0, 0))
  expect_false(file.exists(paste0(r$out, ".sigmaGuard.txt")))
  expect_match(sg_line(r$log), "observation only")
})

test_that("plant 0.1: every test on the band is flagged, listed with the ratio, and moved", {
  skip_unless_level("standard")
  base <- sg_run()
  r <- sg_run(env = c(SAIGE_SIGMA_PLANT_RATIO = "0.1"))
  cnt <- sg_counts(r$log)
  expect_identical(cnt[1], sg_counts(base$log)[1])
  expect_identical(cnt[2], cnt[1])          # all flagged at 0.1 < 0.5
  expect_identical(cnt[3:4], c(0, 0))
  expect_match(sg_line(r$log), "PLANT: solve x 0.1", fixed = TRUE)
  side <- sg_side(r$out)
  expect_identical(nrow(side), as.integer(cnt[2]))
  expect_true(all(side$negative == "FALSE"))
  expect_true(all(side$reported == "grm"))
  # the ratio is the plant times the true ratio, which on the bundled PD Sigma
  # sits near 1: nowhere near zero
  ratio <- as.numeric(side$ratio)
  expect_true(all(ratio > 0.01 & ratio < 0.2))
  expect_true(all(as.numeric(side$var_grm) > 0 & as.numeric(side$var_diag) > 0))
  expect_true(all(side$Kind %in% c("ultra_rare", "burden")))
  # the p-values that moved are flagged cells, and some moved: a cell whose p
  # comes from the exact test (LEDGER #172 routing, MAC <= --max_MAC_for_ER)
  # never reads the variance, so not every flagged cell moves
  g <- sg_results(base$out)
  p <- sg_results(r$out)
  expect_identical(dim(g), dim(p))
  moved <- which(g$Pvalue_Burden != p$Pvalue_Burden & g$Group != "Cauchy")
  expect_gt(length(moved), 0)
  flagged <- paste(side$Region, side$Group, as.numeric(side$max_MAF))
  moved_cells <- paste(g$Region[moved], g$Group[moved], as.numeric(g$max_MAF[moved]))
  expect_true(all(moved_cells %in% flagged))
  # and the sidecar's reported p IS the results table's, cell for cell
  for (k in which(side$Kind == "burden")) {
    row <- sg_row(p, side[k, ])
    expect_identical(nrow(row), 1L)
    expect_identical(row$Pvalue_Burden, side$p_reported[k])
  }
})

test_that("plant -1: every test is negative and reads the #44 collapse in the results", {
  skip_unless_level("standard")
  r <- sg_run(env = c(SAIGE_SIGMA_PLANT_RATIO = "-1"))
  cnt <- sg_counts(r$log)
  expect_identical(cnt[2], cnt[1])
  expect_identical(cnt[3], cnt[1])
  side <- sg_side(r$out)
  expect_true(all(side$negative == "TRUE"))
  expect_true(all(as.numeric(side$var_grm) < 0))
  p <- sg_results(r$out)
  # the #44 collapse on every flagged burden cell: BETA and SE print NA, and the
  # p prints 1 unless the cell's p comes from the exact test (which reads no
  # variance); at least one cell reads the collapsed p
  n_one <- 0L
  for (k in which(side$Kind == "burden")) {
    row <- sg_row(p, side[k, ])
    expect_identical(nrow(row), 1L)
    expect_identical(row$BETA_Burden, "NA")
    expect_identical(row$SE_Burden, "NA")
    expect_identical(row$Pvalue_Burden, side$p_reported[k])
    if (as.numeric(row$Pvalue_Burden) == 1) n_one <- n_one + 1L
  }
  expect_gt(n_one, 0L)
})

test_that("plant 0.1 with APPLY=1: flagged tests take the diagonal variance and report it", {
  skip_unless_level("standard")
  obs <- sg_run(env = c(SAIGE_SIGMA_PLANT_RATIO = "0.1"))
  r <- sg_run(env = c(SAIGE_SIGMA_PLANT_RATIO = "0.1", SAIGE_SIGMA_GUARD_APPLY = "1"))
  cnt <- sg_counts(r$log)
  expect_identical(cnt[2], cnt[1])
  expect_identical(cnt[4], cnt[2])          # every flagged test applied
  expect_match(sg_line(r$log), "SAIGE_SIGMA_GUARD_APPLY=1", fixed = TRUE)
  side <- sg_side(r$out)
  expect_true(all(side$reported == "diagonal"))
  # the reported p is no longer the planted one: it is the diagonal variance's,
  # and the sidecar's p_reported is what the results table carries
  p <- sg_results(r$out)
  o <- sg_results(obs$out)
  cells <- which(p$Group != "Cauchy")
  expect_gt(sum(p$Pvalue_Burden[cells] != o$Pvalue_Burden[cells]), 0)
  for (k in which(side$Kind == "burden")) {
    row <- sg_row(p, side[k, ])
    expect_identical(nrow(row), 1L)
    expect_identical(row$Pvalue_Burden, side$p_reported[k])
    expect_false(row$BETA_Burden == "NA")
  }
})

test_that("kill switch: golden bytes, nothing evaluated, the line says off", {
  skip_unless_level("standard")
  r <- sg_run(env = c(SAIGE_SIGMA_GUARD = "0"))
  expect_output_matches_golden(sg_gold, r$out)
  expect_match(sg_line(r$log), "guard off", fixed = TRUE)
  expect_false(file.exists(paste0(r$out, ".sigmaGuard.txt")))
})
