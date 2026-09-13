# LEDGER #135 FIXED: the BURDEN emitter gated a cell on a FOLDED pseudo-frequency
# of its aggregate allele count -- MAF = min(f, 1 - f) with f = count / 2N,
# clamped to N above f = 1 -- so the MAC that --minGroupMAC_in_BurdenTest, the
# variance-ratio band and the ER predicate read DECREASED as the cell's count
# grew past N. On P65's configuration (the bundled 100 markers at max MAF 0.5,
# collapse cutoff 700, in-sample N = 1,000) the four distinct cells print MAC
# 748, 3064, 1915 and 2997 and were gated as 748, 1000, 85 and 1000. The gating
# MAC is now the count itself, which is what the printed MAC already was.
#
# Two assertions, each red at unpatched HEAD:
#   * the gate reads the PRINTED count: at 1500 exactly the cells whose gate-1
#     row prints MAC >= 1500 survive, with their rows verbatim (before: none
#     survived, the file was two `Cauchy NA` rows);
#   * the two region paths refuse the SAME cells at a gate past N -- the SKAT-O
#     half of #126 always gated on MAC_GroupVec, the honest sum, so the two
#     halves disagreed exactly where the fold bit (before: BURDEN refused all
#     four cells, SKAT-O only GENE1 lof).
# Non-vacuity: the fixture must hold a tested cell with MORE alleles than N and
# one the gate refuses, or the regime the row is about is not reached.
# Tier: standard.

bam_cli <- run_saige_cli
bam_ok <- expect_cli_ok
bam_read <- read_saige_table

bam_run <- function(extra) {
  out <- file.path(tempfile("step2_bam_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- bam_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.5", "--MACCutoff_to_CollapseUltraRare=700",
    "--is_fastTest=FALSE",
    paste0("--SAIGEOutputFile=", out), extra
  ))
  bam_ok(status)
  list(path = out, rows = bam_read(out), pooled = paste0(out, ".pooledTests.txt"))
}
bam_cells <- function(r) r$rows[r$rows$Group != "Cauchy", , drop = FALSE]
bam_key <- function(d) paste(d$Region, d$Group, d$max_MAF, sep = "|")

test_that("LEDGER #135: the BURDEN gate reads the cell's aggregate count, not a folded frequency", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  g1 <- bam_run(c("--r.corr=1", "--minGroupMAC_in_BurdenTest=1"))
  g1500 <- bam_run(c("--r.corr=1", "--minGroupMAC_in_BurdenTest=1500"))
  c1 <- bam_cells(g1)
  c1500 <- bam_cells(g1500)
  mac1 <- as.numeric(c1$MAC)
  # the regime: a cell past N = 1,000 (the fold would gate it below N) and a
  # cell below the gate (so the gate is live)
  expect_gt(sum(mac1 > 1000), 0)
  expect_gt(sum(mac1 < 1500), 0)
  expect_gt(sum(mac1 >= 1500), 0)
  # membership is decided by the printed count ...
  expect_identical(bam_key(c1500), bam_key(c1[mac1 >= 1500, , drop = FALSE]))
  # ... and the surviving rows are the gate-1 rows verbatim: the gate changes
  # membership only, never a tested cell's statistic
  keep <- c1[mac1 >= 1500, , drop = FALSE]
  rownames(keep) <- NULL
  rownames(c1500) <- NULL
  expect_identical(c1500, keep)
  # the #127 sidecar discloses the refusals region by region; a region with
  # nothing to disclose writes no row at all
  expect_true(file.exists(g1500$pooled))
  sc <- bam_read(g1500$pooled)
  refused <- tapply(as.integer(sc$Cells_refused_minGroupMAC), sc$Region, sum)
  expected <- tapply(mac1 < 1500, c1$Region, sum)
  present <- names(expected)[expected > 0]
  expect_setequal(unique(sc$Region), present)
  expect_identical(as.integer(refused[present]), as.integer(expected[present]))
})

test_that("LEDGER #135/#126: both region paths refuse the same cells at a gate past N", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  b <- bam_run(c("--r.corr=1", "--minGroupMAC_in_BurdenTest=1500"))
  s <- bam_run(c("--r.corr=0", "--minGroupMAC_in_BurdenTest=1500"))
  # cell identity only, never the p-values: the two paths do not share a
  # burden estimator (LEDGER #126)
  expect_identical(bam_key(bam_cells(b)), bam_key(bam_cells(s)))
  expect_setequal(unique(bam_cells(b)$Region), c("GENE1", "GENE2"))
})
