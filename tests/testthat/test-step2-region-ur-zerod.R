# LEDGER #129 (characterization, OPEN): the dosage-zeroing filter is applied to
# the WEIGHTED ultra-rare pseudo-marker in weight units. src/Main.cpp passes the
# collapsed `genoURVec` (max of w * G over the set's ultra-rare markers) through
# imputeGenoAndFlip with --dosage_zerod_cutoff (0.2) and
# --dosage_zerod_MAC_cutoff (10); src/UTIL.cpp cleans every entry <= 0.2 when
# the weighted MAC is <= 10. A custom or transform set whose support is all
# ultra-rare and whose weights put every carrier at or below 0.2 therefore loses
# its pseudo-marker and its row -- silently: not refused, not a SKAT failure,
# the pooled row reads NA and the #127 sidecar 0 / 0 / none. The same set with
# weights 0.5 is tested, and at collapse 0 (no pseudo-marker) both are tested
# with the SAME p, because a burden p is invariant to the weight scale.
# Fixture: GENET var rs28 rs32 (two singletons, cell MAC 2), one custom set.
# These pin what the code does TODAY; invert them with the fix, same commit.
# Tier: standard (four step-2 runs of a two-marker region).

urzerod_run <- function(w, collapse) {
  gf <- tempfile(fileext = ".txt")
  writeLines(c("GENET var rs28 rs32", "GENET anno lof lof",
               sprintf("GENET weight:w %s %s", w, w)), gf)
  out <- tempfile(fileext = ".txt")
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE", "--is_fastTest=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
            paste0("--sparseGRMFile=", sg),
            paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            paste0("--groupFile=", gf), "--annotation_in_groupTest=lof",
            "--maxMAF_in_groupTest=0.01", "--r.corr=1",
            paste0("--MACCutoff_to_CollapseUltraRare=", collapse),
            paste0("--SAIGEOutputFile=", out))
  st <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st)
  x <- read_saige_table(out)
  # the #127 sidecar is written only when a region has something to disclose
  sf <- paste0(out, ".pooledTests.txt")
  side <- if (file.exists(sf)) read.table(sf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  list(cells = x[x$Group != "Cauchy", ], pooled = x[x$Group == "Cauchy", ], side = side)
}

test_that("weights 0.01: the all-ultra-rare set vanishes at collapse 10, is tested at collapse 0", {
  skip_unless_level("standard")
  r10 <- urzerod_run(0.01, 10)
  expect_equal(nrow(r10$cells), 0L)
  expect_equal(r10$pooled$Pvalue_Burden, "NA")   # read_saige_table keeps the printed "NA"
  expect_equal(r10$side$Rows, 0L)
  expect_equal(r10$side$Cells_refused_minGroupMAC, 0L)   # not refused: it vanished
  expect_equal(r10$side$Rule, "none")
  r0 <- urzerod_run(0.01, 0)
  expect_equal(nrow(r0$cells), 1L)
  expect_equal(r0$cells$MAC, "2")
  expect_true(is.finite(as.numeric(r0$cells$Pvalue_Burden)))
  expect_null(r0$side)   # nothing to disclose once the set is tested
})

test_that("weights 0.5: the same set is tested at both collapses, with the same burden p", {
  skip_unless_level("standard")
  r10 <- urzerod_run(0.5, 10)
  r0 <- urzerod_run(0.5, 0)
  expect_equal(nrow(r10$cells), 1L)
  expect_equal(nrow(r0$cells), 1L)
  expect_equal(r10$cells$MAC, "2")
  # the weight scale cannot move a burden p: 0.5 at collapse 0 equals 0.01 at
  # collapse 0 (same path, inputs scaled), and the collapsed pseudo-marker of
  # two singletons equals their sum (no sample carries both) -- two different
  # code paths, so compared at printed precision, not for bit-identity
  p001 <- as.numeric(urzerod_run(0.01, 0)$cells$Pvalue_Burden)
  expect_equal(as.numeric(r0$cells$Pvalue_Burden), p001, tolerance = 1e-6)
  expect_equal(as.numeric(r10$cells$Pvalue_Burden), p001, tolerance = 1e-6)
})
