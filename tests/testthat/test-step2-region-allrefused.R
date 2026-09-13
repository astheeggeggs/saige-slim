# LEDGER #128 (FIXED 2026-09-13): on the SKAT-O path a region with TWO or more
# cells, every one of them refused by --minGroupMAC_in_BurdenTest, USED TO KILL
# the run in regionFinish's cross-region rbind ("names do not match previous
# names") -- the pooled row of a region with no surviving cell was built by the
# degenerate `rbind(NULL, cctVec)` branch of region_postprocess as a name-less
# character matrix, which rbind.data.frame refused one region later.
#
# The fix builds NO pooled row in that branch, which is what this path ALREADY
# did for a region with ONE cell (second test below, unchanged and now the
# CONTROL for the first: a two-cell all-refused region must behave exactly like
# the one-cell one). Unreachable at the default gate 1 (#126: a surviving cell
# has MAC >= 1), so the gate is NAMED here, as test-step2-region-threads.R does.
#
# LEDGER #116 (FIXED 2026-09-13) then decided how an empty region is reported:
# BOTH paths write one "Region Cauchy NA ..." row for it, the row the BURDEN
# writer always emitted (third test, unchanged). The first two tests were
# inverted with it: the SKAT-O path now writes that row too, for a two-cell
# and a one-cell all-refused region alike, from R/RegionSchema.R.
#
# Fixture: GENET var rs28 rs32 (AC 1 each, cell MAC 2) beside GENE1 of
# group_new_snpid.txt, on the frozen example_binary_sparseGRM model.
# Tier: standard (three step-2 runs).

allrefused_group_file <- function() {
  gf <- tempfile(fileext = ".txt")
  lines <- readLines(file.path(saige_extdata, "input", "group_new_snpid.txt"))
  gene1 <- grep("^GENE1 ", lines, value = TRUE)
  writeLines(c(gene1, "GENET var rs28 rs32", "GENET anno lof lof"), gf)
  gf
}

allrefused_run <- function(anno, r_corr = NULL, gate = 5) {
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
            paste0("--groupFile=", allrefused_group_file()),
            paste0("--annotation_in_groupTest=", anno), "--maxMAF_in_groupTest=0.01",
            paste0("--minGroupMAC_in_BurdenTest=", gate),
            paste0("--SAIGEOutputFile=", out))
  if (!is.null(r_corr)) args <- c(args, paste0("--r.corr=", r_corr))
  st <- run_saige_cli("step2_SPAtests.R", args)
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}

test_that("SKAT-O: a region whose two cells are both refused by the gate no longer kills the run", {
  skip_unless_level("standard")
  r <- allrefused_run("lof,lof:missense")
  expect_cli_ok(r$status)
  expect_true(any(grepl("Analyzing Region GENET", r$log, fixed = TRUE)))
  # the crash signature is gone, and the run reaches the end
  expect_false(any(grepl("names do not match previous names", r$log, fixed = TRUE)))
  expect_true(any(grepl("Analysis done", r$log, fixed = TRUE)))
  # GENE1 still reports; GENET contributes exactly ONE row, the NA Cauchy
  # row the BURDEN path writes for it (LEDGER #116 FIXED; it used to be
  # absent, and the two paths disagreed)
  x <- read_saige_table(r$out)
  genet <- x[x$Region == "GENET", , drop = FALSE]
  expect_identical(nrow(genet), 1L)
  expect_identical(genet$Group, "Cauchy")
  expect_identical(genet$Pvalue, "NA")
  expect_identical(genet$Pvalue_Burden, "NA")
  # ... and the sidecar is what carries GENET, with BOTH cells recorded refused
  side <- read.table(paste0(r$out, ".pooledTests.txt"), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  expect_equal(side$Cells_refused_minGroupMAC[side$Region == "GENET"], 2L)
  expect_equal(side$Rows[side$Region == "GENET"], 0L)
})

test_that("SKAT-O: the same region with ONE refused cell writes the same NA Cauchy row", {
  skip_unless_level("standard")
  r <- allrefused_run("lof")
  expect_cli_ok(r$status)
  x <- read_saige_table(r$out)
  # a 1 x 1 configuration writes no Cauchy row for a TESTED region (GENE1 has
  # its one cell row only), and one NA Cauchy row for the untested one
  expect_identical(x$Region, c("GENE1", "GENET"))
  expect_identical(x$Group, c("lof", "Cauchy"))
  side <- read.table(paste0(r$out, ".pooledTests.txt"), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  expect_equal(side$Region, "GENET")
  expect_equal(side$Cells_refused_minGroupMAC, 1L)
  expect_equal(side$Rows, 0L)
})

test_that("BURDEN: the two-cell all-refused region survives with an NA Cauchy row", {
  skip_unless_level("standard")
  r <- allrefused_run("lof,lof:missense", r_corr = 1)
  expect_cli_ok(r$status)
  x <- read_saige_table(r$out)
  expect_equal(x$Region[x$Group != "Cauchy"], c("GENE1", "GENE1"))
  genet <- x[x$Region == "GENET", ]
  expect_equal(nrow(genet), 1L)
  expect_equal(genet$Group, "Cauchy")
  expect_equal(genet$Pvalue_Burden, "NA")   # read_saige_table keeps the printed "NA"
  side <- read.table(paste0(r$out, ".pooledTests.txt"), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  expect_equal(side$Cells_refused_minGroupMAC[side$Region == "GENET"], 2L)
})
