# LEDGER #61: a --maxMAC_in_groupTest bin above the largest
# --maxMAF_in_groupTest was emitted as a byte-for-byte DUPLICATE of the bin
# below it, at exit 0, after the log said the cutoff "will be applied".
#
# Not reachable by any golden: no golden passes --maxMAC_in_groupTest at all,
# and the default c(0) skips the whole block. So the runs here are derived from
# golden P2's command with one flag added, and the detector is a second run
# that asks for the SAME bin the legitimate way -- as a MAF -- which is what
# the flag is documented to be equivalent to.
#
# The fixture makes the arithmetic concrete: the model has 1,000 samples, so
# --maxMAC_in_groupTest=100 implies MAF 100/2000 = 0.05, five times the largest
# --maxMAF_in_groupTest=0.01. Every marker in (0.01, 0.05] was therefore
# discarded by the C++ per-marker QC gate, which had been fixed from the
# UN-extended vector before the model was even loaded.

cliMM <- run_saige_cli
BEDMM <- file.path(saige_extdata, "input", "genotype_100markers")
GRPMM <- file.path(saige_extdata, "input", "group_new_snpid.txt")
readMM <- read_saige_table

group_run <- function(extra, tag) {
  out <- file.path(tempfile(paste0("step2_mm_", tag, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- cliMM("step2_SPAtests.R", c(
    paste0("--bedFile=", BEDMM, ".bed"), paste0("--bimFile=", BEDMM, ".bim"),
    paste0("--famFile=", BEDMM, ".fam"), "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
           "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    paste0("--groupFile=", GRPMM),
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--is_fastTest=TRUE", extra,
    paste0("--SAIGEOutputFile=", out)
  ))
  list(status = st, path = out)
}

bin_rows <- function(path, maxMAF) {
  d <- readMM(path)
  d[d$max_MAF == maxMAF, ]
}

test_that("a maxMAC-derived bin is actually applied, not a copy of the bin below", {
  skip_unless_level("standard")

  byMAC <- group_run(c("--maxMAF_in_groupTest=0.0001,0.001,0.01",
                       "--maxMAC_in_groupTest=100"), "bymac")
  expect_cli_ok(byMAC$status)

  # The same bin requested the way that always worked. This is the reference
  # answer, not a second implementation: after the fix both runs hand C++ the
  # identical cutoff vector c(1e-4, 1e-3, 0.01, 0.05).
  byMAF <- group_run("--maxMAF_in_groupTest=0.0001,0.001,0.01,0.05", "bymaf")
  expect_cli_ok(byMAF$status)

  a <- bin_rows(byMAC$path, "0.05")
  b <- bin_rows(byMAF$path, "0.05")
  expect_identical(nrow(a), 6L)          # 2 regions x 3 annotation sets
  expect_identical(nrow(a), nrow(b))
  expect_identical(a$Region, b$Region)
  expect_identical(a$Group, b$Group)
  expect_lt(rel_err(a$Pvalue, b$Pvalue, tol = 1e-10), 1e-10)
  expect_identical(a$MAC, b$MAC)
  expect_identical(a$Number_rare, b$Number_rare)

  # DETECTOR. Before the fix the 0.05 bin of the maxMAC run was every statistic
  # of its own 0.01 bin, unchanged -- so an assertion that only compared the two
  # runs' 0.01 bins would have passed throughout. The 0.05 bin has to contain
  # markers the 0.01 bin does not, or the test cannot fail for the right reason.
  small <- bin_rows(byMAC$path, "0.01")
  expect_identical(nrow(small), nrow(a))
  expect_true(all(as.integer(a$Number_rare) >= as.integer(small$Number_rare)))
  expect_true(any(as.integer(a$Number_rare) > as.integer(small$Number_rare)))
})
