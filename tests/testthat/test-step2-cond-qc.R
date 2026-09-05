# LEDGER #73 FIXED: conditioning markers bypassed the marker QC.
#
# Every other getMarkerPval call site sits behind the missingness/MAF/MAC/
# INFO gate; assign_conditionMarkers_factors had none. A marker the scan
# itself drops could be named to --condition and was used anyway: rs1 of the
# bundled fixture is monomorphic in the analysed samples, and --condition=rs1
# ran to exit 0 with the conditional columns IDENTICAL to the marginal ones on
# 70/70 rows -- a conditional analysis that did nothing while labelled as if
# it had (the only notice was a std::cerr line). A monomorphic marker is
# refused now, by name. One that merely fails the SCAN's thresholds is used
# with a warning: golden P4 conditions on rs13 (MAC 2) under --minMAC=20, and
# a first draft that refused it was caught by the pre-push gate.
# Tier: standard (three marker runs).

cond_run <- function(condition, minMAC) {
  out <- tempfile(fileext = ".txt")
  m <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  args <- c("--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0",
            paste0("--minMAC=", minMAC), "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
            "--is_output_moreDetails=TRUE", paste0("--condition=", condition))
  st <- run_saige_cli("step2_SPAtests.R", args)
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}

test_that("a monomorphic conditioning marker is refused, by name", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- cond_run("rs1", 1)
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("conditioning marker rs1 .*is monomorphic in the analysed samples", r$log)))
  expect_true(any(grepl("It cannot be conditioned on", r$log, fixed = TRUE)))
  expect_false(any(grepl("Analysis done", r$log, fixed = TRUE)))
})

test_that("a marker below the scan's own threshold is used with a WARNING, and silently above it", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # rs28 has MAC 1 in the analysed samples (golden P1: AC_Allele2 / AF_Allele2).
  # Conditioning on a marker rarer than the scan's --minMAC is legitimate --
  # golden P4 conditions on rs13 (MAC 2) under --minMAC=20 -- so it is NOT
  # refused; the scan's thresholds are reported and the marker is used.
  warned <- cond_run("rs28", 5)
  expect_identical(warned$status, 0L)
  expect_true(any(grepl("conditioning marker rs28 .*would be dropped by the scan's own QC.*MAC 1", warned$log)))
  expect_true(any(grepl("^Warning message", warned$log)))
  d5 <- read_saige_table(warned$out)
  expect_true("p.value_c" %in% names(d5))
  expect_true(any(d5$var_c != d5$var))          # the conditioning did something

  # at the threshold P1 uses, the same marker passes without a word
  quiet <- cond_run("rs28", 1)
  expect_identical(quiet$status, 0L)
  expect_false(any(grepl("would be dropped by the scan", quiet$log, fixed = TRUE)))
  expect_false(any(grepl("^Warning message", quiet$log)))
})
