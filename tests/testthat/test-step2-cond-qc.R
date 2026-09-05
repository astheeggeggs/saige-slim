# LEDGER #73 FIXED: conditioning markers bypassed the marker QC.
#
# Every other getMarkerPval call site sits behind the missingness/MAF/MAC/
# INFO gate; assign_conditionMarkers_factors had none. A marker the scan
# itself drops could be named to --condition and was used anyway: rs1 of the
# bundled fixture is monomorphic in the analysed samples, and --condition=rs1
# ran to exit 0 with the conditional columns IDENTICAL to the marginal ones on
# 70/70 rows -- a conditional analysis that did nothing while labelled as if
# it had (the only notice was a std::cerr line). It is refused now, naming the
# marker and the thresholds. Tier: standard (three marker runs).

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
  expect_true(any(grepl("conditioning marker rs1 .*fails the marker QC", r$log)))
  expect_true(any(grepl("it is monomorphic", r$log, fixed = TRUE)))
  expect_false(any(grepl("Analysis done", r$log, fixed = TRUE)))
})

test_that("the gate is the scan's own: a MAC-1 marker is refused at --minMAC=5 and accepted at 1", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # rs28 has MAC 1 in the analysed samples (golden P1: AC_Allele2 / AF_Allele2)
  refused <- cond_run("rs28", 5)
  expect_false(identical(refused$status, 0L))
  expect_true(any(grepl("conditioning marker rs28 .*fails the marker QC.*MAC 1", refused$log)))
  expect_false(any(grepl("monomorphic", refused$log, fixed = TRUE)))

  # POSITIVE CONTROL: the same marker at the threshold the scan uses in P1
  ok <- cond_run("rs28", 1)
  expect_identical(ok$status, 0L)
  d <- read_saige_table(ok$out)
  expect_true("p.value_c" %in% names(d))
  # and the conditioning did something: conditional variance moved somewhere
  expect_true(any(d$var_c != d$var))
})
