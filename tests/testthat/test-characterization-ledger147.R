# LEDGER #147 (FIXED 2026-09-14): --minGroupMAC_in_BurdenTest is validated by
# name, like every other numeric step-2 argument.
#
# Before the fix checkArgsListNumeric took it as a formal and checked nothing:
# optparse (type = "double") warned `double expected, got "abc"`, KEPT the
# string, and the run died ~200 lines later in an Rcpp conversion that named
# no argument (`Not compatible with requested type: [type=character;
# target=double]`, after "group-based test will be performed"); a negative
# value was accepted and the #126 gate went inert. Now a non-numeric or
# negative value stops in checkArgs with the argument's name, before any file
# is opened. Tier: standard (runs the CLI). The inverted characterization.

l147_cli <- run_saige_cli
l147_extdata <- function() get("saige_extdata")

l147_run <- function(value) {
  out <- file.path(tempfile("step2_l147_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  args <- c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--is_fastTest=TRUE",
    paste0("--minGroupMAC_in_BurdenTest=", value),
    paste0("--SAIGEOutputFile=", out)
  )
  status <- l147_cli("step2_SPAtests.R", args, log = log)
  list(status = as.integer(status), log = readLines(log, warn = FALSE), out = out)
}
l147_has <- function(r, pattern) any(grepl(pattern, r$log, fixed = TRUE))

test_that("a non-numeric --minGroupMAC_in_BurdenTest stops by name (LEDGER #147)", {
  skip_unless_level("standard")
  r <- l147_run("abc")
  expect_false(identical(r$status, 0L))
  expect_true(l147_has(r, "minGroupMAC_in_BurdenTest should be a numeric value"))
  # the pre-fix death: an unnamed Rcpp type error after the group test had started
  expect_false(l147_has(r, "Not compatible with requested type"))
  expect_false(l147_has(r, "group-based test will be performed"))
  expect_false(file.exists(r$out))
})

test_that("a negative --minGroupMAC_in_BurdenTest stops by name (LEDGER #147)", {
  skip_unless_level("standard")
  r <- l147_run("-1")
  expect_false(identical(r$status, 0L))
  expect_true(l147_has(r, paste("minGroupMAC_in_BurdenTest should be a numeric value",
                                "greater than or equal to 0")))
  expect_false(file.exists(r$out))
})
