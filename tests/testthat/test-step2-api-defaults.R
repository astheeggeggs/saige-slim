# The R API (SPAGMMATtest) and the command line (extdata/step2_SPAtests.R) are
# two entry points to one function, and every committed golden is a CLI run. A
# default that differs between them is therefore invisible to the goldens while
# every direct R caller gets different behaviour from what the goldens pin.
# LEDGER #18 and #82 were both of that shape. Each test here drives the API from
# a child Rscript with cwd = extdata/ (so the fixture paths are the goldens'
# own), passing ONLY the arguments the corresponding golden command passes, and
# compares against that golden or against the CLI run.
#
# Tier: standard (runs step 2).

api_run <- function(out, args_r, extra_r = "") {
  skip_unless_level("standard")
  skip_if_no_goldens()
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    "sink(tempfile())",
    "call <- list(", paste(args_r, collapse = ",\n"), ")",
    extra_r,
    sprintf('call$SAIGEOutputFile <- "%s"', out),
    "do.call(SPAGMMATtest, call)",
    "sink()"
  ), script)
  log <- tempfile(fileext = ".log")
  owd <- setwd(saige_extdata); on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = log, stderr = log))
  if (!identical(as.integer(st), 0L)) {
    testthat::fail(paste0("API child exited with status ", as.integer(st), "; log tail:\n",
                          paste(utils::tail(readLines(log, warn = FALSE), 15), collapse = "\n")))
  }
  invisible(out)
}

PLINK_R <- c(
  'bedFile = "./input/genotype_100markers.bed"',
  'bimFile = "./input/genotype_100markers.bim"',
  'famFile = "./input/genotype_100markers.fam"',
  'AlleleOrder = "alt-first"'
)
PLINK_CLI <- c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first"
)

# ---- LEDGER #18, second half: is_single_in_groupTest -------------------------
# Golden P5 (plink_group_burden.txt) is a BURDEN-only group run at the CLI
# default, which writes NO .singleAssoc.txt. The API default used to be TRUE, so
# the same call from R wrote one. Both entry points now default to FALSE.

GRP_R <- c(
  PLINK_R,
  'GMMATmodelFile = "./output/example_binary_sparseGRM.rda"',
  'varianceRatioFile = "./output/example_binary_sparseGRM.varianceRatio.txt"',
  'sparseGRMFile = "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"',
  'sparseGRMSampleIDFile = "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"',
  'chrom = "1"', "min_MAF = 0", "min_MAC = 0.5", "LOCO = FALSE",
  'groupFile = "./input/group_new_snpid.txt"',
  'annotation_in_groupTest = c("lof", "missense;lof", "missense;lof;synonymous")',
  "maxMAF_in_groupTest = c(0.0001, 0.001, 0.01)",
  "r.corr = 1", "is_fastTest = FALSE"
)

test_that("LEDGER #18 FIXED: the API's BURDEN-only default writes no single-variant file, like the CLI", {
  out <- file.path(tempfile("api18_"), "burden.txt")
  dir.create(dirname(out), recursive = TRUE)
  api_run(out, GRP_R)
  # the region file is the golden's, so this IS P5's run through the other door
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_group_burden.txt"), out)
  expect_false(file.exists(paste0(out, ".singleAssoc.txt")))
})

test_that("LEDGER #18: is_single_in_groupTest = TRUE still writes the file (positive control)", {
  out <- file.path(tempfile("api18on_"), "burden.txt")
  dir.create(dirname(out), recursive = TRUE)
  api_run(out, GRP_R, "call$is_single_in_groupTest <- TRUE")
  sa <- paste0(out, ".singleAssoc.txt")
  expect_true(file.exists(sa))
  expect_gt(nrow(read_saige_table(sa)), 0)
})
