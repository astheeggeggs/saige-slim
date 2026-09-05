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
  owd <- setwd(saige_extdata)
  on.exit(setwd(owd), add = TRUE)
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
# ---- LEDGER #18, second half: is_single_in_groupTest -------------------------
# Golden P5 (plink_group_burden.txt) is a BURDEN-only group run at the CLI
# default, which writes NO .singleAssoc.txt. The API default used to be TRUE, so
# the same call from R wrote one. Both entry points now default to FALSE.

GRP_R <- c(
  PLINK_R,
  'GMMATmodelFile = "./output/example_binary_sparseGRM.rda"',
  'varianceRatioFile = "./output/example_binary_sparseGRM.varianceRatio.txt"',
  paste0('sparseGRMFile = "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed',
         '.sparseGRM.mtx"'),
  paste0('sparseGRMSampleIDFile = "output/sparseGRM_relatednessCutoff_0.125_1000_',
         'randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"'),
  'chrom = "1"', "min_MAF = 0", "min_MAC = 0.5", "LOCO = FALSE",
  'groupFile = "./input/group_new_snpid.txt"',
  'annotation_in_groupTest = c("lof", "missense;lof", "missense;lof;synonymous")',
  "maxMAF_in_groupTest = c(0.0001, 0.001, 0.01)",
  "r.corr = 1", "is_fastTest = FALSE"
)

test_that("LEDGER #18 FIXED: the API's BURDEN-only default writes no single-variant file", {
  out <- file.path(tempfile("api18_"), "burden.txt")
  dir.create(dirname(out), recursive = TRUE)
  api_run(out, GRP_R)
  # the region file is the golden's, so this IS P5's run through the other door
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_group_burden.txt"), out)
  expect_false(file.exists(paste0(out, ".singleAssoc.txt")))
})

# ---- LEDGER #82: is_noadjCov -------------------------------------------------
# Golden P1 (plink_marker_binary.txt) is the CLI default, is_noadjCov = TRUE,
# the unadjusted score test with the null_noXadj variance ratio. The API used
# to default to FALSE, the covariate-adjusted test, which is golden P30
# (plink_marker_noadjCov.txt): 68 of P1's 70 p-values differ from P30's. So the
# API called with P1's arguments must reproduce P1, and with is_noadjCov =
# FALSE it must reproduce P30 (positive control, true before and after).

MARKER_R <- c(
  PLINK_R,
  sprintf('GMMATmodelFile = "%s"', file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  sprintf('varianceRatioFile = "%s"',
          file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
  'chrom = "1"', "min_MAF = 0", "min_MAC = 1", "LOCO = FALSE",
  "is_output_moreDetails = TRUE"
)

test_that("LEDGER #82 FIXED: the API's is_noadjCov default reproduces golden P1", {
  out <- file.path(tempfile("api82_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  api_run(out, MARKER_R)
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_marker_binary.txt"), out)
})

test_that("LEDGER #82: is_noadjCov = FALSE from the API reproduces golden P30 (control)", {
  out <- file.path(tempfile("api82adj_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  api_run(out, MARKER_R, "call$is_noadjCov <- FALSE")
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_marker_noadjCov.txt"), out)
  # and the two goldens really are different tests, so the pair above has power
  p1 <- read_saige_table(file.path(saige_goldens, "step2", "plink_marker_binary.txt"))
  p30 <- read_saige_table(file.path(saige_goldens, "step2", "plink_marker_noadjCov.txt"))
  expect_gt(sum(p1$p.value != p30$p.value), 60)
})

test_that("LEDGER #18: is_single_in_groupTest = TRUE still writes the file (positive control)", {
  out <- file.path(tempfile("api18on_"), "burden.txt")
  dir.create(dirname(out), recursive = TRUE)
  api_run(out, GRP_R, "call$is_single_in_groupTest <- TRUE")
  sa <- paste0(out, ".singleAssoc.txt")
  expect_true(file.exists(sa))
  expect_gt(nrow(read_saige_table(sa)), 0)
})
