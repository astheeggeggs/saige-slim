# The ER path's entry into the R API (S2 in docs/design/step2-threading-verdict.md).
#
# ComputeExact::Run seeds via SL_setseed(1), which in the R build is
# Rcpp::Environment("package:base") -> set.seed -> GetRNGstate (the LEDGER #2
# fix). It used to be called on EVERY ER marker; it is now called only when a
# resampling arm will actually draw, because the exact recursions consume no
# randomness. That change is invisible to every golden -- the outputs are
# byte-identical by construction -- so without this file nothing in the suite
# would notice it being reverted, which is precisely the failure mode CLAUDE.md
# warns about ("golden-identical is a statement about coverage, not a pass").
#
# Counting is done by tracing base::set.seed in a CHILD process that calls
# SPAGMMATtest directly: the C++ singletons are one-dataset-per-process, so this
# must not run inside the testthat worker.
#
# Tier: standard (two step-2 runs).

seed_count_run <- function(max_mac_for_er, outfile) {
  script <- tempfile(fileext = ".R")
  writeLines(c(
    'args <- commandArgs(trailingOnly = TRUE)',
    'suppressPackageStartupMessages(library(SAIGE))',
    'cnt <- new.env(); cnt$n <- 0L',
    'suppressMessages(trace(base::set.seed,',
    '  tracer = quote(cnt$n <- cnt$n + 1L), print = FALSE))',
    'sink(tempfile())',
    'SPAGMMATtest(bedFile = "./input/genotype_100markers.bed",',
    '  bimFile = "./input/genotype_100markers.bim",',
    '  famFile = "./input/genotype_100markers.fam",',
    '  AlleleOrder = "alt-first",',
    '  GMMATmodelFile = args[3], varianceRatioFile = args[4],',
    '  SAIGEOutputFile = args[2],',
    '  chrom = "1", min_MAF = 0, min_MAC = 1, LOCO = FALSE,',
    '  is_noadjCov = TRUE, is_output_moreDetails = TRUE,',
    '  max_MAC_use_ER = as.numeric(args[1]))',
    'sink()',
    'cat("SEEDCOUNT", cnt$n, "\n")'
  ), script)
  s1 <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  log <- tempfile(fileext = ".log")
  owd <- setwd(saige_extdata); on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(script, max_mac_for_er, outfile,
      paste0(s1, ".rda"), paste0(s1, ".varianceRatio.txt")),
    stdout = log, stderr = log))
  expect_identical(as.integer(st), 0L)
  out <- grep("^SEEDCOUNT ", readLines(log, warn = FALSE), value = TRUE)
  expect_length(out, 1L)
  as.integer(sub("^SEEDCOUNT ", "", out))
}

test_that("ER seeds the R RNG only on markers that actually resample", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  # Default cutoff: 2 ER markers on this fixture, both with k <= 2 carriers, so
  # every one of their IsExact entries is 1 and no arm draws. Before S2 this
  # was 2 calls into the R API; it must now be none at all.
  o1 <- tempfile(fileext = ".txt")
  expect_identical(seed_count_run(4, o1), 0L)

  # --max_MAC_for_ER=200 admits rs38 (MAC 102, ~98 carriers), whose n_total
  # exceeds NResampling so its IsExact entries are not all 1. Exactly one of
  # the three ER markers draws, so exactly one seeding must happen -- and it
  # must still happen, or the resampling draws stop being reproducible.
  o2 <- tempfile(fileext = ".txt")
  expect_identical(seed_count_run(200, o2), 1L)

  # Null control for the whole exercise: both runs must still reproduce their
  # goldens. This also pins that the direct SPAGMMATtest call and the CLI agree.
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_binary.txt"), o1)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_er.txt"), o2)
})
