# The fastTest ESCALATION pass of a non-BURDEN region test (LEDGER #50).
#
# `--is_fastTest=TRUE` runs each region twice when the fast pass looks
# interesting: pass 1 with the cheap variance ratios, and -- if the region's
# Cauchy p is below `pval_cutoff_for_fastTest` -- pass 2 with the sparse GRM,
# which is the answer that gets kept. R signals pass 2 to C++ by
# set_flagSparseGRM_cur_SAIGE(TRUE) before the second mainRegionInCPP call.
#
# The single-variant-in-group rows have TWO possible routes to the output file
# and exactly one fires per region: on an unescalated region C++ parks them in
# `.singleAssoc.txt_temp` and R copies them across itself (copy_singleInGroup);
# on an escalated region C++ must write them directly, which it decides from
# that signal. Before the #50 fix it read the signal off a model member the
# marker and ultra-rare loops overwrite from each marker's MAC, so the escalated
# pass parked the rows in the temp bin too -- and the next region truncated it.
# Every escalated region's single-variant rows were silently lost.
#
# NO GOLDEN REACHED THIS WHEN #50 WAS FIXED, which is why it survived: at the
# default cutoff of 0.05 neither bundled region escalates, so P2 exercises only
# the copy_singleInGroup route and P6 (is_fastTest=FALSE) only the unconditional
# write, and `pval_cutoff_for_fastTest` was a `SPAGMMATtest()` parameter with no
# `--` option here or upstream. It IS a CLI option now, and golden P35 is this
# run frozen -- but the direct-call route below is KEPT, because it is the one
# that does not depend on the CLI plumbing being right, and because it is what
# proves the two agree. Both go through a child process: the C++ singletons are
# one-dataset-per-process (same reason as test-step2-er-seeding.R).
#
# Tier: standard (three step-2 runs).

escalation_run <- function(cutoff, outfile) {
  script <- tempfile(fileext = ".R")
  writeLines(c(
    'args <- commandArgs(trailingOnly = TRUE)',
    'suppressPackageStartupMessages(library(SAIGE))',
    'SPAGMMATtest(bedFile = "./input/genotype_100markers.bed",',
    '  bimFile = "./input/genotype_100markers.bim",',
    '  famFile = "./input/genotype_100markers.fam",',
    '  AlleleOrder = "alt-first",',
    '  GMMATmodelFile = args[3], varianceRatioFile = args[4],',
    '  sparseGRMFile = args[5], sparseGRMSampleIDFile = args[6],',
    '  SAIGEOutputFile = args[2],',
    '  chrom = "1", min_MAF = 0, min_MAC = 0.5, LOCO = FALSE,',
    '  groupFile = "./input/group_new_snpid.txt",',
    '  annotation_in_groupTest = c("lof", "missense;lof",',
    '                              "missense;lof;synonymous"),',
    '  maxMAF_in_groupTest = c(0.0001, 0.001, 0.01),',
    '  is_fastTest = TRUE,',
    '  pval_cutoff_for_fastTest = as.numeric(args[1]))'
  ), script)
  s1 <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  log <- tempfile(fileext = ".log")
  owd <- setwd(saige_extdata); on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(script, cutoff, outfile,
      paste0(s1, ".rda"), paste0(s1, ".varianceRatio.txt"),
      sg, paste0(sg, ".sampleIDs.txt")),
    stdout = log, stderr = log))
  expect_cli_ok(st)
  # How many regions actually took the escalated pass. Asserting this is what
  # separates "the escalated writer works" from "the run quietly never
  # escalated and wrote its rows the other way" -- the LEDGER #30 shape.
  sum(grepl("^Non-fast test is performed",
            readLines(log, warn = FALSE), fixed = FALSE))
}

test_that("an escalated non-BURDEN region still writes its single-variant rows", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  golden <- file.path(saige_goldens, "step2", "plink_group_fastTest.txt")

  # pval_cutoff_for_fastTest = 1 escalates every region whose Cauchy p is below
  # 1, i.e. both of them.
  esc <- tempfile(fileext = ".txt")
  expect_identical(escalation_run(1, esc), 2L)

  # The rows must be there at all. Before the fix this file was the header and
  # nothing else, which is the whole defect.
  esc_single <- paste0(esc, ".singleAssoc.txt")
  expect_true(file.exists(esc_single))
  expect_identical(nrow(read_saige_table(esc_single)), 19L)

  # And they must be the RIGHT rows. Escalated-fastTest and unescalated-fastTest
  # are two different control flows, so this compares labels exactly and
  # statistics at tolerance rather than asserting bit-identity between paths
  # (they do in fact agree byte for byte here -- escalation changes the route,
  # not the arithmetic -- but that is an observation, not something to pin).
  expect_tables_equivalent(paste0(golden, ".singleAssoc.txt"), esc_single,
                           tol = 1e-12)
  expect_tables_equivalent(golden, esc, tol = 1e-12)
})

test_that("an unescalated fastTest region still reproduces golden P2", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  # The negative control for the test above, and the reason the fix is known to
  # be inert everywhere a golden looks: at the default cutoff neither region
  # escalates, the rows come through copy_singleInGroup as before, and the run
  # matches the committed golden in whichever mode the suite is comparing in.
  # It also pins that this direct SPAGMMATtest call and P2's CLI run agree.
  ctrl <- tempfile(fileext = ".txt")
  expect_identical(escalation_run(0.05, ctrl), 0L)

  golden <- file.path(saige_goldens, "step2", "plink_group_fastTest.txt")
  expect_output_matches_golden(golden, ctrl)
  expect_output_matches_golden(paste0(golden, ".singleAssoc.txt"),
                              paste0(ctrl, ".singleAssoc.txt"))
})

test_that("--pval_cutoff_for_fastTest reaches the escalation path from the CLI", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  # What golden P35 buys: the escalated writer is now reachable by the scripts,
  # so the path has a frozen byte reference instead of only a live comparison.
  out <- file.path(tempfile("step2_p35_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- run_saige_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
           "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--is_fastTest=TRUE", "--pval_cutoff_for_fastTest=1",
    paste0("--SAIGEOutputFile=", out)
  ))
  expect_cli_ok(st)

  # The option has to have ARRIVED. Without this the test would pass on a run
  # that silently ignored the flag and never escalated -- the LEDGER #30 shape,
  # and the more so here because P35's bytes equal P2's (LEDGER #51: the region
  # path recomputes the same arithmetic on the second pass).
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_identical(sum(grepl("^Non-fast test is performed", log)), 2L)
  # And the run must have said the flag is pure cost on a group test.
  expect_true(any(grepl("no computational benefit for group tests", log,
                        fixed = TRUE)))

  g35 <- file.path(saige_goldens, "step2", "plink_group_fastTest_escalated.txt")
  expect_output_matches_golden(g35, out)
  expect_output_matches_golden(paste0(g35, ".singleAssoc.txt"),
                               paste0(out, ".singleAssoc.txt"))
})
