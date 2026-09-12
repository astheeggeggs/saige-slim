# Two region runs in ONE R process at DIFFERENT sample counts (2026-09-12).
#
# `CarrierStore` (src/Main.cpp, the 2026-09-12 carrier-space lever) is a FILE
# STATIC, so it survives from one `SPAGMMATtest` call to the next in a process,
# and `begin()` is what resets it. It used to reset the slot map on a changed
# sample count but NOT the three value buffers: the changed-n branch clears
# `carriers`, so the `used = carriers.size() * q` the else arm computes came
# back 0 and nothing was zeroed, while `ur` / `urNw` / `sum` still held the
# previous call's values in exactly the slots the next carriers were about to be
# handed (`slot()` hands out 0, 1, 2, ... again; `resize` only zeroes what it
# ADDS). Both accumulations assume a zero start -- the ultra-rare collapse is
# `u = std::max(u, ...)` and the burden sum is `+=` -- so the stale slots
# survived as maxima and were added in.
#
# It was a WRONG ANSWER, not merely undefined behaviour, and it was measured:
# on the pre-fix build this fixture's second run reported `GENE1 lof` at
# max_MAF 0.001 with MAC 7 and p 7.860290E-02 where the truth is MAC 2 and
# p 3.182528E-01, and `GENE1 lof` at 0.01 with MAC 69 against 65. Inflated MACs
# are the mechanism's own signature.
#
# CLAUDE.md says "one dataset per R process; no reentrancy", so this is at the
# edge of the contract -- but nothing REFUSES it, both calls return normally,
# and the numbers were silently wrong, which is the worst of the three
# possibilities. What this test pins is the property the fix delivers: a second
# run in one process is byte-identical to the same run standing alone. If the
# branch is ever made to refuse instead, this test is what should be inverted.
#
# The two in-process runs happen in a CHILD process (the pattern of
# test-step2-skat-failures.R): doing them here would leave this worker's own
# globals set by a foreign sample set for every later test in the file.

tr_cli <- run_saige_cli
tr_ok <- expect_cli_ok
tr_read <- read_saige_table

# helper-saige.R's objects, bound at top level: lintr's object-usage check
# cannot see another file's bindings from inside a function body
TR_EXT <- saige_extdata
TR_S1 <- file.path(TR_EXT, "output", "example_binary_sparseGRM")
TR_SG <- file.path(TR_EXT, "output",
                   "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")

# the model's own sample IDs -- 1,000 of the fam file's 10,000
tr_model_ids <- function() {
  e <- new.env()
  nm <- load(paste0(TR_S1, ".rda"), envir = e)
  as.character(get(nm[1], envir = e)$sampleID)
}

# One child process, two SPAGMMATtest calls: the full sample set then a subset,
# so CarrierStore::begin() is re-entered with a DIFFERENT n and the SAME cell
# count. Returns the two output paths.
tr_two_in_one <- function(subFile) {
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressMessages(library(SAIGE))",
    "args <- commandArgs(trailingOnly = TRUE)",
    "common <- list(",
    "  bedFile = './input/genotype_100markers.bed',",
    "  bimFile = './input/genotype_100markers.bim',",
    "  famFile = './input/genotype_100markers.fam',",
    "  AlleleOrder = 'alt-first',",
    "  GMMATmodelFile = args[3], varianceRatioFile = args[4],",
    "  sparseGRMFile = args[5], sparseGRMSampleIDFile = args[6],",
    "  chrom = '1', min_MAF = 0, min_MAC = 0.5, LOCO = FALSE,",
    "  groupFile = './input/group_new_snpid.txt',",
    "  annotation_in_groupTest = c('lof', 'missense;lof', 'missense;lof;synonymous'),",
    "  maxMAF_in_groupTest = c(0.0001, 0.001, 0.01),",
    "  is_fastTest = FALSE, r.corr = 1)",
    "# the full sample set first, then the subset: n changes between the calls",
    "do.call(SPAGMMATtest, c(common, list(SAIGEOutputFile = args[1])))",
    "do.call(SPAGMMATtest, c(common, list(SAIGEOutputFile = args[2],",
    "                                     subSampleFile = args[7])))"
  ), script)
  outA <- tempfile(pattern = "tworun_A_", fileext = ".txt")
  outB <- tempfile(pattern = "tworun_B_", fileext = ".txt")
  log <- tempfile(fileext = ".log")
  owd <- setwd(TR_EXT)
  on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(script, outA, outB, paste0(TR_S1, ".rda"), paste0(TR_S1, ".varianceRatio.txt"),
      TR_SG, paste0(TR_SG, ".sampleIDs.txt"), subFile),
    stdout = log, stderr = log
  ))
  list(status = st, A = outA, B = outB, log = readLines(log, warn = FALSE))
}

# the same run standing alone, through the CLI
tr_alone <- function(subFile) {
  out <- tempfile(pattern = "tworun_alone_", fileext = ".txt")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam",
            "--AlleleOrder=alt-first",
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", TR_S1, ".rda"),
            paste0("--varianceRatioFile=", TR_S1, ".varianceRatio.txt"),
            paste0("--sparseGRMFile=", TR_SG),
            paste0("--sparseGRMSampleIDFile=", TR_SG, ".sampleIDs.txt"),
            "--groupFile=./input/group_new_snpid.txt",
            "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
            "--maxMAF_in_groupTest=0.0001,0.001,0.01",
            "--is_fastTest=FALSE", "--r.corr=1",
            paste0("--SAIGEOutputFile=", out))
  if (nzchar(subFile)) args <- c(args, paste0("--subSampleFile=", subFile))
  st <- tr_cli("step2_SPAtests.R", args)
  list(status = st, path = out)
}

test_that("a second region run in one process is byte-identical to running it alone", {
  skip_unless_level("standard")
  ids <- tr_model_ids()
  expect_identical(length(ids), 1000L)
  sub <- tempfile(fileext = ".txt")
  writeLines(head(ids, 950), sub)      # 5% dropped: n changes, the cell grid does not

  two <- tr_two_in_one(sub)
  tr_ok(two$status)
  expect_true(file.exists(two$A))
  expect_true(file.exists(two$B))

  aloneA <- tr_alone("")
  tr_ok(aloneA$status)
  aloneB <- tr_alone(sub)
  tr_ok(aloneB$status)

  # The FIRST in-process run is the control: it shares its starting state with a
  # fresh process, so it was byte-identical even on the pre-fix build. If this
  # ever fails, the fixture or the flags moved and the real assertion below is
  # not measuring what it claims.
  expect_identical(readLines(two$A), readLines(aloneA$path))

  # The SECOND is the regression: it re-enters CarrierStore::begin() with a
  # changed n. Pre-fix this differed at line 2 (MAC 7 against 2).
  expect_identical(readLines(two$B), readLines(aloneB$path))

  # and the numbers are the subset's own, not the full set's -- a run that
  # silently ignored subSampleFile would pass the line above
  a <- tr_read(aloneA$path)
  b <- tr_read(two$B)
  expect_false(identical(readLines(aloneA$path), readLines(aloneB$path)))
  expect_identical(nrow(a), nrow(b))
})
