# LEDGER #149 (FIXED 2026-09-14): getNeff's LOCO windows are computed over the
# QC-retained markers, by the shared helper, after ingest.
#
# The LOCO block used to be a stale inline copy of updateChrStartEndIndexVec
# run BEFORE ingest over the raw bim chromosome column, so its windows were
# in raw-bim space while every consumer of setStartEndIndexVec indexes the
# QC'd genotype store; it also lacked the LEDGER #62 contiguity guard and
# named `plinkFile.bim` in a message reachable from a PGEN run. Harmless
# only because nothing on the getNeff path reads the windows today (LEDGER
# #153: tau is fixed). Now the block is fitNULLGLMM's three lines after each
# setgenoFromInput: getQCdMarkerIndex(), subset the chromosome vector, call
# the helper.
#
# The observable is the printed windows. Their total span,
# sum(end - start + 1) over the chromosomes present, was the raw bim count
# (128,868 on the 22-chromosome fixture) and is now the QC'd count the same
# log prints ("N markers with MAF >= 0.01 and missing rate <= 0"; 111,305).
# `--LOCO` has no CLI flag in extractNglmm.R, so this goes through the R API
# in a child process (getNeff drives the one-dataset-per-process C++ state).
# Tier: standard (one getNeff run, ~30 s).
#
# 2026-09-15, LEDGER #153: getNeff always solves on the sparse GRM, whose
# branch sets LOCO <- FALSE, so the windows are no longer printed by getNeff
# at all; the test below pins that, and the helper is covered by step 1.

l149_extdata <- function() get("saige_extdata")

l149_run <- function() {
  log <- tempfile(fileext = ".log")
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  q <- function(x) paste0("\"", x, "\"")
  code <- paste0(
    "suppressMessages(library(SAIGE)); SAIGE::getNeff(",
    "plinkFile = ", q("./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr"), ", ",
    "sparseGRMFile = ", q(sg), ", sparseGRMSampleIDFile = ", q(paste0(sg, ".sampleIDs.txt")), ", ",
    "phenoFile = ", q("./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"), ", ",
    "phenoCol = ", q("y_binary"), ", covarColList = c(", q("x1"), ", ", q("x2"), "), ",
    "sampleIDColinphenoFile = ", q("IID"), ", traitType = ", q("binary"), ", ",
    "outputPrefix = ", q(tempfile("nglmm_loco_")), ", LOCO = TRUE)"
  )
  rscript <- file.path(R.home("bin"), "Rscript")
  owd <- setwd(l149_extdata())
  on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(rscript, c("-e", shQuote(code)), stdout = log, stderr = log))
  list(status = as.integer(st), log = readLines(log, warn = FALSE))
}

l149_vec <- function(lg, name) {
  ln <- grep(paste0("^", name, ": "), lg, value = TRUE)
  stopifnot(length(ln) == 1L)
  as.numeric(strsplit(trimws(sub(paste0("^", name, ": "), "", ln)), "\\s+")[[1]])
}

test_that("getNeff with LOCO = TRUE drops LOCO and says why (LEDGER #149 -> #153)", {
  skip_unless_level("standard")
  # Since LEDGER #153 (2026-09-15) every getNeff run solves on the sparse GRM,
  # and that branch has always set LOCO <- FALSE, so the QC'd-window block
  # this test was written for is unreachable from getNeff: no window vector is
  # printed. The helper it calls, updateChrStartEndIndexVec, is fitNULLGLMM's
  # and stays covered there (test-step1-lowmem-loco.R, test-step2-loco.R).
  r <- l149_run()
  expect_identical(r$status, 0L)
  expect_true(any(grepl("LOCO is not applied in getNeff", r$log, fixed = TRUE)))
  expect_false(any(grepl("^chromosomeStartIndexVec: ", r$log)))
  expect_false(any(grepl("^chromosomeEndIndexVec: ", r$log)))
  expect_length(grep("^Nglmm ", r$log), 1L)
})
