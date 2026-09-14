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

test_that("getNeff's LOCO windows span the QC'd markers, not the raw bim (LEDGER #149)", {
  skip_unless_level("standard")
  r <- l149_run()
  expect_identical(r$status, 0L)
  starts <- l149_vec(r$log, "chromosomeStartIndexVec")
  ends <- l149_vec(r$log, "chromosomeEndIndexVec")
  expect_length(starts, 22L)
  expect_length(ends, 22L)
  ok <- !is.na(starts) & !is.na(ends)
  span <- sum(ends[ok] - starts[ok] + 1)
  qc_line <- grep("markers with MAF >= 0.01 and missing rate <= 0", r$log, value = TRUE)
  expect_length(qc_line, 1L)
  n_qc <- as.numeric(sub("^\\s*(\\d+) markers.*$", "\\1", qc_line))
  bim <- file.path(l149_extdata(), "input",
                   "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim")
  n_raw <- length(readLines(bim))
  # the fixture really has QC to apply
  expect_lt(n_qc, n_raw)
  # the inverted characterization: before the fix span == n_raw (128,868)
  expect_identical(span, n_qc)
  expect_false(span == n_raw)
  # the windows are the QC'd store's: contiguous, ordered, starting at 0
  expect_identical(starts[ok][1], 0)
  expect_true(all(diff(starts[ok]) > 0))
  expect_identical(starts[ok][-1], ends[ok][-sum(ok)] + 1)
})
