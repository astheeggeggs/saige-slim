# A SECOND setgeno in one process must produce exactly the state a fresh
# process would (LEDGER #55).
#
# No shipped path re-ingests -- createSparseGRM is called with isSetGeno=FALSE
# from step 1, and getNeff runs once per process -- so this was latent until the
# multi-phenotype driver started reusing one process for several fits. The fix
# is therefore byte-neutral for every existing golden BY CONSTRUCTION, and this
# test is the only thing that can see it at all.
#
# Driven in a child process: the genotype store is a C++ global, so a run of
# this in the test session would leave it loaded for every later test.

reingest_probe <- function(lines) {
  out <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('inp <- "%s"; N_FILE <- %d', file.path(saige_extdata, "input"), N_FILE),
    sprintf('out <- "%s"', out),
    "ingest <- function(pfx, n, seed) {",
    "  base <- file.path(inp, pfx)",
    "  set.seed(seed)",
    "  setminMAFforGRM(0.01); setmaxMissingRateforGRM(0.15)",
    "  setminMAC_VarianceRatio(10.5, 20.5, TRUE)",
    # indicatorGenoSamplesWithPheno is indexed by FILE sample, not by kept
    # sample: it must have one entry per row of the .fam, TRUE where that sample
    # is in the fit. A short vector is read out of bounds in C++, which is UB
    # and -- being UB -- reads different memory in a fresh process than in one
    # that has already ingested, so it looks exactly like the leak under test.
    "  idx <- seq_len(n); ind <- c(rep(TRUE, n), rep(FALSE, N_FILE - n))",
    '  setgeno(paste0(base, ".bed"), paste0(base, ".bim"),',
    '          paste0(base, ".fam"), idx, ind, 2, FALSE)',
    "}",
    "snapshot <- function() list(",
    "  M = gettotalMarker(), N = getNnomissingOut(),",
    "  msub = getMsub_MAFge_minMAFtoConstructGRM(),",
    "  qcd = getQCdMarkerIndex(),",
    "  freq = getAlleleFreqVec(), mac = getMACVec(),",
    "  macvr = getMACVec_forVarRatio(), idxvr = getIndexVec_forVarRatio(),",
    "  diag = get_GRMdiagVec(),",
    "  g = sapply(0:9, function(i) Get_OneSNP_Geno(i))",
    ")",
    lines,
    "saveRDS(res, out)"
  ), script)
  log <- tempfile(fileext = ".log")
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = log, stderr = log))
  if (!identical(as.integer(st), 0L)) {
    testthat::fail(paste0(
      "re-ingest probe exited ", as.integer(st), "; log tail:\n",
      paste(utils::tail(readLines(log, warn = FALSE), 15), collapse = "\n")
    ))
  }
  readRDS(out)
}

# ONE fixture, two SAMPLE SUBSETS. The GRM fixture is the only bundled file
# whose markers survive the MAF >= 0.01 GRM gate at all (genotype_100markers is
# entirely rare: 0 markers pass, and the GRM diagonal comes back NaN), so both
# arms use it and differ in how many samples they ingest. That still moves every
# accumulator the bug appends to -- Nnomissing, the QC flags, the allele
# frequencies, the variance-ratio staging -- and it additionally exercises the
# m_DiagStd guard, which keys on the sample count.
FIXTURE <- "plinkforGRM_1000samples_10kMarkers"   # 1000 samples, 10,000 markers
N_FILE <- 1000L   # rows in the .fam
N_ALL <- 1000L
N_SUB <- 900L

test_that("a fresh ingest is what a re-ingest after close produces", {
  skip_unless_level("standard")
  fresh <- reingest_probe(c(
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_SUB),
    "res <- snapshot()"
  ))
  reused <- reingest_probe(c(
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_ALL),
    "closeGenoFile_plink()",
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_SUB),
    "res <- snapshot()"
  ))
  # Every field, not a summary: the failure mode is an APPEND, so the counts are
  # what go wrong first and the values only afterwards.
  expect_identical(reused$M, fresh$M)
  expect_identical(reused$N, fresh$N)
  expect_identical(reused$msub, fresh$msub)
  expect_identical(length(reused$qcd), length(fresh$qcd))
  expect_identical(reused$qcd, fresh$qcd)
  expect_identical(reused$freq, fresh$freq)
  expect_identical(reused$mac, fresh$mac)
  expect_identical(reused$macvr, fresh$macvr)
  expect_identical(reused$idxvr, fresh$idxvr)
  expect_identical(reused$diag, fresh$diag)
  expect_identical(reused$g, fresh$g)
})

test_that("a re-ingest WITHOUT closing first is equally clean", {
  skip_unless_level("standard")
  # The reset lives at the top of setGenoObjFromSource rather than in
  # closeGenoFile_plink on purpose: a caller that forgets to close must not be
  # able to corrupt the next fit.
  fresh <- reingest_probe(c(
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_SUB),
    "res <- snapshot()"
  ))
  reused <- reingest_probe(c(
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_ALL),
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_SUB),
    "res <- snapshot()"
  ))
  expect_identical(length(reused$qcd), length(fresh$qcd))
  expect_identical(reused$diag, fresh$diag)
  expect_identical(reused$g, fresh$g)
  expect_identical(reused$freq, fresh$freq)
})

test_that("two ingests of the SAME sample set agree (the m_DiagStd guard)", {
  skip_unless_level("standard")
  # m_DiagStd is rebuilt only when its length differs from Nnomissing, so a
  # second ingest at the same sample count would have reused the first one's GRM
  # diagonal -- the one case the marker/sample counts cannot reveal.
  fresh <- reingest_probe(c(
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_SUB),
    "res <- snapshot()"
  ))
  twice <- reingest_probe(c(
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_ALL),
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_SUB),
    sprintf('ingest("%s", %d, 1)', FIXTURE, N_SUB),
    "res <- snapshot()"
  ))
  expect_identical(twice$diag, fresh$diag)
  expect_identical(twice$msub, fresh$msub)
  expect_identical(length(twice$qcd), length(fresh$qcd))
})
