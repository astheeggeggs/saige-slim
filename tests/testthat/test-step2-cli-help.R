# LEDGER #94 FIXED: the --nThreads help of step2_SPAtests.R says what the code
# does. It said "Ignored for group tests" (the threads drive the SKAT-O
# post-processing there, R/SkatoParallel.R) and that the one-thread fallback
# covers "PGEN, or PLINK without mmap" (PgenReader::clone() landed in S8; only
# the PLINK pread fallback cannot clone). The README documented the staleness
# in a second place; this pins the string itself so a third copy cannot drift.
#
# Tier: standard (one --help invocation, ~1.5 s).

test_that("LEDGER #94 FIXED: --nThreads help describes group tests and PGEN correctly", {
  skip_unless_level("standard")
  log <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", "--help", log = log)
  expect_cli_ok(st)
  help <- readLines(log, warn = FALSE)
  i <- grep("^\\s*--nThreads=", help)
  expect_length(i, 1L)
  para <- help[i + 1L]
  expect_false(grepl("Ignored for group tests", para, fixed = TRUE))
  expect_false(grepl("PGEN, or PLINK without mmap", para, fixed = TRUE))
  expect_true(grepl("SKAT-O post-processing", para, fixed = TRUE))
  expect_true(grepl("PGEN threads", para, fixed = TRUE))
  expect_true(grepl("could not be memory-mapped", para, fixed = TRUE))
  # 2026-09-06: threads past the physical core count cost CPU and buy no wall
  # time on a binary scan (harness/bench/NOTES.md); the help says so.
  expect_true(grepl("PHYSICAL core count", para, fixed = TRUE))
})
