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

# LEDGER #111 FIXED: the --LOCO help of step2_SPAtests.R stated the OPPOSITE of
# its own default -- declared `default = TRUE`, help ended "[default=FALSE]".
# That is not cosmetic: requesting LOCO against a model that has none is a hard
# error (#112), so a user who believes the help gets a failed run rather than a
# silently different answer. The default itself is unchanged; only the help was
# wrong. Pinned here so the two cannot drift apart again.
#
# Tier: standard (one --help invocation, ~1.5 s).

test_that("LEDGER #111 FIXED: --LOCO help states its actual default, TRUE", {
  skip_unless_level("standard")
  log <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", "--help", log = log)
  expect_cli_ok(st)
  help <- readLines(log, warn = FALSE)
  i <- grep("^\\s*--LOCO=", help)
  expect_length(i, 1L)
  # the paragraph wraps, so take everything up to the next option line
  nxt <- grep("^\\s*--", help)
  nxt <- nxt[nxt > i]
  para <- paste(help[(i + 1L):(if (length(nxt)) nxt[1] - 1L else length(help))],
                collapse = " ")
  expect_false(grepl("[default=FALSE]", para, fixed = TRUE))
  expect_true(grepl("[default=TRUE]", para, fixed = TRUE))
  # the help also names the second requirement the option carries, which the
  # user cannot discover anywhere else: the model must have been fitted with
  # LOCO, or this run is a hard error (#112)
  expect_true(grepl("fitted with LOCO", para, fixed = TRUE))
  # The declared `default = TRUE` in the source is pinned by the fast-tier
  # check in test-step2-loco.R, so the two are not asserted twice here.
})
