# LEDGER #102: setStartEndIndex looped `i` over geno.M -- every marker in the
# plink file -- while indexing geno.alleleFreqVec, which is sized to the
# QC-passing count, so a LOCO fit with any marker under --minMAFforGRM read past
# the end of the allocation. The value it computed (geno.Msub) was dead, so no
# golden ever moved; what it did was crash the process at random, once in eight
# standard-tier runs and never in twelve isolated ones.
#
# A repeat count is therefore not a test of this. The characterization is a
# page-granular allocator: harness/guardmalloc.sh runs the same two-chromosome
# LOCO fit under /usr/lib/libgmalloc.dylib, which gives every allocation its own
# page followed by an unmapped one, so the read faults deterministically on the
# instruction that performs it. Measured 2026-09-07: red at the unpatched build
# (SIGSEGV, traceback frame 1 `setStartEndIndex`), green after.
#
# NOT a bounds-checked armadillo build: the access is alleleFreqVec[i], and
# armadillo's operator[] is unchecked whether or not ARMA_NO_DEBUG is defined.
# Dropping ARMA_NO_DEBUG catches nothing here.
#
# The script's own plant arm carries the detection-power control, and it is not
# optional -- Guard Malloc's standard failure mode is to print its banner while
# guarding nothing (see the header of harness/guardmalloc.sh). Tier: standard.

test_that("a step-1 LOCO fit reads nothing past the end of the heap", {
  skip_unless_level("standard")
  if (Sys.info()[["sysname"]] != "Darwin") skip("Guard Malloc is macOS-only")
  if (!file.exists("/usr/lib/libgmalloc.dylib")) skip("no /usr/lib/libgmalloc.dylib")
  script <- file.path(saige_repo, "harness", "guardmalloc.sh")
  skip_if_not(file.exists(script), "harness/guardmalloc.sh missing")

  work <- tempfile("guardmalloc_")
  # the library this suite actually loaded, not whatever the pixi env holds
  lib <- dirname(system.file(package = "SAIGE"))
  log <- tempfile(fileext = ".log")
  envv <- c(paste0("SAIGE_GUARD_WORK=", work), paste0("SAIGE_GUARD_LIB=", lib))
  status <- suppressWarnings(system2("bash", c(script, "both"), stdout = log,
                                     stderr = log, env = envv))
  txt <- paste(readLines(log, warn = FALSE), collapse = "\n")

  # the three arms, named, so a green run says which claim it is making
  expect_match(txt, "clean fit, guard ON\\s+status=0\\s+PASS")
  expect_match(txt, "planted OOB read, guard OFF\\s+status=0\\s+PASS")
  expect_match(txt, "planted OOB read, guard ON\\s+status=\\d+\\s+PASS")
  expect_identical(as.integer(status), 0L, info = txt)
})
