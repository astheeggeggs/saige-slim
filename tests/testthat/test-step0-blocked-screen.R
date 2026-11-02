# Step 0's opt-in blocked-GEMM relatedness screen (--screenMethod=blocked,
# findIndiceRelatedSampleBlocked in src/SAIGE_fitGLMM_fast.cpp) against the
# default pairwise screen (findIndiceRelatedSample).
#
# The claim under test is SET IDENTITY, by construction: the GEMM value is a
# prefilter with a rigorous float-rounding margin, and every pair that clears
# it is re-scored with the pairwise worker's own two lines (arma::dot, then
# / m in float). So the screened pair list must be the same SET as the
# pairwise one (its ORDER is nondeterministic in both at nThreads > 1, which is
# why sets are compared sorted), and the .mtx -- written from a column-sorted
# Matrix by refineKin's per-pair, order-independent values -- byte-identical.
#
# Controls:
#   - cutoff -1e6 keeps every pair: the blocked tiling must visit exactly the
#     C(N,2) pairs i < j, each once (the bijection the pairwise fold has, #63),
#     for block sizes that do and do not divide N, B = 1, and B > N;
#   - a fixture with planted duplicates and parent-offspring pairs at a real
#     cutoff: both screens keep the planted pairs (positive control) and the
#     same set;
#   - the CLI on the bundled fixture reproduces the committed step-0 golden.

# write_plink_fixture_b, make_related_geno and run_screen_b live in helper-step0.R
# (shared with test-step0-refine-pairmajor.R since 2026-11-02).

test_that("blocked screen visits exactly the C(N,2) pairs i < j, any block size", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  prefix <- file.path(dir, "blk")
  n <- 23L
  fx <- make_related_geno(50L, 20L, 2L, 1L, 7L)
  write_plink_fixture_b(prefix, fx$g)
  expect_identical(ncol(fx$g), n)
  allkeys <- {
    p <- utils::combn(n, 2)
    sort(paste(p[1, ], p[2, ], sep = ":"))
  }
  pw <- run_screen_b(prefix, n, "pairwise")
  expect_identical(sort(pw$keys), allkeys)
  for (b in c(1L, 4L, 5L, 22L, 23L, 64L)) {
    for (t in c(1L, 3L)) {
      bl <- run_screen_b(prefix, n, "blocked", block = b, threads = t)
      expect_true(all(bl$i < bl$j))
      expect_identical(anyDuplicated(bl$keys), 0L)
      expect_identical(sort(bl$keys), allkeys)
      # same pairs, same refined values (refineKin is order-independent per pair)
      expect_identical(bl$k, pw$k)
    }
  }
})

test_that("blocked screen keeps the same pair set as pairwise at a real cutoff", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  prefix <- file.path(dir, "blkrel")
  fx <- make_related_geno(400L, 150L, 6L, 4L, 11L)
  write_plink_fixture_b(prefix, fx$g)
  n <- ncol(fx$g)
  for (cut in c(0.05, 0.125)) {
    pw <- run_screen_b(prefix, n, "pairwise", cutoff = cut)
    # POSITIVE CONTROL: the planted PO (0.5) and duplicate (1.0) pairs are screened in
    expect_true(all(fx$planted %in% pw$keys))
    # and the cutoff actually filters (not every pair kept)
    expect_lt(length(pw$keys), choose(n, 2))
    for (b in c(7L, 32L, 2048L)) {
      bl <- run_screen_b(prefix, n, "blocked", block = b, cutoff = cut, threads = 2L)
      expect_identical(bl$keys, pw$keys)
      expect_identical(bl$k, pw$k)
      expect_true(any(grepl("^Blocked screen: .* pairs kept of", bl$log)))
    }
  }
})

test_that("--screenMethod=blocked reproduces the step-0 golden byte for byte", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  suffix <- "_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
  out <- file.path(tempdir(), "s0_blocked")
  st <- run_saige_cli("createSparseGRM.R", c(
    "--plinkFile=input/plinkforGRM_1000samples_10kMarkers",
    "--nThreads=2", paste0("--outputPrefix=", out),
    "--numRandomMarkerforSparseKin=5000", "--relatednessCutoff=0.05",
    "--screenMethod=blocked", "--screenBlockSize=300"
  ))
  expect_cli_ok(st)
  golden <- file.path(saige_goldens, "step0", "test_sparseGRM")
  expect_files_byte_identical(paste0(golden, suffix), paste0(out, suffix))
  expect_files_byte_identical(paste0(golden, suffix, ".sampleIDs.txt"),
                              paste0(out, suffix, ".sampleIDs.txt"))
})

test_that("screenMethod and screenBlockSize are validated", {
  expect_error(createSparseGRM(screenMethod = "gemm"), "screenMethod must be")
  expect_error(createSparseGRM(screenMethod = "blocked", screenBlockSize = 0),
               "screenBlockSize must be a positive integer")
  expect_error(createSparseGRM(screenMethod = "blocked", screenBlockSize = 2.5),
               "screenBlockSize must be a positive integer")
})
