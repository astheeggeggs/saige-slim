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

write_plink_fixture_b <- function(prefix, geno) {
  m <- nrow(geno)
  n <- ncol(geno)
  writeLines(sprintf("1\tsnp%d\t0\t%d\tA\tT", seq_len(m), seq_len(m)),
             paste0(prefix, ".bim"))
  writeLines(sprintf("f%d\ts%d\t0\t0\t0\t-9", seq_len(n), seq_len(n)),
             paste0(prefix, ".fam"))
  code <- c(3L, 2L, 0L)
  nb <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (i in seq_len(m)) {
    cd <- c(code[geno[i, ] + 1L], rep(0L, nb * 4L - n))
    bytes <- vapply(seq_len(nb), function(b) {
      u <- cd[(4L * (b - 1L) + 1L):(4L * b)]
      as.raw(u[1] + 4L * u[2] + 16L * u[3] + 64L * u[4])
    }, raw(1))
    writeBin(bytes, con)
  }
  invisible(prefix)
}

# Haplotype fixture: founders, then `npo` parent-offspring pairs and `ndup`
# duplicates appended. Every marker polymorphic.
make_related_geno <- function(m, nfounder, npo, ndup, seed) {
  set.seed(seed)
  p <- stats::runif(m, 0.05, 0.5)
  h1 <- matrix(stats::rbinom(m * nfounder, 1, p), m, nfounder)
  h2 <- matrix(stats::rbinom(m * nfounder, 1, p), m, nfounder)
  par <- seq_len(npo)
  pick <- matrix(stats::runif(m * npo) < 0.5, m, npo)
  kid1 <- ifelse(pick, h1[, par, drop = FALSE], h2[, par, drop = FALSE])
  kid2 <- matrix(stats::rbinom(m * npo, 1, p), m, npo)
  dup <- seq_len(ndup) + npo
  g <- cbind(h1 + h2, kid1 + kid2, (h1 + h2)[, dup, drop = FALSE])
  for (i in seq_len(m)) {
    if (sum(g[i, ]) == 0) g[i, 1] <- 1L
    if (sum(g[i, ]) == 2 * ncol(g)) g[i, 1] <- 1L
  }
  list(g = g,
       planted = c(paste(par, nfounder + par, sep = ":"),
                   paste(dup, nfounder + npo + seq_len(ndup), sep = ":")))
}

# One child process per screen: setgeno and both screens write C++ globals.
# Returns the screened pairs as sorted "i:j" keys (1-based, as refineKin
# returns them) plus the log.
run_screen_b <- function(prefix, n, method, block = 2048L, cutoff = -1e6,
                         threads = 1L) {
  out <- tempfile(fileext = ".rds")
  log <- tempfile(fileext = ".log")
  script <- tempfile(fileext = ".R")
  call <- if (method == "blocked") {
    sprintf("findIndiceRelatedSampleBlocked(%dL)", as.integer(block))
  } else {
    "findIndiceRelatedSample()"
  }
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('prefix <- "%s"; out <- "%s"; n <- %d; cutoff <- %g',
            prefix, out, n, cutoff),
    sprintf("RcppParallel::setThreadOptions(numThreads = %d)", as.integer(threads)),
    "setminMAFforGRM(0.01); setmaxMissingRateforGRM(0.15)",
    'setgeno(paste0(prefix, ".bed"), paste0(prefix, ".bim"),',
    '        paste0(prefix, ".fam"), seq_len(n), rep(TRUE, n), 2, FALSE)',
    "mq <- length(getAlleleFreqVec())",
    "setRelatednessCutoff(cutoff)",
    "setSubMarkerIndex(seq_len(mq) - 1L)",
    "Get_MultiMarkersBySample_StdGeno_Mat()",
    call,
    # refineKin at -1e30 returns every screened pair; refineKin itself is
    # unchanged and shared by both screens.
    "kl <- refineKin(-1e30)",
    "saveRDS(list(i = kl$iIndex, j = kl$jIndex, k = kl$kinValue), out)"
  ), script)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = log, stderr = log))
  testthat::expect_identical(as.integer(st), 0L)
  res <- readRDS(out)
  o <- order(res$i, res$j)
  list(i = res$i[o], j = res$j[o], k = res$k[o],
       keys = paste(res$i[o], res$j[o], sep = ":"),
       log = readLines(log, warn = FALSE))
}

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
