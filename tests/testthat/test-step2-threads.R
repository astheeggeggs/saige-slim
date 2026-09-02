# --nThreads for the single-variant marker scan (Phase 4, C2; LEDGER #30/#31).
#
# The claim under test is byte-exactness of a threaded run against the EXISTING
# serial goldens. That is a legitimate byte claim, and not an instance of the
# CLAUDE.md "never assert bit-identity between two different code paths" trap:
# each shard runs the same serial code over a disjoint block of the same
# canonical marker order, no marker's arithmetic depends on another's, and the
# only cross-marker aggregates (chunk progress, numtest, mFirth) go to stdout
# and never reach the file. One code path, re-run over a partition.
#
# make_goldens.sh is deliberately NOT modified: an nThreads > 1 golden would
# defeat the whole gate.
#
# Tier: standard for the CLI runs; the partition/merge unit tests below are
# pure R and run in the fast tier.

# ---- unit: splitMarkerShards ------------------------------------------------

test_that("splitMarkerShards partitions contiguously with no empty block", {
  for (M in c(1L, 2L, 7L, 70L, 128868L)) {
    for (k in c(1L, 2L, 3L, 4L, 7L, 64L)) {
      blocks <- SAIGE:::splitMarkerShards(M, k)
      expect_equal(length(blocks), min(k, M))
      # every index exactly once, in order
      expect_identical(unlist(blocks), seq_len(M))
      # contiguous, non-empty, near-equal
      sizes <- vapply(blocks, length, integer(1))
      expect_true(all(sizes >= 1L))
      expect_lte(max(sizes) - min(sizes), 1L)
    }
  }
})

test_that("splitMarkerShards clamps to M rather than emitting empty blocks", {
  # parallel::splitIndices(3, 64) would give 61 empty blocks; the old CLI's
  # `total %/% nThreads` split gave nThreads + 1 files. Neither happens here.
  expect_equal(length(SAIGE:::splitMarkerShards(3L, 64L)), 3L)
  expect_equal(length(SAIGE:::splitMarkerShards(0L, 4L)), 0L)
  expect_identical(SAIGE:::splitMarkerShards(10L, 4L),
                   list(1:3, 4:6, 7:8, 9:10))
})

# ---- unit: mergeShardOutputs ------------------------------------------------

shard_write <- function(dir, name, lines) {
  p <- file.path(dir, name)
  con <- file(p, open = "wb")
  writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), "\n")), con)
  close(con)
  p
}

test_that("mergeShardOutputs concatenates bodies under one header", {
  d <- tempfile("shardmerge_"); dir.create(d)
  s1 <- shard_write(d, "s1.txt", c("CHR\tPOS", "1\t10", "1\t11"))
  s2 <- shard_write(d, "s2.txt", c("CHR\tPOS", "1\t12"))
  s3 <- shard_write(d, "s3.txt", "CHR\tPOS")            # header-only: legal
  out <- file.path(d, "merged.txt")
  SAIGE:::mergeShardOutputs(c(s1, s2, s3), out)
  expect_identical(readLines(out), c("CHR\tPOS", "1\t10", "1\t11", "1\t12"))
  expect_false(file.exists(paste0(out, ".merging")))
})

test_that("mergeShardOutputs merges identically in one block or many", {
  d <- tempfile("shardmerge2_"); dir.create(d)
  body <- sprintf("1\t%d\trs%d", 1:500, 1:500)
  s1 <- shard_write(d, "s1.txt", c("CHR\tPOS\tID", body[1:250]))
  s2 <- shard_write(d, "s2.txt", c("CHR\tPOS\tID", body[251:500]))
  big <- file.path(d, "big.txt"); small <- file.path(d, "small.txt")
  SAIGE:::mergeShardOutputs(c(s1, s2), big)
  SAIGE:::mergeShardOutputs(c(s1, s2), small, blockSize = 7L)  # forces re-reads
  expect_true(files_byte_identical(big, small))
})

test_that("mergeShardOutputs refuses empty shards and header mismatches", {
  d <- tempfile("shardmerge3_"); dir.create(d)
  s1 <- shard_write(d, "s1.txt", c("CHR\tPOS", "1\t10"))
  empty <- file.path(d, "empty.txt"); file.create(empty)
  other <- shard_write(d, "other.txt", c("CHR\tPOS\tEXTRA", "1\t10\tx"))
  out <- file.path(d, "merged.txt")
  expect_error(SAIGE:::mergeShardOutputs(c(s1, empty), out), "empty")
  expect_error(SAIGE:::mergeShardOutputs(c(s1, other), out), "different header")
  # a failed merge leaves no plausible-looking partial result behind
  expect_false(file.exists(out))
})

# ---- end-to-end: P1 (ER path) and P4 (conditional, quantitative).
# P7/P11 and the repeat-identity check are in
# test-step2-threads-goldens.R; scaffolding in helper-threads.R.

for (key in c("P1", "P4")) {
  local({
    k <- key
    test_that(paste0("threaded marker scan reproduces golden ", k), {
      skip_unless_level("standard")
      skip_if_no_goldens()
      golden <- file.path(saige_goldens, "step2", th_workloads[[k]]$golden)
      for (nt in TH_NTHREADS) {
        out <- th_run(k, nt)
        expect_output_matches_golden(golden, out)
      }
    })
  })
}

test_that("the merged index file marks one completed run", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- th_run("P1", 4, tag = "_idx")
  idx <- readLines(paste0(out, ".index"))
  # Not byte-identical to a serial run's index (there are no chunks to restart
  # from), and deliberately so: compare_goldens.sh excludes *.index.
  expect_true(any(grepl("^This is a Marker level analysis", idx)))
  expect_true(any(grepl("^Have completed the analyses of all chunks", idx)))
  # no shard scratch left behind
  expect_length(list.files(dirname(out), pattern = "^\\.SAIGE_marker_shards_"), 0)
  expect_false(file.exists(paste0(out, ".merging")))
})

test_that("group tests fall back to one thread with a warning", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step2_th_group_"), "group.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  status <- run_saige_cli("step2_SPAtests.R", c(
    TH_PLINK_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", "./output/example_binary_sparseGRM.rda"),
    paste0("--varianceRatioFile=",
           "./output/example_binary_sparseGRM.varianceRatio.txt"),
    paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=",
           "./output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.",
           "sparseGRM.mtx.sampleIDs.txt"),
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--is_fastTest=TRUE",
    paste0("--SAIGEOutputFile=", out),
    "--nThreads=4"), log = log)
  expect_cli_ok(status)
  expect_true(any(grepl("only supported for single-variant tests",
                        readLines(log, warn = FALSE))))
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_group_fastTest.txt"), out)
})

test_that("varianceRatioBinsHaveGap accepts contiguous bins, rejects gapped", {
  # the CLI default, and the harness's 8-category model b: both contiguous
  expect_false(SAIGE:::varianceRatioBinsHaveGap(c(10.5, 20.5), 20.5))
  expect_false(SAIGE:::varianceRatioBinsHaveGap(
    c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)))
  # bin 3 starts above everything covered so far, leaving MAC = 3 uncovered
  expect_true(SAIGE:::varianceRatioBinsHaveGap(
    c(0.5, 1.5, 3.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5)))
  # single variance ratio: no bins consulted, nothing to prove
  expect_false(SAIGE:::varianceRatioBinsHaveGap(10.5, numeric(0)))
  # cannot prove coverage -> refuse
  expect_true(SAIGE:::varianceRatioBinsHaveGap(0.5, c(1.5, 2.5)))
  expect_true(SAIGE:::varianceRatioBinsHaveGap(c(0.5, NA), c(1.5, 2.5)))
})

test_that("gapped variance-ratio bins are refused, serially and threaded", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # LEDGER #32, now FIXED (S3). With this gap at MAC = 3 the uncovered markers
  # used to reuse whichever marker was tested before them: wrong serially, and
  # measurably order-dependent under threading (the serial run disagreed with
  # the threaded one at nThreads 5,6,7,10,12,15,20,30,50). The carry is gone --
  # selectVarianceRatio throws when a MAC falls in no bin -- so BOTH runs must
  # now stop. This assertion was inverted with the fix: it previously required
  # the serial run to succeed, which is what pinned the bug.
  # The two arms stop for different reasons and that is deliberate: the threaded
  # arm is refused up front in R by varianceRatioBinsHaveGap (a clearer message,
  # before a single marker is read), the serial arm at the first marker that
  # lands in the gap.
  gapArgs <- c(
    TH_PLINK_G, "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=",
           file.path(saige_goldens, "step1", "binary_sparseGRM_b.rda")),
    paste0("--varianceRatioFile=",
           file.path(saige_goldens, "step1", "binary_sparseGRM_b.varianceRatio.txt")),
    "--cateVarRatioMinMACVecExclude=0.5,1.5,3.5,3.5,4.5,5.5,10.5,20.5",
    "--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,20.5",
    "--is_output_moreDetails=TRUE")

  out1 <- file.path(tempfile("step2_th_gap1_"), "marker.txt")
  dir.create(dirname(out1), recursive = TRUE)
  log1 <- tempfile(fileext = ".log")
  status1 <- run_saige_cli("step2_SPAtests.R",
                           c(gapArgs, paste0("--SAIGEOutputFile=", out1),
                             "--nThreads=1"), log = log1)
  expect_false(identical(as.integer(status1), 0L))
  expect_true(any(grepl("no categorical variance ratio covers MAC",
                        readLines(log1, warn = FALSE))))

  out2 <- file.path(tempfile("step2_th_gap2_"), "marker.txt")
  dir.create(dirname(out2), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  status <- run_saige_cli("step2_SPAtests.R",
                          c(gapArgs, paste0("--SAIGEOutputFile=", out2),
                            "--nThreads=5"), log = log)
  expect_false(identical(as.integer(status), 0L))
  expect_true(any(grepl("gap in MAC", readLines(log, warn = FALSE))))
})

test_that("restarting is refused at nThreads > 1", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # A serial run interrupted after chunk 1 leaves an index file whose last
  # record is not "all chunks"; re-running it threaded must stop rather than
  # silently produce a partial or duplicated result.
  out <- file.path(tempfile("step2_th_restart_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  file.create(out)
  writeLines(c(paste0("This is the output index file for SAIGE package to record the end ",
                      "point in case users want to restart the analysis. Please do not ",
                      "modify this file."),
               "This is a Marker level analysis.",
               "nEachChunk = 1000",
               "Have completed the analysis of chunk 2"),
             paste0(out, ".index"))
  log <- tempfile(fileext = ".log")
  status <- run_saige_cli("step2_SPAtests.R",
                          c(th_workloads[["P7"]]$args(),
                            paste0("--SAIGEOutputFile=", out),
                            "--is_overwrite_output=FALSE",
                            "--nThreads=4"), log = log)
  expect_false(identical(as.integer(status), 0L))
  expect_true(any(grepl("Restarting an interrupted analysis is not supported",
                        readLines(log, warn = FALSE))))
})
