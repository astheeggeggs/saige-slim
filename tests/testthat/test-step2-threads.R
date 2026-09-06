# --nThreads for the single-variant marker scan.
#
# S7 retired the forked mechanism, so --nThreads now drives the IN-PROCESS
# marker threads and this file tests that one path through the real CLI flag.
# The partition/merge unit tests that used to live here went with
# splitMarkerShards/appendShardBytes/mergeShardOutputs.
#
# Two assertions here were INVERTED by S7 and the reasons are in each block:
# restarting is now supported (the fork path refused it), and a reader that
# cannot be cloned degrades with a warning rather than erroring.
#
# P7/P11 and the repeat-identity check are in test-step2-threads-goldens.R;
# scaffolding in helper-threads.R.

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

test_that("a threaded run writes the same index a serial run does", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # S7: with the forked path gone there are no shards and no merge, so the
  # index is the ordinary per-chunk one -- not the single "all chunks" record
  # the fork path used to leave. Compared to the SERIAL index rather than to a
  # fixed string, because the point is that threading changed nothing here.
  out4 <- th_run("P7", 4, tag = "_idx4")
  out1 <- th_run("P7", 1, tag = "_idx1")
  expect_identical(readLines(paste0(out4, ".index")),
                   readLines(paste0(out1, ".index")))
  # and no scratch left behind from the retired mechanism
  expect_length(list.files(dirname(out4), pattern = "^\\.SAIGE_marker_shards_"), 0)
  expect_false(file.exists(paste0(out4, ".merging")))
})

test_that("group tests thread the marker loop and say so", {
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
  lg <- readLines(log, warn = FALSE)
  # 2026-09-04: --nThreads stopped being DISCARDED on a group run (it drove the
  # SKAT-O post-processing, R/SkatoParallel.R). 2026-09-06: it threads the C++
  # region marker loop as well, and the count is PRINTED so a knob that never
  # arrived cannot pass as green (LEDGER #30). The old "does not thread the
  # group test's marker loop" line is deliberately gone.
  expect_true(any(grepl("Region marker loop will use 4 in-process threads", lg,
                        fixed = TRUE)))
  expect_false(any(grepl("does not thread the group test's marker loop", lg,
                         fixed = TRUE)))
  expect_false(any(grepl("only supported for single-variant tests", lg,
                         fixed = TRUE)))
  # and the threads went somewhere rather than nowhere
  expect_true(any(grepl("SKAT-O post-processing will use 4 workers", lg,
                        fixed = TRUE)))
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
  # lands in the gap. Post-S7 the R guard's message no longer claims the run is
  # unreproducible -- S3 fixed that -- only that the bins leave MAC uncovered.
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

test_that("restarting an interrupted run works, and matches a serial restart", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # INVERTED BY S7. The forked path refused restart outright -- its shards were
  # a partition of the marker list, not of the chunk sequence, so a partially
  # written serial output meant nothing to them. The in-process path is one
  # process running the ordinary chunk loop, so restart is simply the serial
  # behaviour, and that is what is asserted: threaded and serial restarts from
  # the SAME interrupted state must agree byte for byte.
  #
  # This deliberately does NOT assert the restarted file equals a complete run,
  # and the reason CHANGED with the LEDGER #47 fix. It used to be that a restart
  # duplicated rows; now it does not, but the state seeded below is a LEGACY
  # index -- no byte count, which is what an older SAIGE wrote -- so the restart
  # cannot verify what is on disk and says so. Here the seeded output is empty
  # while the index claims two chunks done, so chunks 1-2 are missing from the
  # result either way. The parity claim (both arms agree) is what this test is
  # for; the #47 regression is the test below.
  seed_state <- function(out) {
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    file.create(out)
    writeLines(c(paste0("This is the output index file for SAIGE package to record the ",
                        "end point in case users want to restart the analysis. Please ",
                        "do not modify this file."),
                 "This is a Marker level analysis.",
                 "nEachChunk = 1000",
                 "Have completed the analysis of chunk 2"),
               paste0(out, ".index"))
  }
  run_restart <- function(nt, tag) {
    out <- file.path(tempfile(paste0("step2_th_restart_", tag, "_")), "marker.txt")
    seed_state(out)
    status <- run_saige_cli("step2_SPAtests.R",
                            c(th_workloads[["P7"]]$args(),
                              paste0("--SAIGEOutputFile=", out),
                              "--is_overwrite_output=FALSE",
                              paste0("--nThreads=", nt)))
    expect_cli_ok(status)
    out
  }
  serial <- run_restart(1, "serial")
  threaded <- run_restart(4, "threaded")
  expect_identical(readLines(threaded), readLines(serial))
})

test_that("a restart drops what the interrupted run wrote past its last chunk (LEDGER #47)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # THE #47 REGRESSION. The bug: the `.index` is written once per chunk from R
  # while the output stream flushes on its own ~192 KiB boundary, so a killed
  # run leaves rows belonging to chunks the index never recorded, and appending
  # from the indexed chunk writes them AGAIN (measured at 3,826 duplicated
  # markers on this fixture).
  #
  # The interrupted state is RECONSTRUCTED rather than produced by a kill, on
  # purpose: where the stream buffer happens to sit when a signal lands is
  # timing-dependent, and a flaky test that only sometimes reproduces a
  # duplication bug is worse than none. Truncating a complete run's output to a
  # point PAST a recorded chunk boundary is the same state, exactly and every
  # time.
  complete <- file.path(tempfile("step2_restart_full_"), "marker.txt")
  dir.create(dirname(complete), recursive = TRUE)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(th_workloads[["P7"]]$args(),
                                paste0("--SAIGEOutputFile=", complete))))
  idx <- readLines(paste0(complete, ".index"), warn = FALSE)
  chunkLines <- grep("^Have completed the analysis of chunk ", idx, value = TRUE)
  expect_gt(length(chunkLines), 4)
  # The fix's own record: each chunk line now carries the output size at that
  # point. Without it there is nothing to truncate to, so assert it is there.
  expect_true(all(grepl(" bytes [0-9]+$", chunkLines)),
              label = "every chunk line records a byte count")
  bytesAt <- as.numeric(sub("^.* bytes ([0-9]+)$", "\\1", chunkLines))
  expect_false(is.unsorted(bytesAt))

  # Interrupted state: index stops after chunk 3, but the file on disk holds
  # chunk 3 plus HALF of chunk 4 -- rows that a naive append would duplicate.
  keep <- 3L
  partialBytes <- floor((bytesAt[keep] + bytesAt[keep + 1L]) / 2)
  expect_gt(partialBytes, bytesAt[keep])

  restarted <- file.path(tempfile("step2_restart_partial_"), "marker.txt")
  dir.create(dirname(restarted), recursive = TRUE)
  file.copy(complete, restarted)
  con <- file(restarted, open = "r+b")
  seek(con, where = partialBytes, origin = "start")
  truncate(con)
  close(con)
  writeLines(c(idx[seq_len(3L)], chunkLines[seq_len(keep)]),
             paste0(restarted, ".index"))

  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(th_workloads[["P7"]]$args(),
                                paste0("--SAIGEOutputFile=", restarted),
                                "--is_overwrite_output=FALSE")))
  # The claim: a restart reproduces the file a complete run would have written,
  # byte for byte. Before the fix this was 132,695 rows against 128,869.
  expect_identical(readLines(restarted), readLines(complete))
})

test_that("a reader that cannot be cloned falls back to one thread with a warning", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # S7 decision: DEGRADE, do not error. Erroring would break anyone scripting
  # one --nThreads across a mix of inputs, and it is the same posture group
  # tests already take. (The mix that motivated it was PLINK vs PGEN; PGEN
  # threads since S8, so the pread fallback below is the last refusal.) It is not the silent no-op of LEDGER #30
  # because the reason is printed AND raised as an R warning -- which is why
  # this asserts the message, not just that the run succeeded.
  #
  # SAIGE_PLINK_NO_MMAP forces the pread fallback, which refuses to clone (S5).
  out <- file.path(tempfile("step2_th_noclone_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  status <- run_saige_cli("step2_SPAtests.R",
                          c(th_workloads[["P1"]]$args(),
                            paste0("--SAIGEOutputFile=", out),
                            "--nThreads=4"),
                          log = log,
                          env = c(SAIGE_PLINK_NO_MMAP = "1"))
  expect_cli_ok(status)
  txt <- readLines(log, warn = FALSE)
  expect_true(any(grepl("cannot be cloned", txt)))
  expect_true(any(grepl("run on ONE thread|run on one thread", txt)))
  # and the fallback still produces the right answer
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", th_workloads[["P1"]]$golden), out)
})
