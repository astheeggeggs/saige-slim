# The multi-model marker scan over the 128,868-marker fixture: chunking,
# threading and RESTART, which are the three things a per-model output file
# could plausibly break and a 70-marker single-chunk run cannot see.
#
# Its own file because it is slow, and the tier's wall time is bounded below by
# its slowest single file.

MMC_BED <- "./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr"

mmc_manifest <- function() {
  path <- tempfile(fileext = ".txt")
  bin <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  qnt <- file.path(saige_goldens, "step1", "quant_sparseGRM")
  writeLines(c(
    paste("binary", paste0(bin, ".rda"), paste0(bin, ".varianceRatio.txt")),
    paste("quant", paste0(qnt, ".rda"), paste0(qnt, ".varianceRatio.txt"))
  ), path)
  path
}

# P7's flags, minus --SAIGEOutputFile and the single-model file arguments.
mmc_args <- function(man) {
  c(
    paste0("--bedFile=", MMC_BED, ".bed"),
    paste0("--bimFile=", MMC_BED, ".bim"),
    paste0("--famFile=", MMC_BED, ".fam"),
    "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    "--markers_per_chunk=1000",
    paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)
  )
}

mmc_run <- function(out, extra = character()) {
  st <- run_saige_cli("step2_SPAtests.R",
                      c(mmc_args(mmc_manifest()), extra,
                        paste0("--SAIGEOutputFile=", out)))
  expect_cli_ok(st)
  st
}

test_that("a chunked multi-model scan reproduces the chunked golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step2_mmc_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  mmc_run(out)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_chunked.txt"),
    paste0(out, "_binary")
  )
  # Both files got every chunk's index line, each with its own byte count.
  for (sfx in c("_binary", "_quant")) {
    idx <- readLines(paste0(out, sfx, ".index"), warn = FALSE)
    chunkLines <- grep("^Have completed the analysis of chunk ", idx, value = TRUE)
    expect_gt(length(chunkLines), 4)
    expect_true(all(grepl(" bytes [0-9]+$", chunkLines)))
  }
})

test_that("threading a chunked multi-model scan is byte-identical to serial", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  serial <- file.path(tempfile("step2_mmc_s_"), "marker.txt")
  threaded <- file.path(tempfile("step2_mmc_t_"), "marker.txt")
  dir.create(dirname(serial), recursive = TRUE)
  dir.create(dirname(threaded), recursive = TRUE)
  mmc_run(serial)
  st <- mmc_run(threaded, "--nThreads=4")
  expect_true(any(grepl("Marker scan will use 4 in-process threads",
                        readLines(attr(st, "log"), warn = FALSE))))
  # ONE code path re-run over a disjoint partition of the markers, so
  # byte-identity is the right gate (the same argument --nThreads itself rests
  # on), not a tolerance.
  for (sfx in c("_binary", "_quant")) {
    expect_identical(readLines(paste0(threaded, sfx)),
                     readLines(paste0(serial, sfx)))
  }
})

test_that("a restart with BOTH files mid-chunk reproduces the complete run", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # Reconstructed rather than killed, for the reason test-step2-threads.R gives:
  # where the stream buffer sits when a signal lands is timing-dependent, and a
  # flaky duplication test is worse than none. Truncating past a recorded chunk
  # boundary is the same state, exactly, every time -- and here it must be done
  # to BOTH files, each by its OWN recorded byte count, which is the thing a
  # single shared count would get wrong.
  complete <- file.path(tempfile("step2_mmcr_full_"), "marker.txt")
  dir.create(dirname(complete), recursive = TRUE)
  mmc_run(complete)

  restarted <- file.path(tempfile("step2_mmcr_part_"), "marker.txt")
  dir.create(dirname(restarted), recursive = TRUE)
  keep <- 3L
  for (sfx in c("_binary", "_quant")) {
    idx <- readLines(paste0(complete, sfx, ".index"), warn = FALSE)
    chunkLines <- grep("^Have completed the analysis of chunk ", idx, value = TRUE)
    bytesAt <- as.numeric(sub("^.* bytes ([0-9]+)$", "\\1", chunkLines))
    partialBytes <- floor((bytesAt[keep] + bytesAt[keep + 1L]) / 2)
    expect_gt(partialBytes, bytesAt[keep])

    file.copy(paste0(complete, sfx), paste0(restarted, sfx))
    con <- file(paste0(restarted, sfx), open = "r+b")
    seek(con, where = partialBytes, origin = "start")
    truncate(con)
    close(con)
    writeLines(c(idx[seq_len(3L)], chunkLines[seq_len(keep)]),
               paste0(restarted, sfx, ".index"))
  }

  st <- mmc_run(restarted, "--is_overwrite_output=FALSE")
  for (sfx in c("_binary", "_quant")) {
    expect_identical(readLines(paste0(restarted, sfx)),
                     readLines(paste0(complete, sfx)))
  }
  # LEDGER #49 FIXED: the restart resumes AFTER the last completed chunk, as
  # checkOutputFile's message always claimed ("we restart the analysis from
  # the N + 1 chunk"). It used to resume AT chunk N and recompute it -- the
  # output was right (the truncation target was chunk N - 1 to match), but a
  # chunk of work was repeated on every restart.
  lg <- readLines(attr(st, "log"), warn = FALSE)
  analysed <- as.integer(sub("^.*Analyzing Chunk ([0-9]+)/.*$", "\\1",
                             grep("Analyzing Chunk [0-9]+/", lg, value = TRUE)))
  nChunks <- length(grep("^Have completed the analysis of chunk ",
                         readLines(paste0(complete, "_binary.index"), warn = FALSE)))
  expect_identical(analysed, seq.int(keep + 1L, nChunks))
  expect_true(any(grepl(paste0("restart the analysis from the ", keep + 1L, " chunk"), lg)))
})

test_that("index files that disagree about the restart chunk stop the run", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # v1 refuses rather than resuming at the minimum: a disagreement means the
  # files were not all written by the same interrupted run, and guessing which
  # is authoritative is how rows get silently duplicated.
  complete <- file.path(tempfile("step2_mmcd_"), "marker.txt")
  dir.create(dirname(complete), recursive = TRUE)
  mmc_run(complete)
  idxB <- readLines(paste0(complete, "_binary.index"), warn = FALSE)
  idxQ <- readLines(paste0(complete, "_quant.index"), warn = FALSE)
  cutB <- grep("^Have completed the analysis of chunk ", idxB, value = TRUE)
  cutQ <- grep("^Have completed the analysis of chunk ", idxQ, value = TRUE)
  writeLines(c(idxB[seq_len(3L)], cutB[seq_len(3L)]), paste0(complete, "_binary.index"))
  writeLines(c(idxQ[seq_len(3L)], cutQ[seq_len(5L)]), paste0(complete, "_quant.index"))

  st <- run_saige_cli("step2_SPAtests.R",
                      c(mmc_args(mmc_manifest()),
                        paste0("--SAIGEOutputFile=", complete),
                        "--is_overwrite_output=FALSE"))
  expect_false(identical(as.integer(st), 0L))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("disagree about where to restart", lg)))
  expect_true(any(grepl("_binary.index", lg)) && any(grepl("_quant.index", lg)))
})
