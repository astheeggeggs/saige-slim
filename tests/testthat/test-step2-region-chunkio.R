# LEDGER #57: the chunked VarMat assembly's temporary .bin files.
#
# A region test with more than one chunk writes P1Mat/P2Mat to
# `<SAIGEOutputFile>_P1Mat_Chunk_<i>.bin` and reads them back. Both
# `arma::Mat::save()` and `.load()` return a bool and every return value used to
# be discarded, and Armadillo does not cover for that either: save/load failures
# warn at `arma_warn(3, ...)` while `ARMA_WARN_LEVEL` defaults to 2. A failed
# load left the matrix empty, the assembly skipped that block, SKAT threw for
# every cell of the region, and `get_SKAT_pvalue`'s `try(silent = TRUE)` turned
# the whole gene into NA rows -- exit 0, no warning.
#
# A thing worth knowing before touching this: the bundled fixture DOES take the
# chunked path, which is not obvious from the group file. Each region has only
# three or four non-ultra-rare markers, far under `markers_per_chunk_in_groupTest`,
# so it looks like a single chunk -- but the ultra-rare collapse always appends
# one more, so `nchunks == 2` and the chunked assembly runs on every non-BURDEN
# group golden. `--markers_per_chunk_in_groupTest` therefore changes nothing
# here; do not reach for it to "force" chunking.
#
# Since 2026-09-06 the chunks stay IN MEMORY (RegionChunkStore, src/Main.cpp)
# and the .bin round trip runs only when their bytes exceed
# SAIGE_REGION_CHUNK_MEMORY_MB (default 8192). At N = 50,000 the round trip plus
# the resize/truncation churn around it was 42% of the region C++ time for a
# 1.6% matrix product (harness/bench/NOTES.md). The tests below force the spill
# with a budget of 0 so the disk path and its LEDGER #57 checks stay gated, and
# prove the default path never touches the file system by planting a directory
# where the chunk file would go and expecting the run to SUCCEED.

CHUNK_G <- c("--bedFile=./input/genotype_100markers.bed",
             "--bimFile=./input/genotype_100markers.bim",
             "--famFile=./input/genotype_100markers.fam",
             "--AlleleOrder=alt-first")
CHUNK_MODEL <- c(
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
)
CHUNK_SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)

# Aliases: helper-saige.R is a separate file, so lintr cannot see these.
cli <- run_saige_cli
P6 <- file.path(saige_goldens, "step2", "plink_group_nofast.txt")

# Golden P6's command (group test, SKAT-O, is_fastTest=FALSE).
chunk_run <- function(out, env = character()) {
  cli("step2_SPAtests.R", c(
    CHUNK_G, CHUNK_MODEL, CHUNK_SG,
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=FALSE",
    paste0("--SAIGEOutputFile=", out)
  ), env = env)
}
# A budget of 0 MB forces every chunk to disk: the pre-2026-09-06 path.
SPILL <- c(SAIGE_REGION_CHUNK_MEMORY_MB = "0")

test_that("the bundled group fixture really does take the chunked VarMat path", {
  # Positive control for the failure test below: injecting a write failure at
  # chunk 0 proves nothing unless chunk 0 is written in the first place.
  skip_unless_level("standard")
  skip_if_no_goldens()

  out <- file.path(tempfile("step2_chunkio_ok_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- chunk_run(out)
  expect_cli_ok(st)

  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_gt(sum(grepl("Analyzing chunks (", lg, fixed = TRUE)), 0L)

  # And the checked I/O is byte-neutral where it succeeds, which is here.
  expect_output_matches_golden(P6, out)

  # The temporaries are cleaned up, so a failure to remove them is visible too.
  expect_length(list.files(dirname(out), pattern = "_Chunk_[0-9]+\\.bin$"), 0L)
})

test_that("at the default budget the chunks never touch the file system", {
  # Plant a DIRECTORY where chunk 0 would be written. On the disk path this is
  # the failure the next test pins; on the in-memory path nothing is written,
  # so the run must succeed and reproduce the golden byte for byte.
  skip_unless_level("standard")
  skip_if_no_goldens()

  out <- file.path(tempfile("step2_chunkio_mem_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  dir.create(paste0(out, "_P1Mat_Chunk_0.bin"))

  st <- chunk_run(out)
  expect_cli_ok(st)
  expect_output_matches_golden(P6, out)
})

test_that("the forced spill to disk is byte-neutral and cleans up", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  out <- file.path(tempfile("step2_chunkio_spill_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- chunk_run(out, env = SPILL)
  expect_cli_ok(st)
  expect_output_matches_golden(P6, out)
  expect_length(list.files(dirname(out), pattern = "_Chunk_[0-9]+\\.bin$"), 0L)
})

test_that("a chunk file that cannot be written stops the run instead of NA-ing a gene", {
  skip_unless_level("standard")

  out <- file.path(tempfile("step2_chunkio_fail_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  # A DIRECTORY where the first chunk file must go: portable, deterministic, and
  # it leaves the output file itself perfectly writable, so the run fails at
  # exactly the step under test and nowhere earlier.
  dir.create(paste0(out, "_P1Mat_Chunk_0.bin"))

  st <- chunk_run(out, env = SPILL)
  expect_false(identical(as.integer(st), 0L))

  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("Failed to WRITE the temporary chunk file", lg, fixed = TRUE)))
  # The message has to name the region, or an operator cannot act on it.
  expect_true(any(grepl("region 'GENE1'", lg, fixed = TRUE)))

  # Before the fix this run exited 0 and wrote a full-looking table whose GENE1
  # rows were silently all-NA. Nothing partial should survive now: since
  # LEDGER #116 the header is written when the run starts, so at most that
  # one line is on disk, and no result row.
  if (file.exists(out)) {
    lines <- readLines(out, warn = FALSE)
    expect_lte(length(lines), 1L)
    if (length(lines) == 1L) expect_true(startsWith(lines[1], "Region\t"))
  }
})
