# In-process threading of the step-2 marker loop (S6).
#
# The claim is that a threaded scan reproduces the serial one BYTE FOR BYTE, and
# that claim is sound only because each marker's arithmetic depends on nothing
# but its own genotype, the model (const since S3) and the config; its outputs
# land in its own slot (S4); and its reader is its own handle onto one shared
# mapping (S5). Byte-exactness IS the right standard here, unlike a comparison
# between two DIFFERENT code paths -- this is one code path re-run over a
# disjoint partition, which is the same argument the forked-shard gate makes.
#
# This file exists because nothing else in the suite can see the threaded path:
# nThreadsMarker defaults to 1 and there is no CLI flag, so every other golden
# test runs serially and would stay green if threading corrupted every row.
#
# Kept separate from test-step2-marker-paths.R so its CLI runs go in parallel
# with that file's rather than lengthening the tier's slowest file.
#
# Detection power is drill-verified, not assumed: with one MarkerScratch shared
# across workers instead of one per worker, this comparison fails loudly (16
# failures across the marker goldens at 4 threads), while forcing 1 thread
# THROUGH TBB rather than down the serial shortcut passes -- so the gate
# separates "TBB is involved" from "markers ran concurrently".

TH_PLINK <- c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first"
)
TH_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
               "--is_output_moreDetails=TRUE")
th_model <- function(prefix) c(
  paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", paste0(prefix, ".rda"))),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", paste0(prefix, ".varianceRatio.txt")))
)

# S7: the knob is --nThreads, the real CLI flag. This used to go through
# SAIGE_MARKER_THREADS because S6 deliberately changed no CLI contract while two
# mechanisms coexisted; that env var is gone with the forked path. `env` is now
# only for variables that have no flag at all, such as SAIGE_PLINK_NO_MMAP, and
# it is passed to the CHILD rather than set in this process -- these test files
# run in parallel workers, so mutating the parent's environment would leak.
th_run <- function(golden_base, extra_args, threads, env = character(),
                   expect_threaded = TRUE) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step2_th_"), golden_base)
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli(
    "step2_SPAtests.R",
    c(TH_PLINK, TH_COMMON, extra_args, paste0("--SAIGEOutputFile=", out),
      paste0("--nThreads=", threads)),
    env = env
  )
  expect_cli_ok(status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", golden_base), out
  )
  # PROVE THE RUN WAS ACTUALLY THREADED. Without this the whole file is a
  # potential silent no-op: if --nThreads failed to reach the marker loop,
  # every comparison above would pass by running serially and gate nothing --
  # which is exactly what --nThreads did before LEDGER #30 caught it. The
  # expected thread count is asserted, not just "some threading happened", so a
  # fallback to one thread also fails.
  if (expect_threaded) {
    log_txt <- readLines(attr(status, "log"), warn = FALSE)
    expect_true(
      any(grepl(paste0("use ", threads, " in-process threads"), log_txt, fixed = TRUE)),
      label = paste0("the run reported ", threads, " in-process threads")
    )
  }
  invisible(status)
}

test_that("the default binary marker scan is unchanged at 2 and 4 threads (P1)", {
  # P1 takes the ER branch on 2 ultra-rare markers, so this also covers the
  # main-thread deferral: those two are handed back and re-run serially while
  # the other 68 are tested on workers.
  th_run("plink_marker_binary.txt", th_model("binary_sparseGRM_a"), threads = 2)
  th_run("plink_marker_binary.txt", th_model("binary_sparseGRM_a"), threads = 4)
})

test_that("the quantitative marker scan is unchanged at 4 threads (P32)", {
  # The fully-threaded case: quantitative traits never enter the ER branch
  # (getMarkerPval gates it on traitType == "binary"), so no marker is deferred
  # and every one is tested on a worker.
  th_run("plink_marker_quantitative.txt", th_model("quant_sparseGRM"), threads = 4)
})

test_that("the ER resampling branch is unchanged at 4 threads (P33)", {
  # --max_MAC_for_ER=200 admits rs38 (~98 carriers), the one marker on this
  # fixture that reaches the resampling arm and so R's RNG. Threaded, it is
  # deferred to the main thread and re-seeded there exactly as it would be
  # serially, which is what makes the deferral byte-neutral rather than merely
  # safe.
  th_run("plink_marker_er.txt",
         c(th_model("binary_sparseGRM_a"), "--max_MAC_for_ER=200"), threads = 4)
})

test_that("a reader that cannot clone degrades to one thread instead of failing", {
  # The pread fallback refuses to clone (S5). It is the only refusal left --
  # PGEN cloned from S8 -- and it must be a notice and a serial run, not an
  # error, and NOT the silent no-op of LEDGER #30, so the reason is printed.
  status <- th_run("plink_marker_binary.txt", th_model("binary_sparseGRM_a"),
                   threads = 4, env = c(SAIGE_PLINK_NO_MMAP = "1"),
                   expect_threaded = FALSE)
  log_txt <- readLines(attr(status, "log"), warn = FALSE)
  expect_true(any(grepl("cannot be cloned", log_txt)),
              label = "the degrade-to-one-thread notice is printed")
  expect_true(any(grepl("one thread", log_txt)),
              label = "the notice says what it did")
  # ... and it really did fall back: no threaded-scan line was printed
  expect_false(any(grepl("in-process threads", log_txt)),
               label = "the degraded run did not also claim to be threaded")
})

test_that("runMarkerScan validates nThreads rather than trusting it", {
  # The old test here asserted that the two parallelism mechanisms could not be
  # combined. S7 deleted the forked one, so there is nothing left to combine and
  # that assertion has been replaced rather than kept passing vacuously.
  expect_error(
    SAIGE:::runMarkerScan(
      traitType = "binary", genoType = "plink", idstoIncludeFile = "",
      rangestoIncludeFile = "", genoIndex = c("0", "1"), CHROM = c("1", "1"),
      OutputFile = tempfile(), nMarkersEachChunk = 10000,
      isMoreOutput = FALSE, isImputation = FALSE, isFirth = FALSE,
      LOCO = FALSE, chrom = "1", isCondition = FALSE,
      isOverWriteOutput = TRUE, isAnyInclude = FALSE,
      nThreads = 0
    ),
    "positive integer"
  )
})
