# The carrier-space PGEN read (2026-09-16). A record the pgen WRITER stored as
# a difference list is handed to the marker path as (common code, listed
# samples) and prepared over the listed samples only (src/Main.cpp
# prepareSparseMarker, src/UTIL.cpp imputeSparseAndFlip); every other record,
# and every PLINK record, takes the dense path as before. The claim under test
# is that the two paths are the SAME arithmetic: the sparse run must reproduce
# the committed pgen goldens byte for byte AND the dense run of the same
# command (SAIGE_PGEN_SPARSE=0, the kill switch) byte for byte, and the run
# must REPORT that it took the sparse path -- a knob that fails to arrive would
# otherwise pass every byte gate as green (the LEDGER #30 shape).
#
# This file adds CROSS-RUN comparisons, so it is gated in BOTH compare modes
# (CLAUDE.md, test-only changes). Tier: standard.

PGEN_G <- c("--pgenFile=./input/genotype_100markers.pgen",
            "--pvarFile=./input/genotype_100markers.pvar",
            "--psamFile=./input/genotype_100markers.psam")
PGEN_CHUNK_G <- c(
  "--pgenFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.pgen",
  "--pvarFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.pvar",
  "--psamFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.psam"
)
MISS_PGEN_G <- c("--pgenFile=./input/genotype_100markers_missing.pgen",
                 "--pvarFile=./input/genotype_100markers_missing.pvar",
                 "--psamFile=./input/genotype_100markers_missing.psam")
BIN_MODEL <- function() c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt"))
)
QUANT_MODEL <- function() c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.varianceRatio.txt"))
)
MARKER_FLAGS <- c("--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
                  "--is_output_moreDetails=TRUE")

sparse_run <- function(args, tag, sparse = TRUE) {
  out <- file.path(tempfile(paste0("step2_sparse_", tag, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  env <- if (sparse) character() else c(SAIGE_PGEN_SPARSE = "0")
  status <- run_saige_cli("step2_SPAtests.R",
                          c(args, paste0("--SAIGEOutputFile=", out)), env = env)
  list(status = status, path = out, log = attr(status, "log"))
}

# The per-chunk report lines, summed: how many markers the run read sparsely.
sparse_reads <- function(log) {
  ln <- grep("^Sparse reads: ", readLines(log, warn = FALSE), value = TRUE)
  if (!length(ln)) return(0L)
  sum(as.integer(sub("^Sparse reads: ([0-9]+) of .*$", "\\1", ln)))
}

expect_sparse_equals_dense <- function(args, tag, golden = NULL) {
  s <- sparse_run(args, paste0(tag, "_sparse"), sparse = TRUE)
  d <- sparse_run(args, paste0(tag, "_dense"), sparse = FALSE)
  expect_cli_ok(s$status)
  expect_cli_ok(d$status)
  expect_gt(sparse_reads(s$log), 0L)
  expect_identical(sparse_reads(d$log), 0L)
  expect_files_byte_identical(d$path, s$path)
  if (!is.null(golden)) {
    expect_output_matches_golden(golden, s$path)
    expect_output_matches_golden(golden, d$path)
  }
  invisible(list(sparse = s, dense = d))
}

test_that("sparse read reproduces golden P13 (pgen marker scan) and reports itself", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_sparse_equals_dense(
    c(PGEN_G, MARKER_FLAGS, BIN_MODEL()), "p13",
    golden = file.path(saige_goldens, "step2", "pgen_marker_binary.txt"))
})

test_that("sparse read reproduces golden P20 (chunked pgen scan, 22chr fixture)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- expect_sparse_equals_dense(
    c(PGEN_CHUNK_G, "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
      BIN_MODEL(), "--markers_per_chunk=1000"), "p20",
    golden = file.path(saige_goldens, "step2", "pgen_marker_chunked.txt"))
  # one report line per chunk that had a sparse read; the chunks are 1000 wide
  ln <- grep("^Sparse reads: ", readLines(r$sparse$log, warn = FALSE), value = TRUE)
  expect_gt(length(ln), 1L)
})

test_that("missing genotypes: mean imputation and the zero-dosage clean are the same sum", {
  # The two places the carrier-space arithmetic could drift from the dense
  # vector: a FRACTIONAL imputed value (--impute_method=mean writes 2*AF into
  # every missing call, so altCounts is a sum that depends on association
  # order -- imputeSparseAndFlip reproduces armadillo's parity-split
  # accumulate) and the clean that then zeroes it on MAC <= 10 markers.
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_sparse_equals_dense(
    c(MISS_PGEN_G, MARKER_FLAGS, BIN_MODEL(), "--impute_method=mean"), "miss_mean")
  expect_sparse_equals_dense(
    c(MISS_PGEN_G, MARKER_FLAGS, BIN_MODEL(), "--impute_method=minor"), "miss_minor")
  expect_sparse_equals_dense(
    c(MISS_PGEN_G, MARKER_FLAGS, BIN_MODEL()), "miss_best")
})

test_that("quantitative trait and the ER branch (empty zero index) are the same bytes", {
  # The ER exact test is the one consumer of the zero-genotype INDICES; the
  # marker path hands it an empty vector now and it builds the complement
  # itself. --max_MAC_for_ER=200 is what reaches the resampling arm on this
  # fixture (golden P33, CLAUDE.md Traps).
  skip_unless_level("standard")
  skip_if_no_goldens()
  expect_sparse_equals_dense(
    c(PGEN_G, MARKER_FLAGS, BIN_MODEL(), "--max_MAC_for_ER=200"), "er200")
  expect_sparse_equals_dense(
    c(PGEN_G, MARKER_FLAGS, QUANT_MODEL()), "quant")
})

test_that("four threads over the sparse read are byte-identical to one", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  one <- sparse_run(c(PGEN_G, MARKER_FLAGS, BIN_MODEL(), "--nThreads=1"), "t1")
  four <- sparse_run(c(PGEN_G, MARKER_FLAGS, BIN_MODEL(), "--nThreads=4"), "t4")
  expect_cli_ok(one$status)
  expect_cli_ok(four$status)
  expect_gt(sparse_reads(four$log), 0L)
  expect_files_byte_identical(one$path, four$path)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_marker_binary.txt"), four$path)
})
