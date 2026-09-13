# LEDGER #140: a region run interrupted MID-chunk restarts to the complete
# run's bytes, in both files, and says so when it cannot.
#
# The marker path's #47 fix (S9) records the output size on every index line
# and truncates to it on restart. The region path never got it: its index was
# appended by R at once while the BURDEN rows sat in an unflushed C++ stream,
# so the two could be out of step in either direction -- rows on disk past the
# indexed chunk (duplicated on restart, #47's shape) or an index ahead of the
# file (regions silently MISSING, with no duplicate row to reveal it). The
# existing restart test (test-step2-region-restart.R) reconstructs the
# interruption exactly AT a chunk boundary, the benign case, which is why the
# suite never saw this.
#
# Reconstructed rather than killed, as that test is: the complete run's files
# are cut back to region 1 and then given a TORN tail (one whole GENE2 row and
# half of another, in both the results file and .singleAssoc.txt), the index
# says chunk 1 is done at the sizes BEFORE the tail, and the same command is
# re-run with --is_overwrite_output=FALSE. Before the fix the tail was kept and
# GENE2 appended after it; the restart could not reproduce the complete run.
# Tier: standard (one complete run per path plus the restarts).

rmc_model <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
rmc_grm <- file.path(saige_extdata, "output",
                     "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")

# top-level aliases, as test-step2-region-flip.R does, so the helper functions
# below can be linted (the helpers live in helper-saige.R)
rmc_cli <- run_saige_cli
rmc_ok <- expect_cli_ok

rmc_run <- function(out, extra = character()) {
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0", "--minMAC=0.5",
            "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", rmc_model, ".rda"),
            paste0("--varianceRatioFile=", rmc_model, ".varianceRatio.txt"),
            paste0("--sparseGRMFile=", rmc_grm),
            paste0("--sparseGRMSampleIDFile=", rmc_grm, ".sampleIDs.txt"),
            "--groupFile=./input/group_new_snpid.txt",
            "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
            "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE", extra)
  st <- rmc_cli("step2_SPAtests.R", args)
  rmc_ok(st)
  readLines(attr(st, "log"), warn = FALSE)
}

# Which rows of a results file belong to region 1. The region file keys on the
# region name; .singleAssoc.txt keys on the marker, so region 1's rows there are
# its group-file markers (bim column 2 rsIDs) plus its "GENE1:..." pseudo-marker
# rows. Regions are written one after another, so the split is a prefix.
rmc_gene1_ids <- local({
  gl <- readLines(file.path(saige_extdata, "input", "group_new_snpid.txt"))
  v <- grep("^GENE1[ \t]+var[ \t]", gl, value = TRUE)[1]
  strsplit(v, "[ \t]+")[[1]][-(1:2)]
})
rmc_is_gene1 <- function(lines, single) {
  if (!single) return(grepl("^GENE1\t", lines))
  id <- vapply(strsplit(lines, "\t", fixed = TRUE), `[`, character(1), 3L)
  id %in% rmc_gene1_ids | startsWith(id, "GENE1:")
}

# Write `path` as the header plus region 1's rows of `full`, return its size,
# then append a torn tail: region 2's first row whole and its second cut in half.
rmc_write_torn <- function(path, full, single = FALSE) {
  g1 <- which(rmc_is_gene1(full, single))
  testthat::expect_true(length(g1) > 0 && all(diff(g1) == 1) && g1[1] == 2)
  g2 <- seq(max(g1) + 1, length(full))
  testthat::expect_true(length(g2) >= 2)
  writeLines(full[c(1, g1)], path)
  size <- file.size(path)
  half <- substr(full[g2[2]], 1, nchar(full[g2[2]]) %/% 2)
  cat(full[g2[1]], "\n", half, sep = "", file = path, append = TRUE)
  testthat::expect_gt(file.size(path), size)
  size
}

rmc_restart_case <- function(label, extra, hasSingle) {
  complete <- file.path(tempfile(paste0("step2_rmc_full_", label, "_")), "region.txt")
  dir.create(dirname(complete), recursive = TRUE)
  rmc_run(complete, extra)
  full <- readLines(complete)
  idx <- readLines(paste0(complete, ".index"), warn = FALSE)
  testthat::expect_identical(tail(idx, 1), "Have completed the analyses of all chunks.")
  # the complete run's index carries the sizes (the #140 record itself)
  chunkLines <- grep("^Have completed the analysis of chunk ", idx, value = TRUE)
  testthat::expect_true(all(grepl(" bytes [0-9]+$", chunkLines)))
  if (hasSingle) {
    testthat::expect_true(all(grepl(" single_bytes [0-9]+ bytes", chunkLines)))
  } else {
    testthat::expect_false(any(grepl("single_bytes", chunkLines)))
  }
  singleFull <- if (hasSingle) readLines(paste0(complete, ".singleAssoc.txt")) else NULL
  testthat::expect_true(sum(grepl("^GENE2\t", full)) >= 2)

  restarted <- file.path(tempfile(paste0("step2_rmc_part_", label, "_")), "region.txt")
  dir.create(dirname(restarted), recursive = TRUE)
  bytes <- rmc_write_torn(restarted, full)
  singleBytes <- if (hasSingle) {
    rmc_write_torn(paste0(restarted, ".singleAssoc.txt"), singleFull, single = TRUE)
  }
  line <- paste("Have completed the analysis of chunk 1",
                if (hasSingle) paste("single_bytes", singleBytes) else NULL,
                "bytes", bytes)
  writeLines(c(idx[1:3], line), paste0(restarted, ".index"))
  lg <- rmc_run(restarted, c(extra, "--is_overwrite_output=FALSE"))
  testthat::expect_true(any(grepl("we restart the analysis from the 2 chunk", lg, fixed = TRUE)))
  # LEDGER #140 FIXED: the torn tail is discarded (and the run says so, per
  # file) and the restart reproduces the complete run's bytes in BOTH files.
  # Before the fix the tail stayed and GENE2 was appended after it.
  testthat::expect_true(any(grepl("Discarded", lg, fixed = TRUE)))
  testthat::expect_identical(readLines(restarted), full)
  if (hasSingle) {
    testthat::expect_identical(readLines(paste0(restarted, ".singleAssoc.txt")), singleFull)
  }
  list(full = full, idx = idx, bytes = bytes)
}

test_that("LEDGER #140: a SKAT-O run torn mid-chunk restarts to the complete run, both files", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  rmc_restart_case("skato", character(), hasSingle = TRUE)
})

test_that("LEDGER #140: a BURDEN run (C++ writer) torn mid-chunk restarts to the complete run", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- rmc_restart_case("burden", "--r.corr=1", hasSingle = FALSE)

  # The other direction: the index AHEAD of the file. Cut the results file
  # SHORTER than the size the index records for chunk 1 -- what a kill between
  # R's index append and the C++ stream's flush used to leave -- and the
  # restart must say the file cannot be trusted rather than silently resume
  # with region 1's rows missing.
  short <- file.path(tempfile("step2_rmc_short_"), "region.txt")
  dir.create(dirname(short), recursive = TRUE)
  g1 <- grep("^GENE1\t", r$full)
  writeLines(r$full[c(1, g1[seq_len(length(g1) - 1)])], short)
  expect_lt(file.size(short), r$bytes)
  writeLines(c(r$idx[1:3], paste("Have completed the analysis of chunk 1 bytes", r$bytes)),
             paste0(short, ".index"))
  lg <- rmc_run(short, c("--r.corr=1", "--is_overwrite_output=FALSE"))
  expect_true(any(grepl("shorter than the index claims", lg, fixed = TRUE)))
})
