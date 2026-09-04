# Two malformed-file cases that used to exit 0 with a wrong or unreadable
# result: LEDGER #58 (a .bed shorter than its .bim and .fam imply) and #59 (a
# region results file written with data rows and no header).
#
# Neither is reachable by any golden, and for the same reason in both cases:
# every bundled fixture is internally consistent, and every bundled group run
# has a first chunk that emits rows. That is the same blind spot that let #52's
# FID/IID keying pass 13 PGEN goldens. So the inputs here are DERIVED at run
# time from the committed fixtures, which also keeps the property under test
# visible in the source rather than buried in a binary.

cli3 <- run_saige_cli
BED <- file.path(saige_extdata, "input", "genotype_100markers")
P1g <- file.path(saige_goldens, "step2", "plink_marker_binary.txt")
BINMOD <- c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt"))
)

# Golden P1's command, against a caller-supplied bed/bim/fam prefix.
marker_run <- function(prefix, tag) {
  out <- file.path(tempfile(paste0("step2_int_", tag, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- cli3("step2_SPAtests.R", c(
    paste0("--bedFile=", prefix, ".bed"), paste0("--bimFile=", prefix, ".bim"),
    paste0("--famFile=", prefix, ".fam"), "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    BINMOD, "--is_output_moreDetails=TRUE",
    paste0("--SAIGEOutputFile=", out)
  ))
  list(status = st, path = out, log = readLines(attr(st, "log"), warn = FALSE))
}

# Copy the trio to a temp prefix, optionally appending phantom rows to the .bim
# so it claims more markers than the .bed holds.
stage_bed <- function(nPhantom) {
  pfx <- tempfile("bedint_")
  file.copy(paste0(BED, ".bed"), paste0(pfx, ".bed"))
  file.copy(paste0(BED, ".fam"), paste0(pfx, ".fam"))
  bim <- readLines(paste0(BED, ".bim"))
  if (nPhantom > 0) {
    last <- strsplit(bim[length(bim)], "\t")[[1]]
    extra <- vapply(seq_len(nPhantom), function(k) {
      row <- last
      row[2] <- paste0("rsPHANTOM", k)
      row[4] <- as.character(as.integer(last[4]) + 1000L * k)
      paste(row, collapse = "\t")
    }, character(1))
    bim <- c(bim, extra)
  }
  writeLines(bim, paste0(pfx, ".bim"))
  pfx
}

test_that("an intact bed/bim/fam is unaffected by the size check", {
  # Control for the two below: the check must not cost a valid run anything.
  skip_unless_level("standard")
  skip_if_no_goldens()

  r <- marker_run(stage_bed(0), "intact")
  expect_cli_ok(r$status)
  expect_output_matches_golden(P1g, r$path)
})

test_that("a .bed shorter than its .bim implies is refused, not silently short", {
  skip_unless_level("standard")

  # ONE phantom marker. Before the fix this exited 0 and wrote output
  # BYTE-IDENTICAL to the intact run -- the marker simply did not exist -- and
  # the pread fallback was the only path that noticed. The read lands in the
  # last page's mmap slack, so it does not even fault.
  r <- marker_run(stage_bed(1), "phantom1")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(any(grepl("is too short for the .bim and .fam", r$log, fixed = TRUE)))
  # The message has to carry the numbers, or an operator cannot tell which file
  # is wrong.
  expect_true(any(grepl("101 markers x 10000 samples", r$log, fixed = TRUE)))

  # THREE phantom markers pushed the read past the mapped extent entirely and
  # killed R with a SIGSEGV (exit 139). Same clean error now.
  r3 <- marker_run(stage_bed(3), "phantom3")
  expect_false(identical(as.integer(r3$status), 0L))
  expect_true(any(grepl("is too short for the .bim and .fam", r3$log, fixed = TRUE)))
})

test_that("a region file gets its header even when the first chunk emits no rows", {
  skip_unless_level("standard")

  # rs1, rs2, rs19 and rs100 are monomorphic in this fixture, so a region built
  # only from them is FOUND in the geno file (the chunk is not empty, and the
  # chunk counter advances) but produces no result rows. With
  # --groups_per_chunk=1 it is its own chunk, so before the fix it consumed the
  # one chunk allowed to write column names and GENE_REAL's row landed as the
  # first bytes of the file.
  grp <- tempfile(fileext = ".txt")
  real <- grep("^GENE1 ", readLines(file.path(saige_extdata, "input",
                                              "group_new_snpid.txt")), value = TRUE)
  writeLines(c(
    "GENE_BARREN var rs1 rs2 rs19 rs100",
    "GENE_BARREN anno lof lof lof lof",
    sub("^GENE1", "GENE_REAL", real)
  ), grp)

  out <- file.path(tempfile("step2_hdr_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- cli3("step2_SPAtests.R", c(
    paste0("--bedFile=", BED, ".bed"), paste0("--bimFile=", BED, ".bim"),
    paste0("--famFile=", BED, ".fam"), "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
           "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    paste0("--groupFile=", grp), "--annotation_in_groupTest=lof",
    "--maxMAF_in_groupTest=0.01", "--is_fastTest=FALSE", "--groups_per_chunk=1",
    paste0("--SAIGEOutputFile=", out)
  ))
  expect_cli_ok(st)

  lines <- readLines(out, warn = FALSE)
  expect_true(startsWith(lines[1], "Region\t"))
  # The barren region really was barren -- if these markers ever stop being
  # monomorphic this test would otherwise pass for the wrong reason.
  expect_identical(length(lines), 2L)
  expect_true(startsWith(lines[2], "GENE_REAL\t"))
})

# --- LEDGER #60: the same short .bed, through the step-1 GRM ingest ----------
#
# #58 fixed step 2's reader. Step 1 has its own .bed reader (BedByteSource in
# src/SAIGE_fitGLMM_fast.cpp) and it discarded the seekg/read status entirely,
# so the SAME trio of files that step 2 now refuses was accepted here and the
# GRM was quietly built from fewer markers. Driven through setgeno directly
# rather than a step-1 fit: it is the ingest that is under test, it costs
# seconds instead of minutes, and the marker counts are readable.

# Copy the trio to a temp prefix, truncating the .bed to `nKeep` whole markers
# (plus `extraBytes` of the next one) while leaving the .bim at 100 rows.
stage_short_bed <- function(nKeep, extraBytes = 0) {
  pfx <- tempfile("bedshort_")
  file.copy(paste0(BED, ".bim"), paste0(pfx, ".bim"))
  file.copy(paste0(BED, ".fam"), paste0(pfx, ".fam"))
  nSample <- length(readLines(paste0(BED, ".fam"), warn = FALSE))
  nbyte <- ceiling(nSample / 4)
  con <- file(paste0(BED, ".bed"), "rb")
  on.exit(close(con), add = TRUE)
  keep <- readBin(con, "raw", 3L + nKeep * nbyte + extraBytes)
  writeBin(keep, paste0(pfx, ".bed"))
  pfx
}

# setgeno populates a C++ global, so it is driven in a child process. Returns
# the exit status, the log, and the QC-passing marker count when it got one.
setgeno_run <- function(prefix) {
  out <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('base <- "%s"; out <- "%s"', prefix, out),
    "setminMAFforGRM(0.01); setmaxMissingRateforGRM(0.15)",
    'ids <- SAIGE:::readGenoSampleIDs("plink", famFile = paste0(base, ".fam"))',
    "n <- length(ids)",
    'setgeno(paste0(base, ".bed"), paste0(base, ".bim"), paste0(base, ".fam"),',
    "        seq_len(n), rep(TRUE, n), 2, FALSE)",
    "saveRDS(sum(getQCdMarkerIndex()), out)"
  ), script)
  log <- tempfile(fileext = ".log")
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = log, stderr = log))
  list(status = as.integer(st), log = readLines(log, warn = FALSE),
       nQC = if (file.exists(out)) readRDS(out) else NA_integer_)
}

test_that("the step-1 ingest refuses a .bed shorter than its .bim implies", {
  skip_unless_level("standard")

  # NEGATIVE CONTROL: the intact fixture, so the check costs a valid run
  # nothing and the marker count below has a reference.
  ok <- setgeno_run(stage_short_bed(100))
  expect_identical(ok$status, 0L)
  expect_gt(ok$nQC, 0L)

  # 60 of 100 markers present. Before the fix: exit 0, no warning, and a GRM
  # built from 18 QC-passing markers instead of the intact file's 33 -- the
  # zero-filled buffer decodes as genotype 2 for every sample, so altFreq is
  # 1.0 and the minMAF gate drops it.
  short <- setgeno_run(stage_short_bed(60))
  expect_false(identical(short$status, 0L))
  expect_true(any(grepl("is too short for the .bim", short$log, fixed = TRUE)))
  expect_true(any(grepl("100 markers x 10000 samples", short$log, fixed = TRUE)))

  # Truncating MID-RECORD was the worse case: half a real marker plus
  # fabricated homozygotes PASSED QC and entered the GRM (19 markers, one of
  # which the intact file rejects outright). Same refusal.
  half <- setgeno_run(stage_short_bed(60, extraBytes = 1250))
  expect_false(identical(half$status, 0L))
  expect_true(any(grepl("is too short for the .bim", half$log, fixed = TRUE)))

  # POSITIVE CONTROL that the input really is inconsistent and that the two
  # readers now agree about it: step 2's #58 check refuses the same trio.
  r2 <- marker_run(stage_short_bed(60), "step1short")
  expect_false(identical(as.integer(r2$status), 0L))
  expect_true(any(grepl("is too short for the .bim and .fam", r2$log, fixed = TRUE)))
})
