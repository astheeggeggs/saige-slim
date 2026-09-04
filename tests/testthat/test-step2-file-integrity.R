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
