# LEDGER #116 (FIXED 2026-09-13): a region run that ANALYSES a region and finds
# nothing to test no longer leaves a 0-byte results file at exit 0.
#
# The row's bundled reproduction: GENE1 = rs1 + rs2, both lof, on the frozen
# example_binary_sparseGRM model. Both markers are polymorphic across the
# 10,000 fam rows and MONOMORPHIC in the analysed 1,000 (ALT count 0, measured
# with plink2 --freq --keep), so at --minMAC=0.5 both die at the ingest QC gate
# in src/Main.cpp, nothing reaches the collapse or the cell filters, and the
# region is "analysed" with no cell to test. Before the fix the SKAT-O path
# wrote nothing (no header, no row) and the BURDEN path wrote a header plus an
# NA Cauchy row -- the two paths disagreed about the same region -- and the
# early return in regionFinish skipped the index bookkeeping, so the run's
# .index recorded "chunk 0" and no completion line.
#
# The decision: BOTH paths write a header at start and ONE row per region
# with nothing tested, "Region Cauchy NA ...", the row the BURDEN writer has
# always emitted; the #127 sidecar names the region with Rows = 0. The SKAT-O
# header and row come from R/RegionSchema.R, whose schema is asserted against
# every chunk of real rows a run writes (so every committed golden gates it).
# Tier: standard (five step-2 runs).

re_model <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
re_grm <- file.path(saige_extdata, "output",
                    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
re_cli <- run_saige_cli
re_ok <- expect_cli_ok
re_read <- read_saige_table
re_extdata <- saige_extdata

re_group_file <- function(withReal) {
  gf <- tempfile(fileext = ".txt")
  lines <- c("GENE1 var rs1 rs2", "GENE1 anno lof lof")
  if (withReal) {
    real <- readLines(file.path(re_extdata, "input", "group_new_snpid.txt"))
    lines <- c(lines, grep("^GENE2 ", real, value = TRUE))
  }
  writeLines(lines, gf)
  gf
}

re_run <- function(groupFile, extra = character(), anno = "lof", maf = "0.01") {
  out <- tempfile(fileext = ".txt")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", re_model, ".rda"),
            paste0("--varianceRatioFile=", re_model, ".varianceRatio.txt"),
            paste0("--sparseGRMFile=", re_grm),
            paste0("--sparseGRMSampleIDFile=", re_grm, ".sampleIDs.txt"),
            paste0("--groupFile=", groupFile),
            paste0("--annotation_in_groupTest=", anno),
            paste0("--maxMAF_in_groupTest=", maf), "--is_fastTest=TRUE",
            paste0("--SAIGEOutputFile=", out), extra)
  st <- re_cli("step2_SPAtests.R", args)
  re_ok(st)
  list(out = out, log = readLines(attr(st, "log"), warn = FALSE))
}

re_expect_empty_region_row <- function(out, region) {
  x <- re_read(out)
  rows <- x[x$Region == region, , drop = FALSE]
  testthat::expect_identical(nrow(rows), 1L)
  testthat::expect_identical(rows$Group, "Cauchy")
  other <- setdiff(colnames(rows), c("Region", "Group"))
  testthat::expect_true(all(as.character(unlist(rows[, other])) == "NA"))
  sidePath <- paste0(out, ".pooledTests.txt")
  testthat::expect_true(file.exists(sidePath))
  if (!file.exists(sidePath)) return(invisible(x))
  side <- read.table(sidePath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  testthat::expect_identical(side$Rows[side$Region == region], 0L)
  testthat::expect_identical(side$Rule[side$Region == region], "none")
  invisible(x)
}

test_that("LEDGER #116: a run whose only region has nothing to test writes a header + NA row", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  gf <- re_group_file(withReal = FALSE)
  for (extra in list(character(), "--r.corr=1")) {
    r <- re_run(gf, extra)
    # the row's symptom: analysed, "Analysis done!", and nothing on disk
    expect_true(any(grepl("Analyzing Region GENE1", r$log, fixed = TRUE)))
    expect_true(any(grepl("Analysis done", r$log, fixed = TRUE)))
    expect_gt(file.size(r$out), 0)
    x <- re_expect_empty_region_row(r$out, "GENE1")
    expect_identical(nrow(x), 1L)
    # the ultra-rare tell line is ABSENT: nothing passed QC, so this is the
    # ingest-gate route and not the collapse or the cell filters
    expect_false(any(grepl("markers are ultra-rare and", r$log, fixed = TRUE)))
    # and the index bookkeeping ran: chunk 1 recorded, the run marked complete
    idx <- readLines(paste0(r$out, ".index"), warn = FALSE)
    expect_true(any(grepl("^Have completed the analysis of chunk 1", idx)))
    expect_identical(tail(idx, 1), "Have completed the analyses of all chunks.")
  }
})

test_that("LEDGER #116: the empty region's row sits beside a real region's rows, both paths", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  gf <- re_group_file(withReal = TRUE)
  ref <- re_run(re_group_file(withReal = FALSE), character())   # GENE1 alone, SKAT-O
  for (extra in list(character(), "--r.corr=1")) {
    r <- re_run(gf, extra, anno = "lof,missense:lof", maf = "0.001,0.01")
    x <- re_expect_empty_region_row(r$out, "GENE1")
    g2 <- x[x$Region == "GENE2", , drop = FALSE]
    expect_true(nrow(g2) > 1)
    expect_identical(sum(g2$Group == "Cauchy"), 1L)
    # the real region's numbers are untouched by the NA row in the same
    # chunk: a real p-value is printed as a number, not coerced
    pb <- g2$Pvalue_Burden[g2$Group != "Cauchy"]
    expect_true(any(grepl("^[0-9.]+E-[0-9]+$|^0\\.[0-9]+$", pb)))
  }
  expect_gt(file.size(ref$out), 0)
})

test_that("LEDGER #116: a run with nothing tested under the multi-weight schema, header + NA row", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  gf <- tempfile(fileext = ".txt")
  writeLines(c("GENE1 var rs1 rs2", "GENE1 anno lof lof", "GENE1 weight:CADD 1 2"), gf)
  # K = 1 keeps the legacy schema; the forced column is what turns it on
  r <- re_run(gf, "--is_output_weight_column=TRUE")
  x <- read_saige_table(r$out)
  expect_true("Weight" %in% colnames(x))
  expect_true(all(c("Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(x)))
  re_expect_empty_region_row(r$out, "GENE1")
})
