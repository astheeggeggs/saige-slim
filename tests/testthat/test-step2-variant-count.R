# LEDGER #87 FIXED: the variant file must describe exactly the genotype
# file's marker set.
#
# Marker metadata comes from the text file by row index and the genotype from
# the binary file by the same index, and nothing tied the two counts together.
# A .bim or .pvar missing an INTERIOR row relabelled every marker after it at
# exit 0 -- measured: rs5 deleted, 70 rows, 68 of them carrying the previous
# marker's statistics. PLINK's #58 check was an inequality (a longer .bed was
# allowed), which is exactly the case a missing row produces; PGEN never
# compared the header's raw_variant_ct with the pvar's rows. Both are exact
# now, in step 2 and in step 1's BedByteSource.
#
# Fixtures are DERIVED: the committed .bed/.pgen are reused unchanged and only
# the text file loses its rs5 line. Tier: standard (four short runs).

short_variant_file <- function(src, is_pvar) {
  lines <- readLines(src, warn = FALSE)
  id_col <- if (is_pvar) 3 else 2
  ids <- vapply(strsplit(lines, "\t| "), function(f) if (length(f) >= id_col) f[id_col] else "", "")
  keep <- !(ids == "rs5")
  expect_identical(sum(!keep), 1L)
  out <- tempfile(fileext = if (is_pvar) ".pvar" else ".bim")
  writeLines(lines[keep], out)
  out
}

m1 <- function() file.path(saige_goldens, "step1", "binary_sparseGRM_a")

step2_args <- function(out, m = m1()) {
  c(paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", m, ".rda"), paste0("--varianceRatioFile=", m, ".varianceRatio.txt"))
}

test_that("PLINK step 2: a .bim missing an interior row is refused, not relabelled", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  bim <- short_variant_file(file.path(saige_extdata, "input", "genotype_100markers.bim"), FALSE)
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed", paste0("--bimFile=", bim),
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first", step2_args(out)))
  expect_false(identical(as.integer(st), 0L))
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("is LONGER than the .bim and .fam", log, fixed = TRUE)))
  expect_false(file.exists(out) && file.size(out) > 0)
})

test_that("PGEN step 2: a .pvar missing an interior row is refused by the header count", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  pvar <- short_variant_file(file.path(saige_extdata, "input", "genotype_100markers.pvar"), TRUE)
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", c(
    "--pgenFile=./input/genotype_100markers.pgen", paste0("--pvarFile=", pvar),
    "--psamFile=./input/genotype_100markers.psam", step2_args(out)))
  expect_false(identical(as.integer(st), 0L))
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("holds 100 variants but the .pvar", log, fixed = TRUE)))
  expect_true(any(grepl("describes 99", log, fixed = TRUE)))
})

test_that("step 1: a .bim missing an interior row is refused before any GRM is built", {
  skip_unless_level("standard")
  src <- file.path(saige_extdata, "input", "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")
  lines <- readLines(paste0(src, ".bim"), warn = FALSE)
  dir <- tempfile("s1short_"); dir.create(dir)
  prefix <- file.path(dir, "geno")
  writeLines(lines[-5], paste0(prefix, ".bim"))
  file.symlink(normalizePath(paste0(src, ".bed")), paste0(prefix, ".bed"))
  file.symlink(normalizePath(paste0(src, ".fam")), paste0(prefix, ".fam"))
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  st <- run_saige_cli("step1_fitNULLGLMM.R", c(
    paste0("--plinkFile=", prefix), paste0("--sparseGRMFile=", sg),
    paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"), "--useSparseGRMtoFitNULL=TRUE",
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    "--phenoCol=y_binary", "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID", "--traitType=binary", "--IsOverwriteVarianceRatioFile=TRUE",
    "--isCateVarianceRatio=TRUE", paste0("--outputPrefix=", file.path(dir, "fit"))))
  expect_false(identical(as.integer(st), 0L))
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("is LONGER than the .bim", log, fixed = TRUE)))
})

test_that("POSITIVE CONTROL: the intact trio still runs and matches golden P1", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first", step2_args(out),
    "--is_output_moreDetails=TRUE"))
  expect_cli_ok(st)
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_marker_binary.txt"), out)
})
