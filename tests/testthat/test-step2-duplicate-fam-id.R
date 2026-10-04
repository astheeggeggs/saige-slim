# A duplicated IID in the genotype sample file (LEDGER #141).
#
# A repeated IID is legal PLINK1 (the key is FID+IID) but SAIGE discards FID
# everywhere, so two fam rows sharing an IID are two individuals SAIGE cannot
# tell apart. Before the fix: steps 0/1 read the fam/psam IDs in R
# (readGenoSampleIDs, R/Geno.R) with no uniqueness test, so the phenotype merge
# (by IID) emitted that phenotype row twice with two different IndexGeno and
# both fam rows were loaded as samples of the model; step 2's PLINK reader
# (readFamFile, src/PLINK.cpp) had no test either, and Rcpp::match mapped BOTH
# model slots to the FIRST fam row -- never NA, so the only guard did not fire
# -- so step 1 and step 2 used different genotype rows for the same model row.
# The PGEN reader already refused a duplicate psam IID at step 2
# (src/PGEN.cpp) but not at steps 0/1, which never build the IID map.
#
# LEDGER #67 is the same defect in the PHENOTYPE file and is refused by name;
# these tests pin the same refusal for the genotype sample file, through both
# readers, at step 1 (R) and at step 2 (C++ PLINK). Each is a fixture copied
# to a temp dir with ONE IID rewritten to another row's, so the file is
# otherwise byte-identical to a bundled one that every golden runs on.
# Tier: standard (CLI runs; each refusal happens before any fit or scan).

dfi_copy_fixture <- function(src, dst, exts) {
  for (e in exts) {
    stopifnot(file.copy(paste0(src, e), paste0(dst, e), overwrite = TRUE))
  }
  dst
}

# Rewrite the IID (column 2) of fam/psam row `row` to `iid`; the header line of
# a psam is left alone.
dfi_plant_duplicate <- function(path, row, iid) {
  lines <- readLines(path)
  isHeader <- grepl("^#", lines)
  dataRows <- which(!isHeader)
  target <- dataRows[row]
  f <- strsplit(lines[target], "[ \t]+")[[1]]
  f[2] <- iid
  sep <- if (grepl("\t", lines[target])) "\t" else " "
  lines[target] <- paste(f, collapse = sep)
  writeLines(lines, path)
  invisible(path)
}

test_that("LEDGER #141: step 2 refuses a fam with a repeated IID (it read the wrong row)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- tempfile("dupfam_")
  dir.create(d)
  src <- file.path(saige_extdata, "input", "genotype_100markers")
  pfx <- dfi_copy_fixture(src, file.path(d, "geno"), c(".bed", ".bim", ".fam"))
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sgName <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  sg <- file.path(saige_extdata, "output", sgName)
  # a GROUP run (golden P2's shape): this model's variance-ratio file has no
  # null_noXadj row, so the marker path would refuse it before the fam is read
  args <- function(out) {
    c(paste0("--bedFile=", pfx, ".bed"), paste0("--bimFile=", pfx, ".bim"),
      paste0("--famFile=", pfx, ".fam"), "--AlleleOrder=alt-first",
      "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
      paste0("--GMMATmodelFile=", m, ".rda"),
      paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
      paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
      "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE",
      paste0("--SAIGEOutputFile=", out))
  }

  # NEGATIVE CONTROL: the untouched copy runs to completion, so a refusal
  # below is the planted duplicate and not the copying.
  st0 <- run_saige_cli("step2_SPAtests.R", args(file.path(d, "ok.txt")))
  expect_cli_ok(st0)
  expect_true(file.size(file.path(d, "ok.txt")) > 0)

  # Row 5000 of the 10,000-row fam is not one of the model's 1,000 samples;
  # give it the IID of the model's first sample. Two fam rows now claim 1a1.
  load(paste0(m, ".rda"))
  dupID <- modglmm$sampleID[1]
  expect_identical(dupID, "1a1")
  dfi_plant_duplicate(paste0(pfx, ".fam"), 5000L, dupID)
  ids <- vapply(strsplit(readLines(paste0(pfx, ".fam")), "[ \t]+"), `[`, character(1), 2L)
  expect_identical(sum(ids == dupID), 2L)

  st <- run_saige_cli("step2_SPAtests.R", args(file.path(d, "dup.txt")))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  # LEDGER #141 FIXED: refused by name, before any marker is read. It used to
  # exit 0 and scan with both model slots pointing at the first fam row.
  expect_false(identical(as.integer(st), 0L))
  expect_true(any(grepl(paste0("fam file has a duplicate IID (", dupID, ")"), lg, fixed = TRUE)))
  expect_false(any(grepl("Analysis done!", lg, fixed = TRUE)))
})

dfi_step1_src <- file.path(saige_extdata, "input",
                           "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")

dfi_step1_args <- function(genoFlag, d, tag) {
  grm <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  c(genoFlag,
    paste0("--sparseGRMFile=", grm), paste0("--sparseGRMSampleIDFile=", grm, ".sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    "--phenoCol=y_binary", "--covarColList=x1,x2", "--sampleIDColinphenoFile=IID",
    "--traitType=binary", "--LOCO=FALSE", "--IsOverwriteVarianceRatioFile=TRUE",
    paste0("--outputPrefix=", file.path(d, paste0("out_", tag))))
}

test_that("LEDGER #141: step 1 refuses a fam with a repeated IID through the PLINK reader", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- tempfile("dupfam1_")
  dir.create(d)
  pfx <- dfi_copy_fixture(dfi_step1_src, file.path(d, "g"), c(".bed", ".bim", ".fam"))
  # row 2 takes row 1's IID: both are phenotyped, so before the fix the merge
  # carried 1a1 twice and the fit went ahead on 1,001 sample slots
  dfi_plant_duplicate(paste0(pfx, ".fam"), 2L, "1a1")
  st <- run_saige_cli("step1_fitNULLGLMM.R",
                      dfi_step1_args(paste0("--plinkFile=", pfx), d, "plink"))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_false(identical(as.integer(st), 0L))
  expect_true(any(grepl("1 sample ID(s) appear more than once in the fam file", lg, fixed = TRUE)))
  expect_true(any(grepl(": 1a1.", lg, fixed = TRUE)))
  # refused before the count that #67's provenance names as the only tell,
  # and before the fit. fitNULLGLMM pre-creates the model file EMPTY as its
  # writability check; LEDGER #171 FIXED (2026-10-04): a fit that does not
  # complete removes it again, so the refusal leaves NO model file.
  expect_false(any(grepl("samples will be used for analysis", lg, fixed = TRUE)))
  expect_false(any(grepl("Start fitting the NULL GLMM", lg, fixed = TRUE)))
  expect_false(file.exists(file.path(d, "out_plink.rda")))
})

test_that("LEDGER #141: step 1 refuses a psam with a repeated IID through the PGEN reader", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- tempfile("duppsam1_")
  dir.create(d)
  pfx <- dfi_copy_fixture(dfi_step1_src, file.path(d, "g"), c(".pgen", ".pvar", ".psam"))
  dfi_plant_duplicate(paste0(pfx, ".psam"), 2L, "1a1")
  st <- run_saige_cli("step1_fitNULLGLMM.R",
                      dfi_step1_args(paste0("--plink2File=", pfx), d, "pgen"))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_false(identical(as.integer(st), 0L))
  expect_true(any(grepl("1 sample ID(s) appear more than once in the psam file", lg, fixed = TRUE)))
  expect_false(any(grepl("samples will be used for analysis", lg, fixed = TRUE)))
  expect_false(any(grepl("Start fitting the NULL GLMM", lg, fixed = TRUE)))
  expect_false(file.exists(file.path(d, "out_pgen.rda")))   # LEDGER #171 FIXED
})
