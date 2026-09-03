# PGEN (plink2) input for step 0 and step 1.
#
# The claim under test is BYTE-IDENTITY against the PLINK goldens, which is
# ordinarily forbidden here (see CLAUDE.md, "Never assert bit-identity between
# two different code paths"). It is legitimate in this one case because the two
# readers are not two code paths: the pgen front-end decodes a variant with
# pgenlib and RE-ENCODES it as a PLINK-1 .bed record, which is then handed to
# the same Get_OneSNP_Geno_atBeginning as the .bed reader's own bytes. The
# unpack, the MAF/missingness QC, the round(2*altFreq) imputation, the 2-bit
# repack and every float downstream are literally one code path with two
# front-ends. If any of these ever needs a tolerance, the re-encode has broken
# -- that is the finding; do not relax the assertion.
#
# What makes the two readers comparable at all was verified when the fixtures
# were made: pvar ALT == bim column 5 and pvar REF == bim column 6 for every
# marker (10,000 in the step-0 fixture, 128,868 in the step-1 one), and psam
# IID order == fam IID order. Step 1 counts copies of bim column 5; pgenlib's
# alleleIdx = 1 reads ALT. Same quantity.

test_that("psam and fam yield the same sample IDs, in the same order", {
  for (pfx in c("plinkforGRM_1000samples_10kMarkers",
                "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")) {
    base <- file.path(saige_extdata, "input", pfx)
    fam <- SAIGE:::readGenoSampleIDs("plink", famFile = paste0(base, ".fam"))
    psam <- SAIGE:::readGenoSampleIDs("pgen", psamFile = paste0(base, ".psam"))
    expect_identical(psam, fam)
  }
})

test_that("pvar and bim yield the same chromosome vector, in the same order", {
  for (pfx in c("plinkforGRM_1000samples_10kMarkers",
                "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")) {
    base <- file.path(saige_extdata, "input", pfx)
    bim <- SAIGE:::readVariantChrom("plink", bimFile = paste0(base, ".bim"))
    pvar <- SAIGE:::readVariantChrom("pgen", pvarFile = paste0(base, ".pvar"))
    expect_identical(as.character(pvar), as.character(bim))
  }
})

test_that("checkGenoInputStep1 resolves prefixes and refuses mixed input", {
  base <- file.path(saige_extdata, "input", "plinkforGRM_1000samples_10kMarkers")

  p <- SAIGE:::checkGenoInputStep1(plinkFile = base)
  expect_identical(p$genoType, "plink")
  expect_identical(p$bedFile, paste0(base, ".bed"))

  g <- SAIGE:::checkGenoInputStep1(plink2File = base)
  expect_identical(g$genoType, "pgen")
  expect_identical(g$pgenFile, paste0(base, ".pgen"))
  expect_identical(g$psamFile, paste0(base, ".psam"))

  # pgenFile alone derives its two companions.
  g2 <- SAIGE:::checkGenoInputStep1(pgenFile = paste0(base, ".pgen"))
  expect_identical(g2$pvarFile, paste0(base, ".pvar"))

  expect_error(
    SAIGE:::checkGenoInputStep1(plinkFile = base, plink2File = base),
    "mutually exclusive"
  )
  expect_error(SAIGE:::checkGenoInputStep1(), "required")
  # Step 1 may legitimately be given no genotype file at all.
  expect_identical(SAIGE:::checkGenoInputStep1(allowNone = TRUE)$genoType, "")
})

test_that("the pgen reader loads a GRM genotype store identical to the bed one", {
  skip_unless_level("standard")
  base <- file.path(saige_extdata, "input", "plinkforGRM_1000samples_10kMarkers")

  # Driven in child processes: setgeno populates a C++ GLOBAL, so loading both
  # formats in one session would have the second overwrite the first.
  dump <- function(mode) {
    out <- tempfile(fileext = ".rds")
    script <- tempfile(fileext = ".R")
    writeLines(c(
      "suppressPackageStartupMessages(library(SAIGE))",
      sprintf('base <- "%s"; mode <- "%s"; out <- "%s"', base, mode, out),
      "setminMAFforGRM(0.01); setmaxMissingRateforGRM(0.15)",
      'ids <- if (mode == "plink") {',
      '  SAIGE:::readGenoSampleIDs("plink", famFile = paste0(base, ".fam"))',
      "} else {",
      '  SAIGE:::readGenoSampleIDs("pgen", psamFile = paste0(base, ".psam"))',
      "}",
      "n <- length(ids); idx <- seq_len(n); ind <- rep(TRUE, n)",
      'if (mode == "plink") {',
      '  setgeno(paste0(base, ".bed"), paste0(base, ".bim"),',
      '          paste0(base, ".fam"), idx, ind, 2, FALSE)',
      "} else {",
      '  setgeno_pgen(paste0(base, ".pgen"), paste0(base, ".pvar"),',
      '               paste0(base, ".psam"), idx, ind, 2, FALSE)',
      "}",
      "saveRDS(list(M = gettotalMarker(), N = getNnomissingOut(),",
      "             freq = getAlleleFreqVec(), mac = getMACVec(),",
      "             qcd = getQCdMarkerIndex(),",
      "             g = sapply(0:49, function(i) Get_OneSNP_Geno(i))), out)"
    ), script)
    st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                   stdout = FALSE, stderr = FALSE))
    expect_identical(as.integer(st), 0L)
    readRDS(out)
  }

  a <- dump("plink")
  b <- dump("pgen")

  expect_identical(a$M, b$M)
  expect_identical(a$N, b$N)
  # Bit-for-bit, not toleranced: these are the inputs to everything else.
  expect_identical(a$freq, b$freq)
  expect_identical(a$mac, b$mac)
  expect_identical(a$qcd, b$qcd)
  expect_identical(a$g, b$g)
})

test_that("step 0 through the pgen reader is byte-identical to the plink golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  suffix <- "_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
  out <- file.path(tempdir(), "s0_pgen")
  # set.seed(1) lives in extdata/createSparseGRM.R, and it is load-bearing:
  # 5,000 markers are sampled out of the 9,588 that pass MAF, so two UNSEEDED
  # runs differ from each other and this comparison would be measuring the RNG.
  # Driving the CLI (not createSparseGRM directly) is what keeps the seed in.
  st <- run_saige_cli("createSparseGRM.R", c(
    paste0("--plink2File=input/plinkforGRM_1000samples_10kMarkers"),
    "--nThreads=1", paste0("--outputPrefix=", out),
    "--numRandomMarkerforSparseKin=5000", "--relatednessCutoff=0.05"
  ))
  expect_cli_ok(st)

  golden <- file.path(saige_goldens, "step0", "test_sparseGRM")
  expect_files_byte_identical(paste0(golden, suffix), paste0(out, suffix))
  expect_files_byte_identical(paste0(golden, suffix, ".sampleIDs.txt"),
                              paste0(out, suffix, ".sampleIDs.txt"))
})

test_that("step 1 through the pgen reader is byte-identical to the plink goldens", {
  skip_unless_level("full")
  skip_if_no_goldens()

  pfx <- "./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr"
  grm <- c(
    "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE"
  )
  pheno <- "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"
  common <- c("--sampleIDColinphenoFile=IID", "--IsOverwriteVarianceRatioFile=TRUE")

  # The five step-1 golden configurations. binary_denseGRM is the only one on
  # the DENSE GRM path, so dropping it would leave the TBB kernel reached by no
  # pgen test at all.
  cfgs <- list(
    binary_sparseGRM_a = c(grm, pheno, "--phenoCol=y_binary",
                           "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
                           "--traitType=binary", "--isCateVarianceRatio=TRUE", common),
    binary_sparseGRM_b = c(grm, pheno, "--phenoCol=y_binary",
                           "--covarColList=x1,x2", "--qCovarColList=x2",
                           "--cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5",
                           "--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,20.5",
                           "--traitType=binary", "--isCateVarianceRatio=TRUE", common),
    quant_sparseGRM = c(grm, pheno, "--phenoCol=y_quantitative",
                        "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
                        "--traitType=quantitative", "--invNormalize=TRUE",
                        "--isCateVarianceRatio=TRUE", common),
    binary_denseGRM = c(pheno, "--phenoCol=y_binary",
                        "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
                        "--traitType=binary", "--nThreads=1", "--LOCO=FALSE",
                        "--minMAFforGRM=0.01", common),
    binary_sparseGRM_vr = c(grm, "--useSparseGRMforVarRatio=TRUE", pheno,
                            "--phenoCol=y_binary", "--covarColList=x1,x2,a9,a10",
                            "--qCovarColList=a9,a10", "--traitType=binary",
                            "--isCateVarianceRatio=TRUE", common)
  )

  for (nm in names(cfgs)) {
    out <- file.path(tempdir(), paste0("s1_pgen_", nm))
    st <- run_saige_cli("step1_fitNULLGLMM.R",
                        c(paste0("--plink2File=", pfx), cfgs[[nm]],
                          paste0("--outputPrefix=", out)))
    expect_cli_ok(st)
    golden <- file.path(saige_goldens, "step1", nm)
    expect_files_byte_identical(paste0(golden, ".rda"), paste0(out, ".rda"))
    expect_files_byte_identical(paste0(golden, ".varianceRatio.txt"),
                                paste0(out, ".varianceRatio.txt"))
  }
})

test_that("a pgen carrying dosages is hard-called, and says so once", {
  # extdata/input/dosage_20markers_100samples.{pgen,pvar,psam}: 20 markers x
  # 100 samples, built with
  #   plink2 --gen <oxford .gen> ref-first --sample <.sample> --make-pgen
  # from probability triples chosen so the ALT dosage cycles 0.1, 0.2, ... 1.0
  # across samples. 1,800 of the 2,000 genotypes are fractional and 200 (the
  # dosage-1.0 ones) are exactly integral -- so the fixture also pins that an
  # integral dosage is NOT miscounted as hard-called.
  base <- file.path(saige_extdata, "input", "dosage_20markers_100samples")
  skip_if_not(file.exists(paste0(base, ".pgen")), "dosage fixture not present")

  out <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('base <- "%s"; out <- "%s"', base, out),
    "setminMAFforGRM(0); setmaxMissingRateforGRM(1)",
    'ids <- SAIGE:::readGenoSampleIDs("pgen", psamFile = paste0(base, ".psam"))',
    "n <- length(ids); msg <- NULL",
    "withCallingHandlers(",
    '  setgeno_pgen(paste0(base, ".pgen"), paste0(base, ".pvar"), paste0(base, ".psam"),',
    "               seq_len(n), rep(TRUE, n), 2, FALSE),",
    "  warning = function(w) { msg <<- c(msg, conditionMessage(w));",
    '                          invokeRestart("muffleWarning") })',
    "saveRDS(list(msg = msg, M = gettotalMarker(), n = n,",
    "             g = Get_OneSNP_Geno(0)), out)"
  ), script)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = FALSE, stderr = FALSE))
  expect_identical(as.integer(st), 0L)
  r <- readRDS(out)

  # Exactly one warning, naming the count of hard-called genotypes.
  hits <- grep("hard-called", r$msg, value = TRUE)
  expect_length(hits, 1L)
  expect_match(hits, "1800 genotype")

  # Hard calls are the nearest integer: dosages 0.1-0.4 -> 0, 0.5-1.0 -> 1,
  # i.e. 40 zeros and 60 ones per marker on this fixture.
  expect_identical(sort(unique(as.integer(r$g))), c(0L, 1L))
  expect_identical(as.integer(table(as.integer(r$g))), c(40L, 60L))
})
