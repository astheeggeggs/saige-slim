# LEDGER #100 FIXED: getNeff (extractNglmm.R) builds its phenotype frame with
# preparePhenotypeFrame, the function fitNULLGLMM uses, instead of an older
# copy of that block that accepted a duplicated sample ID (#67), never turned
# qCovarCol into factors (its roxygen said "NOT work yet") and had no
# SampleIDIncludeFile handling. Two arms:
#   1. a phenotype file with one IID repeated: extractNglmm must now refuse
#      it by name. At unpatched HEAD it ran to exit 0 on 1001 rows.
#   2. the numeric-covariate binary run of cmd_extractNeff.sh: Nglmm must be
#      what the old block produced, 142.2884 (measured at 33a470c5, before
#      the change), over 1000 samples -- the shared frame carries an extra
#      intercept column and the same covariates, so nothing downstream moves.
#
# Tier: standard (two extractNglmm runs, ~20 s each).

np_args <- function(phenoFile, prefix) {
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  c("--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    paste0("--phenoFile=", phenoFile),
    "--phenoCol=y_binary", "--covarColList=x1,x2", "--sampleIDColinphenoFile=IID",
    "--traitType=binary", paste0("--outputPrefix=", prefix))
}

np_pheno <- file.path(saige_extdata, "input",
                      "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")

test_that("LEDGER #100 FIXED: extractNglmm refuses a duplicated phenotype ID", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  ph <- readLines(np_pheno)
  dup <- tempfile("pheno_dup_", fileext = ".txt")
  writeLines(c(ph, ph[2]), dup) # the first data row twice
  st <- run_saige_cli("extractNglmm.R", np_args(dup, tempfile("nglmm_dup_")))
  expect_false(identical(as.integer(st), 0L))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("1 sample ID(s) appear more than once in the phenoFile column IID",
                        lg, fixed = TRUE)))
})

test_that("LEDGER #100: the numeric-covariate Nglmm is unchanged by the shared frame", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  st <- run_saige_cli("extractNglmm.R", np_args(np_pheno, tempfile("nglmm_")))
  expect_cli_ok(st)
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("^1000  samples will be used for analysis", lg)))
  n <- grep("^Nglmm ", lg, value = TRUE)
  expect_length(n, 1L)
  expect_equal(as.numeric(sub("^Nglmm +", "", n)), 142.2884, tolerance = 1e-6)
})
