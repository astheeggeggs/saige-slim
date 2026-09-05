# LEDGER #77 FIXED: --FemaleOnly / --MaleOnly could not work unless the sex
# column happened to be a covariate, and the error blamed the wrong file.
#
# The phenotype frame was subset to (phenoCol, covarColList, sampleID) BEFORE
# the sex check, which then tested the SUBSET and reported that the sex column
# "does not exist in the phenoFile" -- for a file that contained it. The fix
# checks the file and keeps the column through the subset.
#
# The bundled phenotype file has no sex column, and its only 0/1 covariate
# with a sensible split (a1) has 20 ones -- too few to fit. So a sex column is
# DERIVED at run time: 1 for even rows, 2 for odd, 500 each, and the model is
# the committed binary sparse-GRM configuration otherwise unchanged. Tier:
# standard (four step-1 fits, ~4 s each on the 1,000-sample fixture).

sex_pheno <- local({
  path <- NULL
  function() {
    if (is.null(path)) {
      src <- file.path(saige_extdata, "input",
                       "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
      d <- utils::read.table(src, header = TRUE, stringsAsFactors = FALSE)
      d$sex <- ifelse(seq_len(nrow(d)) %% 2 == 0, 1L, 2L)
      path <<- tempfile("pheno_sex_", fileext = ".txt")
      utils::write.table(d, path, quote = FALSE, row.names = FALSE)
    }
    path
  }
})

sex_run <- function(prefix, extra = character(), covar = "x1,x2,a9,a10") {
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  args <- c("--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
            paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            "--useSparseGRMtoFitNULL=TRUE", paste0("--phenoFile=", sex_pheno()),
            "--phenoCol=y_binary", paste0("--covarColList=", covar), "--qCovarColList=a9,a10",
            "--sampleIDColinphenoFile=IID", "--traitType=binary",
            "--IsOverwriteVarianceRatioFile=TRUE", "--isCateVarianceRatio=TRUE",
            paste0("--outputPrefix=", prefix), extra)
  st <- run_saige_cli("step1_fitNULLGLMM.R", args)
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE))
}

model_ids <- function(rda) {
  e <- new.env(); load(rda, envir = e)
  as.character(e$modglmm$sampleID)
}

expected_ids <- function(code) {
  d <- utils::read.table(sex_pheno(), header = TRUE, stringsAsFactors = FALSE)
  as.character(d$IID[d$sex == code])
}

test_that("--FemaleOnly runs with a sex column that is NOT a covariate, on exactly the females", {
  skip_unless_level("standard")
  prefix <- tempfile("s1_female_")
  r <- sex_run(prefix, c("--FemaleOnly=TRUE", "--sexCol=sex", "--FemaleCode=1"))
  expect_identical(r$status, 0L)
  expect_true(any(grepl("^500  samples will be used for analysis", r$log)))
  # LEDGER #90 FIXED: the driver's _FemaleOnly / _MaleOnly suffix now reaches
  # the output names (it used to be appended after they were built).
  rda <- paste0(prefix, "_FemaleOnly.rda")
  expect_false(file.exists(paste0(prefix, ".rda")))
  expect_true(file.exists(rda))
  expect_setequal(model_ids(rda), expected_ids(1L))
})

test_that("--MaleOnly picks the complementary half", {
  skip_unless_level("standard")
  prefix <- tempfile("s1_male_")
  r <- sex_run(prefix, c("--MaleOnly=TRUE", "--sexCol=sex", "--MaleCode=2"))
  expect_identical(r$status, 0L)
  expect_true(any(grepl("^500  samples will be used for analysis", r$log)))
  expect_setequal(model_ids(paste0(prefix, "_MaleOnly.rda")), expected_ids(2L))
  expect_true(file.exists(paste0(prefix, "_MaleOnly.varianceRatio.txt")))
})

test_that("NEGATIVE CONTROL: a sex column that really is absent is still refused, blaming the file", {
  skip_unless_level("standard")
  r <- sex_run(tempfile("s1_nosex_"), c("--FemaleOnly=TRUE", "--sexCol=nosuchcol", "--FemaleCode=1"))
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("column for sex nosuchcol does not exist in the phenoFile", r$log, fixed = TRUE)))
})

test_that("POSITIVE CONTROL: without a sex restriction the same file fits all 1,000", {
  skip_unless_level("standard")
  r <- sex_run(tempfile("s1_all_"))
  expect_identical(r$status, 0L)
  expect_true(any(grepl("^1000  samples will be used for analysis", r$log)))
})
