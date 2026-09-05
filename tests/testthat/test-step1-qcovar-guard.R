# LEDGER #78 FIXED: a column named in --qCovarColList was made a factor with no
# cardinality check, so a continuous column became one dummy per distinct
# value (upstream: 43,790 levels, 216 GB, no output). Now: every sample its
# own value -> refused by name; more than 20 levels -> accepted with a warning
# saying how many dummies. Tier: standard (two ~4 s fits; the refusal is
# instant).

qc_fit <- function(pheno, covar, qcovar) {
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  st <- run_saige_cli("step1_fitNULLGLMM.R", c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE", paste0("--phenoFile=", pheno), "--phenoCol=y_binary",
    paste0("--covarColList=", covar), paste0("--qCovarColList=", qcovar),
    "--sampleIDColinphenoFile=IID", "--traitType=binary", "--IsOverwriteVarianceRatioFile=TRUE",
    "--isCateVarianceRatio=TRUE", paste0("--outputPrefix=", tempfile("s1qc_"))))
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE))
}

test_that("a 30-level column is accepted with a warning that counts its dummies", {
  skip_unless_level("standard")
  src <- file.path(saige_extdata, "input", "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
  d <- utils::read.table(src, header = TRUE, stringsAsFactors = FALSE)
  d$site <- rep_len(1:30, nrow(d))
  p30 <- tempfile("pheno30_", fileext = ".txt")
  utils::write.table(d, p30, quote = FALSE, row.names = FALSE)
  w <- qc_fit(p30, "x1,x2,site", "site")
  expect_identical(w$status, 0L)
  expect_true(any(grepl("covariate site is categorical with 30 levels, i.e. 29 dummy columns", w$log, fixed = TRUE)))
  expect_true(any(grepl("^Warning message", w$log)))
})

test_that("a continuous column named categorical is refused by name", {
  # DO NOT drill this block at unpatched HEAD: without the guard it asks for a
  # fit with 1,000 dummy columns, which is the pathological run the row is
  # about and does not finish in useful time. The block above is the drill.
  skip_unless_level("standard")
  if (nzchar(Sys.getenv("SAIGE_DRILL_HEAD"))) skip("pathological at unpatched HEAD; the warning block is the drill")
  src <- file.path(saige_extdata, "input", "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
  r <- qc_fit(src, "x1,x2,a9,a10", "x1,a9,a10")   # x1 is continuous: 1,000 distinct values
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("covariate x1 is named in qCovarColList as categorical but has 1000 distinct values", r$log, fixed = TRUE)))
})
