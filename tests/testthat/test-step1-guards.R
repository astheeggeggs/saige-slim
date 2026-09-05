# Three step-1 guards from the type-polymorphic-builtin audit, FIXED:
#   LEDGER #35  one MAC category: `1:(numCate - 1)` was 1:0 = c(1, 0), and the
#               second pass indexed listOfMarkersForVarRatio[[0]] -> crash.
#   LEDGER #36  Covariate_Transform dropped collinear columns without
#               drop = FALSE; with one column left, nrow() of a bare vector
#               was NULL and the fit died in "non-numeric argument".
#   LEDGER #37  seq(1, nrow(frame), by = 1) on an empty frame errored with
#               "wrong sign in 'by' argument", so the intended "No
#               individuals left in the analysis" could never fire.
# Each was red at unpatched HEAD in the form the row predicts. Tier: standard
# (two ~4 s fits and two instant refusals).

g_fit <- function(extra, pheno = NULL, covar = "x1,x2,a9,a10", qcovar = "a9,a10", trait = "binary", phenoCol = "y_binary") {
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  if (is.null(pheno)) pheno <- "./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"
  args <- c("--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
            paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            "--useSparseGRMtoFitNULL=TRUE", paste0("--phenoFile=", pheno), paste0("--phenoCol=", phenoCol),
            paste0("--covarColList=", covar), "--sampleIDColinphenoFile=IID", paste0("--traitType=", trait),
            "--IsOverwriteVarianceRatioFile=TRUE", paste0("--outputPrefix=", tempfile("s1g_")), extra)
  if (nzchar(qcovar)) args <- c(args, paste0("--qCovarColList=", qcovar))
  st <- run_saige_cli("step1_fitNULLGLMM.R", args)
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE))
}

test_that("LEDGER #35: a single MAC category estimates one categorical ratio instead of crashing", {
  skip_unless_level("standard")
  r <- g_fit(c("--isCateVarianceRatio=TRUE", "--cateVarRatioMinMACVecExclude=20.5", "--cateVarRatioMaxMACVecInclude="))
  expect_identical(r$status, 0L)
  expect_false(any(grepl("attempt to select less than one element", r$log, fixed = TRUE)))
})

test_that("LEDGER #37: a SampleIDIncludeFile disjoint from the phenotypes gets the intended message", {
  skip_unless_level("standard")
  ids <- tempfile(); writeLines(c("nobody1", "nobody2"), ids)
  r <- g_fit(paste0("--SampleIDIncludeFile=", ids))
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("No individuals left in the analysis", r$log, fixed = TRUE)))
  expect_false(any(grepl("wrong sign in 'by' argument", r$log, fixed = TRUE)))
})

test_that("LEDGER #36: a quantitative fit whose only covariate is constant drops it and continues", {
  skip_unless_level("standard")
  src <- file.path(saige_extdata, "input", "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
  d <- utils::read.table(src, header = TRUE, stringsAsFactors = FALSE)
  d$const <- 1
  p <- tempfile("pheno_const_", fileext = ".txt"); utils::write.table(d, p, quote = FALSE, row.names = FALSE)
  r <- g_fit(character(), pheno = p, covar = "const", qcovar = "", trait = "quantitative", phenoCol = "y_quantitative")
  expect_identical(r$status, 0L)
  expect_true(any(grepl("multi collinearity is detected in covariates", r$log, fixed = TRUE)))
  expect_false(any(grepl("non-numeric argument to mathematical function", r$log, fixed = TRUE)))
})
