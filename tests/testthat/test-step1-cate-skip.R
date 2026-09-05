# LEDGER #83 FIXED: a MAC category the caller asks NOT to estimate
# (cateVarRatioIndexVec entry 0, R API only -- there is no CLI option) used to
# decide whether to emit a `sparse` row from the PREVIOUS category's vector,
# or from an unbound symbol at k = 1. The decision is now the run's: a sparse
# row (value 1) is written for a skipped category exactly when sparse-GRM
# ratios are being estimated. Driven through a child Rscript calling
# fitNULLGLMM directly, since no CLI reaches the argument. Tier: standard.

cate_fit <- function(sparseVR) {
  prefix <- tempfile("s1cate_")
  script <- tempfile(fileext = ".R")
  writeLines(c(
    'suppressPackageStartupMessages(library(SAIGE))',
    'args <- commandArgs(trailingOnly = TRUE)',
    'sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"',
    'fitNULLGLMM(plinkFile = "./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",',
    '  sparseGRMFile = sg, sparseGRMSampleIDFile = paste0(sg, ".sampleIDs.txt"),',
    '  useSparseGRMtoFitNULL = TRUE, useSparseGRMforVarRatio = as.logical(args[2]),',
    '  phenoFile = "./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",',
    '  phenoCol = "y_binary", covarColList = c("x1", "x2", "a9", "a10"), qCovarCol = c("a9", "a10"),',
    '  sampleIDColinphenoFile = "IID", traitType = "binary", IsOverwriteVarianceRatioFile = TRUE,',
    '  isCateVarianceRatio = TRUE, cateVarRatioIndexVec = c(0, 1), outputPrefix = args[1])'
  ), script)
  log <- tempfile(fileext = ".log")
  owd <- setwd(saige_extdata); on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(script, prefix, if (sparseVR) "TRUE" else "FALSE"), stdout = log, stderr = log))
  list(status = as.integer(st), log = readLines(log, warn = FALSE),
       vr = if (file.exists(paste0(prefix, ".varianceRatio.txt")))
              utils::read.table(paste0(prefix, ".varianceRatio.txt"), stringsAsFactors = FALSE) else NULL)
}

test_that("a skipped first category writes ratio-1 rows, with a sparse row iff sparse ratios are estimated", {
  skip_unless_level("standard")
  for (sparseVR in c(FALSE, TRUE)) {
    r <- cate_fit(sparseVR)
    expect_identical(r$status, 0L)
    expect_false(any(grepl("object 'varRatio_sparseGRM_vec' not found", r$log, fixed = TRUE)))
    vr <- r$vr
    expect_false(is.null(vr))
    k1 <- vr[vr$V3 == 1, ]
    expect_identical(sort(k1$V2), sort(c("null", "null_noXadj", if (sparseVR) "sparse")))
    expect_true(all(k1$V1 == 1))
    k2 <- vr[vr$V3 == 2, ]
    expect_true(all(c("null", "null_noXadj") %in% k2$V2))
    expect_identical("sparse" %in% k2$V2, sparseVR)
    expect_true(all(k2$V1 != 1))
  }
})
