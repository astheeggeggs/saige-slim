# LEDGER #97 FIXED: the `sparse` variance-ratio rows are scale-free and must
# not be rescaled with the `null` rows.
#
# The #65 standardisation fits a quantitative trait on y / s and multiplies
# every variance-ratio row by 1 / s^2 on the way out. That is right for the
# `null` and `null_noXadj` rows: their numerator G' Sigma^-1 G / AC is in
# 1 / y^2 units and their denominator g'g carries no y. It is wrong for the
# `sparse` row: its denominator g' Sigma_sparse^-1 g ALSO goes through Sigma^-1,
# so multiplying y by f multiplies both sides by 1 / f^2 and the ratio does not
# move. Applied anyway, the sparse rows came out 1 / s^2 too high -- 11-12% on
# the N = 50,000 fixture, ~1.5x on the bundled one (adjusted sd ~ 0.8) -- and
# step 2 consumes the sparse row on the region path and the sparse marker test.
#
# Two detections, derived before running:
#   1. y * 100 leaves every `sparse` row unchanged (to the rounding of the
#      scaling itself, 1e-8 relative, as in test-step1-pheno-scale.R) while
#      every `null` / `null_noXadj` row scales by exactly 1e-4. At unpatched
#      HEAD the sparse rows scaled by 1e-4 too.
#   2. With the null model fitted on the SAME sparse GRM the ratio compares
#      against, numerator and denominator are the same quadratic form up to the
#      marker sample, so the row is ~1 (binary_sparseGRM_vr's are 0.99998). At
#      unpatched HEAD the unscaled fixture's rows were ~1 / s^2 ~ 1.5.
# No committed golden has quantitative `sparse` rows (quant_sparseGRM is fitted
# without --useSparseGRMforVarRatio), which is why nothing caught it.
#
# Tier: standard. Two step-1 fits, ~10 s.

vs_fit <- function(scale) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  prefix <- file.path(tempfile(sprintf("step1_vrsparse%g_", scale)), "fit")
  dir.create(dirname(prefix), recursive = TRUE)
  st <- run_saige_cli("step1_fitNULLGLMM.R", c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed",
           ".sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    "--useSparseGRMforVarRatio=TRUE",
    paste0("--phenoFile=", ps_pheno_file(scale)),
    "--phenoCol=y_quantitative", "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID", "--traitType=quantitative",
    "--IsOverwriteVarianceRatioFile=TRUE", "--isCateVarianceRatio=TRUE",
    paste0("--outputPrefix=", prefix)
  ))
  expect_cli_ok(st)
  read.table(paste0(prefix, ".varianceRatio.txt"), header = FALSE, stringsAsFactors = FALSE)
}

vs_cache <- new.env(parent = emptyenv())
vs_run <- function(scale) {
  key <- format(scale)
  if (is.null(vs_cache[[key]])) vs_cache[[key]] <- vs_fit(scale)
  vs_cache[[key]]
}

test_that("LEDGER #97 FIXED: the sparse rows are invariant to y * 100; the null rows scale", {
  base <- vs_run(1)
  x100 <- vs_run(100)
  # the configuration reaches the sparse rows at all, in every category
  expect_true(any(base[[2]] == "sparse"))
  expect_setequal(unique(base[[3]][base[[2]] == "sparse"]), unique(base[[3]]))
  # same rows, same labels
  expect_identical(dim(x100), dim(base))
  expect_identical(x100[[2]], base[[2]])
  expect_identical(x100[[3]], base[[3]])
  sp <- base[[2]] == "sparse"
  expect_lt(rel_err(x100[[1]][sp], base[[1]][sp]), 1e-8)
  expect_lt(rel_err(x100[[1]][!sp] * 100^2, base[[1]][!sp]), 1e-8)
})

test_that("LEDGER #97 FIXED: a sparse-GRM fit's sparse ratio is ~1 in phenotype units", {
  base <- vs_run(1)
  sp <- base[[2]] == "sparse"
  expect_lt(max(abs(base[[1]][sp] - 1)), 0.05)
})
