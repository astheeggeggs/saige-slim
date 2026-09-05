# LEDGER #43 FIXED: extractNglmm.R forwards --noEstFixedEff to getNeff(). It
# was declared, printed back, and replaced by a hardcoded FALSE at the call,
# so the fixed effects were always estimated. The only observable is the
# path getNeff takes, which it announces; Nglmm itself does not depend on the
# coefficients (see the row). Tier: standard (one extractNglmm run).

test_that("--noEstFixedEff=TRUE reaches getNeff", {
  skip_unless_level("standard")
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  st <- run_saige_cli("extractNglmm.R", c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    "--phenoCol=y_binary", "--covarColList=x1,x2", "--sampleIDColinphenoFile=IID",
    "--traitType=binary", paste0("--outputPrefix=", tempfile("nglmm_")),
    "--noEstFixedEff=TRUE"))
  expect_cli_ok(st)
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("noEstFixedEff=TRUE, so fixed effects coeff", lg)))
})
