# LEDGER #148 (FIXED 2026-09-14): a quantitative getNeff run with
# --noEstFixedEff=TRUE no longer dies on an unassigned `covoffset`.
#
# Before the fix the flag's binary arm assigned the covariate offset and the
# quantitative arm residualised the phenotype instead and assigned nothing,
# while the quantitative fit0 read `offset = covoffset` -- so the run died
# after the GRM and phenotype preparation with `object 'covoffset' not
# found`. Both arms now build the offset the same way (one glm, family by
# trait, offset = its linear predictor), which is fitNULLGLMM's shape.
#
# Nglmm itself cannot see the choice: both arms solve 1' Sigma^-1 1 at the
# fixed tau = (0, 1) and read fit0 only for length(fit0$y) (see the row, and
# the Nglmm-magnitude observation recorded beside it). So the observable is
# that the run completes, announces the path, and prints the SAME Nglmm as
# the run without the flag. Tier: standard (two extractNglmm runs).

# helper-saige.R supplies the CLI runner; fetched with get() so the linter can
# see the binding (the idiom of helper-step2-flexrv.R)
l148_cli <- function(...) get("run_saige_cli")(...)

l148_run <- function(flag) {
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  args <- c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    "--phenoCol=y_quantitative", "--covarColList=x1,x2", "--sampleIDColinphenoFile=IID",
    "--traitType=quantitative", paste0("--outputPrefix=", tempfile("nglmm_q_")),
    if (flag) "--noEstFixedEff=TRUE" else character()
  )
  st <- l148_cli("extractNglmm.R", args)
  lg <- readLines(attr(st, "log"), warn = FALSE)
  ng <- grep("^Nglmm ", lg, value = TRUE)
  list(status = as.integer(st), log = lg,
       nglmm = if (length(ng) == 1) as.numeric(sub("^Nglmm ", "", ng)) else NA_real_)
}

test_that("quantitative getNeff with --noEstFixedEff=TRUE completes (LEDGER #148)", {
  skip_unless_level("standard")
  a <- l148_run(TRUE)
  expect_identical(a$status, 0L)
  expect_true(any(grepl("noEstFixedEff=TRUE, so fixed effects coeff", a$log)))
  expect_false(any(grepl("object 'covoffset' not found", a$log)))
  expect_false(is.na(a$nglmm))
  # the control: same Nglmm without the flag (fit0 is read only for its length)
  b <- l148_run(FALSE)
  expect_identical(b$status, 0L)
  expect_equal(a$nglmm, b$nglmm, tolerance = 1e-10)
})
