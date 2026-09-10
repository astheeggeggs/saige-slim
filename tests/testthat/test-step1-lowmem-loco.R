# --isLowMemLOCO writes the same LOCO results as the in-memory layout, for
# both traits.
#
# A LOCO fit refits the fixed effects per left-out chromosome with
# Get_Coef_LOCO. In the default layout that happens inside the fitter
# (glmmkin.ai_PCG_Rcpp_{Binary,Quantitative}) and lands in
# modglmm$LOCOResult[[j]]; under --isLowMemLOCO fitNULLGLMM strips the model,
# runs the same refit itself and saves one <prefix>_chr<j>.rda per chromosome.
# The two layouts should agree per chromosome: same coefficients, fitted
# values and cov. This guards the block fitNULLGLMM runs for BOTH traits,
# which no golden reaches (the sparse-GRM goldens have LOCO forced off, the
# dense one is LOCO=FALSE) and which the 2026-09-05 fold of the
# binary/quantitative halves rewrote. It is not a bit-identity claim; the
# fold itself was drilled by cmp of the model files against the pre-fold
# build (identical for both traits and both layouts).
#
# Fixture: the two-chromosome dense-GRM subset from test-step1-stale-model-file.R
# (helper-pheno-scale.R), ~3 s per fit. Tier: standard.

ll_fit <- function(traitType, prefix, plink, lowmem) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  st <- run_saige_cli("step1_fitNULLGLMM.R", c(
    paste0("--plinkFile=", plink),
    "--nThreads=1", "--minMAFforGRM=0.01", "--LOCO=TRUE",
    paste0("--isLowMemLOCO=", if (lowmem) "TRUE" else "FALSE"),
    "--skipVarianceRatioEstimation=TRUE",
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    paste0("--phenoCol=", if (traitType == "binary") "y_binary" else "y_quantitative"),
    "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID", paste0("--traitType=", traitType),
    paste0("--outputPrefix=", prefix)
  ))
  expect_cli_ok(st)
  invisible(st)
}

ll_load <- function(path) {
  e <- new.env()
  load(path, e)
  e$modglmm
}

# `tol` is the per-chromosome agreement the trait has TODAY. Quantitative: the
# two layouts agree to 1e-10 on every chromosome, obj.noK included (the
# gaussian refit is one linear step, so it cannot depend on where it starts).
# Binary: chr1 agrees to 1e-10, but chr2 only to ~1e-5 (coefficients 8.5e-6,
# fitted values 2.4e-5, cov 7.7e-6 on this fixture, measured 2026-09-05 at
# c0688cc9 and identical after the fold) and its obj.noK$S_a -- the score at
# the refit, near zero when converged -- not at all (rel_err 1.0). The cause
# is in the code: the fitter's own LOCO loop (glmmkin.ai_PCG_Rcpp_Binary)
# passes `alpha`/`eta` that it OVERWRITES with each chromosome's result, so
# chr2 warm-starts from chr1's refit, while the low-memory loop starts every
# chromosome from the full model's alpha0/eta0; at tol.coef = 0.02 the two
# starts stop at different points. Pinned as found (see docs/state/next.md)
# so that a fix inverts the chr2 assertions.
ll_check_trait <- function(traitType, plink, chr_tol, sa_agrees) {
  d <- tempfile(paste0("lowmem_", traitType, "_"))
  dir.create(d)
  ll_fit(traitType, file.path(d, "mem"), plink, lowmem = FALSE)
  ll_fit(traitType, file.path(d, "low"), plink, lowmem = TRUE)
  mem <- ll_load(file.path(d, "mem.rda"))
  expect_true(isTRUE(mem$LOCO))
  low0 <- ll_load(file.path(d, "low_noLOCO.rda"))
  expect_lt(rel_err(low0$theta, mem$theta), 1e-10)
  expect_lt(rel_err(low0$coefficients, mem$coefficients), 1e-10)
  chr <- which(vapply(mem$LOCOResult, function(r) isTRUE(r$isLOCO), logical(1)))
  expect_identical(chr, 1:2)
  for (j in chr) {
    r_mem <- mem$LOCOResult[[j]]
    r_low <- ll_load(file.path(d, sprintf("low_chr%d.rda", j)))$LOCOResult[[j]]
    expect_true(isTRUE(r_low$isLOCO), info = j)
    expect_lt(rel_err(r_low$coefficients, r_mem$coefficients), chr_tol[j])
    expect_lt(rel_err(r_low$fitted.values, r_mem$fitted.values), chr_tol[j])
    expect_lt(rel_err(r_low$cov, r_mem$cov), chr_tol[j])
    if (sa_agrees[j]) {
      expect_lt(rel_err(r_low$obj.noK$S_a, r_mem$obj.noK$S_a), chr_tol[j])
    } else {
      # the characterization: not merely loose, but a different stopping point
      expect_gt(rel_err(r_low$coefficients, r_mem$coefficients), 1e-7)
      expect_gt(rel_err(r_low$obj.noK$S_a, r_mem$obj.noK$S_a), 0.5)
    }
  }
}

test_that("binary: the low-memory per-chromosome files match the in-memory LOCO results", {
  skip_unless_level("standard")
  plink <- file.path(tempfile("twochr_"), "sub")
  dir.create(dirname(plink), recursive = TRUE)
  sm_two_chr_fixture(plink)
  # chromosome 2's refit used to stop at a different point in the two layouts
  # (the in-memory path warm-starts from chromosome 1's alpha, the low-memory
  # path from the full model's, and at tol.coef = 0.02 that was the answer);
  # with the inner pass at tolCoef = 1e-5 (LEDGER #125) both reach the same
  # point: coefficients identical, S_a within 5e-4 (measured 2026-09-10)
  ll_check_trait("binary", plink, chr_tol = c(1e-10, 1e-3), sa_agrees = c(TRUE, TRUE))
})

test_that("quantitative: the low-memory per-chromosome files match the in-memory LOCO results", {
  skip_unless_level("standard")
  plink <- file.path(tempfile("twochr_"), "sub")
  dir.create(dirname(plink), recursive = TRUE)
  sm_two_chr_fixture(plink)
  ll_check_trait("quantitative", plink, chr_tol = c(1e-10, 1e-10), sa_agrees = c(TRUE, TRUE))
})
