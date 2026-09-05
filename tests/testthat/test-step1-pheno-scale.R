# LEDGER #65 FIXED: a quantitative null model is invariant to the unit of the
# phenotype. Three AI-REML rules in glmmkin.ai_PCG_Rcpp_Quantitative were
# absolute in phenotype units -- a stop() when var(Y) < 0.1, the `+ tol` in the
# convergence rule's denominator, and the divergence cap max(tau) > tol^-2 =
# 2500 -- and tau started at (1, 0) whatever the scale, so y in millimetres and
# y in metres took different paths: on this fixture y * 100 hit the cap
# ("model not converged") and y / 100 hit the stop(). The fit now runs on
# y / sd(Y) and scales back, so the three runs below must agree exactly up to
# the scale factor.
#
# Derived fixture: the committed phenotype file with y_quantitative multiplied
# by 100 and by 1/100, everything else untouched; the same .bed, sparse GRM and
# covariates as the quant_sparseGRM golden, WITHOUT invNormalize (which would
# undo the scaling and hide the bug). Each fit is ~5 s at N = 1,000.
#
# Expected relations, derived before running: theta scales by s^2,
# coefficients, linear predictors, fitted values and residuals by s, and the
# variance ratio by 1/s^2 -- its numerator variance is in 1/y^2 units while its
# denominator is in genotype units only (measured at y * 100 before the fix:
# ratios 1e-4 x the unscaled ones to 5e-8 relative, step-2 p-values identical
# to 1.5e-7; at y / 100 the two absolute 1e-4 floors on Sigma's diagonal made
# the ratio 33% wrong and step-2 p-values 57% wrong, which is why the variance
# ratios are estimated on the standardised trait too). The three fits are the
# SAME arithmetic on y/sd(y) up to the rounding of the scaling itself, so the
# tolerance is loose only against that rounding, 1e-8 relative. The end-to-end
# check is step 2: the single-variant p-values from the three models must agree.
#
# Tier: standard.

ps_pheno_file <- function(scale) {
  src <- file.path(saige_extdata, "input",
                   "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
  d <- read.table(src, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  d$y_quantitative <- d$y_quantitative * scale
  out <- tempfile(sprintf("pheno_x%g_", scale), fileext = ".txt")
  write.table(d, out, quote = FALSE, row.names = FALSE, sep = " ")
  out
}

ps_fit <- function(scale) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  prefix <- file.path(tempfile(sprintf("step1_scale%g_", scale)), "fit")
  dir.create(dirname(prefix), recursive = TRUE)
  st <- run_saige_cli("step1_fitNULLGLMM.R", c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed",
           ".sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    paste0("--phenoFile=", ps_pheno_file(scale)),
    "--phenoCol=y_quantitative", "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID", "--traitType=quantitative",
    "--IsOverwriteVarianceRatioFile=TRUE", "--isCateVarianceRatio=TRUE",
    paste0("--outputPrefix=", prefix)
  ))
  log <- readLines(attr(st, "log"), warn = FALSE)
  if (!identical(as.integer(st), 0L)) {
    return(list(status = as.integer(st), log = log))
  }
  e <- new.env()
  load(paste0(prefix, ".rda"), e)
  vr <- read.table(paste0(prefix, ".varianceRatio.txt"), header = FALSE, stringsAsFactors = FALSE)
  # step 2 on this model: golden P32's command (single-variant, quantitative)
  s2 <- paste0(prefix, "_step2.txt")
  st2 <- run_saige_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first", "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", prefix, ".rda"),
    paste0("--varianceRatioFile=", prefix, ".varianceRatio.txt"),
    "--is_output_moreDetails=TRUE",
    paste0("--SAIGEOutputFile=", s2)
  ))
  expect_cli_ok(st2)
  list(status = 0L, log = log, m = e$modglmm, vr = vr, s2 = read_saige_table(s2))
}

ps_cache <- new.env(parent = emptyenv())
ps_run <- function(scale) {
  key <- format(scale)
  if (is.null(ps_cache[[key]])) ps_cache[[key]] <- ps_fit(scale)
  ps_cache[[key]]
}

expect_scaled_fit <- function(base, other, s) {
  expect_identical(other$status, 0L)
  expect_true(isTRUE(other$m$converged))
  expect_lt(rel_err(other$m$theta, base$m$theta * s^2), 1e-8)
  expect_lt(rel_err(other$m$coefficients, base$m$coefficients * s), 1e-8)
  expect_lt(rel_err(other$m$linear.predictors, base$m$linear.predictors * s), 1e-8)
  expect_lt(rel_err(other$m$fitted.values, base$m$fitted.values * s), 1e-8)
  expect_lt(rel_err(other$m$residuals, base$m$residuals * s), 1e-8)
  # the variance ratios: same rows, same labels, values scaled by 1/s^2
  expect_identical(dim(other$vr), dim(base$vr))
  expect_identical(other$vr[[2]], base$vr[[2]])
  expect_identical(other$vr[[3]], base$vr[[3]])
  expect_lt(rel_err(other$vr[[1]] * s^2, base$vr[[1]]), 1e-8)
  # end to end: the same markers, the same p-values, BETA and SE in phenotype units
  expect_identical(other$s2$MarkerID, base$s2$MarkerID)
  expect_lt(rel_err(as.numeric(other$s2$p.value), as.numeric(base$s2$p.value)), 1e-6)
  expect_lt(rel_err(as.numeric(other$s2$BETA), as.numeric(base$s2$BETA) * s), 1e-6)
  expect_lt(rel_err(as.numeric(other$s2$SE), as.numeric(base$s2$SE) * s), 1e-6)
}

test_that("LEDGER #65 FIXED: multiplying the phenotype by 100 rescales the fit and nothing else", {
  base <- ps_run(1)
  expect_identical(base$status, 0L)
  expect_true(isTRUE(base$m$converged))
  # at unpatched HEAD this arm hit max(tau) > 2500 and reported 'not converged'
  expect_scaled_fit(base, ps_run(100), 100)
})

test_that("LEDGER #65 FIXED: dividing the phenotype by 100 rescales the fit and nothing else", {
  base <- ps_run(1)
  # at unpatched HEAD this arm was a stop(): 'variance of the phenotype is much
  # smaller than 1. Please consider invNormalize=T'
  expect_scaled_fit(base, ps_run(0.01), 0.01)
})

test_that("LEDGER #65: the unscaled fit is the committed golden within the harness tolerances", {
  # the standardisation must not have moved the answer beyond rounding on a
  # phenotype whose adjusted sd is already ~0.8 (the golden's own, minus its
  # invNormalize, which on this fixture is a near no-op: same theta to 6 digits)
  base <- ps_run(1)
  g <- new.env()
  load(file.path(saige_goldens, "step1", "quant_sparseGRM.rda"), g)
  expect_lt(rel_err(base$m$theta, g$modglmm$theta), 1e-2)
  expect_lt(rel_err(base$m$coefficients, g$modglmm$coefficients), 1e-3)
  # and the log says what it did
  expect_true(any(grepl("^Quantitative trait: the null model and the variance ratios are estimated",
                        base$log)))
})
