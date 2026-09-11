# FlexRV under --GMMATmodel_varianceRatio_multiTraits_File: several null
# models per region in ONE genotype pass, the transform sets expanded once
# and shared. The claim is the region multi-model contract's
# (test-step2-region-multimodel.R): EACH model's region file and single-in-
# group file are byte-identical to what that model's own single-model run
# writes. One code path checked against itself with one model, not two paths
# compared. Scaffolding (FXM_MODELS, fxm_check) in helper-step2-flexrv.R.
# Tier: standard.

test_that("two models, 192-set BURDEN with single-in-group: each model's files reproduce its run", {
  fx_skip("standard")
  skip_if_no_goldens()
  man <- fxm_check(c(FX_FLEX, FX_BURDEN, "--is_single_in_groupTest=TRUE"), "burden")
  # the expansion (the R line) and the memory floor (the C++ line) are each
  # printed once for the pass, not once per model
  expect_identical(sum(grepl("FlexRV: 192 transform sets (16 score x 12 MAF)", man$log,
                             fixed = TRUE)), 1L)
  expect_identical(sum(grepl("dense per-cell working memory is at least", man$log,
                             fixed = TRUE)), 1L)
  for (nm in names(FXM_MODELS)) {
    d <- fx_read(paste0(man$path, "_", nm))
    expect_identical(sum(d$Group != "Cauchy"), 384L, label = nm)
    expect_identical(sum(d$Group == "Cauchy"), 2L, label = nm)
  }
})

test_that("two models, reduced-spec SKAT-O: each model's files reproduce its own run", {
  fx_skip("standard")
  skip_if_no_goldens()
  man <- fxm_check(c("--flexRV_score=REVEL", "--flexRV_maxMAF=0.02",
                     "--annotation_in_groupTest=lof:missense",
                     "--flexRV_scoreTransforms=unw;sbeta:2", "--flexRV_MAFTransforms=unw;mafbeta:2",
                     "--is_fastTest=FALSE"), "skato")
  l <- grep("^SKAT-O cells: ", man$log, value = TRUE)
  expect_length(l, 1L)
  n <- as.integer(regmatches(l, gregexpr("[0-9]+", l))[[1]])
  expect_identical(n[1], 16L)   # 2 models x 2 genes x 4 sets
  expect_identical(n[2] + n[3], n[1])
})
