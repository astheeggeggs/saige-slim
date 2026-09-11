# FlexRV with a spec that selects ONE transform set: K = 1, the legacy output
# schema, and the one place the engine's implicit beta(1,25) set could have
# reached a FlexRV run. It never does on the BURDEN path (C++ composes the
# implicit beta set only when NO other set exists), and on the SKAT-O path
# region_postprocess's K = 1 fall-through used to STOP on a transform set
# (LEDGER #119, fixed: the set takes the custom-set branch). These pin both:
# a single unw x unw set is the No_Weights run and NOT the beta(1,25) run; a
# single non-unit set is its own external K = 1 custom-weight run; and
# --weights.beta is inert on both. Scaffolding (fxk_flex, fxk_plain) in
# helper-step2-flexrv.R. Tier: standard.

FXK_SPEC <- c("--flexRV_scoreTransforms=unw", "--flexRV_MAFTransforms=unw")

test_that("a single unw x unw set is the No_Weights run and not the beta(1,25) run: BURDEN", {
  fx_skip("standard")
  one <- fxk_flex(c(FXK_SPEC, FX_BURDEN), "k1-burden")
  fx_ok(one)
  expect_false("Weight" %in% colnames(one$region))   # K = 1: the legacy schema
  expect_identical(nrow(fx_rows(one)), 2L)
  now <- fxk_plain(c("--is_no_weight_in_groupTest=TRUE", FX_BURDEN), "k1-noweight-burden")
  fx_ok(now)
  fx_expect_rows_equal(fx_rows(one), fx_rows(now), tol = 1e-5, label = "unw_unw vs No_Weights")
  fx_expect_accounting_identical(fx_rows(one), fx_rows(now))
  # positive control: the implicit beta(1,25) run is a DIFFERENT burden
  beta <- fxk_plain(FX_BURDEN, "k1-beta-burden")
  fx_ok(beta)
  expect_false(isTRUE(all.equal(as.numeric(fx_rows(one)$Pvalue_Burden),
                                as.numeric(fx_rows(beta)$Pvalue_Burden), tolerance = 1e-5)))
  # and --weights.beta cannot reach the set
  b11 <- fxk_flex(c(FXK_SPEC, FX_BURDEN, "--weights.beta=1,1"), "k1-burden-beta11")
  fx_ok(b11)
  expect_identical(readLines(b11$path), readLines(one$path))
})

test_that("a single unw x unw set is the No_Weights run and not the beta(1,25) run: SKAT-O", {
  # LEDGER #119: before the fix this run STOPPED on the implicit-beta arm's
  # guard, which is what kept the beta set from being applied silently
  fx_skip("standard")
  one <- fxk_flex(c(FXK_SPEC, "--is_fastTest=FALSE"), "k1-skato")
  fx_ok(one)
  expect_false(fx_log_has(one, "identical(ws_kind"))
  expect_false("Weight" %in% colnames(one$region))
  expect_identical(colnames(one$region)[1:6],
                   c("Region", "Group", "max_MAF", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT"))
  now <- fxk_plain(c("--is_no_weight_in_groupTest=TRUE", "--is_fastTest=FALSE"),
                   "k1-noweight-skato")
  fx_ok(now)
  fx_expect_rows_equal(fx_rows(one), fx_rows(now), tol = 1e-5,
                       label = "SKAT-O unw_unw vs No_Weights")
  fx_expect_accounting_identical(fx_rows(one), fx_rows(now))
  beta <- fxk_plain("--is_fastTest=FALSE", "k1-beta-skato")
  fx_ok(beta)
  for (cc in c("Pvalue_Burden", "Pvalue_SKAT")) {
    expect_false(isTRUE(all.equal(as.numeric(fx_rows(one)[[cc]]),
                                  as.numeric(fx_rows(beta)[[cc]]), tolerance = 1e-5)), label = cc)
  }
  b11 <- fxk_flex(c(FXK_SPEC, "--is_fastTest=FALSE", "--weights.beta=1,1"), "k1-skato-beta11")
  fx_ok(b11)
  expect_identical(readLines(b11$path), readLines(one$path))
})

test_that("a single non-unit set (unw x mafbeta2) is its own external K = 1 custom run: SKAT-O", {
  fx_skip("standard")
  skip_if_no_goldens()
  one <- fxk_flex(c("--flexRV_scoreTransforms=unw", "--flexRV_MAFTransforms=mafbeta:2",
                    "--is_fastTest=FALSE"), "k1-mafbeta-skato")
  fx_ok(one)
  expect_false("Weight" %in% colnames(one$region))
  gf <- fx_external_groupfile(fx_mac_from_golden(), scoreSpec = "unw", mafSpec = "mafbeta:2")
  expect_identical(gf$names, "REVEL_unw_mafbeta2@0.01")
  ext <- fx_run(c(FX_G, fx_quant(), "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                  paste0("--groupFile=", gf$path), "--annotation_in_groupTest=lof:missense",
                  "--maxMAF_in_groupTest=0.01", "--MACCutoff_to_CollapseUltraRare=0",
                  "--is_fastTest=FALSE"), "k1-mafbeta-ext")
  fx_ok(ext)
  expect_false("Weight" %in% colnames(ext$region))
  fx_expect_rows_equal(fx_rows(one), fx_rows(ext), tol = 1e-5, label = "K=1 mafbeta2 vs external")
  fx_expect_accounting_identical(fx_rows(one), fx_rows(ext))
  # the set is live: it is not the unweighted run
  now <- fxk_plain(c("--is_no_weight_in_groupTest=TRUE", "--is_fastTest=FALSE"),
                   "k1-noweight-skato")
  fx_ok(now)
  expect_false(isTRUE(all.equal(as.numeric(fx_rows(one)$Pvalue_SKAT),
                                as.numeric(fx_rows(now)$Pvalue_SKAT), tolerance = 1e-5)))
})
