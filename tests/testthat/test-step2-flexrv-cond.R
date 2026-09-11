# FlexRV with --condition (the conditional reconciliation the F2-F4 handoff
# left open; only the missing --weights_for_condition error was tested). The
# transform sets take the custom-set route through the conditional export
# (`w0 = 1` on the pseudo-marker, LEDGER #17), so the internal expansion must
# give the rows the same weights give as external `weight:NAME` lines under
# the same --condition/--weights_for_condition. Scaffolding in
# helper-step2-flexrv.R. Tier: standard.

FX_COND <- c("--condition=rs30,rs79", "--weights_for_condition=1,1")

test_that("--condition: internal transform sets == the external weight-line route, cond columns", {
  fx_skip("standard")
  skip_if_no_goldens()
  int <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN, "--flexRV_MAFTransforms=unw",
                  FX_COND), "cond-int")
  fx_ok(int)
  d <- int$region
  expect_identical(colnames(d)[1:10],
                   c("Region", "Group", "max_MAF", "Weight", "Pvalue_Burden", "BETA_Burden",
                     "SE_Burden", "Pvalue_Burden_c", "Beta_Burden_c", "seBeta_Burden_c"))
  gf <- fx_external_groupfile(fx_mac_from_golden(), mafSpec = "unw")
  ext <- fx_run(c(FX_G, FX_BIN, FX_SG, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                  paste0("--groupFile=", gf$path), "--annotation_in_groupTest=lof:missense",
                  "--maxMAF_in_groupTest=0.01", FX_BURDEN, FX_COND), "cond-ext")
  fx_ok(ext)
  a <- fx_rows(int)
  b <- fx_rows(ext)
  expect_identical(nrow(a), 32L)
  expect_identical(a$Weight, rep(gf$names, 2))
  expect_identical(a$Weight, b$Weight)
  # unconditional AND conditional columns, at the printed precision
  fx_expect_rows_equal(a, b, tol = 1e-5, label = "cond internal vs external")
  fx_expect_accounting_identical(a, b)
  # the pooled row: every p column is the CCT over the gene's rows, the
  # conditional one included
  ca <- fx_cauchy(int)
  cb <- fx_cauchy(ext)
  for (g in c("GENE1", "GENE2")) {
    for (cc in c("Pvalue_Burden", "Pvalue_Burden_c")) {
      p <- as.numeric(a[[cc]][a$Region == g])
      expect_true(all(is.finite(p)), label = paste(g, cc))
      expect_equal(as.numeric(ca[[cc]][ca$Region == g]), get_CCT_pvalue(p), tolerance = 1e-5,
                   label = paste("internal", g, cc))
      expect_equal(as.numeric(cb[[cc]][cb$Region == g]), get_CCT_pvalue(p), tolerance = 1e-5,
                   label = paste("external", g, cc))
    }
  }
  # positive control: conditioning does something, and set-specifically
  expect_false(identical(a$Pvalue_Burden_c, a$Pvalue_Burden))
  expect_gt(length(unique(a$Pvalue_Burden_c[a$Region == "GENE1"])), 10L)
})

test_that("--condition leaves the unconditional columns exactly as the unconditional run wrote", {
  # one path, one knob: the conditional run's shared columns are the
  # unconditional 192-set run's rows for the same 16 sets (P2-vs-P3 shape)
  fx_skip("standard")
  int <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN, "--flexRV_MAFTransforms=unw",
                  FX_COND), "cond-int")
  base <- fx_burden_bin()
  fx_ok(int)
  fx_ok(base)
  a <- fx_rows(int)
  b <- fx_rows(base)
  b <- b[b$Weight %in% a$Weight, ]
  rownames(a) <- NULL
  rownames(b) <- NULL
  shared <- intersect(colnames(a), colnames(b))
  expect_identical(a[shared], b[shared])
})

test_that("--weights.beta is inert on the conditional FlexRV run too", {
  # the pair's last legacy role is the conditioning-weight DEFAULT, which
  # --weights_for_condition (mandatory here) replaces; nothing else on the
  # conditional path may read it
  fx_skip("standard")
  int <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN, "--flexRV_MAFTransforms=unw",
                  FX_COND), "cond-int")
  b11 <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN, "--flexRV_MAFTransforms=unw",
                  FX_COND, "--weights.beta=1,1"), "cond-int-beta11")
  fx_ok(int)
  fx_ok(b11)
  expect_identical(readLines(b11$path), readLines(int$path))
})
