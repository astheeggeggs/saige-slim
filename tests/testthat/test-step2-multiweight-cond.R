# Multi-weight (K > 1) region tests: CONDITIONAL ANALYSIS (C11, design D7).
#
# At K > 1 the C++ side exports the G1-raw G1tilde_P_G2tilde matrix plus
# VarInvMat_G2_cond_w0 / Tstat_G2_cond_w0 and skips its baked
# TstatAdjCond/VarMatAdjCond, and R applies set r's G1-side weight per cell;
# at q_weight == 1 the baked C++ path runs verbatim. These reconcile the two
# and pin the schema the conditional run emits.
#
# The reconciliation goes through the BETA route, not a custom-weighted K = 1
# run: the latter reads an uninitialised w0 on ultra-rare rows (LEDGER #17)
# and cannot serve as a baseline.
#
# Split out of test-step2-multiweight.R; shared scaffolding is in
# helper-multiweight.R. Tier: standard.

test_that("conditional analysis composes with K>1: schema, per-set adjustment, pooled Cauchy", {
  d <- cond_multi_run()$region
  # full schema: Weight + omnibus columns + the conditional columns, in order
  expect_identical(colnames(d)[1:9],
                   c("Region", "Group", "max_MAF", "Weight", "Pvalue",
                     "Pvalue_Burden", "Pvalue_SKAT", "Pvalue_ACATV", "Pvalue_ACATO"))
  expect_identical(colnames(d)[12:16],
                   c("Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond",
                     "BETA_Burden_cond", "SE_Burden_cond"))
  r <- d[d$Group != "Cauchy", ]
  # weight-innermost row order, exactly one pooled Cauchy row per region
  expect_identical(unique(r$Weight), WSETS)
  expect_identical(r$Weight, rep(WSETS, times = nrow(r) / 2))
  expect_identical(as.integer(table(d$Group[d$Group == "Cauchy"])), 2L)
  # positive control: the per-set adjustment is set-specific -- the two sets'
  # conditional p-values differ (their unconditional ones already do)
  cadd <- r[r$Weight == "CADD", ]
  revel <- r[r$Weight == "REVEL", ]
  expect_false(identical(cadd$Pvalue_cond, revel$Pvalue_cond))
  # conditioning must do something: cond and uncond differ somewhere per set
  expect_false(identical(cadd$Pvalue_cond, cadd$Pvalue))
  # the pooled Cauchy row combines BOTH sets' rows, conditional columns too
  cz <- d[d$Group == "Cauchy", ]
  for (gene in unique(r$Region)) {
    rows <- r[r$Region == gene, ]
    expect_equal(pnum(cz$Pvalue_cond[cz$Region == gene]),
                 acat_p(pnum(rows$Pvalue_cond)), tolerance = 1e-12)
    expect_equal(pnum(cz$Pvalue_Burden_cond[cz$Region == gene]),
                 acat_p(pnum(rows$Pvalue_Burden_cond)), tolerance = 1e-12)
  }
  # conditioning leaves the unconditional columns untouched: shared columns
  # equal the unconditional K=2 run's, column for column (same computation,
  # same platform, same run configuration -- the P2-vs-P3 golden invariant)
  u <- multi_run()$region
  shared <- intersect(colnames(u), colnames(d))
  expect_identical(d[shared], u[shared])
})

test_that("the K>1 conditional adjustment reconciles with the K=1 baked path (beta route)", {
  # Two-ways check per the header comment: the composed run's Beta_1_25 rows
  # against the plain conditional run (the implicit beta(1,25) K = 1 path,
  # C++-baked w0). Labels exact, statistics at 1e-12.
  a <- cond_beta2_run()$region
  b <- cond_plain_run()$region
  a25 <- a[a$Weight == "Beta_1_25" & a$Group != "Cauchy", ]
  bb <- b[b$Group != "Cauchy", ]
  expect_identical(a25$Region, bb$Region)
  expect_identical(a25$Group, bb$Group)
  expect_identical(a25$max_MAF, bb$max_MAF)
  expect_stats_close(a25, bb)
  # the second set is a genuinely different weighting, not a copy
  a11 <- a[a$Weight == "Beta_1_1" & a$Group != "Cauchy", ]
  expect_false(identical(a11$Pvalue_cond, a25$Pvalue_cond))
})

test_that("a forced Weight column composes with conditioning at K=1", {
  # is_output_weight_column at a single set: the C++ side is the literal
  # q_weight == 1 baked path (w0 baked into the exported matrix), and R's
  # multi-weight branch must NOT re-apply the set weight to it. For this
  # binary trait the baked matrix flows through the legacy scaled adjustment;
  # the sibling quantitative arm (baked TstatAdjCond/VarMatAdjCond exports) is
  # covered by test-unit-region-cond.R. Dropping the schema columns must
  # recover the plain conditional run.
  a <- cond_forced_run()$region
  b <- cond_plain_run()$region
  expect_identical(colnames(a)[4], "Weight")
  a1 <- a[a$Group != "Cauchy", ]
  bb <- b[b$Group != "Cauchy", ]
  expect_identical(unique(a1$Weight), "Beta")
  expect_identical(a1$Region, bb$Region)
  expect_identical(a1$Group, bb$Group)
  expect_identical(a1$max_MAF, bb$max_MAF)
  expect_stats_close(a1, bb)
})

test_that("BURDEN conditional at K>1 reconciles with a K=1 single-set run", {
  # The BURDEN conditional path is pure C++ (per-jmr genoSumVec through
  # getMarkerPval's internal conditional machinery; no G1tilde matrix, no w0,
  # so LEDGER #17 does not pollute this comparison and the CUSTOM fixture can
  # be reconciled directly): the K=2 run's CADD rows against a CADD-only K=1
  # file derived from the same fixture.
  skip_unless_level("standard")
  src <- file.path(saige_extdata, "input", "group_new_snpid_multiWeights.txt")
  caddfile <- tempfile(fileext = ".txt")
  writeLines(grep("weight,REVEL", readLines(src),
                  fixed = TRUE, invert = TRUE, value = TRUE), caddfile)
  run_burden_cond <- function(groupFile) {
    out <- file.path(tempfile("step2_mwbcond_"), "region.txt")
    dir.create(dirname(out), recursive = TRUE)
    args <- c(MW_ARGS("group_new_snpid_multiWeights.txt", out),
              "--r.corr=1", "--condition=rs30,rs79")
    # is_fastTest=TRUE + BURDEN drops all regions but the last (LEDGER #27);
    # run the correct non-fast path so both regions are compared
    args[grep("^--is_fastTest=", args)] <- "--is_fastTest=FALSE"
    args[grep("^--groupFile=", args)] <- paste0("--groupFile=", groupFile)
    expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
    read_saige_table(out)
  }
  a <- run_burden_cond(file.path(saige_extdata, "input",
                                 "group_new_snpid_multiWeights.txt"))
  b <- run_burden_cond(caddfile)
  # conditional BURDEN columns are present in the C++-written schema
  expect_true(all(c("Pvalue_Burden_c", "Beta_Burden_c", "seBeta_Burden_c")
                  %in% colnames(a)))
  aC <- a[a$Group != "Cauchy" & a$Weight == "CADD", ]
  bb <- b[b$Group != "Cauchy", ]
  expect_identical(aC$Region, bb$Region)
  expect_identical(aC$Group, bb$Group)
  expect_identical(aC$max_MAF, bb$max_MAF)
  for (cc in c("Pvalue_Burden", "BETA_Burden", "SE_Burden",
               "Pvalue_Burden_c", "Beta_Burden_c", "seBeta_Burden_c",
               "MAC", "Number_rare", "Number_ultra_rare")) {
    if (cc %in% colnames(a) && cc %in% colnames(b)) {
      expect_equal(pnum(aC[[cc]]), pnum(bb[[cc]]), tolerance = 1e-12, label = cc)
    }
  }
  # positive control: REVEL's conditional burden differs from CADD's
  aR <- a[a$Group != "Cauchy" & a$Weight == "REVEL", ]
  expect_false(identical(aR$Pvalue_Burden_c, aC$Pvalue_Burden_c))
})

test_that("a single-weight conditional run is still accepted", {
  # Baseline control: the plain conditional run without any multi-weight
  # schema keeps the legacy columns, so everything the tests above observe is
  # attributable to the weight configuration and not to --condition itself.
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_mwcond_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(MW_ARGS("group_new_snpid_withWeights.txt", out),
                            "--condition=rs30,rs79"))
  expect_cli_ok(status)
  res <- read_saige_table(out)
  expect_false("Weight" %in% colnames(res))
  # the SKAT-O region writer names the conditional columns with a _cond suffix
  expect_true(all(c("Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond")
                  %in% colnames(res)))
})
