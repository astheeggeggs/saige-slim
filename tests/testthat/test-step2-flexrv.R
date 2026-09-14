# FlexRV (F3/F4): the transform sets end to end through the step-2 CLI, on
# extdata/input/group_new_snpid_flexrv.txt (fixture facts in helper-step2-flexrv.R).
#
# Every cross-run comparison here is NUMERIC at ~1e-10, never byte, because
# each pairs two different code paths (internal transform sets against the
# external weight-line route, a transform set against a legacy single-set
# run). The one byte-identity claim -- --nThreads and --weights.beta must not
# move a FlexRV file -- is about ONE path re-run with an inert knob.
#
# Tier: standard.

fx_names <- function(scoreSpec = SAIGE::flexrvDefaultScoreSpec,
                     mafSpec = SAIGE::flexrvDefaultMAFSpec) {
  # N = 1000: the sets a run composes at the models' in-sample count, i.e.
  # without the MAC rungs at or above 2*N*maxMAF = 20 (LEDGER #142)
  SAIGE::flexrvTransformTable(scoreSpec, mafSpec, "REVEL", maxMAF = 0.01, N = 1000)$names
}
test_that("the default run writes 112 sets per gene, in the normative order, plus one Cauchy row", {
  r <- fx_burden_bin()
  fx_ok(r)
  d <- r$region
  expect_identical(colnames(d)[1:5], c("Region", "Group", "max_MAF", "Weight", "Pvalue_Burden"))
  rows <- fx_rows(r)
  for (g in c("GENE1", "GENE2")) {
    w <- rows$Weight[rows$Region == g]
    expect_identical(length(w), 112L, label = g)
    expect_identical(w, fx_names(), label = g)
    expect_true(all(rows$Group[rows$Region == g] == "lof;missense"))
    expect_true(all(rows$max_MAF[rows$Region == g] == "0.01"))
  }
  expect_identical(nrow(fx_cauchy(r)), 2L)
  expect_true(fx_log_has(r, paste0("FlexRV: 112 transform sets (16 score x 7 MAF) from score ",
                                   "'REVEL' at max MAF 0.01")))
  # the memory line, printed once
  expect_identical(sum(grepl("dense per-cell working memory is at least", r$log, fixed = TRUE)),
                   1L)
  expect_true(fx_log_has(r, paste0("112 transform sets over 1 (annotation x max_MAF) cell(s) = ",
                                   "112 weighted cells")))
  # lever 5 (2026-09-12): identical BURDEN cells are evaluated once. The
  # tally is printed once and adds up. served is a measurement: since LEDGER
  # #142 dropped the rungs that duplicated unw, no twin is guaranteed on this
  # fixture (one is served on the box that made the goldens), so only the
  # accounting is asserted.
  tally <- grep("^BURDEN cells: ", r$log, value = TRUE)
  expect_identical(length(tally), 1L)
  nums <- as.integer(regmatches(tally, gregexpr("[0-9]+", tally))[[1]])
  expect_identical(nums[1], 224L)
  expect_identical(nums[2] + nums[3], nums[1])
  # the twin bound: the one printed-equal pair left (GENE2 lof_MAC5 / lof_MAC10)
  expect_lte(nums[3], 1L)
})

test_that("the pooled Cauchy row is the CCT over the gene's set rows (p_FlexRV)", {
  r <- fx_burden_bin()
  rows <- fx_rows(r)
  cau <- fx_cauchy(r)
  for (g in c("GENE1", "GENE2")) {
    p <- as.numeric(rows$Pvalue_Burden[rows$Region == g])
    expect_true(all(is.finite(p)))
    # from the PRINTED p-values, so the agreement is bounded by their 7 digits
    expect_equal(as.numeric(cau$Pvalue_Burden[cau$Region == g]), get_CCT_pvalue(p),
                 tolerance = 1e-5, label = g)
  }
})

test_that("a MAC threshold above every MAC is never composed; the ladder below it is live", {
  # 2N = 2,000 and max MAF 0.01: no admitted marker has MAC > 20, so
  # MAC25..MAC500 would select everything, exactly as unw does. Until LEDGER
  # #142 (2026-09-13) those rows were emitted, byte-identical to their score
  # transform's unw row, and pooled; now they are dropped before composition
  # (test-step2-flexrv-rungs.R has the count and the warning).
  r <- fx_burden_bin()
  rows <- fx_rows(r)
  strip <- function(d) {
    d$Weight <- NULL
    rownames(d) <- NULL
    d
  }
  expect_false(any(grepl("_MAC(25|50|100|200|500)$", rows$Weight)))
  for (sc in c("unw", "lof", "cbrt0.5", "sbeta20")) {
    base <- strip(rows[rows$Weight == paste0("REVEL_", sc, "_unw"), ])
    expect_identical(nrow(base), 2L, label = sc)
    # the ladder is NOT trivially degenerate: MAC10 (below the bound) differs from unw
    expect_false(identical(strip(rows[rows$Weight == paste0("REVEL_", sc, "_MAC10"), ]), base))
  }
})

test_that("internal expansion == the external weight-line route (score transforms, MAF unw)", {
  # The reconciliation control for the whole machinery: the 16 score
  # transforms written out as 16 `weight:NAME` lines at full precision, run as
  # ordinary custom sets, must give the transform sets' rows. MAF unw keeps
  # the weights free of in-sample quantities, which a file cannot carry.
  # (test-step2-flexrv-reconcile.R does the same for all 192 sets, with the
  # in-sample MAC from two sources; this one needs no MAC at all)
  fx_skip("standard")
  gf <- fx_external_groupfile(numeric(0), mafSpec = "unw")
  ext <- fx_run(c(FX_G, FX_BIN, FX_SG, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                  paste0("--groupFile=", gf$path), "--annotation_in_groupTest=lof:missense",
                  "--maxMAF_in_groupTest=0.01", FX_BURDEN), "external")
  fx_ok(ext)
  int <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                  "--flexRV_MAFTransforms=unw"), "internal-unw")
  fx_ok(int)
  a <- fx_rows(int)
  b <- fx_rows(ext)
  expect_identical(nrow(a), 32L)
  expect_identical(a$Weight, rep(gf$names, 2))
  expect_identical(a$Weight, b$Weight)
  fx_expect_rows_equal(a, b, label = "internal vs external")
  # the allele accounting is weight-independent and must match exactly
  fx_expect_accounting_identical(a, b)
  expect_equal(as.numeric(fx_cauchy(int)$Pvalue_Burden), as.numeric(fx_cauchy(ext)$Pvalue_Burden),
               tolerance = 1e-10)
  # positive control on the comparison: the internal run's OTHER MAF
  # transforms are not the external rows
  full <- fx_rows(fx_burden_bin())
  expect_false(isTRUE(all.equal(as.numeric(full$Pvalue_Burden[full$Weight == "REVEL_sbeta5_MAC5"]),
                                as.numeric(b$Pvalue_Burden[b$Weight == "REVEL_sbeta5_unw"]),
                                tolerance = 1e-10)))
})

test_that("unw_unw is the legacy all-ones run, and lof_unw the LoF-only run (quantitative)", {
  # no ultra-rare collapse on either side, so the comparison is between two
  # burdens over the same real markers; quantitative, so no ER path can gate
  # the two sets differently
  int <- fx_burden_quant()
  fx_ok(int)
  rows <- fx_rows(int)
  plain <- fx_run(c(FX_G, fx_quant(), FX_COMMON, "--annotation_in_groupTest=lof:missense",
                    "--maxMAF_in_groupTest=0.01", FX_BURDEN, "--is_no_weight_in_groupTest=TRUE",
                    "--MACCutoff_to_CollapseUltraRare=0"), "quant-noweight")
  fx_ok(plain)
  fx_expect_rows_equal(rows[rows$Weight == "REVEL_unw_unw", ], fx_rows(plain),
                       label = "unw_unw vs No_Weights")
  lofonly <- fx_run(c(FX_G, fx_quant(), FX_COMMON, "--annotation_in_groupTest=lof",
                      "--maxMAF_in_groupTest=0.01", FX_BURDEN, "--is_no_weight_in_groupTest=TRUE",
                      "--MACCutoff_to_CollapseUltraRare=0"), "quant-lofonly")
  fx_ok(lofonly)
  fx_expect_rows_equal(rows[rows$Weight == "REVEL_lof_unw", ], fx_rows(lofonly),
                       label = "lof_unw vs lof-only")
})

test_that("--weights.beta is inert on the transform sets", {
  a <- fx_burden_bin()
  b <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN, "--weights.beta=1,1"),
              "beta11")
  fx_ok(b)
  expect_identical(readLines(b$path), readLines(a$path))
})

test_that("a set with no qualifying marker is omitted, and the Cauchy pool is over the rest", {
  # --minMAC=1.5 removes every singleton at QC, so the MAC1 sets have an
  # all-zero weight vector: BURDEN omits their 16 rows per gene
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, "--chrom=1", "--minMAF=0", "--minMAC=1.5", "--LOCO=FALSE",
                "--groupFile=./input/group_new_snpid_flexrv.txt", FX_FLEX, FX_BURDEN), "minmac15")
  fx_ok(r)
  rows <- fx_rows(r)
  cau <- fx_cauchy(r)
  for (g in c("GENE1", "GENE2")) {
    w <- rows$Weight[rows$Region == g]
    expect_identical(length(w), 96L, label = g)
    expect_false(any(grepl("_MAC1$", w)))
    expect_identical(w, grep("_MAC1$", fx_names(), value = TRUE, invert = TRUE))
    p <- as.numeric(rows$Pvalue_Burden[rows$Region == g])
    expect_equal(as.numeric(cau$Pvalue_Burden[cau$Region == g]), get_CCT_pvalue(p),
                 tolerance = 1e-5, label = g)
  }
})

test_that("SKAT-O on a reduced spec: schema, rows, and unw_unw against the legacy all-ones run", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, "--is_fastTest=FALSE",
                "--flexRV_scoreTransforms=unw;sbeta:2", "--flexRV_MAFTransforms=unw;mac:5"),
              "skato")
  fx_ok(r)
  d <- r$region
  expect_identical(colnames(d)[1:5], c("Region", "Group", "max_MAF", "Weight", "Pvalue"))
  expect_true(all(c("Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(d)))
  rows <- fx_rows(r)
  expect_identical(rows$Weight, rep(fx_names("unw;sbeta:2", "unw;mac:5"), 2))
  expect_identical(nrow(fx_cauchy(r)), 2L)
  plain <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, "--annotation_in_groupTest=lof:missense",
                    "--maxMAF_in_groupTest=0.01", "--is_fastTest=FALSE",
                    "--is_no_weight_in_groupTest=TRUE"), "skato-noweight")
  fx_ok(plain)
  fx_expect_rows_equal(rows[rows$Weight == "REVEL_unw_unw", ], fx_rows(plain),
                       label = "SKAT-O unw_unw vs No_Weights")
  # the sbeta2 sets differ from unw (the score axis is live on this fixture)
  expect_false(identical(rows$Pvalue[rows$Weight == "REVEL_unw_unw"],
                         rows$Pvalue[rows$Weight == "REVEL_sbeta2_unw"]))
  # LEDGER #117 FIXED 2026-09-11; this block is the INVERTED characterization.
  # Every real marker in the cell has MAC > 10, so the MAC5 sets weight them
  # exactly 0. That used to leave their zero-variance rows in the weighted Phi,
  # and the binary SKAT-O kernel's adjustment computed Score^2 / 0 -- the cell's
  # three p-values came back NA through the #57 failure path and the gene's
  # Cauchy row silently combined only the survivors (#127 route (a)).
  # A zero weight now means the marker is not in the set, so what is left of a
  # MAC5 cell is its collapsed ultra-rare pseudo-marker, which a custom or
  # transform set weights 1: a one-marker cell, which SKAT-O reduces to burden.
  mac5 <- rows[grepl("_MAC5$", rows$Weight), ]
  expect_identical(nrow(mac5), 4L)
  for (cc in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")) {
    expect_true(all(mac5[[cc]] != "NA"), label = cc)
  }
  expect_true(all(mac5$Pvalue_ACATV != "NA"))
  # a one-marker SKAT-O cell IS its burden test, which is the check that what
  # survived is the pseudo-marker alone rather than something arbitrary
  expect_identical(mac5$Pvalue, mac5$Pvalue_Burden)
  # nothing failed, so no sidecar and the summary says zero
  expect_false(file.exists(paste0(r$path, ".skatFailures.txt")))
  expect_true(fx_log_has(r, "SKAT failures: 0 cells returned NA"))
})

test_that("--nThreads=4 is byte-identical to one thread on the default BURDEN run, and says so", {
  a <- fx_burden_bin()
  b <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN, "--nThreads=4"), "threads4")
  fx_ok(b)
  expect_true(fx_log_has(b, "Region marker loop will use 4 in-process threads"))
  expect_identical(readLines(b$path), readLines(a$path))
})
