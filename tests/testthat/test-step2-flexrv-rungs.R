# The LEDGER #142 characterization (decided 2026-09-13, the owner's option b).
# Region ingest admits a marker only if its folded MAF is at or below
# --flexRV_maxMAF, so every admitted MAC is at most 2*N*maxMAF; a MAC rung at
# or above that bound weights every marker 1 and is `unw` under another name.
# Until the fix those sets were composed, evaluated, emitted byte-identical to
# their score transform's unw row and pooled into p_FlexRV, which was therefore
# 1.45x (GENE1) and 0.966x (GENE2) the value over distinct sets on this fixture.
# They are now dropped before composition and named in a warning.
#
# The count the row specified: sets dropped = (rungs at or above 2*N*maxMAF) x
# (score transforms) = 5 x 16 = 80 on the bundled fixture (N = 1,000, max MAF
# 0.01, bound 20), and 0 once --flexRV_maxMAF lifts the bound above the top
# rung (0.26: bound 520 > 500). Both arms are red at unpatched HEAD, where the
# table builder has no N argument (and, behind that, 192 rows were emitted and
# no warning printed).
# Tier: standard.

fr_names <- function(N = NULL, maxMAF = 0.01) {
  SAIGE::flexrvTransformTable(SAIGE::flexrvDefaultScoreSpec, SAIGE::flexrvDefaultMAFSpec,
                              "REVEL", maxMAF = maxMAF, N = N)
}
fr_drop_line <- "dropped before the sets are composed (LEDGER #142)"

test_that("LEDGER #142: the rungs at or above 2*N*maxMAF are dropped, named, and absent", {
  fx_skip("standard")
  skip_if_no_goldens()
  full <- fr_names()
  kept <- fr_names(N = 1000)
  degenerate <- c(25, 50, 100, 200, 500)
  # the arithmetic the row specified, from the table builder alone
  expect_equal(kept$bound, 20)
  expect_identical(kept$dropped, paste0("MAC", degenerate))
  expect_identical(length(full$names) - length(kept$names), length(degenerate) * 16L)
  expect_identical(length(kept$names), 112L)
  # ... and from the run: 112 rows per gene, the dropped names absent, the
  # rung below the bound present, the warning naming N, the bound and the rungs
  r <- fx_burden_bin()
  fx_ok(r)
  rows <- fx_rows(r)
  for (g in c("GENE1", "GENE2")) {
    w <- rows$Weight[rows$Region == g]
    expect_identical(w, kept$names, label = g)
    expect_false(any(grepl("_MAC(25|50|100|200|500)$", w)), label = g)
    expect_identical(sum(grepl("_MAC10$", w)), 16L, label = g)
  }
  expect_true(fx_log_has(r, paste0("FlexRV: at N = 1000 and --flexRV_maxMAF=0.01 every admitted ",
                                   "marker has MAC <= 20, so the MAC rungs 25, 50, 100, 200, 500 ",
                                   "would weight every marker 1 and duplicate 'unw'")))
  expect_true(fx_log_has(r, paste0(fr_drop_line, ": 112 sets (16 score x 7 MAF), not 192")))
  # the pooled row is the CCT over the 112 surviving rows and nothing else
  cau <- fx_cauchy(r)
  for (g in c("GENE1", "GENE2")) {
    p <- as.numeric(rows$Pvalue_Burden[rows$Region == g])
    expect_identical(length(p), 112L, label = g)
    expect_equal(as.numeric(cau$Pvalue_Burden[cau$Region == g]), get_CCT_pvalue(p),
                 tolerance = 1e-5, label = g)
  }
})

test_that("LEDGER #142: with the bound lifted above the top rung nothing is dropped", {
  fx_skip("standard")
  skip_if_no_goldens()
  # 2 * 1000 * 0.26 = 520 > 500: every rung can select something
  kept <- fr_names(N = 1000, maxMAF = 0.26)
  expect_identical(kept$dropped, character(0))
  expect_identical(length(kept$names), 192L)
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, "--flexRV_score=REVEL", "--flexRV_maxMAF=0.26",
                "--annotation_in_groupTest=lof:missense", FX_BURDEN), "rungs-nodrop")
  fx_ok(r)
  rows <- fx_rows(r)
  for (g in c("GENE1", "GENE2")) {
    w <- rows$Weight[rows$Region == g]
    expect_identical(w, kept$names, label = g)
    expect_identical(sum(grepl("_MAC500$", w)), 16L, label = g)
  }
  expect_false(fx_log_has(r, fr_drop_line))
  expect_true(fx_log_has(r, "FlexRV: 192 transform sets (16 score x 12 MAF)"))
})
