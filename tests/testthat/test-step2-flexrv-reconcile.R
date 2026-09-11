# FlexRV end to end, reconciled against SECOND computations of the same
# quantities (helper-step2-flexrv.R has the fixture arithmetic and the
# builders). Three second ways:
#
#  1. the whole 192-set table written out as 192 `weight:NAME` lines and run
#     as custom sets, twice: binary with the ultra-rare collapse on (MAC from
#     P1's golden) and quantitative with it off (MAC from the run's OWN
#     .singleAssoc.txt) -- two MAC sources, two collapse regimes, every MAF
#     transform live;
#  2. the MAC ladder as `--maxMAC_in_groupTest` cells of a plain unweighted
#     run, which is the same marker subset reached by a different flag;
#  3. a K = 1 custom 0/1 mask, which says whether what (2) finds is FlexRV's
#     or the engine's.
#
# (2) is where LEDGER #118 lived: a burden cell's variance-ratio category, ER
# eligibility and --minGroupMAC_in_BurdenTest gate used to read the CELL's
# unweighted MAC, shared by every set, so a set whose weights zeroed most of the
# cell was tested as if it were the whole cell. FIXED 2026-09-11 -- the count is
# per set now -- and the block that used to be marked CHARACTERIZATION is
# INVERTED below: the two parameterizations must agree.
#
# Cross-path comparisons are at 1e-5 (the printed precision); the accounting
# columns exactly. Tier: standard.

test_that("all 192 sets == 192 external weight lines: binary, collapse on, MAC from P1's golden", {
  fx_skip("standard")
  skip_if_no_goldens()
  int <- fx_burden_bin()
  fx_ok(int)
  ext <- fx_external_run(fx_external_groupfile(fx_mac_from_golden())$path, c(FX_BIN, FX_SG),
                         character(0), "ext192-bin")
  fx_ok(ext)
  fx_expect_int_eq_ext(int, ext, "binary collapse-on")
  # with the collapse on, the low MAC sets are carried by the pseudo-marker
  # (design doc, "What enters the test"): Number_ultra_rare > 0 on every row
  expect_true(all(as.numeric(fx_rows(int)$Number_ultra_rare) > 0))
})

test_that("all 192 sets == 192 external lines: quantitative, collapse off, MAC from the run", {
  fx_skip("standard")
  skip_if_no_goldens()
  q <- fx_quant_single()
  fx_ok(q)
  mac <- fx_mac_from_single(q)
  # the fixture arithmetic: 21 + 18 polymorphic markers at MAC <= 20 in the cell
  expect_identical(length(mac), 39L)
  expect_true(all(mac >= 1 & mac <= 20))
  # ... and the two MAC sources agree marker for marker (a control on both)
  expect_identical(mac, fx_mac_from_golden()[names(mac)])
  # --is_single_in_groupTest does not move the region file (one path, one knob)
  expect_identical(readLines(q$path), readLines(fx_burden_quant()$path))
  ext <- fx_external_run(fx_external_groupfile(mac)$path, fx_quant(),
                         "--MACCutoff_to_CollapseUltraRare=0", "ext192-quant")
  fx_ok(ext)
  fx_expect_int_eq_ext(q, ext, "quantitative collapse-off")
  expect_true(all(fx_rows(q)$Number_ultra_rare == "0"))
})

# ---------------------------------------------------------------------------
# The MAC ladder as --maxMAC_in_groupTest cells. --maxMAC_in_groupTest=t is
# converted to a max MAF of t / (2N) and unioned with --maxMAF_in_groupTest,
# so the plain unweighted run tests the SAME marker subset the MAC<=t set
# weights 1 -- one statistic, two flags. With the collapse off and one `null`
# variance ratio per MAC band, the two agree exactly wherever both burdens
# fall in the same band; where they do not, LEDGER #118.
# ---------------------------------------------------------------------------

test_that("the MAC ladder == --maxMAC_in_groupTest cells where both burdens share a VR band", {
  fx_skip("standard")
  skip_if_no_goldens()
  q <- fx_quant_single()
  fx_ok(q)
  plain <- fx_maxmac_run()
  fx_ok(plain)
  expect_true(fx_log_has(plain, "maxMAC_in_groupTest:  1 3 5 10  is specified"))
  mac <- fx_mac_from_single(q)
  gene_of <- function(ids) ifelse(as.integer(sub("^rs", "", ids)) <= 50, "GENE1", "GENE2")
  # the cells the plain run wrote: 1, 3, 5, 10 and the 0.01 cell. The MAC<=1
  # cell is the singletons alone -- 4 in GENE1, 3 in GENE2 -- and it IS written
  # now: LEDGER #126 moved --minGroupMAC_in_BurdenTest's default from 5 to 1,
  # where it is inert by construction. It used to be refused, which is what the
  # old characterization below turned on.
  expect_equal(sort(unique(as.numeric(fx_rows(plain)$max_MAF))),
               c(1 / 2000, 3 / 2000, 5 / 2000, 10 / 2000, 0.01), tolerance = 1e-12)
  for (g in c("GENE1", "GENE2")) {
    expect_identical(sum(mac[gene_of(names(mac)) == g & mac <= 1]), if (g == "GENE1") 4 else 3)
  }
  # LEDGER #118 FIXED: the two parameterizations agree at EVERY rung, not only
  # where both burdens happened to land in the same variance-ratio band. The
  # rungs 5 and 10 were green before the fix too (they were the control); 1 and
  # 3 are the ones that moved, and they are asserted in the block below.
  for (g in c("GENE1", "GENE2")) {
    for (t in c(5, 10)) {
      a <- fx_ladder_row(q, g, t)
      b <- fx_cell_row(plain, g, t)
      expect_identical(nrow(b), 1L, label = paste(g, t))
      fx_expect_rows_equal(a, b, tol = 1e-5, label = paste(g, "MAC", t))
      expect_identical(as.numeric(b$MAC), sum(mac[gene_of(names(mac)) == g & mac <= t]),
                       label = paste(g, t))
      expect_identical(as.integer(b$Number_rare), sum(gene_of(names(mac)) == g & mac <= t),
                       label = paste(g, t))
    }
  }
  a <- fx_ladder_row(q, "GENE2", 3)
  b <- fx_cell_row(plain, "GENE2", 3)
  expect_identical(as.numeric(b$MAC), 21)   # > 20.5 on both sides
  fx_expect_rows_equal(a, b, tol = 1e-5, label = "GENE2 MAC 3")
})

test_that("LEDGER #118 FIXED: a set's VR band, MAC gate and MAC column are the set's own", {
  # GENE1's MAC<=3 burden is 10 markers with MAC 19. Before the fix it was
  # tested with the CELL's unweighted MAC, 121 (band 2) where the plain 0.0015
  # cell selects band 1, and printed MAC 121 and Number_rare 21; the two rows
  # differed by exactly the ratio of the two `null` variance ratios in BETA and
  # in SE^2. They are the same row now. (Inverted from the CHARACTERIZATION
  # block in the same commit as the fix; the ratio identity below is kept as the
  # measurement of what the bug was worth, now asserted to be ABSENT.)
  fx_skip("standard")
  skip_if_no_goldens()
  q <- fx_quant_single()
  plain <- fx_maxmac_run()
  a <- fx_ladder_row(q, "GENE1", 3)
  b <- fx_cell_row(plain, "GENE1", 3)
  expect_identical(as.numeric(b$MAC), 19)
  expect_identical(as.numeric(b$Number_rare), 10)
  expect_identical(a$MAC, "19")
  expect_identical(a$Number_rare, "10")
  fx_expect_rows_equal(a, b, tol = 1e-5, label = "GENE1 MAC 3")
  # the bands really are distinguishable on this model, so the equality above
  # is a statement about the code and not about a degenerate variance-ratio
  # file (the control the old characterization needed, kept)
  vr <- fx_quant_vr()
  expect_identical(names(vr), c("1", "2"))
  expect_gt(abs(vr[["1"]] / vr[["2"]] - 1), 1e-3)

  # the MAC<=1 set is 4 singletons, MAC 4. At the default gate (1 since LEDGER
  # #126) both parameterizations write it; before, the FlexRV set was written
  # because the gate saw 121 while the plain cell was refused at 5.
  a1 <- fx_ladder_row(q, "GENE1", 1)
  b1 <- fx_cell_row(plain, "GENE1", 1)
  expect_identical(nrow(a1), 1L)
  expect_identical(nrow(b1), 1L)
  expect_identical(a1$MAC, "4")
  expect_identical(as.numeric(b1$MAC), 4)
  fx_expect_rows_equal(a1, b1, tol = 1e-5, label = "GENE1 MAC 1")

  # ... and it is the ENGINE's semantics, not FlexRV's, in both directions: a
  # K = 1 custom 0/1 mask over the same markers gives the same row as well.
  mask <- fx_weight_groupfile(fx_mac_from_single(q), "MASK3",
                              function(s, isLof, maf, m) matrix(as.numeric(m <= 3), ncol = 1))
  k1 <- fx_external_run(mask, fx_quant(), "--MACCutoff_to_CollapseUltraRare=0", "mask3")
  fx_ok(k1)
  m1 <- fx_rows(k1)
  expect_identical(nrow(m1), 2L)
  expect_false("Weight" %in% colnames(m1))   # K = 1, legacy schema
  m1g <- m1[m1$Region == "GENE1", ]
  expect_identical(m1g$MAC, "19")
  fx_expect_rows_equal(a, m1g, tol = 1e-5, label = "K=1 mask vs FlexRV MAC3")
  fx_expect_rows_equal(b, m1g, tol = 1e-5, label = "K=1 mask vs plain cell")
})
