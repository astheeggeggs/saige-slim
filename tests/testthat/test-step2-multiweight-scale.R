# LEDGER #134: the ultra-rare pseudo-marker's variance-ratio band and sparse-GRM
# choice used to be decided by the WEIGHTED collapse sum, whose units are
# weight x genotype, so the SCALE of a group file's weights decided which
# statistical treatment a set of alleles received. They now read MAC_nw, the
# UNWEIGHTED collapse's allele count -- already computed in the same scope and
# already what the reported group MAC uses (LEDGER #118).
#
# THE ASSERTION, and why it is not vacuous.
# A uniform rescale of every weight in a group file is a no-op for a burden
# p-value: it scales BETA and SE together and cancels. It is byte-exactly a
# no-op for Pvalue, Pvalue_Burden and Pvalue_SKAT -- that is pinned in
# test-step2-multiweight-behaviour.R for a factor of 2. It must therefore be a
# no-op for EVERY factor, and before this fix it was not.
#
# The factor here is 0.25 and it is chosen, not arbitrary. With the CLI defaults
# --cateVarRatioMinMACVecExclude=10,20.5 and --cateVarRatioMaxMACVecInclude=20.5,
# selectVarianceRatio (src/SAIGE_test.cpp) maps MAC <= 10 AND 10 < MAC <= 20.5
# both to varRatio(0), and only MAC > 20.5 to varRatio(1). So the ONLY boundary
# that changes the selected value is 20.5, and a factor must carry a cell across
# THAT to be detectable. Read off the committed golden
# harness/goldens/step2/plink_group_multiweight.txt.singleAssoc.txt, the six
# CADD UR-only cells have weighted sums 38.4, 43.5, 43.5, 22.1, 65.4, 65.4 --
# all above 20.5 -- and x0.25 takes ALL SIX below it, while the six REVEL cells
# (3.34 .. 9.51) stay below throughout. A factor of 2 would NOT do: it crosses
# the 10 boundary for two REVEL cells and nothing else, and both sides of 10
# select the same ratio, which is exactly why the doubling measurement recorded
# on the row moved no p-value.
#
# 0.25 is also a NEGATIVE POWER OF TWO, and that is load-bearing for the
# expect_identical assertions below: scaling by 2^-2 is exact in binary floating
# point, so a weight that was exactly representable still is, and any difference
# in the printed p-values is the variance-ratio band and not a 1-ULP rounding
# lottery. Do not "simplify" it to 0.2 or 0.3.
#
# Tier: standard. FIVE step-2 CLI runs: mws_run_scaled has NO cache (it is not
# mw_run) and each test_that block calls it afresh.

# Rescale every weight on every `weight:NAME` / `weight,NAME` line of a group
# file by a constant, leaving var / anno lines untouched.
mws_rescale_group_file <- function(src, factor) {
  lines <- readLines(file.path(saige_extdata, "input", src))
  out <- tempfile(fileext = ".txt")
  scaled <- vapply(lines, function(l) {
    tok <- strsplit(l, "[[:space:]]+")[[1]]
    if (length(tok) < 3 || !grepl("^weight", tok[2])) return(l)
    v <- suppressWarnings(as.numeric(tok[-(1:2)]))
    if (anyNA(v)) return(l)
    paste(c(tok[1], tok[2], format(v * factor, trim = TRUE, scientific = FALSE)),
          collapse = " ")
  }, character(1), USE.NAMES = FALSE)
  writeLines(scaled, out)
  out
}

mws_run_scaled <- function(factor) {
  skip_unless_level("standard")
  gf <- mws_rescale_group_file("group_new_snpid_multiWeights.txt", factor)
  out <- file.path(tempfile("step2_mws_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  # MW_ARGS pastes ./input/ in front of the group file, so pass an absolute path
  # by overriding the argument it builds rather than reusing it verbatim.
  args <- MW_ARGS("group_new_snpid_multiWeights.txt", out)
  args[grepl("^--groupFile=", args)] <- paste0("--groupFile=", gf)
  status <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(status)
  sa <- paste0(out, ".singleAssoc.txt")
  list(region = read_saige_table(out),
       single = if (file.exists(sa)) read_saige_table(sa) else NULL)
}

# The pseudo-marker's WEIGHTED collapse sum, per cell: the UR rows of the
# single-assoc file carry it as AC_Allele2, keyed by the cell label in MarkerID.
mws_ur_weighted <- function(run) {
  sg <- run$single
  ur <- sg[sg$CHR == "UR", ]
  stats::setNames(pnum(ur$AC_Allele2), ur$MarkerID)
}

test_that("a uniform rescale of a group file's weights moves no p-value", {
  skip_unless_level("standard")
  base <- mws_run_scaled(1)$region
  quarter <- mws_run_scaled(0.25)$region

  # the two runs must describe the same grid, or the comparison is meaningless
  expect_identical(base$Region, quarter$Region)
  expect_identical(base$Group, quarter$Group)
  expect_identical(base$max_MAF, quarter$max_MAF)
  expect_identical(base$Weight, quarter$Weight)

  # THE ASSERTION. The three p-value columns that a uniform rescale cannot
  # legitimately touch are byte-identical. At unpatched HEAD the six CADD UR-only
  # cells cross 20.5, pick a different variance ratio, and these differ.
  #
  # Pvalue_ACATV and Pvalue_ACATO are p-value columns in this schema and are
  # DELIBERATELY EXCLUDED: they are scale-dependent through a different
  # mechanism, the hardcoded 1 the pseudo-marker gets in the ACAT-V combination
  # vector, which is LEDGER #132 and is still OPEN. Asserting them here would
  # fail for a reason this row does not fix. When #132 lands, add them.
  for (col in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT")) {
    expect_identical(base[[col]], quarter[[col]],
                     info = paste("column", col, "is not rescale-invariant"))
  }

  # the allele counts never depended on the weights and must not start to
  for (col in c("MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")) {
    expect_identical(base[[col]], quarter[[col]],
                     info = paste("column", col, "moved under a pure rescale"))
  }
})

test_that("the rescale really does straddle the live variance-ratio boundary", {
  skip_unless_level("standard")
  # ANTI-VACUITY, and it is the assertion that keeps this file honest. The test
  # above passes trivially if the rescale never carries a cell across the only
  # boundary that changes the selected variance ratio -- which would happen
  # silently if the fixture's weights, the model's variance-ratio bins, or the
  # --cateVarRatioMinMACVecExclude default ever changed. So assert the regime.
  #
  # The live boundary is 20.5: MAC <= 10 and 10 < MAC <= 20.5 both select
  # varRatio(0), only > 20.5 differs (src/SAIGE_test.cpp). The quantity compared
  # against it before the fix is the pseudo-marker's WEIGHTED sum, which is its
  # AC_Allele2 on the UR rows of the single-assoc file.
  LIVE_EDGE <- 20.5
  b <- mws_run_scaled(1)
  q <- mws_run_scaled(0.25)
  wb <- mws_ur_weighted(b)
  wq <- mws_ur_weighted(q)
  expect_gt(length(wb), 0)
  expect_identical(names(wb), names(wq))

  # the weighted sums really are quartered -- proves the rescaled file was used
  expect_equal(unname(wq), unname(wb) / 4, tolerance = 1e-6)

  # THE STRADDLE: at least one cell above the edge at scale 1 and below it at
  # 0.25. Without this, the p-value test above could pass vacuously.
  straddlers <- names(wb)[wb > LIVE_EDGE & wq <= LIVE_EDGE]
  expect_gt(length(straddlers), 0)

  # and the UNWEIGHTED MAC -- what the fixed code reads -- is scale-invariant,
  # which is why the p-values no longer move
  ur_b <- b$region[b$region$Number_rare == "0" & b$region$Number_ultra_rare != "0", ]
  ur_q <- q$region[q$region$Number_rare == "0" & q$region$Number_ultra_rare != "0", ]
  expect_gt(nrow(ur_b), 0)
  expect_identical(ur_b$MAC, ur_q$MAC)
  expect_true(any(grepl("CADD", ur_b$Weight)))
  expect_true(any(grepl("REVEL", ur_b$Weight)))
})

test_that("BETA and SE do scale with the weights, so the rescale was real", {
  skip_unless_level("standard")
  # The positive control for the assertion above: a rescale that changed NOTHING
  # at all would also make the p-value test pass. BETA_Burden is in units of the
  # burden, so quartering every weight must quadruple it -- that proves the
  # rescaled file really was used.
  base <- mws_run_scaled(1)$region
  quarter <- mws_run_scaled(0.25)$region
  keep <- base$Group != "Cauchy" & base$BETA_Burden != "NA" & quarter$BETA_Burden != "NA"
  b <- pnum(base$BETA_Burden[keep])
  q <- pnum(quarter$BETA_Burden[keep])
  expect_gt(length(b), 0)
  expect_equal(q, b * 4, tolerance = 1e-6)
})
