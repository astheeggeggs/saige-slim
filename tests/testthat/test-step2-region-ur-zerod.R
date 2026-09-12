# LEDGER #129 (FIXED 2026-09-12, inverted in the fixing commit): the
# dosage-zeroing filter used to be applied to the WEIGHTED ultra-rare
# pseudo-marker, in weight units. src/Main.cpp passed the collapsed `genoURVec`
# (max of w * G over the set's ultra-rare markers) through the filter with
# --dosage_zerod_cutoff (0.2) and --dosage_zerod_MAC_cutoff (10), cleaning every
# entry <= 0.2 when the WEIGHTED MAC was <= 10. A custom or transform set whose
# support is all ultra-rare and whose weights put every carrier at or below 0.2
# therefore lost its pseudo-marker and its row -- silently: not refused, not a
# SKAT failure, the pooled row NA and the #127 sidecar 0 / 0 / none. So whether
# a set was TESTED depended on the scale of its weights, which a burden p does
# not depend on at all. The filter is no longer applied to the pseudo-marker
# (owner's decision, fix (b) of the two on the row): it is not a dosage, and
# each constituent marker has already been through the same filter as a dosage
# -- at the defaults the two cutoffs are co-extensive, both MAC 10.
#
# What these three configurations now assert is the invariance itself, three
# ways over one fixture: weights 0.01 and weights 0.5, at collapse 10 and at
# collapse 0, all give the SAME burden p, and the BETAs stand in the 50x ratio
# of the weights. Golden P64 is the collapse-10 / weights-0.01 corner.
# Fixture: GENET var rs28 rs32 (two singletons, cell MAC 2), one custom set.
# Tier: standard -- FOUR step-2 runs of a two-marker region, which is what the
# (weight, collapse) cache below makes true; the assertions read each of the
# four configurations more than once.

# One CLI run per (weight, collapse) pair, cached per process: the assertions
# below read the same four configurations several times each, and without this
# the file spent 7 runs on 4 answers while its header claimed four.
urzerod_cache <- new.env(parent = emptyenv())
urzerod_run <- function(w, collapse) {
  key <- paste(w, collapse)
  if (!is.null(urzerod_cache[[key]])) return(urzerod_cache[[key]])
  res <- urzerod_run_uncached(w, collapse)
  urzerod_cache[[key]] <- res
  res
}
urzerod_run_uncached <- function(w, collapse) {
  gf <- tempfile(fileext = ".txt")
  writeLines(c("GENET var rs28 rs32", "GENET anno lof lof",
               sprintf("GENET weight:w %s %s", w, w)), gf)
  out <- tempfile(fileext = ".txt")
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE", "--is_fastTest=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
            paste0("--sparseGRMFile=", sg),
            paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            paste0("--groupFile=", gf), "--annotation_in_groupTest=lof",
            "--maxMAF_in_groupTest=0.01", "--r.corr=1",
            paste0("--MACCutoff_to_CollapseUltraRare=", collapse),
            paste0("--SAIGEOutputFile=", out))
  st <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st)
  x <- read_saige_table(out)
  # the #127 sidecar is written only when a region has something to disclose
  sf <- paste0(out, ".pooledTests.txt")
  side <- if (file.exists(sf)) read.table(sf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  list(cells = x[x$Group != "Cauchy", ], pooled = x[x$Group == "Cauchy", ], side = side)
}

test_that("weights 0.01: the all-ultra-rare set is tested at collapse 10 and at collapse 0", {
  skip_unless_level("standard")
  r10 <- urzerod_run(0.01, 10)
  # INVERTED. Was: nrow(cells) == 0, pooled "NA", sidecar 0 / 0 / "none" --
  # the set silently untested because 0.01 <= --dosage_zerod_cutoff.
  expect_equal(nrow(r10$cells), 1L)
  expect_equal(r10$cells$MAC, "2")
  # numeric, NOT the printed string: a p-value is a computed quantity and the
  # suite has a compare mode for exactly this (CLAUDE.md's byte-vs-numeric
  # trap). The string form of this assertion passed locally and would have been
  # a platform bet on CI.
  expect_equal(as.numeric(r10$cells$Pvalue_Burden), 0.2522646, tolerance = 1e-6)
  expect_equal(as.numeric(r10$pooled$Pvalue_Burden), 0.2522646,   # CCT of one p is that p
               tolerance = 1e-5)
  expect_null(r10$side)   # nothing to disclose once the set is tested
  # ANTI-VACUITY (the LEDGER #30 shape): without these, every assertion in this
  # file is also satisfied by a build where the ultra-rare collapse never
  # happens at all -- the two singletons would then be tested as rare markers,
  # giving one cell at MAC 2 with the same p. These two lines are what make the
  # file a test of the COLLAPSE.
  expect_equal(r10$cells$Number_ultra_rare, "2")
  expect_equal(r10$cells$Number_rare, "0")
  # and the collapse-0 arm, which was ALWAYS tested, is unchanged -- it is the
  # control that says the fix moved the broken corner and nothing else
  r0 <- urzerod_run(0.01, 0)
  expect_equal(nrow(r0$cells), 1L)
  expect_equal(r0$cells$MAC, "2")
  expect_true(is.finite(as.numeric(r0$cells$Pvalue_Burden)))
  expect_null(r0$side)
  # the other half of the anti-vacuity pair: at collapse 0 the SAME markers are
  # rare, so the two arms differ in which code path produced the identical p
  expect_equal(r0$cells$Number_ultra_rare, "0")
  expect_equal(r0$cells$Number_rare, "2")
})

test_that("the burden p is invariant to the weight scale and to the collapse cutoff", {
  skip_unless_level("standard")
  r10 <- urzerod_run(0.5, 10)
  r0 <- urzerod_run(0.5, 0)
  expect_equal(nrow(r10$cells), 1L)
  expect_equal(nrow(r0$cells), 1L)
  expect_equal(r10$cells$MAC, "2")
  # the weight scale cannot move a burden p: 0.5 at collapse 0 equals 0.01 at
  # collapse 0 (same path, inputs scaled), and the collapsed pseudo-marker of
  # two singletons equals their sum (no sample carries both) -- two different
  # code paths, so compared at printed precision, not for bit-identity
  p001 <- as.numeric(urzerod_run(0.01, 0)$cells$Pvalue_Burden)
  expect_equal(as.numeric(r0$cells$Pvalue_Burden), p001, tolerance = 1e-6)
  expect_equal(as.numeric(r10$cells$Pvalue_Burden), p001, tolerance = 1e-6)
  # ADDED with the #129 fix: the fourth corner, which used to have no p at all.
  # All four now agree, so the invariance holds over the whole 2 x 2 -- and it
  # is the weights-0.01-at-collapse-10 corner that the bug was hiding in.
  p001_10 <- as.numeric(urzerod_run(0.01, 10)$cells$Pvalue_Burden)
  expect_equal(p001_10, p001, tolerance = 1e-6)
  # the BETAs, by contrast, DO scale with the weights -- 1/w on a burden
  # coefficient -- which is what says the p-invariance is not two runs of the
  # same arithmetic: 232.633 at 0.01 against 4.65266 at 0.5, a factor of 50
  b001 <- as.numeric(urzerod_run(0.01, 10)$cells$BETA_Burden)
  b05 <- as.numeric(r10$cells$BETA_Burden)
  expect_equal(b001 / b05, 50, tolerance = 1e-4)
})
