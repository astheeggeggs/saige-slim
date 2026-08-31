# Characterization tests for the two AF > 0.5 bugs, LEDGER #28 and #29.
#
# These assert what the code does TODAY, including the broken behaviour. Read
# the LEDGER entries before "correcting" one; both flip in the dedicated fix
# phase, together with golden P11 (plink_marker_reffirst.txt).
#
# The construction: the bundled bed has no marker with AF_Allele2 > 0.5 in its
# native orientation (max 0.498), so every marker is minor-coded and neither
# bug can fire. Reading the SAME bed with --AlleleOrder=ref-first swaps which
# allele is counted, giving a mirror-image dataset in which EVERY marker has
# AF > 0.5. That is a proxy for "ALT is the major allele", not the real thing
# (a native fixture would need an edited .bed, not just a .bim) — but it is a
# sound one, and the first test below is the control that proves it: under
# --is_noadjCov=FALSE the flip is EXACTLY invariant, so the reader, the
# genotype flip in imputeGenoAndFlip, and the score test are all correct and
# any difference seen elsewhere is attributable to the bug under test.
#
# Tier: standard for the CLI runs; the gate-1 arithmetic reconciliation at the
# bottom is pure R and runs in the fast tier.

AO_BASE <- function(alleleOrder) c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  paste0("--AlleleOrder=", alleleOrder),
  "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
  "--is_output_moreDetails=TRUE"
)

# One CLI run per (allele order, extra args), cached across the blocks below.
ao_cache <- new.env(parent = emptyenv())
ao_run <- function(alleleOrder, extra = character(0)) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  key <- paste(c(alleleOrder, extra), collapse = " ")
  if (!is.null(ao_cache[[key]])) return(ao_cache[[key]])
  out <- file.path(tempfile("step2_ao_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(AO_BASE(alleleOrder), extra,
                            paste0("--SAIGEOutputFile=", out)))
  expect_cli_ok(status)
  res <- list(path = out, tab = read_saige_table(out))
  ao_cache[[key]] <- res
  res
}

# The alt-first / default-args run is byte-identical to golden P1, so read that
# rather than paying for a fourth CLI run.
alt_default <- function() {
  skip_if_no_goldens()
  read_saige_table(file.path(saige_goldens, "step2", "plink_marker_binary.txt"))
}

# read_saige_table keeps every column as the literal printed string (so that
# golden comparison sees exactly what was written). Numeric assertions below
# parse explicitly; string identity of a printed column is the STRONGER claim
# and is used wherever the two runs should agree exactly.
ao_num <- function(x) suppressWarnings(as.numeric(x))

# Pair the two orientations on MarkerID. Only markers surviving BOTH runs are
# compared (the ref-first runs lose 7 to LEDGER #29, asserted separately).
ao_pair <- function(alt, ref) {
  m <- merge(as.data.frame(alt), as.data.frame(ref),
             by = "MarkerID", suffixes = c(".alt", ".ref"))
  m[order(m$MarkerID), ]
}

test_that("CONTROL: with --is_noadjCov=FALSE the allele flip is exactly invariant", {
  # The negative control for both bugs. The covariate-adjusted score test
  # builds gtilde by projecting out the covariates (intercept included), so it
  # centres the genotype empirically and never consults altFreq. Flipping the
  # allele order must therefore change NOTHING except the sign of the effect.
  # Exact equality, not a tolerance: the tested genotype vector is bit-identical
  # in the two runs, so every downstream double must be too.
  alt <- ao_run("alt-first", "--is_noadjCov=FALSE")$tab
  ref <- ao_run("ref-first", "--is_noadjCov=FALSE")$tab
  m <- ao_pair(alt, ref)
  expect_gt(nrow(m), 60)
  expect_true(all(ao_num(m$AF_Allele2.ref) > 0.5))   # the mirror is AF > 0.5

  # Printed-string identity: the strongest form of "nothing changed".
  expect_identical(m$SE.alt, m$SE.ref)
  expect_identical(m$var.alt, m$var.ref)
  expect_identical(m$p.value.alt, m$p.value.ref)
  # Sign convention: BETA/Tstat are negated exactly (IEEE negation is exact).
  expect_identical(ao_num(m$BETA.alt), -ao_num(m$BETA.ref))
  expect_identical(ao_num(m$Tstat.alt), -ao_num(m$Tstat.ref))
})

test_that("LEDGER #28: the default (noadjCov) path inflates variance at AF > 0.5", {
  # scoreTestFast_noadjCov centres on 2*altFreq with altFreq in the ALT
  # orientation while the genotype is minor-coded, so the centre is off by
  # 2*(1 - 2*MAF) and the variance picks up sum(mu2)*4*(1-2*MAF)^2.
  # is_fastTest=FALSE here: without a sparse GRM this path is taken either way,
  # which is itself the point — --is_fastTest=FALSE alone does NOT avoid it.
  alt <- alt_default()
  ref <- ao_run("ref-first")$tab
  m <- ao_pair(alt, ref)

  ratio <- ao_num(m$var.ref) / ao_num(m$var.alt)
  expect_gt(median(ratio), 100)         # observed ~178x at the time of writing
  expect_gt(max(ratio), 1000)           # observed ~5408x
  expect_true(all(ratio > 0.99))        # inflation only, never deflation

  # The score is mis-centred toward zero the same way, so p is driven to 1 and
  # every signal in the fixture is erased.
  expect_lt(min(ao_num(alt$p.value)), 0.05)
  expect_gt(min(ao_num(ref$p.value)), 0.05)
  expect_gt(median(ao_num(ref$p.value)), 0.9)

  # Self-concealing under --is_fastTest: the escalation at Main.cpp:440 fires
  # only when pass-1 p < pval_cutoff_for_fastTest (0.05), and the inflation
  # pushes every marker above it, so the exact re-test never runs and the
  # fastTest output is identical to the non-fast output.
  fast <- ao_run("ref-first", "--is_fastTest=TRUE")$tab
  mf <- ao_pair(ref, fast)
  expect_identical(mf$var.alt, mf$var.ref)
  expect_identical(mf$p.value.alt, mf$p.value.ref)
})

test_that("LEDGER #29: markers at exactly MAC = minMAC are dropped when ALT is major", {
  alt <- alt_default()
  ref <- ao_run("ref-first")$tab
  dropped <- setdiff(alt$MarkerID, ref$MarkerID)

  # Exactly the MAC = 1 singletons, and nothing else, under --minMAC=1.
  expect_length(dropped, 7)
  ac <- ao_num(alt$AC_Allele2)
  nn <- ao_num(alt$N_case) + ao_num(alt$N_ctrl)
  mac <- pmin(ac, 2 * nn - ac)
  expect_true(all(mac[alt$MarkerID %in% dropped] == 1))
  expect_true(all(mac[!alt$MarkerID %in% dropped] != 1))
  expect_length(setdiff(ref$MarkerID, alt$MarkerID), 0L)
})

test_that("LEDGER #29: the gate-1 MAC round trip is reproducible in exact arithmetic", {
  # Pure R, fast tier: reproduce src/Main.cpp:328-329 on the counts of the
  # dropped singletons and show the value lands below the cutoff it equals.
  # PLINK.cpp:303-307 (ref-first) and a major ALT under alt-first both reach
  # gate 1 through the same `1 - f` subtraction, which is what loses precision;
  # altCounts itself stays exact.
  n <- 1000
  gate1_mac <- function(altFreq, n) {
    MAF <- min(altFreq, 1 - altFreq)
    MAF * n * (1 - 0) * 2                 # missingRate = 0 in this fixture
  }
  # minor-coded (what every golden marker is): exact integers
  for (ac in c(1, 2, 10, 100)) {
    expect_identical(gate1_mac(ac / n / 2, n), as.numeric(ac))
  }
  # ALT major: the singleton lands 1.1e-13 short of 1 and fails `MAC < cutoff`
  mac1 <- gate1_mac(1 - 1 / n / 2, n)
  expect_lt(mac1, 1)
  expect_equal(mac1, 1, tolerance = 1e-12)
  # ...while neighbouring counts happen to round high and survive. The error is
  # an ulp coin-flip, not systematically signed — hence "data-dependent".
  expect_gt(gate1_mac(1 - 2 / n / 2, n), 2)
  expect_gt(gate1_mac(1 - 10 / n / 2, n), 10)

  # Gate 2 (Main.cpp:355) is exact, but runs only after gate 1 has dropped it.
  expect_identical(min(1999, 2 * n - 1999), 1)
})

test_that("golden P11: the ref-first characterization run is stable", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  ref <- ao_run("ref-first")
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_reffirst.txt"),
    ref$path
  )
})
