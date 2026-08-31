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

test_that("LEDGER #28 FIXED: the default (noadjCov) path is flip-invariant", {
  # Was: scoreTestFast_noadjCov centred on 2*altFreq with altFreq in the ALT
  # orientation while the genotype is minor-coded, so the centre was off by
  # 2*(1 - 2*MAF) and the variance picked up sum(mu2)*4*(1-2*MAF)^2 (observed
  # median inflation 178x, max 5408x, every signal erased).
  # Now: the flip flag is threaded down from the Main.cpp marker loop and the
  # centre is the frequency of the allele t_GVec actually counts, so this path
  # matches the covariate-adjusted CONTROL arm above EXACTLY.
  # is_fastTest=FALSE here: without a sparse GRM this path is taken either way,
  # which was the point — --is_fastTest=FALSE alone did NOT avoid the bug.
  alt <- alt_default()
  ref <- ao_run("ref-first")$tab
  m <- ao_pair(alt, ref)
  expect_gt(nrow(m), 60)
  expect_true(all(ao_num(m$AF_Allele2.ref) > 0.5))   # the mirror is AF > 0.5

  # The same printed-string identity the CONTROL arm asserts. This is now the
  # strong claim: the two orientations agree bit-for-bit, not to a tolerance.
  expect_identical(m$SE.alt, m$SE.ref)
  expect_identical(m$var.alt, m$var.ref)
  expect_identical(m$p.value.alt, m$p.value.ref)
  expect_identical(ao_num(m$BETA.alt), -ao_num(m$BETA.ref))
  expect_identical(ao_num(m$Tstat.alt), -ao_num(m$Tstat.ref))

  # Signal is preserved rather than erased: the ref-first run recovers the same
  # minimum p-value as alt-first (it was driven above 0.05 by the inflation).
  expect_lt(min(ao_num(ref$p.value)), 0.05)

  # The self-concealment is gone, and this is the sharpest evidence of it.
  # BEFORE the fix the inflation drove every marker's pass-1 p above
  # pval_cutoff_for_fastTest (0.05), so the escalation at src/Main.cpp:440 never
  # fired and the fastTest output was identical to the non-fast output for all
  # 63 markers. AFTER the fix, signal is preserved, so exactly the markers below
  # the cutoff get re-tested by the exact second pass.
  #
  # Assert the partition exactly: differing <=> pass-1 p < cutoff. Not
  # bit-identity between fast and non-fast for the escalated rows -- those ARE
  # two different code paths (CLAUDE.md's standing trap), so they are compared
  # numerically and only required to agree to within the refinement they make.
  fast <- ao_run("ref-first", "--is_fastTest=TRUE")$tab
  mf <- ao_pair(ref, fast)
  escalated <- ao_num(mf$p.value.alt) < 0.05
  expect_gt(sum(escalated), 0)          # the re-test is now REACHABLE (was not)

  # Below the cutoff: pass 1 stands, byte for byte.
  expect_identical(mf$var.alt[!escalated], mf$var.ref[!escalated])
  expect_identical(mf$p.value.alt[!escalated], mf$p.value.ref[!escalated])

  # Above it: the exact re-test ran and moved the answer, but only slightly --
  # it is a refinement of the same statistic, not a different one. Note the
  # moved set is a SUBSET of the escalated set, not equal to it: at the time of
  # writing 4 markers sit below the cutoff and 3 move, because rs56 (MAC 2)
  # re-tests to the same 6-significant-digit string. So assert containment plus
  # reachability, never "all escalated markers differ".
  moved <- mf$p.value.alt != mf$p.value.ref
  expect_true(all(escalated[moved]))    # nothing outside the cutoff can move
  expect_gt(sum(moved), 0)              # and the second pass demonstrably runs
  expect_lt(max(abs(ao_num(mf$var.alt[escalated]) - ao_num(mf$var.ref[escalated])) /
                ao_num(mf$var.ref[escalated])), 0.05)
})

test_that("LEDGER #29 FIXED: markers at exactly MAC = minMAC survive when ALT is major", {
  # Was: gate 1 derived MAC from the FREQUENCY, so under --minMAC=1 the seven
  # MAC = 1 singletons (rs28, rs32, rs42, rs44, rs55, rs62, rs89) landed
  # ~1.1e-13 short of 1 and were dropped by the strictly-less-than gate.
  # Now: gate 1 derives MAC from the counts, and PLINK.cpp takes the ref-first
  # complement by exact integer arithmetic instead of rebuilding it from the
  # frequency, so both orientations keep all 70 markers.
  alt <- alt_default()
  ref <- ao_run("ref-first")$tab

  expect_length(setdiff(alt$MarkerID, ref$MarkerID), 0L)   # nothing dropped now
  expect_length(setdiff(ref$MarkerID, alt$MarkerID), 0L)
  expect_identical(nrow(alt), nrow(ref))

  # The seven are back, and they really are the MAC = 1 singletons.
  ac <- ao_num(alt$AC_Allele2)
  nn <- ao_num(alt$N_case) + ao_num(alt$N_ctrl)
  mac <- pmin(ac, 2 * nn - ac)
  expect_identical(sum(mac == 1), 7L)
  expect_true(all(alt$MarkerID[mac == 1] %in% ref$MarkerID))

  # And the recovered rows are not junk: the whole 70-marker pairing is exactly
  # flip-invariant, singletons included. AC is exactly complementary.
  m <- ao_pair(alt, ref)
  expect_identical(nrow(m), 70L)
  expect_identical(m$var.alt, m$var.ref)
  expect_identical(m$p.value.alt, m$p.value.ref)
  expect_identical(ao_num(m$BETA.alt), -ao_num(m$BETA.ref))
  expect_identical(ao_num(m$AC_Allele2.ref),
                   2 * (ao_num(m$N_case.alt) + ao_num(m$N_ctrl.alt)) -
                     ao_num(m$AC_Allele2.alt))
})

test_that("LEDGER #29: the OLD gate-1 MAC round trip, kept as a record of the bug", {
  # Pure R, fast tier. This reproduces the FORMER src/Main.cpp gate-1 expression
  # (MAF * n * (1-missingRate) * 2) to document why the seven singletons were
  # dropped. It is deliberately retained after the fix: it pins the float
  # behaviour that motivated the change, and it is the reason the replacement
  # had to come from the counts rather than the frequency.
  # PLINK.cpp (ref-first) and a major ALT under alt-first both reached gate 1
  # through the same `1 - f` subtraction, which is what loses precision.
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

  # Gate 2 (Main.cpp:355) was always exact, but ran only after gate 1 had
  # dropped the marker. The fix makes gate 1 use the same count-based form, so
  # the two gates now agree by construction -- exactly, for every count.
  gate1_fixed <- function(altCounts, nNonMissing) {
    min(altCounts, 2 * nNonMissing - altCounts)
  }
  for (ac in c(1, 2, 10, 100, 1900, 1990, 1998, 1999)) {
    expect_identical(gate1_fixed(ac, n), min(ac, 2 * n - ac))   # gate 1 == gate 2
  }
  # The singleton that used to fail `MAC < 1` now lands on exactly 1.
  expect_identical(gate1_fixed(1999, n), 1)
  expect_false(gate1_fixed(1999, n) < 1)
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
