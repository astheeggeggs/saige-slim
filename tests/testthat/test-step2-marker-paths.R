# The step-2 MARKER-path branches that no golden reached before S0.
#
# Every pre-existing single-variant golden (P1/P4/P7/P11/P12) runs the same one
# arm of the getMarkerPval cascade -- scoreTestFast_noadjCov -- so a change to
# any other arm left every golden byte-identical while being completely
# untested. That is the coverage gap recorded as F8 in
# docs/design/step2-threading-verdict.md, and it is what made the C3 allocation
# work stop short of the sparse scoreTest and its PCG.
#
# Each case here is P1's command with ONE thing changed, so a diff against
# plink_marker_binary.txt isolates the branch. The branch counts quoted in the
# comments were measured with an instrumented build (counters on each arm, run
# at 61ac4be9 and discarded); they are recorded in harness/make_goldens.sh
# beside the runs that produce these goldens.
#
# Tier: standard (runs the CLI). Kept out of test-step2-goldens.R so the two
# files run in parallel rather than serialising ~14 + 7 CLI runs into one.

MP_PLINK <- c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first"
)
MP_GRM <- paste0("output/sparseGRM_relatednessCutoff_0.125_1000_",
                 "randomMarkersUsed.sparseGRM.mtx")
MP_SG <- c(
  paste0("--sparseGRMFile=", MP_GRM),
  paste0("--sparseGRMSampleIDFile=", MP_GRM, ".sampleIDs.txt")
)
MP_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
               "--is_output_moreDetails=TRUE")
mp_model <- function(prefix) c(
  paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", paste0(prefix, ".rda"))),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", paste0(prefix, ".varianceRatio.txt")))
)

# The two bundled markers that ARE the model's categorical covariates a9/a10
# (R2 = 1.000000; CLAUDE.md records the identification). On any path that
# regresses the covariates out of the genotype, the adjusted vector is exactly
# zero, so Tstat and var are sums that cancel to nothing: in these goldens
# |Tstat| ~ 3e-15 and var is 0 or ~3e-14, against [0.021, 14.5] and
# [0.020, 39.9] for the 68 ordinary markers -- thirteen orders of magnitude
# below the smallest real value in the same column.
#
# BETA is then Tstat/var, a RATIO of two noise quantities, which lands at O(0.1)
# while carrying no information at all: -0.110908 here, something else entirely
# under a different summation order. That is why these rows cannot be compared
# across platforms at any relative tolerance, and why an absolute floor on the
# value would not have caught the BETA column -- the failure is not that the
# number is small, it is that the number is arbitrary. CI found it as
# rel_err 1.2e-01 [BETA] / 1.6e-01 [Tstat] on Linux against macOS-generated
# goldens.
#
# LEDGER #44 IS NOW FULLY FIXED and this exclusion is RETIRED. Both halves:
#
#   P27 (sparse scoreTest): var is 0 exactly, so `var1 <= DBL_MIN` fires.
#   P30 (scoreTestFast):    var is ~2.7e-14, ~294 orders of magnitude ABOVE
#       DBL_MIN, so that guard never could. The two "fast" paths build var2 as a
#       sum of signed terms, so the same degeneracy appears as CANCELLATION, and
#       they now ask whether |var2| is at the noise floor of its own summation
#       (n*eps times the sum of the term magnitudes).
#
# Both now print BETA/SE as the literal "NA". What is left on those rows is
# Tstat ~3e-15 and var ~2.7e-14 -- genuinely-near-zero values, handled in the
# comparator by rel_err()'s column-scaled absolute floor rather than by excluding
# the rows. So nothing is skipped any more: all 70 markers are compared in both
# modes, and the characterization test at the bottom of this file pins the NAs.

mp_case <- function(golden_base, extra_args, drop_markers = character()) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step2_mp_"), golden_base)
  dir.create(dirname(out), recursive = TRUE)
  expect_cli_ok(run_saige_cli(
    "step2_SPAtests.R",
    c(MP_PLINK, MP_COMMON, extra_args, paste0("--SAIGEOutputFile=", out))
  ))
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", golden_base), out,
    drop_markers = drop_markers
  )
  invisible(out)
}

test_that("marker-path sparse scoreTest + PCG matches golden (P27)", {
  # 70/70 markers take the sparse branch; 70 PCG solves, 94 iterations. Needs
  # the binary_sparseGRM_vr model: sparse VR rows AND null_noXadj rows in one
  # file, which no other committed variance-ratio file has.
  mp_case("plink_marker_sparse.txt",
          c(mp_model("binary_sparseGRM_vr"), MP_SG, "--is_fastTest=FALSE"))
})

test_that("marker fastTest re-test matches golden (P28)", {
  # Main.cpp:436-475. 73 getMarkerPval calls for 70 markers: 3 re-tests.
  mp_case("plink_marker_fastTest.txt",
          c(mp_model("binary_sparseGRM_vr"), MP_SG, "--is_fastTest=TRUE"))
})

test_that("marker fastTest re-test, sparse arm, matches golden (P29)", {
  # max_MAC_for_ER=0 turns ER off, which is the only way the two ultra-rare
  # markers reach the re-test (it is gated on MAC > MACCutoffforER). 75 calls:
  # 2 re-tests take the sparse branch, 3 scoreTestFast; ER runs zero times.
  mp_case("plink_marker_fastTest_sparse.txt",
          c(mp_model("binary_sparseGRM_vr"), MP_SG,
            "--is_fastTest=TRUE", "--max_MAC_for_ER=0"))
})

test_that("marker scoreTestFast (is_noadjCov=FALSE) matches golden (P30)", {
  # The third arm of SAIGE_test.cpp:514-528, on all 70 markers.
  mp_case("plink_marker_noadjCov.txt",
          c(mp_model("binary_sparseGRM_a"), "--is_noadjCov=FALSE"))
})

test_that("marker Firth path matches golden (P31)", {
  # 5 of 70 markers qualify at 0.05: 3 via the SPA path, 2 via the ER path.
  # Do not raise the cutoff (LEDGER #26).
  mp_case("plink_marker_firth.txt",
          c(mp_model("binary_sparseGRM_a"),
            "--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"))
})

test_that("quantitative marker path, unconditional, matches golden (P32)", {
  mp_case("plink_marker_quantitative.txt", mp_model("quant_sparseGRM"))
})

test_that("ER resampling branch matches golden (P33)", {
  # Raising the ER cap to 200 admits rs38 (MAC 102, ~98 carriers), the only
  # marker on this fixture that can reach n_total = 2^k > 2e6 and so the
  # resampling arm of ComputeExact::Run -- the arm that draws from R's RNG
  # after SL_setseed(1). Byte-equality with the golden is therefore also the
  # determinism check on that draw.
  mp_case("plink_marker_er.txt",
          c(mp_model("binary_sparseGRM_a"), "--max_MAC_for_ER=200"))
})

test_that("the pread fallback reproduces a marker golden end to end (S5)", {
  # SAIGE can read a .bed two ways -- mmap, or a pread fallback when the map
  # fails -- and until S5 the fallback was reachable only by an mmap failure,
  # which does not happen for a small local file. So one of the two decode paths
  # in this package was exercised by nothing at all, on any platform, in any
  # golden. SAIGE_PLINK_NO_MMAP forces it.
  #
  # test-unit-genoreader.R compares the two paths' decoded genotype vectors for
  # five markers, including the LEDGER #16 mid-stream row-0 pattern. This is the
  # same claim over a whole run: every marker, every column, through the CLI.
  # Byte-exactness IS the right standard here -- same file, same offsets, same
  # decode -- unlike a comparison between two different code paths.
  #
  # The ER golden is the one chosen because it has no covariate-collinear
  # degenerate rows (the default noadjCov path does not adjust), so it compares
  # in full in both modes.
  Sys.setenv(SAIGE_PLINK_NO_MMAP = "1")
  on.exit(Sys.unsetenv("SAIGE_PLINK_NO_MMAP"), add = TRUE)
  # This case is about the READ path. It needs no thread pinning any more: S7
  # made --nThreads the only knob and it defaults to 1, so mp_case() is serial
  # unless it asks otherwise. (Before S7 this had to unset SAIGE_MARKER_THREADS,
  # because the pread fallback refuses to clone and a threaded run would have
  # been answered by the degrade-to-one-thread notice, testing something else.)
  mp_case("plink_marker_er.txt",
          c(mp_model("binary_sparseGRM_a"), "--max_MAC_for_ER=200"))
})

test_that("every new marker golden differs from P1 in the intended direction", {
  # Guard against a golden that gates nothing: --max_MAC_for_ER=25 was proposed
  # for the ER case and is byte-identical to P1 on this fixture. Each golden
  # here must differ from P1 where its branch says it should, and agree with it
  # on the allele accounting, which no branch in the cascade may touch.
  skip_if_no_goldens()
  g <- function(f) read_saige_table(file.path(saige_goldens, "step2", f))
  p1 <- g("plink_marker_binary.txt")
  # branches that change the score test itself -> the p-value must move
  for (f in c("plink_marker_sparse.txt", "plink_marker_fastTest.txt",
              "plink_marker_fastTest_sparse.txt", "plink_marker_noadjCov.txt",
              "plink_marker_er.txt")) {
    d <- g(f)
    expect_identical(d$MarkerID, p1$MarkerID, label = paste("same markers:", f))
    expect_identical(d$AC_Allele2, p1$AC_Allele2, label = paste("same AC:", f))
    expect_false(identical(d$p.value, p1$p.value), label = paste("p differs:", f))
  }
  # Firth is the exception BY CONSTRUCTION: it refines the effect size after the
  # p-value is settled, so p.value must be untouched and BETA/SE must move. That
  # is the same claim test-step2-firth.R makes about the live run; here it is
  # made about the frozen golden.
  fi <- g("plink_marker_firth.txt")
  expect_identical(fi$MarkerID, p1$MarkerID)
  expect_identical(fi$p.value, p1$p.value)
  expect_false(identical(fi$BETA, p1$BETA))
  expect_false(identical(fi$SE, p1$SE))
  # 5 markers qualify at pCutoffforFirth=0.05, and only those change
  changed <- which(fi$BETA != p1$BETA)
  expect_identical(changed, which(suppressWarnings(as.numeric(p1$p.value)) <= 0.05))
  expect_length(changed, 5L)
})

test_that("LEDGER #44: both degeneracies now print NA", {
  # rs20 and rs22 ARE the categorical covariates a9 and a10 (verified: R^2 = 1
  # against the model matrix, next highest marker 0.022), so every
  # covariate-adjusting path tests a genotype that is exactly explained away.
  #
  # #44 is FIXED on BOTH paths, by two different predicates, because the
  # degeneracy has two different shapes:
  #
  #   P27 (sparse scoreTest): var is 0 EXACTLY -> `var1 <= DBL_MIN` fires.
  #   P30 (scoreTestFast):    var is ~2.7e-14, ~294 orders ABOVE DBL_MIN, so the
  #       absolute guard cannot fire. var2 there is a sum of signed terms, so the
  #       degeneracy is CANCELLATION, caught by comparing |var2| to the noise
  #       floor of its own summation. Measured across all 70 markers on that
  #       configuration: these two give 1.7e-16 and 4.2e-16, every other marker
  #       >= 0.237 -- 14.8 orders of separation.
  #
  # Tstat and var stay as they are, deliberately: the score really is ~0 here and
  # that is information. They are portable because rel_err() has a column-scaled
  # absolute floor, not because these rows are skipped -- nothing is skipped now.
  skip_if_no_goldens()
  g <- function(f) read_saige_table(file.path(saige_goldens, "step2", f))
  sp <- g("plink_marker_sparse.txt")
  i <- match(c("rs20", "rs22"), sp$MarkerID)
  expect_false(anyNA(i))
  # FIXED half: literal "NA", which is also the point -- a NaN's printed
  # spelling is the C library's (glibc "-nan" vs macOS "nan"), so the old value
  # was not byte-portable. "NA" says the same thing on every platform.
  expect_identical(sp$BETA[i], c("NA", "NA"))
  expect_identical(sp$SE[i], c("NA", "NA"))
  expect_identical(sp$p.value[i], c("1.000000E+00", "1.000000E+00"))
  expect_true(all(abs(suppressWarnings(as.numeric(sp$Tstat[i]))) < 1e-12))
  # The other adjusting path, now caught by the cancellation test rather than
  # left printing BETA=-0.11 SE=6.1e6:
  na <- g("plink_marker_noadjCov.txt")
  expect_identical(na$BETA[i], c("NA", "NA"))
  expect_identical(na$SE[i], c("NA", "NA"))
  expect_true(all(suppressWarnings(as.numeric(na$var[i])) < 1e-13))
  # and the NON-adjusting default path is untouched -- it never degenerates,
  # which is why P1 reports these two as ordinary null markers
  p1 <- g("plink_marker_binary.txt")
  k <- match(c("rs20", "rs22"), p1$MarkerID)
  expect_true(all(is.finite(suppressWarnings(as.numeric(p1$BETA[k])))))
  # and nothing else on either path is degenerate
  expect_true(all(is.finite(suppressWarnings(as.numeric(sp$BETA[-i])))))
  expect_true(all(is.finite(suppressWarnings(as.numeric(na$SE[-i])))))
})
