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

mp_case <- function(golden_base, extra_args) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step2_mp_"), golden_base)
  dir.create(dirname(out), recursive = TRUE)
  expect_cli_ok(run_saige_cli(
    "step2_SPAtests.R",
    c(MP_PLINK, MP_COMMON, extra_args, paste0("--SAIGEOutputFile=", out))
  ))
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", golden_base), out
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

test_that("the covariate-collinear fixture markers degenerate as LEDGER #44 records", {
  # rs20 and rs22 ARE the categorical covariates a9 and a10 (verified: R^2 = 1
  # against the model matrix, next highest marker 0.022). Every covariate-
  # adjusting path therefore tests a genotype that is exactly explained away,
  # and the variance underflows. The p-value guard holds (p = 1) but BETA/SE
  # are unguarded, so the sparse path prints -inf/nan. This is pinned, not
  # fixed, per the branch policy -- if it changes, update LEDGER #44 with it.
  skip_if_no_goldens()
  g <- function(f) read_saige_table(file.path(saige_goldens, "step2", f))
  sp <- g("plink_marker_sparse.txt")
  i <- match(c("rs20", "rs22"), sp$MarkerID)
  expect_false(anyNA(i))
  expect_identical(sp$BETA[i], c("-inf", "-inf"))
  # NaN's printed spelling is the C library's: glibc writes "-nan" where macOS
  # writes "nan" for the same sign-bit-set NaN. The committed golden carries
  # whichever platform generated it, so accept either.
  expect_true(all(grepl("^-?nan$", sp$SE[i])))
  expect_identical(sp$p.value[i], c("1.000000E+00", "1.000000E+00"))
  # the same two markers on the other covariate-adjusting path: finite, but the
  # variance has underflowed to ~1e-14 and the SE is ~1e6
  na <- g("plink_marker_noadjCov.txt")
  expect_true(all(suppressWarnings(as.numeric(na$var[i])) < 1e-13))
  expect_true(all(suppressWarnings(as.numeric(na$SE[i])) > 1e6))
  # and nothing else on either path is degenerate
  expect_true(all(is.finite(suppressWarnings(as.numeric(sp$BETA[-i])))))
  expect_true(all(is.finite(suppressWarnings(as.numeric(na$SE[-i])))))
})
