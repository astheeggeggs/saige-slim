# The LEDGER #139 fix. On chromosome X with the male rewrite on, a
# marker with ANY missing call got two different silent wrong answers depending
# on the reader. applyMaleXnonPAR recounted altCounts/altFreq from a vector that
# still held the reader's MISSING SENTINEL (-1 from PLINK's genoMaps, NA_REAL
# from pgenlib), and the garbage count drove the QC gates and the imputed value:
# PLINK's rows VANISHED as ordinary QC (6 of P42's 48 markers survived, exactly
# the six with no missing call), PGEN's went NaN and either wrote `BETA=NA SE=NA`
# beside `nan` frequencies (single-ratio model) or killed the run (categorical
# ratios: "no categorical variance ratio covers MAC nan"), and --condition on
# such a marker conditioned on a garbage vector. The recount is now over the
# non-missing entries only, with the reader's own denominator.
#
# Fixture: genotype_100markers_missing (60 of its 100 markers carry missing
# calls, 42 of the 48 that pass the missing-rate gate) with the bim/pvar
# chromosome column rewritten to X (extdata/input/genotype_100markers_missing_chrX
# .{bim,pvar}; the .bed/.pgen are row-indexed and shared), the bundled male list
# (500 of the model's 1,000) and the categorical-ratio golden model. Every test
# here is red at unpatched HEAD: the PGEN arm exits non-zero and the PLINK arm
# tests 6 markers. Goldens P67/P68 pin the bytes; this file pins the mechanism.
# Tier: standard.

xm_cli <- run_saige_cli
xm_ok <- expect_cli_ok
xm_read <- read_saige_table

xm_model <- function() {
  m <- file.path(get("saige_goldens"), "step1", "binary_sparseGRM_a")
  c(paste0("--GMMATmodelFile=", m, ".rda"),
    paste0("--varianceRatioFile=", m, ".varianceRatio.txt"))
}
XM_COMMON <- c("--minMAF=0", "--minMAC=1", "--LOCO=FALSE", "--is_output_moreDetails=TRUE")
XM_X <- c("--chrom=X", "--is_rewrite_XnonPAR_forMales=TRUE",
          "--sampleFile_male=./input/sampleid_males.txt",
          "--X_PARregion=60001-2699520,154931044-155260560")
xm_geno <- function(reader, chrX) {
  suf <- if (chrX) "_chrX" else ""
  if (reader == "plink") {
    c("--bedFile=./input/genotype_100markers_missing.bed",
      paste0("--bimFile=./input/genotype_100markers_missing", suf, ".bim"),
      "--famFile=./input/genotype_100markers_missing.fam", "--AlleleOrder=alt-first")
  } else {
    c("--pgenFile=./input/genotype_100markers_missing.pgen",
      paste0("--pvarFile=./input/genotype_100markers_missing", suf, ".pvar"),
      "--psamFile=./input/genotype_100markers_missing.psam", "--AlleleOrder=alt-first")
  }
}
xm_run <- function(reader, extra = character()) {
  out <- tempfile(paste0("xm_", reader, "_"), fileext = ".txt")
  st <- xm_cli("step2_SPAtests.R", c(
    xm_geno(reader, chrX = TRUE), xm_model(), XM_COMMON, XM_X,
    paste0("--SAIGEOutputFile=", out), extra
  ))
  xm_ok(st)
  list(out = out, rows = xm_read(out),
       log = readLines(attr(st, "log"), warn = FALSE))
}
xm_num <- function(x) suppressWarnings(as.numeric(x))

test_that("LEDGER #139: chrX doubling with missing calls -- both readers agree, nothing vanishes", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  p <- xm_run("plink")
  g <- xm_run("pgen")
  expect_true(any(grepl("males are found in the test", p$log, fixed = TRUE)))
  expect_true(any(grepl("males are found in the test", g$log, fixed = TRUE)))
  # the invariant the row is about: same data, same flags, one answer
  expect_identical(readLines(p$out), readLines(g$out))
  # nothing vanished: the marker set of the autosomal run of the same file
  # (P42's golden). Doubling raises a MAC; it never fails a gate.
  a <- xm_read(file.path(get("saige_goldens"), "step2", "plink_marker_missing.txt"))
  expect_identical(p$rows$MarkerID, a$MarkerID)
  # the regime: most tested rows carry missing calls
  expect_gt(sum(xm_num(p$rows$MissingRate) > 0), 30)
  # every number is a number: no sentinel reached a count
  expect_false(any(grepl("nan", readLines(p$out), fixed = TRUE)))
  expect_true(all(is.finite(xm_num(p$rows$AC_Allele2))))
  expect_true(all(is.finite(xm_num(p$rows$AF_Allele2))))
  # the count is the frequency's: AC = AF * 2N after imputation (N = 1,000)
  expect_equal(xm_num(p$rows$AC_Allele2), xm_num(p$rows$AF_Allele2) * 2000, tolerance = 1e-3)
  # and the doubling was live: males carry alleles, so X counts exceed autosomal ones
  expect_gt(sum(xm_num(p$rows$AC_Allele2) > xm_num(a$AC_Allele2)), 30)
  # the corrected stdout line (LEDGER #150): neither the #139 sentence ("won't
  # be applied") nor its first correction ("replaces --minInfo and --minMAF",
  # "minMAF to 1/(2N)") may come back -- the flag changes neither gate
  expect_true(any(grepl("--minInfo and --minMAF apply to the X chromosome markers at the values given",
                        p$log, fixed = TRUE)))
  expect_false(any(grepl("won't be applied", p$log, fixed = TRUE)))
  expect_false(any(grepl("minMAF to 1/(2N)", p$log, fixed = TRUE)))
})

test_that("LEDGER #139: conditioning on a chrX marker with missing calls -- both readers agree", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # rs23: 609 ALT alleles, 5 missing calls in the analysed samples (P42)
  p <- xm_run("plink", extra = "--condition=rs23")
  g <- xm_run("pgen", extra = "--condition=rs23")
  expect_identical(readLines(p$out), readLines(g$out))
  expect_true(all(c("BETA_c", "SE_c", "p.value_c") %in% colnames(p$rows)))
  expect_false(any(grepl("nan", readLines(p$out), fixed = TRUE)))
  others <- p$rows$MarkerID != "rs23"
  expect_true(all(is.finite(xm_num(p$rows$p.value_c[others]))))
  # the adjustment is live: the conditional column is not the marginal one
  expect_false(identical(p$rows$p.value[others], p$rows$p.value_c[others]))
})
