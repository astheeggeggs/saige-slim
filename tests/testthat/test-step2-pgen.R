# PGEN (plink2 native format) reader, Phase 7. New scope, not a ledger fix:
# fresh code implementing GenoReader, carrying the CURRENT (post-#28/#29-fix)
# marker-path semantics rather than main's pre-refactor bugs.
#
# genotype_100markers.{pgen,pvar,psam} was generated once via
# `plink2 --bfile genotype_100markers --make-pgen --out genotype_100markers`;
# its pvar REF/ALT matches the alt-first bim's REF/ALT exactly (both label the
# same allele as ALT), so a PGEN run and the matching PLINK run over the SAME
# model/phenotype should reconcile at float tolerance -- that reconciliation
# is the core correctness gate (per CLAUDE.md: never assert bit-identity
# between two different code paths, since the decode paths genuinely differ).
#
# Tier: standard (CLI runs against goldens + the positive control).

PGEN_G <- c("--pgenFile=./input/genotype_100markers.pgen",
            "--pvarFile=./input/genotype_100markers.pvar",
            "--psamFile=./input/genotype_100markers.psam")
PLINK_G <- c("--bedFile=./input/genotype_100markers.bed",
             "--bimFile=./input/genotype_100markers.bim",
             "--famFile=./input/genotype_100markers.fam",
             "--AlleleOrder=alt-first")

BIN_MODEL <- function() c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt"))
)
GRP_MODEL <- c(
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
)
SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)

QUANT_MODEL <- function() c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.varianceRatio.txt"))
)
PGEN_CHUNK_G <- c(
  "--pgenFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.pgen",
  "--pvarFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.pvar",
  "--psamFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.psam"
)
PLINK_CHUNK_G <- c(
  "--bedFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bed",
  "--bimFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim",
  "--famFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.fam",
  "--AlleleOrder=alt-first"
)

pgen_run <- function(args, tag = "pgen") {
  out <- file.path(tempfile(paste0("step2_", tag, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(args, paste0("--SAIGEOutputFile=", out)))
  list(status = status, path = out)
}

pnum <- function(x) suppressWarnings(as.numeric(x))
rel_err_vec <- function(a, b) {
  a <- pnum(a); b <- pnum(b)
  abs(a - b) / (abs(a) + abs(b) + .Machine$double.eps)
}

# Reconcile a pgen run's table against its plink-golden mirror at float
# tolerance, merging on `key` columns (per CLAUDE.md: never assert bit-identity
# between two different code paths -- the decode paths genuinely differ).
# Region-test tables carry a pooled Cauchy summary row per region whose
# per-set columns (BETA_Burden, MAC, ...) are the literal string "NA" -- both
# sides must agree on WHICH rows are "NA" (checked below), then the tolerance
# comparison runs only over the rows where both sides report a number.
expect_reconciles_with_plink <- function(pgen_path, plink_golden, key, cols, tol = 1e-8) {
  pgen <- read_saige_table(pgen_path)
  plink <- read_saige_table(plink_golden)
  expect_identical(nrow(pgen), nrow(plink))
  m <- merge(as.data.frame(pgen), as.data.frame(plink), by = key,
            suffixes = c(".pgen", ".plink"))
  expect_identical(nrow(m), nrow(pgen))
  for (col in cols) {
    a <- m[[paste0(col, ".pgen")]]; b <- m[[paste0(col, ".plink")]]
    if (is.null(a) || is.null(b)) next
    expect_identical(a == "NA", b == "NA", label = paste("NA pattern:", col))
    keep <- a != "NA" & b != "NA"
    if (!any(keep)) next
    err <- max(rel_err_vec(a[keep], b[keep]), 0)
    expect_lt(err, tol)
  }
}

# ---- golden CLI runs: marker / group fastTest / group conditional ---------

test_that("golden P13: pgen single-variant marker scan", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
                 BIN_MODEL(), "--is_output_moreDetails=TRUE"), "marker")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_marker_binary.txt"), r$path)
})

test_that("golden P14: pgen group test, binary, fastTest", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--is_fastTest=TRUE"), "group")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_fastTest.txt"), r$path)
})

test_that("golden P15: pgen group test conditional, binary", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--condition=rs30,rs79", "--is_fastTest=TRUE"), "groupcond")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_cond_binary.txt"), r$path)
})

# ---- golden CLI runs P16-P24: extended path coverage (Phase 7 follow-up) --
# Mirrors of every DISTINCT code path the PLINK goldens exercise beyond
# P1-P3/P13-P15, except P11 (--AlleleOrder=ref-first is bim-specific; pvar has
# fixed REF/ALT columns, so there is no ref-first/alt-first ambiguity to
# mirror). Each test also reconciles the pgen run against its PLINK golden
# mirror at float tolerance -- the same positive-control logic as P13-P15's
# check below, extended to the paths that had no pgen coverage at all
# (quantitative traits, BURDEN, non-fastTest SKAT-O, chunking, multiweight,
# composed sets, conditional multiweight).

GROUP_KEY <- c("Region", "Group", "max_MAF")
BURDEN_COLS <- c("Pvalue_Burden", "BETA_Burden", "SE_Burden", "MAC",
                 "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")
SKATO_COLS <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "MAC",
                "Number_rare", "Number_ultra_rare")

test_that("golden P16: pgen single-variant conditional, quantitative", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=20", "--LOCO=FALSE",
                 QUANT_MODEL(), "--is_output_moreDetails=TRUE",
                 "--condition=rs13,rs79"), "marker_cond_quant")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_marker_cond_quantitative.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_marker_cond_quantitative.txt"),
    key = "MarkerID",
    cols = c("AF_Allele2", "AC_Allele2", "BETA", "SE", "Tstat", "var",
             "Tstat_cond", "var_cond", "BETA_cond", "SE_cond"))
})

test_that("golden P17: pgen single-variant conditional, binary (LEDGER #40 path)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
                 BIN_MODEL(), "--is_output_moreDetails=TRUE",
                 "--condition=rs13,rs79"), "marker_cond_binary")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_marker_cond_binary.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_marker_cond_binary.txt"),
    key = "MarkerID",
    cols = c("AF_Allele2", "AC_Allele2", "BETA", "SE", "Tstat", "var",
             "Tstat_cond", "var_cond", "BETA_cond", "SE_cond"))
})

test_that("golden P18: pgen group test, BURDEN only (r.corr=1, is_fastTest=FALSE)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--r.corr=1", "--is_fastTest=FALSE"), "group_burden")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_burden.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_group_burden.txt"),
    key = GROUP_KEY, cols = BURDEN_COLS)
})

test_that("golden P19: pgen group test, SKAT-O, is_fastTest=FALSE", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--is_fastTest=FALSE"), "group_nofast")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_nofast.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_group_nofast.txt"),
    key = GROUP_KEY, cols = SKATO_COLS)
})

test_that("golden P20: pgen single-variant, multi-chunk (different, larger fixture)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_CHUNK_G, "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
                 BIN_MODEL(), "--markers_per_chunk=1000"), "marker_chunked")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_marker_chunked.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_marker_chunked.txt"),
    key = "MarkerID", cols = c("AF_Allele2", "AC_Allele2", "BETA", "SE", "Tstat", "var"))
})

test_that("golden P21: pgen group test, SKAT-O, K=2 custom weight sets, is_fastTest=TRUE", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid_multiWeights.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--is_fastTest=TRUE"), "group_multiweight")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_multiweight.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_group_multiweight.txt"),
    key = c(GROUP_KEY, "Weight"), cols = SKATO_COLS)
})

test_that("golden P22: pgen group test, multiweight BURDEN (r.corr=1, is_fastTest=FALSE)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid_multiWeights.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--r.corr=1", "--is_fastTest=FALSE"), "group_multiweight_burden")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_multiweight_burden.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_group_multiweight_burden.txt"),
    key = c(GROUP_KEY, "Weight"), cols = BURDEN_COLS)
})

test_that("golden P23: pgen group test, composed set (beta weights + No_Weights)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 paste0("--weights.beta=", shQuote("1,25;1,1")),
                 "--is_no_weight_in_groupTest=TRUE",
                 "--is_fastTest=TRUE"), "group_composed")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_composed.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_group_composed.txt"),
    key = c(GROUP_KEY, "Weight"), cols = SKATO_COLS)
})

test_that("golden P24: pgen group test, conditional multi-weight", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid_multiWeights.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--condition=rs30,rs79", "--is_fastTest=TRUE"), "group_cond_multiweight")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_cond_multiweight.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_group_cond_multiweight.txt"),
    key = c(GROUP_KEY, "Weight"),
    cols = c(SKATO_COLS, "Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond",
             "BETA_Burden_cond", "SE_Burden_cond"))
})

# ---- positive control: PGEN vs PLINK reconcile over the SAME genotypes ----

test_that("POSITIVE CONTROL: pgen and plink marker scans reconcile at float tolerance", {
  # The core correctness gate (CLAUDE.md: run a positive control before
  # trusting any pgen-specific golden). Hard-call dosages are exactly 0/1/2 in
  # both formats and the pvar/bim REF/ALT labelling agrees, so the score-test
  # statistics should match to float tolerance -- NOT byte-identity, since the
  # decode paths genuinely differ (pgenlib vs the .bed bit-unpacker).
  skip_unless_level("standard")
  skip_if_no_goldens()
  pgen <- read_saige_table(file.path(saige_goldens, "step2", "pgen_marker_binary.txt"))
  plink <- read_saige_table(file.path(saige_goldens, "step2", "plink_marker_binary.txt"))

  expect_identical(sort(pgen$MarkerID), sort(plink$MarkerID))
  expect_identical(nrow(pgen), nrow(plink))

  m <- merge(as.data.frame(pgen), as.data.frame(plink),
            by = "MarkerID", suffixes = c(".pgen", ".plink"))
  expect_identical(nrow(m), nrow(pgen))

  num <- function(x) suppressWarnings(as.numeric(x))
  rel_err_vec <- function(a, b) {
    a <- num(a); b <- num(b)
    abs(a - b) / (abs(a) + abs(b) + .Machine$double.eps)
  }
  for (col in c("AF_Allele2", "AC_Allele2", "BETA", "SE", "Tstat", "var")) {
    err <- max(rel_err_vec(m[[paste0(col, ".pgen")]], m[[paste0(col, ".plink")]]), 0)
    expect_lt(err, 1e-8)
  }
  # p-values are pre-rounded strings (%.6E); compare parsed at a looser but
  # still tight tolerance, since rounding of two independently-derived
  # z-statistics can differ in the last printed digit.
  perr <- max(rel_err_vec(m$p.value.pgen, m$p.value.plink), 0)
  expect_lt(perr, 1e-4)
})

# ---- --nThreads > 1 refusal for PGEN (v1: not fork-safe) -------------------

# ---- live Firth path through the pgen reader -------------------------------
# Mirror of tests/testthat/test-step2-firth.R (the last arm of the
# getMarkerPval cascade: score -> SPA -> ER -> Firth), swapping the PLINK
# genotype flags for the pgen ones. Firth's inputs (decoded dosages, allele
# counts) come from the SAME reader path already exercised by P13/P17/P20
# above; this confirms the cascade's most expensive arm also gets correctly
# decoded genotypes when reached through pgen. As in test-step2-firth.R, do
# NOT raise pCutoffforFirth above 0.05 -- LEDGER #26.

PGEN_FIRTH_ARGS <- function(out, extra = character(0)) c(
  PGEN_G,
  "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
  BIN_MODEL(),
  "--is_output_moreDetails=TRUE",
  extra,
  paste0("--SAIGEOutputFile=", out)
)

pgen_firth_cache <- new.env(parent = emptyenv())
pgen_firth_run <- function(extra = character(0)) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  key <- if (length(extra)) paste(extra, collapse = " ") else "<none>"
  if (!is.null(pgen_firth_cache[[key]])) return(pgen_firth_cache[[key]])
  out <- file.path(tempfile("step2_pgen_firth_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", PGEN_FIRTH_ARGS(out, extra)))
  res <- read_saige_table(out)
  pgen_firth_cache[[key]] <- res
  res
}

test_that("pgen: is_Firth_beta at the default cutoff changes nothing on this fixture", {
  base <- pgen_firth_run()
  firth <- pgen_firth_run("--is_Firth_beta=TRUE")
  expect_identical(base, firth)
  expect_gt(min(pnum(base$p.value)), 0.01)
})

test_that("pgen: pCutoffforFirth selects exactly the markers at or below the cutoff", {
  base <- pgen_firth_run()
  firth <- pgen_firth_run(c("--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"))
  changed <- which(base$BETA != firth$BETA)
  qualifying <- which(pnum(base$p.value) <= 0.05)
  expect_identical(changed, qualifying)
  expect_gt(length(qualifying), 3)
  expect_identical(base$BETA[-qualifying], firth$BETA[-qualifying])
  expect_identical(base$SE[-qualifying], firth$SE[-qualifying])
})

test_that("pgen: Firth touches only the effect-size columns, never the score test", {
  base <- pgen_firth_run()
  firth <- pgen_firth_run(c("--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"))
  for (col in c("MarkerID", "AC_Allele2", "AF_Allele2", "Tstat", "var",
                "p.value", "N_case", "N_ctrl", "N_case_het", "N_ctrl_het")) {
    expect_identical(base[[col]], firth[[col]], label = paste("unchanged:", col))
  }
  expect_false(identical(base$BETA, firth$BETA))
  expect_false(identical(base$SE, firth$SE))
})

test_that("pgen: live Firth estimates reconcile with the plink Firth mirror at float tolerance", {
  # Positive control across readers: same fixture genotypes, same penalised
  # fit, decoded through a different reader -- BETA/SE for the qualifying
  # markers should agree with the PLINK live-Firth run to float tolerance.
  pgen_base <- pgen_firth_run()
  pgen_firth <- pgen_firth_run(c("--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"))
  out <- file.path(tempfile("step2_plinkfirth_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  plink_args <- c(PLINK_G, "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
                  BIN_MODEL(), "--is_output_moreDetails=TRUE",
                  "--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05",
                  paste0("--SAIGEOutputFile=", out))
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", plink_args))
  plink_firth <- read_saige_table(out)

  m <- merge(as.data.frame(pgen_firth), as.data.frame(plink_firth),
            by = "MarkerID", suffixes = c(".pgen", ".plink"))
  expect_identical(nrow(m), nrow(pgen_firth))
  qualifying <- which(pnum(pgen_base$p.value) <= 0.05)
  expect_gt(length(qualifying), 3)
  qids <- pgen_base$MarkerID[qualifying]
  mq <- m[m$MarkerID %in% qids, ]
  for (col in c("BETA", "SE")) {
    err <- max(rel_err_vec(mq[[paste0(col, ".pgen")]], mq[[paste0(col, ".plink")]]), 0)
    expect_lt(err, 1e-6)
  }
})

test_that("golden P26: pgen group test, BURDEN conditional (r.corr=1, rs30,rs79)", {
  # PGEN mirror of P25 (the LEDGER #41 coverage golden).
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- pgen_run(c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 GRP_MODEL, SG,
                 "--groupFile=./input/group_new_snpid.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01",
                 "--r.corr=1", "--is_fastTest=FALSE",
                 "--condition=rs30,rs79"), "group_cond_burden")
  expect_cli_ok(r$status)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_group_cond_burden.txt"), r$path)
  expect_reconciles_with_plink(
    r$path, file.path(saige_goldens, "step2", "plink_group_cond_burden.txt"),
    key = GROUP_KEY,
    cols = c(BURDEN_COLS, "Pvalue_Burden_c", "Beta_Burden_c", "seBeta_Burden_c"))
})

test_that("--nThreads > 1 on a pgen input degrades to one thread, loudly", {
  # INVERTED BY S7, deliberately. This used to assert a hard ERROR, because
  # --nThreads forked and PgenReader::isForkSafe() is always false. S7 routed
  # --nThreads to the in-process path, where the question is whether the reader
  # can CLONE -- PgenReader::clone() still throws (v1 decision unchanged) -- and
  # the agreed behaviour became: degrade to one thread rather than fail.
  #
  # The reasoning, so this is not re-inverted on taste: erroring would break
  # anyone scripting a single --nThreads across a mix of PLINK and PGEN inputs,
  # and degrading is the posture group tests already take. It is NOT the silent
  # no-op of LEDGER #30 -- the reason is both printed and raised as an R
  # warning -- which is why the message is asserted here and not merely the
  # exit status. PGEN gains real parallelism only when PgenReader::clone() is
  # implemented; until then this is the honest answer.
  skip_unless_level("standard")
  skip_if_no_goldens()
  log <- tempfile(fileext = ".log")
  out <- file.path(tempfile("step2_pgen_nthreads_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(PGEN_G, "--chrom=1", "--minMAF=0", "--minMAC=1",
                            "--LOCO=FALSE", BIN_MODEL(),
                            "--is_output_moreDetails=TRUE", "--nThreads=4",
                            paste0("--SAIGEOutputFile=", out)),
                          log = log)
  expect_cli_ok(status)
  txt <- readLines(log, warn = FALSE)
  expect_true(any(grepl("cannot be cloned", txt)))
  expect_true(any(grepl("ONE thread", txt)))
  # It really did fall back -- no threaded-scan line was printed ...
  expect_false(any(grepl("in-process threads", txt)))
  # ... and the degraded run still produced the right answer.
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_marker_binary.txt"), out)
})
