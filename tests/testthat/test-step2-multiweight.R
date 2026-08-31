# Multi-weight (K > 1) region tests, end-to-end through the step-2 CLI.
#
# Two runs of the P2 group-fastTest configuration (same model, sparse GRM,
# annotations and max-MAF cutoffs as test-step2-goldens.R's "group fastTest"
# case), differing only in the group file:
#   - group_new_snpid_multiWeights.txt  (K = 2: weight:CADD, weight,REVEL)
#   - group_new_snpid_withWeights.txt   (K = 1: one bare weight line)
#
# What is pinned here:
#   K = 2 -> rows are (anno x max_MAF x weight set) with the weight innermost,
#            a Weight column right after max_MAF, additive Pvalue_ACATV /
#            Pvalue_ACATO columns (Pvalue stays SKAT-O), and exactly ONE
#            pooled Cauchy row per region.
#   K = 1 -> the legacy schema, byte-for-byte the same columns as the
#            committed unweighted golden: no Weight/ACATV/ACATO.
#   Either way a real marker is score-tested ONCE, no matter how many weight
#   sets are in play; only the collapsed ultra-rare pseudo-markers are per set.
#
# Tier: standard (runs the CLI), like the golden tests.

MW_ARGS <- function(groupFile, out) c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first",
  "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
  "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
  "--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
  paste0("--groupFile=./input/", groupFile),
  "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
  "--maxMAF_in_groupTest=0.0001,0.001,0.01",
  "--is_fastTest=TRUE",
  paste0("--SAIGEOutputFile=", out)
)

# One CLI run per (group file, extra args) pair, cached across the test_that
# blocks below.
mw_cache <- new.env(parent = emptyenv())
mw_run <- function(groupFile, extra = character(0)) {
  skip_unless_level("standard")
  key <- paste(c(groupFile, extra), collapse = " ")
  if (!is.null(mw_cache[[key]])) return(mw_cache[[key]])
  out <- file.path(tempfile("step2_mw_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R", c(MW_ARGS(groupFile, out), extra))
  expect_cli_ok(status)
  # BURDEN-only runs do not emit the single-variant file
  single_path <- paste0(out, ".singleAssoc.txt")
  res <- list(region = read_saige_table(out),
              single = if (file.exists(single_path)) read_saige_table(single_path) else NULL)
  mw_cache[[key]] <- res
  res
}

multi_run <- function() mw_run("group_new_snpid_multiWeights.txt")
single_run <- function() mw_run("group_new_snpid_withWeights.txt")
# BURDEN-only (--r.corr=1): the region file is written from C++
# (openOutfile/writeOutfile_BURDEN), not the R SKAT-O post-processing
burden_multi_run <- function() mw_run("group_new_snpid_multiWeights.txt", "--r.corr=1")
burden_single_run <- function() mw_run("group_new_snpid_withWeights.txt", "--r.corr=1")
# C8 set composition (D10): explicit beta sets from a weights.beta pair list
# on an UNWEIGHTED group file (the semicolon needs shell quoting), the
# appended No_Weights set on a weighted file, and the forced Weight column at
# a single set (SKAT-O and BURDEN writers)
beta_run <- function() mw_run("group_new_snpid.txt",
                              paste0("--weights.beta=", shQuote("1,25;1,1")))
noweights_run <- function() mw_run("group_new_snpid_multiWeights.txt",
                                   "--is_no_weight_in_groupTest=TRUE")
forced_run <- function() mw_run("group_new_snpid_withWeights.txt",
                                "--is_output_weight_column=TRUE")
burden_forced_run <- function() mw_run("group_new_snpid_withWeights.txt",
                                       c("--r.corr=1", "--is_output_weight_column=TRUE"))

WSETS <- c("CADD", "REVEL")
pnum <- function(x) suppressWarnings(as.numeric(x))
acat_p <- function(p) 0.5 - atan(mean(tan((0.5 - p) * pi))) / pi

test_that("K>1 emits the multi-weight schema: Weight after max_MAF, ACATV/ACATO added", {
  d <- multi_run()$region
  expect_identical(colnames(d)[1:9],
                   c("Region", "Group", "max_MAF", "Weight", "Pvalue",
                     "Pvalue_Burden", "Pvalue_SKAT", "Pvalue_ACATV", "Pvalue_ACATO"))
  # Pvalue keeps its SKAT-O meaning: same sandwich bounds the unweighted
  # goldens are held to in test-invariants-goldens.R (>= the smaller component
  # p-value, and within 2x of it after the rho-grid multiplicity correction)
  r <- d[d$Group != "Cauchy", ]
  lo <- pmin(pnum(r$Pvalue_Burden), pnum(r$Pvalue_SKAT))
  expect_true(all(pnum(r$Pvalue) >= lo * (1 - 1e-12)))
  expect_true(all(pnum(r$Pvalue) <= 2 * lo * (1 + 1e-12)))
})

test_that("K>1 row layout: weight-innermost, K rows per (anno, max_MAF), one Cauchy row", {
  d <- multi_run()$region
  s <- single_run()$region
  for (region in unique(d$Region)) {
    rows <- d[d$Region == region & d$Group != "Cauchy", ]
    srows <- s[s$Region == region & s$Group != "Cauchy", ]
    # exactly K rows per (annotation, max_MAF) cell, and the same set of cells
    # as the single-weight run
    expect_identical(nrow(rows), 2L * nrow(srows))
    expect_identical(unique(paste(rows$Group, rows$max_MAF)),
                     unique(paste(srows$Group, srows$max_MAF)))
    expect_identical(nrow(rows), length(unique(paste(rows$Group, rows$max_MAF))) * 2L)
    # weight fastest: the Weight column cycles CADD, REVEL within each cell
    expect_identical(rows$Weight, rep(WSETS, nrow(rows) / 2L))
    expect_identical(rows$Group, rep(srows$Group, each = 2L))
    expect_identical(rows$max_MAF, rep(srows$max_MAF, each = 2L))
    # exactly one pooled Cauchy row per region, named and un-labelled by weight
    cauchy <- d[d$Region == region & d$Group == "Cauchy", ]
    expect_identical(nrow(cauchy), 1L)
    expect_identical(cauchy$Weight, "NA")
    expect_identical(cauchy$max_MAF, "NA")
  }
  expect_identical(sum(d$Group == "Cauchy"), length(unique(d$Region)))
})

test_that("K>1 Pvalue_ACATO reconciles with CCT(ACATV, Burden, SKAT) of its own row", {
  d <- multi_run()$region
  r <- d[d$Group != "Cauchy", ]
  for (i in seq_len(nrow(r))) {
    trio <- pnum(c(r$Pvalue_ACATV[i], r$Pvalue_Burden[i], r$Pvalue_SKAT[i]))
    expect_lt(rel_err(pnum(r$Pvalue_ACATO[i]), CCT(trio)), 1e-12)
  }
})

test_that("K>1 pooled Cauchy row is ACAT over every row of the region, per p-column", {
  d <- multi_run()$region
  pcols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "Pvalue_ACATV", "Pvalue_ACATO")
  for (region in unique(d$Region)) {
    rows <- d[d$Region == region & d$Group != "Cauchy", ]
    cauchy <- d[d$Region == region & d$Group == "Cauchy", ]
    for (col in pcols) {
      expect_lt(rel_err(pnum(cauchy[[col]]), acat_p(pnum(rows[[col]]))), 1e-12)
    }
  }
})

test_that("K>1 p-values are all in (0, 1]", {
  d <- multi_run()$region
  for (col in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                "Pvalue_ACATV", "Pvalue_ACATO")) {
    p <- pnum(d[[col]])
    expect_false(anyNA(p))
    expect_true(all(p > 0 & p <= 1))
  }
})

test_that("K=1 keeps the legacy schema: no Weight/ACATV/ACATO columns", {
  skip_if_no_goldens()
  s <- single_run()$region
  expect_false(any(c("Weight", "Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(s)))
  # same columns as the committed unweighted golden of this configuration
  g <- read_saige_table(file.path(saige_goldens, "step2", "plink_group_fastTest.txt"))
  expect_identical(colnames(s), colnames(g))
  expect_identical(s$Region, g$Region)
  expect_identical(s$Group, g$Group)
  expect_identical(s$max_MAF, g$max_MAF)
})

# Legacy BURDEN header written by openOutfile for a binary trait, no
# conditioning (Region.txt is C++-written when --r.corr=1)
BURDEN_HDR_K1 <- c("Region", "Group", "max_MAF", "Pvalue_Burden", "BETA_Burden",
                   "SE_Burden", "MAC", "MAC_case", "MAC_control",
                   "Number_rare", "Number_ultra_rare")

test_that("BURDEN K>1: Weight column after max_MAF, every set named, weight-innermost rows", {
  d <- burden_multi_run()$region
  # additive schema: the legacy columns with Weight inserted right after max_MAF
  expect_identical(colnames(d), append(BURDEN_HDR_K1, "Weight", after = 3L))
  r <- d[d$Group != "Cauchy", ]
  s <- burden_single_run()$region
  sr <- s[s$Group != "Cauchy", ]
  # every set name appears, and rows cycle weight-innermost within each
  # (annotation, max_MAF) cell
  expect_true(all(WSETS %in% r$Weight))
  expect_identical(r$Weight, rep(WSETS, nrow(r) / 2L))
  expect_identical(nrow(r), 2L * nrow(sr))
  expect_identical(r$Group, rep(sr$Group, each = 2L))
  expect_identical(r$max_MAF, rep(sr$max_MAF, each = 2L))
  # exactly one pooled Cauchy row per region, not tied to any weight set
  cauchy <- d[d$Group == "Cauchy", ]
  expect_identical(nrow(cauchy), length(unique(d$Region)))
  expect_true(all(cauchy$Weight == "NA"))
  expect_true(all(cauchy$max_MAF == "NA"))
  # sanity on the numbers: burden p-values in (0, 1]
  p <- pnum(d$Pvalue_Burden)
  expect_false(anyNA(p))
  expect_true(all(p > 0 & p <= 1))
})

test_that("BURDEN K=1 keeps the legacy header exactly", {
  s <- burden_single_run()$region
  expect_identical(colnames(s), BURDEN_HDR_K1)
  expect_false(any(c("Weight", "Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(s)))
})

# ---- C8: weights.beta pair list -> explicit beta sets ----------------------

BETA_SETS <- c("Beta_1_25", "Beta_1_1")

test_that("unweighted file + weights.beta pair list: one beta set per pair, in given order", {
  d <- beta_run()$region
  expect_identical(colnames(d)[1:9],
                   c("Region", "Group", "max_MAF", "Weight", "Pvalue",
                     "Pvalue_Burden", "Pvalue_SKAT", "Pvalue_ACATV", "Pvalue_ACATO"))
  r <- d[d$Group != "Cauchy", ]
  # weight-innermost rows cycling the pair-derived set names, K rows per cell
  expect_identical(r$Weight, rep(BETA_SETS, nrow(r) / 2L))
  expect_identical(nrow(r), length(unique(paste(r$Region, r$Group, r$max_MAF))) * 2L)
  # one pooled Cauchy row per region
  expect_identical(sum(d$Group == "Cauchy"), length(unique(d$Region)))
  # positive control: the two beta parameterizations must actually differ
  p25 <- pnum(r$Pvalue[r$Weight == "Beta_1_25"])
  p11 <- pnum(r$Pvalue[r$Weight == "Beta_1_1"])
  expect_true(any(abs(p25 - p11) > 0))
  # per-row ACAT-O still reconciles on the beta-set rows
  for (i in seq_len(nrow(r))) {
    trio <- pnum(c(r$Pvalue_ACATV[i], r$Pvalue_Burden[i], r$Pvalue_SKAT[i]))
    expect_lt(rel_err(pnum(r$Pvalue_ACATO[i]), CCT(trio)), 1e-12)
  }
  # UR pseudo-markers are per set, suffixed with the beta set names
  ur <- beta_run()$single
  ur <- ur[ur$CHR == "UR", ]
  expect_gt(nrow(ur), 0)
  expect_identical(sub("^.*:", "", ur$MarkerID), rep(BETA_SETS, nrow(ur) / 2L))
})

# ---- C8: is_no_weight_in_groupTest appends No_Weights to a weighted file ----

test_that("weighted file + is_no_weight_in_groupTest: No_Weights appended after CADD/REVEL", {
  d <- noweights_run()$region
  m <- multi_run()$region
  r <- d[d$Group != "Cauchy", ]
  rm_ <- m[m$Group != "Cauchy", ]
  sets3 <- c(WSETS, "No_Weights")
  # the No_Weights set is appended LAST: rows cycle CADD, REVEL, No_Weights
  expect_identical(r$Weight, rep(sets3, nrow(r) / 3L))
  # same (anno, max_MAF) cells as the 2-set run, one extra row per cell
  expect_identical(2L * nrow(r), 3L * nrow(rm_))
  expect_identical(r$Group, rep(rm_$Group[seq(1, nrow(rm_), by = 2L)], each = 3L))
  # still exactly one pooled Cauchy row per region
  expect_identical(sum(d$Group == "Cauchy"), length(unique(d$Region)))
  # sanity: the all-ones set is not a copy of a custom set
  pn <- pnum(r$Pvalue[r$Weight == "No_Weights"])
  pc <- pnum(r$Pvalue[r$Weight == "CADD"])
  expect_true(any(abs(pn - pc) > 0))
  # its UR pseudo-markers carry the No_Weights suffix
  ur <- noweights_run()$single
  ur <- ur[ur$CHR == "UR", ]
  expect_true(any(grepl(":No_Weights$", ur$MarkerID)))
})

# ---- C8: --is_output_weight_column forces the schema at K = 1 ---------------

test_that("forced Weight column at K=1 (SKAT-O): set name labels rows; excising the new columns reproduces the plain run", {
  f <- forced_run()$region
  s <- single_run()$region
  newcols <- c("Weight", "Pvalue_ACATV", "Pvalue_ACATO")
  expect_true(all(newcols %in% colnames(f)))
  expect_identical(colnames(f)[1:5],
                   c("Region", "Group", "max_MAF", "Weight", "Pvalue"))
  # the single set keeps its group-file name (bare 'weight' line -> "WEIGHT")
  expect_identical(unique(f$Weight[f$Group != "Cauchy"]), "WEIGHT")
  expect_identical(unique(f$Weight[f$Group == "Cauchy"]), "NA")
  # excising the additive columns reproduces the un-flagged run byte-for-byte
  f2 <- f[, setdiff(colnames(f), newcols)]
  rownames(f2) <- NULL
  expect_identical(colnames(f2), colnames(s))
  expect_identical(f2, s)
  # the single-variant output is schema-independent
  expect_identical(forced_run()$single, single_run()$single)
})

test_that("forced Weight column at K=1 (BURDEN writer): header and rows gain only the Weight column", {
  f <- burden_forced_run()$region
  s <- burden_single_run()$region
  expect_identical(colnames(f), append(BURDEN_HDR_K1, "Weight", after = 3L))
  expect_identical(unique(f$Weight[f$Group != "Cauchy"]), "WEIGHT")
  expect_identical(unique(f$Weight[f$Group == "Cauchy"]), "NA")
  f2 <- f[, setdiff(colnames(f), "Weight")]
  rownames(f2) <- NULL
  expect_identical(colnames(f2), colnames(s))
  expect_identical(f2, s)
})

# ---- C8: weights.beta pair-list parser (fast tier: pure R, no CLI run) ------

test_that("parse_weights_beta accepts the legacy and pair-list forms", {
  pw <- SAIGE:::parse_weights_beta
  expect_identical(pw("1,25"), matrix(c(1, 25), nrow = 1, ncol = 2, byrow = TRUE))
  expect_identical(pw("1,25;1,1"),
                   matrix(c(1, 25, 1, 1), nrow = 2, ncol = 2, byrow = TRUE))
  expect_identical(pw("0.5,0.5;1,25;2,3"),
                   matrix(c(0.5, 0.5, 1, 25, 2, 3), nrow = 3, ncol = 2, byrow = TRUE))
  expect_identical(pw(c(1, 25)), matrix(c(1, 25), nrow = 1, ncol = 2))
  m <- matrix(c(1, 25, 0.5, 0.5), nrow = 2, ncol = 2, byrow = TRUE)
  expect_identical(pw(m), m)
})

test_that("parse_weights_beta rejects malformed pairs, naming the offender", {
  pw <- SAIGE:::parse_weights_beta
  expect_error(pw("1,25;2"), "offending pair: '2'", fixed = TRUE)
  expect_error(pw("1,2,3"), "offending pair: '1,2,3'", fixed = TRUE)
  expect_error(pw("a,b"), "offending pair: 'a,b'", fixed = TRUE)
  expect_error(pw("1,25;;1,1"), "offending pair")
  expect_error(pw(""), "weights.beta")
  expect_error(pw(c(1, 25, 3)), "exactly 2 elements")
  expect_error(pw(TRUE), "weights.beta")
})

test_that("parse_weights_beta rejects non-positive or non-finite parameters", {
  pw <- SAIGE:::parse_weights_beta
  expect_error(pw("0,25"), "offending pair 1: '0,25'", fixed = TRUE)
  expect_error(pw("1,-1"), "finite and > 0")
  expect_error(pw("1,25;1,Inf"), "offending pair 2", fixed = TRUE)
  expect_error(pw(c(-1, 25)), "finite and > 0")
})

test_that("real markers are tested once regardless of K; only UR pseudo-markers are per set", {
  dm <- multi_run()$single
  ds <- single_run()$single
  is_ur <- function(d) d$CHR == "UR"
  # real markers: identical rows in both runs, each appearing exactly once
  drop_rn <- function(d) { rownames(d) <- NULL; d }
  real_m <- drop_rn(dm[!is_ur(dm), ])
  real_s <- drop_rn(ds[!is_ur(ds), ])
  expect_gt(nrow(real_m), 0)
  expect_false(anyDuplicated(real_m$MarkerID) > 0)
  expect_identical(real_m, real_s)
  # ultra-rare collapsed pseudo-markers: one per (anno, max_MAF, weight set),
  # tagged with the weight-set name; twice as many as the single-set run
  ur_m <- dm[is_ur(dm), ]
  ur_s <- ds[is_ur(ds), ]
  expect_identical(nrow(ur_m), 2L * nrow(ur_s))
  expect_false(anyDuplicated(ur_m$MarkerID) > 0)
  suffix <- sub("^.*:", "", ur_m$MarkerID)
  expect_identical(suffix, rep(WSETS, nrow(ur_m) / 2L))
  # dropping the suffix recovers exactly the single-set pseudo-marker names
  expect_identical(unique(sub(":[^:]+$", "", ur_m$MarkerID)), ur_s$MarkerID)
})

test_that("negative weights in a K=1 file are rejected under the multi-weight schema", {
  skip_unless_level("standard")
  # checkGroupFile's parse-time sign check is K>1-file-scoped, so a K=1
  # weighted file with a negative weight parses fine -- but the multi-weight
  # SCHEMA (forced here via is_output_weight_column; composition via
  # is_no_weight_in_groupTest or multiple weights.beta pairs is the same
  # path) computes ACAT-V from those weights, and CCT stop()s mid-run on a
  # negative. The per-region guard in SAIGE.Region must fail fast instead,
  # with a message naming the region.
  src <- file.path(saige_extdata, "input", "group_new_snpid_withWeights.txt")
  lines <- readLines(src)
  wIdx <- which(vapply(strsplit(lines, "[ \t]+"),
                       function(p) length(p) >= 3 && p[2] == "weight",
                       logical(1)))[1]
  parts <- strsplit(lines[wIdx], "[ \t]+")[[1]]
  parts[3] <- paste0("-", parts[3])
  lines[wIdx] <- paste(parts, collapse = " ")
  negfile <- tempfile(fileext = ".txt")
  writeLines(lines, negfile)
  out <- file.path(tempfile("step2_mwneg_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- c(MW_ARGS("group_new_snpid.txt", out), "--is_output_weight_column=TRUE")
  args[grep("^--groupFile=", args)] <- paste0("--groupFile=", negfile)
  status <- run_saige_cli("step2_SPAtests.R", args)
  expect_false(identical(as.integer(status), 0L))
  log <- readLines(attr(status, "log"), warn = FALSE)
  expect_true(any(grepl("negative marker weights", log)))
})

# ---- guards that must fail fast (C8/C9) ------------------------------------
# Runs expected to exit non-zero; returns the CLI log so the exact message can
# be asserted. The messages are user-facing contract, not incidental text.
mw_fail_log <- function(groupFile, extra = character(0)) {
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_mwfail_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R", c(MW_ARGS(groupFile, out), extra))
  expect_false(identical(as.integer(status), 0L))
  readLines(attr(status, "log"), warn = FALSE)
}

test_that("conditional analysis with K>1 is refused, naming the option and the way out", {
  # The per-set conditional adjustment (design D7) is not implemented; until it
  # lands the combination must stop rather than silently condition every set on
  # the same single-weight vector.
  log <- mw_fail_log("group_new_snpid_multiWeights.txt", "--condition=rs30,rs79")
  expect_true(any(grepl("not yet supported together with conditional analysis", log)))
  expect_true(any(grepl("--condition", log, fixed = TRUE)))
  expect_true(any(grepl("single weight set", log)))
})

test_that("conditional analysis with a forced Weight column at K=1 is refused too", {
  # is_output_weight_column switches a single-set run onto the multi-weight
  # schema, which the conditional writer cannot emit either.
  log <- mw_fail_log("group_new_snpid_withWeights.txt",
                     c("--condition=rs30,rs79", "--is_output_weight_column=TRUE"))
  expect_true(any(grepl("not yet supported together with conditional analysis", log)))
  expect_true(any(grepl("is_output_weight_column", log, fixed = TRUE)))
})

test_that("a single-weight conditional run is still accepted", {
  # Negative control for the two guards above: the same conditional run
  # without any multi-weight schema must succeed, so the stop() is attributable
  # to the weight configuration and not to --condition itself.
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_mwcond_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(MW_ARGS("group_new_snpid_withWeights.txt", out),
                            "--condition=rs30,rs79"))
  expect_cli_ok(status)
  res <- read_saige_table(out)
  expect_false("Weight" %in% colnames(res))
  # the SKAT-O region writer names the conditional columns with a _cond suffix
  expect_true(all(c("Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond")
                  %in% colnames(res)))
})

test_that("a region run cannot be restarted, so the append schema guard never fires", {
  # SAIGE.Region writes "nEachChunk = <numberRegionsInChunk>" into the index,
  # but numberRegionsInChunk is still 0 when nEachChunk is captured
  # (SAIGE_SPATest_Region.R:881,890) while checkOutputFile is called with a
  # hardcoded 1 (:631). The two never match, so --is_overwrite_output=FALSE
  # always fails validation on an untouched index -- even restarting the very
  # same command. Consequence: the C6 append-restart schema guard at :770
  # (multi-weight columns vs legacy header) is unreachable through the CLI.
  # LEDGER #25: recorded, not fixed.
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_mwrestart_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- MW_ARGS("group_new_snpid_withWeights.txt", out)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
  idx <- readLines(paste0(out, ".index"), warn = FALSE)
  expect_true(any(grepl("^nEachChunk = 0$", idx)))   # writer's value

  status <- run_saige_cli("step2_SPAtests.R", c(args, "--is_overwrite_output=FALSE"))
  expect_false(identical(as.integer(status), 0L))
  log <- readLines(attr(status, "log"), warn = FALSE)
  expect_true(any(grepl("is not as expected", log)))
})

# ---- the implicit_beta set (LEDGER #19) -------------------------------------
# An UNWEIGHTED group file with --is_output_weight_column=TRUE composes no set
# at all, so region_postprocess falls to kind "implicit_beta": the literal
# legacy dbeta(MAFVec, 1, 25). This is the only way to observe that arm, and
# it is where the R/C++ weight asymmetry becomes visible in the output.
implicit_beta_run <- function(beta) {
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_implbeta_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- c(MW_ARGS("group_new_snpid.txt", out),
            "--is_output_weight_column=TRUE",
            paste0("--weights.beta=", beta))
  # single anno/max_MAF keeps the run small; the asymmetry is per row
  args[grep("^--annotation_in_groupTest=", args)] <- "--annotation_in_groupTest=lof"
  args[grep("^--maxMAF_in_groupTest=", args)] <- "--maxMAF_in_groupTest=0.01"
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
  read_saige_table(out)
}

test_that("an unweighted file with a forced Weight column labels its single set Beta", {
  res <- implicit_beta_run("1,25")
  expect_true("Weight" %in% colnames(res))
  expect_identical(which(colnames(res) == "Weight"), 4L)   # right after max_MAF
  lab <- unique(res$Weight[res$Group != "Cauchy"])
  expect_identical(lab, "Beta")
  expect_identical(unique(res$Weight[res$Group == "Cauchy"]), "NA")
})

test_that("weights.beta moves the C++ statistics but not the R-side ACAT-V (LEDGER #19)", {
  # region_postprocess hardcodes dbeta(MAF, 1, 25) for the implicit_beta set
  # while C++ honours g_cfg.weights_beta, so changing --weights.beta shifts
  # every C++-computed column and leaves Pvalue_ACATV exactly where it was.
  # Pvalue_ACATO moves because it is CCT(ACATV, Burden, SKAT) and two of its
  # three inputs moved. Preserved bug-for-bug: this test documents the
  # asymmetry, it does not endorse it.
  a <- implicit_beta_run("1,25")
  b <- implicit_beta_run("1,1")
  expect_identical(colnames(a), colnames(b))

  # frozen in R, therefore identical across the two runs
  expect_identical(a$Pvalue_ACATV, b$Pvalue_ACATV)
  # computed in C++ from the CLI value, therefore different
  for (cn in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")) {
    expect_false(identical(a[[cn]], b[[cn]]),
                 label = paste0("weights.beta changes ", cn))
  }
  # allele accounting is weight-independent
  for (cn in c("MAC", "Number_rare", "Number_ultra_rare")) {
    expect_identical(a[[cn]], b[[cn]])
  }
})

# ---- C10 invariants --------------------------------------------------------

test_that("BURDEN + is_fastTest=TRUE emits only the last region (LEDGER #27)", {
  # Silent data loss, independent of the weight dimension: every earlier region
  # is analysed (the log announces it) but never written. RECORDED, NOT FIXED
  # -- this test asserts the defect so the fix phase trips it. is_fastTest=FALSE
  # is the correct-output control immediately below.
  skip_unless_level("standard")
  for (gf in c("group_new_snpid.txt", "group_new_snpid_multiWeights.txt")) {
    out <- file.path(tempfile("step2_burdenfast_"), "region.txt")
    dir.create(dirname(out), recursive = TRUE)
    args <- MW_ARGS(gf, out)
    args <- c(args[!grepl("^--is_fastTest=", args)], "--r.corr=1", "--is_fastTest=TRUE")
    expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
    d <- read_saige_table(out)
    expect_identical(unique(d$Region), "GENE2")   # GENE1 analysed, never written
  }
})

test_that("BURDEN + is_fastTest=FALSE keeps every region", {
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_burdenslow_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- MW_ARGS("group_new_snpid_multiWeights.txt", out)
  args <- c(args[!grepl("^--is_fastTest=", args)], "--r.corr=1", "--is_fastTest=FALSE")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
  d <- read_saige_table(out)
  expect_identical(sort(unique(d$Region)), c("GENE1", "GENE2"))
})

test_that("burden linearity: scaling every weight by 2 is exact without the UR collapse", {
  # A uniform rescaling of user-supplied weights is a reparameterisation of the
  # burden test: the score scales by c, its variance by c^2, so the p-value is
  # unchanged and BETA/SE both scale by 1/c. Checked with the ultra-rare
  # collapse DISABLED, where the identity is exact to machine precision.
  skip_unless_level("standard")
  doubled <- local({
    src <- file.path(saige_extdata, "input", "group_new_snpid_multiWeights.txt")
    lines <- readLines(src)
    out <- tempfile(fileext = ".txt")
    writeLines(vapply(lines, function(l) {
      p <- strsplit(l, "[ \t]+")[[1]]
      if (length(p) >= 3 && grepl("^weight", p[2])) {
        p[-(1:2)] <- format(as.numeric(p[-(1:2)]) * 2, trim = TRUE, scientific = FALSE)
        paste(p, collapse = " ")
      } else l
    }, "", USE.NAMES = FALSE), out)
    out
  })
  run_one <- function(groupFile) {
    out <- file.path(tempfile("step2_lin_"), "region.txt")
    dir.create(dirname(out), recursive = TRUE)
    args <- MW_ARGS("group_new_snpid_multiWeights.txt", out)
    args[grep("^--groupFile=", args)] <- paste0("--groupFile=", groupFile)
    args <- c(args, "--MACCutoff_to_CollapseUltraRare=0")
    expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
    read_saige_table(out)
  }
  a <- run_one(file.path(saige_extdata, "input", "group_new_snpid_multiWeights.txt"))
  b <- run_one(doubled)
  ra <- a[a$Group != "Cauchy", ]; rb <- b[b$Group != "Cauchy", ]
  expect_identical(nrow(ra), nrow(rb))
  expect_true(all(pnum(ra$Number_ultra_rare) == 0))    # collapse really is off
  # scale-invariant: the p-values are unchanged, to the printed precision
  expect_identical(ra$Pvalue_Burden, rb$Pvalue_Burden)
  expect_identical(ra$Pvalue_SKAT, rb$Pvalue_SKAT)
  # equivariant: BETA and SE both halve
  expect_lt(max(abs(pnum(rb$BETA_Burden) / pnum(ra$BETA_Burden) - 0.5)), 1e-9)
  expect_lt(max(abs(pnum(rb$SE_Burden) / pnum(ra$SE_Burden) - 0.5)), 1e-9)
})

test_that("the ultra-rare collapse makes custom weights only approximately scale-free", {
  # With the collapse ON the same doubling stops being exact: the collapsed
  # pseudo-genotype is max(w_r*G), so its apparent allele count doubles too
  # (AC 38.40 -> 76.80 in the singleAssoc output) while its weight in
  # AnnoWeights stays 1, and the MAC-dependent machinery downstream therefore
  # sees a different marker. So a uniform rescaling of a group file -- which a
  # user would reasonably expect to be a no-op -- shifts the answer slightly.
  # Bounded, not asserted equal: a real weight-handling regression would move
  # these by orders of magnitude, while this residual is ~1e-2.
  skip_unless_level("standard")
  doubled <- local({
    lines <- readLines(file.path(saige_extdata, "input",
                                 "group_new_snpid_multiWeights.txt"))
    out <- tempfile(fileext = ".txt")
    writeLines(vapply(lines, function(l) {
      p <- strsplit(l, "[ \t]+")[[1]]
      if (length(p) >= 3 && grepl("^weight", p[2])) {
        p[-(1:2)] <- format(as.numeric(p[-(1:2)]) * 2, trim = TRUE, scientific = FALSE)
        paste(p, collapse = " ")
      } else l
    }, "", USE.NAMES = FALSE), out)
    out
  })
  out <- file.path(tempfile("step2_lin_on_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- MW_ARGS("group_new_snpid_multiWeights.txt", out)
  args[grep("^--groupFile=", args)] <- paste0("--groupFile=", doubled)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))

  a <- multi_run()$region                      # collapse on, original weights
  b <- read_saige_table(out)                   # collapse on, doubled weights
  ra <- a[a$Group != "Cauchy", ]; rb <- b[b$Group != "Cauchy", ]
  expect_true(all(pnum(ra$Number_ultra_rare) > 0))   # collapse really is active
  rel <- function(x, y) max(abs(x - y) / (abs(x) + abs(y) + 1e-300))
  # not exact any more ...
  expect_false(identical(ra$Pvalue_Burden, rb$Pvalue_Burden))
  # ... but still close, and BETA/SE still halve to within a couple of percent
  expect_lt(rel(pnum(ra$Pvalue_Burden), pnum(rb$Pvalue_Burden)), 0.05)
  expect_lt(max(abs(pnum(rb$BETA_Burden) / pnum(ra$BETA_Burden) - 0.5)), 0.05)
  expect_lt(max(abs(pnum(rb$SE_Burden) / pnum(ra$SE_Burden) - 0.5)), 0.05)
})

test_that("No_Weights rows reproduce a plain all-ones run exactly", {
  # The appended all-ones set must be the same analysis as running the
  # unweighted group file in no-weight mode: same cells, same numbers, column
  # for column. This is what pins design D2's "noweight = all ones with an
  # unweighted UR collapse" against the C++ side.
  skip_unless_level("standard")
  plain_out <- file.path(tempfile("step2_nwplain_"), "region.txt")
  dir.create(dirname(plain_out), recursive = TRUE)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
    c(MW_ARGS("group_new_snpid.txt", plain_out), "--is_no_weight_in_groupTest=TRUE")))
  plain <- read_saige_table(plain_out)

  d <- noweights_run()$region                  # weighted file + No_Weights set
  nw <- d[d$Weight == "No_Weights" & d$Group != "Cauchy", ]
  pr <- plain[plain$Group != "Cauchy", ]
  expect_identical(nrow(nw), nrow(pr))
  expect_identical(nw$Group, pr$Group)
  expect_identical(nw$max_MAF, pr$max_MAF)
  for (col in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden",
                "SE_Burden", "MAC", "MAC_case", "MAC_control",
                "Number_rare", "Number_ultra_rare")) {
    expect_identical(nw[[col]], pr[[col]], label = paste("No_Weights ==", col))
  }
})
