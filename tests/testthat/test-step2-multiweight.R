# Multi-weight (K > 1) region tests: SCHEMA AND LAYOUT.
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
# Conditional analysis lives in test-step2-multiweight-cond.R; weighting
# behaviour and the ledger cases in test-step2-multiweight-behaviour.R.
# Shared scaffolding is in helper-multiweight.R.
#
# Tier: standard (runs the CLI), like the golden tests.

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
  drop_rn <- function(d) {
    rownames(d) <- NULL
    d
  }
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
