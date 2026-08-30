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

# One CLI run per group file, cached across the test_that blocks below.
mw_cache <- new.env(parent = emptyenv())
mw_run <- function(groupFile) {
  skip_unless_level("standard")
  key <- groupFile
  if (!is.null(mw_cache[[key]])) return(mw_cache[[key]])
  out <- file.path(tempfile("step2_mw_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R", MW_ARGS(groupFile, out))
  expect_cli_ok(status)
  res <- list(region = read_saige_table(out),
              single = read_saige_table(paste0(out, ".singleAssoc.txt")))
  mw_cache[[key]] <- res
  res
}

multi_run <- function() mw_run("group_new_snpid_multiWeights.txt")
single_run <- function() mw_run("group_new_snpid_withWeights.txt")

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
