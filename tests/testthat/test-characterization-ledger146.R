# LEDGER #146 (half (a) FIXED 2026-09-14): the BURDEN writer's Weight column
# can no longer disagree with the header it sits under.
#
# R composes the weight sets (K custom sets from the group file, beta sets,
# No_Weights, FlexRV transforms) and builds the header for that count; C++
# composes them again from the same inputs. The one place the two could
# diverge: a region whose ENTIRE custom weight matrix is zero makes C++ drop
# the custom sets (`hasCustomWeights` is false) and fall through to ONE
# implicit beta set, while R's header was built for K. With K >= 2 that was a
# header WITH a Weight column over rows WITHOUT it, at exit 0 -- a
# column-shifted results file -- and the `"W" + r` label fallbacks existed
# only because the code tolerated that size disagreement. C++ now asserts
# `q_weight == t_weightNames.size()` once the sets are composed and the two
# fallbacks are gone, so the same input is a loud error naming the counts.
# Half (b), FIXED 2026-09-15 (owner's decision): an all-zero matrix composes
# its K sets like any other, so every GENE1 cell is empty under the region
# path's existing rule (no row for a cell with nothing to test -- the same
# thing a gene whose annotation is absent gets), GENE1's pooled row is NA
# with no inputs, and GENE2 is untouched. The stop of half (a) stays as an
# invariant that can no longer fire on this input.
#
# Fixture derived at run time: the K = 2 multi-weight group file with GENE1's
# two weight lines zeroed (GENE2 untouched, so the run reaches GENE1 with a
# valid composition on record). Tier: standard (one BURDEN run).

l146_cli <- function(...) get("run_saige_cli")(...)
l146_extdata <- function() get("saige_extdata")

l146_zero_fixture <- function(dir) {
  src <- file.path(l146_extdata(), "input", "group_new_snpid_multiWeights.txt")
  lines <- readLines(src)
  is_w <- grepl("^GENE1 weight", lines)
  stopifnot(sum(is_w) == 2L)
  lines[is_w] <- vapply(lines[is_w], function(l) {
    f <- strsplit(l, " ", fixed = TRUE)[[1]]
    paste(c(f[1:2], rep("0", length(f) - 2L)), collapse = " ")
  }, character(1))
  path <- file.path(dir, "group_zero_gene1.txt")
  writeLines(lines, path)
  path
}

test_that("an all-zero K=2 region is empty, not a stop and not a column shift (LEDGER #146)", {
  skip_unless_level("standard")
  dir <- tempfile("step2_l146_")
  dir.create(dir, recursive = TRUE)
  grp <- l146_zero_fixture(dir)
  out <- file.path(dir, "burden.txt")
  common <- c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--r.corr=1", "--is_fastTest=TRUE"
  )
  st <- l146_cli("step2_SPAtests.R", c(common, paste0("--groupFile=", grp),
                                       paste0("--SAIGEOutputFile=", out)))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  # inverted 2026-09-15: before (b) this stopped naming "2 weight set(s)"
  expect_identical(as.integer(st), 0L)
  expect_false(any(grepl("LEDGER #146", lg, fixed = TRUE)))
  tab <- read.delim(out, stringsAsFactors = FALSE, check.names = FALSE)
  # every row has the header's field count (the (a) defect was a column shift)
  nf <- lengths(strsplit(readLines(out, warn = FALSE), "\t", fixed = TRUE))
  expect_true(all(nf == nf[1]))
  # GENE1: no set rows (every cell empty), one pooled row with nothing in it
  g1 <- tab[tab$Region == "GENE1", ]
  expect_identical(nrow(g1), 1L)
  expect_identical(g1$Group, "Cauchy")
  expect_true(is.na(g1$Pvalue_Burden))
  # GENE2: byte-identical rows to the untouched K = 2 fixture's run
  ref <- file.path(dir, "burden_ref.txt")
  st2 <- l146_cli("step2_SPAtests.R", c(common,
    "--groupFile=./input/group_new_snpid_multiWeights.txt", paste0("--SAIGEOutputFile=", ref)))
  expect_identical(as.integer(st2), 0L)
  g2 <- grep("^GENE2\t", readLines(out, warn = FALSE), value = TRUE)
  g2_ref <- grep("^GENE2\t", readLines(ref, warn = FALSE), value = TRUE)
  expect_true(length(g2) > 1L)
  expect_identical(g2, g2_ref)
  # the disclosure: the pooled sidecar says GENE1 combined nothing
  pooled <- paste0(out, ".pooledTests.txt")
  expect_true(file.exists(pooled))
  pt <- read.delim(pooled, stringsAsFactors = FALSE, check.names = FALSE)
  expect_true("GENE1" %in% pt$Region)
})
