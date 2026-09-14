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
# Half (b), making C++ honour R's layout on an all-zero matrix, is a
# behaviour change on a reachable input and stays open for the owner.
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

test_that("an all-zero K=2 region is a loud error, not a column-shifted file (LEDGER #146)", {
  skip_unless_level("standard")
  dir <- tempfile("step2_l146_")
  dir.create(dir, recursive = TRUE)
  grp <- l146_zero_fixture(dir)
  out <- file.path(dir, "burden.txt")
  args <- c(
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
    paste0("--groupFile=", grp),
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--r.corr=1", "--is_fastTest=TRUE",
    paste0("--SAIGEOutputFile=", out)
  )
  st <- l146_cli("step2_SPAtests.R", args)
  lg <- readLines(attr(st, "log"), warn = FALSE)
  # the inverted characterization: before the fix this exited 0 and the
  # header (built for 2 sets, with a Weight column) had one more field than
  # GENE1's rows (written for 1 implicit set, without it)
  expect_false(identical(as.integer(st), 0L))
  expect_true(any(grepl("LEDGER #146", lg, fixed = TRUE)))
  expect_true(any(grepl("2 weight set(s)", lg, fixed = TRUE)))
  if (file.exists(out)) {
    rows <- readLines(out, warn = FALSE)
    nf <- lengths(strsplit(rows, "\t", fixed = TRUE))
    expect_true(all(nf == nf[1]))
  }
})
