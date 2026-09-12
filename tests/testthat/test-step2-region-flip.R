# The FLIPPED ultra-rare pseudo-marker on the region path (2026-09-12).
#
# No golden reaches this. `--maxMAF_in_groupTest=0.5` admits every marker of
# the two bundled genes and `--MACCutoff_to_CollapseUltraRare=2000` is above
# 2N, so EVERY marker is ultra-rare and each (region x annotation) cell is one
# collapsed pseudo-marker `max(G)` over 14-35 markers. FIVE of the six cells
# then carry an alt frequency above 0.5 and `imputeGenoAndFlip` (UTIL.cpp)
# flips them, which makes the pseudo-marker DENSE: every non-carrier holds 2.
# Since the carrier-space rewrite of the region path (2026-09-12, `d73a9981`)
# that density is a per-cell FILL (`CarrierStore::urFill`) rather than a stored
# N-vector, and a flipped cell keeps the dense routines rather than the
# reductions over the carriers -- a branch nothing else in the suite executes
# (sim50k never flips, and neither does any bundled golden).
#
# This is a characterization test: it pins what the code does TODAY, not a
# published expectation. What it holds, on BOTH region paths:
#
#  - the run succeeds, and the single-variant file holds the six pseudo-markers
#    and NOTHING else. That is the guard against a vacuous pass: if a change
#    stopped collapsing here, the frequency assertions below would go on
#    passing over ordinary un-flipped markers;
#  - `GENE1 missense;lof` reports AC 1556 of 2N = 2000, i.e. AF 0.778. The
#    report is the UNFLIPPED frequency by construction -- imputeGenoAndFlip
#    ends with `if(flip){ t_altFreq = 1 - t_altFreq; ... }` -- so a flipped row
#    is exactly one whose PRINTED AF exceeds 0.5;
#  - `GENE1 lof` at AF 0.374 is the ONE cell of the six that is NOT flipped,
#    and it is the in-fixture negative control: without it "AF > 0.5" could be
#    read as a property of every UR row, and a build that flipped everything
#    (or nothing) would still be green. The five/one split is also the count
#    that must not be confused with the memo line's "4 evaluated" below --
#    that 4 is 6 cells less 2 twins, and has nothing to do with the flip;
#  - the region table's MAC is the FOLDED count of the same pseudo-marker,
#    min(AC, 2N - AC) -- one quantity on the two sides of the R/C++ bridge;
#  - each arm REPORTS which of the two per-cell memos it drove (LEDGER #30's
#    shape: a path that silently stopped being taken must not pass as green).
#    `served` is a measurement and only its bound is asserted (CLAUDE.md);
#  - LEDGER #130, found writing this file: on the BURDEN arm the region row's
#    `Pvalue_Burden` equals the run's OWN single-marker p exactly on the
#    un-flipped cell and differs on all five flipped ones. Each cell holds ONE
#    marker, so its burden IS that marker's test and the two numbers are the
#    same quantity computed twice -- the un-flipped cell is the positive
#    control that says so, and the five flipped ones are the bug.

FL_G <- c("--bedFile=./input/genotype_100markers.bed",
          "--bimFile=./input/genotype_100markers.bim",
          "--famFile=./input/genotype_100markers.fam",
          "--AlleleOrder=alt-first")
# the frozen example model, P6's model flags
FL_MODEL <- c(
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
# P6's group flags with the two that reach the flip: max MAF 0.5 admits every
# marker, and a collapse cutoff above 2N makes every one of them ultra-rare.
FL_GROUP <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
              "--groupFile=./input/group_new_snpid.txt",
              "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
              "--maxMAF_in_groupTest=0.5", "--MACCutoff_to_CollapseUltraRare=2000",
              "--is_fastTest=FALSE", "--is_single_in_groupTest=TRUE")
# The models' 1,000 samples, NOT the fam file's 10,000 (helper-step2-flexrv.R).
FL_2N <- 2000
# the six cells, in output order: the annotation groups are printed in the
# group file's own annotation order, so `lof,missense` prints "missense;lof"
FL_MARKERS <- paste0(rep(c("GENE1", "GENE2"), each = 3),
                     c(":lof:0.5", ":missense;lof:0.5", ":missense;lof;synonymous:0.5"))
FL_FLIPPED <- "GENE1:missense;lof:0.5"     # AC 1556 of 2000, flipped
FL_UNFLIPPED <- "GENE1:lof:0.5"            # AC 748, not flipped

# helper-saige.R's functions, bound at top level: lintr cannot see them from
# inside a function body (the aliasing test-step2-region-threads.R uses).
fl_skip <- skip_unless_level
fl_cli <- run_saige_cli
fl_ok <- expect_cli_ok
fl_read <- read_saige_table

# One CLI run per arm, cached per worker process.
fl_cache <- new.env(parent = emptyenv())
fl_run <- function(extra, tag) {
  fl_skip("standard")
  if (!is.null(fl_cache[[tag]])) return(fl_cache[[tag]])
  out <- file.path(tempfile(paste0("step2_flip_", tag, "_")), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  status <- fl_cli("step2_SPAtests.R",
                   c(FL_G, FL_MODEL, FL_GROUP, extra, paste0("--SAIGEOutputFile=", out)),
                   log = log)
  res <- list(status = status, path = out, log = readLines(log, warn = FALSE),
              region = NULL, single = NULL)
  if (identical(as.integer(status), 0L) && file.exists(out)) {
    res$region <- fl_read(out)
    res$single <- fl_read(paste0(out, ".singleAssoc.txt"))
  }
  fl_cache[[tag]] <- res
  res
}
fl_skato <- function() fl_run(character(0), "skato")
# the run's own single-marker p per cell, in the region table's row order
fl_single_p <- function(res) as.numeric(res$single[["p.value"]])
fl_region_p <- function(res, col) {
  r <- res$region[res$region$Group != "Cauchy", ]
  as.numeric(r[[col]])
}
fl_burden <- function() fl_run("--r.corr=1", "burden")

# `<kind> cells: A asked for, E evaluated, S served ...`, printed once per run
# whether or not the memo fired (R/BurdenMemo.R, R/SkatoParallel.R).
fl_memo <- function(res, kind) {
  line <- grep(paste0("^", kind, " cells: "), res$log, value = TRUE)
  testthat::expect_identical(length(line), 1L, label = paste(kind, "memo line"))
  n <- as.integer(regmatches(line, gregexpr("[0-9]+", line))[[1]])
  testthat::expect_identical(n[2] + n[3], n[1], label = paste(kind, "memo adds up"))
  n
}

# Everything both arms must show: the collapse, the flip, and the folded MAC.
fl_expect_flip <- function(res, label) {
  fl_ok(res$status)
  s <- res$single
  # nothing but the six pseudo-markers survives the collapse at cutoff 2000
  testthat::expect_identical(s$MarkerID, FL_MARKERS, label = paste(label, "pseudo-markers"))
  testthat::expect_identical(unique(s$CHR), "UR", label = paste(label, "CHR"))
  ac <- as.numeric(s$AC_Allele2)
  af <- as.numeric(s$AF_Allele2)
  # the printed AF is the count's own, un-folded: that is what makes "> 0.5"
  # readable as "this cell was flipped"
  testthat::expect_equal(af, ac / FL_2N, tolerance = 1e-6, label = paste(label, "AF"))
  flipped <- af > 0.5
  testthat::expect_identical(sum(flipped), 5L, label = paste(label, "flipped count"))
  # the drill's cell, and the in-fixture negative control
  testthat::expect_equal(ac[s$MarkerID == FL_FLIPPED], 1556, label = paste(label, "AC"))
  testthat::expect_equal(af[s$MarkerID == FL_FLIPPED], 0.778, tolerance = 1e-6,
                         label = paste(label, "flipped AF"))
  testthat::expect_true(flipped[s$MarkerID == FL_FLIPPED])
  testthat::expect_false(flipped[s$MarkerID == FL_UNFLIPPED])
  testthat::expect_equal(af[s$MarkerID == FL_UNFLIPPED], 0.374, tolerance = 1e-6,
                         label = paste(label, "unflipped AF"))
  # the region table: the same six cells, every marker ultra-rare, and the
  # group MAC the folded allele count of the pseudo-marker
  r <- res$region[res$region$Group != "Cauchy", ]
  testthat::expect_identical(paste0(r$Region, ":", r$Group, ":", r$max_MAF), FL_MARKERS,
                             label = paste(label, "region rows"))
  testthat::expect_identical(unique(r$Number_rare), "0", label = paste(label, "Number_rare"))
  testthat::expect_true(all(as.numeric(r$Number_ultra_rare) > 0))
  testthat::expect_equal(as.numeric(r$MAC), pmin(ac, FL_2N - ac), label = paste(label, "MAC"))
}

test_that("SKAT-O: every marker collapses and five of six pseudo-markers flip", {
  res <- fl_skato()
  fl_expect_flip(res, "skato")
  # the cell is tested singly on this path, so every cell has a p-value
  r <- res$region[res$region$Group != "Cauchy", ]
  expect_true(all(as.numeric(r$Pvalue) > 0))
  expect_true(all(as.numeric(r$Pvalue_Burden) > 0))
  # the drill's flipped cell, pinned by value. Pvalue == Pvalue_SKAT ==
  # Pvalue_Burden here: a cell of one marker, tested through the per-marker
  # exports rather than the BURDEN emitter (LEDGER #130).
  i <- which(r$Region == "GENE1" & r$Group == "missense;lof")
  expect_equal(as.numeric(r$Pvalue[i]), 0.8334812, tolerance = 1e-5)
  expect_equal(as.numeric(r$Pvalue_SKAT[i]), 0.8334812, tolerance = 1e-5)
  expect_equal(as.numeric(r$Pvalue_Burden[i]), 0.8334812, tolerance = 1e-5)
  # the same single-marker file both arms write, byte for byte: the ultra-rare
  # prepare and its single tests are ONE path re-run, not two paths compared.
  expect_identical(readLines(paste0(res$path, ".singleAssoc.txt")),
                   readLines(paste0(fl_burden()$path, ".singleAssoc.txt")))
  # this arm drove the SKAT-O memo and not the BURDEN one
  expect_identical(fl_memo(res, "SKAT-O")[1], 6L)
  expect_identical(fl_memo(res, "BURDEN")[1], 0L)
  # GENE1/GENE2 `missense;lof` and `missense;lof;synonymous` are the same
  # collapse (no marker of the fixture is annotated synonymous), so the memo
  # has twins to serve; how many it serves is a measurement (CLAUDE.md).
  expect_gte(fl_memo(res, "SKAT-O")[3], 1L)
  expect_lte(fl_memo(res, "SKAT-O")[3], 2L)
})

test_that("LEDGER #130: the BURDEN emitter's flipped cells disagree with their own single test", {
  res <- fl_burden()
  fl_ok(res$status)
  single <- fl_single_p(res)
  region <- fl_region_p(res, "Pvalue_Burden")
  flipped <- as.numeric(res$single$AF_Allele2) > 0.5
  # the positive control FIRST: on the one cell that does not flip, the burden
  # of a one-marker cell is that marker's own test, to every printed digit.
  expect_equal(region[!flipped], single[!flipped], tolerance = 1e-6)
  # and on every flipped cell it is not -- the bug, pinned by value. Both
  # numbers are plain normal score tests over the SAME single marker (no SPA,
  # no ER, one variance-ratio band), so this is a contradiction, not a
  # tolerance: see the row's PROVENANCE for the two causes measured out.
  expect_true(all(abs(region[flipped] - single[flipped]) > 1e-3))
  expect_equal(single, c(0.2516356, 0.8334812, 0.8334812, 0.1663821, 0.7175196, 0.7175196),
               tolerance = 1e-5)
  expect_equal(region, c(0.2516356, 0.5758082, 0.5758082, 0.2634589, 0.7021056, 0.7021056),
               tolerance = 1e-5)
  # the sign correction `Beta * (1 - 2*flip)` is applied on the single path and
  # not by the emitter, so the same one-marker cell reports opposite signs
  r <- res$region[res$region$Group != "Cauchy", ]
  cell <- r$Region == "GENE1" & r$Group == "missense;lof"
  expect_lt(as.numeric(res$single$BETA[res$single$MarkerID == FL_FLIPPED]), 0)
  expect_gt(as.numeric(r$BETA_Burden[cell]), 0)
  expect_equal(as.numeric(r$BETA_Burden[cell]), 2.34563, tolerance = 1e-5)
  expect_equal(as.numeric(r$SE_Burden[cell]), 4.19223, tolerance = 1e-5)
})

test_that("BURDEN: the flipped pseudo-marker reaches the burden emitter too", {
  res <- fl_burden()
  fl_expect_flip(res, "burden")
  r <- res$region[res$region$Group != "Cauchy", ]
  expect_true(all(as.numeric(r$Pvalue_Burden) > 0))
  expect_false("Pvalue_SKAT" %in% colnames(r))
  # and this arm drove the BURDEN memo and not the SKAT-O one
  expect_identical(fl_memo(res, "BURDEN")[1], 6L)
  expect_identical(fl_memo(res, "SKAT-O")[1], 0L)
  expect_gte(fl_memo(res, "BURDEN")[3], 1L)
  expect_lte(fl_memo(res, "BURDEN")[3], 2L)
})
