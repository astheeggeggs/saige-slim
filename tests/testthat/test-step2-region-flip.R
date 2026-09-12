# The high-alt-frequency ultra-rare pseudo-marker on the region path
# (2026-09-12; LEDGER #130 FIXED in the same commit as this inversion).
#
# `--maxMAF_in_groupTest=0.5` admits every marker of the two bundled genes and
# `--MACCutoff_to_CollapseUltraRare=2000` is above 2N, so EVERY marker is
# ultra-rare and each (region x annotation) cell is one collapsed pseudo-marker
# `max(G)` over 14-35 markers. FIVE of the six cells then carry an alt
# frequency above 0.5. Until #130 was fixed `imputeGenoAndFlip` (UTIL.cpp)
# folded those five to the minor allele -- replacing the vector with `2 - G`,
# which made it DENSE (every non-carrier holding 2) -- while the BURDEN
# emitter read the POST-fold values at the PRE-fold carrier positions, so one
# collapsed marker was tested as two different vectors and the emitted column
# was neither. The pseudo-marker is now never folded: the group file's
# annotation is a property of the ALT allele, and this marker IS the collapse
# of ALT-annotated ultra-rare variants, so its burden counts the ALT allele
# (src/Main.cpp, the note at the collapse). Goldens P62/P63.
#
# What it holds, on BOTH region paths:
#
#  - the run succeeds, and the single-variant file holds the six pseudo-markers
#    and NOTHING else. That is the guard against a vacuous pass: if a change
#    stopped collapsing here, the frequency assertions below would go on
#    passing over ordinary un-flipped markers;
#  - `GENE1 missense;lof` reports AC 1556 of 2N = 2000, i.e. AF 0.778. The
#    report is the ALT frequency, which is what makes "AF > 0.5" readable as
#    "this cell is one #130 used to fold". The printed AC/AF/MAC did not move
#    with the fix: they never came from the folded vector (imputeGenoAndFlip
#    ended with `if(flip){ t_altFreq = 1 - t_altFreq; ... }`, undoing its own
#    fold for the report), and MAF/MAC are still min(AF, 1 - AF);
#  - `GENE1 lof` at AF 0.374 is the ONE cell of the six BELOW 0.5, and it is
#    the in-fixture negative control: it is the cell the old code and the new
#    agree on byte for byte, so it is what says the reconciliation below is a
#    property of the statistic and not of the arithmetic. The five/one split is
#    also the count that must not be confused with the memo line's "4
#    evaluated" below -- that 4 is 6 cells less 2 twins;
#  - the region table's MAC is the FOLDED count of the same pseudo-marker,
#    min(AC, 2N - AC) -- one quantity on the two sides of the R/C++ bridge;
#  - each arm REPORTS which of the two per-cell memos it drove (LEDGER #30's
#    shape: a path that silently stopped being taken must not pass as green).
#    `served` is a measurement and only its bound is asserted (CLAUDE.md);
#  - LEDGER #130, found writing this file and FIXED with this inversion: each
#    cell holds ONE marker, so its burden IS that marker's test, and the region
#    row's `Pvalue_Burden` now equals the run's OWN single-marker p on all six
#    cells and on BOTH arms -- one quantity computed three ways (the BURDEN
#    emitter's summed column in C++, the SKAT-O path's Score/Phi in R, and the
#    per-marker test). Before the fix only `GENE1 lof` agreed. The sign goes
#    with it: `BETA_Burden` is negative wherever the single row's `BETA` is,
#    because both now count the ALT allele.

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
FL_FLIPPED <- "GENE1:missense;lof:0.5"     # AC 1556 of 2000, AF 0.778 (was folded)
FL_UNFLIPPED <- "GENE1:lof:0.5"            # AC 748, AF 0.374 (never folded)

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
  # the printed AF is the count's own ALT frequency: that is what makes "> 0.5"
  # readable as "this cell is one #130 used to fold"
  testthat::expect_equal(af, ac / FL_2N, tolerance = 1e-6, label = paste(label, "AF"))
  flipped <- af > 0.5
  testthat::expect_identical(sum(flipped), 5L, label = paste(label, "AF > 0.5 count"))
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

test_that("SKAT-O: every marker collapses and five of six cells sit above AF 0.5", {
  res <- fl_skato()
  fl_expect_flip(res, "skato")
  # the cell is tested singly on this path, so every cell has a p-value
  r <- res$region[res$region$Group != "Cauchy", ]
  expect_true(all(as.numeric(r$Pvalue) > 0))
  expect_true(all(as.numeric(r$Pvalue_Burden) > 0))
  # the drill's cell, pinned by value. Pvalue == Pvalue_SKAT == Pvalue_Burden
  # here: a cell of one marker, tested through the per-marker exports.
  i <- which(r$Region == "GENE1" & r$Group == "missense;lof")
  expect_equal(as.numeric(r$Pvalue[i]), 0.8334812, tolerance = 1e-5)
  expect_equal(as.numeric(r$Pvalue_SKAT[i]), 0.8334812, tolerance = 1e-5)
  expect_equal(as.numeric(r$Pvalue_Burden[i]), 0.8334812, tolerance = 1e-5)
  # LEDGER #130 moved this arm too, by a route the row did not name: on a
  # BINARY trait region_postprocess rescales Phi through `g.sum`, which IS the
  # burden column the emitter builds (get_newPhi_scaleFactor_traitType,
  # R/SAIGE_SPATest_Region.R). With the column fixed, this arm's cells also
  # reproduce their own single-marker p exactly -- `GENE2 lof` went
  # 0.166847187 -> 0.1663821, which is the single row's p (golden P63).
  expect_equal(as.numeric(r$Pvalue_Burden), fl_single_p(res), tolerance = 1e-6)
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

test_that("LEDGER #130 FIXED: every one-marker cell's burden IS its own single test", {
  res <- fl_burden()
  fl_ok(res$status)
  single <- fl_single_p(res)
  region <- fl_region_p(res, "Pvalue_Burden")
  flipped <- as.numeric(res$single$AF_Allele2) > 0.5
  # The positive control FIRST, and it is the one cell the fix cannot have
  # moved: below AF 0.5 nothing was ever folded, and the burden of a
  # one-marker cell is that marker's own test to every printed digit.
  expect_equal(region[!flipped], single[!flipped], tolerance = 1e-6)
  # INVERTED (was: differ by > 1e-3 on all five). The five cells above AF 0.5
  # now agree as well. Both numbers are plain normal score tests over the SAME
  # single marker (no SPA, no ER, one variance-ratio band), so agreement here
  # is the reconciliation the row asked for, not a tolerance being generous.
  expect_equal(region[flipped], single[flipped], tolerance = 1e-6)
  expect_equal(single, c(0.2516356, 0.8334812, 0.8334812, 0.1663821, 0.7175196, 0.7175196),
               tolerance = 1e-5)
  # was c(0.2516356, 0.5758082, 0.5758082, 0.2634589, 0.7021056, 0.7021056)
  expect_equal(region, single, tolerance = 1e-6)
  # The SIGN. Nothing is folded any more, so no `Beta * (1 - 2*flip)`
  # correction is in play on either path and both count the ALT allele: the
  # emitter's BETA_Burden is negative wherever the single row's BETA is. Its
  # MAGNITUDE is the single Beta divided by the cell's weight -- the burden
  # column is w0 * G with w0 = dbeta(MAF, 1, 25) = 0.060565 at MAF 0.222, and
  # 0.05205 / 0.060565 = 0.859 against the printed 0.860935. The last 0.2 % is
  # the binary Phi rescale (get_newPhi_scaleFactor_traitType), so this is a
  # two-figure reconciliation of the WEIGHT, not an identity -- the identity in
  # this file is the p-value one above. (Was +2.34563, which was neither
  # reading.)
  r <- res$region[res$region$Group != "Cauchy", ]
  cell <- r$Region == "GENE1" & r$Group == "missense;lof"
  expect_lt(as.numeric(res$single$BETA[res$single$MarkerID == FL_FLIPPED]), 0)
  expect_lt(as.numeric(r$BETA_Burden[cell]), 0)
  expect_equal(as.numeric(r$BETA_Burden[cell]), -0.860935, tolerance = 1e-5)
  expect_equal(as.numeric(r$SE_Burden[cell]), 4.09503, tolerance = 1e-5)
  # and the two region paths now report the SAME burden coefficient for the
  # same cell -- sum(Score)/sum(Phi) in R against getMarkerPval on the summed
  # column in C++, one quantity computed two ways (printed precision, not bits)
  rs <- fl_skato()$region
  cs <- rs$Region == "GENE1" & rs$Group == "missense;lof"
  expect_equal(as.numeric(r$BETA_Burden[cell]), as.numeric(rs$BETA_Burden[cs]),
               tolerance = 1e-5)
})

test_that("BURDEN: the high-AF pseudo-marker reaches the burden emitter too", {
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
