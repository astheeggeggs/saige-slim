# The variance-ratio marker loop's cursor (LEDGER #75, FIXED).
#
# R/SAIGE_fitGLMM_fast.R advanced `indexInMarkerList` once unconditionally after
# reading a marker, and a SECOND time on each of the two skip branches (the
# chromosome parses to NA, or it is non-autosomal and
# includeNonautoMarkersforVarRatio is FALSE). The commented-out increment still
# sitting inside the accepted branch is the record of what happened: the
# increment was hoisted to the top of the iteration and the two copies were left
# behind. So a skipped marker consumed the marker AFTER it, which was never read
# at all; and because the exhaustion guard tested EQUALITY, a double step from
# the last position jumped straight over it and the next read was an NA that
# `if (genoInd == 0)` could not evaluate.
#
# Every bundled step-1 fixture is autosome-only, so none of this is reachable
# from a committed input. The fixtures here are DERIVED AT RUN TIME: a small
# plink file over the step-1 fixture's own 1,000 fam rows, every marker at
# MAC 40 so all of them enter the single (non-categorical) variance-ratio pool,
# with the chromosome column set per test. There is no golden -- the assertions
# are on exit status, on the number of markers the loop READS and TESTS, and on
# the files left behind.
#
# The counts are the point, and they are derived rather than observed. The
# cursor now advances by exactly one per iteration, so a pool of M markers is
# read in M reads and every usable marker in it is tested. At HEAD the cursor
# visited 1, 3, 5, ... on a fully-skipped pool and so read ceil(M/2) of them.
#
# NOTE the marker list is ordered by MAC, not by file order, so which markers a
# given cursor pattern lands on is NOT derivable -- only the totals are. Every
# assertion below is order-independent for that reason.

vr_mac <- 40L
vr_nmarkers <- 30L                  # --numRandomMarkerforVarianceRatio default
vr_step1_src <- file.path(saige_extdata, "input",
                          "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")
vr_cli <- run_saige_cli
vr_grm_base <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
vr_grm <- c(paste0("--sparseGRMFile=", vr_grm_base),
            paste0("--sparseGRMSampleIDFile=", vr_grm_base, ".sampleIDs.txt"))

# One marker per entry of `chr`, each at MAC vr_mac.
build_vr_fixture <- function(prefix, chr) {
  fam <- readLines(paste0(vr_step1_src, ".fam"))
  n <- length(fam)
  m <- length(chr)
  set.seed(7)
  writeLines(sprintf("%d\tvr%d\t0\t%d\tA\tT", chr, seq_len(m), seq_len(m)),
             paste0(prefix, ".bim"))
  file.copy(paste0(vr_step1_src, ".fam"), paste0(prefix, ".fam"), overwrite = TRUE)
  code <- c(3L, 2L, 0L)
  nb <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (i in seq_len(m)) {
    g <- integer(n)
    g[sample.int(n, vr_mac)] <- 1L
    cd <- c(code[g + 1L], rep(0L, nb * 4L - n))
    mm <- matrix(cd, nrow = 4L)
    writeBin(as.raw(mm[1, ] + 4L * mm[2, ] + 16L * mm[3, ] + 64L * mm[4, ]), con)
  }
  prefix
}

# Returns the exit status, the log, how many markers the loop READ, and how many
# of those were autosomal -- which is the same as how many were TESTED, since
# the chromosome is printed before the accept/skip branch and every autosomal
# marker here is accepted.
run_vr_step1 <- function(dir, chr, tag) {
  prefix <- build_vr_fixture(file.path(dir, tag), chr)
  out <- file.path(dir, paste0("out_", tag))
  args <- c(paste0("--plinkFile=", prefix), vr_grm, "--useSparseGRMtoFitNULL=TRUE",
            "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
            "--phenoCol=y_binary", "--covarColList=x1,x2",
            "--sampleIDColinphenoFile=IID", "--traitType=binary", "--LOCO=FALSE",
            "--IsOverwriteVarianceRatioFile=TRUE", paste0("--outputPrefix=", out))
  st <- vr_cli("step1_fitNULLGLMM.R", args)
  lg <- readLines(attr(st, "log"), warn = FALSE)
  list(status = as.integer(st), log = lg, prefix = out,
       read = sum(grepl("th marker in geno", lg, fixed = TRUE)),
       tested = sum(grepl("^CHR  1 $", lg)),
       skipped = sum(grepl("^CHR  23 $", lg)))
}

vr_na_error <- function(r) {
  any(grepl("missing value where TRUE/FALSE needed", r$log, fixed = TRUE))
}

test_that("NEGATIVE CONTROL: an all-autosomal pool estimates the ratio normally", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  r <- run_vr_step1(dir, rep(1L, 40L), "auto40")
  expect_identical(r$status, 0L)
  # nothing is skipped, so no skip branch runs and this arm is UNCHANGED by the
  # fix: the loop stops as soon as vr_nmarkers markers have been tested
  expect_identical(r$read, vr_nmarkers)
  expect_identical(r$tested, vr_nmarkers)
  expect_identical(r$skipped, 0L)
  vr <- utils::read.table(paste0(r$prefix, ".varianceRatio.txt"),
                          header = FALSE, stringsAsFactors = FALSE)
  expect_setequal(vr[[2]], c("null", "null_noXadj"))
  expect_true(all(is.finite(vr[[1]])))
})

test_that("LEDGER #75 FIXED: a fully-skipped pool is read once through, not halved", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  m <- 40L                                   # EVEN pool, every marker skipped
  r <- run_vr_step1(dir, rep(23L, m), "skip40")
  # THE SIGNATURE, inverted: every marker is read exactly once. At HEAD the
  # cursor visited 1, 3, 5, ... and read only 20 of the 40.
  expect_identical(r$read, m)
  expect_identical(r$skipped, m)
  expect_identical(r$tested, 0L)
  # The fix reads the whole pool; it does NOT invent data. With no usable marker
  # there is still no variance ratio to estimate, and the run still fails --
  # leaving the model on disk and the variance-ratio file empty, which is the
  # pair step 2 cannot use.
  expect_false(identical(r$status, 0L))
  expect_true(file.exists(paste0(r$prefix, ".rda")))
  expect_identical(file.info(paste0(r$prefix, ".varianceRatio.txt"))$size, 0)
})

test_that("LEDGER #75 FIXED: an odd pool no longer steps over the exhaustion guard", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  m <- 41L                                   # ODD pool, every marker skipped
  r <- run_vr_step1(dir, rep(23L, m), "skip41")
  expect_identical(r$read, m)
  expect_identical(r$skipped, m)
  # THE POINT: at HEAD the odd pool behaved DIFFERENTLY from the even one --
  # the cursor landed on M+2 rather than M+1, the equality guard did not fire,
  # and the next read was an NA that `if (genoInd == 0)` could not evaluate.
  # The asymmetry is gone: this now ends exactly as the even pool does.
  expect_false(vr_na_error(r))
  expect_false(any(grepl("genoInd == 0", r$log, fixed = TRUE)))
  expect_false(identical(r$status, 0L))
  expect_identical(file.info(paste0(r$prefix, ".varianceRatio.txt"))$size, 0)
})

test_that("LEDGER #75 FIXED: a mixed pool tests every marker it was asked for", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  # 30 usable markers -- exactly --numRandomMarkerforVarianceRatio -- among 30
  # skipped ones. At HEAD each skip consumed a list entry that was never read,
  # so the loop exhausted the pool having tested only 24 of the 30 requested
  # and returned a ratio estimated from a silently thinned marker set, at
  # exit 0 and with no warning. That is the consequence with no crash attached
  # to it, and the one a real run is most likely to meet.
  r <- run_vr_step1(dir, rep(c(23L, 1L), vr_nmarkers), "mix60")
  expect_identical(r$status, 0L)
  # DERIVED: exactly vr_nmarkers usable markers exist and the loop cannot stop
  # before it has tested that many, so it must find all of them. This is the
  # assertion with teeth -- HEAD tested 24.
  expect_identical(r$tested, vr_nmarkers)
  # NOT asserted exactly: how many entries were READ depends on where the
  # permutation places the last usable marker, which is a property of one
  # seeded sample() and not of the cursor. (Measured here at 60, i.e. the whole
  # list, but the CV never triggers a marker-count increase -- 1.2e-04 against
  # a 1e-03 cutoff -- so that 60 is incidental, not forced.) The derived
  # bounds are: at least vr_nmarkers reads to test vr_nmarkers, and at most the
  # length of the list.
  expect_gte(r$read, vr_nmarkers)
  expect_lte(r$read, 2L * vr_nmarkers)
  expect_identical(r$read, r$tested + r$skipped)
  vr <- utils::read.table(paste0(r$prefix, ".varianceRatio.txt"),
                          header = FALSE, stringsAsFactors = FALSE)
  expect_setequal(vr[[2]], c("null", "null_noXadj"))
  expect_true(all(is.finite(vr[[1]])))
})
