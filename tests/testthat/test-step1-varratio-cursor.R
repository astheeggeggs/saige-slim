# The variance-ratio marker loop's cursor (LEDGER #75), and the arithmetic
# signature of its double increment.
#
# R/SAIGE_fitGLMM_fast.R advances `indexInMarkerList` once unconditionally after
# reading a marker, and a SECOND time on each of the two skip branches (the
# chromosome parses to NA, or it is non-autosomal and
# includeNonautoMarkersforVarRatio is FALSE). So a skipped marker consumes the
# marker AFTER it as well, which is never read at all; and because the
# exhaustion guard tests EQUALITY (`indexInMarkerList - 1 == length(list)`), a
# double step from the last position jumps straight over it and the next read is
# an NA that `if (genoInd == 0)` cannot evaluate.
#
# Every bundled step-1 fixture is autosome-only, so none of this is reachable
# from a committed input. The fixtures here are DERIVED AT RUN TIME: a small
# plink file over the step-1 fixture's own 1,000 fam rows, every marker at
# MAC 40 so all of them enter the single (non-categorical) variance-ratio pool,
# with the chromosome column set per test. There is no golden -- the assertions
# are on exit status, on the number of markers the loop READS, and on the files
# left behind.
#
# The read counts are the point. With all M markers skipped the cursor visits
# 1, 3, 5, ... so it reads ceil(M/2) of them; a single-increment loop would read
# all M. That is what separates this from "the pool merely ran out".

vr_mac <- 40L
vr_step1_src <- file.path(saige_extdata, "input",
                          "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")
vr_cli <- run_saige_cli
vr_grm_base <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
vr_grm <- c(paste0("--sparseGRMFile=", vr_grm_base),
            paste0("--sparseGRMSampleIDFile=", vr_grm_base, ".sampleIDs.txt"))

# M markers at MAC vr_mac; the last n_nonauto of them on chromosome 23.
build_vr_fixture <- function(prefix, m, n_nonauto) {
  fam <- readLines(paste0(vr_step1_src, ".fam"))
  n <- length(fam)
  set.seed(7)
  chr <- c(rep(1L, m - n_nonauto), rep(23L, n_nonauto))
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

# Returns the exit status, the log lines, and how many markers the loop READ.
run_vr_step1 <- function(dir, m, n_nonauto) {
  prefix <- build_vr_fixture(file.path(dir, sprintf("vr%d_%d", m, n_nonauto)),
                             m, n_nonauto)
  out <- file.path(dir, sprintf("out%d_%d", m, n_nonauto))
  args <- c(paste0("--plinkFile=", prefix), vr_grm, "--useSparseGRMtoFitNULL=TRUE",
            "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
            "--phenoCol=y_binary", "--covarColList=x1,x2",
            "--sampleIDColinphenoFile=IID", "--traitType=binary", "--LOCO=FALSE",
            "--IsOverwriteVarianceRatioFile=TRUE", paste0("--outputPrefix=", out))
  st <- vr_cli("step1_fitNULLGLMM.R", args)
  lg <- readLines(attr(st, "log"), warn = FALSE)
  list(status = as.integer(st), log = lg, prefix = out,
       read = sum(grepl("th marker in geno", lg, fixed = TRUE)),
       accepted = sum(grepl("^CHR  1 $", lg)),
       skipped = sum(grepl("^CHR  23 $", lg)))
}

test_that("NEGATIVE CONTROL: an all-autosomal pool estimates the ratio normally", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  r <- run_vr_step1(dir, 40L, 0L)
  expect_identical(r$status, 0L)
  # nothing is skipped, so the cursor advances by one and stops as soon as
  # --numRandomMarkerforVarianceRatio (default 30) markers have been tested
  expect_identical(r$read, 30L)
  expect_identical(r$skipped, 0L)
  vr <- utils::read.table(paste0(r$prefix, ".varianceRatio.txt"),
                          header = FALSE, stringsAsFactors = FALSE)
  expect_setequal(vr[[2]], c("null", "null_noXadj"))
  expect_true(all(is.finite(vr[[1]])))
})

test_that("LEDGER #75: a skipped marker consumes the next one, so half the pool is never read", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  m <- 40L                                   # EVEN pool, every marker skipped
  r <- run_vr_step1(dir, m, m)
  # THE SIGNATURE: the cursor visits 1, 3, 5, ... so exactly half the pool is
  # examined. A single-increment loop would have read all 40.
  expect_identical(r$read, m %/% 2L)
  expect_identical(r$accepted, 0L)
  expect_identical(r$skipped, m %/% 2L)
  # An even pool lands the cursor exactly on the guard, so the loop exits
  # cleanly -- and then fails downstream instead, with no marker accepted.
  expect_false(identical(r$status, 0L))
  # the model is already on disk; the variance-ratio file is created and EMPTY,
  # which is the pair step 2 cannot use
  expect_true(file.exists(paste0(r$prefix, ".rda")))
  expect_identical(file.info(paste0(r$prefix, ".varianceRatio.txt"))$size, 0)
})

test_that("LEDGER #75: an odd pool steps over the exhaustion guard and reads an NA", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  m <- 41L                                   # ODD pool, every marker skipped
  r <- run_vr_step1(dir, m, m)
  # One marker more than the even case, and the failure MODE changes: the
  # cursor lands on M+2 rather than M+1, the equality guard does not fire, and
  # the next read is off the end of the list.
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("missing value where TRUE/FALSE needed", r$log, fixed = TRUE)))
  expect_true(any(grepl("genoInd == 0", r$log, fixed = TRUE)))
  # ceil(M/2) reads plus the one that read past the end
  expect_identical(r$read, as.integer(ceiling(m / 2) + 1L))
  expect_true(file.exists(paste0(r$prefix, ".rda")))
  expect_identical(file.info(paste0(r$prefix, ".varianceRatio.txt"))$size, 0)
})
