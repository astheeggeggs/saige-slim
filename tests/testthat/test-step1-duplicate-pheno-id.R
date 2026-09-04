# A duplicated sample ID in the phenotype file (LEDGER #67).
#
# R/SAIGE_fitGLMM_fast.R merges the phenotype frame onto the genotype sample list
# with `merge(mmat_nomissing, sampleListwithGeno, by.x = "IID", by.y = "IIDgeno")`
# and NOTHING anywhere in the step-1 phenotype path checks for duplicate IDs --
# the only `duplicated()` call in that file is on SampleIDIncludeFile. A repeated
# ID therefore multiplies that sample's rows through the merge, and the inflated
# frame becomes `subPheno`, whose `IndexGeno` column is the vector `setgeno` uses
# to choose genotype columns.
#
# The consequence is an internal inconsistency rather than a crash, and it is
# arithmetic rather than a guess. In setGenoObjFromSource,
# `Nnomissing = subSampleInGeno.size()` counts the DUPLICATED list, while the
# allele-count accumulation is gated on `indicatorGenoSamplesWithPheno_in[ind]`
# with `ind` running over the FILE's samples -- so each distinct sample is
# counted at most once. `altFreq = alleleCount / float(Nnomissing * 2)` then puts
# the duplicated sample in the DENOMINATOR but not in the NUMERATOR.

dup_src <- file.path(saige_extdata, "input",
                     "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")
dup_pheno <- file.path(saige_extdata, "input",
                       "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
dup_cli <- run_saige_cli
dup_grm_base <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"

# M markers at MAC 40; the FIRST fam row carries the ALT allele in markers
# 1..M/2 and not in the rest, so one run exercises both arms.
build_carrier_fixture <- function(prefix, m = 40L) {
  fam <- readLines(paste0(dup_src, ".fam"))
  n <- length(fam)
  set.seed(11)
  geno <- t(vapply(seq_len(m), function(i) {
    g <- integer(n)
    g[sample(2:n, 39)] <- 1L
    g[1] <- if (i <= m %/% 2L) 1L else 0L
    g
  }, integer(n)))
  writeLines(sprintf("1\tc%d\t0\t%d\tA\tT", seq_len(m), seq_len(m)), paste0(prefix, ".bim"))
  file.copy(paste0(dup_src, ".fam"), paste0(prefix, ".fam"), overwrite = TRUE)
  code <- c(3L, 2L, 0L)
  nb <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (i in seq_len(m)) {
    cd <- c(code[geno[i, ] + 1L], rep(0L, nb * 4L - n))
    mm <- matrix(cd, nrow = 4L)
    writeBin(as.raw(mm[1, ] + 4L * mm[2, ] + 16L * mm[3, ] + 64L * mm[4, ]), con)
  }
  list(prefix = prefix, n = n, m = m, g_dup = geno[, 1])
}

test_that("LEDGER #67: a duplicated phenotype ID inflates the analysed sample count", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- tempfile()
  dir.create(d)
  fx <- build_carrier_fixture(file.path(d, "carrier"))

  ph <- readLines(dup_pheno)
  writeLines(c(ph, ph[2]), file.path(d, "dup.txt"))     # one row repeated verbatim
  writeLines(ph, file.path(d, "orig.txt"))
  iid <- vapply(strsplit(utils::tail(readLines(file.path(d, "dup.txt")), -1), "[ \t]+"),
                `[`, character(1), 5L)
  expect_identical(length(iid), 1001L)                  # 1001 rows ...
  expect_identical(length(unique(iid)), 1000L)          # ... 1000 distinct people

  runs <- lapply(c("orig", "dup"), function(tag) {
    args <- c(paste0("--plinkFile=", fx$prefix),
              paste0("--sparseGRMFile=", dup_grm_base),
              paste0("--sparseGRMSampleIDFile=", dup_grm_base, ".sampleIDs.txt"),
              "--useSparseGRMtoFitNULL=TRUE",
              paste0("--phenoFile=", file.path(d, paste0(tag, ".txt"))),
              "--phenoCol=y_binary", "--covarColList=x1,x2",
              "--sampleIDColinphenoFile=IID", "--traitType=binary", "--LOCO=FALSE",
              "--IsOverwriteVarianceRatioFile=TRUE",
              paste0("--outputPrefix=", file.path(d, paste0("out_", tag))))
    st <- dup_cli("step1_fitNULLGLMM.R", args)
    lg <- readLines(attr(st, "log"), warn = FALSE)
    n <- sub("^\\s*([0-9]+).*", "\\1",
             grep("samples will be used for analysis", lg, value = TRUE)[1])
    list(status = as.integer(st), n = as.integer(n))
  })
  # NEGATIVE CONTROL: the unmodified file reports the true count
  expect_identical(runs[[1]]$status, 0L)
  expect_identical(runs[[1]]$n, 1000L)
  # and the duplicated one carries the repeated individual as a second sample,
  # at EXIT 0 and with no warning
  expect_identical(runs[[2]]$status, 0L)
  expect_identical(runs[[2]]$n, 1001L)
})

test_that("LEDGER #67: the duplicate lands in the frequency denominator, not the numerator", {
  skip_unless_level("standard")
  d <- tempfile()
  dir.create(d)
  fx <- build_carrier_fixture(file.path(d, "carrier"))

  # setgeno writes C++ globals, so it is driven in a child process. The two
  # arguments are exactly what the R side passes after a merge with one
  # duplicated ID: the index vector gains a repeat, the indicator keeps one
  # entry per FILE sample.
  out <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('pre <- "%s"; out <- "%s"; n <- %dL; m <- %dL', fx$prefix, out, fx$n, fx$m),
    "setminMAFforGRM(0); setmaxMissingRateforGRM(1)",
    "sub <- c(seq_len(n), 1L)",
    "ind <- rep(TRUE, n)",
    'setgeno(paste0(pre, ".bed"), paste0(pre, ".bim"), paste0(pre, ".fam"),',
    "        sub, ind, 2, FALSE)",
    "N <- getNnomissingOut()",
    "saveRDS(list(N = N, freq = getAlleleFreqVec(),",
    "             store = vapply(seq_len(m) - 1L, function(i) Get_OneSNP_Geno(i),",
    "                            numeric(N))), out)"
  ), script)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = FALSE, stderr = FALSE))
  expect_identical(as.integer(st), 0L)
  s <- readRDS(out)

  # the store is sized by the DUPLICATED list, not by the file
  expect_identical(s$N, fx$n + 1L)

  # getAlleleFreqVec() comes back as a one-column matrix (arma::fvec)
  reported_ac <- as.vector(s$freq) * 2 * s$N  # the numerator the frequency implies
  store_ac <- colSums(s$store)           # what the store actually holds
  # THE POINT: the gap is exactly the duplicated sample's genotype -- 1 where it
  # carries the ALT allele, 0 where it does not. Compared at float32 tolerance
  # because altFreq is a float.
  expect_equal(store_ac - reported_ac, as.numeric(fx$g_dup), tolerance = 1e-5)
  # both arms are present, so this is not vacuous
  expect_gt(sum(fx$g_dup == 1L), 0L)
  expect_gt(sum(fx$g_dup == 0L), 0L)
  # and stated as the two arms: carriers disagree by one allele, non-carriers agree
  expect_true(all(abs((store_ac - reported_ac)[fx$g_dup == 1L] - 1) < 1e-4))
  expect_true(all(abs((store_ac - reported_ac)[fx$g_dup == 0L]) < 1e-4))
})
