# preparePhenotypeFrame (R/PhenoPrep.R): the phenotype-preparation block of
# fitNULLGLMM as one function, exercised on a synthetic phenotype file with a
# synthetic genotype sample list, no genotypes needed. Pins what the block
# has always done (the messages, the merge order, the dummy expansion) and
# the guards the Phase 3 sweep put in it (#67 duplicate IDs, #78 categorical
# cardinality, #37 an empty include list), so that a caller other than
# fitNULLGLMM (getNeff, LEDGER #100) gets the same behaviour.
#
# Tier: fast.

pp_pheno <- function(n = 30, dup = FALSE) {
  set.seed(7)
  d <- data.frame(
    IID = sprintf("S%02d", seq_len(n)),
    y = rnorm(n), x1 = rnorm(n),
    cat3 = sample(c(1, 2, 3), n, replace = TRUE),
    sex = rep(c(1, 0), length.out = n),
    stringsAsFactors = FALSE
  )
  d$y[3] <- NA # one incomplete row
  if (dup) d <- rbind(d, d[5, ])
  f <- tempfile("pheno_", fileext = ".txt")
  write.table(d, f, quote = FALSE, row.names = FALSE, sep = " ")
  f
}

# genotype sample list: the phenotype IDs in a scrambled order, plus two
# samples with genotypes and no phenotype
pp_geno <- function(n = 30) {
  ids <- c(sprintf("S%02d", rev(seq_len(n))), "G01", "G02")
  data.frame(IIDgeno = ids, IndexGeno = seq_along(ids), stringsAsFactors = FALSE)
}

pp_run <- function(phenoFile = pp_pheno(), covarColList = c("x1", "cat3"),
                   qCovarCol = "cat3", invNormalize = FALSE, ...) {
  suppressMessages(capture.output(
    out <- preparePhenotypeFrame(
      phenoFile, "y", covarColList, qCovarCol, "IID", pp_geno(),
      traitType = "quantitative", invNormalize = invNormalize, ...
    )
  ))
  out
}

test_that("the frame is the genotype-ordered merge with the categorical covariate expanded", {
  out <- pp_run()
  # 29 complete phenotype rows all have genotypes; the two genotype-only samples do not
  expect_identical(nrow(out$data), 29L)
  expect_identical(out$data$IndexGeno, sort(out$data$IndexGeno))
  expect_identical(out$data$IID, rev(sprintf("S%02d", 1:30))[rev(sprintf("S%02d", 1:30)) != "S03"])
  # cat3 became two dummy columns and the formula names them
  expect_setequal(out$covarColList, c("x1", "cat32", "cat33"))
  expect_true(all(c("cat32", "cat33") %in% colnames(out$data)))
  expect_identical(deparse(out$formula.null), "y ~ x1 + cat32 + cat33")
  expect_true(out$hasCovariate)
  # indicator over the genotype list's rows: TRUE except S03 and the two G samples
  g <- pp_geno()
  expect_identical(out$indicatorGenoSamplesWithPheno, !(g$IIDgeno %in% c("S03", "G01", "G02")))
})

test_that("no covariates: intercept-only formula, hasCovariate FALSE", {
  out <- pp_run(covarColList = NULL, qCovarCol = NULL)
  expect_identical(deparse(out$formula.null), "y ~ 1")
  expect_false(out$hasCovariate)
  expect_identical(nrow(out$data), 29L)
})

test_that("LEDGER #67: a duplicated sample ID is refused by name", {
  expect_error(pp_run(phenoFile = pp_pheno(dup = TRUE)),
               "1 sample ID\\(s\\) appear more than once .* S05")
})

test_that("LEDGER #78: a categorical column with one level per sample is refused", {
  f <- pp_pheno()
  d <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  d$cat3 <- seq_len(nrow(d))
  write.table(d, f, quote = FALSE, row.names = FALSE, sep = " ")
  expect_error(pp_run(phenoFile = f), "every sample its own level")
})

test_that("LEDGER #37: an include list matching nobody stops with the intended message", {
  inc <- tempfile("include_", fileext = ".txt")
  writeLines(c("Z1", "Z2"), inc)
  expect_error(pp_run(SampleIDIncludeFile = inc), "No individuals left in the analysis")
})

test_that("a sex-specific subset keeps only that code, and the sex column must exist", {
  out <- pp_run(sexCol = "sex", FemaleOnly = TRUE, FemaleCode = 1)
  expect_identical(nrow(out$data), 14L) # 15 odd-numbered females minus the incomplete S03
  expect_true(all(out$data$sex == 1))
  expect_error(pp_run(sexCol = "nosuch", MaleOnly = TRUE), "column for sex nosuch does not exist")
})

test_that("invNormalize replaces the phenotype by its normal scores", {
  out <- pp_run(invNormalize = TRUE)
  y <- out$data$y
  expect_equal(sort(y), qnorm((seq_along(y) - 0.5) / length(y)))
})

test_that("restricting to the sparse GRM's samples happens before the merge", {
  ids <- tempfile("grm_", fileext = ".txt")
  writeLines(sprintf("S%02d", 1:10), ids)
  out <- pp_run(sparseGRMSampleIDFile = ids, restrictToSparseGRM = TRUE)
  expect_identical(nrow(out$data), 9L) # S01..S10 minus the incomplete S03
})
