# LEDGER #68 FIXED: covariate names are backtick-quoted in every formula
# built from them. Bare, a name containing formula syntax was REINTERPRETED:
# `x1^2` vanished from y ~ x1 + x2 + x1^2 (term labels x1, x2), `x1*x2` became
# an interaction, and the --covarColList membership check had already
# confirmed the column present -- so a covariate the user named was silently
# dropped from the null model. Tier: standard (two ~4 s fits + two marker runs).
#
# Fixture: the bundled phenotype file plus one fresh covariate, written twice
# under two names -- `x1^2` (the offending one) and `x3` (a plain one). Same
# numbers, different name, so the two fits must agree on everything but the
# name: that is the reconciliation, and the plain-name run is the control.

covar_pheno <- local({
  paths <- NULL
  function() {
    if (is.null(paths)) {
      src <- file.path(saige_extdata, "input",
                       "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
      d <- utils::read.table(src, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      set.seed(68L)
      extra <- stats::rnorm(nrow(d))
      p_odd <- tempfile("pheno_odd_", fileext = ".txt")
      p_plain <- tempfile("pheno_plain_", fileext = ".txt")
      d$`x1^2` <- extra
      utils::write.table(d, p_odd, quote = FALSE, row.names = FALSE)
      names(d)[names(d) == "x1^2"] <- "x3"
      utils::write.table(d, p_plain, quote = FALSE, row.names = FALSE)
      paths <<- list(odd = p_odd, plain = p_plain)
    }
    paths
  }
})

covar_fit <- function(pheno, covar, prefix) {
  sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  args <- c("--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
            paste0("--sparseGRMFile=", sg), paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            "--useSparseGRMtoFitNULL=TRUE", paste0("--phenoFile=", pheno),
            "--phenoCol=y_binary", paste0("--covarColList=", covar), "--qCovarColList=a9,a10",
            "--sampleIDColinphenoFile=IID", "--traitType=binary",
            "--IsOverwriteVarianceRatioFile=TRUE", "--isCateVarianceRatio=TRUE",
            paste0("--outputPrefix=", prefix))
  st <- run_saige_cli("step1_fitNULLGLMM.R", args)
  expect_cli_ok(st)
  e <- new.env(); load(paste0(prefix, ".rda"), envir = e)
  e$modglmm
}

covar_scan <- function(prefix) {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
    paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", prefix, ".rda"),
    paste0("--varianceRatioFile=", prefix, ".varianceRatio.txt"),
    "--is_noadjCov=FALSE"))   # the covariate-adjusting path reads X
  expect_cli_ok(st)
  out
}

test_that("a covariate named x1^2 is fitted as a covariate, not parsed as a power", {
  skip_unless_level("standard")
  ph <- covar_pheno()
  p_odd <- tempfile("s1_odd_"); p_plain <- tempfile("s1_plain_")
  odd <- covar_fit(ph$odd, "x1,x2,a9,a10,x1^2", p_odd)
  plain <- covar_fit(ph$plain, "x1,x2,a9,a10,x3", p_plain)
  # the named column is in the design: intercept + x1 + x2 + a9 (2 dummies) +
  # a10 (2 dummies) + the extra one = 8 columns, and it was 7 when it vanished
  expect_identical(ncol(odd$X), ncol(plain$X))
  expect_identical(ncol(odd$X), 8L)
  # The stored design's column is called X.x1.2., not x1^2: Covariate_Transform
  # passes the model matrix through data.frame(), whose check.names mangles it.
  # Cosmetic -- the VALUES are the covariate's -- and pinned as such.
  expect_identical(unname(odd$X[, 8]), unname(plain$X[, 8]))
  # same numbers under a different name: the fits agree exactly
  expect_identical(unname(odd$theta), unname(plain$theta))
  expect_identical(unname(odd$fitted.values), unname(plain$fitted.values))
  expect_identical(unname(odd$coefficients), unname(plain$coefficients))
  # and step 2 consumes the odd-named model on the covariate-adjusting path,
  # giving the same bytes as the plain-named one (names never reach the output)
  o_out <- covar_scan(p_odd); p_out <- covar_scan(p_plain)
  expect_identical(readLines(o_out), readLines(p_out))
})
