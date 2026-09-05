# LEDGER #98 FIXED: the #65 rescale touches exactly the model files this run
# wrote, not whatever `<prefix>_chr<j>.rda` happens to be on disk.
#
# The end of fitNULLGLMM's quantitative branch scales the saved model(s) from
# the standardised fit back to phenotype units. It used to find them with a
# glob -- c(modelOut, paste0(outputPrefix, "_chr", 1:22, ".rda")) filtered by
# file.exists -- so a per-chromosome file left by an EARLIER run under the same
# prefix (an --isLowMemLOCO run followed by a plain one, or any leftover) was
# multiplied by this run's s and re-saved, silently. Now the three save sites
# record what they wrote and only those files are rescaled.
#
# Arm 1 (the defect): plant the committed quant_sparseGRM.rda as
#   <prefix>_chr5.rda, run a plain quantitative fit under <prefix>. The planted
#   file must come back byte-identical. At unpatched HEAD it came back with its
#   theta multiplied by s^2 and its coefficients by s.
# Arm 2 (the coverage the fix must not lose): an --isLowMemLOCO fit at y and at
#   y * 100. Its _noLOCO.rda and its _chr<j>.rda files are the ones the old glob
#   DID legitimately reach, so their contents must still scale (theta by 1e4,
#   the per-chromosome LOCO coefficients by 100), which proves the tracked list
#   covers the low-memory layout. No other test drives --isLowMemLOCO. A
#   sparse-GRM fit cannot serve here (fitNULLGLMM forces LOCO off for it), and a
#   dense-GRM fit over the full 22-chromosome fixture costs ~15 s per
#   chromosome, so the arm derives a two-chromosome fixture at run time: the
#   .bed is SNP-major, so the first 1,000 records of chr 1 and of chr 2 (two
#   contiguous byte slices) with the matching bim lines and the same fam is a
#   valid PLINK set. Dense GRM, LOCO over two chromosomes, variance-ratio
#   estimation skipped (not what this arm tests): ~3 s per fit.
#
# Tier: standard. Three step-1 fits, ~15 s.

sm_common_args <- function(scale, prefix) {
  c(paste0("--phenoFile=", ps_pheno_file(scale)),
    "--phenoCol=y_quantitative", "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID", "--traitType=quantitative",
    paste0("--outputPrefix=", prefix))
}

sm_run <- function(args) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  st <- run_saige_cli("step1_fitNULLGLMM.R", args)
  expect_cli_ok(st)
  readLines(attr(st, "log"), warn = FALSE)
}

# the quant_sparseGRM golden's configuration (sparse GRM, so no LOCO)
sm_fit_sparse <- function(scale, prefix) {
  sm_run(c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed",
           ".sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    "--IsOverwriteVarianceRatioFile=TRUE", "--isCateVarianceRatio=TRUE",
    sm_common_args(scale, prefix)
  ))
}

# sm_two_chr_fixture (the two-chromosome PLINK subset) is in helper-pheno-scale.R,
# shared with test-step1-lowmem-loco.R.

# dense GRM over the two-chromosome fixture, LOCO written per chromosome
sm_fit_lowmem <- function(scale, prefix, plink) {
  sm_run(c(
    paste0("--plinkFile=", plink),
    "--nThreads=1", "--minMAFforGRM=0.01", "--LOCO=TRUE", "--isLowMemLOCO=TRUE",
    "--skipVarianceRatioEstimation=TRUE",
    sm_common_args(scale, prefix)
  ))
}

sm_load <- function(path) {
  e <- new.env()
  load(path, e)
  e$modglmm
}

test_that("LEDGER #98 FIXED: a stale _chr<j>.rda under the prefix is left byte-identical", {
  prefix <- file.path(tempfile("step1_stale_"), "fit")
  dir.create(dirname(prefix), recursive = TRUE)
  golden <- file.path(saige_goldens, "step1", "quant_sparseGRM.rda")
  planted <- paste0(prefix, "_chr5.rda")
  file.copy(golden, planted)
  log <- sm_fit_sparse(1, prefix)
  expect_true(file.exists(paste0(prefix, ".rda")))
  # the run rescaled its own file and names it ...
  expect_true(any(grepl("^Saved model\\(s\\) scaled back to phenotype units .*: fit\\.rda$", log)))
  # ... and did not touch the planted one
  expect_files_byte_identical(golden, planted)
})

test_that("LEDGER #98: the low-memory LOCO files are still rescaled (tracked, not globbed)", {
  skip_unless_level("standard")
  plink <- file.path(tempfile("twochr_"), "sub")
  dir.create(dirname(plink), recursive = TRUE)
  sm_two_chr_fixture(plink)
  fits <- lapply(c(1, 100), function(scale) {
    prefix <- file.path(tempfile(sprintf("step1_lowmem%g_", scale)), "fit")
    dir.create(dirname(prefix), recursive = TRUE)
    log <- sm_fit_lowmem(scale, prefix, plink)
    chr <- sort(list.files(dirname(prefix), "^fit_chr[0-9]+\\.rda$", full.names = TRUE))
    expect_length(chr, 2L)
    expect_true(file.exists(paste0(prefix, "_noLOCO.rda")))
    expect_false(file.exists(paste0(prefix, ".rda")))
    # every file the run wrote is named in its rescale line
    line <- grep("^Saved model\\(s\\) scaled back to phenotype units", log, value = TRUE)
    expect_length(line, 1L)
    for (f in c(paste0(prefix, "_noLOCO.rda"), chr)) {
      expect_true(grepl(basename(f), line, fixed = TRUE), info = basename(f))
    }
    list(noLOCO = sm_load(paste0(prefix, "_noLOCO.rda")),
         chr = setNames(lapply(chr, sm_load), basename(chr)))
  })
  base <- fits[[1]]
  x100 <- fits[[2]]
  expect_lt(rel_err(x100$noLOCO$theta, base$noLOCO$theta * 100^2), 1e-8)
  expect_lt(rel_err(x100$noLOCO$coefficients, base$noLOCO$coefficients * 100), 1e-8)
  expect_identical(names(x100$chr), names(base$chr))
  for (nm in names(base$chr)) {
    j <- as.integer(sub("^fit_chr([0-9]+)\\.rda$", "\\1", nm))
    r0 <- base$chr[[nm]]$LOCOResult[[j]]
    r1 <- x100$chr[[nm]]$LOCOResult[[j]]
    expect_true(isTRUE(r0$isLOCO), info = nm)
    expect_lt(rel_err(r1$coefficients, r0$coefficients * 100), 1e-8)
    expect_lt(rel_err(r1$fitted.values, r0$fitted.values * 100), 1e-8)
    expect_lt(rel_err(r1$cov, r0$cov * 100^2), 1e-8)
    # LEDGER #99 FIXED: the LOCO entry's obj.noK is rebuilt from the model's y
    # and theta in phenotype units, ONCE. Closed form within each file --
    # S_a = X'(y - mu_LOCO) and V = 1/theta[1] -- and across the two scales
    # (S_a by 100, V by 1e-4). At unpatched HEAD y was scaled a second time on
    # its way into the entry, so S_a was X'(s y - mu) instead of X'(y - mu).
    for (m in list(base$chr[[nm]], x100$chr[[nm]])) {
      r <- m$LOCOResult[[j]]
      res <- as.vector(m$y - r$fitted.values) # y is stored as a one-column matrix
      expect_lt(rel_err(r$obj.noK$S_a, colSums(m$X * res)), 1e-8)
      expect_lt(rel_err(r$obj.noK$V, rep(1 / m$theta[1], length(m$y))), 1e-10)
    }
    expect_lt(rel_err(r1$obj.noK$S_a, r0$obj.noK$S_a * 100), 1e-8)
    expect_lt(rel_err(r1$obj.noK$V, r0$obj.noK$V / 100^2), 1e-8)
  }
})
