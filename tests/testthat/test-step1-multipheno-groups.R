# Two phenotypes with DIFFERENT complete-case sample sets: two ingests, and each
# phenotype still byte-identical to its own standalone run.
#
# Its own file because its four step-1 fits are slow and the tier's wall time is
# bounded below by its slowest file.

MPG_GRM <- paste0("output/sparseGRM_relatednessCutoff_0.125_1000_",
                  "randomMarkersUsed.sparseGRM.mtx")

# The bundled phenotype file with y_quantitative set to NA on every hundredth
# row, so that phenotype has 990 complete cases against y_binary's 1000.
mpg_phenofile <- function() {
  src <- file.path(saige_extdata, "input",
                   "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
  ped <- read.table(src, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE)
  ped$y_quantitative[seq(100, 1000, by = 100)] <- NA
  path <- tempfile(fileext = ".txt")
  write.table(ped, path, row.names = FALSE, quote = FALSE, sep = "\t")
  path
}

mpg_common <- function(phenoFile) {
  c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=", MPG_GRM),
    paste0("--sparseGRMSampleIDFile=", MPG_GRM, ".sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    paste0("--phenoFile=", phenoFile),
    "--covarColList=x1,x2,a9,a10",
    "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID",
    "--IsOverwriteVarianceRatioFile=TRUE",
    "--isCateVarianceRatio=TRUE"
  )
}

mpg_run <- function(phenoFile, extra, prefix) {
  st <- run_saige_cli("step1_fitNULLGLMM.R",
                      c(mpg_common(phenoFile), extra, paste0("--outputPrefix=", prefix)))
  expect_cli_ok(st)
  st
}

test_that("two sample sets mean two ingests, and both fits are standalone-identical", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  phenoFile <- mpg_phenofile()

  multi <- file.path(tempfile("step1_mpg_"), "multi")
  dir.create(dirname(multi), recursive = TRUE)
  st <- mpg_run(phenoFile, c(
    "--phenoCol=y_binary,y_quantitative",
    "--traitType=binary,quantitative",
    "--invNormalize=FALSE,TRUE"
  ), multi)

  lg <- readLines(attr(st, "log"), warn = FALSE)
  # EXACTLY two: one per sample-set group. One would mean the second fit reused
  # a store built over the wrong samples; three would mean the sharing is not
  # working at all.
  expect_identical(sum(grepl("setgeno mark2", lg, fixed = TRUE)), 2L)
  expect_true(any(grepl("2 phenotypes in 2 sample-set group", lg)))

  # Each phenotype against its OWN standalone run on the same file.
  for (spec in list(list(p = "y_binary", args = c("--phenoCol=y_binary", "--traitType=binary")),
                    list(p = "y_quantitative",
                         args = c("--phenoCol=y_quantitative", "--traitType=quantitative",
                                  "--invNormalize=TRUE")))) {
    solo <- file.path(tempfile(paste0("step1_mpg_solo_", spec$p, "_")), "solo")
    dir.create(dirname(solo), recursive = TRUE)
    stSolo <- mpg_run(phenoFile, spec$args, solo)
    # The standalone y_binary run doubles as the legacy-naming check.
    expect_true(file.exists(paste0(solo, ".rda")))
    expect_false(file.exists(paste0(solo, ".multiTraits.txt")))
    expect_identical(sum(grepl("setgeno mark2", readLines(attr(stSolo, "log"),
                                                          warn = FALSE),
                               fixed = TRUE)), 1L)
    for (ext in c(".rda", ".varianceRatio.txt")) {
      expect_step1_matches(paste0(solo, ext),
                           paste0(multi, "_", spec$p, ext))
    }
  }
})
