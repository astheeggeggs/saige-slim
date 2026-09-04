# Step 1 for several phenotypes in ONE process, sharing the genotype ingest.
#
# The claim is byte-identity: each phenotype's `.rda` and `.varianceRatio.txt`
# are exactly what a standalone run writes. That is legitimate here for the same
# reason it is for `--nThreads`: this is ONE code path (`fitNULLGLMM`) re-run per
# phenotype, over a genotype store whose contents do not depend on which
# phenotype asked for it. The only RNG-dependent input to the ingest is the
# variance-ratio marker pool, and the driver re-seeds before every fit.
#
# The two step-1 goldens used here were fitted on the same plink file, the same
# sparse GRM, the same covariates and the same sample ID column, so ONE
# multi-phenotype command reproduces both.

MP1_GRM <- paste0("output/sparseGRM_relatednessCutoff_0.125_1000_",
                  "randomMarkersUsed.sparseGRM.mtx")
MP1_COMMON <- function() {
  c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=", MP1_GRM),
    paste0("--sparseGRMSampleIDFile=", MP1_GRM, ".sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    "--covarColList=x1,x2,a9,a10",
    "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID",
    "--IsOverwriteVarianceRatioFile=TRUE",
    "--isCateVarianceRatio=TRUE"
  )
}

mp1_run <- function(extra, prefix) {
  st <- run_saige_cli("step1_fitNULLGLMM.R",
                      c(MP1_COMMON(), extra, paste0("--outputPrefix=", prefix)))
  expect_cli_ok(st)
  st
}

mp1_log <- function(status) readLines(attr(status, "log"), warn = FALSE)

# The per-ingest marker setgeno prints. Counting it is what separates "the
# outputs are right" from "the ingest was actually shared" -- exactly the
# LEDGER #30 distinction, and drill D4's target.
mp1_ingests <- function(status) sum(grepl("setgeno mark2", mp1_log(status), fixed = TRUE))

test_that("two phenotypes in one run reproduce both step-1 goldens", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step1_mp_"), "multi")
  dir.create(dirname(out), recursive = TRUE)
  st <- mp1_run(c(
    "--phenoCol=y_binary,y_quantitative",
    "--traitType=binary,quantitative",
    "--invNormalize=FALSE,TRUE"
  ), out)

  for (pair in list(c("y_binary", "binary_sparseGRM_a"),
                    c("y_quantitative", "quant_sparseGRM"))) {
    for (ext in c(".rda", ".varianceRatio.txt")) {
      expect_step1_matches(
        file.path(saige_goldens, "step1", paste0(pair[2], ext)),
        paste0(out, "_", pair[1], ext)
      )
    }
  }

  # ONE ingest for both phenotypes: they have the same complete-case sample set.
  expect_identical(mp1_ingests(st), 1L)
  expect_true(any(grepl("already loaded; reusing the ingest", mp1_log(st))))

  # The manifest, in --phenoCol order, ready for step 2.
  man <- read.table(paste0(out, ".multiTraits.txt"), header = FALSE,
                    stringsAsFactors = FALSE, colClasses = "character")
  expect_identical(nrow(man), 2L)
  expect_identical(man[[1]], c("y_binary", "y_quantitative"))
  expect_identical(man[[2]], paste0(out, "_", c("y_binary", "y_quantitative"), ".rda"))
  expect_true(all(file.exists(man[[3]])))
})

test_that("the phenotype ORDER does not change any phenotype's output", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # The gate for "no state leaks between fits": if anything survived one fit
  # into the next, running them the other way round would move something.
  out <- file.path(tempfile("step1_mprev_"), "rev")
  dir.create(dirname(out), recursive = TRUE)
  mp1_run(c(
    "--phenoCol=y_quantitative,y_binary",
    "--traitType=quantitative,binary",
    "--invNormalize=TRUE,FALSE"
  ), out)

  for (pair in list(c("y_binary", "binary_sparseGRM_a"),
                    c("y_quantitative", "quant_sparseGRM"))) {
    for (ext in c(".rda", ".varianceRatio.txt")) {
      expect_step1_matches(
        file.path(saige_goldens, "step1", paste0(pair[2], ext)),
        paste0(out, "_", pair[1], ext)
      )
    }
  }
  man <- read.table(paste0(out, ".multiTraits.txt"), header = FALSE,
                    stringsAsFactors = FALSE, colClasses = "character")
  expect_identical(man[[1]], c("y_quantitative", "y_binary"))
})

test_that("a single phenotype keeps the legacy names and writes no manifest", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # The negative control for the whole feature: same driver, one phenotype,
  # unsuffixed output, no manifest.
  out <- file.path(tempfile("step1_mp1_"), "one")
  dir.create(dirname(out), recursive = TRUE)
  st <- mp1_run(c("--phenoCol=y_binary", "--traitType=binary"), out)
  expect_true(file.exists(paste0(out, ".rda")))
  expect_false(file.exists(paste0(out, "_y_binary.rda")))
  expect_false(file.exists(paste0(out, ".multiTraits.txt")))
  expect_identical(mp1_ingests(st), 1L)
  for (ext in c(".rda", ".varianceRatio.txt")) {
    expect_step1_matches(
      file.path(saige_goldens, "step1", paste0("binary_sparseGRM_a", ext)),
      paste0(out, ext)
    )
  }
})

test_that("recyclePerPheno and the driver's validation", {
  expect_identical(SAIGE:::recyclePerPheno(c("a", "b"), 2L, "x"), c("a", "b"))
  expect_identical(SAIGE:::recyclePerPheno("a", 3L, "x"), rep("a", 3))
  expect_error(SAIGE:::recyclePerPheno(c("a", "b"), 3L, "traitType"),
               "traitType has 2 values but there are 3 phenotypes")

  expect_error(SAIGE::fitNULLGLMM_multi(phenoCol = character(0)),
               "at least one phenotype")
  expect_error(SAIGE::fitNULLGLMM_multi(phenoCol = c("a", "a")),
               "repeated names")
  expect_error(SAIGE::fitNULLGLMM_multi(phenoCol = c("a", "")),
               "non-empty")
  expect_error(SAIGE::fitNULLGLMM_multi(phenoCol = c("a", "b"), traitType = "survival"),
               "must be binary or quantitative")
})
