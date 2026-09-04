# Step 1 for several phenotypes in one process, sharing the genotype ingest.
#
# The ingest (`setgeno`) reads and bit-packs the whole GRM marker set into a C++
# global. It is O(M x N) and it depends on NOTHING about the phenotype except
# which SAMPLES are in the fit -- so two phenotypes with the same complete-case
# sample set can share one, and two with different sets cannot.
#
# The fits themselves stay SEQUENTIAL and untouched: each phenotype runs the
# same `fitNULLGLMM` body it would run alone, over the already-loaded store.
# Cross-phenotype batched multi-RHS PCG -- the win on the DENSE-GRM path, where
# the fit is 300-450 GRM passes against one for the ingest -- is deliberately
# out of scope here.
#
# **Why the outputs are byte-identical to standalone runs.** The only
# RNG-dependent input to the ingest is the `arma::randi(1000)` variance-ratio
# marker pool, drawn from R's RNG inside setgeno. The driver calls
# `set.seed(seed)` before every fit, so the FIRST fit of a group draws exactly
# the pool a standalone run would; later fits in the group reuse that same
# store, so they get the same pool by construction. Nothing between the seed and
# the ingest consumes RNG (`fread`, `glm`, `checkPerfectSep`,
# `Covariate_Transform`, `getsubGRM`, and `invNormalize`'s rank(ties="average")
# are all deterministic), `GetTrace` re-seeds 200 on entry, and
# `extractVarianceRatio` re-seeds seed+1 before its own `sample()`s.

# The loaded genotype store. One per process, because the C++ side has one.
.step1Geno <- new.env(parent = emptyenv())
.step1Geno$loaded <- FALSE
.step1Geno$key <- NULL
.step1Geno$nIngests <- 0L

# What makes two ingests interchangeable: the same files, the same sample
# subset, and the same QC/variance-ratio gates -- the gates because
# setminMAFforGRM / setmaxMissingRateforGRM / setminMAC_VarianceRatio are pushed
# BEFORE setgeno and decide which markers it keeps and which it stages for the
# variance ratio.
step1IngestKey <- function(genoIn, subSampleInGeno, indicatorGenoSamplesWithPheno,
                           memoryChunk, isDiagofKinSetAsOne, ingestGates) {
  list(
    genoIn = genoIn,
    subSampleInGeno = subSampleInGeno,
    indicator = indicatorGenoSamplesWithPheno,
    memoryChunk = memoryChunk,
    isDiagofKinSetAsOne = isDiagofKinSetAsOne,
    gates = ingestGates
  )
}

# Load the genotype store, or say that the one already loaded will do.
#
# `identical()` on the whole key, deliberately: anything that is not provably
# the same ingest gets a fresh one, so reuse can only ever be a saving and never
# a difference. A re-ingest is safe because setGenoObjFromSource resets every
# accumulator on entry (LEDGER #55); before that fix this function would have
# been silently wrong on its second call.
ensureGenoLoaded <- function(genoIn, subSampleInGeno, indicatorGenoSamplesWithPheno,
                             memoryChunk, isDiagofKinSetAsOne, ingestGates = NULL) {
  key <- step1IngestKey(
    genoIn, subSampleInGeno, indicatorGenoSamplesWithPheno,
    memoryChunk, isDiagofKinSetAsOne, ingestGates
  )
  if (.step1Geno$loaded && identical(key, .step1Geno$key)) {
    cat("The genotypes for this sample set are already loaded; reusing the ingest\n")
    return(invisible(FALSE))
  }
  if (.step1Geno$loaded) {
    closeGenoFile_plink()
    .step1Geno$loaded <- FALSE
  }
  re1 <- system.time({
    setgenoFromInput(
      genoIn, subSampleInGeno, indicatorGenoSamplesWithPheno,
      memoryChunk, isDiagofKinSetAsOne
    )
  })
  cat("genotype ingest (setgeno) took\n")
  print(re1)
  .step1Geno$loaded <- TRUE
  .step1Geno$key <- key
  .step1Geno$nIngests <- .step1Geno$nIngests + 1L
  invisible(TRUE)
}

releaseGeno <- function() {
  if (.step1Geno$loaded) {
    closeGenoFile_plink()
    .step1Geno$loaded <- FALSE
    .step1Geno$key <- NULL
  }
  invisible(NULL)
}
