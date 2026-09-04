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
# Per-phenotype option handling: length P, or length 1 to mean "the same for
# every phenotype". Anything else is a user error worth naming.
recyclePerPheno <- function(x, nPheno, name) {
  if (length(x) == nPheno) {
    return(x)
  }
  if (length(x) == 1L) {
    return(rep(x, nPheno))
  }
  stop(
    name, " has ", length(x), " values but there are ", nPheno,
    " phenotypes; give one value per phenotype, or one for all of them\n"
  )
}

# Which phenotypes share a complete-case sample set?
#
# ORDERING ONLY. This decides the order the fits run in, so that phenotypes over
# the same samples are adjacent and the ingest between them can be reused.
# Correctness does not rest on it: `ensureGenoLoaded` compares the actual ingest
# key, so a wrong guess here costs an extra ingest and changes no output.
step1PhenoGroups <- function(phenoFile, phenoCol, covarColList, sampleIDColinphenoFile,
                             SampleIDIncludeFile = "", sexCol = "",
                             FemaleOnly = FALSE, MaleOnly = FALSE,
                             FemaleCode = 1, MaleCode = 0) {
  ped <- data.table::fread(phenoFile, header = TRUE, stringsAsFactors = FALSE,
                           data.table = FALSE)
  needed <- unique(c(covarColList, sampleIDColinphenoFile,
                     if (nzchar(sexCol)) sexCol else NULL))
  needed <- needed[nzchar(needed)]
  missingCols <- setdiff(needed, colnames(ped))
  if (length(missingCols) > 0) {
    stop("Column(s) ", paste(missingCols, collapse = ", "), " are not in ", phenoFile, "\n")
  }
  keepIDs <- NULL
  if (nzchar(SampleIDIncludeFile)) {
    keepIDs <- data.table::fread(SampleIDIncludeFile, header = FALSE,
                                 colClasses = "character", data.table = FALSE)[, 1]
  }
  sexKeep <- rep(TRUE, nrow(ped))
  if (nzchar(sexCol)) {
    if (FemaleOnly) sexKeep <- ped[[sexCol]] == FemaleCode
    if (MaleOnly) sexKeep <- ped[[sexCol]] == MaleCode
    sexKeep[is.na(sexKeep)] <- FALSE
  }
  ids <- lapply(phenoCol, function(p) {
    if (!(p %in% colnames(ped))) {
      stop("Phenotype column ", p, " is not in ", phenoFile, "\n")
    }
    ok <- stats::complete.cases(ped[, c(p, needed), drop = FALSE]) & sexKeep
    keep <- as.character(ped[[sampleIDColinphenoFile]][ok])
    if (!is.null(keepIDs)) keep <- keep[keep %in% keepIDs]
    sort(keep)
  })
  groupOf <- integer(length(ids))
  reps <- list()
  for (k in seq_along(ids)) {
    g <- 0L
    for (j in seq_along(reps)) {
      if (identical(reps[[j]], ids[[k]])) {
        g <- j
        break
      }
    }
    if (g == 0L) {
      reps[[length(reps) + 1L]] <- ids[[k]]
      g <- length(reps)
    }
    groupOf[k] <- g
  }
  list(groupOf = groupOf, nGroups = length(reps), sampleIDs = ids)
}

#' Fit null GLMMs for one or more phenotypes, sharing the genotype ingest
#'
#' A thin driver over \code{fitNULLGLMM}. With a single phenotype it IS
#' \code{fitNULLGLMM}, called once with the legacy output names and no manifest,
#' so nothing about a single-phenotype run changes. With several it groups them
#' by complete-case sample set, runs the fits sequentially with one
#' \code{setgeno} per group, names each phenotype's outputs
#' \code{<outputPrefix>_<phenoCol>}, and writes a manifest
#' \code{<outputPrefix>.multiTraits.txt} that step 2's
#' \code{--GMMATmodel_varianceRatio_multiTraits_File} consumes directly.
#'
#' @param phenoCol character vector of phenotype column names.
#' @param traitType character, one per phenotype or one for all.
#' @param invNormalize logical, one per phenotype or one for all.
#' @param outputPrefix character. Per-phenotype outputs are suffixed with the
#'   phenotype name when there is more than one.
#' @param outputPrefix_varRatio character. Suffixed the same way.
#' @param seed integer. Re-applied before EVERY fit, which is what makes each
#'   phenotype's variance-ratio marker pool identical to a standalone run's.
#' @param ... every other \code{fitNULLGLMM} argument, shared by all phenotypes.
#' @return invisibly, a data frame of the manifest rows.
#' @export
fitNULLGLMM_multi <- function(phenoCol,
                              traitType = "binary",
                              invNormalize = FALSE,
                              outputPrefix = "",
                              outputPrefix_varRatio = NULL,
                              seed = 1,
                              ...) {
  phenoCol <- as.character(phenoCol)
  nPheno <- length(phenoCol)
  if (nPheno < 1L) {
    stop("phenoCol must name at least one phenotype\n")
  }
  if (any(is.na(phenoCol)) || any(!nzchar(phenoCol))) {
    stop("Every phenotype name must be non-empty\n")
  }
  if (anyDuplicated(phenoCol) != 0) {
    dup <- unique(phenoCol[duplicated(phenoCol)])
    stop("phenoCol has repeated names: ", paste(dup, collapse = ", "), "\n")
  }
  traitType <- recyclePerPheno(as.character(traitType), nPheno, "traitType")
  bad <- setdiff(unique(traitType), c("binary", "quantitative"))
  if (length(bad) > 0) {
    stop("traitType must be binary or quantitative; got ",
         paste(bad, collapse = ", "), "\n")
  }
  invNormalize <- recyclePerPheno(as.logical(invNormalize), nPheno, "invNormalize")
  if (anyNA(invNormalize)) {
    stop("invNormalize must be TRUE or FALSE for every phenotype\n")
  }
  if (is.null(outputPrefix_varRatio) || !nzchar(outputPrefix_varRatio[1])) {
    outputPrefix_varRatio <- outputPrefix
  }

  dots <- list(...)

  runOne <- function(k, prefix, prefixVR, releaseAtEnd) {
    set.seed(seed)
    do.call(fitNULLGLMM, c(
      list(
        phenoCol = phenoCol[k], traitType = traitType[k],
        invNormalize = invNormalize[k],
        outputPrefix = prefix, outputPrefix_varRatio = prefixVR,
        seed = seed, releaseGenoAtEnd = releaseAtEnd
      ),
      dots
    ))
  }

  if (nPheno == 1L) {
    # Legacy naming, no manifest, one fit: byte- and name-identical to calling
    # fitNULLGLMM directly, which is what the CLI used to do.
    runOne(1L, outputPrefix, outputPrefix_varRatio, TRUE)
    return(invisible(NULL))
  }

  order <- seq_len(nPheno)
  nGroups <- NA_integer_
  groupOf <- rep(1L, nPheno)
  if (!is.null(dots$phenoFile) && !is.null(dots$covarColList) &&
        !is.null(dots$sampleIDColinphenoFile)) {
    groups <- step1PhenoGroups(
      dots$phenoFile, phenoCol, dots$covarColList, dots$sampleIDColinphenoFile,
      SampleIDIncludeFile = if (is.null(dots$SampleIDIncludeFile)) "" else dots$SampleIDIncludeFile,
      sexCol = if (is.null(dots$sexCol)) "" else dots$sexCol,
      FemaleOnly = isTRUE(dots$FemaleOnly), MaleOnly = isTRUE(dots$MaleOnly),
      FemaleCode = if (is.null(dots$FemaleCode)) 1 else dots$FemaleCode,
      MaleCode = if (is.null(dots$MaleCode)) 0 else dots$MaleCode
    )
    groupOf <- groups$groupOf
    nGroups <- groups$nGroups
    order <- base::order(groupOf, seq_len(nPheno))
  }
  cat(nPheno, "phenotypes in", nGroups, "sample-set group(s)\n")

  results <- vector("list", nPheno)
  for (pos in seq_along(order)) {
    k <- order[pos]
    # Release only after the LAST phenotype of a group, so the next group's
    # ingest replaces it and the last one leaves the process clean.
    isLastOfGroup <- (pos == length(order)) || (groupOf[order[pos + 1L]] != groupOf[k])
    cat("---- step 1 for phenotype", k, phenoCol[k], "----\n")
    results[[k]] <- runOne(
      k, paste0(outputPrefix, "_", phenoCol[k]),
      paste0(outputPrefix_varRatio, "_", phenoCol[k]), isLastOfGroup
    )
  }
  releaseGeno()

  # The manifest, in --phenoCol order rather than in the order the fits ran, so
  # it reads like the command line. Paths come from what fitNULLGLMM actually
  # wrote (the _noLOCO.rda reassignment under --isLowMemLOCO, and the in-place
  # _FemaleOnly/_MaleOnly prefix mutation, both happen inside it).
  manifest <- data.frame(
    name = phenoCol,
    modelFile = vapply(results, function(r) r$modelOut, character(1)),
    varianceRatioFile = vapply(
      results, function(r) if (is.null(r$varRatioFile)) "NA" else r$varRatioFile,
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  manifestFile <- paste0(outputPrefix, ".multiTraits.txt")
  utils::write.table(manifest, manifestFile,
                     row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Wrote the step-2 model manifest to", manifestFile, "\n")
  invisible(manifest)
}
