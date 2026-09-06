# Multi-model single-variant scanning: the manifest, the sample-set grouping,
# and the per-model null-model preparation that SPAGMMATtest used to do inline
# for the one model it could have.
#
# The design in one line: models whose `sampleID` vectors are `identical()`
# share ONE genotype reader, ONE decode and ONE QC pass per marker; models with
# different sample sets get their own reader. Nothing is unioned and nothing is
# gathered, so each model's rows are byte-identical to the rows a single-model
# run would have written -- which is what makes the goldens the gate for this
# feature rather than a tolerance.
#
# Scope, decided and not to be re-opened without new evidence: the marker path
# only. Region tests, the male-X rewrite and `--subSampleFile` stop for P > 1
# (see SPAGMMATtest's validation block), because each of them keeps per-run
# state that is model-order-dependent.

# Read the 3-column model manifest: name, .rda, variance-ratio file. No header,
# whitespace-separated, everything character (a phenotype called `NA` or `T`
# must not become a logical, and a numeric-looking name must not lose its
# formatting).
readModelManifest <- function(path) {
  if (!file.exists(path)) {
    stop("GMMATmodel_varianceRatio_multiTraits_File ", path, " does not exist\n")
  }
  # na.strings = character(0) as well as colClasses = "character": without it a
  # phenotype literally named NA becomes a missing value rather than a name, and
  # so does the literal NA that step 1 writes in column 3 when it skipped
  # variance-ratio estimation -- which has its own message below.
  manifest <- utils::read.table(path,
    header = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", comment.char = "", na.strings = character(0)
  )
  if (ncol(manifest) != 3) {
    stop(
      "GMMATmodel_varianceRatio_multiTraits_File ", path, " has ", ncol(manifest),
      " columns; it must have exactly 3: name, model .rda, variance ratio file\n"
    )
  }
  if (nrow(manifest) < 1) {
    stop("GMMATmodel_varianceRatio_multiTraits_File ", path, " has no rows\n")
  }
  names(manifest) <- c("name", "modelFile", "varianceRatioFile")

  if (any(is.na(manifest$name)) || any(manifest$name == "")) {
    stop("Every row of ", path, " must have a non-empty name in column 1\n")
  }
  if (anyDuplicated(manifest$name) != 0) {
    dup <- unique(manifest$name[duplicated(manifest$name)])
    stop(
      "Names in ", path, " must be unique; repeated: ",
      paste(dup, collapse = ", "), "\n"
    )
  }
  # The name becomes a filename suffix, so it must not contain a path separator.
  bad <- grep("/", manifest$name, fixed = TRUE)
  if (length(bad) > 0) {
    stop(
      "Names in ", path, " become output-file suffixes and must not contain '/'; ",
      "row ", paste(bad, collapse = ", "), "\n"
    )
  }

  for (i in seq_len(nrow(manifest))) {
    if (is.na(manifest$modelFile[i]) || manifest$modelFile[i] == "") {
      stop("Row ", i, " of ", path, " has no model file in column 2\n")
    }
    if (!file.exists(manifest$modelFile[i])) {
      stop(
        "Row ", i, " of ", path, " (name ", manifest$name[i],
        "): model file ", manifest$modelFile[i], " does not exist\n"
      )
    }
    if (is.na(manifest$varianceRatioFile[i]) || manifest$varianceRatioFile[i] == "") {
      stop("Row ", i, " of ", path, " has no variance ratio file in column 3\n")
    }
    # Step 1 writes a literal NA into the manifest when it skipped variance-ratio
    # estimation. That model cannot be single-variant tested at all, so say why
    # rather than letting Get_Variance_Ratio fail on a missing file.
    if (manifest$varianceRatioFile[i] == "NA") {
      stop(
        "Row ", i, " of ", path, " (name ", manifest$name[i],
        "): step 1 skipped variance-ratio estimation; single-variant tests need one\n"
      )
    }
    if (!file.exists(manifest$varianceRatioFile[i])) {
      stop(
        "Row ", i, " of ", path, " (name ", manifest$name[i],
        "): variance ratio file ", manifest$varianceRatioFile[i], " does not exist\n"
      )
    }
  }
  manifest
}

# Turn the two mutually exclusive ways of naming models into ONE list of specs.
#
# The suffix rule is upstream's and is deliberate: a manifest with a SINGLE row
# produces an unsuffixed output file, so a 1-row manifest run is byte- AND
# name-identical to the equivalent --GMMATmodelFile run.
resolveModelSpecs <- function(GMMATmodelFile, varianceRatioFile, manifestFile) {
  if (manifestFile != "" && GMMATmodelFile != "") {
    stop(
      "--GMMATmodel_varianceRatio_multiTraits_File and --GMMATmodelFile are ",
      "mutually exclusive: give the manifest, or the single model, not both\n"
    )
  }
  if (manifestFile == "") {
    return(list(list(
      name = "", modelFile = GMMATmodelFile,
      varianceRatioFile = varianceRatioFile, suffix = ""
    )))
  }
  manifest <- readModelManifest(manifestFile)
  nModels <- nrow(manifest)
  lapply(seq_len(nModels), function(i) {
    list(
      name = manifest$name[i],
      modelFile = manifest$modelFile[i],
      varianceRatioFile = manifest$varianceRatioFile[i],
      suffix = if (nModels > 1) paste0("_", manifest$name[i]) else ""
    )
  })
}

# Which models share a sample set? Group k against the representatives seen so
# far, in manifest order, so group indices are stable and group 1 is model 1's.
# identical() rather than setequal(): the ORDER of sampleID is the reader's
# output order, so two models agreeing as sets but not as vectors are two
# groups, and each still gets bytes identical to its standalone run.
groupModelsBySampleSet <- function(sampleIDs) {
  groupOf <- integer(length(sampleIDs))
  reps <- list()
  for (k in seq_along(sampleIDs)) {
    g <- 0L
    for (j in seq_along(reps)) {
      if (identical(reps[[j]], sampleIDs[[k]])) {
        g <- j
        break
      }
    }
    if (g == 0L) {
      reps[[length(reps) + 1L]] <- sampleIDs[[k]]
      g <- length(reps)
    }
    groupOf[k] <- g
  }
  list(groupOf = groupOf, representatives = reps)
}

# Build the reader for a sample-set group after the first. The first group's
# reader comes from setGenoInput exactly as it always did; this reuses the paths
# setGenoInput resolved so a bim/fam default derived there is not re-derived
# here (and cannot drift from it).
addGenoReaderForGroup <- function(objGeno, sampleInModel) {
  files <- objGeno$files
  if (is.null(files)) {
    stop("setGenoInput did not record the resolved genotype file paths\n")
  }
  if (objGeno$genoType == "plink") {
    spec <- list(
      bimFile = files$bimFile, famFile = files$famFile,
      bedFile = files$bedFile, AlleleOrder = files$AlleleOrder
    )
  } else if (objGeno$genoType == "pgen") {
    spec <- list(
      pgenFile = files$pgenFile, pvarFile = files$pvarFile,
      psamFile = files$psamFile
    )
  } else {
    stop("Unsupported genotype type for a second sample-set group: ", objGeno$genoType, "\n")
  }
  addGenoReaderInCPP(objGeno$genoType, spec, sampleInModel)
}

# All the models' restart state must agree, because there is one marker loop and
# it resumes at one chunk. v1 stops rather than resuming at the minimum: a
# disagreement means the files were not all written by the same interrupted run,
# and guessing which one is authoritative is how rows get silently duplicated
# (LEDGER #47's shape).
restartChunkForModels <- function(outChecks, OutputFileIndex) {
  idx <- vapply(outChecks, function(x) as.integer(x$indexChunk), integer(1))
  if (length(unique(idx)) != 1L) {
    stop(
      "The models' progress index files disagree about where to restart, so ",
      "the scan cannot resume: ",
      paste(paste0(OutputFileIndex, " -> chunk ", idx), collapse = "; "),
      ". Delete them and re-run, or set --is_overwrite_output=TRUE\n"
    )
  }
  idx[1]
}

# One model's null-model preparation: everything SPAGMMATtest used to do inline
# between ReadModel and setSAIGEobjInCPP, per model. The body is that code
# verbatim; what is new is that it returns rather than assigning into the
# caller's frame, and that it says which model it is working on so a stop is
# attributable to a manifest row.
#
# Note the sparse Sigma is built PER MODEL and there are P copies of it. That is
# deliberate and not an oversight: Sigma depends on the model's theta and mu2,
# so it is not shareable even between two models over the same samples.
prepareNullModel <- function(spec, k, chrom, LOCO, is_Firth_beta, subSampleFile,
                             sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff,
                             cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude,
                             isGroupTest, is_fastTest, pval_cutoff_for_fastTest,
                             is_noadjCov, max_MAC_use_ER, nThreads) {
  label <- if (spec$name == "") "" else paste0(" ", spec$name)
  cat("---- null model ", k, label, " ----\n", sep = "")

  if (subSampleFile == "") {
    obj.model <- ReadModel(spec$modelFile, chrom, LOCO, is_Firth_beta) # readInGLMM.R
  } else {
    # readInGLMM.R
    obj.model <- ReadModel_subsample(
      spec$modelFile, chrom, LOCO, is_Firth_beta, subSampleFile
    )
  }

  if (obj.model$traitType == "binary") {
    if (max_MAC_use_ER > 0) {
      cat(
        "P-values of genetic variants with MAC <= ", max_MAC_use_ER,
        " will be calculated via effecient resampling.\n"
      )
      if (max_MAC_use_ER > 4) {
        cat("WARNING: Efficient resampling may not work well for MAC > 4!")
      }
    }
  }

  isSparseGRM <- TRUE
  if (sparseGRMFile != "") {
    sparseSigmaRList <- setSparseSigma_new(
      sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff,
      obj.model$sampleID, obj.model$theta, obj.model$mu2, obj.model$traitType
    )
    isSparseGRM <- TRUE
  } else {
    sparseSigmaRList <- list(
      nSubj = 0, locations = matrix(0, nrow = 2, ncol = 2), values = rep(0, 2)
    )
    isSparseGRM <- FALSE
  }

  cat("isSparseGRM is ", isSparseGRM, "\n")
  # readInGLMM.R
  ratioVecList <- Get_Variance_Ratio(
    spec$varianceRatioFile, cateVarRatioMinMACVecExclude,
    cateVarRatioMaxMACVecInclude, isGroupTest, isSparseGRM
  )

  # LEDGER #10 fix: single-variant analysis dereferences the null_noXadj
  # variance ratios (the no-covariate-adjustment fast pass); a file without
  # such rows used to crash in C++ mid-analysis (Mat index out of bounds).
  # `&` not `&&`, and the same below: these conditions are moved VERBATIM from
  # SPAGMMATtest, where they have always been vectorised operators on
  # length-1 values. Changing them here would be a behaviour change smuggled in
  # as a style fix, which the bug-for-bug rule forbids outside a fix commit.
  noXadjMissing <- length(ratioVecList$ratioVec_null_noXadj) == 0
  if (!isGroupTest & noXadjMissing) { # nolint: vector_logic_linter.
    stop(
      "varianceRatioFile ", spec$varianceRatioFile,
      " contains no 'null_noXadj' rows, ",
      "which single-variant tests require. This file was produced by an older ",
      "SAIGE version or configuration; please re-run Step 1.\n"
    )
  }

  # A gap in the categorical variance-ratio bins means some MAC falls in no
  # category. This guard was written against LEDGER #32, where such a marker
  # silently reused the PREVIOUS marker's ratio and so made a parallel run
  # differ from a serial one; S3 FIXED that -- selectVarianceRatio throws now --
  # so the divergence is gone and both arms stop at the same marker. It is kept
  # only because it fires before any marker is read and names the offending
  # flags, where the C++ throw names a MAC. See R/MarkerScan.R for why it stays
  # gated on a parallel scan rather than firing always.
  hasGap <- varianceRatioBinsHaveGap(
    cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude
  )
  manyRatios <- length(ratioVecList$ratioVec_null) > 1
  if (nThreads > 1 & manyRatios & hasGap) { # nolint: vector_logic_linter.
    stop(
      "cateVarRatioMinMACVecExclude/cateVarRatioMaxMACVecInclude leave a gap in MAC ",
      "coverage, so some markers fall in no variance-ratio category at all. ",
      "Such a marker has no ratio to use and the test will stop when it is ",
      "reached, whatever --nThreads is. Give contiguous MAC categories.\n"
    )
  }

  if (is_fastTest) {
    if (!file.exists(spec$varianceRatioFile)) {
      is_fastTest <- FALSE
      cat("No variance ratio file is specified, so is_fastTest is not working.\n")
    }
    if (pval_cutoff_for_fastTest == 0) {
      is_fastTest <- FALSE
      cat("pval_cutoff_for_fastTest == 0, so is_fastTest is not working.\n")
    }
  }

  if (!isSparseGRM & !is_noadjCov) { # nolint: vector_logic_linter.
    is_fastTest <- FALSE
    cat("isSparseGRM is FALSE and is_noadjCov is FALSE,, so is_fastTest is not working.\n")
  }

  if (is_fastTest) {
    if (isSparseGRM) {
      cat("is_fastTest is TRUE.\n")
      if (ratioVecList$ratioVec_null[1] == -1) {
        stop(
          "Variance ratios estimated without GRM are not found, so the fast ",
          "tests can't be performed.\n Please set is_fastTest=FALSE or re-run ",
          "the variance ratio estimation in Step 1 with --skipModelFitting=TRUE ",
          "using the most recent version of the program.\n"
        )
      } else {
        cat("The fast tests will be performed (when p-values >= ", pval_cutoff_for_fastTest, ").\n")
      }
    }
  }

  nsample <- length(obj.model$y)

  list(
    obj.model = obj.model,
    sparseSigmaRList = sparseSigmaRList,
    ratioVecList = ratioVecList,
    isSparseGRM = isSparseGRM,
    is_fastTest = is_fastTest,
    # The trailing sample size is the last categorical bin's upper bound and is
    # therefore per model, not per run.
    cateVarRatioMaxMACVecInclude = c(cateVarRatioMaxMACVecInclude, nsample)
  )
}

# Push the prepared models into C++: setSAIGEobjInCPP for the first (which
# clears the previous run's models) and addSAIGEobjInCPP for the rest, each with
# the zero-based index of its sample-set group. For P == 1 this is exactly the
# single setSAIGEobjInCPP call SPAGMMATtest always made.
pushNullModelsToCPP <- function(prepared, groupOf, SPAcutoff, impute_method,
                                is_noadjCov, pval_cutoff_for_fastTest,
                                isCondition, condition_genoIndex,
                                is_Firth_beta, pCutoffforFirth,
                                firth_solver,
                                cateVarRatioMinMACVecExclude) {
  for (k in seq_along(prepared)) {
    p <- prepared[[k]]
    obj.model <- p$obj.model
    args <- list(
      t_XVX = obj.model$obj.noK$XVX,
      t_XXVX_inv = obj.model$obj.noK$XXVX_inv,
      t_XV = obj.model$obj.noK$XV,
      t_XVX_inv_XV = obj.model$obj.noK$XVX_inv_XV,
      t_Sigma_iXXSigma_iX = obj.model$Sigma_iXXSigma_iX,
      t_X = obj.model$X,
      t_S_a = obj.model$obj.noK$S_a,
      t_res = obj.model$residuals,
      t_mu2 = obj.model$mu2,
      t_mu = obj.model$mu,
      t_varRatio_sparse = as.vector(p$ratioVecList$ratioVec_sparse),
      t_varRatio_null = as.vector(p$ratioVecList$ratioVec_null),
      t_varRatio_null_noXadj = as.vector(p$ratioVecList$ratioVec_null_noXadj),
      t_cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
      t_cateVarRatioMaxMACVecInclude = p$cateVarRatioMaxMACVecInclude,
      t_SPA_Cutoff = SPAcutoff,
      t_tauvec = obj.model$theta,
      t_traitType = obj.model$traitType,
      t_y = obj.model$y,
      t_impute_method = impute_method,
      t_flagSparseGRM = p$isSparseGRM,
      t_isFastTest = p$is_fastTest,
      t_isnoadjCov = is_noadjCov,
      t_pval_cutoff_for_fastTest = pval_cutoff_for_fastTest,
      t_locationMat = as.matrix(p$sparseSigmaRList$locations),
      t_valueVec = p$sparseSigmaRList$values,
      t_dimNum = p$sparseSigmaRList$nSubj,
      t_isCondition = isCondition,
      t_condition_genoIndex = condition_genoIndex,
      t_is_Firth_beta = is_Firth_beta,
      t_pCutoffforFirth = pCutoffforFirth,
      t_firthSolver = firth_solver,
      t_offset = obj.model$offset,
      t_resout = as.integer(obj.model$obj_cc$res.out)
    )
    if (k == 1L) {
      do.call(setSAIGEobjInCPP, args)
    } else {
      do.call(addSAIGEobjInCPP, c(args, list(t_groupIndex = as.integer(groupOf[k] - 1L))))
    }
  }
  invisible(NULL)
}
