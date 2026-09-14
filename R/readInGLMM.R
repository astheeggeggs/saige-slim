getChromNumber <- function(chrom = "") {
  if (chrom == "") {
    stop("chrom is not specified\n")
  } else {
    chrom_v2 <- as.character(chrom)
    chrom_v2 <- gsub("CHR", "", chrom_v2, ignore.case = TRUE)
    chrom_v3 <- as.numeric(gsub("[^0-9.]", "", chrom_v2))

    if (!is.na(chrom_v3)) {
      if (chrom_v3 <= 22 & chrom_v3 >= 1) {
        # stop("chromosome ", chrom, " is out of the range of null model LOCO results\n")
        # }else {
        cat("Leave chromosome ", chrom_v3, " out will be applied\n")
      } else {
        chrom_v3 <- chrom
      }
    } else {
      chrom_v3 <- chrom
    }
  }

  chrom_v3
}

removeLOCOResult <- function(chromList, obj.glmm.null) {
  for (chr in 1:22) {
    if (chr %in% chromList) {
      obj.glmm.null$LOCOResult[chr] <- list(NULL)
      cat("chromosome ", chr, " model results are removed to save memory\n")
      gc()
    }
  }
  obj.glmm.null
}


ReadModel <- function(GMMATmodelFile = "", chrom = "", LOCO = TRUE, is_Firth_beta = FALSE) {
  # Check file existence
  Check_File_Exist(GMMATmodelFile, "GMMATmodelFile")
  if (!LOCO %in% c(TRUE, FALSE)) {
    stop("LOCO should be TRUE or FALSE.")
  }
  # load GMMATmodelFile
  load(GMMATmodelFile)
  # if(!modglmm$LOCO){
  #     if(LOCO){
  #             stop("LOCO is TRUE but the null model was not fit with LOCO=TRUE\n")
  #             LOCO=FALSE
  #     }
  # }

  modglmm$Y <- NULL
  # modglmm$offset = modglmm$linear.predictors - modglmm$coefficients[1]
  modglmm$linear.predictors <- NULL
  modglmm$coefficients <- NULL
  modglmm$cov <- NULL
  # rm(modglmm)
  gc()
  # traitType = modglmm$traitType
  # y = modglmm$y
  # X = modglmm$X
  # N = length(y)
  # tauVec = modglmm$theta
  # indChromCheck = FALSE
  chrom_v3 <- NULL

  if (!LOCO) {
    print("Leave-one-chromosome-out is not applied")
    if (modglmm$LOCO) {
      for (chr in 1:22) {
        modglmm$LOCOResult[chr] <- list(NULL)
        cat("chromosome ", chr, " model results are removed to save memory\n")
        gc()
      }
    }
  } else {
    if (!modglmm$LOCO) {
      stop(
        "LOCO is TRUE but this null model has no LOCO results. Step 1 stores them only ",
        "when it was run with LOCO=TRUE against a full GRM. A model fitted with ",
        "--useSparseGRMtoFitNULL=TRUE never has them: a sparse GRM has no per-chromosome ",
        "structure to leave out, so Step 1 refuses that combination. Set LOCO=FALSE in ",
        "this step (Step 2), or refit the null model in Step 1 with ",
        "--useSparseGRMtoFitNULL=FALSE --LOCO=TRUE.\n"
      )
    } else {
      if (chrom == "") {
        stop(
          "--chrom must be specified when LOCO is TRUE: the null model holds one fit ",
          "per chromosome and --chrom selects which to use. This applies to every ",
          "step-2 test, single-variant as well as region-based. Either pass --chrom, ",
          "or set LOCO=FALSE in this step (Step 2).\n"
        )
      } else {
        chrom_v3 <- getChromNumber(chrom)
      }
    }

    if (is.numeric(chrom_v3)) {
      if (chrom_v3 >= 1 & chrom_v3 <= 22) {
        chromList <- c(1:22)
        chromList <- chromList[which(chromList != chrom_v3)]
        modglmm <- removeLOCOResult(chromList, modglmm)
        if (!is.null(modglmm$LOCOResult[[chrom_v3]])) {
          modglmm$fitted.values <- modglmm$LOCOResult[[chrom_v3]]$fitted.values
          modglmm$residuals <- modglmm$LOCOResult[[chrom_v3]]$residuals
          modglmm$obj.noK <- modglmm$LOCOResult[[chrom_v3]]$obj.noK
          if (is_Firth_beta) {
            if (!is.null(modglmm$LOCOResult[[chrom_v3]]$offset)) {
              modglmm$offset <- modglmm$LOCOResult[[chrom_v3]]$offset
            }
          }
        }
        # modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$linear.predictors -  modglmm$LOCOResult[[chrom_v3]]$coefficients[1]
        modglmm$LOCOResult[chrom_v3] <- list(NULL)
      } else { # if(chrom_v3 >= 1 & chrom_v3 <= 22){
        chromList <- c(1:22)
        modglmm <- removeLOCOResult(chromList, modglmm)
      }
    } else {
      chromList <- c(1:22)
      modglmm <- removeLOCOResult(chromList, modglmm)
    }
    gc()
  }

  modglmm$mu <- as.vector(modglmm$fitted.values)
  tau <- modglmm$theta
  N <- length(modglmm$mu)
  if (modglmm$traitType == "binary") {
    modglmm$mu2 <- (modglmm$mu) * (1 - modglmm$mu)
    modglmm$obj_cc <- SKAT::SKAT_Null_Model(modglmm$y ~ modglmm$X - 1, out_type = "D", Adjustment = FALSE)
    modglmm$obj_cc$mu <- modglmm$mu
    modglmm$obj_cc$res <- modglmm$res
    modglmm$obj_cc$pi_1 <- modglmm$mu2
  } else if (modglmm$traitType == "quantitative") {
    modglmm$mu2 <- (1 / tau[1]) * rep(1, N)
  }
  # if(FALSE){

  if (is_Firth_beta) {
    if (modglmm$traitType == "binary") {
      if (is.null(modglmm$offset)) {
        # if(FALSE){
        # print("WARNING. is_Firth_beta = TRUE, but offset was not computed in Step 1. Please re-run Step 1 using the more rencet version of SAIGE/SAIGE-GENE.")
        cat("Applying is_Firth_beta = TRUE and note the results would be more accurate withe Step 1 results using the more rencet version of SAIGE/SAIGE-GENE. \n")

        if (ncol(modglmm$X) == 1) {
          covoffset <- rep(0, nrow(modglmm$X))
        } else {
          formulastr <- paste0("y ~ ", paste(quoteCovarNames(colnames(modglmm$X)[-1]), collapse = "+"))
          modglmm$X <- cbind(modglmm$X, as.vector(modglmm$y))
          colnames(modglmm$X)[ncol(modglmm$X)] <- "y"
          modglmm$X <- as.data.frame(modglmm$X, check.names = FALSE) # LEDGER #68: keep names the quoted formula can find
          formula.new <- as.formula(formulastr)
          modwitcov <- glm(formula.new, data = modglmm$X, family = binomial)
          modglmm$X <- modglmm$X[, -ncol(modglmm$X)]
          covoffset <- as.matrix(modglmm$X[, -1]) %*% modwitcov$coefficients[-1]
          modglmm$X <- as.matrix(modglmm$X)
        }
        modglmm$offset <- covoffset
        # }
      }
    }
  }


  if (is.null(modglmm$offset)) {
    modglmm$offset <- rep(0, nrow(modglmm$X))
  }


  if (is.null(modglmm$Sigma_iXXSigma_iX)) {
    modglmm$Sigma_iXXSigma_iX <- matrix(0, nrow = 1, ncol = 1)
  }
  modglmm
}


ReadModel_subsample <- function(GMMATmodelFile = "", chrom = "", LOCO = TRUE, is_Firth_beta = FALSE, subSampleFile = "") {
  # Check file existence
  Check_File_Exist(GMMATmodelFile, "GMMATmodelFile")
  if (!LOCO %in% c(TRUE, FALSE)) {
    stop("LOCO should be TRUE or FALSE.")
  }

  if (!file.exists(subSampleFile)) {
    stop("ERROR! ", subSampleFile, " does not exist\n")
  } else {
    cat("subSampleFile is specified. This option is used when any sample included in Step 1 but does not have dosages/genotypes for Step 2. Please make sure it contains one column of sample IDs that will be used for subsetting samples from the Step 1 results for Step 2 jobs\n")
    subSampleID0 <- data.table::fread(subSampleFile, header = FALSE, colClasses = c("character"), data.table = FALSE)
    subSampleID0 <- subSampleID0[, 1]
  }

  load(GMMATmodelFile)

  includeIndex <- which(modglmm$sampleID %in% subSampleID0)
  if (length(includeIndex) == 0) {
    stop("No samples analyzed in Step 1 are included in ", subSampleFile, "\n")
  } else {
    cat(length(includeIndex), " samples analyzed in Step 1 will be included for Step 2\n")

    percNum <- length(includeIndex) / (length(modglmm$sampleID))
    percNum <- percNum * 100
    if (percNum < 90) {
      cat("WARNING: ", percNum, "% samples analyzed in Step 1 will be included for Step 2. This proportion is small and using subsetted Step 1 results could lead to incorrect Step 2 results. Please re-run Step 1 for the sample subset\n")
    }
  }


  modglmm$Y <- NULL
  # modglmm$offset = modglmm$linear.predictors - modglmm$coefficients[1]
  modglmm$linear.predictors <- NULL
  modglmm$coefficients <- NULL
  modglmm$cov <- NULL
  modglmm$obj.glm.null <- NULL
  modglmm$fitted.values <- modglmm$fitted.values[includeIndex]
  modglmm$Y <- modglmm$Y[includeIndex]
  modglmm$residuals <- modglmm$residuals[includeIndex]
  modglmm$sampleID <- modglmm$sampleID[includeIndex]
  modglmm$y <- modglmm$y[includeIndex]
  modglmm$X <- modglmm$X[includeIndex, , drop = FALSE]


  # rm(modglmm)
  gc()
  # traitType = modglmm$traitType
  # y = modglmm$y
  # X = modglmm$X
  # N = length(y)
  # tauVec = modglmm$theta
  # indChromCheck = FALSE
  chrom_v3 <- NULL

  if (!LOCO) {
    print("Leave-one-chromosome-out is not applied")
    if (modglmm$LOCO) {
      for (chr in 1:22) {
        modglmm$LOCOResult[chr] <- list(NULL)
        cat("chromosome ", chr, " model results are removed to save memory\n")
        gc()
      }
    }
  } else {
    if (!modglmm$LOCO) {
      stop(
        "LOCO is TRUE but this null model has no LOCO results. Step 1 stores them only ",
        "when it was run with LOCO=TRUE against a full GRM. A model fitted with ",
        "--useSparseGRMtoFitNULL=TRUE never has them: a sparse GRM has no per-chromosome ",
        "structure to leave out, so Step 1 refuses that combination. Set LOCO=FALSE in ",
        "this step (Step 2), or refit the null model in Step 1 with ",
        "--useSparseGRMtoFitNULL=FALSE --LOCO=TRUE.\n"
      )
    } else {
      if (chrom == "") {
        stop(
          "--chrom must be specified when LOCO is TRUE: the null model holds one fit ",
          "per chromosome and --chrom selects which to use. This applies to every ",
          "step-2 test, single-variant as well as region-based. Either pass --chrom, ",
          "or set LOCO=FALSE in this step (Step 2).\n"
        )
      } else {
        chrom_v3 <- getChromNumber(chrom)
      }
    }

    if (is.numeric(chrom_v3)) {
      if (chrom_v3 >= 1 & chrom_v3 <= 22) {
        chromList <- c(1:22)
        chromList <- chromList[which(chromList != chrom_v3)]
        modglmm <- removeLOCOResult(chromList, modglmm)
        if (!is.null(modglmm$LOCOResult[[chrom_v3]])) {
          # LEDGER #14 FIXED: the LOCO replacement is subset here, like the
          # non-LOCO values were at load, so mu below is built from a vector
          # that is ALREADY the subset in both branches.
          modglmm$fitted.values <- modglmm$LOCOResult[[chrom_v3]]$fitted.values[includeIndex]
          modglmm$residuals <- modglmm$LOCOResult[[chrom_v3]]$residuals
          modglmm$residuals <- modglmm$residuals[includeIndex]
          modglmm$obj.noK <- modglmm$LOCOResult[[chrom_v3]]$obj.noK


          if (is_Firth_beta) {
            if (!is.null(modglmm$LOCOResult[[chrom_v3]]$offset)) {
              modglmm$offset <- modglmm$LOCOResult[[chrom_v3]]$offset
            }
          }
        }
        # modglmm$offset = modglmm$LOCOResult[[chrom_v3]]$linear.predictors -  modglmm$LOCOResult[[chrom_v3]]$coefficients[1]
        modglmm$LOCOResult[chrom_v3] <- list(NULL)
      } else { # if(chrom_v3 >= 1 & chrom_v3 <= 22){
        chromList <- c(1:22)
        modglmm <- removeLOCOResult(chromList, modglmm)
      }
    } else {
      chromList <- c(1:22)
      modglmm <- removeLOCOResult(chromList, modglmm)
    }
    gc()
  }

  # LEDGER #14 FIXED: fitted.values is already the subset (at load for the
  # non-LOCO path, above for LOCO), so it was being indexed TWICE here -- a
  # no-op for a prefix subset, NAs and wrong elements for any other.
  modglmm$mu <- as.vector(modglmm$fitted.values)


  tau <- modglmm$theta
  N <- length(modglmm$mu)
  if (modglmm$traitType == "binary") {
    modglmm$mu2 <- (modglmm$mu) * (1 - modglmm$mu)
    modglmm$obj_cc <- SKAT::SKAT_Null_Model(modglmm$y ~ modglmm$X - 1, out_type = "D", Adjustment = FALSE)
    modglmm$obj_cc$mu <- modglmm$mu

    if (!is.null(modglmm$res)) {
      modglmm$res <- modglmm$res[includeIndex]
    }
    modglmm$obj_cc$res <- modglmm$res
    modglmm$obj_cc$pi_1 <- modglmm$mu2
  } else if (modglmm$traitType == "quantitative") {
    modglmm$mu2 <- (1 / tau[1]) * rep(1, N)
  }
  # if(FALSE){

  if (is_Firth_beta) {
    if (modglmm$traitType == "binary") {
      if (is.null(modglmm$offset)) {
        # if(FALSE){
        # print("WARNING. is_Firth_beta = TRUE, but offset was not computed in Step 1. Please re-run Step 1 using the more rencet version of SAIGE/SAIGE-GENE.")
        cat("Applying is_Firth_beta = TRUE and note the results would be more accurate withe Step 1 results using the more rencet version of SAIGE/SAIGE-GENE. \n")

        if (ncol(modglmm$X) == 1) {
          covoffset <- rep(0, nrow(modglmm$X))
        } else {
          formulastr <- paste0("y ~ ", paste(quoteCovarNames(colnames(modglmm$X)[-1]), collapse = "+"))
          modglmm$X <- cbind(modglmm$X, as.vector(modglmm$y))
          colnames(modglmm$X)[ncol(modglmm$X)] <- "y"
          modglmm$X <- as.data.frame(modglmm$X, check.names = FALSE) # LEDGER #68: keep names the quoted formula can find
          formula.new <- as.formula(formulastr)
          modwitcov <- glm(formula.new, data = modglmm$X, family = binomial)
          modglmm$X <- modglmm$X[, -ncol(modglmm$X)]
          covoffset <- as.matrix(modglmm$X[, -1]) %*% modwitcov$coefficients[-1]
          modglmm$X <- as.matrix(modglmm$X)
        }
        modglmm$offset <- covoffset
        # }
      }
    }
  }


  obj.noK <- ScoreTest_NULL_Model(modglmm$mu, modglmm$mu2, modglmm$y, modglmm$X)
  modglmm$obj.noK <- obj.noK
  rm(obj.noK)


  if (is.null(modglmm$offset)) {
    modglmm$offset <- rep(0, nrow(modglmm$X))
  } else {
    if (length(modglmm$offset) > N) {
      modglmm$offset <- modglmm$offset[includeIndex]
    }
  }


  if (is.null(modglmm$Sigma_iXXSigma_iX)) {
    modglmm$Sigma_iXXSigma_iX <- matrix(0, nrow = 1, ncol = 1)
  } else {
    modglmm$Sigma_iXXSigma_iX <- modglmm$Sigma_iXXSigma_iX[includeIndex, , drop = FALSE]
  }
  modglmm
}


# isGroupTest and useSparseGRMtoFitNULL were formals until 2026-09-14 (C6):
# neither appeared in the body, and the production caller passed five
# arguments, so useSparseGRMtoFitNULL was a never-forced missing promise (the
# LEDGER #85 shape). test-unit-varianceratio.R pins the formals.
Get_Variance_Ratio <- function(varianceRatioFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isSparseGRM) {
  iscateVR <- FALSE
  # check variance ratio
  if (!file.exists(varianceRatioFile)) {
    if (varianceRatioFile != "") {
      stop("varianceRatioFile is specified but the file ", varianceRatioFile, " does not exist\n")
    } else {
      # LEDGER #70 FIXED: this used to be a bare cat(). A missing ratio scales
      # var by 1/r, BETA by r and SE by sqrt(r) in EVERY row (5-26% on the
      # bundled fixtures, anti-conservative where r > 1), so it is raised as an
      # R warning -- visible at the end of an Rscript run and in the log -- and
      # says what is substituted. Still not a stop(): a ratio of 1 is what a
      # step-1 run that skipped variance-ratio estimation intends.
      warning("varianceRatioFile is not specified, so a variance ratio of 1 is used ",
              "for every marker: var is scaled by 1/r, BETA by r and SE by sqrt(r) ",
              "relative to a run with the step-1 ratio. Pass --varianceRatioFile ",
              "unless step 1 deliberately skipped the ratio.", call. = FALSE)
      cat("varianceRatioFile is not specified so variance ratio won't be used\n")
    }
    if (isSparseGRM) {
      cat("WARNING: Sparse GRM is specified. Please make sure the null model was fit using the same sparse GRM in Step 1.\n")
      ratioVec_sparse <- c(1)
      ratioVec_null <- c(-1)
      ratioVec_null_noXadj <- c(-1)
    } else {
      ratioVec_sparse <- c(-1)
      ratioVec_null <- c(1)
      ratioVec_null_noXadj <- c(1)
    }
  } else {
    varRatioData <- data.frame(data.table:::fread(varianceRatioFile, header = FALSE, stringsAsFactors = FALSE))
    if (ncol(varRatioData) == 3) {
      spindex <- which(varRatioData[, 2] == "sparse")
      if (length(spindex) > 0) {
        ratioVec_sparse <- varRatioData[which(varRatioData[, 2] == "sparse"), 1]
        ratioVec_sparse <- as.numeric(ratioVec_sparse)
        # if(!isSparseGRM & sum(ratioVec_sparse > 1.0001 | ratioVec_sparse < 0.9999) > 0){
        #       	stop("sparse GRM is not specified but it was used for estimating variance ratios in Step 1. Please specify --sparseGRMFile and --sparseGRMSampleIDFile\n")
        # }
      } else {
        ratioVec_sparse <- c(-1)
        if (isSparseGRM) {
          stop("sparse GRM is specified but the variance ratio for sparse GRM was not estimatedin Step 1. Pleae remove --sparseGRMFile and --sparseGRMSampleIDFile\n")
        }
      }
      ratioVec_null <- varRatioData[which(varRatioData[, 2] == "null"), 1]
      ratioVec_null <- as.numeric(ratioVec_null)
      # Refuse a 3-column file with no `null` row (audit E1, 2026-09-14): every
      # step-2 test selects its ratio from the null rows for a common marker
      # even under a sparse GRM, and the build is -DARMA_NO_DEBUG, so an empty
      # vector in C++ is read off the end silently, not caught.
      if (length(ratioVec_null) == 0) {
        stop("varianceRatioFile ", varianceRatioFile, " has no 'null' row (types present: ",
             paste(unique(varRatioData[, 2]), collapse = ","), "). Every step-2 test reads ",
             "the null variance ratio for at least the common markers; re-run Step 1, or ",
             "restore the null rows.\n")
      }
      cat("variance Ratio null is ", ratioVec_null, "\n")
      if (length(ratioVec_null) > 1) {
        iscateVR <- TRUE
        nrv <- length(ratioVec_null)
      } else {
        nrv <- 1
      }

      ratioVec_null_noXadj <- varRatioData[which(varRatioData[, 2] == "null_noXadj"), 1]
      # LEDGER #39 FIXED: the coercion and the cardinality guard sat in the
      # EMPTY branch, so they never ran; they run on the populated one now,
      # like the sibling ratios above.
      if (length(ratioVec_null_noXadj) > 0) {
        ratioVec_null_noXadj <- as.numeric(ratioVec_null_noXadj)
        cat("variance Ratio null_noXadj is ", ratioVec_null_noXadj, "\n")
        if (length(ratioVec_null_noXadj) == 1 & iscateVR) {
          stop("categorical variance Ratio null_noXadj are needed\n")
        }
      }
    } else {
      # LEDGER #13 fix: this legacy branch could never complete
      # (ratioVec_null_noXadj was never assigned). Reject cleanly instead.
      stop(
        "varianceRatioFile ", varianceRatioFile, " is in the pre-1.0.6 legacy format (",
        ncol(varRatioData), " column(s)); this version requires the 3-column ",
        "'value type category' format. Please re-run Step 1.\n"
      )
      # The 19 lines that followed this stop() -- two WARNING strings describing
      # the legacy format the branch refuses -- were unreachable since the
      # LEDGER #13 fix and were deleted 2026-09-14 (C6).
    }


    if (iscateVR) {
      ln <- length(cateVarRatioMinMACVecExclude)
      hn <- length(cateVarRatioMaxMACVecInclude)

      if (nrv != ln) {
        stop("ERROR! The number of variance ratios are different from the length of cateVarRatioMinMACVecExclude\n")
      }
      if (ln != (hn + 1)) {
        stop("ERROR! The length of cateVarRatioMaxMACVecInclude does not match with the lenght of cateVarRatioMinMACVecExclude (-1)\n")
      }
    }
  }

  print("ratioVec_sparse")
  print(ratioVec_sparse)
  print("ratioVec_null")
  print(ratioVec_null)
  print("ratioVec_null_noXadj")
  print(ratioVec_null_noXadj)

  list(ratioVec_sparse = ratioVec_sparse, ratioVec_null = ratioVec_null, ratioVec_null_noXadj = ratioVec_null_noXadj)
}
