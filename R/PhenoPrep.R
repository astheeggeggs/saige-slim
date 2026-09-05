# The phenotype frame a null-model fit is built on, shared by fitNULLGLMM and
# getNeff (which used to carry an older copy of this block: no categorical
# covariates, no duplicate-ID check). Reads the phenotype file, checks the
# named columns exist, refuses a duplicated sample ID (LEDGER #67), keeps
# complete cases, applies SampleIDIncludeFile, turns qCovarCol into factors
# after a cardinality check (LEDGER #78), applies a sex-specific subset,
# expands the design with model.matrix (so qCovarCol dummies become the
# covariates), merges with the genotype sample list `sampleListwithGeno`
# (columns IIDgeno, IndexGeno) in genotype order, optionally restricts to the
# sparse GRM's samples first, and inverse-normalises a quantitative
# phenotype on request. The messages are the ones step 1 has always printed.
#
# Returns data (the merged frame, genotype order), formula.null, covarColList
# (rewritten to the dummy names when qCovarCol is used), hasCovariate and
# indicatorGenoSamplesWithPheno (over sampleListwithGeno's rows).
preparePhenotypeFrame <- function(phenoFile, phenoCol, covarColList, qCovarCol,
                                  sampleIDColinphenoFile, sampleListwithGeno,
                                  traitType, invNormalize,
                                  sexCol = "", FemaleOnly = FALSE, FemaleCode = "1",
                                  MaleOnly = FALSE, MaleCode = "0",
                                  SampleIDIncludeFile = "",
                                  sparseGRMSampleIDFile = NULL,
                                  restrictToSparseGRM = FALSE) {
  if (!file.exists(phenoFile)) {
    stop("ERROR! phenoFile ", phenoFile, " does not exsit\n")
  }
  if (grepl(".gz$", phenoFile) | grepl(".bgz$", phenoFile)) {
    data <- data.table:::fread(
      cmd = paste0(
        "gunzip -c ",
        phenoFile
      ), header = TRUE, stringsAsFactors = FALSE,
      colClasses = list(character = sampleIDColinphenoFile), data.table = FALSE
    )
  } else {
    data <- data.table:::fread(phenoFile,
      header = TRUE, stringsAsFactors = FALSE,
      colClasses = list(character = sampleIDColinphenoFile), data.table = FALSE
    )
  }
  # data = data.frame(ydat)
  for (i in c(phenoCol, covarColList, sampleIDColinphenoFile)) {
    if (!(i %in% colnames(data))) {
      stop("ERROR! column for ", i, " does not exist in the phenoFile \n")
    }
  }
  # LEDGER #67 FIXED: a repeated sample ID multiplied that sample's rows
  # through the genotype merge (model matrix, phenotype vector and the
  # genotype index all doubled; counted in the allele-frequency denominator
  # but not the numerator) and "N+1 samples will be used" was printed.
  dupIDs <- unique(data[[sampleIDColinphenoFile]][duplicated(data[[sampleIDColinphenoFile]])])
  if (length(dupIDs) > 0) {
    stop("ERROR! ", length(dupIDs), " sample ID(s) appear more than once in the phenoFile column ",
         sampleIDColinphenoFile, ": ", paste(head(dupIDs, 10), collapse = ", "),
         if (length(dupIDs) > 10) ", ..." else "", ". Each sample must have exactly one row.\n")
  }


  # LEDGER #77 FIXED: the sex column is kept through this subset when a
  # sex-specific run asks for it, and its existence is checked against the
  # FILE here, before the subset -- the check below used to test the subset,
  # so --FemaleOnly/--MaleOnly failed unless the sex column happened to be a
  # covariate, blaming a phenotype file that was correct.
  keepCols <- c(phenoCol, covarColList, sampleIDColinphenoFile)
  if (FemaleOnly | MaleOnly) {
    if (!sexCol %in% colnames(data)) {
      stop("ERROR! column for sex ", sexCol, " does not exist in the phenoFile \n")
    }
    keepCols <- c(keepCols, sexCol)
  }
  data <- data[, which(colnames(data) %in% keepCols), drop = FALSE]

  data <- data[complete.cases(data), , drop = FALSE]

  if (SampleIDIncludeFile != "") {
    if (!file.exists(SampleIDIncludeFile)) {
      stop("ERROR! SampleIDIncludeFile ", SampleIDIncludeFile, " does not exsit\n")
    } else {
      sampleIDInclude <- data.table:::fread(SampleIDIncludeFile,
        header = FALSE, stringsAsFactors = FALSE,
        colClasses = c("character"), data.table = FALSE
      )
      sampleIDInclude <- as.vector(sampleIDInclude[!duplicated(sampleIDInclude), ])
      cat(length(sampleIDInclude), " non-duplicated sample IDs were found in SampleIDIncludeFile\n")
      idCol <- as.vector(data[, which(colnames(data) == sampleIDColinphenoFile)])
      data <- data[which(idCol %in% sampleIDInclude), , drop = FALSE]
      cat(nrow(data), " samples in sampleIDInclude have non-missing phenotypes and covariates\n")
      # LEDGER #37 FIXED: say so HERE. Downstream, model.matrix on a 0-row
      # frame with factor covariates died in `contrasts<-` and the seq()
      # before the intended stop() errored with "wrong sign in 'by'".
      if (nrow(data) == 0) {
        stop("No individuals left in the analysis. Exiting... (no sample in ", SampleIDIncludeFile,
             " has non-missing phenotypes and covariates)\n")
      }
    }
  }


  if (length(qCovarCol) > 0) {
    cat(qCovarCol, "are categorical covariates\n")
    if (!all(qCovarCol %in% covarColList)) {
      stop("ERROR! all covariates in qCovarCol must be in covarColList\n")
    } else {
      for (q in qCovarCol) {
        # LEDGER #78 FIXED: a column named categorical was made a factor with
        # no look at its cardinality, so a continuous column became one dummy
        # column per distinct value (upstream: 43,790 levels, 216 GB, no
        # output). A column in which every sample has its own value cannot be
        # a category and is refused; one with more than 20 levels is
        # accepted with a warning that says how many dummies it will cost.
        nLevels <- length(unique(data[, q]))
        if (nLevels >= nrow(data) && nrow(data) > 1) {
          stop("ERROR! covariate ", q, " is named in qCovarColList as categorical but has ",
               nLevels, " distinct values in ", nrow(data), " samples -- every sample its own ",
               "level. It is not a category; remove it from qCovarColList to treat it as ",
               "continuous.\n")
        }
        if (nLevels > 20) {
          warning("covariate ", q, " is categorical with ", nLevels, " levels, i.e. ",
                  nLevels - 1, " dummy columns in the null model. Check that it is really ",
                  "categorical.", call. = FALSE)
        }
        data[, q] <- as.factor(data[, q])
      }
    }
  }


  if (FemaleOnly | MaleOnly) {
    if (!sexCol %in% colnames(data)) {
      stop("ERROR! column for sex ", sexCol, " does not exist in the phenoFile \n")
    } else {
      if (FemaleOnly) {
        data <- data[which(data[, which(colnames(data) ==
          sexCol)] == FemaleCode), ]
        if (nrow(data) == 0) {
          stop(
            "ERROR! no samples in the phenotype are coded as ",
            FemaleCode, " in the column ", sexCol,
            "\n"
          )
        }
      } else if (MaleOnly) {
        data <- data[which(data[, which(colnames(data) ==
          sexCol)] == MaleCode), ]
        if (nrow(data) == 0) {
          stop(
            "ERROR! no samples in the phenotype are coded as ",
            MaleCode, " in the column ", sexCol, "\n"
          )
        }
      }
    }
  }


  if (length(covarColList) > 0) {
    formula <- paste0(phenoCol, "~", paste0(quoteCovarNames(covarColList), collapse = "+"))
    hasCovariate <- TRUE
  } else {
    formula <- paste0(phenoCol, "~ 1")
    hasCovariate <- FALSE
  }

  cat("formula is ", formula, "\n")
  formula.null <- as.formula(formula)
  mmat <- model.matrix(formula.null, data, na.action = NULL)
  mmat <- data.frame(mmat)
  mmat <- cbind(mmat, data[, which(colnames(data) == phenoCol), drop = FALSE])
  colnames(mmat)[ncol(mmat)] <- phenoCol


  if (length(covarColList) > 0) {
    if (length(qCovarCol) > 0) {
      covarColList <- colnames(mmat)[2:(ncol(mmat) - 1)]
      formula <- paste0(phenoCol, "~", paste0(quoteCovarNames(covarColList), collapse = "+"))
      formula.null <- as.formula(formula)
    }
  }

  mmat$IID <- data[, which(sampleIDColinphenoFile == colnames(data))]
  mmat_nomissing <- mmat[complete.cases(mmat), , drop = FALSE]
  # LEDGER #37 FIXED: seq(1, 0, by = 1) errored here before the stop() below could
  mmat_nomissing$IndexPheno <- seq_len(nrow(mmat_nomissing))
  cat(nrow(mmat_nomissing), " samples have non-missing phenotypes\n")
  if (nrow(mmat_nomissing) == 0) {
    stop("No individuals left in the analysis. Exiting...\n")
  }

  if (restrictToSparseGRM) {
    sampleListwithGenov0 <- data.table:::fread(sparseGRMSampleIDFile,
      header = FALSE, , colClasses = c("character"), data.table = FALSE
    )
    colnames(sampleListwithGenov0) <- c("IIDgeno")
    cat(length(sampleListwithGenov0$IIDgeno), " samples are in the sparse GRM\n")
    mmat_nomissing <- mmat_nomissing[which(mmat_nomissing$IID %in% sampleListwithGenov0$IIDgeno), ]
    cat(nrow(mmat_nomissing),
        " samples who have non-missing phenotypes are also in the sparse GRM\n")
  }


  dataMerge <- merge(mmat_nomissing, sampleListwithGeno,
    by.x = "IID", by.y = "IIDgeno"
  )

  dataMerge_sort <- dataMerge[with(dataMerge, order(IndexGeno)), ]

  rm(mmat)
  rm(mmat_nomissing)
  gc()
  indicatorGenoSamplesWithPheno <- (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)

  if (nrow(dataMerge_sort) < nrow(sampleListwithGeno)) {
    cat(
      nrow(sampleListwithGeno) - nrow(dataMerge_sort),
      " samples in geno file do not have phenotypes\n"
    )
  }
  cat(nrow(dataMerge_sort), " samples will be used for analysis\n")

  if (traitType == "quantitative" & invNormalize) {
    cat(
      "Perform the inverse nomalization for ", phenoCol,
      "\n"
    )
    invPheno <- qnorm((rank(dataMerge_sort[, which(colnames(dataMerge_sort) ==
      phenoCol)], na.last = "keep") - 0.5) / sum(!is.na(dataMerge_sort[
      ,
      which(colnames(dataMerge_sort) == phenoCol)
    ])))
    dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] <- invPheno
  }
  list(
    data = dataMerge_sort,
    formula.null = formula.null,
    covarColList = covarColList,
    hasCovariate = hasCovariate,
    indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno
  )
}
