# Parse/canonicalize the weights.beta argument into the Kb x 2 numeric matrix
# of beta(a, b) parameter pairs (one row per pair). Accepted inputs:
#  - a single string "a1,b1;a2,b2;..." -- semicolon-separated pairs, comma
#    within a pair (the CLI format; "1,25" is the single-pair legacy form)
#  - a numeric vector of length 2 (the legacy R API form, e.g. c(1, 25))
#  - a Kb x 2 numeric matrix (passed through after validation)
# Each pair must be exactly 2 finite numerics with a > 0 and b > 0; the error
# names the offending pair. Row 1 keeps the legacy roles (single-variant tests,
# conditioning-weight default). More than one pair adds explicit beta weight
# sets to the group test (see SAIGE.Region).
parse_weights_beta <- function(weights.beta) {
  if (is.character(weights.beta)) {
    if (length(weights.beta) != 1 || is.na(weights.beta) || !nzchar(weights.beta)) {
      stop("weights.beta must be a single string of semicolon-separated 'a,b' pairs, e.g. '1,25' or '1,25;1,1'\n")
    }
    pairStrings <- strsplit(weights.beta, ";", fixed = TRUE)[[1]]
    if (length(pairStrings) == 0) {
      stop("weights.beta must contain at least one 'a,b' pair, e.g. '1,25'\n")
    }
    pairList <- lapply(pairStrings, function(p) {
      vals <- suppressWarnings(as.numeric(strsplit(p, ",", fixed = TRUE)[[1]]))
      if (length(vals) != 2 || anyNA(vals)) {
        stop("weights.beta: each pair must be exactly 2 numeric values 'a,b'; offending pair: '", p, "'\n")
      }
      vals
    })
    mat <- do.call(rbind, pairList)
  } else if (is.matrix(weights.beta)) {
    if (!is.numeric(weights.beta) || ncol(weights.beta) != 2 || nrow(weights.beta) < 1) {
      stop("weights.beta given as a matrix must be a numeric Kb x 2 matrix of beta(a, b) parameter pairs\n")
    }
    mat <- matrix(as.numeric(weights.beta), nrow = nrow(weights.beta), ncol = 2)
  } else if (is.numeric(weights.beta)) {
    if (length(weights.beta) != 2) {
      stop("weights.beta given as a numeric vector must have exactly 2 elements c(a, b), e.g. c(1, 25)\n")
    }
    mat <- matrix(as.numeric(weights.beta), nrow = 1, ncol = 2)
  } else {
    stop("weights.beta must be a pair-list string ('a1,b1;a2,b2;...'), a numeric vector c(a, b), or a Kb x 2 numeric matrix\n")
  }
  for (k in seq_len(nrow(mat))) {
    if (!all(is.finite(mat[k, ])) || any(mat[k, ] <= 0)) {
      stop(
        "weights.beta: beta distribution parameters must be finite and > 0; offending pair ",
        k, ": '", mat[k, 1], ",", mat[k, 2], "'\n"
      )
    }
  }
  dimnames(mat) <- NULL
  mat
}


#' Run single variant or gene- or region-based score tests with SPA based on the linear/logistic mixed model.
#'
#' @param bedFile character. Path to bed file (PLINK)
#' @param bimFile character. Path to bim file (PLINK)
#' @param famFile character. Path to fam file (PLINK)
#' @param pgenFile character. Path to pgen file (PGEN/plink2). Mutually exclusive with bedFile
#' @param pvarFile character. Path to pvar file (PGEN/plink2)
#' @param psamFile character. Path to psam file (PGEN/plink2)
#' @param AlleleOrder character. alt-first or ref-first for bgen or PLINK files. By default, alt-first
#' @param idstoIncludeFile character. Path to a file containing variant ids to be included from the dosage file. The file does not have a header and each line is for a marker ID. Variant ids are in the format chr:pos_ref/alt
#' @param rangestoIncludeFile character. Path to a file containing genome regions to be included from the dosage file. The file contains three columns for chromosome, start, and end respectively with no header. Note for vcf and sav files, only the first line in the file will be used.
#' @param chrom character. If LOCO is specified, chrom is required. chrom is also required for VCF/BCF/SAV input. Note: the string needs to exactly match the chromosome string in the vcf/sav file. For example, 1 does not match chr1.
#' @param is_imputed_data logical. Whether the dosages/genotypes imputed are imputed. If TRUE, the program will output the imputed info score. By default, FALSE.
#' @param minMAC numeric. Minimum minor allele count of markers to test. By default, 0.5. The higher threshold between minMAC and minMAF will be used
#' @param minMAF numeric. Minimum minor allele frequency of markers to test. By default 0. The higher threshold between minMAC and minMAF will be used
#' @param minInfo numeric. Minimum imputation info of markers to test. By default, 0.
#' @param maxMissing numeric. Maximum missing rate for markers to be tested. By default, 0.15
#' @param impute_method character. Imputation method for missing dosages. best_guess, mean or minor. best_guess: missing dosages imputed as best guessed genotyes round(2*allele frequency). mean: missing dosages are imputed as mean (2*allele frequency). minor: missing dosages are imputed as minor allele homozygotes. By default, minor
#' @param LOCO logical. Whether to apply the leave-one-chromosome-out option. If TRUE, --chrom is required. By default, TRUE
#' @param GMMATmodelFile character. Path to the input file containing the glmm model, which is output from previous step. Will be used by load()
#' @param varianceRatioFile character. Path to the input file containing the variance ratio, which is output from the previous step
#' @param SAIGEOutputFile character. Prefix of the output files containing assoc test results
#' @param markers_per_chunk character. Number of markers to be tested and output in each chunk in the single-variant assoc tests. By default, 10000
#' @param groups_per_chunk character. Number of groups/sets to be read in and tested in each chunk in the set-based assoc tests. By default, 100
#' @param is_output_moreDetails logical. Whether to output heterozygous and homozygous counts in cases and controls. By default, FALSE. If True, the columns homN_Allele2_cases, hetN_Allelelogical2_cases, homN_Allele2_ctrls, hetN_Allele2_ctrls will be output. By default, FALSE
#' @param is_overwrite_output logical. Whether to overwrite the output file if it exists. If FALSE, the program will continue the unfinished analysis instead of starting over from the beginining. By default, TRUE
#' @param maxMAF_in_groupTest. vector of numeric. Max MAF for markers tested in group test seperated by comma. e.g. c(0.0001,0.001,0.01). By default, c(0.01)
#' @param maxMAC_in_groupTest. vector of numeric. Max MAC for markers tested in group test seperated by comma. This vector will be combined with maxMAF_in_groupTest. e.g. c(1) to only test singletons. By default, c(0) and no Max MAC cutoffs are applied.
#' @param minGroupMAC_in_BurdenTest numeric. Only applied when only Burden tests are performed (r.corr=1). Minimum minor allele count in the Burden test for the psueodo marker. By default, 5
#' @param  annotation_in_groupTest. vector of character. annotations of markers to be tested in the set-based tests. using ; to combine multiple annotations in the same test. e.g. c("lof","missense;lof","missense;lof;synonymous")  will test lof variants only, missense+lof variants, and missense+lof+synonymous variants. By default:  c("lof","missense;lof","missense;lof;synonymous")
#' @param groupFile character. Path to the file containing the group information for gene-based tests. Each gene/set has 2 + K lines in the group file (K >= 0 weight lines). The first element of each line is the gene/set name. The second element indicates whether the line contains variant IDs (var), annotations (anno), or weights. A weight line's second element is 'weight' (weight set named "WEIGHT"), or 'weight:NAME' / 'weight,NAME' to name the weight set NAME. Weight lines are optional; with none, the default weights are generated based on beta(MAF, 1, 25) (use weights.beta to change the parameters). With more than one weight line per gene/set (each named, same names and order in every gene/set), each weight line defines one weight set and results are output per (annotation, max MAF, weight set) with a Weight column labeling the set (plus Pvalue_ACATV/Pvalue_ACATO columns for SKAT-O). The variant ids must match the genotype input's marker IDs (bim column 2 for PLINK). Elements are seperated by tab or space.
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file that was used in Step 1
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to sample IDs in the sparse GRM
#' @param relatednessCutoff float. The threshold for coefficient of relatedness to treat two samples as unrelated in the sparse GRM. By default, 0
#' @param MACCutoff_to_CollapseUltraRare numeric. MAC cutoff to collpase the ultra rare variants (<= MACCutoff_to_CollapseUltraRare) in the set-based association tests. By default, 10.
#' @param cateVarRatioMinMACVecExclude vector of float. Lower bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation. By default, c(10.5,20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param weights.beta parameters for the beta distribution(s) to weight genetic markers in gene-based tests: a numeric vector with two elements c(a, b), a Kb x 2 numeric matrix with one (a, b) pair per row, or a string "a1,b1;a2,b2;..." (semicolon-separated pairs, comma within a pair). By default, c(1,25). The first pair keeps the legacy roles (single-variant weighting and the conditioning-weight default). When more than one pair is given, each pair becomes one beta(MAF; a, b) weight set in the group test, named "Beta_a_b" (e.g. "Beta_1_25"), appended after any group-file weight sets; results are then output per weight set with a Weight column. With a single pair, group-file weights (if present) take precedence, exactly as before.
#' @param r.corr numeric. bewteen 0 and 1. parameters for gene-based tests. If r.corr = 1, only Burden tests will be performed. If r.corr = 0, SKAT-O tests will be performed and results for Burden tests and SKAT tests will be output too.  By default, 0.
#' @param markers_per_chunk_in_groupTest numeric. Number of markers in each chunk when calculating the variance covariance matrix in the set/group-based tests. By default, 100.
#' @param condition character. For conditional analysis. Variant ids are in the format chr:pos_ref/alt and seperated by by comma. e.g."chr3:101651171:C:T,chr3:101651186:G:A".
#' @param weights_for_condition. vector of numeric. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, e.g. c(1,2,3). If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use weights.beta to change the parameters for the Beta distribution.
#' @param SPAcutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation)
#' @param dosage_zerod_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 0.2
#' @param dosage_zerod_MAC_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 10
#' @param is_single_in_groupTest logical.  Whether to output single-variant assoc test results when perform group tests. Note, single-variant assoc test results will always be output when SKAT and SKAT-O tests are conducted with r.corr=0. This parameter should only be used when only Burden tests are condcuted with r.corr=1. By default, TRUE
#' @param is_no_weight_in_groupTest logical. Whether no weights (all-ones) are used in the group test. If FALSE, weights will be calcuated based on MAF from the Beta distribution with paraemters weights.beta or weights will be extracted from the group File if available. If TRUE and the run has no other weight sets (unweighted group file, single weights.beta pair), the all-ones weights simply replace the beta weights, exactly as before. If TRUE and the run does have other weight sets (weight lines in the group file, or more than one weights.beta pair), an all-ones weight set named "No_Weights" is appended as the last weight set (this replaces the previous hard error for weighted group files). By default, FALSE
#' @param is_output_weight_column logical. Whether to force the multi-weight output schema even when only a single weight set is in play: the Weight column (labeled with the set's name) right after max_MAF, and for SKAT-O the additive Pvalue_ACATV/Pvalue_ACATO columns and an always-emitted Cauchy row. Useful for a stable schema across runs with different weight configurations. With more than one weight set the schema is always on and this flag is redundant. By default, FALSE
#' @param is_output_markerList_in_groupTest logical. Whether to output the marker lists included in the set-based tests for each mask. By default, FALSE
#' @param is_Firth_beta logical. Whether to estimate effect sizes using approx Firth, only for binary traits. By default, FALSE
#' @param pCutoffforFirth numeric. p-value cutoff to use approx Firth to estiamte the effect sizes. Only for binary traits. The effect sizes of markers with p-value <= pCutoffforFirth will be estimated using approx Firth. By default, 0.01.
#' @param X_PARregion character. ranges of (pseudoautosomal) PAR region on chromosome X, which are seperated by comma and in the format start:end. By default: '60001-2699520,154931044-155260560' in the UCSC build hg19. For males, there are two X alleles in the PAR region, so PAR regions are treated the same as autosomes. In the NON-PAR regions (outside the specified PAR regions on chromosome X), for males, there is only one X allele. If is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in the NON-PAR regions on chromosome X will be multiplied by 2 (Not activated).
#' @param is_rewrite_XnonPAR_forMales logical. Whether to rewrite gentoypes or dosages of variants in the NON-PAR regions on chromosome X for males (multiply by 2). By default, FALSE. Note, only use is_rewrite_XnonPAR_forMales=TRUE when the specified VCF or Bgen file only has variants on chromosome X. When is_rewrite_XnonPAR_forMales=TRUE, the program does not check the chromosome value by assuming all variants are on chromosome X (Not activated)
#' @param sampleFile_male character. Path to the file containing one column for IDs of MALE samples in the bgen or vcf file with NO header. Order does not matter
#' @param nThreads integer. Number of threads for the single-variant marker scan. The C++ marker loop is threaded in place, so there is one process, one copy of the null model and one output file; the result is byte-identical to a one-thread run, which is the gate in tests/testthat/test-step2-threads-inprocess.R. Restarting an interrupted run works exactly as it does serially. Group tests always run on one thread, and a genotype reader that cannot be cloned (PLINK when the .bed could not be memory mapped) falls back to one thread with a warning. PGEN threads too, since S8. By default, 1
#' @return SAIGEOutputFile
#' @export
SPAGMMATtest <- function(bedFile = "",
                         bimFile = "",
                         famFile = "",
                         pgenFile = "",
                         pvarFile = "",
                         psamFile = "",
                         AlleleOrder = "alt-first", # new
                         idstoIncludeFile = "",
                         rangestoIncludeFile = "",
                         chrom = "", # for vcf file
                         max_missing = 0.15, # new
                         impute_method = "best_guess", # "mean", "minor", "best_guess"     #new
                         min_MAC = 0.5,
                         min_MAF = 0,
                         min_Info = 0,
                         is_imputed_data = FALSE, # new
                         GMMATmodelFile = "",
                         LOCO = TRUE,
                         varianceRatioFile = "",
                         cateVarRatioMinMACVecExclude = c(10.5, 20.5),
                         cateVarRatioMaxMACVecInclude = c(20.5),
                         SPAcutoff = 2,
                         SAIGEOutputFile = "",
                         markers_per_chunk = 10000,
                         groups_per_chunk = 100,
                         markers_per_chunk_in_groupTest = 100, # new
                         condition = "",
                         sparseGRMFile = "",
                         sparseGRMSampleIDFile = "",
                         relatednessCutoff = 0,
                         groupFile = "",
                         weights.beta = c(1, 25),
                         weights_for_condition = NULL,
                         r.corr = 0,
                         dosage_zerod_cutoff = 0.2,
                         dosage_zerod_MAC_cutoff = 10,
                         is_output_moreDetails = FALSE, # new
                         X_PARregion = "60001-2699520,154931044-155270560",
                         is_rewrite_XnonPAR_forMales = FALSE,
                         sampleFile_male = "",
                         MACCutoff_to_CollapseUltraRare = 10,
                         annotation_in_groupTest = c("lof", "missense;lof", "missense;lof;synonymous"), # new
                         maxMAF_in_groupTest = c(0.01),
                         maxMAC_in_groupTest = c(0),
                         minGroupMAC_in_BurdenTest = 5,
                         is_Firth_beta = FALSE,
                         pCutoffforFirth = 0.01,
                         is_overwrite_output = TRUE,
                         is_single_in_groupTest = TRUE,
                         is_no_weight_in_groupTest = FALSE,
                         is_output_markerList_in_groupTest = FALSE,
                         is_fastTest = FALSE,
                         pval_cutoff_for_fastTest = 0.05,
                         max_MAC_use_ER = 4,
                         subSampleFile = "",
                         is_noadjCov = FALSE,
                         is_output_weight_column = FALSE,
                         nThreads = 1) {
  cat("is_noadjCov ", is_noadjCov, "\n")
  cat("is_fastTest ", is_fastTest, "\n")

  # time_0 = proc.time()

  # cat("r.corr is ", r.corr, "\n")
  if (!(impute_method %in% c("best_guess", "mean", "minor"))) {
    stop("impute_method should be 'best_guess', 'mean' or 'minor'.")
  }

  nThreads <- as.integer(nThreads)
  if (is.na(nThreads) || nThreads < 1) {
    stop("nThreads should be a positive integer.")
  }
  if (nThreads > 1 && groupFile != "") {
    cat("WARNING: nThreads > 1 is only supported for single-variant tests; the group test will run on one thread\n")
    nThreads <- 1L
  }


  checkArgsListBool(
    is_imputed_data = is_imputed_data,
    LOCO = LOCO,
    is_output_moreDetails = is_output_moreDetails,
    is_overwrite_output = is_overwrite_output
  )
  # is_rewrite_XnonPAR_forMales = is_rewrite_XnonPAR_forMales)
  cat("dosage_zerod_cutoff ", dosage_zerod_cutoff, "\n")
  checkArgsListNumeric(
    start = 1,
    end = 250000000,
    max_missing = max_missing,
    min_MAC = min_MAC,
    min_MAF = min_MAF,
    min_Info = min_Info,
    SPAcutoff = SPAcutoff,
    dosage_zerod_cutoff = dosage_zerod_cutoff,
    dosage_zerod_MAC_cutoff = dosage_zerod_MAC_cutoff,
    markers_per_chunk = markers_per_chunk,
    groups_per_chunk = groups_per_chunk,
    minGroupMAC_in_BurdenTest = minGroupMAC_in_BurdenTest,
    max_MAC_use_ER = max_MAC_use_ER
  )

  # time_1 = proc.time()
  # if(file.exists(SAIGEOutputFile)) {print("ok -2 file exist")}


  ## check and create the output file
  # Check_OutputFile_Create(SAIGEOutputFile)
  OutputFile <- SAIGEOutputFile
  OutputFileIndex <- NULL
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  ## canonicalize weights.beta into the Kb x 2 matrix of beta(a, b) pairs
  ## (design D10). Row 1 keeps the legacy roles below; extra rows become
  ## explicit beta weight sets in the group test (SAIGE.Region).
  weights.beta <- parse_weights_beta(weights.beta)

  ## check the variance ratio file and extract the variance ratio vector
  setAssocTest_GlobalVarsInCPP(
    impute_method,
    max_missing,
    min_MAF,
    min_MAC,
    min_Info,
    dosage_zerod_cutoff,
    dosage_zerod_MAC_cutoff,
    as.numeric(weights.beta[1, ]),
    OutputFile,
    max_MAC_use_ER
  )

  # time_2 = proc.time()

  if (groupFile == "") {
    isGroupTest <- FALSE
    cat("single-variant association test will be performed\n")

    setMarker_GlobalVarsInCPP(
      is_output_moreDetails,
      markers_per_chunk
    )
  } else {
    isGroupTest <- TRUE
    Check_File_Exist(groupFile, "groupFile")
    cat("group-based test will be performed\n")
    # checkArgsList_for_Region(method_to_CollapseUltraRare,
    # order the max MAF from lowest to highest
    maxMAF_in_groupTest <- maxMAF_in_groupTest[order(maxMAF_in_groupTest)]
    maxMAC_in_groupTest <- maxMAC_in_groupTest[order(maxMAC_in_groupTest)]

    checkArgsList_for_Region(
      MACCutoff_to_CollapseUltraRare,
      # DosageCutoff_for_UltraRarePresence,
      maxMAF_in_groupTest = maxMAF_in_groupTest,
      maxMAC_in_groupTest = maxMAC_in_groupTest,
      markers_per_chunk_in_groupTest = markers_per_chunk_in_groupTest
    )


    # if(file.exists(SAIGEOutputFile)) {print("ok -1 file exist")}

    # IsOutputlogPforSingle = FALSE   #to check
    # OUT_Filename_Single<-sprintf("%s.single",SAIGEOutputFile)
    # Check_OutputFile_Create(OUT_Filename_Single)
    # if (sum(weights.beta.rare != weights.beta.common) > 0) {
    #  cat("WARNING:The option for weights.beta.common is not fully developed\n")
    #  cat("weights.beta.common is set to be equal to weights.beta.rare\n")
    #  weights.beta.common = weights.beta.rare
    # }

    # method_to_CollapseUltraRare,
    # DosageCutoff_for_UltraRarePresence,
    setRegion_GlobalVarsInCPP(
      maxMAF_in_groupTest,
      markers_per_chunk_in_groupTest,
      MACCutoff_to_CollapseUltraRare,
      minGroupMAC_in_BurdenTest
    )
    # cat("dosage_zerod_cutoff is ", dosage_zerod_cutoff, "\n")
    # cat("dosage_zerod_MAC_cutoff is ", dosage_zerod_MAC_cutoff, "\n")
  }


  # time_3 = proc.time()

  if (subSampleFile == "") {
    obj.model <- ReadModel(GMMATmodelFile, chrom, LOCO, is_Firth_beta) # readInGLMM.R
  } else {
    obj.model <- ReadModel_subsample(GMMATmodelFile, chrom, LOCO, is_Firth_beta, subSampleFile) # readInGLMM.R
  }


  # time_4 = proc.time()


  if (is_rewrite_XnonPAR_forMales) {
    cat("is_rewrite_XnonPAR_forMales is TRUE, so genotypes/dosages in the non-PAR regions of X chromosome for males will be multiplied by 2\n")
    if (!file.exists(sampleFile_male)) {
      stop(
        "ERROR! The sample file for male IDs ", sampleFile_male,
        " does not exist\n"
      )
    } else {
      sampleList_male <- data.frame(data.table:::fread(sampleFile_male,
        header = FALSE, stringsAsFactors = FALSE, colClasses = c("character"),
        data.table = FALSE
      ))
      colnames(sampleList_male) <- c("sampleID_male")
      cat(
        nrow(sampleList_male), " sample IDs are found in ",
        sampleFile_male, "\n"
      )
      indexInModel_male <- which(obj.model$sampleID %in% (sampleList_male$sampleID_male))
      cat(length(indexInModel_male), " males are found in the test\n")
      if (length(indexInModel_male) == 0) {
        is_rewrite_XnonPAR_forMales <- FALSE
        if (nrow(sampleList_male) > 0) {
          cat(
            "WARNING: no male IDs specified in the ",
            sampleFile_male, " are found sample IDs used to fit in the null model in Step 1\n"
          )
        }
      } else {
        cat("is_rewrite_XnonPAR_forMales=TRUE and minInfo and minMAF won't be applied to all X chromosome variants\n")
        minInfo <- 0
        minMAF <- 1 / (2 * length(obj.model$sampleID))

        setAssocTest_GlobalVarsInCPP_indexInModel_male(indexInModel_male - 1)
      }
    }
    X_PARregion_list <- unlist(strsplit(X_PARregion, split = ","))
    X_PARregion_mat <- NULL
    if (length(X_PARregion_list) > 0) {
      for (lxp in 1:length(X_PARregion_list)) {
        X_PARregion_list_sub <- as.numeric(unlist(strsplit(X_PARregion_list[lxp],
          split = "-"
        )))
        X_PARregion_mat <- rbind(X_PARregion_mat, X_PARregion_list_sub)
      }
      setAssocTest_GlobalVarsInCPP_X_PARregion_mat(X_PARregion_mat)
    } else {
      cat("PAR region on X chromosome is not specified\n")
    }
  }


  # time_5 = proc.time()


  if (obj.model$traitType == "binary") {
    if (max_MAC_use_ER > 0) {
      cat("P-values of genetic variants with MAC <= ", max_MAC_use_ER, " will be calculated via effecient resampling.\n")
      if (max_MAC_use_ER > 4) {
        cat("WARNING: Efficient resampling may not work well for MAC > 4!")
      }
    }
  } else {
    max_MAC_use_ER <- 0
  }

  if (!LOCO) {
    # 	LOCO = FALSE
    print("LOCO = FASLE and leave-one-chromosome-out is not applied")
  }

  sparseSigmaRList <- list()
  isSparseGRM <- TRUE
  if (sparseGRMFile != "") {
    # sparseSigmaRList = setSparseSigma(sparseSigmaFile)
    sparseSigmaRList <- setSparseSigma_new(sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff, obj.model$sampleID, obj.model$theta, obj.model$mu2, obj.model$traitType)
    isSparseGRM <- TRUE
  } else {
    # if(!is.null(obj.model$useSparseGRMforVarRatio)){
    # 	if(obj.model$useSparseGRMforVarRatio == TRUE){
    # 		stop("sparse GRM is not specified but it was used in Step 1.\n")
    # 	}
    # }
    sparseSigmaRList <- list(nSubj = 0, locations = matrix(0, nrow = 2, ncol = 2), values = rep(0, 2))
    isSparseGRM <- FALSE
  }

  cat("isSparseGRM is ", isSparseGRM, "\n")
  ratioVecList <- Get_Variance_Ratio(varianceRatioFile, cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude, isGroupTest, isSparseGRM) # readInGLMM.R

  # LEDGER #10 fix: single-variant analysis dereferences the null_noXadj
  # variance ratios (the no-covariate-adjustment fast pass); a file without
  # such rows used to crash in C++ mid-analysis (Mat index out of bounds).
  if (!isGroupTest & length(ratioVecList$ratioVec_null_noXadj) == 0) {
    stop(
      "varianceRatioFile ", varianceRatioFile, " contains no 'null_noXadj' rows, ",
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
  if (nThreads > 1 & length(ratioVecList$ratioVec_null) > 1 &
    varianceRatioBinsHaveGap(cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude)) {
    stop(
      "cateVarRatioMinMACVecExclude/cateVarRatioMaxMACVecInclude leave a gap in MAC ",
      "coverage, so some markers fall in no variance-ratio category at all. ",
      "Such a marker has no ratio to use and the test will stop when it is ",
      "reached, whatever --nThreads is. Give contiguous MAC categories.\n"
    )
  }


  # time_6 = proc.time()


  if (is_fastTest) {
    if (!file.exists(varianceRatioFile)) {
      is_fastTest <- FALSE
      cat("No variance ratio file is specified, so is_fastTest is not working.\n")
    }
    if (pval_cutoff_for_fastTest == 0) {
      is_fastTest <- FALSE
      cat("pval_cutoff_for_fastTest == 0, so is_fastTest is not working.\n")
    }
  }

  if (!isSparseGRM & !is_noadjCov) {
    is_fastTest <- FALSE
    cat("isSparseGRM is FALSE and is_noadjCov is FALSE,, so is_fastTest is not working.\n")
  }


  if (is_fastTest) {
    if (isSparseGRM) {
      cat("is_fastTest is TRUE.\n")
      if (ratioVecList$ratioVec_null[1] == -1) {
        stop("Variance ratios estimated without GRM are not found, so the fast tests can't be performed.\n Please set is_fastTest=FALSE or re-run the variance ratio estimation in Step 1 with --skipModelFitting=TRUE using the most recent version of the program.\n")
      } else {
        cat("The fast tests will be performed (when p-values >= ", pval_cutoff_for_fastTest, ").\n")
      }
    }

    # else{
    #     is_fastTest = FALSE
    #  cat("No sparse GRM is specified, so is_fastTest is not working.\n")
    # }
  }

  nsample <- length(obj.model$y)
  cateVarRatioMaxMACVecInclude <- c(cateVarRatioMaxMACVecInclude, nsample)

  # time_6 = proc.time()


  # in Geno.R
  objGeno <- setGenoInput(
    bedFile = bedFile,
    bimFile = bimFile,
    famFile = famFile,
    pgenFile = pgenFile,
    pvarFile = pvarFile,
    psamFile = psamFile,
    idstoIncludeFile = idstoIncludeFile,
    rangestoIncludeFile = rangestoIncludeFile,
    chrom = chrom,
    AlleleOrder = AlleleOrder,
    sampleInModel = obj.model$sampleID
  )
  # time_7 = proc.time()


  genoType <- objGeno$genoType
  if (condition != "") {
    isCondition <- TRUE
    if (is_fastTest) {
      is_fastTest <- FALSE
      cat("is_fastTest is not working for conditional analysis.\n")
    }
  } else {
    isCondition <- FALSE
  }

  condition_genoIndex <- c(-1)
  if (isCondition) {
    cat("Conducting conditional analysis. Please specify the conditioning markers in the order as they are store in the genotype/dosage file.\n")
  }
  # set up the SAIGE object based on the null model results
  setSAIGEobjInCPP(
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
    t_varRatio_sparse = as.vector(ratioVecList$ratioVec_sparse),
    t_varRatio_null = as.vector(ratioVecList$ratioVec_null),
    t_varRatio_null_noXadj = as.vector(ratioVecList$ratioVec_null_noXadj),
    t_cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
    t_cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
    t_SPA_Cutoff = SPAcutoff,
    t_tauvec = obj.model$theta,
    t_traitType = obj.model$traitType,
    t_y = obj.model$y,
    t_impute_method = impute_method,
    t_flagSparseGRM = isSparseGRM,
    t_isFastTest = is_fastTest,
    t_isnoadjCov = is_noadjCov,
    t_pval_cutoff_for_fastTest = pval_cutoff_for_fastTest,
    t_locationMat = as.matrix(sparseSigmaRList$locations),
    t_valueVec = sparseSigmaRList$values,
    t_dimNum = sparseSigmaRList$nSubj,
    t_isCondition = isCondition,
    t_condition_genoIndex = condition_genoIndex,
    t_is_Firth_beta = is_Firth_beta,
    t_pCutoffforFirth = pCutoffforFirth,
    t_offset = obj.model$offset,
    t_resout = as.integer(obj.model$obj_cc$res.out)
  )
  rm(sparseSigmaRList)
  gc()

  # time_8 = proc.time()
  # process condition
  if (isCondition) {
    n <- length(obj.model$y) # sample size

    ## re-order the conditioning markers
    ## condition_original = unlist(strsplit(condition, ","))
    condition_genoIndex <- extract_genoIndex_condition(condition, objGeno$markerInfo, genoType)
    if (!is.null(weights_for_condition)) {
      condition_weights <- weights_for_condition
      # print(condition_weights)
      # print(condition_genoIndex$cond_genoIndex)
      # condition_weights = as.numeric(unlist(strsplit(weights_for_condition, ",")))
      if (length(condition_weights) != length(condition_genoIndex$cond_genoIndex)) {
        stop("The length of the provided weights for conditioning markers is not equal to the number of conditioning markers\n")
      }
    } else {
      condition_weights <- rep(0, length(condition_genoIndex$cond_genoIndex))
    }

    condition_genoIndex_a <- as.character(format(condition_genoIndex$cond_genoIndex, scientific = FALSE))
    assign_conditionMarkers_factors(genoType, condition_genoIndex_a, n, condition_weights)
    if (obj.model$traitType == "binary" & isGroupTest) {
      outG2cond <- RegionSetUpConditional_binary_InCPP(condition_weights)


      outG2cond$pval_G2_cond <- unlist(lapply(outG2cond$pval_G2_cond, convert_str_to_log))

      G2condList <- get_newPhi_scaleFactor(q.sum = outG2cond$qsum_G2_cond, mu.a = obj.model$mu, g.sum = outG2cond$gsum_G2_cond, p.new = outG2cond$pval_G2_cond, Score = outG2cond$Score_G2_cond, Phi = outG2cond$VarMat_G2_cond, "SKAT-O")
      # print(G2condList)
      scaleFactorVec <- as.vector(G2condList$scaleFactor)
      # print(scaleFactorVec)
      assign_conditionMarkers_factors_binary_region(scaleFactorVec)
    }
  } else {
    condition_weights <- c(0)
  }

  traitType <- obj.model$traitType
  mu <- obj.model$mu
  rm(obj.model)
  gc()
  # print(gc(v=T))
  # if(file.exists(SAIGEOutputFile)) {print("ok 0 file exist")}


  # cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  # cat("Number of markers in each chunk:\t", numLinesOutput, "\n")
  # cat("Number of chunks for all markers:\t", nChunks, "\n")
  # }

  # time_9 = proc.time()


  if (!isGroupTest) {
    OutputFile <- SAIGEOutputFile

    # if(file.exists(SAIGEOutputFile)) {print("ok 2 file exist")}
    if (!is.null(objGeno$markerInfo$CHROM)) {
      setorderv(objGeno$markerInfo, col = c("CHROM", "POS"))
    }

    runMarkerScan(
      traitType,
      genoType,
      idstoIncludeFile,
      rangestoIncludeFile,
      objGeno$markerInfo$genoIndex,
      objGeno$markerInfo$CHROM,
      OutputFile,
      OutputFileIndex,
      markers_per_chunk,
      is_output_moreDetails,
      is_imputed_data,
      is_Firth_beta,
      LOCO,
      chrom,
      isCondition,
      is_overwrite_output,
      objGeno$anyInclude,
      nThreads
    )
  } else {
    maxMACbinind <- which(maxMAC_in_groupTest > 0)
    if (length(maxMACbinind) > 0) {
      maxMAC_in_groupTest_to_MAF <- (maxMAC_in_groupTest[maxMACbinind]) / (2 * length(mu))
      cat("maxMAC_in_groupTest: ", maxMAC_in_groupTest, " is specified, corresponding to max MAF ", maxMAC_in_groupTest_to_MAF, "\n")
      for (i in 1:length(maxMACbinind)) {
        checkArgNumeric(maxMAC_in_groupTest_to_MAF[i], deparse(substitute(maxMAC_in_groupTest_to_MAF[i])), 0, 0.5, FALSE, TRUE)
      }
      maxMAF_in_groupTest <- unique(c(maxMAF_in_groupTest, maxMAC_in_groupTest_to_MAF))
      maxMAF_in_groupTest <- maxMAF_in_groupTest[order(maxMAF_in_groupTest)]
      cat("max MAF cutoff ", maxMAF_in_groupTest, "will be applied\n")
    }
    # method_to_CollapseUltraRare,
    # DosageCutoff_for_UltraRarePresence,
    if (r.corr == 0 & is_Firth_beta) {
      print("WARNING: Note that the Firth correction has not been implemented for Burden test effect sizes when SKAT-O test is conducted. If the corrected Burden test effect sizes is needed, please use r.corr=1 to only conduct Burden test.")
      if (is_single_in_groupTest) {
        print("Firth correction will be used for effect sizes of single variant tests")
      }
    }

    SAIGE.Region(mu,
      OutputFile,
      MACCutoff_to_CollapseUltraRare,
      groupFile,
      annotation_in_groupTest,
      maxMAF_in_groupTest,
      markers_per_chunk_in_groupTest,
      genoType,
      objGeno$markerInfo,
      traitType,
      is_imputed_data,
      isCondition,
      condition_weights,
      groups_per_chunk,
      r.corr,
      is_overwrite_output,
      is_single_in_groupTest,
      is_no_weight_in_groupTest,
      is_output_markerList_in_groupTest,
      chrom,
      is_fastTest,
      pval_cutoff_for_fastTest,
      is_output_moreDetails,
      weights_beta_mat = weights.beta,
      is_output_weight_column = is_output_weight_column
    )
  }
}
