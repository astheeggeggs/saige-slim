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

# One "<start>-<end>" token of X_PARregion as c(start, end). LEDGER #38 FIXED:
# a token that is not of that form used to be recycled by rbind into a row (a
# lone "5" became c(5, 5)) and handed to C++ as a PAR boundary with no warning.
parse_X_PARregion_token <- function(token) {
  bounds <- suppressWarnings(as.numeric(unlist(strsplit(token, split = "-"))))
  if (length(bounds) != 2 || anyNA(bounds) || bounds[1] > bounds[2]) {
    stop("X_PARregion token '", token, "' is not of the form start-end ",
         "with start <= end (the option is a comma-separated list of such ranges)\n")
  }
  bounds
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
#' @param cateVarRatioMinMACVecExclude vector of float. Exclusive lower bound of MAC for MAC categories (a marker with MAC equal to the bound is in the category BELOW it). The length equals to the number of MAC categories for variance ratio estimation. By default, c(10,20.5), the same as step 1 and as both CLI drivers. If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). If groupFile="", only one variance ratio corresponding to MAC >= 20 is used
#' @param weights.beta parameters for the beta distribution(s) to weight genetic markers in gene-based tests: a numeric vector with two elements c(a, b), a Kb x 2 numeric matrix with one (a, b) pair per row, or a string "a1,b1;a2,b2;..." (semicolon-separated pairs, comma within a pair). By default, c(1,25). The first pair keeps the legacy roles (single-variant weighting and the conditioning-weight default). When more than one pair is given, each pair becomes one beta(MAF; a, b) weight set in the group test, named "Beta_a_b" (e.g. "Beta_1_25"), appended after any group-file weight sets; results are then output per weight set with a Weight column. With a single pair, group-file weights (if present) take precedence, exactly as before.
#' @param r.corr numeric. bewteen 0 and 1. parameters for gene-based tests. If r.corr = 1, only Burden tests will be performed. If r.corr = 0, SKAT-O tests will be performed and results for Burden tests and SKAT tests will be output too.  By default, 0. The per-region Cauchy row is the ACAT of the cells' SKAT-O Pvalue, except that when any cell's Pvalue is exactly 1 (the SKAT package saturates its SKAT-O integral to 1 when the smallest per-rho p exceeds ~0.78) the row is min(1, n_cells * min(Pvalue)), a Bonferroni under the Cauchy label; see LEDGER #86.
#' @param markers_per_chunk_in_groupTest numeric. Number of markers in each chunk when calculating the variance covariance matrix in the set/group-based tests. By default, 100.
#' @param condition character. For conditional analysis. Variant ids are in the format chr:pos_ref/alt and seperated by by comma. e.g."chr3:101651171:C:T,chr3:101651186:G:A".
#' @param weights_for_condition. vector of numeric. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, e.g. c(1,2,3). If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use weights.beta to change the parameters for the Beta distribution.
#' @param SPAcutoff by default = 2 (SPA test would be used when p value < 0.05 under the normal approximation). The output column Is.SPA is true whenever the reported SE was back-derived from the reported p-value rather than from 1/sqrt(var): after the saddlepoint approximation and also after the exact resampling test on MAC <= max_MAC_use_ER markers, which has no saddlepoint (LEDGER #80).
#' @param dosage_zerod_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 0.2
#' @param dosage_zerod_MAC_cutoff numeric. If is_imputed_data = TRUE, For variants with MAC <= dosage_zerod_MAC_cutoff, dosages <= dosageZerodCutoff with be set to 0. By derault, 10
#' @param is_single_in_groupTest logical.  Whether to output single-variant assoc test results when perform group tests. Note, single-variant assoc test results will always be output when SKAT and SKAT-O tests are conducted with r.corr=0. This parameter should only be used when only Burden tests are condcuted with r.corr=1. By default, FALSE (the same as the command line since LEDGER #18; it was TRUE before 2026-09-05)
#' @param is_no_weight_in_groupTest logical. Whether no weights (all-ones) are used in the group test. If FALSE, weights will be calcuated based on MAF from the Beta distribution with paraemters weights.beta or weights will be extracted from the group File if available. If TRUE and the run has no other weight sets (unweighted group file, single weights.beta pair), the all-ones weights simply replace the beta weights, exactly as before. If TRUE and the run does have other weight sets (weight lines in the group file, or more than one weights.beta pair), an all-ones weight set named "No_Weights" is appended as the last weight set (this replaces the previous hard error for weighted group files). By default, FALSE
#' @param is_output_weight_column logical. Whether to force the multi-weight output schema even when only a single weight set is in play: the Weight column (labeled with the set's name) right after max_MAF, and for SKAT-O the additive Pvalue_ACATV/Pvalue_ACATO columns and an always-emitted Cauchy row. Useful for a stable schema across runs with different weight configurations. With more than one weight set the schema is always on and this flag is redundant. By default, FALSE
#' @param is_output_markerList_in_groupTest logical. Whether to output the marker lists included in the set-based tests for each mask. By default, FALSE
#' @param is_Firth_beta logical. Whether to estimate effect sizes using approx Firth, only for binary traits. By default, FALSE
#' @param pCutoffforFirth numeric. p-value cutoff to use approx Firth to estiamte the effect sizes. Only for binary traits. The effect sizes of markers with p-value <= pCutoffforFirth will be estimated using approx Firth. By default, 0.01.
#' @param firth_solver character. "closed" (default: 2x2 closed-form Fisher matrix and
#'   hat values, no QR; ~2x faster per Firth fit) or "qr" (the QR-based hat values of
#'   SAIGE 1.x). Same estimator; the two agree to <= 8e-13 relative, which does not reach
#'   the printed digits. Only for binary traits with is_Firth_beta = TRUE.
#' @param X_PARregion character. ranges of (pseudoautosomal) PAR region on chromosome X, which are seperated by comma and in the format start:end. By default: '60001-2699520,154931044-155260560' in the UCSC build hg19. For males, there are two X alleles in the PAR region, so PAR regions are treated the same as autosomes. In the NON-PAR regions (outside the specified PAR regions on chromosome X), for males, there is only one X allele. If is_rewrite_XnonPAR_forMales=TRUE, genotypes/dosages of all variants in the NON-PAR regions on chromosome X will be multiplied by 2 (Not activated).
#' @param is_rewrite_XnonPAR_forMales logical. Whether to rewrite gentoypes or dosages of variants in the NON-PAR regions on chromosome X for males (multiply by 2). By default, FALSE. Only variants whose chromosome is X (X, x, 23, chrX, chr23) are rewritten; other chromosomes in the same file are left alone (LEDGER #69).
#' @param sampleFile_male character. Path to the file containing one column for IDs of MALE samples in the bgen or vcf file with NO header. Order does not matter
#' @param nThreads integer. Number of threads for the single-variant marker scan. The C++ marker loop is threaded in place, so there is one process, one copy of the null model and one output file; the result is byte-identical to a one-thread run, which is the gate in tests/testthat/test-step2-threads-inprocess.R. Restarting an interrupted run works exactly as it does serially. For a GROUP test the marker loop stays on one thread, but `nThreads` is not wasted there: it drives the SKAT-O post-processing, which Rprof puts at 45% of a group run at N = 50,000 (measured 1.29x at 4 workers, byte-identical -- see harness/bench/NOTES.md). Unix only; on other platforms a group run stays serial. A genotype reader that cannot be cloned (PLINK when the .bed could not be memory mapped) falls back to one thread with a warning. PGEN threads too, since S8. By default, 1
#' @param GMMATmodel_varianceRatio_multiTraits_File character. Path to a manifest of null models to test in ONE marker scan, sharing the genotype read. Three whitespace-separated columns, no header: a name, the step-1 .rda, and that model's variance ratio file. Mutually exclusive with GMMATmodelFile/varianceRatioFile. Output goes to <SAIGEOutputFile>_<name> per model, each with its own .index -- except for a single-row manifest, which writes the unsuffixed <SAIGEOutputFile> and is byte- and name-identical to the equivalent GMMATmodelFile run. Models whose step-1 sample IDs are identical share one genotype reader and one decode per marker; models with different sample sets each get their own. Single-variant tests only: a groupFile with more than one model is an error. By default, ""
#' @param is_noadjCov logical. Whether to SKIP the covariate adjustment of the tested genotype in the single-variant score test. TRUE (the default, matching the command line) uses the unadjusted, allele-frequency-centred genotype with the null_noXadj variance ratio; FALSE regresses the covariates out of the genotype first and uses the null variance ratio. The two are different score tests and give different p-values (LEDGER #82: the API used to default to FALSE while the command line defaulted to TRUE). Ignored for group tests and for the sparse-GRM marker test, which always adjust. Under condition the conditioning terms are always built from the ADJUSTED genotype, so at TRUE the _c columns mix an unadjusted marginal score with adjusted conditioning terms (a marker in the conditioning set keeps a non-zero var_c and a finite p.value_c; a correlated marker's conditional variance is slightly overstated); use FALSE for exact conditional columns (LEDGER #82, settled). By default, TRUE
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
                         cateVarRatioMinMACVecExclude = c(10, 20.5), # LEDGER #71: was c(10.5, 20.5); step 1 and both CLIs use 10
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
                         firth_solver = "closed",
                         is_overwrite_output = TRUE,
                         is_single_in_groupTest = FALSE,
                         is_no_weight_in_groupTest = FALSE,
                         is_output_markerList_in_groupTest = FALSE,
                         is_fastTest = FALSE,
                         pval_cutoff_for_fastTest = 0.05,
                         max_MAC_use_ER = 4,
                         subSampleFile = "",
                         is_noadjCov = TRUE,
                         is_output_weight_column = FALSE,
                         nThreads = 1,
                         GMMATmodel_varianceRatio_multiTraits_File = "") {
  cat("is_noadjCov ", is_noadjCov, "\n")
  cat("is_fastTest ", is_fastTest, "\n")

  # time_0 = proc.time()

  # cat("r.corr is ", r.corr, "\n")
  if (!(impute_method %in% c("best_guess", "mean", "minor"))) {
    stop("impute_method should be 'best_guess', 'mean' or 'minor'.")
  }

  # Which null models is this run testing? One (the legacy --GMMATmodelFile
  # path) or a manifest of them. Resolved and validated HERE, before any state
  # has been pushed into C++, so a rejected run leaves nothing behind.
  modelSpecs <- resolveModelSpecs(
    GMMATmodelFile, varianceRatioFile,
    GMMATmodel_varianceRatio_multiTraits_File
  )
  nModels <- length(modelSpecs)
  if (nModels > 1) {
    # Three deliberate refusals, each because the feature keeps per-RUN state
    # that is model-order-dependent (see docs/state/multimodel-step2.md):
    if (groupFile != "") {
      stop(
        "A model manifest with ", nModels, " models is supported for ",
        "single-variant tests only; --groupFile names a region test. Run the ",
        "region tests one model at a time.\n"
      )
    }
    if (is_rewrite_XnonPAR_forMales) {
      stop(
        "is_rewrite_XnonPAR_forMales keeps ONE set of male sample indices in ",
        "model order, so it cannot be combined with a model manifest.\n"
      )
    }
    if (subSampleFile != "") {
      stop(
        "subSampleFile applies one subsample to one model, so it cannot be ",
        "combined with a model manifest.\n"
      )
    }
  }

  nThreads <- as.integer(nThreads)
  if (is.na(nThreads) || nThreads < 1) {
    stop("nThreads should be a positive integer.")
  }
  # Group tests keep the C++ region loop single-threaded, but --nThreads is no
  # longer thrown away there: it drives the SKAT-O post-processing, which Rprof
  # puts at 45% of a group run (harness/bench/NOTES.md). Carried separately so
  # nothing downstream of `nThreads` can be surprised by a value it never saw
  # before.
  nThreadsSKATO <- 1L
  if (nThreads > 1 && groupFile != "") {
    nThreadsSKATO <- nThreads
    cat("nThreads > 1 does not thread the group test's marker loop; it will be used for",
        "the SKAT-O post-processing instead\n")
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
  # One output file per model, and one .index each. The suffix is "" for a
  # single model (and for a one-row manifest), so this is the same single path
  # it always was in every existing configuration. The region path below keeps
  # using the scalar OutputFile.
  OutputFiles <- paste0(SAIGEOutputFile, vapply(modelSpecs, function(x) x$suffix, character(1)))
  OutputFileIndexes <- paste0(OutputFiles, ".index")
  OutputFileIndex <- OutputFileIndexes[1]

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

  # setAssocTest_GlobalVarsInCPP installed exactly one output prefix (the run
  # prefix). Only a genuine multi-model run needs more, and only it calls this,
  # so a single-model run's C++ call sequence is unchanged.
  if (nModels > 1) {
    setSingleOutputPrefixesInCPP(OutputFiles)
  }

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

  # Every model's null-model preparation, in manifest order. For nModels == 1
  # this is the same sequence of calls the inline block made, in the same order
  # -- only the male-X block below has moved after it, and that block reads
  # nothing the preparation writes (its two setters are C++ globals the marker
  # loop reads, and its minMAF/minInfo assignments are dead locals).
  prepared <- lapply(seq_len(nModels), function(k) {
    prepareNullModel(
      modelSpecs[[k]], k, chrom, LOCO, is_Firth_beta, subSampleFile,
      sparseGRMFile, sparseGRMSampleIDFile, relatednessCutoff,
      cateVarRatioMinMACVecExclude, cateVarRatioMaxMACVecInclude,
      isGroupTest, is_fastTest, pval_cutoff_for_fastTest,
      is_noadjCov, max_MAC_use_ER, nThreads
    )
  })
  obj.model <- prepared[[1]]$obj.model
  # The region path and the group-test setup still read these as scalars, and
  # they only ever run with one model.
  isSparseGRM <- prepared[[1]]$isSparseGRM
  is_fastTest <- prepared[[1]]$is_fastTest

  # Which models can share a genotype reader? Identical sampleID vectors mean
  # one reader, one decode and one QC pass per marker for all of them.
  grouping <- groupModelsBySampleSet(lapply(prepared, function(p) p$obj.model$sampleID))
  groupOf <- grouping$groupOf
  nGroups <- length(grouping$representatives)
  cat(
    "Single-variant tests will run", nModels, "null model(s) in",
    nGroups, "sample-set group(s)\n"
  )

  if (obj.model$traitType != "binary") {
    max_MAC_use_ER <- 0
  }

  if (!LOCO) {
    # 	LOCO = FALSE
    print("LOCO = FASLE and leave-one-chromosome-out is not applied")
  }

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
      X_PARregion_mat <- do.call(rbind, lapply(X_PARregion_list, parse_X_PARregion_token))
      setAssocTest_GlobalVarsInCPP_X_PARregion_mat(X_PARregion_mat)
    } else {
      cat("PAR region on X chromosome is not specified\n")
    }
  }


  # time_5 = proc.time()



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
    sampleInModel = grouping$representatives[[1]]
  )

  # Group 1's reader came from setGenoInput; the rest are appended here, in
  # group order, so g_readers[g] in C++ is group g. A single-model run has one
  # group and never enters this loop.
  if (nGroups > 1) {
    for (g in 2:nGroups) {
      addGenoReaderForGroup(objGeno, grouping$representatives[[g]])
    }
  }
  # time_7 = proc.time()


  genoType <- objGeno$genoType
  if (condition != "") {
    isCondition <- TRUE
    # Every model, not just the first: the flag is baked into each SAIGEClass at
    # push time, so disabling it on one and not the others would give the same
    # marker two different test paths in one scan.
    if (any(vapply(prepared, function(p) p$is_fastTest, logical(1)))) {
      cat("is_fastTest is not working for conditional analysis.\n")
    }
    for (k in seq_along(prepared)) {
      prepared[[k]]$is_fastTest <- FALSE
    }
    is_fastTest <- FALSE
  } else {
    isCondition <- FALSE
  }

  condition_genoIndex <- c(-1)
  if (isCondition) {
    cat("Conducting conditional analysis. Please specify the conditioning markers in the order as they are store in the genotype/dosage file.\n")
  }
  # set up the SAIGE objects based on the null model results, one per model
  pushNullModelsToCPP(
    prepared, groupOf, SPAcutoff, impute_method, is_noadjCov,
    pval_cutoff_for_fastTest, isCondition, condition_genoIndex,
    is_Firth_beta, pCutoffforFirth, firth_solver, cateVarRatioMinMACVecExclude
  )
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
    # Per model: the conditioning markers are read through THAT model's
    # sample-set group and the factors are written onto THAT model. For one
    # model this is the single call it always was, with index 0.
    for (k in seq_len(nModels)) {
      assign_conditionMarkers_factors(
        genoType, condition_genoIndex_a,
        length(prepared[[k]]$obj.model$y), condition_weights, as.integer(k - 1L)
      )
    }
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

  # Model 1's, for the region path and the group-test setup below; the marker
  # path takes each model's trait type from the model itself in C++.
  traitType <- obj.model$traitType
  mu <- obj.model$mu
  rm(obj.model)
  rm(prepared)
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
      genoType,
      idstoIncludeFile,
      rangestoIncludeFile,
      objGeno$markerInfo$genoIndex,
      objGeno$markerInfo$CHROM,
      OutputFiles,
      OutputFileIndexes,
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
      # LEDGER #61. The vector is only extended HERE, because converting a MAC
      # cutoff into a MAF one needs `mu`, which is not loaded until :551 -- but
      # setRegion_GlobalVarsInCPP was already called above with the UN-extended
      # vector, and it is what fixes g_cfg.maxMAFLimit, the per-marker QC gate in
      # mainRegionInCPP. So any bin the extension adds above the largest
      # --maxMAF_in_groupTest had every marker in its new range discarded at
      # ingest, and its rows came out a byte-for-byte DUPLICATE of the previous
      # bin, at exit 0, after the log said the cutoff "will be applied".
      # Re-issued rather than moved: this runs only when the extension actually
      # happened, so a run without --maxMAC_in_groupTest makes exactly the calls
      # it always made.
      setRegion_GlobalVarsInCPP(
        maxMAF_in_groupTest,
        markers_per_chunk_in_groupTest,
        MACCutoff_to_CollapseUltraRare,
        minGroupMAC_in_BurdenTest
      )
    }
    # method_to_CollapseUltraRare,
    # DosageCutoff_for_UltraRarePresence,
    # LEDGER #51. fastTest is a two-stage design -- screen with the cheap
    # variance ratio, re-run the interesting units accurately -- and the REGION
    # path only ever had the second stage. All three of its variance-ratio sites
    # (src/Main.cpp:1435, :1826, :2164) choose from the marker's MAC and never
    # consult m_isFastTest, so the "screen" runs the IDENTICAL arithmetic as the
    # confirmation, and an escalated region recomputes numbers it already had.
    # Measured at N = 50,000 over 210 genes: +16.3% for BYTE-IDENTICAL output.
    # Inherited verbatim from upstream (reported as saigegit/SAIGE#214), so this
    # warns rather than implementing the missing half -- doing that would
    # preserve no golden (a real screen uses a different variance ratio AND a
    # different score test, so genes near the cutoff could fail to escalate) and
    # the whole prize is bounded at 8%.
    if (is_fastTest) {
      cat("WARNING: --is_fastTest has no computational benefit for group tests.",
          "The region path does not use the fast variance ratio, so the screening",
          "pass costs a full extra pass over every escalated region for the same",
          "answer (+16.3% measured at N = 50,000). The results are unaffected;",
          "only pval_cutoff_for_fastTest changes which regions pay it.
")
    }
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
      is_output_weight_column = is_output_weight_column,
      nThreadsSKATO = nThreadsSKATO
    )
  }
}
