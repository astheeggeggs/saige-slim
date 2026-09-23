#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

options(stringsAsFactors = FALSE)

## load R libraries
library(SAIGE)
library(optparse)

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option(
    "--seed",
    type = "integer",
    default = 1L,
    help = paste0(
      "Random seed for the run; also decouples variance-ratio marker selection ",
      "from model-fitting convergence (default 1) [LEDGER #8]"
    )
  ),
  make_option(
    "--plinkFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to plink file for creating the genetic relationship matrix (GRM). ",
      "minMAFforGRM can be used to specify the minimum MAF and maxMissingRate can ",
      "be used to specify the maximum missing rates  of markers in the plink file ",
      "to be used for constructing GRM. Genetic markers are also randomly ",
      "selected from the plink file to estimate the variance ratios"
    )
  ),
  make_option(
    "--bedFile",
    type = "character",
    default = "",
    help = "Path to bed file. If plinkFile is specified, 'plinkFile'.bed will be used"
  ),
  make_option(
    "--bimFile",
    type = "character",
    default = "",
    help = "Path to bim file. If plinkFile is specified, 'plinkFile'.bim will be used"
  ),
  make_option(
    "--famFile",
    type = "character",
    default = "",
    help = "Path to fam file. If plinkFile is specified, 'plinkFile'.fam will be used"
  ),
  make_option(
    "--plink2File",
    type = "character",
    default = "",
    help = paste0(
      "Prefix for plink2 pgen/pvar/psam files, used in place of plinkFile. ",
      "Mutually exclusive with the PLINK input"
    )
  ),
  make_option(
    "--pgenFile",
    type = "character",
    default = "",
    help = "Path to pgen file. If plink2File is specified, 'plink2File'.pgen will be used"
  ),
  make_option(
    "--pvarFile",
    type = "character",
    default = "",
    help = "Path to pvar file. If plink2File is specified, 'plink2File'.pvar will be used"
  ),
  make_option(
    "--psamFile",
    type = "character",
    default = "",
    help = "Path to psam file. If plink2File is specified, 'plink2File'.psam will be used"
  ),
  make_option(
    "--phenoFile",
    type = "character",
    default = "",
    help = paste0(
      "Required. Path to the phenotype file. The file can be either tab or space ",
      "delimited. The phenotype file has a header and contains at least two ",
      "columns. One column is for phentoype and the other column is for sample ",
      "IDs. Additional columns can be included in the phenotype file for ",
      "covariates in the null model. Please specify the names of the covariates ",
      "using the argument covarColList and specify categorical covariates using ",
      "the argument qCovarColList. All categorical covariates must also be ",
      "included in covarColList."
    )
  ),
  make_option(
    "--phenoCol",
    type = "character",
    default = "",
    help = paste0(
      "Required. Column name for phenotype to be tested in the phenotype file, ",
      "e.g CAD. A COMMA-SEPARATED LIST fits several phenotypes in one process, ",
      "sharing the genotype ingest: outputs are named ",
      "<outputPrefix>_<phenoCol> and a manifest <outputPrefix>.multiTraits.txt ",
      "is written for step 2's ",
      "--GMMATmodel_varianceRatio_multiTraits_File. Each phenotype's outputs ",
      "are identical to a standalone run's"
    )
  ),
  make_option(
    "--traitType",
    type = "character",
    default = "binary",
    help = paste0(
      "Required. binary or quantitative [default=binary]. With a ",
      "comma-separated --phenoCol, give one value per phenotype or one for ",
      "all of them"
    )
  ),
  make_option(
    "--invNormalize",
    type = "character",
    default = "FALSE",
    help = paste0(
      "Optional. Only for quantitative. Whether to perform the inverse ",
      "normalization for the phenotype [default='FALSE']. With a ",
      "comma-separated --phenoCol, give one value per phenotype or one for ",
      "all of them. Not needed for scale: the null model is fitted on the ",
      "phenotype standardised by the sd of its covariate-adjusted values and ",
      "scaled back, so a trait in any unit fits the same way (LEDGER #65)"
    )
  ),
  make_option(
    "--covarColList",
    type = "character",
    default = "",
    help = "List of covariates (comma separated)"
  ),
  make_option(
    "--qCovarColList",
    type = "character",
    default = "",
    help = paste0(
      "List of categorical covariates (comma separated). All categorical ",
      "covariates must also be in covarColList"
    )
  ),
  make_option(
    "--sampleIDColinphenoFile",
    type = "character",
    default = "IID",
    help = "Required. Column name of sample IDs in the phenotype file, e.g. IID"
  ),
  make_option(
    "--tol",
    type = "double",
    default = 0.02,
    help = "Optional. Tolerance for fitting the null GLMM to converge [default=0.02]."
  ),
  make_option(
    "--tolCoef",
    type = "double",
    default = 1e-5,
    help = paste0("Optional. Tolerance of the inner coefficient pass at each variance-component ",
                  "iterate of a binary fit (LEDGER #125; --tolCoef=0.02 reproduces the old ",
                  "behaviour, where it was tied to --tol) [default=1e-5].")
  ),
  make_option(
    "--maxiter",
    type = "integer",
    default = 20L,
    help = paste0(
      "Optional. Maximum number of iterations used to fit the null GLMM ",
      "[default=20]."
    )
  ),
  make_option(
    "--tolPCG",
    type = "double",
    default = 1e-05,
    help = "Optional. Tolerance for PCG to converge [default=1e-5]."
  ),
  make_option(
    "--maxiterPCG",
    type = "integer",
    default = 500L,
    help = "Optional. Maximum number of iterations for PCG [default=500]."
  ),
  make_option(
    "--nThreads",
    type = "integer",
    default = 1L,
    help = "Optional. Number of threads (CPUs) to use [default=1]."
  ),
  make_option(
    "--SPAcutoff",
    type = "double",
    default = 2,
    help = paste0(
      "Optional. Cutoff for the deviation of score test statistics from mean in ",
      "the unit of sd to perform SPA [default=2]."
    )
  ),
  make_option(
    "--numRandomMarkerforVarianceRatio",
    type = "integer",
    default = 30L,
    help = paste0(
      "Optional. An integer greater than 0. Number of markers to be randomly ",
      "selected for estimating the variance ratio. The number will be ",
      "automatically added by 10 until the coefficient of variantion (CV) for the ",
      "variance ratio estimate is below ratioCVcutoff [default=30]."
    )
  ),
  make_option(
    "--skipModelFitting",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to skip model fitting and only to estimate the variance ",
      "ratio. If TRUE, the file outputPrefix.rda is required [default='FALSE']"
    )
  ),
  make_option(
    "--skipVarianceRatioEstimation",
    type = "logical",
    default = FALSE,
    help = "Optional. Whether to skip the variance ratio estimation. [default='FALSE']"
  ),
  make_option(
    "--memoryChunk",
    type = "double",
    default = 2,
    help = "Optional. Size (Gb) for each memory chunk [default=2]"
  ),
  make_option(
    "--tauInit",
    type = "character",
    default = "0,0",
    help = paste0(
      "Optional. Initial values for tau, in phenotype units for a quantitative ",
      "trait (the fit itself runs on the phenotype standardised by the sd of its ",
      "covariate-adjusted values and is scaled back; LEDGER #65). [default=0,0]"
    )
  ),
  make_option(
    "--LOCO",
    type = "logical",
    default = NULL,
    help = paste0(
      "Whether to apply the leave-one-chromosome-out (LOCO) approach when fitting ",
      "the null model using the full GRM. LOCO requires a full GRM: it is refused ",
      "with --useSparseGRMtoFitNULL=TRUE, because a sparse GRM has no per-chromosome ",
      "structure to leave out. If not supplied, defaults to FALSE for a sparse-GRM ",
      "fit and TRUE otherwise [default=TRUE, except sparse-GRM fits]."
    )
  ),
  make_option(
    "--isLowMemLOCO",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whehter to output the model file by chromosome when LOCO=TRUE. If TRUE, ",
      "the memory usage in Step 1 and Step 2 will be lower [default=FALSE]"
    )
  ),
  make_option(
    "--traceCVcutoff",
    type = "double",
    default = 0.0025,
    help = paste0(
      "Optional. Threshold for coefficient of variation (CV) for the trace ",
      "estimator. Number of runs for trace estimation will be increased until the ",
      "CV is below the threshold [default=0.0025]."
    )
  ),
  make_option(
    "--nrun",
    type = "double",
    default = 30,
    help = "Number of rums in trace estimation. [default=30]"
  ),
  make_option(
    "--ratioCVcutoff",
    type = "double",
    default = 0.001,
    help = paste0(
      "Optional. Threshold for coefficient of variation (CV) for estimating the ",
      "variance ratio. The number of randomly selected markers will be increased ",
      "until the CV is below the threshold [default=0.001]"
    )
  ),
  make_option(
    "--outputPrefix",
    type = "character",
    default = "",
    help = paste0(
      "Required. Path and prefix of the output files. Give a path on a mounted volume ",
      "in a container; ~ resolves to the container's own home [default='']"
    )
  ),
  make_option(
    "--outputPrefix_varRatio",
    type = "character",
    default = "",
    help = paste0(
      "Optional. Path and prefix of the output the variance ratio file. if not ",
      "specified, it will be the same as the outputPrefix"
    )
  ),
  make_option(
    "--IsOverwriteVarianceRatioFile",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to overwrite the variance ratio file if the file ",
      "exist.[default='FALSE']"
    )
  ),
  make_option(
    "--sparseGRMFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the pre-calculated sparse GRM file. If not specified and  ",
      "IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]"
    )
  ),
  make_option(
    "--sparseGRMSampleIDFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the sample ID file for the pre-calculated sparse GRM. No header is ",
      "included. The order of sample IDs is corresponding to sample IDs in the ",
      "sparse GRM [default=NULL]"
    )
  ),
  make_option(
    "--isCateVarianceRatio",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Required. Whether to estimate variance ratio based on different MAC ",
      "categories. If yes, variance ratio will be estiamted for multiple MAC ",
      "categories corresponding to cateVarRatioMinMACVecExclude and ",
      "cateVarRatioMaxMACVecInclude. Combines freely with --LOCO (nothing forces ",
      "LOCO=FALSE) [default=FALSE]"
    )
  ),
  make_option(
    "--relatednessCutoff",
    type = "double",
    default = 0,
    help = paste0(
      "Optional. Threshold (minimum relatedness coefficient) to treat two samples ",
      "as unrelated when the sparse GRM is used [default=0]"
    )
  ),
  make_option(
    "--cateVarRatioMinMACVecExclude",
    type = "character",
    default = "10,20.5",
    help = paste0(
      "Optional. vector of float. Lower bound for MAC categories. The length ",
      "equals to the number of MAC categories for variance ratio estimation. ",
      "[default='10,20.5']"
    )
  ),
  make_option(
    "--cateVarRatioMaxMACVecInclude",
    type = "character",
    default = "20.5",
    help = paste0(
      "Optional. vector of float. Higher bound for MAC categories. The length ",
      "equals to the number of MAC categories for variance ratio estimation minus ",
      "1. [default='20.5']"
    )
  ),
  make_option(
    "--isCovariateTransform",
    type = "logical",
    default = TRUE,
    help = "Optional. Whether use qr transformation on covariates [default='TRUE']."
  ),
  make_option(
    "--isDiagofKinSetAsOne",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to set the diagnal elements in GRM to be 1 ",
      "[default='FALSE']."
    )
  ),
  make_option(
    "--useSparseGRMtoFitNULL",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to use sparse GRM to fit the null model ",
      "[default='FALSE']."
    )
  ),
  make_option(
    "--usePCGwithSparseGRM",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to use PCG for inverse of the sparse GRM if ",
      "--useSparseGRMtoFitNULL=TRUE. [default='FALSE']."
    )
  ),
  make_option(
    "--useSparseGRMforVarRatio",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to use sparse GRM to estimate the variance Ratios. If ",
      "TRUE, the variance ratios will be estimated using the full GRM (numerator) ",
      "and the sparse GRM (denominator). By default, FALSE"
    )
  ),
  make_option(
    "--minMAFforGRM",
    type = "double",
    default = 0.01,
    help = "Optional. Minimum MAF of markers used for GRM"
  ),
  make_option(
    "--maxMissingRateforGRM",
    type = "double",
    default = 0.15,
    help = "Optional. Maximum missing rate of markers used for GRM"
  ),
  make_option(
    "--minCovariateCount",
    type = "double",
    default = -1,
    help = paste0(
      "Optional. Binary covariates with a count less than minCovariateCount will ",
      "be excluded from the model to avoid convergence issues [default=-1] (no ",
      "covariates will be excluded)."
    )
  ),
  make_option(
    "--minCaseCount",
    type = "integer",
    default = 100,
    help = paste0(
      "Optional. Binary traits only. Refuse to fit if the final sample set has ",
      "fewer than this many cases. Below a few hundred cases the step-2 SKAT-O ",
      "null tail is liberal (measured: 2.7x nominal at p=1e-4 with 100 cases, ",
      "1.5x at 1e-3 with 500) [default=100]. 0 disables the gate; a fit with ",
      "fewer than 500 cases is warned about and still run."
    )
  ),
  make_option(
    "--minSampleCount",
    type = "integer",
    default = 0,
    help = paste0(
      "Optional. Refuse to fit if the final sample set has fewer than this many ",
      "samples with a non-missing phenotype. Both trait types [default=0, no gate]."
    )
  ),
  make_option(
    "--minCovariateCellCount",
    type = "integer",
    default = 10,
    help = paste0(
      "Optional. Refuse to fit if any level of a categorical covariate has fewer ",
      "than this many cases or controls (samples, for a quantitative trait). A ",
      "level with zero cases is complete separation [default=10]. 0 disables the ",
      "gate. Distinct from --minCovariateCount, which EXCLUDES a low-count ",
      "covariate instead of refusing; the exclusion runs first."
    )
  ),
  make_option(
    "--separationCheck",
    type = "character",
    default = "error",
    help = paste0(
      "Optional. What to do when the covariate GLM shows separation (fails to ",
      "converge, or a fitted probability reaches the 0/1 boundary): 'error', ",
      "'warn' or 'none' [default='error']. Catches separation on a continuous ",
      "covariate, which a cell count cannot. Binary traits only."
    )
  ),
  make_option(
    "--convergenceCheck",
    type = "character",
    default = "error",
    help = paste0(
      "Optional. What to do when the null fit did not converge -- the AI-REML ",
      "loop ran out without meeting its criterion, so the variance components ",
      "and the variance ratio derived from them never settled: 'error', ",
      "'warn' or 'none' [default='error']. Both trait types."
    )
  ),
  make_option(
    "--separationMaxFraction",
    type = "double",
    default = 0.02,
    help = paste0(
      "Optional. The fraction of samples whose fitted probability may sit at ",
      "the 0/1 boundary before --separationCheck treats the model as ",
      "separated. A pinned sample contributes ~0 to the score variance, so it ",
      "carries no information: 0.5%% pinned is a 0.5%% loss of effective ",
      "sample size, 100%% is total collapse. Above zero but at or below this ",
      "is a warning and the run continues [default=0.02]."
    )
  ),
  make_option(
    "--includeNonautoMarkersforVarRatio",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to allow for non-autosomal markers for variance ratio. ",
      "[default, 'FALSE']"
    )
  ),
  make_option(
    "--FemaleOnly",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to run Step 1 for females only [default=FALSE]. if TRUE, ",
      "--sexCol and --FemaleCode need to be specified"
    )
  ),
  make_option(
    "--MaleOnly",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Optional. Whether to run Step 1 for males only [default=FALSE]. if TRUE, ",
      "--sexCol and --MaleCode need to be specified"
    )
  ),
  make_option(
    "--sexCol",
    type = "character",
    default = "",
    help = "Optional. Column name for sex in the phenotype file, e.g Sex"
  ),
  make_option(
    "--FemaleCode",
    type = "character",
    default = "1",
    help = paste0(
      "Optional. Values in the column for sex in the phenotype file are used for ",
      "females [default, '1']"
    )
  ),
  make_option(
    "--MaleCode",
    type = "character",
    default = "0",
    help = paste0(
      "Optional. Values in the column for sex in the phenotype file are used for ",
      "males [default, '0']"
    )
  ),
  make_option(
    "--isCovariateOffset",
    type = "logical",
    default = TRUE,
    help = paste0(
      "Optional. TRUE: covariate effects are estimated once by a fixed-effects GLM ",
      "and entered as an offset; FALSE: re-estimated inside the mixed model [default=TRUE]"
    )
  ),
  make_option(
    "--SampleIDIncludeFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to the file that contains one column for IDs of samples who will be ",
      "include for null model fitting."
    )
  )
)

## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

covars <- strsplit(opt$covarColList, ",")[[1]]
qcovars <- strsplit(opt$qCovarColList, ",")[[1]]

convertToNumeric <- function(x, stringOutput) {
  y <- tryCatch(expr = as.numeric(x), warning = function(w) NULL)
  if (is.null(y)) {
    stop(stringOutput, " is not numeric\n")
  }
  cat(stringOutput, " is ", y, "\n")
  y
}

tauInit <- convertToNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")
cateVarRatioMinMACVecExclude <- convertToNumeric(
  strsplit(opt$cateVarRatioMinMACVecExclude, ",")[[1]],
  "cateVarRatioMinMACVecExclude"
)
cateVarRatioMaxMACVecInclude <- convertToNumeric(
  strsplit(opt$cateVarRatioMaxMACVecInclude, ",")[[1]],
  "cateVarRatioMaxMACVecInclude"
)

set.seed(opt$seed)

# --phenoCol, --traitType and --invNormalize may each be a comma-separated list.
# One phenotype is the legacy path exactly: fitNULLGLMM_multi calls fitNULLGLMM
# once with the unsuffixed output names and writes no manifest.
phenoCols <- trimws(unlist(strsplit(opt$phenoCol, ",", fixed = TRUE)))
traitTypes <- trimws(unlist(strsplit(opt$traitType, ",", fixed = TRUE)))
invNormalizes <- as.logical(trimws(unlist(strsplit(opt$invNormalize, ",", fixed = TRUE))))
if (anyNA(invNormalizes)) {
  stop("--invNormalize must be TRUE or FALSE (comma-separated for several phenotypes)\n")
}

# The upstream script dispatched on packageVersion("SAIGE") across four
# near-duplicate fitNULLGLMM calls; the in-tree package version is fixed
# (DESCRIPTION: 1.5.2), so only the >= 1.4.5 call is live and the dead
# branches were removed. The argument list below is that branch's, verbatim.
fitNULLGLMM_multi(
  plinkFile = opt$plinkFile,
  seed = opt$seed,
  bedFile = opt$bedFile,
  bimFile = opt$bimFile,
  famFile = opt$famFile,
  plink2File = opt$plink2File,
  pgenFile = opt$pgenFile,
  pvarFile = opt$pvarFile,
  psamFile = opt$psamFile,
  useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL,
  usePCGwithSparseGRM = opt$usePCGwithSparseGRM,
  sparseGRMFile = opt$sparseGRMFile,
  sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
  phenoFile = opt$phenoFile,
  phenoCol = phenoCols,
  sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
  traitType = traitTypes,
  outputPrefix = opt$outputPrefix,
  isCovariateOffset = opt$isCovariateOffset,
  nThreads = opt$nThreads,
  useSparseGRMforVarRatio = opt$useSparseGRMforVarRatio,
  invNormalize = invNormalizes,
  covarColList = covars,
  qCovarCol = qcovars,
  tol = opt$tol,
  tolCoef = opt$tolCoef,
  maxiter = opt$maxiter,
  tolPCG = opt$tolPCG,
  maxiterPCG = opt$maxiterPCG,
  SPAcutoff = opt$SPAcutoff,
  numMarkersForVarRatio = opt$numRandomMarkerforVarianceRatio,
  skipModelFitting = opt$skipModelFitting,
  skipVarianceRatioEstimation = opt$skipVarianceRatioEstimation,
  memoryChunk = opt$memoryChunk,
  tauInit = tauInit,
  LOCO = opt$LOCO,
  isLowMemLOCO = opt$isLowMemLOCO,
  traceCVcutoff = opt$traceCVcutoff,
  nrun = opt$nrun,
  ratioCVcutoff = opt$ratioCVcutoff,
  outputPrefix_varRatio = opt$outputPrefix_varRatio,
  IsOverwriteVarianceRatioFile = opt$IsOverwriteVarianceRatioFile,
  relatednessCutoff = opt$relatednessCutoff,
  isCateVarianceRatio = opt$isCateVarianceRatio,
  cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
  cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
  isCovariateTransform = opt$isCovariateTransform,
  isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
  minMAFforGRM = opt$minMAFforGRM,
  maxMissingRateforGRM = opt$maxMissingRateforGRM,
  minCovariateCount = opt$minCovariateCount,
  minCaseCount = opt[["minCaseCount"]],
  minSampleCount = opt[["minSampleCount"]],
  minCovariateCellCount = opt[["minCovariateCellCount"]],
  separationCheck = opt[["separationCheck"]],
  convergenceCheck = opt[["convergenceCheck"]],
  separationMaxFraction = opt[["separationMaxFraction"]],
  includeNonautoMarkersforVarRatio = opt$includeNonautoMarkersforVarRatio,
  sexCol = opt$sexCol,
  FemaleCode = opt$FemaleCode,
  FemaleOnly = opt$FemaleOnly,
  MaleCode = opt$MaleCode,
  MaleOnly = opt$MaleOnly,
  SampleIDIncludeFile = opt$SampleIDIncludeFile
)
