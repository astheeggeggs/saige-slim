#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

# Estimate the effective sample size (Nglmm) of a null GLMM.
#
# This is a diagnostic companion to step 1, NOT a step of the pipeline: it fits
# the same null model as step1_fitNULLGLMM.R (via getNeff() in
# R/SAIGE_extractNeff.R) but stops after the variance-component fit and reports
# the effective sample size implied by it.
#
#   binary:       Nglmm = 4 * Pn * (1 - Pn) * t(1) Rinv 1   (Pn = case fraction)
#   quantitative: Nglmm = t(1) Rinv 1
#
# It writes NO output files -- Nglmm is printed to stdout on a line reading
# "Nglmm  <value>", so redirect the run to a log and read the value off the end
# of it. `outputPrefix` / `outputPrefix_varRatio` are accepted for signature
# compatibility with step 1 and are not used to write anything here.
#
# Typical use is against a sparse GRM built by createSparseGRM.R
# (--useSparseGRMtoFitNULL=TRUE with --sparseGRMFile / --sparseGRMSampleIDFile);
# see extdata/cmd_extractNeff.sh for a worked binary + quantitative example.
#
# The option help strings below are kept verbatim from step1_fitNULLGLMM.R so
# that the shared arguments document identically in both scripts.

options(stringsAsFactors = FALSE)

## load R libraries
library(SAIGE)
library(optparse)

print(sessionInfo())

## set list of cmd line arguments
option_list <- list(
  make_option(
    "--plinkFile",
    type = "character",
    default = "",
    help = paste0(
      "Path to plink file for creating the genetic relationship matrix (GRM). ",
      "minMAFforGRM can be used to specify the minimum MAF of markers in he ",
      "plink file to be used for constructing GRM. Genetic markers are also ",
      "randomly selected from the plink file to estimate the variance ratios"
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
      "Path to the phenotype file. The file can be either tab or space ",
      "delimited. The phenotype file has a header and contains at least two ",
      "columns. One column is for phentoype and the other column is for sample ",
      "IDs. Additional columns can be included in the phenotype file for ",
      "covariates in the null GLMM. Please note that covariates to be used in ",
      "the NULL GLMM need to specified using the argument covarColList"
    )
  ),
  make_option(
    "--phenoCol",
    type = "character",
    default = "",
    help = "Coloumn name for phenotype to be tested in the phenotype file, e.g CAD"
  ),
  make_option(
    "--traitType",
    type = "character",
    default = "binary",
    help = "binary/quantitative [default=binary]"
  ),
  make_option(
    "--invNormalize",
    type = "logical",
    default = FALSE,
    help = paste0(
      "if quantitative, whether asking SAIGE to perform the inverse ",
      "normalization for the phenotype [default='FALSE']"
    )
  ),
  make_option(
    "--covarColList",
    type = "character",
    default = "",
    help = "List of covariates (comma separated)"
  ),
  make_option(
    "--sampleIDColinphenoFile",
    type = "character",
    default = "IID",
    help = "Column name of sample IDs in the phenotype file, e.g. IID"
  ),
  make_option(
    "--tol",
    type = "numeric",
    default = 0.02,
    help = "Tolerance for fitting the null GLMM to converge [default=0.02]."
  ),
  make_option(
    "--maxiter",
    type = "integer",
    default = 20,
    help = "Maximum number of iterations used to fit the null GLMM [default=20]."
  ),
  make_option(
    "--tolPCG",
    type = "numeric",
    default = 1e-5,
    help = "Tolerance for PCG to converge [default=1e-5]."
  ),
  make_option(
    "--maxiterPCG",
    type = "integer",
    default = 500,
    help = "Maximum number of iterations for PCG [default=500]."
  ),
  make_option(
    "--nThreads",
    type = "integer",
    default = 1,
    help = "Number of threads (CPUs) to use [default=1]."
  ),
  make_option(
    "--memoryChunk",
    type = "numeric",
    default = 2,
    help = "Size (Gb) for each memory chunk [default=2]"
  ),
  make_option(
    "--tauInit",
    type = "character",
    default = "0,0",
    help = "Initial values for tau. [default=0,0]"
  ),
  make_option(
    "--outputPrefix",
    type = "character",
    default = "~/",
    help = "Path and prefix of the output files [default='~/']"
  ),
  make_option(
    "--outputPrefix_varRatio",
    type = "character",
    default = NULL,
    help = paste0(
      "Path and prefix of the output the variance ratio file [default=NULL]. ",
      "if NULL, it will be the same as the outputPrefix"
    )
  ),
  make_option(
    "--IsSparseKin",
    type = "logical",
    default = FALSE,
    help = "Whether to use sparse kinship for association test [default='FALSE']"
  ),
  make_option(
    "--sparseGRMFile",
    type = "character",
    default = NULL,
    help = paste0(
      "Path to the pre-calculated sparse GRM file. If not specified and  ",
      "IsSparseKin=TRUE, sparse GRM will be computed [default=NULL]"
    )
  ),
  make_option(
    "--sparseGRMSampleIDFile",
    type = "character",
    default = NULL,
    help = paste0(
      "Path to the sample ID file for the pre-calculated sparse GRM. No header ",
      "is included. The order of sample IDs is corresponding to sample IDs in ",
      "the sparse GRM [default=NULL]"
    )
  ),
  make_option(
    "--numRandomMarkerforSparseKin",
    type = "integer",
    default = 2000,
    help = paste0(
      "Number of randomly selected markers to be used to identify related ",
      "samples for sparse GRM [default=2000]"
    )
  ),
  make_option(
    "--relatednessCutoff",
    type = "numeric",
    default = 0.125,
    help = paste0(
      "Threshold to treat two samples as unrelated if IsSparseKin is TRUE ",
      "[default=0.125]"
    )
  ),
  make_option(
    "--isCovariateTransform",
    type = "logical",
    default = TRUE,
    help = "Whether use qr transformation on non-genetic covariates [default='TRUE']."
  ),
  make_option(
    "--useSparseSigmaConditionerforPCG",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to sparse GRM to speed up the PCG. Current this option is ",
      "deactivated. [default='FALSE']."
    )
  ),
  make_option(
    "--useSparseSigmaforInitTau",
    type = "logical",
    default = FALSE,
    help = "Whether to use sparse Sigma to estiamte initial tau [default='FALSE']."
  ),
  make_option(
    "--useSparseGRMtoFitNULL",
    type = "logical",
    default = FALSE,
    help = "Whether to use sparse GRM to fit the null model [default='FALSE']."
  ),
  make_option(
    "--minMAFforGRM",
    type = "numeric",
    default = 0.01,
    help = "Minimum MAF of markers used for GRM"
  ),
  make_option(
    "--minCovariateCount",
    type = "numeric",
    default = -1,
    help = paste0(
      "If binary covariates have a count less than this, they will be excluded ",
      "from the model to avoid convergence issues [default=-1] (no covariates ",
      "will be excluded)."
    )
  ),
  make_option(
    "--FemaleOnly",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to run Step 1 for females only [default=FALSE]. if TRUE, ",
      "--sexCol and --FemaleCode need to be specified"
    )
  ),
  make_option(
    "--MaleOnly",
    type = "logical",
    default = FALSE,
    help = paste0(
      "Whether to run Step 1 for males only [default=FALSE]. if TRUE, --sexCol ",
      "and --MaleCode need to be specified"
    )
  ),
  make_option(
    "--sexCol",
    type = "character",
    default = "",
    help = "Coloumn name for sex in the phenotype file, e.g Sex"
  ),
  make_option(
    "--FemaleCode",
    type = "character",
    default = "1",
    help = paste0(
      "Values in the column for sex in the phenotype file are used for females ",
      "[default, '1']"
    )
  ),
  make_option(
    "--MaleCode",
    type = "character",
    default = "0",
    help = paste0(
      "Values in the column for sex in the phenotype file are used for males ",
      "[default, '0']"
    )
  ),
  make_option(
    "--noEstFixedEff",
    type = "logical",
    default = FALSE,
    help = "Whether to estimate fixed effect coeffciets. [default, 'FALSE']"
  )
)

## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

covars <- strsplit(opt$covarColList, ",")[[1]]

convertToNumeric <- function(x, stringOutput) {
  y <- tryCatch(expr = as.numeric(x), warning = function(w) NULL)
  if (is.null(y)) {
    stop(stringOutput, " is not numeric\n")
  }
  cat(stringOutput, " is ", y, "\n")
  y
}

tauInit <- convertToNumeric(strsplit(opt$tauInit, ",")[[1]], "tauInit")

# Unlike step1_fitNULLGLMM.R (which takes --seed), the seed is fixed here: this
# script reports a single summary statistic and is not a golden-generating step.
set.seed(1)

getNeff(
  plinkFile = opt$plinkFile,
  bedFile = opt$bedFile,
  bimFile = opt$bimFile,
  famFile = opt$famFile,
  plink2File = opt$plink2File,
  pgenFile = opt$pgenFile,
  pvarFile = opt$pvarFile,
  psamFile = opt$psamFile,
  phenoFile = opt$phenoFile,
  phenoCol = opt$phenoCol,
  traitType = opt$traitType,
  invNormalize = opt$invNormalize,
  covarColList = covars,
  qCovarCol = NULL,
  sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
  tol = opt$tol,
  maxiter = opt$maxiter,
  tolPCG = opt$tolPCG,
  maxiterPCG = opt$maxiterPCG,
  nThreads = opt$nThreads,
  memoryChunk = opt$memoryChunk,
  tauInit = tauInit,
  outputPrefix = opt$outputPrefix,
  outputPrefix_varRatio = opt$outputPrefix_varRatio,
  IsSparseKin = opt$IsSparseKin,
  sparseGRMFile = opt$sparseGRMFile,
  sparseGRMSampleIDFile = opt$sparseGRMSampleIDFile,
  numRandomMarkerforSparseKin = opt$numRandomMarkerforSparseKin,
  relatednessCutoff = opt$relatednessCutoff,
  isCovariateTransform = opt$isCovariateTransform,
  useSparseSigmaConditionerforPCG = opt$useSparseSigmaConditionerforPCG,
  useSparseSigmaforInitTau = opt$useSparseSigmaforInitTau,
  minMAFforGRM = opt$minMAFforGRM,
  minCovariateCount = opt$minCovariateCount,
  sexCol = opt$sexCol,
  FemaleCode = opt$FemaleCode,
  FemaleOnly = opt$FemaleOnly,
  MaleCode = opt$MaleCode,
  MaleOnly = opt$MaleOnly,
  # LEDGER #43 FIXED: forwarded like every other option (it was a hardcoded FALSE).
  noEstFixedEff = opt$noEstFixedEff,
  useSparseGRMtoFitNULL = opt$useSparseGRMtoFitNULL
)
