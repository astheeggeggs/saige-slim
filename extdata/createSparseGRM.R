#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript

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
    help = "path to plink file for creating the genetic relationship matrix (GRM)"
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
    "--nThreads",
    type = "integer",
    default = 1L,
    help = "Number of threads (CPUs) to use [default=1]"
  ),
  make_option(
    "--memoryChunk",
    type = "numeric",
    default = 2,
    help = "The size (Gb) for each memory chunk. By default, 2"
  ),
  make_option(
    "--outputPrefix",
    type = "character",
    default = "",
    help = paste0(
      "path and prefix to the output files. Give a path on a mounted volume ",
      "in a container; ~ resolves to the container's own home [default='']"
    )
  ),
  make_option(
    "--numRandomMarkerforSparseKin",
    type = "integer",
    default = 2000L,
    help = paste0(
      "number of randomly selected markers to be used to identify related ",
      "samples for sparse GRM [default=2000]"
    )
  ),
  make_option(
    "--relatednessCutoff",
    type = "numeric",
    default = 0.125,
    help = paste0(
      "The threshold to treat two samples as unrelated if IsSparseKin is TRUE ",
      "[default=0.125]"
    )
  ),
  make_option(
    "--isDiagofKinSetAsOne",
    type = "logical",
    default = FALSE,
    help = "Whether to set the diagnal elements in GRM to be 1 [default='FALSE']."
  ),
  make_option(
    "--minMAFforGRM",
    type = "numeric",
    default = 0.01,
    help = "minimum MAF of markers used for GRM"
  ),
  make_option(
    "--maxMissingRateforGRM",
    type = "numeric",
    default = 0.15,
    help = "Optional. Maximum missing rate of markers used for GRM"
  ),
  make_option(
    "--screenMethod",
    type = "character",
    default = "pairwise",
    help = paste0(
      "How the relatedness screen over all sample pairs is computed: pairwise ",
      "(one dot product per pair) or blocked (one single-threaded sgemm per ",
      "block pair of samples, parallel over block pairs; pairs near the cutoff ",
      "are re-scored with the pairwise arithmetic, so the output is the same) ",
      "[default=pairwise]"
    )
  ),
  make_option(
    "--screenBlockSize",
    type = "integer",
    default = 2048L,
    help = paste0(
      "Samples per block for --screenMethod=blocked; each thread holds a ",
      "screenBlockSize^2 float matrix (16 MB at 2048) [default=2048]"
    )
  )
)

## list of options
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

set.seed(1)

createSparseGRM(
  plinkFile = opt$plinkFile,
  bedFile = opt$bedFile,
  bimFile = opt$bimFile,
  famFile = opt$famFile,
  plink2File = opt$plink2File,
  pgenFile = opt$pgenFile,
  pvarFile = opt$pvarFile,
  psamFile = opt$psamFile,
  outputPrefix = opt$outputPrefix,
  numRandomMarkerforSparseKin = opt$numRandomMarkerforSparseKin,
  relatednessCutoff = opt$relatednessCutoff,
  memoryChunk = opt$memoryChunk,
  isDiagofKinSetAsOne = opt$isDiagofKinSetAsOne,
  nThreads = opt$nThreads,
  minMAFforGRM = opt$minMAFforGRM,
  maxMissingRateforGRM = opt$maxMissingRateforGRM,
  screenMethod = opt$screenMethod,
  screenBlockSize = opt$screenBlockSize
)
