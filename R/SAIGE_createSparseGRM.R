#' Construct a sparse GRM for a given data set
#'
#' @param bedFile character. Path to plink file (bed) to be used for calculating the sparse GRM
#' @param bimFile character. Path to plink file (bim) to be used for calculating the sparse GRM
#' @param famFile character. Path to plink file (fam) to be used for calculating the sparse GRM
#' @param plinkFile character. Prefix for the plink bed/bim/fam files. Mutually
#'   exclusive with the PGEN input
#' @param plink2File character. Prefix for the plink2 pgen/pvar/psam files.
#'   Mutually exclusive with the PLINK input
#' @param pgenFile character. Path to plink2 file (pgen) to be used for
#'   calculating the sparse GRM
#' @param pvarFile character. Path to plink2 file (pvar) to be used for
#'   calculating the sparse GRM
#' @param psamFile character. Path to plink2 file (psam) to be used for
#'   calculating the sparse GRM
#' @param outputPrefix character. Path to the output files with prefix
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples for sparse GRM. By default, 1000
#' @param relatednessCutoff float. The threshold to treat two samples as unrelated if IsSparseKin is TRUE. By default, 0.125
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param isDiagofKinSetAsOne  logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param nThreads integer. Number of threads to be used. By default, 1
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @return a file ended with sampleIDs.txt that contains sample IDs for the sparse GRM and a file ended with .sparseGRM.mtx that contains the sparse GRM
#' @export
createSparseGRM <- function(bedFile = "",
                            bimFile = "",
                            famFile = "",
                            plinkFile = "",
                            plink2File = "",
                            pgenFile = "",
                            pvarFile = "",
                            psamFile = "",
                            outputPrefix = "",
                            numRandomMarkerforSparseKin = 1000,
                            relatednessCutoff = 0.125,
                            memoryChunk = 2,
                            isDiagofKinSetAsOne = FALSE,
                            nThreads = 1,
                            minMAFforGRM = 0.01,
                            maxMissingRateforGRM = 0.15,
                            isSetGeno = TRUE,
                            isWritetoFiles = TRUE) {
  # LEDGER #7 fix: always pin the TBB pool — nThreads=1 previously left
  # RcppParallel on all hardware threads (nondeterministic float reductions).
  RcppParallel:::setThreadOptions(numThreads = nThreads)
  cat(nThreads, " threads are set to be used ", "\n")
  # cat("numRandomMarkerforSparseKin is ", numRandomMarkerforSparseKin, "\n")
  cat("sparse GRM will be created\n")
  setminMAFforGRM(minMAFforGRM)

  if (minMAFforGRM > 0) {
    cat(
      "Markers in the Plink file with MAF < ", minMAFforGRM,
      " will be removed before constructing GRM\n"
    )
  }
  setmaxMissingRateforGRM(maxMissingRateforGRM)
  if (maxMissingRateforGRM > 0) {
    cat("Markers in the Plink file with missing rate > ", maxMissingRateforGRM, " will be removed before constructing GRM\n")
  }

  genoIn <- checkGenoInputStep1(
    plinkFile = plinkFile, bedFile = bedFile, bimFile = bimFile, famFile = famFile,
    plink2File = plink2File, pgenFile = pgenFile, pvarFile = pvarFile, psamFile = psamFile
  )
  genoType <- genoIn$genoType
  bedFile <- genoIn$bedFile
  bimFile <- genoIn$bimFile
  famFile <- genoIn$famFile
  pgenFile <- genoIn$pgenFile
  pvarFile <- genoIn$pvarFile
  psamFile <- genoIn$psamFile

  # Sample IDs in genotype-file row order -- fam column 2, or the psam IID
  # column. That order is the index space setgeno's subSampleInGeno lives in.
  sparseGRMSampleID <- readGenoSampleIDs(genoType, famFile = famFile, psamFile = psamFile)
  nSample <- length(sparseGRMSampleID)
  sparseGRMSampleIDFile <- paste0(outputPrefix, "_relatednessCutoff_", relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")

  if (isWritetoFiles) {
    cat("write sample IDs for the sparse GRM to ", sparseGRMSampleIDFile, "\n")
    write.table(sparseGRMSampleID, sparseGRMSampleIDFile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  }

  genoSampleIndex <- seq(1, nSample)

  cat("isDiagofKinSetAsOne ", isDiagofKinSetAsOne, "\n")

  indicatorGenoSamplesWithPheno <- rep(TRUE, nSample)

  if (isSetGeno) {
    setgenoFromInput(
      genoIn, genoSampleIndex, indicatorGenoSamplesWithPheno, memoryChunk,
      isDiagofKinSetAsOne
    )
  }
  freqVec <- getAlleleFreqVec()
  if (minMAFforGRM > 0) {
    MAFindex <- which(freqVec >= minMAFforGRM & freqVec <= 1 - minMAFforGRM)
    cat(length(MAFindex), " markers have MAF >= ", minMAFforGRM, "\n")
  } else {
    MAFindex <- which(freqVec > 0 & freqVec < 1)
    cat(length(MAFindex), " markers have MAF > ", minMAFforGRM, "\n")
  }

  cat(numRandomMarkerforSparseKin, " genetic markers are randomly selected to decide which samples are related\n")
  if (length(MAFindex) < numRandomMarkerforSparseKin) {
    if (minMAFforGRM > 0) {
      stop("ERROR! not enough genetic markers with MAF >= ", minMAFforGRM, " to detect which samples are related\n", "Try include at least ", numRandomMarkerforSparseKin, " genetic markers with MAF >= ", minMAFforGRM, " in the plink file\n")
    } else {
      stop("ERROR! not enough genetic markers with MAF > ", minMAFforGRM, " to detect which samples are related\n", "Try include at least ", numRandomMarkerforSparseKin, " genetic markers with MAF > ", minMAFforGRM, " in the plink file\n")
    }
  }

  markerIndexforSparseM <- sample(MAFindex, size = numRandomMarkerforSparseKin, replace = FALSE)

  cat("Start detecting related samples for the sparse GRM\n")
  ta <- proc.time()
  setSubMarkerIndex(markerIndexforSparseM - 1)
  tb <- proc.time()
  cat("tb-ta\n")
  print(tb - ta)


  cat("Start creating sparse GRM\n")
  ta <- proc.time()
  sparseMList <- createSparseKinParallel(nblocks = nThreads, ncore = nThreads, relatednessCutoff)
  tb <- proc.time()
  cat("tb-ta\n")
  print(tb - ta)


  cat("length(sparseMList$iIndex): ", length(sparseMList$iIndex), "\n")
  print(sparseMList$iIndex[1:102])
  cat("length(sparseMList$jIndex): ", length(sparseMList$jIndex), "\n")
  print(sparseMList$jIndex[1:102])
  cat("length(sparseMList$kinValue): ", length(sparseMList$kinValue), "\n")
  print(sparseMList$kinValue[1:102])
  sparseGRM <- Matrix:::sparseMatrix(i = as.vector(sparseMList$iIndex), j = as.vector(sparseMList$jIndex), x = as.vector(sparseMList$kinValue), symmetric = TRUE)
  cat("nrow(sparseGRM): ", nrow(sparseGRM), "\n")
  cat("ncol(sparseGRM): ", ncol(sparseGRM), "\n")
  cat("ncol(sparseGRM): ", sum(sparseGRM != 0), "\n")

  tc <- proc.time()
  cat("tc-tb\n")
  print(tc - tb)


  if (isWritetoFiles) {
    sparseGRMFile <- paste0(outputPrefix, "_relatednessCutoff_", relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx")
    cat("write sparse GRM to ", sparseGRMFile, "\n")
    Matrix:::writeMM(sparseGRM, sparseGRMFile)
  }
  list(sparseGRMSampleID = sparseGRMSampleID, sparseGRM = sparseGRM)
}
