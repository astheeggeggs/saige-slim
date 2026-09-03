SAIGE.Marker <- function(traitType,
                         genoType,
                         idstoIncludeFile,
                         rangestoIncludeFile,
                         genoIndex,
                         CHROM,
                         OutputFile,
                         OutputFileIndex = NULL,
                         nMarkersEachChunk,
                         isMoreOutput,
                         isImputation,
                         isFirth,
                         LOCO,
                         chrom,
                         isCondition,
                         isOverWriteOutput,
                         isAnyInclude,
                         # In-process threads for the marker scan; this is what
                         # the CLI's --nThreads becomes (S7 retired the forked
                         # path). At 1 the C++ side runs a literally serial loop
                         # that touches no threading machinery at all, so the
                         # default configuration cannot be affected by any of it.
                         nThreadsMarker = 1) {
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # genoType = objGeno$genoType

  outCheck <- checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific = FALSE), isOverWriteOutput) # this function is in 'Util.R'
  outIndex <- outCheck$indexChunk
  isappend <- FALSE
  if (outIndex != 1) {
    cat("Restart the analysis from chunk:\t", outIndex, "\n")
    isappend <- TRUE
    # LEDGER #47. The interrupted run may have flushed rows belonging to chunks
    # the index never recorded; appending would write them twice. Drop them
    # before the output is reopened (truncateOutputToIndex is in 'Util.R').
    #
    # The target is the size at the end of chunk outIndex - 1, NOT outIndex,
    # because the loop below RESUMES AT outIndex and re-runs it. That off-by-one
    # is pre-existing and contradicts checkOutputFile's own message ("we restart
    # the analysis from the indexChunk + 1 chunk"); it costs one chunk of
    # recomputation and, with this truncation in place, produces exactly the
    # file a complete run would. It was found by the #47 regression test, which
    # failed with a whole duplicated chunk when the target was outIndex.
    truncateOutputToIndex(OutputFile,
                          outCheck$chunkBytes[as.character(outIndex - 1)])
  }

  isOpenOutFile_single <- openOutfile_single(traitType, isImputation, isappend, isMoreOutput)

  if (!isOpenOutFile_single) {
    stop("Output file ", OutputFile, " can't be opened\n")
  }

  ## set up an object for genotype
  if (genoType %in% c("plink", "pgen")) {
    # markerInfo = objGeno$markerInfo
    if (LOCO) {
      genoIndex <- genoIndex[which(CHROM == chrom)]
      CHROM <- CHROM[which(CHROM == chrom)]
      # markerInfo = markerInfo[which(markerInfo$CHROM == chrom),]
    }
    # CHROM = markerInfo$CHROM
    # genoIndex = markerInfo$genoIndex
    ## only for one chrom
    # all markers were split into multiple chunks,
    # print(markerInfo[1:10,])
    # print(genoIndex[1:10])

    genoIndexList <- splitMarker(genoIndex, nMarkersEachChunk, CHROM)
    nChunks <- length(genoIndexList)

    if (nChunks == 0) {
      stop("No markers on chrom ", chrom, " are found\n")
    }


    cat("Number of all markers to test:\t", length(genoIndex), "\n")
    cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
    cat("Number of chunks for all markers:\t", nChunks, "\n")
    if (outIndex > nChunks) {
      cat("The analysis has been finished! Please delete ", OutputFileIndex, " if the analysis needs to be run again or set --is_overwrite_output=TRUE\n")
      is_marker_test <- FALSE
    } else {
      is_marker_test <- TRUE
      i <- outIndex
    }
  }
  # chrom = "InitialChunk"
  # set_flagSparseGRM_cur_SAIGE_org()
  while (is_marker_test) {
    # for(i in outIndex:nChunks)
    # {
    # time_left = system.time({


    if (genoType %in% c("plink", "pgen")) {
      tempList <- genoIndexList[[i]]
      genoIndex <- as.character(format(tempList$genoIndex, scientific = FALSE))
      tempChrom <- tempList$chrom
    }


    # print("tempList")
    # print(tempList)
    # print(tempList$genoIndex)
    # })
    # print("time_left")
    # print(time_left)
    # print("genoIndex here")
    # print(genoIndex)
    # set up objects that do not change for different variants
    # if(tempChrom != chrom){
    #  setMarker("SAIGE", objNull, control, chrom, Group, ifOutGroup)
    #  chrom = tempChrom
    # }
    # ptm <- proc.time()
    # print(ptm)
    # print("gc()")
    # print(gc())
    cat(paste0("(", Sys.time(), ") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom, " ---- \n"))
    # main function to calculate summary statistics for markers in one chunk

    # resMarker = as.data.frame(mainMarkerInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation))
    # resMarker = resMarker[which(!is.na(resMarker$BETA)), ]

    mainMarkerInCPP(genoType, traitType, genoIndex, isMoreOutput, isImputation, isFirth,
                    as.integer(nThreadsMarker))

    # timeoutput=system.time({writeOutputFile(Output = list(resMarker),
    # if(nrow(resMarker) > 0){

    isEnd_Output <- (i == nChunks)

    # }

    # LEDGER #47. Flush the C++ stream FIRST, then record how big the output is,
    # then write the index line: that order is what makes the recorded count
    # mean "everything on disk belongs to chunks <= i". Without the flush the
    # stream would still be holding this chunk's tail, and the count would be an
    # underestimate that a restart would use to discard good rows.
    flushOutfile_singleInCPP()
    outputBytes <- if (file.exists(OutputFile)) file.size(OutputFile) else NA_real_

    writeOutputFileIndex(
      OutputFileIndex = OutputFileIndex,
      AnalysisType = "Marker",
      nEachChunk = format(nMarkersEachChunk, scientific = FALSE),
      indexChunk = i,
      Start = (i == 1),
      End = isEnd_Output,
      outputBytes = outputBytes
    )

    # writeOutputFile(Output = list(resMarker),
    #                  OutputFile = list(OutputFile),
    #                  OutputFileIndex = OutputFileIndex,
    #                  AnalysisType = "Marker",
    #                  nEachChunk = format(nMarkersEachChunk, scientific=F),
    #                  indexChunk = i,
    #                  Start = (i==1),
    #                  End = isEnd_Output)

    # End = (i==nChunks))})
    # print("timeoutput")
    # print(timeoutput)
    ptm <- proc.time()
    print(ptm)
    # (per-chunk gc() removed in Phase 4: it was the top frame of the marker-loop profile)
    # rm(resMarker)


    i <- i + 1
    if (i > nChunks) {
      is_marker_test <- FALSE
    }
  } # while(is_marker_test){

  # Flush and close g_out.single explicitly. This is byte-neutral -- it only
  # moves the flush that static destruction would have done at exit -- and is
  # kept now that the forked path is gone, because it makes the stream's
  # lifetime explicit. It was MANDATORY under fork (experiment E1: mcexit ends a
  # child with _exit(), no static destructor runs, and the child silently lost
  # everything past the 192 KiB buffer boundary while still returning success).
  # That same buffer is the mechanism behind LEDGER #47. See R/MarkerScan.R.
  closeOutfile_single()

  # information to users
  output <- paste0("Analysis done! The results have been saved to '", OutputFile, "'.")

  output
}


splitMarker <- function(genoIndex, nMarkersEachChunk, CHROM) {
  genoIndexList <- list()
  iTot <- 1

  uCHROM <- unique(CHROM)
  for (chrom in uCHROM) {
    pos <- which(CHROM == chrom)
    gIdx <- genoIndex[pos]
    M <- length(gIdx)

    idxStart <- seq(1, M, nMarkersEachChunk)
    idxEnd <- idxStart + nMarkersEachChunk - 1

    nChunks <- length(idxStart)
    idxEnd[nChunks] <- M

    for (i in 1:nChunks) {
      idxMarker <- idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] <- list(
        chrom = chrom,
        genoIndex = gIdx[idxMarker]
      )

      iTot <- iTot + 1
    }
  }

  genoIndexList
}
