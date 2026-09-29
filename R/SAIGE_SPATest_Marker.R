# The single-variant marker scan's chunk loop.
#
# M1/M2: this drives P >= 1 null models over ONE genotype read. Everything that
# was a scalar path -- the output file, its .index, the restart check, the C++
# stream -- is now one per MODEL, and everything that was a per-chunk C++ call
# is still exactly one per chunk, because the models are tested inside
# mainMarkerInCPP rather than by calling it P times. For P == 1 this is the
# same call sequence in the same order as before, which is what makes the
# marker goldens the gate for the change.
SAIGE.Marker <- function(genoType, # nolint: object_name_linter.
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
  if (length(OutputFileIndex) != length(OutputFile)) {
    stop("Each output file needs exactly one progress index file\n")
  }
  nModels <- length(OutputFile)

  # checkOutputFile is in 'Util.R'
  outChecks <- lapply(seq_len(nModels), function(k) {
    checkOutputFile(
      OutputFile[k], OutputFileIndex[k], "Marker",
      format(nMarkersEachChunk, scientific = FALSE), isOverWriteOutput
    )
  })
  # There is ONE marker loop, so there is one chunk to resume at. If the models
  # disagree the files were not all written by the same interrupted run, and
  # picking one is how rows get silently duplicated (LEDGER #47's shape).
  outIndex <- restartChunkForModels(outChecks, OutputFileIndex)
  isappend <- FALSE
  if (outIndex != 1) {
    cat("Restart the analysis from chunk:\t", outIndex, "\n")
    isappend <- TRUE
    # LEDGER #47. The interrupted run may have flushed rows belonging to chunks
    # the index never recorded; appending would write them twice. Drop them
    # before the output is reopened (truncateOutputToIndex is in 'Util.R').
    #
    # outIndex is the chunk the loop below STARTS AT (LEDGER #49 FIXED: on a
    # restart checkOutputFile now returns the completed chunk plus one, as its
    # message always said), so the target is the size at the end of the chunk
    # before it -- the last one the index recorded as complete. Before #49 the
    # loop resumed AT the completed chunk and recomputed it; the target was the
    # chunk before THAT, so the output was right and one chunk of work was
    # repeated per restart. The #47 regression test is what pins this.
    #
    # Per file, by ITS OWN recorded byte count: two models write different
    # numbers of bytes per chunk, so one shared count would be wrong for both.
    for (k in seq_len(nModels)) {
      truncateOutputToIndex(
        OutputFile[k],
        outChecks[[k]]$chunkBytes[as.character(outIndex - 1)]
      )
    }
  }

  for (k in seq_len(nModels)) {
    # The header comes from null model k in C++ (trait type and the conditional
    # columns are per model), so this is addressed by model index.
    isOpenOutFile_single <- openOutfile_single(
      as.integer(k - 1L), isImputation, isappend, isMoreOutput
    )
    if (!isOpenOutFile_single) {
      stop("Output file ", OutputFile[k], " can't be opened\n")
    }
  }

  ## set up an object for genotype
  if (genoType %in% c("plink", "pgen")) {
    if (LOCO) {
      genoIndex <- genoIndex[which(CHROM == chrom)]
      CHROM <- CHROM[which(CHROM == chrom)]
    }

    genoIndexList <- splitMarker(genoIndex, nMarkersEachChunk, CHROM)
    nChunks <- length(genoIndexList)

    if (nChunks == 0) {
      stop("No markers on chrom ", chrom, " are found\n")
    }
    resetSpaFallbackRecordsInCPP() # LEDGER #164: per-run accounting, R/SpaFallbacks.R
    resetSpaPinnedRecordsInCPP() # LEDGER #174: per-run accounting, R/SpaPinned.R


    cat("Number of all markers to test:\t", length(genoIndex), "\n")
    cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
    cat("Number of chunks for all markers:\t", nChunks, "\n")
    if (outIndex > nChunks) {
      cat(
        "The analysis has been finished! Please delete ", OutputFileIndex[1],
        " if the analysis needs to be run again or set --is_overwrite_output=TRUE\n"
      )
      is_marker_test <- FALSE
    } else {
      is_marker_test <- TRUE
      i <- outIndex
    }
  }
  while (is_marker_test) {
    if (genoType %in% c("plink", "pgen")) {
      tempList <- genoIndexList[[i]]
      genoIndex <- as.character(format(tempList$genoIndex, scientific = FALSE))
    }

    cat(paste0(
      "(", Sys.time(), ") ---- Analyzing Chunk ", i, "/", nChunks,
      ": chrom ", chrom, " ---- \n"
    ))
    # main function to calculate summary statistics for markers in one chunk.
    # ONE call however many models: the models are the inner loop, inside the
    # marker loop, which is the whole point.

    mainMarkerInCPP(genoType, genoIndex, isMoreOutput, isImputation, isFirth,
                    as.integer(nThreadsMarker))

    isEnd_Output <- (i == nChunks)

    # LEDGER 47. Flush the C++ streams FIRST, then record how big each output
    # is, then write the index lines: that order is what makes the recorded
    # count mean "everything on disk belongs to chunks <= i". Without the flush
    # a stream would still be holding this chunk's tail, and the count would be
    # an underestimate that a restart would use to discard good rows.
    flushOutfile_singleInCPP()
    for (k in seq_len(nModels)) {
      outputBytes <- if (file.exists(OutputFile[k])) file.size(OutputFile[k]) else NA_real_

      writeOutputFileIndex(
        OutputFileIndex = OutputFileIndex[k],
        AnalysisType = "Marker",
        nEachChunk = format(nMarkersEachChunk, scientific = FALSE),
        indexChunk = i,
        Start = (i == 1),
        End = isEnd_Output,
        outputBytes = outputBytes
      )
    }

    ptm <- proc.time()
    print(ptm)
    # (per-chunk gc() removed in Phase 4: it was the top frame of the marker-loop profile)

    i <- i + 1
    if (i > nChunks) {
      is_marker_test <- FALSE
    }
  } # while(is_marker_test){
  spaFallbackSummary(OutputFile) # LEDGER #164: once per run, R/SpaFallbacks.R
  spaPinnedSummary(OutputFile) # LEDGER #174

  # Flush and close every g_out.singles stream explicitly. This is byte-neutral
  # -- it only moves the flush that static destruction would have done at exit
  # -- and is kept now that the forked path is gone, because it makes the
  # streams' lifetime explicit. It was MANDATORY under fork (experiment E1:
  # mcexit ends a child with _exit(), no static destructor runs, and the child
  # silently lost everything past the 192 KiB buffer boundary while still
  # returning success). That same buffer is the mechanism behind LEDGER #47.
  # See R/MarkerScan.R.
  closeOutfile_single()

  # information to users
  output <- paste0(
    "Analysis done! The results have been saved to '",
    paste(OutputFile, collapse = "', '"), "'."
  )

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
