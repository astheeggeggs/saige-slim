SAIGE.Marker = function(traitType,
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
			isAnyInclude)
{

  if(is.null(OutputFileIndex))
    OutputFileIndex = paste0(OutputFile, ".index")
  
  #genoType = objGeno$genoType

  outIndex = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F), isOverWriteOutput)    # this function is in 'Util.R'
  outIndex = outIndex$indexChunk
  isappend = FALSE
  if(outIndex != 1){
    cat("Restart the analysis from chunk:\t", outIndex, "\n")
    isappend = TRUE
  }  

  isOpenOutFile_single = openOutfile_single(traitType, isImputation, isappend, isMoreOutput)

  if(!isOpenOutFile_single){
    stop("Output file ", OutputFile, " can't be opened\n")
  }

  ## set up an object for genotype
  if(genoType == "plink"){
      #markerInfo = objGeno$markerInfo
    if(LOCO){
      genoIndex = genoIndex[which(CHROM == chrom)]
      CHROM = CHROM[which(CHROM == chrom)]  
      #markerInfo = markerInfo[which(markerInfo$CHROM == chrom),]  
    }
    #CHROM = markerInfo$CHROM
    #genoIndex = markerInfo$genoIndex
    ##only for one chrom
    # all markers were split into multiple chunks,
    #print(markerInfo[1:10,])
    #print(genoIndex[1:10])    

    genoIndexList = splitMarker(genoIndex, nMarkersEachChunk, CHROM);
    nChunks = length(genoIndexList)

    if(nChunks == 0){
          stop("No markers on chrom ", chrom, " are found\n")
    }


    cat("Number of all markers to test:\t", length(genoIndex), "\n")
    cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
    cat("Number of chunks for all markers:\t", nChunks, "\n")
    if(outIndex > nChunks){
      cat("The analysis has been finished! Please delete ", OutputFileIndex, " if the analysis needs to be run again or set --is_overwrite_output=TRUE\n")
      is_marker_test = FALSE 
    }else{
      is_marker_test = TRUE
      i = outIndex    
    }
    
  }
  #chrom = "InitialChunk"
  #set_flagSparseGRM_cur_SAIGE_org()
  while(is_marker_test){
  #for(i in outIndex:nChunks)
  #{
#time_left = system.time({


    if(genoType == "plink"){	
      tempList = genoIndexList[[i]]
      genoIndex = as.character(format(tempList$genoIndex, scientific = FALSE))
      tempChrom = tempList$chrom
    }



    #print("tempList")
    #print(tempList)
    #print(tempList$genoIndex)
#})
#print("time_left")
#print(time_left)
    #print("genoIndex here")
    #print(genoIndex)
    # set up objects that do not change for different variants
    #if(tempChrom != chrom){
    #  setMarker("SAIGE", objNull, control, chrom, Group, ifOutGroup)
    #  chrom = tempChrom
    #}
    #ptm <- proc.time()
    #print(ptm)
    #print("gc()")
    #print(gc())
    cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \n"))
    # main function to calculate summary statistics for markers in one chunk
   
   #resMarker = as.data.frame(mainMarkerInCPP(genoType, traitType, genoIndex_prev, genoIndex, isMoreOutput, isImputation)) 
   #resMarker = resMarker[which(!is.na(resMarker$BETA)), ]

  mainMarkerInCPP(genoType, traitType, genoIndex, isMoreOutput, isImputation, isFirth)

    #timeoutput=system.time({writeOutputFile(Output = list(resMarker),
  #if(nrow(resMarker) > 0){

  isEnd_Output = (i==nChunks)
  
  #}

  writeOutputFileIndex(OutputFileIndex = OutputFileIndex,
			AnalysisType = "Marker",
			nEachChunk = format(nMarkersEachChunk, scientific=F),
			indexChunk = i,
			Start = (i==1),
			End = isEnd_Output)

  #writeOutputFile(Output = list(resMarker),
  #                  OutputFile = list(OutputFile),
  #                  OutputFileIndex = OutputFileIndex,
  #                  AnalysisType = "Marker",
  #                  nEachChunk = format(nMarkersEachChunk, scientific=F),
  #                  indexChunk = i,
  #                  Start = (i==1),
  #                  End = isEnd_Output)

                    #End = (i==nChunks))})
    #print("timeoutput")
    #print(timeoutput)
    ptm <- proc.time()
    print(ptm)
    # (per-chunk gc() removed in Phase 4: it was the top frame of the marker-loop profile)
    #rm(resMarker)


    
    i = i + 1
    if(i > nChunks){
      is_marker_test = FALSE
    }
	  
  } #while(is_marker_test){

  # information to users
  output = paste0("Analysis done! The results have been saved to '", OutputFile,"'.")

  return(output)
}


splitMarker = function(genoIndex, nMarkersEachChunk, CHROM)
{
  genoIndexList = list()
  iTot = 1;

  uCHROM = unique(CHROM)
  for(chrom in uCHROM){
    pos = which(CHROM == chrom)
    gIdx = genoIndex[pos]
    M = length(gIdx)

    idxStart = seq(1, M, nMarkersEachChunk)
    idxEnd = idxStart + nMarkersEachChunk - 1

    nChunks = length(idxStart)
    idxEnd[nChunks] = M

    for(i in 1:nChunks){
      idxMarker = idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] = list(chrom = chrom,
                                   genoIndex = gIdx[idxMarker])

      iTot = iTot + 1;
    }
  }

  return(genoIndexList)
}
