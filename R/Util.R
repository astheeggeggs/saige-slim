checkArgsBool = function(arg0, arg0name){
  if(!arg0 %in% c(TRUE, FALSE))
    stop(arg0name, " should be TRUE or FALSE.")
}

checkArgsNumeric = function(arg0, arg0name, minVal, maxVal){
  if(!is.numeric(arg0) | arg0 < minVal | arg0 > maxVal)
    stop(arg0name, " should be a numeric value ranging from 0 to 0.5.")
}



Check_File_Exist<-function(file, filetype){
    if (!file.exists(file)) {
        stop("ERROR! ", filetype , " ", file, " does not exsit\n")
    }
}

Check_OutputFile_Create<-function(file){
  if (file.exists(file)) {
    file.remove(file)
  }
  if (!file.exists(file)) {
    file.create(file, showWarnings = TRUE)
  }
}	

###Step 1
updateChrStartEndIndexVec <- function(chrVec){
 chromosomeStartIndexVec = NULL
 chromosomeEndIndexVec = NULL
 LOCO = TRUE 
 #cat("WARNING: leave-one-chromosome-out is activated! Note this option will only be applied to autosomal variants\n")
 cat("WARNING: Genetic variants needs to be ordered by chromosome and position in the Plink file\n")
 #bimData = data.table:::fread(paste0(plinkFile,".bim"), header = F)
 for (i in 1:22) {
    if (length(which(chrVec == i)) > 0) {
    	chromosomeStartIndexVec = c(chromosomeStartIndexVec,min(which(chrVec == i)) - 1)
        chromosomeEndIndexVec = c(chromosomeEndIndexVec,max(which(chrVec == i)) - 1)
        if (i > 1) {
       	  if (!is.na(chromosomeStartIndexVec[i -1])) {
            if (chromosomeStartIndexVec[i] <= chromosomeStartIndexVec[i - 1] | chromosomeEndIndexVec[i] <= chromosomeEndIndexVec[i -1]) {
                          stop(paste0("ERROR! chromosomes need to be ordered from 1 to 22 in the bim file.\n"))
            }
          }
        }
     }else{
        chromosomeStartIndexVec = c(chromosomeStartIndexVec, NA)
        chromosomeEndIndexVec = c(chromosomeEndIndexVec, NA)
     }
 }
 cat("chromosomeStartIndexVec: ", chromosomeStartIndexVec, "\n")
 cat("chromosomeEndIndexVec: ", chromosomeEndIndexVec,"\n")
 if (sum(!is.na(chromosomeStartIndexVec)) <= 1 | sum(!is.na(chromosomeEndIndexVec)) <= 1) {
   cat("WARNING: The number of autosomal chromosomes is less than 2 and leave-one-chromosome-out can't be conducted! \n")
   LOCO = FALSE
 }
 chromosomeStartIndexVec_forcpp = chromosomeStartIndexVec
 chromosomeStartIndexVec_forcpp[is.na(chromosomeStartIndexVec_forcpp)] = -1
 chromosomeEndIndexVec_forcpp = chromosomeEndIndexVec
 chromosomeEndIndexVec_forcpp[is.na(chromosomeEndIndexVec_forcpp)] = -1

 setStartEndIndexVec(chromosomeStartIndexVec_forcpp, chromosomeEndIndexVec_forcpp)
 return(list(LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec, chromosomeEndIndexVec = chromosomeEndIndexVec))
}
checkOutputFile = function(OutputFile,
                           OutputFileIndex,
                           AnalysisType,      ## "Marker" or "Region"
                           nEachChunk,
			   isOverwriteOutput)
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  # message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."

  ## an R list of output
  if(OutputFile == "")
    stop("Argument of 'OutputFile' is required.")

if(isOverwriteOutput){
    Check_OutputFile_Create(OutputFile)
    #Check_OutputFile_Create(OutputFileIndex) 
    Start = TRUE;
    End = FALSE;
    indexChunk = 0;
    if(AnalysisType == "Marker"){
      indexChunk = 1;
    }
  }else{



  if(file.exists(OutputFile)){
    if(!file.exists(OutputFileIndex)){
      stop(paste0("'OutputFile' of '", OutputFile,"' has existed.
                  Please use another 'OutputFile' or specify --is_overwrite_output=TRUE to overwrite the OutputFile."))
    }else{
      outIndexData = read.table(OutputFileIndex, header = F, sep="\t")

      if(outIndexData[1,1] != message1 | outIndexData[2,1] != message2 | outIndexData[3,1] != message3)
        stop(paste0("'OutputFileIndex' of '", OutputFileIndex, "' is not as expected.
                    Probably, it has been modified by user, which is not permitted.
                    Please remove the existing files of 'OutputFile' and 'OutputFileIndex' or specify --is_overwrite_output=TRUE to overwrite the OutputFile."))

      lastMessage = outIndexData[nrow(outIndexData), 1]
      if(lastMessage == message5){
        End = TRUE
        indexChunk = outIndexData[nrow(outIndexData)-1, 1];
        indexChunk = as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        cat("Based on 'OutputFile' and 'OutputFileIndex', the analysis has been completed for the toal", indexChunk, "chunks.\n")
      }else{
        End = FALSE;
        indexChunk = lastMessage;
        indexChunk = as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        cat("Based on 'OutputFile' and 'OutputFileIndex', we restart the analysis from the", indexChunk+1, "chunk.\n")
      }
    }
    Start = FALSE
  }else{
    Check_OutputFile_Create(OutputFile)	  
    #Check_OutputFile_Create(OutputFileIndex)	  
    Start = TRUE;
    End = FALSE;
    indexChunk = 0;
    if(AnalysisType == "Marker"){
      indexChunk = 1;
    }	    
  }
}


  returnList = list(Start = Start, End = End, indexChunk = indexChunk)
  return(returnList)
}



writeOutputFile = function(Output,
                           OutputFile,
                           OutputFileIndex,
                           AnalysisType,
                           nEachChunk,
                           indexChunk,
                           Start,   # TRUE or FALSE, to indicate is the 'Output' is the first one to save into 'OutputFile'
                           End)     # TRUE or FALSE, to indicate is the 'Output' is the last one to save into 'OutputFile'
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."

  n1 = length(Output)
  n2 = length(OutputFile)
  cat("write to output\n")
  #cat("n1 is ", n1, "\n")
  #cat("n2 is ", n2, "\n")
  if(n1 != n2)
    stop("length(Output) != length(OutputFile)")

  if(n1 != 0){
    for(i in 1:n1){
      if(Start){
	if(!is.null(Output[[i]])){
          write.table(Output[[i]], OutputFile[[i]], quote = F, sep = "\t", append = F, col.names = T, row.names = F)
	}
      }else{
	if(!is.null(Output[[i]])){
          write.table(Output[[i]], OutputFile[[i]], quote = F, sep = "\t", append = T, col.names = F, row.names = F)
	}
        #write.table(Output, OutputFile, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
    }
  }
  #print("write Output 2")
  if(Start)
    write.table(c(message1, message2, message3), OutputFileIndex,
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)
  write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)

  if(End)
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
}


isContainValue = function(x, val){
	return(val %in% x)
}


writeOutputFileIndex = function(
                           OutputFileIndex,
                           AnalysisType,
                           nEachChunk,
                           indexChunk,
                           Start,   # TRUE or FALSE, to indicate is the 'Output' is the first one to save into 'OutputFile'
                           End)     # TRUE or FALSE, to indicate is the 'Output' is the last one to save into 'OutputFile'
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."

  cat("write to output\n")

  #print("write Output 2")
  if(Start)
    write.table(c(message1, message2, message3), OutputFileIndex,
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)
  write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)

  if(End)
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
}


convert_str_to_log = function(a){
	if(a != "NA"){
		b = as.numeric(unlist(strsplit(a, split="E")))
		b1 = b[2]+log10(b[1])
		d = b1 * log(10)
	}else{
		d = NA	
	}	
}	
