checkArgsBool <- function(arg0, arg0name) {
  if (!arg0 %in% c(TRUE, FALSE)) {
    stop(arg0name, " should be TRUE or FALSE.")
  }
}

checkArgsNumeric <- function(arg0, arg0name, minVal, maxVal) {
  if (!is.numeric(arg0) | arg0 < minVal | arg0 > maxVal) {
    stop(arg0name, " should be a numeric value ranging from ", minVal, " to ", maxVal, ".")
  }
}


Check_File_Exist <- function(file, filetype) {
  if (!file.exists(file)) {
    stop("ERROR! ", filetype, " ", file, " does not exsit\n")
  }
}

Check_OutputFile_Create <- function(file) {
  if (file.exists(file)) {
    file.remove(file)
  }
  if (!file.exists(file)) {
    file.create(file, showWarnings = TRUE)
  }
}

### Step 1
updateChrStartEndIndexVec <- function(chrVec) {
  chromosomeStartIndexVec <- NULL
  chromosomeEndIndexVec <- NULL
  LOCO <- TRUE
  # cat("WARNING: leave-one-chromosome-out is activated! Note this option will only be applied to autosomal variants\n")
  cat("WARNING: Genetic variants needs to be ordered by chromosome and position in the Plink file\n")
  # bimData = data.table:::fread(paste0(plinkFile,".bim"), header = F)
  for (i in 1:22) {
    if (length(which(chrVec == i)) > 0) {
      chromosomeStartIndexVec <- c(chromosomeStartIndexVec, min(which(chrVec == i)) - 1)
      chromosomeEndIndexVec <- c(chromosomeEndIndexVec, max(which(chrVec == i)) - 1)
      # LEDGER #62. Each chromosome is carried downstream as the single window
      # [min, max] of its marker indices, so the windows are only the
      # chromosomes if every chromosome's markers are CONTIGUOUS. The check
      # below tests monotonicity of consecutive starts and ends, which catches
      # a marker displaced past a whole chromosome but not one displaced INTO
      # the middle of the next: chrVec = c(rep(1, 10), rep(2, 5), 1, rep(2, 5),
      # rep(3, 10)) passes it, and chr1's window then spans 16 markers for its
      # 11, swallowing five chr2 markers. The leave-one-out property survives
      # -- it is over-exclusion, never contamination -- but those five markers
      # are silently in neither LOCO GRM. Contiguity is the invariant the
      # windows actually need, so it is what is asserted.
      nOnChr <- length(which(chrVec == i))
      spanOnChr <- max(which(chrVec == i)) - min(which(chrVec == i)) + 1
      if (spanOnChr != nOnChr) {
        stop(paste0(
          "ERROR! the markers on chromosome ", i, " are not contiguous in the bim file: ",
          nOnChr, " markers span ", spanOnChr, " positions, so ", spanOnChr - nOnChr,
          " marker(s) from other chromosomes lie between the first and last of them.",
          " Leave-one-chromosome-out uses one [start, end] window per chromosome, so those",
          " markers would be dropped from more than one LOCO GRM. Sort the bim file by",
          " chromosome and position.\n"
        ))
      }
      if (i > 1) {
        if (!is.na(chromosomeStartIndexVec[i - 1])) {
          if (chromosomeStartIndexVec[i] <= chromosomeStartIndexVec[i - 1] | chromosomeEndIndexVec[i] <= chromosomeEndIndexVec[i - 1]) {
            stop(paste0("ERROR! chromosomes need to be ordered from 1 to 22 in the bim file.\n"))
          }
        }
      }
    } else {
      chromosomeStartIndexVec <- c(chromosomeStartIndexVec, NA)
      chromosomeEndIndexVec <- c(chromosomeEndIndexVec, NA)
    }
  }
  cat("chromosomeStartIndexVec: ", chromosomeStartIndexVec, "\n")
  cat("chromosomeEndIndexVec: ", chromosomeEndIndexVec, "\n")
  if (sum(!is.na(chromosomeStartIndexVec)) <= 1 | sum(!is.na(chromosomeEndIndexVec)) <= 1) {
    cat("WARNING: The number of autosomal chromosomes is less than 2 and leave-one-chromosome-out can't be conducted! \n")
    LOCO <- FALSE
  }
  chromosomeStartIndexVec_forcpp <- chromosomeStartIndexVec
  chromosomeStartIndexVec_forcpp[is.na(chromosomeStartIndexVec_forcpp)] <- -1
  chromosomeEndIndexVec_forcpp <- chromosomeEndIndexVec
  chromosomeEndIndexVec_forcpp[is.na(chromosomeEndIndexVec_forcpp)] <- -1

  setStartEndIndexVec(chromosomeStartIndexVec_forcpp, chromosomeEndIndexVec_forcpp)
  list(LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec, chromosomeEndIndexVec = chromosomeEndIndexVec)
}
# The per-chunk progress line of a marker `.index` file, parsed.
#
# Two formats exist and both must read. The original is
#   "Have completed the analysis of chunk 7"
# and carries no way to tell how much of the output file belongs to chunk 7.
# Since the LEDGER #47 fix it is
#   "Have completed the analysis of chunk 7 bytes 123456"
# where the byte count is the size of the output file at the moment chunk 7's
# rows had all been flushed. An index written by an older version parses with
# bytes = NA, which restartOutputBytes() below treats as "cannot verify".
#
# LEDGER #140 (region path): a region run has TWO files a restart must cut back,
# the results file and `.singleAssoc.txt`, so its line carries both counts:
#   "Have completed the analysis of chunk 7 single_bytes 98765 bytes 123456"
# `single_bytes` sits BEFORE `bytes` so the `bytes N$` anchor above is unchanged;
# a line without it (a run that wrote no single-variant file, or an older
# version) parses with single = NA. indexChunkLine() below is the one writer of
# this format and this function its one reader; keep them together.
parseIndexChunkLine <- function(line) {
  chunk <- suppressWarnings(as.numeric(
    sub("^Have completed the analysis of chunk ([0-9]+).*$", "\\1", line)
  ))
  bytes <- NA_real_
  if (grepl(" bytes [0-9]+$", line)) {
    bytes <- suppressWarnings(as.numeric(sub("^.* bytes ([0-9]+)$", "\\1", line)))
  }
  single <- NA_real_
  if (grepl(" single_bytes [0-9]+ ", line)) {
    single <- suppressWarnings(as.numeric(sub("^.* single_bytes ([0-9]+) .*$", "\\1", line)))
  }
  list(chunk = chunk, bytes = bytes, single = single)
}

# The per-chunk index line. NULL or NA counts are omitted, which is the legacy
# line (no byte count) that parseIndexChunkLine() reads as "cannot verify".
indexChunkLine <- function(indexChunk, outputBytes = NULL, singleBytes = NULL) {
  line <- paste("Have completed the analysis of chunk", indexChunk)
  hasBytes <- !is.null(outputBytes) && !is.na(outputBytes)
  if (hasBytes && !is.null(singleBytes) && !is.na(singleBytes)) {
    line <- paste(line, "single_bytes", format(singleBytes, scientific = FALSE))
  }
  if (hasBytes) {
    line <- paste(line, "bytes", format(outputBytes, scientific = FALSE))
  }
  line
}

# Drop anything the interrupted run wrote past the last COMPLETED chunk.
#
# LEDGER #47: the output stream flushes on its own ~192 KiB boundary and the
# index once per chunk, so a killed run can leave rows on disk that belong to
# chunks the index never recorded. Appending from the indexed chunk then writes
# those rows a second time -- duplication, never loss, which is benign for a
# scan you re-run and corrupting for anything that joins or counts rows.
#
# Truncating to the recorded byte count is what makes a restart produce the same
# file a complete run would. An index without a byte count (written by an older
# version) cannot be checked: say so and carry on with the old behaviour rather
# than refusing to restart at all.
truncateOutputToIndex <- function(OutputFile, outputBytes) {
  if (!file.exists(OutputFile)) {
    return(invisible(FALSE))
  }
  onDisk <- file.size(OutputFile)
  if (is.null(outputBytes) || is.na(outputBytes)) {
    if (onDisk > 0) {
      warning(
        "'", OutputFile, "' is being appended to from an index that records no ",
        "byte count (written by an older version of SAIGE). Any rows the ",
        "interrupted run wrote past its last completed chunk will be ",
        "DUPLICATED. Re-run from the start with --is_overwrite_output=TRUE to ",
        "be sure of the output.",
        call. = FALSE
      )
    }
    return(invisible(FALSE))
  }
  if (onDisk < outputBytes) {
    warning(
      "'", OutputFile, "' is ", onDisk, " bytes but its index records ",
      outputBytes, " for the last completed chunk. The output file is shorter ",
      "than the index claims, so the restart cannot be trusted; re-run with ",
      "--is_overwrite_output=TRUE.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  if (onDisk == outputBytes) {
    return(invisible(FALSE))
  }
  con <- file(OutputFile, open = "r+b")
  on.exit(close(con))
  seek(con, where = outputBytes, origin = "start")
  truncate(con)
  cat("Discarded", onDisk - outputBytes, "bytes written past the last completed chunk of",
      OutputFile, "\n")
  invisible(TRUE)
}

checkOutputFile <- function(OutputFile,
                            OutputFileIndex,
                            AnalysisType, ## "Marker" or "Region"
                            nEachChunk,
                            isOverwriteOutput) {
  ## The following messages are for 'OutputFileIndex'
  message1 <- "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 <- paste("This is a", AnalysisType, "level analysis.")
  message3 <- paste("nEachChunk =", nEachChunk)
  # message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 <- "Have completed the analyses of all chunks."

  ## an R list of output
  if (OutputFile == "") {
    stop("Argument of 'OutputFile' is required.")
  }

  chunkBytes <- numeric(0)
  chunkSingleBytes <- numeric(0)
  if (isOverwriteOutput) {
    Check_OutputFile_Create(OutputFile)
    # LEDGER #145 FIXED: the region sidecars are written only when a run has
    # something to disclose (.pooledTests.txt, .skatFailures.txt) or asks for
    # one (.markerList.txt), so a previous run's copies survived an overwrite
    # that produced none and described the wrong run. Removed with the
    # results file; unlink() of an absent file is a no-op.
    if (AnalysisType == "Region") {
      unlink(paste0(OutputFile, c(".pooledTests.txt", ".skatFailures.txt", ".markerList.txt")))
    }
    # Check_OutputFile_Create(OutputFileIndex)
    Start <- TRUE
    End <- FALSE
    indexChunk <- 0
    if (AnalysisType == "Marker") {
      indexChunk <- 1
    }
  } else {
    if (file.exists(OutputFile)) {
      if (!file.exists(OutputFileIndex)) {
        stop(paste0("'OutputFile' of '", OutputFile, "' has existed.
                  Please use another 'OutputFile' or specify --is_overwrite_output=TRUE to overwrite the OutputFile."))
      } else {
        outIndexData <- read.table(OutputFileIndex, header = FALSE, sep = "\t")

        if (outIndexData[1, 1] != message1 | outIndexData[2, 1] != message2 | outIndexData[3, 1] != message3) {
          stop(paste0("'OutputFileIndex' of '", OutputFileIndex, "' is not as expected.
                    Probably, it has been modified by user, which is not permitted.
                    Please remove the existing files of 'OutputFile' and 'OutputFileIndex' or specify --is_overwrite_output=TRUE to overwrite the OutputFile."))
        }

        # Every per-chunk line, not just the last: a restart needs the size the
        # output file had at the END OF THE CHUNK BEFORE the one it resumes at
        # (see SAIGE.Marker -- it resumes AT indexChunk, re-running it), so the
        # caller has to be able to look any chunk up.
        chunkLines <- grep("^Have completed the analysis of chunk ",
                           outIndexData[, 1], value = TRUE)
        parsedAll <- lapply(chunkLines, parseIndexChunkLine)
        chunkBytes <- vapply(parsedAll, function(x) x$bytes, numeric(1))
        names(chunkBytes) <- vapply(parsedAll, function(x) as.character(x$chunk),
                                    character(1))
        # LEDGER #140: the region path's second file, `.singleAssoc.txt`
        chunkSingleBytes <- vapply(parsedAll, function(x) x$single, numeric(1))
        names(chunkSingleBytes) <- names(chunkBytes)

        lastMessage <- outIndexData[nrow(outIndexData), 1]
        if (lastMessage == message5) {
          End <- TRUE
          indexChunk <- parseIndexChunkLine(outIndexData[nrow(outIndexData) - 1, 1])$chunk
          cat("Based on 'OutputFile' and 'OutputFileIndex', the analysis has been completed for the toal", indexChunk, "chunks.\n")
        } else {
          End <- FALSE
          indexChunk <- parseIndexChunkLine(lastMessage)$chunk
          cat("Based on 'OutputFile' and 'OutputFileIndex', we restart the analysis from the", indexChunk + 1, "chunk.\n")
          # LEDGER #49 FIXED. For a MARKER run indexChunk is "the chunk to start
          # at" (a fresh run returns 1), so on a restart it is the completed
          # chunk PLUS ONE, as the message above has always said. It used to
          # return the completed chunk itself, and SAIGE.Marker resumed AT it,
          # recomputing one whole chunk per restart (correct output, wasted
          # work). The region path keeps its own convention (LEDGER #81/#25).
          if (AnalysisType == "Marker") {
            indexChunk <- indexChunk + 1
          }
        }
      }
      Start <- FALSE
    } else {
      Check_OutputFile_Create(OutputFile)
      # Check_OutputFile_Create(OutputFileIndex)
      Start <- TRUE
      End <- FALSE
      indexChunk <- 0
      if (AnalysisType == "Marker") {
        indexChunk <- 1
      }
    }
  }


  # chunkBytes maps chunk number -> size the output file had when that chunk
  # finished (LEDGER #47), as a named numeric vector; empty for a fresh run and
  # all-NA for an index written before that record existed.
  returnList <- list(Start = Start, End = End, indexChunk = indexChunk,
                     chunkBytes = chunkBytes, chunkSingleBytes = chunkSingleBytes)
  returnList
}


writeOutputFile <- function(Output,
                            OutputFile,
                            OutputFileIndex,
                            AnalysisType,
                            nEachChunk,
                            indexChunk,
                            Start, # TRUE or FALSE, to indicate is the 'Output' is the first one to save into 'OutputFile'
                            End) { # TRUE or FALSE, to indicate is the 'Output' is the last one to save into 'OutputFile'
  ## The following messages are for 'OutputFileIndex'
  message1 <- "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 <- paste("This is a", AnalysisType, "level analysis.")
  message3 <- paste("nEachChunk =", nEachChunk)
  message4 <- paste("Have completed the analysis of chunk", indexChunk)
  message5 <- "Have completed the analyses of all chunks."

  n1 <- length(Output)
  n2 <- length(OutputFile)
  cat("write to output\n")
  # cat("n1 is ", n1, "\n")
  # cat("n2 is ", n2, "\n")
  if (n1 != n2) {
    stop("length(Output) != length(OutputFile)")
  }

  if (n1 != 0) {
    for (i in 1:n1) {
      if (Start) {
        if (!is.null(Output[[i]])) {
          write.table(Output[[i]], OutputFile[[i]], quote = FALSE, sep = "\t", append = FALSE, col.names = TRUE, row.names = FALSE)
        }
      } else {
        if (!is.null(Output[[i]])) {
          write.table(Output[[i]], OutputFile[[i]], quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
        }
        # write.table(Output, OutputFile, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
    }
  }
  # print("write Output 2")
  if (Start) {
    write.table(c(message1, message2, message3), OutputFileIndex,
      quote = FALSE, sep = "\t", append = FALSE, col.names = FALSE, row.names = FALSE
    )
  }
  write.table(message4, OutputFileIndex, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

  if (End) {
    write.table(message5, OutputFileIndex, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
  }
}


isContainValue <- function(x, val) {
  val %in% x
}


writeOutputFileIndex <- function(OutputFileIndex,
                                 AnalysisType,
                                 nEachChunk,
                                 indexChunk,
                                 Start, # TRUE or FALSE, to indicate is the 'Output' is the first one to save into 'OutputFile'
                                 End, # TRUE or FALSE, to indicate is the 'Output' is the last one to save into 'OutputFile'
                                 # Size of the output file with this chunk's rows all flushed. Recorded
                                 # so a restart can truncate to it (LEDGER #47); NULL keeps the legacy
                                 # line, which is what the region path has always written.
                                 outputBytes = NULL,
                                 # LEDGER #140: the region path's `.singleAssoc.txt`
                                 # size at the same moment; NULL for the marker path
                                 singleBytes = NULL) {
  ## The following messages are for 'OutputFileIndex'
  message1 <- "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 <- paste("This is a", AnalysisType, "level analysis.")
  message3 <- paste("nEachChunk =", nEachChunk)
  message4 <- indexChunkLine(indexChunk, outputBytes, singleBytes)
  message5 <- "Have completed the analyses of all chunks."

  cat("write to output\n")

  # print("write Output 2")
  if (Start) {
    write.table(c(message1, message2, message3), OutputFileIndex,
      quote = FALSE, sep = "\t", append = FALSE, col.names = FALSE, row.names = FALSE
    )
  }
  write.table(message4, OutputFileIndex, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

  if (End) {
    write.table(message5, OutputFileIndex, quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)
  }
}


# LEDGER #23 FIXED. Reads every form a p-value reaches R in: the C++ writer's
# "<mantissa>E<exponent>" (kept on the log scale, which is the point -- the
# exponent can be far below what a double holds), a plain decimal, the string
# "NA", and R's own NA. It used to split on "E" only, so "0.5" became
# NA + log10(0.5) = NA with no warning, and NA_character_ errored in the if.
convert_str_to_log <- function(a) {
  if (is.na(a) || a == "NA") {
    return(NA) # the plain NA the original returned for "NA"; callers test is.na
  }
  if (grepl("E", a, fixed = TRUE)) {
    b <- as.numeric(unlist(strsplit(a, split = "E")))
    b1 <- b[2] + log10(b[1])
    return(b1 * log(10))
  }
  log(as.numeric(a))
}


# LEDGER #68 FIXED. Covariate names were pasted into a formula bare, so a name
# containing formula syntax was silently REINTERPRETED: `x1^2` vanished from
# y ~ x1 + x2 + x1^2 (term labels x1, x2 -- checked), `x1*x2` became an
# interaction, and the --covarColList membership check had already confirmed
# the column present. Backtick-quoting makes every name a single term. For an
# ordinary name the model-matrix column names are unchanged, so nothing moves.
quoteCovarNames <- function(x) {
  if (length(x) == 0) return(x)
  paste0("`", gsub("`", "\\`", x), "`")
}
