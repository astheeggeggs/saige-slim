checkGenoInput <- function(bedFile = "",
                           bimFile = "",
                           famFile = "",
                           pgenFile = "",
                           pvarFile = "",
                           psamFile = "",
                           sampleInModel = NULL) {
  if (is.null(sampleInModel)) {
    stop("sampleInModel is not specified.")
  }
  # genotype data
  if (bedFile != "" && pgenFile != "") {
    stop(
      "bedFile and pgenFile are mutually exclusive: specify PLINK ",
      "bed/bim/fam OR PGEN pgen/pvar/psam, not both\n"
    )
  } else if (bedFile != "") {
    Check_File_Exist(bedFile, "bedFile")
    if (bimFile == "") {
      bimFile <- gsub("bed$", "bim", bedFile)
    }
    if (famFile == "") {
      famFile <- gsub("bed$", "fam", bedFile)
    }
    Check_File_Exist(bimFile, "bimFile")
    Check_File_Exist(famFile, "famFile")
    dosageFileType <- "plink"
  } else if (pgenFile != "") {
    Check_File_Exist(pgenFile, "pgenFile")
    if (pvarFile == "") {
      pvarFile <- gsub("pgen$", "pvar", pgenFile)
    }
    if (psamFile == "") {
      psamFile <- gsub("pgen$", "psam", pgenFile)
    }
    Check_File_Exist(pvarFile, "pvarFile")
    Check_File_Exist(psamFile, "psamFile")
    dosageFileType <- "pgen"
  } else {
    stop(
      "Either bedFile (PLINK bed/bim/fam) or pgenFile (PGEN pgen/pvar/psam) ",
      "is required as genotype input on this branch\n"
    )
  }

  cat("dosageFile type is ", dosageFileType, "\n")

  dosageFileType
}

setGenoInput <- function(bedFile = "",
                         bimFile = "",
                         famFile = "",
                         pgenFile = "",
                         pvarFile = "",
                         psamFile = "",
                         idstoIncludeFile = "",
                         rangestoIncludeFile = "",
                         chrom = "",
                         AlleleOrder = NULL,
                         sampleInModel = NULL) {
  # time_geno_0 = proc.time()

  dosageFileType <- checkGenoInput(
    bedFile = bedFile,
    bimFile = bimFile,
    famFile = famFile,
    pgenFile = pgenFile,
    pvarFile = pvarFile,
    psamFile = psamFile,
    sampleInModel = sampleInModel
  )

  # time_geno_1 = proc.time()

  # cat("time_geno_1 - time_geno_0\n")
  # print(time_geno_1 - time_geno_0)

  ########## ----------  Plink format ---------- ##########

  if (dosageFileType == "plink") {
    if (is.null(AlleleOrder)) AlleleOrder <- "alt-first"

    cat("allele order in the plink file is ", AlleleOrder, ".\n")

    if (bimFile == "") {
      bimFile <- gsub("bed$", "bim", bedFile)
    }
    if (famFile == "") {
      famFile <- gsub("bed$", "fam", bedFile)
    }
    markerInfo <- data.table::fread(bimFile, header = FALSE, select = c(1, 2, 4, 5, 6))
    # markerInfo = as.data.frame(markerInfo)
    # markerInfo is a data.table
    if (AlleleOrder == "alt-first") {
      names(markerInfo) <- c("CHROM", "ID", "POS", "ALT", "REF")
    }
    if (AlleleOrder == "ref-first") {
      names(markerInfo) <- c("CHROM", "ID", "POS", "REF", "ALT")
    }
    markerInfo$genoIndex <- 1:nrow(markerInfo) - 1 # -1 is to convert 'R' to 'C++'
    if (chrom != "") {
      markerInfo <- markerInfo[which(markerInfo[, 1] == chrom), ]
    }
    markerInfo[, REF := NULL]
    markerInfo[, ALT := NULL]
    setkeyv(markerInfo, c("ID"))
    setPLINKobjInCPP(bimFile, famFile, bedFile, sampleInModel, AlleleOrder)
  }

  ########## ----------  PGEN format ---------- ##########

  if (dosageFileType == "pgen") {
    if (pvarFile == "") {
      pvarFile <- gsub("pgen$", "pvar", pgenFile)
    }
    if (psamFile == "") {
      psamFile <- gsub("pgen$", "psam", pgenFile)
    }
    # pvar is always REF/ALT-ordered -- unlike bim, there is no
    # AlleleOrder-dependent column swap to apply here.
    markerInfo <- data.table::fread(pvarFile, header = TRUE, skip = "#CHROM")
    names(markerInfo)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex <- 1:nrow(markerInfo) - 1 # -1 is to convert 'R' to 'C++'
    if (chrom != "") {
      markerInfo <- markerInfo[which(markerInfo[, 1] == chrom), ]
    }
    markerInfo[, REF := NULL]
    markerInfo[, ALT := NULL]
    # Match the PLINK path's column order (CHROM, ID, POS, genoIndex) so
    # downstream code sees the same markerInfo shape regardless of genoType.
    setcolorder(markerInfo, c("CHROM", "ID", "POS", "genoIndex"))
    setkeyv(markerInfo, c("ID"))
    setPGENobjInCPP(pgenFile, pvarFile, psamFile, sampleInModel)
  }

  # time_geno_2 = proc.time()

  # cat("time_geno_2 - time_geno_1\n")
  #  print(time_geno_2 - time_geno_1)

  # Files = c("idstoIncludeFile", "idstoExcludeFile", "rangestoIncludeFile", "rangestoExcludeFile")

  anyInclude <- FALSE
  # anyExclude = FALSE

  markersInclude <- c()
  # markersExclude = c()
  IDsToInclude <- NULL
  RangesToInclude <- NULL
  if (idstoIncludeFile != "") {
    IDsToInclude <- data.table::fread(idstoIncludeFile, header = FALSE, colClasses = c("character"), data.table = FALSE)
    if (ncol(IDsToInclude) != 1) {
      stop("'idstoIncludeFile' of ", idstoIncludeFile, " should only include one column.")
    } else {
      IDsToInclude <- IDsToInclude[, 1]
      cat(length(IDsToInclude), " marker IDs are found in the idstoIncludeFile file\n")
    }

    if (dosageFileType %in% c("plink", "pgen")) {
      if (!is.null(markerInfo$ID2)) {
        posRows <- which((markerInfo$ID %in% IDsToInclude) | (markerInfo$ID2 %in% IDsToInclude))
      } else {
        posRows <- which((markerInfo$ID %in% IDsToInclude))
      }
      if (length(posRows) != 0) {
        cat(length(posRows), " markers in idstoIncludeFile are found in the geno/dosage file.\n")
        markersInclude <- c(markersInclude, markerInfo$ID[posRows])
      }
    }

    anyInclude <- TRUE
  }

  if (rangestoIncludeFile != "") {
    RangesToInclude <- data.table::fread(rangestoIncludeFile, header = FALSE, colClasses = c("character", "numeric", "numeric"), data.table = FALSE)
    if (ncol(RangesToInclude) != 3) {
      stop("rangestoIncludeFile should only include three columns.")
    } else {
      cat(nrow(RangesToInclude), " ranges are in rangestoIncludeFile\n")
    }

    colnames(RangesToInclude) <- c("CHROM", "START", "END")
    # LEDGER #11 fix: the ranges were previously never applied on the PLINK
    # path (anyInclude was set with an empty include list, filtering
    # everything out). Populate markersInclude from the ranges.
    for (i in 1:nrow(RangesToInclude)) {
      CHROM1 <- RangesToInclude$CHROM[i]
      START <- RangesToInclude$START[i]
      END <- RangesToInclude$END[i]
      posRows <- with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
      if (length(posRows) != 0) {
        markersInclude <- c(markersInclude, markerInfo$ID[posRows])
      }
      cat(length(posRows), " markers in range ", i, " are found in the geno file.\n")
    }
    anyInclude <- TRUE
  }

  if (FALSE) {
    if (!is.null(idstoExcludeFile)) {
      if (anyInclude) {
        stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
      }
      IDsToExclude <- data.table::fread(idstoExcludeFile, header = FALSE, colClasses = c("character"))
      if (ncol(IDsToExclude) != 1) {
        stop("idstoExcludeFile should only include one column.")
      }
      IDsToExclude <- IDsToExclude[, 1]

      posRows <- which(markerInfo$ID %in% IDsToExclude)
      if (length(posRows) != 0) {
        markersExclude <- c(markersExclude, markerInfo$ID[posRows])
      }
      anyExclude <- TRUE
    }

    if (!is.null(rangestoExcludeFile)) {
      if (anyInclude) {
        stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
      }
      RangesToExclude <- data.table::fread(rangestoExcludeFile, header = FALSE, colClasses = c("character", "numeric", "numeric"))
      if (ncol(RangesToExclude) != 3) {
        stop("rangestoExcludeFile should only include three columns.")
      }

      colnames(RangesToExclude) <- c("CHROM", "START", "END")

      for (i in 1:nrow(RangesToExclude)) {
        CHROM1 <- RangesToExclude$CHROM[i]
        START <- RangesToExclude$START[i]
        END <- RangesToExclude$END[i]
        posRows <- with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
        if (length(posRows) != 0) {
          markersExclude <- c(markersExclude, markerInfo$ID[posRows])
        }
      }
      anyExclude <- TRUE
    }
  }

  markersInclude <- unique(markersInclude)
  #  markersExclude = unique(markersExclude)

  # return genotype
  # cat("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data.\n")

  if (anyInclude) {
    if (dosageFileType %in% c("plink", "pgen")) {
      markerInfo <- subset(markerInfo, ID %in% markersInclude)
    }
  }
  #  if(anyExclude)
  #    markerInfo = subset(markerInfo, !ID %in% markersExclude)

  #  anyQueue = anyInclude | anyExclude

  # genoList = list(genoType = genoType, markerInfo = markerInfo, SampleIDs = SampleIDs, AlleleOrder = AlleleOrder, GenoFile = GenoFile, GenoFileIndex = GenoFileIndex, anyQueue = anyQueue)
  # genoList = list(dosageFileType = dosageFileType, markerInfo = markerInfo, anyQueue = anyQueue, genoType = dosageFileType)
  # genoList = list(dosageFileType = dosageFileType, markerInfo = markerInfo, genoType = dosageFileType)
  # if(nrow(markerInfo) != 0){
  #   markerInfo[,POS:=NULL]
  # }

  genoList <- list(markerInfo = markerInfo, genoType = dosageFileType, anyInclude = anyInclude)
  genoList
}


extract_genoIndex_condition <- function(condition, markerInfo, genoType) {
  if (condition != "") {
    condition_original <- unlist(strsplit(condition, ","))
    # if(!is.null(weight_cond)){
    # 	weight_original = unlist(strsplit(weight_cond, ","))
    # }
    conditionDat <- data.frame(SNP = condition_original, condIndex = seq(1, length(condition_original)))
    conditionDat <- setDT(conditionDat)
    # print("head(markerInfo)")
    # print(head(markerInfo))
    # CHROM POS    ID genoIndex       ID2 genoIndex_prev

    markerInfo_conditionDat <- merge(conditionDat, markerInfo, by.x = "SNP", by.y = "ID", sort = FALSE, all.x = TRUE)
    # print(markerInfo_conditionDat)
    # print("markerInfo_conditionDat")
    # SNP condIndex CHROM POS genoIndex  ID2 genoIndex_prev
    if (!is.null(markerInfo$ID2)) {
      # setnames(markerInfo, "genoIndex", "genoIndex2")
      # colnames(markerInfo)[which(colnames(markerInfo) == "genoIndex")] = "genoIndex2"

      markerInfo_conditionDat <- merge(markerInfo_conditionDat, markerInfo, by.x = "SNP", by.y = "ID2", sort = FALSE, all.x = TRUE)
      # print("markerInfo_conditionDat")
      # print(markerInfo_conditionDat)

      # SNP condIndex CHROM.x POS.x genoIndex.x  ID2 genoIndex_prev.x CHROM.y
      # 1: 1:33:A:C         1    <NA>    NA          NA <NA>               NA       1
      # POS.y   ID genoIndex.y genoIndex_prev.y
      # 1:    33 rs33       22895            23399
      setnames(markerInfo_conditionDat, "genoIndex.x", "genoIndex")
      setnames(markerInfo_conditionDat, "genoIndex.y", "genoIndex2")
      posNA <- which(is.na(markerInfo_conditionDat$genoIndex) & !is.na(markerInfo_conditionDat$genoIndex2))
      if (length(posNA) != 0) {
        markerInfo_conditionDat$genoIndex[posNA] <- markerInfo_conditionDat$genoIndex2[posNA]
      }
      markerInfo_conditionDat$genoIndex2 <- NULL
    }
    markerInfo_conditionDat <- markerInfo_conditionDat[which(!is.na(markerInfo_conditionDat$genoIndex)), , drop = FALSE]

    markerInfo_conditionDat <- markerInfo_conditionDat[order(markerInfo_conditionDat$condIndex), ]
    # markerInfo_conditionDat = markerInfo_conditionDat[order(markerInfo_conditionDat$genoIndex), ]

    # markerInfo_conditionDat = markerInfo_conditionDat[which()]
    # posInd = which(markerInfo_conditionDat %in% condition_original)
    # if(length(posInd) == length(condition_original)){
    if (nrow(markerInfo_conditionDat) == length(condition_original)) {
      cond_genoIndex <- markerInfo_conditionDat$genoIndex
      # cond_genoIndex = genoIndex[posInd]
    } else {
      stop(length(condition_original) - nrow(markerInfo_conditionDat), " conditioning markers are not found in the geno file. Please Check.\n")
    }
  } else {
    stop("condition is empty!")
  }
  list(cond_genoIndex = cond_genoIndex)
}


# ---------------------------------------------------------------------------
# Genotype input for steps 0/1 and extractNglmm.
#
# These three entry points differ from step 2 in two ways, which is why they do
# not simply call checkGenoInput()/setGenoInput():
#   - they accept a PREFIX (plinkFile / plink2File) as well as explicit paths;
#   - step 1 permits NO genotype file at all, when a sparse GRM is used to fit
#     the null model and variance-ratio estimation is skipped.
# ---------------------------------------------------------------------------

#' Resolve and validate the genotype input for steps 0/1.
#'
#' @return A list with `genoType` ("plink", "pgen" or ""), and the resolved
#'   `bedFile`/`bimFile`/`famFile`/`pgenFile`/`pvarFile`/`psamFile`.
#' @noRd
checkGenoInputStep1 <- function(plinkFile = "",
                                bedFile = "",
                                bimFile = "",
                                famFile = "",
                                plink2File = "",
                                pgenFile = "",
                                pvarFile = "",
                                psamFile = "",
                                allowNone = FALSE) {
  # A prefix fills in only what was not given explicitly.
  if (plinkFile != "") {
    if (bedFile == "") bedFile <- paste0(plinkFile, ".bed")
    if (bimFile == "") bimFile <- paste0(plinkFile, ".bim")
    if (famFile == "") famFile <- paste0(plinkFile, ".fam")
  }
  if (plink2File != "") {
    if (pgenFile == "") pgenFile <- paste0(plink2File, ".pgen")
    if (pvarFile == "") pvarFile <- paste0(plink2File, ".pvar")
    if (psamFile == "") psamFile <- paste0(plink2File, ".psam")
  }

  anyPlink <- bedFile != "" || bimFile != "" || famFile != ""
  anyPgen <- pgenFile != "" || pvarFile != "" || psamFile != ""

  if (anyPlink && anyPgen) {
    stop(
      "PLINK and PGEN genotype input are mutually exclusive: specify ",
      "plinkFile/bedFile/bimFile/famFile OR plink2File/pgenFile/pvarFile/",
      "psamFile, not both\n"
    )
  }

  if (anyPlink) {
    if (bedFile == "") bedFile <- gsub("bim$", "bed", bimFile)
    if (bimFile == "") bimFile <- gsub("bed$", "bim", bedFile)
    if (famFile == "") famFile <- gsub("bed$", "fam", bedFile)
    Check_File_Exist(bedFile, "bedFile")
    Check_File_Exist(bimFile, "bimFile")
    Check_File_Exist(famFile, "famFile")
    genoType <- "plink"
  } else if (anyPgen) {
    if (pgenFile == "") pgenFile <- gsub("pvar$", "pgen", pvarFile)
    if (pvarFile == "") pvarFile <- gsub("pgen$", "pvar", pgenFile)
    if (psamFile == "") psamFile <- gsub("pgen$", "psam", pgenFile)
    Check_File_Exist(pgenFile, "pgenFile")
    Check_File_Exist(pvarFile, "pvarFile")
    Check_File_Exist(psamFile, "psamFile")
    genoType <- "pgen"
  } else {
    if (!allowNone) {
      stop(
        "Either PLINK (plinkFile/bedFile) or PGEN (plink2File/pgenFile) ",
        "genotype input is required\n"
      )
    }
    genoType <- ""
  }

  list(
    genoType = genoType,
    bedFile = bedFile, bimFile = bimFile, famFile = famFile,
    pgenFile = pgenFile, pvarFile = pvarFile, psamFile = psamFile
  )
}

#' Sample IDs (IIDs) in genotype-file row order.
#'
#' The row order is load-bearing: it is the index space `ptrsubSampleInGeno`
#' and `indicatorGenoSamplesWithPheno` live in on the C++ side.
#' @noRd
readGenoSampleIDs <- function(genoType, famFile = "", psamFile = "") {
  if (genoType == "plink") {
    fam <- data.table::fread(
      famFile,
      header = FALSE, colClasses = list(character = 1:4), data.table = FALSE
    )
    return(as.character(fam[, 2]))
  }
  if (genoType == "pgen") {
    # psam carries a '#'-prefixed column header. Take IID BY NAME: column 1 is
    # FID whenever the header reads "#FID IID ...", and steps 0/1 match
    # phenotypes on IID. '##' lines, if any, are metadata and are not it.
    head <- readLines(psamFile, n = 1000L)
    hdrPos <- grep("^#[^#]", head)
    if (length(hdrPos) == 0) {
      stop("psam file ", psamFile, " has no '#'-prefixed column header\n")
    }
    hdrPos <- hdrPos[1]
    hdr <- sub("^#", "", head[hdrPos])
    hdr <- strsplit(trimws(hdr), "[ \t]+")[[1]]
    iidCol <- which(hdr == "IID")
    if (length(iidCol) != 1) {
      stop("psam file ", psamFile, " has no IID column\n")
    }
    psam <- data.table::fread(
      psamFile,
      header = FALSE, skip = hdrPos, select = iidCol,
      colClasses = "character", data.table = FALSE
    )
    return(as.character(psam[, 1]))
  }
  stop("readGenoSampleIDs: unsupported genoType '", genoType, "'\n")
}

#' Per-variant chromosome codes in genotype-file row order.
#'
#' Positionally aligned with the variant file, which is what makes it safe to
#' index with `getQCdMarkerIndex()` and with the C++ variance-ratio marker
#' indices.
#' @noRd
readVariantChrom <- function(genoType, bimFile = "", pvarFile = "") {
  if (genoType == "plink") {
    bim <- data.table::fread(bimFile, header = FALSE, select = 1, data.table = FALSE)
    return(bim[, 1])
  }
  if (genoType == "pgen") {
    pvar <- data.table::fread(pvarFile, header = TRUE, skip = "#CHROM", data.table = FALSE)
    return(pvar[, 1])
  }
  stop("readVariantChrom: unsupported genoType '", genoType, "'\n")
}

#' Default genotype-input list for the PLINK path.
#'
#' Lets functions that still take bedFile/bimFile/famFile positionally keep
#' working unchanged when no resolved `genoIn` is threaded down to them.
#' @noRd
genoInputPlink <- function(bedFile = "", bimFile = "", famFile = "") {
  list(
    genoType = "plink",
    bedFile = bedFile, bimFile = bimFile, famFile = famFile,
    pgenFile = "", pvarFile = "", psamFile = ""
  )
}

#' Load the GRM genotypes, dispatching on format.
#' @noRd
setgenoFromInput <- function(genoIn, subSampleInGeno, indicatorGenoSamplesWithPheno,
                             memoryChunk, isDiagofKinSetAsOne) {
  if (genoIn$genoType == "pgen") {
    setgeno_pgen(
      genoIn$pgenFile, genoIn$pvarFile, genoIn$psamFile,
      subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne
    )
  } else {
    setgeno(
      genoIn$bedFile, genoIn$bimFile, genoIn$famFile,
      subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk, isDiagofKinSetAsOne
    )
  }
}

#' Per-variant chromosome codes for a resolved genotype input.
#' @noRd
readVariantChromFromInput <- function(genoIn) {
  readVariantChrom(genoIn$genoType, bimFile = genoIn$bimFile, pvarFile = genoIn$pvarFile)
}

#' Sample IDs for a resolved genotype input.
#' @noRd
readGenoSampleIDsFromInput <- function(genoIn) {
  readGenoSampleIDs(genoIn$genoType, famFile = genoIn$famFile, psamFile = genoIn$psamFile)
}
