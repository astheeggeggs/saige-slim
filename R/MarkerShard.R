# Sharded single-variant marker scan (Phase 4, C2).
#
# `--nThreads > 1` used to be a silent no-op in extdata/step2_SPAtests.R: it
# forked only when an --idstoIncludeFile was also given, and otherwise exited 0
# having written nothing (LEDGER #30). The partitioning also belongs here
# rather than in the CLI script, for two reasons:
#
#   * the marker order that must be reproduced is the CHROM/POS-sorted
#     genoIndex vector built inside SPAGMMATtest, not raw bim order; and
#   * the fork has to happen AFTER the null model is loaded and pushed into the
#     C++ singletons, so the children share one copy-on-write model instead of
#     each re-reading it.
#
# Sharding is by INPUT index over that canonical vector, so concatenating the
# shards in block order reproduces the serial file byte for byte: markers
# dropped by QC only make a shard shorter, splitMarker() preserves input order
# within a shard, writeOutfile_single() streams its rows in order, and chunk
# boundaries cannot reach the file (markers_per_chunk only drives a stdout
# progress line). Verified empirically as experiment E0 before this was
# written: golden P7 is byte-identical at markers_per_chunk in
# {1000, 1500, 2000, 3000, 10000}.
#
# fork() via mclapply is the only mechanism available: the C++ state is a set
# of file-scope singletons, so there is one dataset per process.

# Contiguous, near-equal partition of seq_len(M) into min(nShards, M) blocks.
# Deliberately not `M %/% nShards` (which leaves a short trailing block, the
# bug behind the old shell-split path's nThreads+1 output files) and not
# parallel::splitIndices (which returns empty blocks when M < nShards).
splitMarkerShards <- function(M, nShards) {
  if (M < 1) {
    return(list())
  }
  k <- min(as.integer(nShards), as.integer(M))
  sizes <- (M %/% k) + as.integer(seq_len(k) <= (M %% k))
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(seq_len(k), function(b) starts[b]:ends[b])
}

# Append one shard file to an open binary connection, skipping its header line
# unless it is the first shard. Returns the shard's header bytes so the caller
# can check that every shard agrees.
appendShardBytes <- function(shardFile, con, writeHeader, blockSize) {
  inCon <- file(shardFile, open = "rb")
  on.exit(close(inCon), add = TRUE)

  block <- readBin(inCon, what = "raw", n = blockSize)
  newline <- which(block == as.raw(10L))
  while (length(newline) == 0) {
    # A header longer than one block: keep reading until the line ends.
    more <- readBin(inCon, what = "raw", n = blockSize)
    if (length(more) == 0) {
      break
    }
    block <- c(block, more)
    newline <- which(block == as.raw(10L))
  }
  if (length(newline) == 0) {
    stop("Marker shard '", shardFile, "' contains no complete header line.")
  }
  headerEnd <- newline[1]
  header <- block[seq_len(headerEnd)]
  if (writeHeader) {
    writeBin(header, con)
  }
  if (length(block) > headerEnd) {
    writeBin(block[(headerEnd + 1L):length(block)], con)
  }
  repeat {
    block <- readBin(inCon, what = "raw", n = blockSize)
    if (length(block) == 0) {
      break
    }
    writeBin(block, con)
  }
  header
}

# Concatenate shard outputs into a single result file. Raw bytes throughout --
# nothing is parsed or reformatted, so the merged file is byte-identical to
# what a serial run would have written. Header-only shards (every marker
# dropped by QC) are legal; an empty or truncated shard is not.
mergeShardOutputs <- function(shardFiles, OutputFile,
                              blockSize = 64L * 1024L * 1024L) {
  for (b in seq_along(shardFiles)) {
    size <- file.size(shardFiles[b])
    if (is.na(size) || size == 0) {
      stop(
        "Marker shard ", b, " ('", shardFiles[b],
        "') is missing or empty; refusing to merge a truncated result."
      )
    }
  }

  # Merge to a side file and rename, so a failure part-way through cannot
  # leave a plausible-looking partial result at OutputFile.
  mergingFile <- paste0(OutputFile, ".merging")
  con <- file(mergingFile, open = "wb")
  header <- NULL
  tryCatch(
    {
      for (b in seq_along(shardFiles)) {
        thisHeader <- appendShardBytes(shardFiles[b], con,
          writeHeader = is.null(header),
          blockSize = blockSize
        )
        if (is.null(header)) {
          header <- thisHeader
        } else if (!identical(header, thisHeader)) {
          stop(
            "Marker shard ", b, " ('", shardFiles[b],
            "') has a different header from shard 1."
          )
        }
      }
    },
    finally = close(con)
  )

  if (!file.rename(mergingFile, OutputFile)) {
    unlink(mergingFile)
    stop(
      "Could not move the merged marker results into place at '",
      OutputFile, "'."
    )
  }
  invisible(OutputFile)
}

# Do the categorical variance-ratio bins leave some MAC uncovered?
#
# WHY THIS STILL EXISTS, which is not the reason it was written. It was written
# against LEDGER #32: a MAC in no bin left m_varRatioVal holding the PREVIOUS
# marker's value, so a gapped specification made a threaded run differ from a
# serial one (measured on the bundled fixture with a hand-made gap at MAC = 3:
# nThreads in {5, 6, 7, 10, 12, 15, 20, 30, 50} all diverged from nThreads = 1).
#
# **#32 was FIXED in S3** at the landmine: SAIGEClass::selectVarianceRatio now
# THROWS on an uncovered MAC, naming the MAC and printing the bins, and the
# marker path builds a TestContext per marker rather than writing the member --
# so there is nothing left to carry and both arms now stop. The divergence this
# guard was built to prevent can no longer happen.
#
# It is kept because the R-side message arrives EARLIER and reads better than
# the C++ throw, which is what the (inverted) characterization test in
# test-step2-threads.R now pins. That is a user-experience reason, not a
# correctness one -- so if the fork path is retired (see the retirement
# criterion in docs/design/step2-threading-verdict.md), this can go with it, or
# move to wherever the surviving mechanism validates its arguments.
#
# Coverage is (-Inf, minMACExclude[1]] then each (minMACExclude[i],
# maxMACInclude[i]] then (maxMACInclude[k], Inf); a gap is a bin whose lower
# bound starts above everything covered so far.
varianceRatioBinsHaveGap <- function(minMACExclude, maxMACInclude) {
  k <- length(maxMACInclude)
  if (k == 0) {
    return(FALSE)
  }
  if (length(minMACExclude) < k || anyNA(minMACExclude) || anyNA(maxMACInclude)) {
    return(TRUE) # cannot prove coverage; refuse rather than guess
  }
  upper <- minMACExclude[1]
  for (i in seq_len(k)) {
    if (minMACExclude[i] > upper) {
      return(TRUE)
    }
    upper <- max(upper, maxMACInclude[i])
  }
  FALSE
}

# Run the single-variant marker scan, serially or over nThreads forked shards.
# The serial path is exactly the old SAIGE.Marker call.
runMarkerScan <- function(traitType,
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
                          nThreads = 1,
                          nThreadsMarker = 1) {
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # "Serial" here means one PROCESS. nThreadsMarker > 1 makes the C++ marker
  # loop run in-process threads inside it (S6) -- the second, temporarily
  # coexisting parallelism mechanism. The two are not combined: the fork path
  # below passes nThreadsMarker = 1 to its shards.
  runSerial <- function(markerThreads = nThreadsMarker) {
    SAIGE.Marker(
      traitType, genoType, idstoIncludeFile, rangestoIncludeFile,
      genoIndex, CHROM, OutputFile, OutputFileIndex,
      nMarkersEachChunk, isMoreOutput, isImputation, isFirth,
      LOCO, chrom, isCondition, isOverWriteOutput, isAnyInclude,
      nThreadsMarker = markerThreads
    )
  }

  nThreads <- as.integer(nThreads)
  nThreadsMarker <- as.integer(nThreadsMarker)
  if (is.na(nThreadsMarker) || nThreadsMarker < 1L) {
    nThreadsMarker <- 1L
  }
  # The two parallelism mechanisms are not composed. Forking N shards that each
  # spawn M threads is N*M workers with none of the byte-identity reasoning done
  # for it, and it is not what either gate measures. One of these is retiring
  # anyway -- see the criterion in docs/design/step2-threading-verdict.md.
  if (nThreads > 1L && nThreadsMarker > 1L) {
    stop("nThreads (forked shards) and nThreadsMarker (in-process threads) ",
         "cannot both exceed 1; pick one mechanism.")
  }
  if (is.na(nThreads) || nThreads < 1L) {
    stop("nThreads must be a positive integer.")
  }
  if (nThreads == 1L) {
    return(runSerial())
  }

  # SAIGE.Marker applies the LOCO filter itself; applying it here as well
  # (idempotently) keeps the shards balanced over the markers that will
  # actually be tested rather than over the whole genome.
  if (LOCO && genoType %in% c("plink", "pgen")) {
    keep <- which(CHROM == chrom)
    genoIndex <- genoIndex[keep]
    CHROM <- CHROM[keep]
  }
  M <- length(genoIndex)
  if (M < 2) {
    # Nothing to split, and the serial path owns the "no markers" error.
    return(runSerial())
  }

  if (!isGenoReaderForkSafeInCPP()) {
    stop(
      "--nThreads > 1 requires a fork-safe genotype reader. For PLINK, the ",
      ".bed file could not be memory mapped and the buffered fallback is in ",
      "use; for PGEN, the reader is not fork-safe in this version. Forked ",
      "workers would share one file offset and silently decode each other's ",
      "records. Re-run with --nThreads=1."
    )
  }

  # Restart is not supported here: the shards are a partition of the marker
  # list, not of the chunk sequence, so a partially written serial output has
  # no meaning for them.
  outIndex <- checkOutputFile(
    OutputFile, OutputFileIndex, "Marker",
    format(nMarkersEachChunk, scientific = FALSE),
    isOverWriteOutput
  )$indexChunk
  if (outIndex != 1) {
    stop(
      "Restarting an interrupted analysis is not supported with ",
      "--nThreads > 1. Re-run with --nThreads=1 to resume, or delete '",
      OutputFileIndex, "' and re-run with --is_overwrite_output=TRUE."
    )
  }

  blocks <- splitMarkerShards(M, nThreads)
  nShards <- length(blocks)

  # Shards live beside the output file, not in tempdir(): TMPDIR is often a
  # small tmpfs on HPC and these can be large.
  shardDir <- tempfile(
    pattern = ".SAIGE_marker_shards_",
    tmpdir = dirname(OutputFile)
  )
  if (!dir.create(shardDir)) {
    stop("Could not create the shard directory '", shardDir, "'.")
  }
  on.exit(unlink(shardDir, recursive = TRUE), add = TRUE)
  shardFiles <- file.path(shardDir, sprintf("shard_%04d.txt", seq_len(nShards)))

  cat("Marker analysis: ", M, " markers in ", nShards,
    " shards over ", nThreads, " threads\n",
    sep = ""
  )

  # mc.set.seed = FALSE so every worker inherits one fixed parent RNG state
  # rather than a pid/time-derived one; the ER path re-seeds inside C++ at
  # every ComputeExact::Run entry (LEDGER #2), so results do not depend on
  # which worker ran a marker.
  set.seed(1)
  shardResults <- parallel::mclapply(seq_len(nShards), function(b) {
    setSingleOutputPrefixInCPP(shardFiles[b])
    SAIGE.Marker(
      traitType, genoType, idstoIncludeFile, rangestoIncludeFile,
      genoIndex[blocks[[b]]], CHROM[blocks[[b]]], shardFiles[b],
      paste0(shardFiles[b], ".index"), nMarkersEachChunk,
      isMoreOutput, isImputation, isFirth, LOCO, chrom,
      isCondition, TRUE, isAnyInclude,
      nThreadsMarker = 1
    )
  }, mc.cores = nThreads, mc.preschedule = TRUE, mc.set.seed = FALSE)

  # Check every shard before merging anything: mclapply reports a worker that
  # died as a try-error or NULL element, and the old shell-split path threw
  # these away and concatenated whatever files happened to exist.
  for (b in seq_len(nShards)) {
    if (inherits(shardResults[[b]], "try-error")) {
      stop(
        "Marker shard ", b, " failed: ",
        conditionMessage(attr(shardResults[[b]], "condition"))
      )
    }
    if (is.null(shardResults[[b]]) || !is.character(shardResults[[b]])) {
      stop(
        "Marker shard ", b, " returned no result; the worker process ",
        "died. Re-run with --nThreads=1 to see the error."
      )
    }
  }

  mergeShardOutputs(shardFiles, OutputFile)

  # One index record for the whole run. This is deliberately not the serial
  # run's per-chunk index: with nThreads > 1 there is nothing to restart from.
  writeOutputFileIndex(
    OutputFileIndex = OutputFileIndex,
    AnalysisType = "Marker",
    nEachChunk = format(nMarkersEachChunk,
      scientific = FALSE
    ),
    indexChunk = nShards,
    Start = TRUE,
    End = TRUE
  )

  paste0("Analysis done! The results have been saved to '", OutputFile, "'.")
}
