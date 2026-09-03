# Single-variant marker scan: argument validation and dispatch (Phase 4 C2,
# rewritten by S7).
#
# HISTORY, because the shape of this file is otherwise puzzling. `--nThreads`
# was once a silent no-op in extdata/step2_SPAtests.R -- it forked only when an
# --idstoIncludeFile was also given, and otherwise exited 0 having written
# nothing (LEDGER #30). C2 made it real by forking `mclapply` shards from here;
# S6 added a second, in-process mechanism alongside it; S7 measured the two and
# retired the fork. `--nThreads` now drives the in-process marker threads, and
# `splitMarkerShards`, `appendShardBytes` and `mergeShardOutputs` are gone.
#
# WHY THE FORK LOST, in one line each -- the full evidence is in
# docs/design/step2-threading-verdict.md and harness/bench/NOTES.md:
#   * equal throughput at every ER-eligible fraction f, once Step B let the ER
#     branch run on a worker (before that the in-process path was SLOWER than
#     serial on a rare-variant binary scan, and fork was the only option);
#   * 80-98 MB less resident memory at 4 threads, growing with N and thread
#     count, because there is one model rather than N copy-on-write copies;
#   * restart works. The fork path refused it outright.
#
# TWO FINDINGS FROM THE FORK ERA THAT MUST NOT BE LOST WITH IT. Both are about
# THIS codebase, not about fork, and both are recorded in CLAUDE.md:
#   * E0 -- chunk boundaries cannot reach the output bytes. `markers_per_chunk`
#     only drives a stdout progress counter; golden P7 is byte-identical at
#     every chunk size tried, from one thousand markers up to ten thousand.
#   * E1 -- a forked child had to close `g_out.single` explicitly, because
#     `mcexit` calls `_exit()`, no static destructor runs, and the child
#     silently lost everything past the stream's 192 KiB buffer boundary WHILE
#     STILL RETURNING SUCCESS. `closeOutfile_single()` in SAIGE.Marker is
#     byte-neutral without the fork and is kept: it makes the stream's lifetime
#     explicit, and the same buffer is the mechanism behind LEDGER #47.

# Do the categorical variance-ratio bins leave some MAC uncovered?
#
# WHY THIS STILL EXISTS, which is not the reason it was written. It was written
# against LEDGER #32: a MAC in no bin left `m_varRatioVal` holding the PREVIOUS
# marker's value, so a gapped specification made a parallel run differ from a
# serial one.
#
# **#32 was FIXED in S3** at the landmine: `SAIGEClass::selectVarianceRatio` now
# THROWS on an uncovered MAC, naming the MAC and printing the bins, and the
# marker path builds a `TestContext` per marker rather than writing the member.
# There is nothing left to carry, and both a serial and a threaded run now stop
# at the same marker. **The divergence this guard was built to prevent can no
# longer happen, and the guard is retained purely because its message is better**
# -- it fires before a single marker is read and names the offending flags,
# where the C++ throw names a MAC. That is a user-experience reason, and the
# message below no longer claims otherwise.
#
# It stays gated on a parallel scan rather than firing always, because making it
# fire at nThreads == 1 would change where a serial run fails and what it says,
# for no gain.
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

# Run the single-variant marker scan. One process either way: nThreads > 1 makes
# the C++ marker loop run that many in-process threads, and at 1 it runs a
# literally serial loop that touches no threading machinery at all.
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
                          nThreads = 1) {
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }
  nThreads <- as.integer(nThreads)
  if (is.na(nThreads) || nThreads < 1L) {
    stop("nThreads must be a positive integer.")
  }

  # A reader that cannot clone cannot be threaded. The C++ loop degrades to one
  # thread on its own and prints why, which is not the silent no-op of LEDGER
  # #30 -- but a line of stdout is easy to lose in a cluster log, so ask first
  # and raise a real R warning too. Degrade rather than error, deliberately: it
  # is the same posture as a group test, and erroring would break anyone
  # scripting one --nThreads across a mix of PLINK and PGEN inputs.
  if (nThreads > 1L && !isGenoReaderCloneableInCPP()) {
    warning(
      "--nThreads=", nThreads, " was requested but this genotype reader cannot ",
      "be cloned, so the marker scan will run on ONE thread. For PGEN the ",
      "reader is single-threaded in this version; for PLINK the .bed could not ",
      "be memory mapped and the buffered fallback is in use.",
      call. = FALSE
    )
  }

  SAIGE.Marker(
    traitType, genoType, idstoIncludeFile, rangestoIncludeFile,
    genoIndex, CHROM, OutputFile, OutputFileIndex,
    nMarkersEachChunk, isMoreOutput, isImputation, isFirth,
    LOCO, chrom, isCondition, isOverWriteOutput, isAnyInclude,
    nThreadsMarker = nThreads
  )
}
