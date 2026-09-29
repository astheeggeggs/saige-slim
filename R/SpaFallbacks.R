# End-of-run accounting for saddlepoint failures (LEDGER #164).
#
# On a binary trait every marker or cell whose |Z| exceeds --SPAcutoff goes
# through the saddlepoint approximation, and the saddlepoint can fail: a root
# the Newton iteration does not reach (LEDGER #163 before its fix), a mirrored
# q outside the attainable range, a NaN step, the iteration cap, or a p that
# underflows. Up to 2026-09-26 the fallback was the normal-approximation
# score p, written to p.value / Pvalue_Burden like any other p with only the
# single-variant sidecar's Is.SPA column to tell, which nothing downstream
# reads. The policy chosen by the owner (2026-09-26) is BOTH halves: the
# fallback is now the exact resampling test at any MAC and for any weight
# set (getMarkerPval, src/SAIGE_test.cpp), and every fallback is counted and
# named -- one line at the end of every step-2 run, and a sidecar
# `<SAIGEOutputFile>.spaFallbacks.txt` with the cell, the normal p that was
# NOT reported and the exact p that was. The sidecar is written only when the
# count is not zero; the line is always printed, so a clean run says "0".
#
# The records are kept in C++ (src/Main.cpp, beside the pooled records): the
# test sets a per-worker flag and the four call sites that know their cell's
# labels record it. The count can exceed the recorded rows, because a
# conditioning-marker test counts but has no cell to name; the line says so
# when that happens. A BURDEN twin served from the memo (R/BurdenMemo.R)
# inherits its primary's p and is not recorded again.
#
# LEDGER #173 (decided 2026-09-29) rides on the same summary: on the marker
# path under --is_noadjCov=TRUE the fast kernel's Gaussian remainder is the
# centred genotype's variance minus the projected carriers' part, which is not
# a variance and could in principle go negative. If it does the dense kernel
# runs instead and the run says how many times.

#' Print the end-of-run summary and, if anything fell back, write the sidecar.
#' @noRd
spaFallbackSummary <- function(OutputFile) {
  # LEDGER #173's safety net: said only when it happened (it never has)
  n_neg <- as.integer(getSpaNegRemainderCountInCPP())
  if (n_neg > 0) {
    cat("SPA remainder: ", n_neg, " test(s) had a fast-kernel Gaussian remainder that was not a ",
        "variance and ran the dense kernel instead (LEDGER #173).\n", sep = "")
  }
  rec <- getSpaFallbackRecordsInCPP()
  n <- as.integer(rec$count)
  if (n == 0) {
    cat("SPA fallbacks: 0 tests fell back from the saddlepoint to the exact test.\n")
    return(invisible(NULL))
  }
  tab <- data.frame(
    Kind = rec$kind, Region = rec$region, Group = rec$group, max_MAF = rec$maxMAF,
    Weight = rec$weight, MarkerID = rec$marker, MAC = rec$MAC,
    p_normal_not_reported = rec$pvalNA, p_exact_reported = rec$pval, Note = rec$note,
    stringsAsFactors = FALSE
  )
  path <- paste0(OutputFile, ".spaFallbacks.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  by_kind <- sort(table(tab$Kind), decreasing = TRUE)
  cat("SPA fallbacks: ", n, " test(s) fell back from the saddlepoint to the exact ",
      "test (LEDGER #164); the normal-approximation p was NOT reported for them. ",
      nrow(tab), " listed with their cell in '", path, "'",
      if (nrow(tab) < n) " (the rest were conditioning-marker tests, which have no cell)",
      ". By kind:\n", sep = "")
  for (k in seq_along(by_kind)) {
    cat("  ", by_kind[[k]], " x ", names(by_kind)[k], "\n", sep = "")
  }
  n_na <- sum(tab$p_exact_reported == "NA")
  if (n_na > 0) {
    cat("  of which ", n_na, " could not run the exact test either and carry p.value NA ",
        "(the Note column says why the saddlepoint failed, and why the exact test could not run).\n", sep = "")
  }
  invisible(tab)
}
