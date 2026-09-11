# LEDGER #127: provenance for the pooled per-region row, the one whose Group
# column reads Cauchy.
#
# The row's ARITHMETIC is not in question and is not touched. Its Bonferroni
# fallback -- when some cell p-value is exactly 1, CCT returns
# min(1, n * min(p)) -- is SAIGE-GENE+'s own prescription (Zhou et al.,
# Nat Genet 2022, PMC9534766; see LEDGER #86), and the Cauchy combination is
# the paper's method. What was missing is DISCLOSURE: three different things
# all produced a well-formed pooled p-value at exit 0 that meant something
# other than it said --
#   (a) a cell whose kernel failed (LEDGER #57) simply vanishes from the
#       combination, so with two weight sets and one failure the pooled row is
#       the surviving set's row verbatim (LEDGER #117);
#   (b) when the min-group-MAC gate refuses every cell of a BURDEN region the
#       row still read Pvalue_Burden = 0.5, which is the CCT of no inputs;
#   (c) the per-set MAC of LEDGER #118 refuses more cells than the per-cell one
#       did, shrinking the denominator the same silent way.
# The fix is one mechanism for all three: record what was asked, tested and
# refused and WHICH CCT rule produced the number, and refuse to print a pooled
# p-value when nothing was tested ((b) is now NA, handled in src/Main.cpp for
# the burden path and by get_CCT_pvalue's own NA for the SKAT-O path).
#
# The counts are a SIDECAR, `<SAIGEOutputFile>.pooledTests.txt`, in the shape
# of `.skatFailures.txt` (LEDGER #57) and for the same reason: no result file
# moves a byte. It is written only when there is something to disclose -- a
# refused cell, an NA that dropped out of the combination, a rule other than
# the plain Cauchy one, or a region in which nothing was tested -- so a clean
# run produces no file, and the one-line summary always prints so that the
# negative control is visible too.
#
# Two sources, one table: the SKAT-O path records here, from
# region_postprocess where the cells are tested; the BURDEN path records in
# C++, where its pooled row is written, and R reads those with
# getPooledRecordsInCPP(). A region that takes the fastTest ESCALATION pass is
# post-processed twice and appears twice, exactly as in `.skatFailures.txt`.

pooled_row_env <- new.env(parent = emptyenv())
pooled_row_env$rows <- list()

#' Forget the previous run's records (one dataset per process, but the tests
#' drive the region scan repeatedly).
#' @noRd
pooledReset <- function() {
  pooled_row_env$rows <- list()
  resetPooledRecordsInCPP()
  invisible(NULL)
}

#' Which rule CCT() will take on this vector of p-values. The predicates and
#' their order are CCT's own (R/CCT_modified.R; src/CCT.cpp is the same test in
#' the same order), so this names the branch rather than guessing at it.
#' @noRd
pooledRule <- function(pvalue) {
  p <- pvalue[!is.na(pvalue)]
  if (length(p) == 0) {
    return("none")
  }
  if (any(p == 0)) {
    return("zero")
  }
  if (any(p == 1)) {
    return("min-p Bonferroni")
  }
  "Cauchy"
}

#' Record one region's pooled row (SKAT-O path).
#' @noRd
pooledRecord <- function(region, model, rows, refused, combined, rule, pvalue) {
  pooled_row_env$rows[[length(pooled_row_env$rows) + 1L]] <- data.frame(
    Region = as.character(region),
    Model = as.integer(model),
    Path = "SKAT-O",
    Rows = as.integer(rows),
    Pvalues_combined = as.integer(combined),
    Cells_refused_minGroupMAC = as.integer(refused),
    Rule = as.character(rule),
    Pvalue_pooled = as.numeric(pvalue),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

#' Every record of this run: the SKAT-O rows made here plus the BURDEN rows
#' C++ kept, in that order.
#' @noRd
pooledTable <- function() {
  empty <- data.frame(
    Region = character(), Model = integer(), Path = character(),
    Rows = integer(), Pvalues_combined = integer(),
    Cells_refused_minGroupMAC = integer(),
    Rule = character(), Pvalue_pooled = numeric(),
    stringsAsFactors = FALSE
  )
  skato <- if (length(pooled_row_env$rows) == 0) empty else do.call(rbind, pooled_row_env$rows)
  cpp <- as.data.frame(getPooledRecordsInCPP(), stringsAsFactors = FALSE)
  if (nrow(cpp) > 0) {
    cpp$Path <- "BURDEN"
    cpp <- cpp[, colnames(empty), drop = FALSE]
  } else {
    cpp <- empty
  }
  rbind(skato, cpp)
}

#' Print the end-of-run summary and, when a region has something to disclose,
#' write the sidecar. A region is worth disclosing when the gate refused one of
#' its cells, when a p-value NA'd out of the combination, when the rule was not
#' the plain Cauchy one, or when nothing was tested at all.
#' @noRd
pooledSummary <- function(OutputFile) {
  tab <- pooledTable()
  if (nrow(tab) == 0) {
    return(invisible(tab))
  }
  interesting <- tab$Cells_refused_minGroupMAC > 0 |
    tab$Pvalues_combined < tab$Rows |
    tab$Rows == 0 |
    tab$Rule != "Cauchy"
  n <- sum(interesting)
  if (n == 0) {
    cat("Pooled (Cauchy) rows: ", nrow(tab),
        " region row(s), every one a plain Cauchy combination of every cell tested.\n",
        sep = "")
    return(invisible(tab))
  }
  path <- paste0(OutputFile, ".pooledTests.txt")
  out <- tab[interesting, , drop = FALSE]
  utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
  nRefused <- sum(out$Cells_refused_minGroupMAC)
  nDropped <- sum(pmax(out$Rows - out$Pvalues_combined, 0))
  nNone <- sum(out$Rule == "none")
  cat("Pooled (Cauchy) rows: ", nrow(tab), " region row(s), ", n,
      " of which do not combine every cell of the region by the plain Cauchy ",
      "rule (", nRefused, " cell(s) refused by --minGroupMAC_in_BurdenTest, ",
      nDropped, " p-value(s) absent from a combination, ", nNone,
      " region(s) with nothing tested, whose pooled p-value is NA). Each is ",
      "listed in '", path, "'. By rule:\n", sep = "")
  # The per-rule tally deliberately does NOT use skatFailureSummary's
  # "  <n> x <message>" shape: `test-step2-skat-failures.R` scrapes the log with
  # `grep("^  [0-9]+ x ")` and would add these lines to ITS count, which is how
  # this was found (8 NA cells, a tally that summed to 10).
  by_rule <- sort(table(out$Rule), decreasing = TRUE)
  for (k in seq_along(by_rule)) {
    cat("  rule ", names(by_rule)[k], ": ", by_rule[[k]], " region row(s)\n", sep = "")
  }
  invisible(tab)
}
