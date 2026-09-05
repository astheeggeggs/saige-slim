# End-of-run accounting for SKAT kernel failures (LEDGER #57, the half left
# open by the chunk-file fix).
#
# `get_SKAT_pvalue` wraps `SKAT:::Met_SKAT_Get_Pvalue` in `try(silent = TRUE)`
# and turns any error into NA for Pvalue, Pvalue_Burden and Pvalue_SKAT. That
# used to be the END of the information: the `error.code` it set was never read,
# and a gene could come back as all-NA rows at exit 0 with nothing in the log.
# The chosen policy (2026-09-05, by the user, over warn-and-NA and abort) is to
# KEEP the NA rows -- so every result file is byte-identical to before -- and
# to count and summarise: one line at the end of the run, and a sidecar file
# `<SAIGEOutputFile>.skatFailures.txt` naming every failed cell with the
# kernel's own error message, so that what tends to drive the NAs can be read
# off rather than guessed. The sidecar is written only when something failed.
#
# The record is made in `region_postprocess`, where the cell's labels are in
# scope, from the `error` field `get_SKAT_pvalue` now returns. Under
# `--nThreads > 1` the kernel runs on a fork pool (R/SkatoParallel.R) and its
# answer list, error field included, is served back into the same call site,
# so the accounting happens in the main process either way. A region that
# takes the fastTest ESCALATION pass runs `region_postprocess` twice and a
# persistent failure is recorded twice; the sidecar shows both.

skat_failure_env <- new.env(parent = emptyenv())
skat_failure_env$rows <- list()

#' Forget the failures of a previous run (one dataset per process, but tests
#' call the driver repeatedly).
#' @noRd
skatFailureReset <- function() {
  skat_failure_env$rows <- list()
  invisible(NULL)
}

#' Record one failed SKAT cell.
#' @noRd
skatFailureRecord <- function(region, group, maxMAF, weight, test, error) {
  skat_failure_env$rows[[length(skat_failure_env$rows) + 1L]] <- data.frame(
    Region = as.character(region),
    Group = as.character(group),
    max_MAF = as.character(maxMAF),
    Weight = if (is.null(weight) || length(weight) == 0) NA_character_ else as.character(weight),
    Test = test,
    Error = as.character(error),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

#' The failures recorded so far, one row per failed cell.
#' @noRd
skatFailureTable <- function() {
  if (length(skat_failure_env$rows) == 0) {
    return(data.frame(Region = character(), Group = character(), max_MAF = character(),
                      Weight = character(), Test = character(), Error = character(),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, skat_failure_env$rows)
}

#' Print the end-of-run summary and, if anything failed, write the sidecar.
#' Always prints, so a clean run says "0" rather than nothing -- the negative
#' control has to be visible too.
#' @noRd
skatFailureSummary <- function(OutputFile) {
  tab <- skatFailureTable()
  n <- nrow(tab)
  if (n == 0) {
    cat("SKAT failures: 0 cells returned NA from the SKAT kernel.\n")
    return(invisible(tab))
  }
  path <- paste0(OutputFile, ".skatFailures.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  by_msg <- sort(table(tab$Error), decreasing = TRUE)
  cat("SKAT failures: ", n, " cell(s) in ", length(unique(tab$Region)),
      " region(s) returned NA p-values (Pvalue, Pvalue_Burden, Pvalue_SKAT, or ",
      "their _cond columns for a conditional cell) because the SKAT kernel ",
      "errored. Each is listed with its message in '", path, "'. By message:\n",
      sep = "")
  for (k in seq_along(by_msg)) {
    cat("  ", by_msg[[k]], " x ", names(by_msg)[k], "\n", sep = "")
  }
  invisible(tab)
}
