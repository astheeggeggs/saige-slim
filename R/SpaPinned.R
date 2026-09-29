# End-of-run accounting for LEDGER #174 (2026-09-29): PINNED samples and the
# saddlepoint.
#
# A sample is pinned when the null model holds its outcome certain: the
# probability of the outcome it does NOT have (its flip probability, mu for a
# control and 1 - mu for a case) is at or below a threshold. Separation
# produces them: a covariate that isolates a group with no cases drives that
# group's fitted values towards 0 for as long as the fit keeps iterating.
# Such a sample adds nothing to the score or to its variance, but its term in
# the saddlepoint's cumulant generating function, mu (exp(g t) - 1), is inert
# while g t is below -log(mu) and explosive beyond, so where the rest of the
# vector needs a root of that order the finder stops on the phantom and the
# p-value is evaluated at a root that means nothing (LEDGER #172's mechanism B).
#
# The binary saddlepoint therefore leaves pinned samples out of the vector.
# The distributions with and without them differ by at most the sum of the
# flip probabilities left out, so a test does it only while that sum is within
# a budget; over the budget the saddlepoint runs on the full vector as before
# and the test is counted as refused. Settings, each read once per process:
#   SAIGE_SPA_PINNED_EPS      the threshold (default 1e-8, the fit gate's)
#   SAIGE_SPA_PINNED_BUDGET   the budget on the mass one test leaves out (1e-6)
#   SAIGE_SPA_PINNED_GUARD=0  the kill switch
# One line is printed at the end of every step-2 run, and when a test was
# touched `<SAIGEOutputFile>.spaPinned.txt` names it.

#' Record what the saddlepoint R runs itself did with the pinned samples
#' (get_newPhi_scaleFactor_traitType's `$pinned`); a no-op when it did nothing.
#' @noRd
spaPinnedRecord <- function(kind, region, group, maxMAF, weight, model_index, done) {
  if (is.null(done) || (done$n_dropped == 0 && !isTRUE(done$refused))) {
    return(invisible(NULL))
  }
  recordSpaPinnedInCPP(kind, region, group, as.character(maxMAF), as.character(weight),
                       as.integer(model_index), as.integer(done$n_dropped), done$mass,
                       isTRUE(done$refused))
  invisible(NULL)
}

#' Print the end-of-run summary and, if any test was touched, write the sidecar.
#' @noRd
spaPinnedSummary <- function(OutputFile) {
  rec <- getSpaPinnedRecordsInCPP()
  n <- as.integer(rec$count)
  n_refused <- as.integer(rec$refused_count)
  if (rec$threshold == 0) {
    cat("Pinned samples: the guard is OFF (SAIGE_SPA_PINNED_GUARD=0); the saddlepoint saw ",
        "every sample (LEDGER #174).\n", sep = "")
    return(invisible(NULL))
  }
  pinned <- as.integer(rec$model_pinned)
  if (sum(pinned) == 0) {
    cat("Pinned samples: 0 (no fitted value within ", format(rec$threshold),
        " of the observed outcome's boundary; LEDGER #174).\n", sep = "")
    return(invisible(NULL))
  }
  per_model <- paste0(pinned, " of ", as.integer(rec$model_n), collapse = "; ")
  cat("Pinned samples: ", per_model, " sample(s) of the null model(s) have a flip probability at ",
      "or below ", format(rec$threshold), " (total ",
      paste(format(rec$model_mass, digits = 3), collapse = "; "),
      "). They were left out of the saddlepoint in ", n, " test(s)",
      if (n_refused > 0) {
        paste0(" and NOT left out in ", n_refused, " whose dropped mass was over the budget of ",
               format(rec$budget))
      },
      " (LEDGER #174).\n", sep = "")
  if (n + n_refused == 0 || length(rec$kind) == 0) {
    return(invisible(NULL))
  }
  tab <- data.frame(
    Model = as.integer(rec$model), Kind = rec$kind, Region = rec$region, Group = rec$group,
    max_MAF = rec$maxMAF, Weight = rec$weight, MarkerID = rec$marker,
    n_dropped = as.integer(rec$n_dropped), dropped_mass = rec$mass,
    guard = ifelse(rec$refused, "refused", "applied"),
    stringsAsFactors = FALSE
  )
  path <- paste0(OutputFile[1], ".spaPinned.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  ", nrow(tab), " listed with their cell in '", path, "'.\n", sep = "")
  invisible(tab)
}
