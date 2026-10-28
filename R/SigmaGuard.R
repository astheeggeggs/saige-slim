# End-of-run accounting for LEDGER #183's sparse-GRM variance guard (2026-10-28).
#
# On the sparse-GRM band step 2 solves against Sigma = tau1 I + tau2 K, where
# K is the sparse GRM with every entry below --relatednessCutoff DELETED, and a
# thinned GRM is not positive definite in general. The conjugate-gradient
# solve checks only its residual, so on an indefinite Sigma a vector along a
# negative direction gets a NEGATIVE g' Sigma^-1 g -- which the LEDGER #44
# guard prints as p 1 with BETA NA -- and a vector near a null direction gets
# a tiny one, an inflated Z printed as a valid p. Found on amr Height (the
# random-score control, a pseudo-marker of eleven singletons): the variance
# walked from 10.7 through 0.39 to below zero across the weight ladder with
# the solve converged, p 1.4e-73 one step before the collapse, while the same
# region without the GRM was smooth at p >= 3e-6.
#
# The guard computes the diagonal (no-GRM) variance of the SAME projected
# vector beside the solve's, each under its own band's variance ratio (the
# two numbers as they would enter a reported p; raw they sit ~1/tau1 apart on
# a quantitative trait), and flags a test whose solve is not positive or is
# below the threshold times the diagonal value. By default it OBSERVES: every
# run prints one `Sigma variance:` line and names the flagged tests in
# `<SAIGEOutputFile>.sigmaGuard.txt`; with SAIGE_SIGMA_GUARD_APPLY=1 a flagged
# test takes the diagonal variance under the band's null ratio instead, and
# the sidecar says so. Settings, each read once per process:
#   SAIGE_SIGMA_GUARD=0            evaluate, print and write nothing
#   SAIGE_SIGMA_GUARD_APPLY=1      flagged tests take the diagonal variance
#   SAIGE_SIGMA_GUARD_THRESHOLD    flag below this ratio (0.5)
#   SAIGE_SIGMA_PLANT_RATIO=<r>    test plant: the solve's value times r
# The sidecar names genes and markers; it carries variances, a ratio and the
# reported p, no carrier quantity.

#' Print the end-of-run summary for the sparse-GRM variance guard and, if any
#' test was flagged, write the sidecar.
#' @noRd
sigmaGuardSummary <- function(OutputFile) {
  rec <- getSigmaGuardRecordsInCPP()
  if (!isTRUE(rec$enabled)) {
    cat("Sigma variance: guard off (SAIGE_SIGMA_GUARD=0).\n")
    return(invisible(NULL))
  }
  n_eval <- as.integer(rec$counts[1])
  n_flag <- as.integer(rec$counts[2])
  n_neg <- as.integer(rec$counts[3])
  n_app <- as.integer(rec$counts[4])
  cat("Sigma variance: ", n_eval, " test(s) on the sparse-GRM band; ", n_flag,
      " flagged (solve below ", signif(rec$threshold, 3),
      " x the diagonal variance, or not positive), ", n_neg, " not positive, ",
      n_app, " took the diagonal variance",
      if (isTRUE(rec$apply_mode)) " (SAIGE_SIGMA_GUARD_APPLY=1)" else " (observation only)",
      if (rec$plant_ratio != 1) paste0(" [PLANT: solve x ", rec$plant_ratio, "]"),
      " (LEDGER #183).\n", sep = "")
  if (n_flag == 0) {
    return(invisible(NULL))
  }
  tab <- data.frame(
    Model = as.integer(rec$model), Kind = rec$kind, Region = rec$region, Group = rec$group,
    max_MAF = rec$maxMAF, Weight = rec$weight, MarkerID = rec$marker, MAC = rec$MAC,
    var_grm = signif(rec$var_grm, 6), var_diag = signif(rec$var_diag, 6),
    ratio = signif(rec$ratio, 4), negative = rec$negative,
    reported = ifelse(rec$applied, "diagonal", "grm"), p_reported = rec$p_reported,
    stringsAsFactors = FALSE
  )
  path <- paste0(OutputFile[1], ".sigmaGuard.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  ", nrow(tab), " listed in '", path, "'.\n", sep = "")
  invisible(tab)
}
