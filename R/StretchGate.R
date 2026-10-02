# End-of-run accounting for LEDGER #176's stretch gate (2026-10-02).
# OBSERVATION ONLY: nothing here, and nothing it reads, reaches a result file.
#
# The binary saddlepoint is evaluated at q = Tstat / sqrt(varianceRatio) + m1
# against a cumulant generating function built on the unscaled vector, so the
# variance ratio STRETCHES the observed score. Where everything but the
# heaviest carrier's term varies by much more than the stretch, that is the
# intended relatedness correction. Where it does not, the saddlepoint is asked
# for a tail the cell can only reach by a freak excursion of the rest, and it
# converges to a confident wrong answer. The gate is the ratio
#   stretch / sd(everything but the heaviest carrier's term)
# evaluated for every weighted set's binary cell whose test reaches the
# saddlepoint. Settings, each read once per process:
#   SAIGE_STRETCH_GATE=0         evaluate, print and write nothing
#   SAIGE_STRETCH_GATE_FLOOR     the gate value from which a cell is listed (0.1)
#   SAIGE_EXACT_CONV_SHADOW=1    also compute the exact convolution at the
#                                UNSTRETCHED score for every evaluated cell and
#                                list every one, with its bracket and its cost
# The sidecar names genes, so it carries no carrier count and no per-carrier
# quantity: only the gate, the score's z, and p-values.

#' Print the end-of-run summary for the stretch gate and, if any cell is
#' listed, write the sidecar.
#' @noRd
stretchGateSummary <- function(OutputFile) {
  rec <- getStretchGateRecordsInCPP()
  if (!isTRUE(rec$enabled)) {
    cat("Stretch gate: OFF (SAIGE_STRETCH_GATE=0; LEDGER #176).\n")
    return(invisible(NULL))
  }
  bands <- as.numeric(rec$bands)
  cat("Stretch gate: ", bands[1], " weighted binary cell test(s) reached the saddlepoint; ",
      bands[3], " at or above ", format(rec$threshold), " (", bands[2], " at or above 0.1, ",
      bands[4], " at or above 1, ", bands[5], " at or above 3). Observation only: no p-value ",
      "was changed (LEDGER #176).\n", sep = "")
  if (length(rec$kind) == 0) {
    return(invisible(NULL))
  }
  tab <- data.frame(
    Model = as.integer(rec$model), Kind = rec$kind, Region = rec$region, Group = rec$group,
    max_MAF = rec$maxMAF, Weight = rec$weight,
    gate = signif(rec$gate, 6), z = signif(rec$z, 6), p_spa = signif(rec$p_spa, 6),
    stringsAsFactors = FALSE
  )
  if (isTRUE(rec$shadow_mode)) {
    tab$p_exact <- signif(rec$p_exact, 6)
    tab$p_exact_lo <- signif(rec$lo, 6)
    tab$p_exact_hi <- signif(rec$hi, 6)
    tab$ms <- signif(rec$ms, 3)
  }
  path <- paste0(OutputFile[1], ".stretchGate.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  ", nrow(tab), " listed in '", path, "'",
      if (isTRUE(rec$shadow_mode)) {
        paste0(" with the exact convolution beside the saddlepoint (shadow mode, ",
               signif(sum(rec$ms) / 1000, 3), " s in total)")
      },
      ".\n", sep = "")
  invisible(tab)
}
