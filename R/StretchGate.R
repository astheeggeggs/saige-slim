# End-of-run accounting for LEDGER #176's stretch gate (2026-10-02; FIXED
# 2026-10-03).
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
# saddlepoint. A cell at or above the threshold reports the exact two-sided
# tail at the UNSTRETCHED score (src/ExactConv.hpp) in place of the
# saddlepoint's p -- the exact test's model, with the tail computed rather than
# approximated, and like a routed cell with no relatedness correction.
# Settings, each read once per process:
#   SAIGE_STRETCH_GATE=0           evaluate, print and write nothing
#   SAIGE_STRETCH_GATE_APPLY=0     evaluate and report, but keep the saddlepoint
#   SAIGE_STRETCH_GATE_THRESHOLD   the threshold (0.5)
#   SAIGE_STRETCH_GATE_FLOOR       the gate value from which a cell is listed (0.1)
#   SAIGE_EXACT_CONV_SHADOW=1      also compute the exact convolution, with its
#                                  bracket, for EVERY evaluated cell, and list
#                                  every one (slow: ~140 ms a cell on AoU)
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
  thr <- format(rec$threshold)
  cat("Stretch gate: ", bands[1], " weighted binary cell test(s) reached the saddlepoint; ",
      bands[3], " at or above ", thr, " (", bands[2], " at or above 0.1, ",
      bands[4], " at or above 1, ", bands[5], " at or above 3). ",
      if (isTRUE(rec$apply_mode)) {
        paste0(as.numeric(rec$applied_count), " of them report the exact convolution at the ",
               "unstretched score in place of the saddlepoint")
      } else {
        "Observation only (SAIGE_STRETCH_GATE_APPLY=0): no p-value was changed"
      },
      " (LEDGER #176).\n", sep = "")
  if (length(rec$kind) == 0) {
    return(invisible(NULL))
  }
  tab <- data.frame(
    Model = as.integer(rec$model), Kind = rec$kind, Region = rec$region, Group = rec$group,
    max_MAF = rec$maxMAF, Weight = rec$weight,
    gate = signif(rec$gate, 6), z = signif(rec$z, 6), p_spa = signif(rec$p_spa, 6),
    stringsAsFactors = FALSE
  )
  # Always the same columns: the wrapper that delivers this file concatenates
  # one per chromosome under the first file's header, and when these two were
  # written only for a chromosome with a gated cell, a chunk whose first
  # chromosome had none carried a 9-column header over 11-field rows (the
  # splice control and every reader of 2026-10-04 had to work around it).
  has_exact <- rec$applied | rec$shadow
  tab$p_exact <- ifelse(has_exact, signif(rec$p_exact, 6), NA)
  tab$reported <- ifelse(rec$applied, "exact", "saddlepoint")
  if (isTRUE(rec$shadow_mode)) {
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
