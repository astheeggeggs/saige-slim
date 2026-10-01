# End-of-run accounting for LEDGER #172 (2026-09-29): cells of a WEIGHTED set
# (a group-file weight line or a FlexRV transform) whose weight one atom
# carries -- Kish's effective carrier count on the vector the test sees at or
# below --max_MAC_for_ER, with at most 13 material carriers -- take the exact
# test although their allele count is above the cutoff. The saddlepoint on such
# a vector converged to 1e-53..1e-93 on All of Us where ~1e-3 was right (two
# mechanisms on the ledger row: the variance-ratio rescaling crossing the atom,
# and a pinned sample's phantom atom). The exact test enumerates the coin flips
# and ignores the variance ratio, as it already did for tiny unweighted sets.
# Every such cell is counted and named: one line at the end of every region
# run and, when the count is not zero, `<SAIGEOutputFile>.exactByWeight.txt`
# (the AoU splice control reads it as its acceptance list). The records are
# kept in C++ (src/Main.cpp), one per cell, model-independent.

#' Print the end-of-run summary and, if any cell was routed, write the sidecar.
#' @noRd
exactByWeightSummary <- function(OutputFile) {
  rec <- getExactByWeightRecordsInCPP()
  n <- as.integer(rec$count)
  if (n == 0) {
    cat("Exact test by weight: 0 cells (no weighted set whose weight one atom carries; ",
        "LEDGER #172).\n", sep = "")
    return(invisible(NULL))
  }
  tab <- data.frame(
    Kind = rec$kind, Region = rec$region, Group = rec$group, max_MAF = rec$maxMAF,
    Weight = rec$weight, MAC = rec$MAC, n_eff = rec$n_eff, n_material = rec$n_material,
    stringsAsFactors = FALSE
  )
  path <- paste0(OutputFile, ".exactByWeight.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Exact test by weight: ", n, " cell(s) of weighted sets took the exact test because one ",
      "atom carries their weight (effective carriers at or below --max_MAC_for_ER, material ",
      "carriers at most 13; LEDGER #172), listed in '", path, "'.\n", sep = "")
  invisible(tab)
}

# LEDGER #176 (2026-10-01): the cells the routing above WANTED and could not
# take. A cell is concentrated when Kish's effective carrier count on the vector
# the test sees is at or below --max_MAC_for_ER; the exact test can only
# ENUMERATE it when it has at most 13 material carriers. Above that the cell
# falls back to the saddlepoint -- the kernel #172 found converging to a
# confident wrong answer in exactly this regime -- and before this counter
# nothing said so anywhere: the All of Us run carrying #176's unexplained
# BETA/SE contradiction reported zero SPA fallbacks, passed every stage check,
# and had its output accepted by the splice control. It was found by a human
# putting two ancestries side by side, which does not scale to a thousand
# traits.
#
# This is OBSERVATION ONLY. No test, no p-value and no routing decision changes;
# the line and the sidecar are the whole of it. A non-zero count is NOT itself a
# defect -- it names the cells whose p-value rests on the approximation in the
# regime where it is least trusted, which is where to look first.

#' Print the end-of-run summary for concentrated cells above the enumeration
#' cap and, if there are any, write the sidecar.
#' @noRd
exactByWeightAboveCapSummary <- function(OutputFile) {
  rec <- getExactByWeightAboveCapRecordsInCPP()
  n <- as.integer(rec$count)
  if (n == 0) {
    cat("Concentrated cells above the exact test's enumeration cap: 0 ",
        "(LEDGER #176).\n", sep = "")
    return(invisible(NULL))
  }
  tab <- data.frame(
    Kind = rec$kind, Region = rec$region, Group = rec$group, max_MAF = rec$maxMAF,
    Weight = rec$weight, MAC = rec$MAC, n_eff = rec$n_eff, n_material = rec$n_material,
    stringsAsFactors = FALSE
  )
  path <- paste0(OutputFile, ".exactByWeightAboveCap.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Concentrated cells above the exact test's enumeration cap: ", n, " cell(s) had one atom ",
      "carrying their weight (effective carriers at or below --max_MAC_for_ER) but more than 13 ",
      "material carriers, so the exact test could not enumerate them and their p-value is the ",
      "saddlepoint's (LEDGER #176), listed in '", path, "'.\n", sep = "")
  invisible(tab)
}
