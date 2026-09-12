# Identical BURDEN cells are evaluated once (lever 5 of the FlexRV RAM/CPU
# work, 2026-09-12). The C++ emitter (src/Main.cpp, the post-pass of D1) keys
# every tested cell on every input its test reads -- the bytes of its burden
# vector over the region's carriers, the bytes of its weight-independent
# allele count and its set's ER eligibility -- and serves a later cell with
# the same key the earlier cell's result verbatim, per model. The per-run
# counts are printed once at the end of the run, like the SKAT-O memo's
# (R/SkatoParallel.R): always, so a run where the memo never fired says "0"
# rather than nothing. served is a measurement, not a number a test may
# assert (CLAUDE.md, the SKAT-O line).

#' @noRd
burdenMemoReset <- function() {
  resetBurdenMemoInCPP()
  invisible(NULL)
}

#' @noRd
burdenMemoSummary <- function() {
  n <- getBurdenMemoInCPP()
  cat("BURDEN cells: ", n[["asked"]], " asked for, ", n[["evaluated"]], " evaluated, ",
      n[["asked"]] - n[["evaluated"]], " served from an identical cell.\n", sep = "")
  invisible(n)
}
