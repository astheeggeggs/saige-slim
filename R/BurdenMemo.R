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
#
# THE UNIT IS CELLS, and deliberately. The dedupe key is model-independent (a
# cell's burden column, its weight-independent count, its set's ER
# eligibility), so ONE decision serves all P models and the natural unit for
# the decision is the cell; what it saves is P tests, so on a multi-model run
# the tests avoided are P times the number this line prints. The SKAT-O line
# accumulates per region_postprocess drive, which IS per model, so above one
# model the two lines do not share a unit. Reviewed 2026-09-12 and left as is:
# each line is self-consistent, this one is the more conservative of the two,
# and re-basing it on tests would move a printed count that a test asserts
# (test-step2-flexrv.R pins asked == 224 = 2 genes x 112 sets at P = 1, since
# LEDGER #142 drops the five MAC rungs at or above 2*N*maxMAF at N = 1,000).

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
