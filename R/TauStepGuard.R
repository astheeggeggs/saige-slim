# LEDGER #122: the guard around the binary AI-REML step on tau2.
#
# The AI step is a Newton step, score / AI, on the one-dimensional REML score
# of tau2 (tau1 is fixed at 1 for a binary trait). Where the curvature is
# small, low prevalence with few relatives, the step overshoots and the plain
# iteration cycles, so the reported tau2 is a function of the iteration budget
# (the 50k FlexRV trait: 1.31, 0 or 0.19 for a fixed point of 0.970).
#
# The guard damps the step: once the iteration has reversed direction twice,
# every further reversal halves the fraction of the step that is applied, and
# the step may then grow at most `growCap`-fold from a positive tau2. It is
# armed by the SECOND reversal on purpose: a monotone approach that overshoots
# once near the end and comes back (every bundled binary fit, one reversal)
# takes exactly the path it took before this file existed, so the guard is
# inert on a healthy fit by construction. Damping was chosen over a bracket
# on the score's signs because the score at a tau2 depends on the inner
# (alpha, b) state, which Get_Coef converges only to tol.coef: after a large
# jump the sign is stale, and a bracket built from stale signs excludes the
# root (measured on the bundled fit's own trace and on the 50k trait, see the
# ledger row). Damping needs no sign to be right.
#
# Pure functions over a small state list so the rules are unit-tested on
# synthetic score functions (tests/testthat/test-unit-tau-guard.R).

tauGuardInit <- function() {
  list(prevDir = 0, reversals = 0L, damp = 1, armed = FALSE, fired = FALSE)
}

# Apply the guard to the proposal `step` made from `tau0`, given the raw
# (unfloored, unhalved) Newton proposal `raw`, whose direction is what counts
# as a reversal. Returns the proposal to apply and the updated state.
tauGuardApply <- function(state, tau0, raw, step, growCap = 4, armAt = 2L) {
  dir <- sign(raw - tau0)
  if (dir != 0) {
    if (state$prevDir != 0 && dir != state$prevDir) {
      state$reversals <- state$reversals + 1L
      if (state$reversals >= armAt) {
        state$armed <- TRUE
        if (state$reversals > armAt) state$damp <- state$damp / 2
      }
    }
    state$prevDir <- dir
  }
  p <- step
  if (state$armed) {
    p <- tau0 + state$damp * (step - tau0)
    if (tau0 > 0 && p > growCap * tau0) p <- growCap * tau0
    if (p != step) state$fired <- TRUE
  }
  list(p = p, state = state)
}
