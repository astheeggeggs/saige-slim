# Handoff 2026-09-10 (late): the binary null-model fitter, LEDGER #121-#124, and the FlexRV study's state

Read this first, then `docs/state/next.md`, then `python3 harness/tools/ledger.py open`.
The previous handoff, `handoff-2026-09-09-flexrv-e2e.md`, carries the day's
full narrative in dated sections; this file is the consolidated state.

## Commits this session (all local on `strip`, NOTHING pushed)

| commit | row | what |
|---|---|---|
| `5281de8d` | #121 | binary tau update extracted into `tauStepBinary`, old rule pinned |
| `d94146a8` | #121 | the `tol` floor no longer pre-empts the step-halving (GMMAT's rule) |
| `3421a0fa` | #121 | null series behind the owner's drop-to-logistic question |
| `3e775aaa` | #122, #123 | rows written; fix design decided |
| `92ace970` | #124 | AI term's covariate projection applied to the solved vector: OPEN |
| `74c237d3` | #122 | adaptive-damping guard (`R/TauStepGuard.R`); goldens byte-identical |
| `5dc679c3` | #123 | `converged` records the exit taken |

Every commit went through `harness/prepush.sh` green in both compare modes
(last: byte FAIL 0 / PASS 4210, numeric FAIL 0 / PASS 4410, the same 2
warnings and 4 skips as before the session). The working tree holds the
PREVIOUS session's uncommitted FlexRV work (17 modified files, the untracked
flexrv tests, `LEDGER.md` rows 118-120) exactly as found, plus this
session's doc edits (`docs/state/*`, this file). Each ledger commit staged
its own row alone via `git update-index --cacheinfo` on a copy of HEAD's
ledger, so rows 118-120 remain the owner's to commit.

## What the fitter work established (the rows are authoritative)

- **#121 (fixed).** One AI overshoot ended a binary fit at tau2 = 0 because
  the `tol` floor ran before the halving loop. Upstream `main` (365deeff)
  identical. Proved data-not-code both ways.
- **#122 (fixed).** Where the curvature is small the plain Newton step
  overshoots and the iteration cycles; the 2% single-pair criterion
  accepted a point on the cycle (50k trait: 1.31 / 0 / 0.19 by budget,
  fixed point 0.970). The decided bracket design was REFUTED by its first
  drill (score signs from an unconverged inner pass are stale; a bracket
  would have bisected the bundled fit's fourth step). What landed:
  adaptive damping armed on the SECOND direction reversal, growth cap 4x
  when armed, zero as an iterate (score at 0 decides), convergence on the
  raw proposal, never on iteration 1, twice running once a step has been
  damped (keyed on a damped step, not on arming: keyed on arming the
  dense-GRM golden took a 4th iteration). Measured: 30/30 goldens
  byte-identical, traces and iteration counts identical; 50k `y_binary`
  0.968-0.973 from starts 0 / 0.3 / 1.45 / 2.5; the wobbling 50k 1% null
  9 -> 11 iterations, 0.331 -> 0.342, the one change on a non-cycling trait.
- **#123 (fixed).** `converged` from the exit, not `i < maxiter`.
- **#124 (OPEN).** `getAIScore`/`getAIScore_q` project `PAPY_1` instead of
  `APY`; AI off by 0.1-0.8% here, step lengths only. One token per site
  (`APY` for `PAPY_1`, `A0PY` for `PA0PY_1`); the step-1 goldens are the
  gate and are predicted to stay inside 1e-2; characterize by exporting the
  AI term or via the planted-trait reference.
- **Estimator tests (designed, not built)**: SAIGE step 1 vs dense exact
  references on the bundled N = 1,000 fixture agrees up to the seeded
  Hutchinson trace noise (se(tau2) computable exactly from tr((PG)^2);
  19% of a binary tau2 of 0.13 at the default 30 probes; `set_seed(200)`
  makes it repeatable, so a same-sign offset across traits is the trace,
  not a bug). `--isCovariateOffset=TRUE` (CLI default, a speed feature)
  = intercept-only mixed model; the reference must match it. Artefacts
  copied out of the scratchpad to `harness/work/flexrv_study/est121/`
  (planted phenotypes, `pql_ref.R`, before/after golden logs, the step-0/1
  golden script with REPO pinned). Design: previous handoff, "Estimator
  tests for step 1".
- Quantitative fitter: no cycle found (50k weak-signal series converges
  cold and warm to the same values; an "expected" table was an assumption
  and was retracted).

## The FlexRV study (`docs/design/flexrv-sim-study.md`; scripts `harness/sim/flexrv_study/`)

Done: final 50k rebuild gated (age changes under FRAC_CAUSAL by rng
ordering, located, invariants file annotated); scores (informative 0.52,
realistic 0.36, uninformative 0.005); collapse bound 0.28% of carriers;
`grid` comparator added (SKAT-O over high,high:medium x 1e-4/1e-3/1e-2,
Cauchy row = the standard practice); the four null models refitted under
the FINAL install (`step1-guard122.log`): `y_null`, `y_quant` byte-identical
to the pre-#121 files; `y_binary` tau2 = 0.9728, converged, VRs 0.986;
`y_null_bin` tau2 = 0.3416526 (was 0.331; the 1% null the guard damps).
Kept for the inflation measurement: `$D/collapsed/` (tau2 = 0 model),
`$D/pre122/` (tau2 = 1.31 model), `$D/prefix121/` (pre-#121 unaffected
models). `D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k`.

NEXT: `bash harness/sim/flexrv_study/study.sh 3` (144 runs, resumable,
refuses at load > 2.5 per run: check `uptime`, Spotlight's `mds_stores` sat
at 240% CPU for hours today), then `Rscript harness/sim/flexrv_study/analyse.R`
(controls first; an OUTSIDE flag on a null line halts). Then the
collapsed-vs-fixed inflation pair on `y_binary` (point a run at
`$D/collapsed/step1_y_binary.*`). Results go to `harness/bench/NOTES.md`,
the design doc's status, next.md.

## Open questions and items, in priority order

1. The sweep (above). 2. #124's one-token fix with its golden prediction.
3. Build the estimator test (standard tier, ~1 min) from the design.
4. The 5% planted trait's start-0.3 run stops on a single pair at 0.671
   (others 0.626/0.634): today's 2% criterion on a monotone tail, left on
   purpose; a tighter default costs every fit an iteration.
5. Mild excess at alpha 0.05 for SKAT-O on every 1% model incl. the
   permuted control (0.076, CI [0.044, 0.121], n = 210).
6. Upstream: #121/#122/#124 are in saigegit/SAIGE main verbatim; an issue
   with the two-commit patch and the kernel test is a drafting job.

## Traps learned today (each cost real time)

- `optparse` keeps the LAST value of a repeated option (measured on
  `--maxiter` both orders); the memory note said first. An override placed
  before a helper's default is silently lost: the #123 red drill was green
  for that reason. Verify the printed `$option`.
- A scratch drill must check the CLI's exit status; grepping the printed
  tau hid a crash on the boundary exit that the test then caught.
- Never `R CMD INSTALL` while a background drill or the gate is running:
  each Rscript loads the installed package at start, so a sequence of fits
  silently straddles two versions; and the staged install briefly removes
  `SAIGE.so`, which made `prepush.sh` refuse once.
- A one-AI-step "score sign" from a cold start is not the score at the
  fixed point (1.45 vs 0.97 on 50k; 2.0 on the planted trait): the inner
  pass is unconverged after one call. Any bracketing on signs needs it
  converged.
- Fit-time comparisons on this box need an idle box (the goldens' times
  doubled while the gate ran alongside).

Scratchpad drills (`s1/`, `s1fix/`, `est/`, `goldA..E`) are gone on clear;
what matters is under `harness/work/flexrv_study/{nullseries,est121}/` and
in the ledger rows.
