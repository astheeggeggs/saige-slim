# Step-1/step-2 eligibility gates, and what the AoU v2 refit found (2026-09-23)

Two commits landed and one LEDGER row opened. The row is the important part:
**a strip-branch regression that makes 21 real binary fits fail where upstream
1.4.2 succeeds**, found because the new gates refused them instead of saving
them silently.

## What landed

`db5bb558` -- `R/FitGates.R` (456 lines) + `tests/testthat/test-step1-fit-gates.R`
(434 lines), five options threaded through `fitNULLGLMM` and `SPAGMMATtest`:

| option | step 1 | step 2 | default |
|---|---|---|---|
| `minCaseCount` | yes | yes | 100 (warns below 500) |
| `minCovariateCellCount` | yes | -- | 10 |
| `separationCheck` + `separationMaxFraction` | covariate GLM | post-fit `mu` | `error`, 0.02 |
| `convergenceCheck` | refuses to SAVE | refuses to TEST | `error` |
| `minSampleCount` | yes | -- | 0 (off) |

`3e61864a` -- the `--separationMaxFraction` help rendered a literal `%%` in
both `man/` and the CLI (Rd wants `\%`, optparse wants `%`; reworded to avoid
the character). Found by reading the BUILT IMAGE's `--help`, not the source.

Gate on both: `FAIL 0 | WARN 3 | SKIP 4 | PASS 6009`, both compare modes, no
golden moved.

## Things that are measured, not assumed

- **`SAIGE_FIT_GATES=0` in `harness/env.sh`** disables all of them and PRINTS a
  warning naming them. That one line replaced an opt-out flag on ~150 command
  lines (16 step-1 golden invocations, 69 step-2, ~38 test files, ~17 harness
  scripts). The bundled fixture fails two gates BY CONSTRUCTION -- 98 cases,
  and `a10` level 2 has 7 -- which is correct: it is a mechanical fixture, not
  a valid analysis.
- **A per-dummy covariate scan is BLIND to a zero-case REFERENCE level**,
  because `model.matrix` gives the reference level no column. Measured on a
  3-level factor whose reference had 0 cases: dummy scan FALSE, factor scan
  TRUE. `smallCovariateCells` therefore reconstructs the factor, searching
  `covarColList` and NOT `names(data)` (the frame also holds the phenotype,
  IID and IndexGeno, and nothing outside the model's covariates may be
  counted).
- **Separation is a FRACTION, not any pinned sample.** The rule this replaced
  refused a real AoU model over fewer than 20 samples in 2,737. Owner set 0.02.
- **`converged = FALSE` covers TWO exits** the flag cannot distinguish: running
  out at `--maxiter`, and the bail-out at `max(tau) > tol^(-2)`. `--tol` sets
  both, so 0.02 aborts above tau 2500 and 0.2 above tau 25.

## LEDGER #161 -- the regression this exposed

`R/SAIGE_fitGLMM_fast.R:207`, the boundary exit added by #122's fix:

```
strip:    if (tau0[2] == 0 && score <= 0) { didConverge <- TRUE; break }
upstream: if (tau[2] == 0) break            # unconditional, AFTER the update
```

Upstream's update is `max(0, tau0 + tau0^2 * score/n)`, for which **zero is
absorbing**; the strip branch's is a Newton step `raw = tau0 + score/AI`
(`src/SAIGE_fitGLMM_fast.cpp:3758`), which is not. So when the score is
POSITIVE at tau2 = 0 the exit does not fire, the guarded step lands back on 0,
`tau` never moves, and the convergence test at `:253` -- judged on the UNDAMPED
proposal -- never fires either. The loop runs the full `--maxiter` with
`tau = (1, 0)` and reports not-converged.

**Controlled against upstream** `gcr.io/.../saige:1.4.2` with identical inputs,
covariates, GRM and `--tol 0.02`: upstream ran ONE iteration (`Iteration 1 1 0`,
`tau: 1 0`, `Final 1 0`); the strip branch ran 20. Not a duplicate of #122 --
that is the iteration CYCLING; this is non-termination at the boundary,
introduced by #122's own fix. In fairness the #122 condition is the better
statistics (upstream calls zero optimal even when the score says increase);
what is missing is a STALL case.

**Not yet measured: the sign of `AI` at the boundary.** `AI < 0` would explain
a negative Newton step with a positive score, hence the clamp. That is the next
datum.

## What the AoU v2 refit found

166 jobs (43 traits x 5 ancestries, after a per-ancestry eligibility screen),
image `1.5.2-dev-3e61864a`, `--tol 0.02` for BOTH trait types.

- **The quantitative `tol` fix is real.** `HbA1c` EUR converged in 3 iterations
  where the BRaVa wrapper's `tol=1e-5` ran it out at 20. That wrapper used
  `0.2` for binary (mislabelled "SAIGE DEFAULT") and `1e-5` for quantitative;
  upstream and the All-by-All pipeline both use 0.02 for both, and All-by-All
  never sets `--tol` at all (`github.com/atgu/aou_gwas`, public).
- **134 ok, 30 refused, 2 no-status, 49 never submitted.** Of the 30: 9 are
  legitimate case-count refusals, 21 are #161.
- **The pre-fit screen overstates case counts by 22-45%**, because it counts in
  the phenotype file while the gate counts in the FITTED set (after
  intersecting genotypes, GRM and complete cases over 26 covariates). e.g. afr
  AMD 160 -> 97, amr CervCanc 142 -> 78. `MIN_CASES = 100` on the file needs to
  be ~170, or the screen should use the intersected set.
- **Cost/size facts**: a fit is 2-4 minutes, peak RSS 703 MB (quant, N 38.5k) to
  1,617 MB (binary, N 170.6k) -- far below the 60 G All-by-All requests. `.rda`
  is 45-140 MB, so ~166 models is 15-20 GB: screen them with a dsub job, not by
  downloading. A sparse-GRM fit is SINGLE-THREADED whatever `--nThreads` says
  (`R/SAIGE_fitGLMM_fast.R:918`), so size the machine for RAM, not cores.

## Open

1. **#161** -- fix the stall detector (terminate when the guarded tau equals
   tau0, reporting boundary-convergence or a distinct stalled status), or
   accept upstream's no-random-effect model for those traits. Blast radius is
   every binary sparse-GRM fit whose tau2 sits at the boundary; check the
   goldens. All-by-All's published answer to their own non-convergence was to
   EXCLUDE the phenotypes.
2. **The screen's `MIN_CASES`** calibration, above.
3. **A trap I set and removed the same day**: `checkConvergence`'s message
   quoted SAIGE's own "Large variance estimate observed in the iterations"
   warning verbatim, so a log classifier grepping for that phrase matched the
   GATE'S message and reported a bail-out for all 21 refusals. Fixed; a comment
   at the site says why. Do not re-introduce it.
