# Handoff, 2026-09-09 (evening): FlexRV F2-F4 landed

Previous handoff: `docs/state/handoff-2026-09-09-flexrv-f0.md`. Plan:
`docs/design/flexrv-plan.md` (status header updated); table and decisions:
`docs/design/flexrv-design.md`, whose new "Implementation notes" section is
the detailed account. This file says what a later session needs.

## What landed

Four commits on `strip`, in this order, each gated (fast tier after F2, the
standard tier in both compare modes after F3/F4 via `harness/prepush.sh`):

1. **F2 kernels.** `src/FlexRV.{hpp,cpp}` (namespace `FlexRV`; the test
   export `flexrvWeightsInCPP`), `R/FlexRV.R` (spec parser
   `parseFlexRVTransformSpec`, `flexrvTransformTable`, the R route
   `flexrvWeightsR` for tests only, `parseFlexRVLofAnno`),
   `tests/testthat/test-unit-flexrv.R` (fast tier: two routes to every
   weight at 1e-12 on the boundary grid, the 16 x 12 ordering, a positive
   control).
2. **F3 engine + CLI.** `mainRegionInCPP` gains four arguments (`t_score`,
   `t_isLof`, `t_transformTable`, `t_flexMaxMAF`) and a set kind
   `SET_TRANSFORM` composed after the custom sets; `region_postprocess` gets
   descriptor kind `transform` over the exported `transformWeightMat`;
   `SPAGMMATtest` gets the five `flexRV_*` parameters and the guards;
   `extdata/step2_SPAtests.R` the five flags. Tests
   `test-step2-flexrv.R`, `test-step2-flexrv-guards.R`,
   `helper-step2-flexrv.R`.
3. **F4 fixture + goldens.** `extdata/input/group_new_snpid_flexrv.txt`,
   P58-P61 in `harness/make_goldens.sh` (binary BURDEN 192 sets; the paper's
   quantitative run with the collapse off; SKAT-O on a 2 x 2 spec; the pgen
   mirror), each generated twice into separate trees and `cmp`-identical
   before commit, invariants in `test-invariants-goldens.R`, golden tests in
   `test-step2-goldens.R`.
4. **LEDGER #117 + docs.** The ledger row, `man/` regenerated with roxygen2
   (which also caught up `fitNULLGLMM.Rd`'s stale `LOCO` text -- a real
   drift, not mine), CLAUDE.md's FlexRV section, this handoff, next.md.

## The four findings

- **The fixture's in-sample N is 1,000.** The plan and the F0 handoff said
  10,000 (the fam's row count); step 2 subsets to the model's 1,000 samples.
  Everything about which cells are live changes by a decade; the plan
  paragraph is corrected in place, the arithmetic is in
  `helper-step2-flexrv.R`.
- **LEDGER #117.** A weight of exactly 0 inside a binary SKAT-O cell makes
  `SPA_ER_kernel_related_Phiadj_fast_new` divide by that marker's zero
  variance and the cell returns NA through the SKAT-failure path.
  Pre-existing (reproduced with a zero custom weight on the multiWeights
  fixture), reached systematically by FlexRV's `lof` mask and MAC ladder.
  Recorded, not fixed; characterized in `test-step2-flexrv.R`.
- **`opt$flexRV_score` partial-matched `--flexRV_scoreTransforms`** when the
  flag was absent, turning every plain group run into a FlexRV run. Found by
  the reconciliation tests. The CLI uses `[[ ]]`; memory note written.
- **N = 1,000 cannot measure the cost.** Both the 192-set and the 1-set
  BURDEN run take ~13-15 s wall at ~219 MB peak RSS. The printed memory line
  is a formula (8.8 MB here, 4.6 GB at N = 500k for BURDEN).

## The checkout

Another session's uncommitted edit to `harness/sim/README.md` was present
throughout and was left alone. `flexrv.pdf` (5.7 MB) sits untracked at the
repo root, as before.

## NEXT UP

**F5 on the 50k fixture**, before anything else FlexRV-shaped: wall time and
peak RSS for the 192-set BURDEN run at one cell against a one-set run, on a
quiet box, against a scratch install (the shared-checkout trap). That number
decides whether the two memory levers in the design doc's implementation notes
get pulled. Then the controls: null KS on `p_FlexRV`, the fraction of
individuals carrying >= 2 qualifying ultra-rare variants per gene. The
`a_lat` workaround for the score line is still the plan's.
