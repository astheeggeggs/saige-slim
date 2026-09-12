# Handoff 2026-09-12 (evening): the review fixes, the CCT answer, and the gate

Everything here is PUSHED (`slim/strip` at `b16fa83c`) and gated. Read this,
then `docs/state/next.md`, then `python3 harness/tools/ledger.py open`.

## 1. Three bugs fixed in the FlexRV lever work, one of which was live

An adversarial review of the then-unpushed lever diff (`48ea91a8..2f1c1960`,
682 changed lines in `src/Main.cpp`) reported 12 findings; 10 survived a 3-lens
adversarial pass and collapsed to 5 distinct defects.

- **`e72e4186` -- a data race this branch had already fixed once.** Lever 2 put
  two cell prepares on the thread pool while they still shared the frame-scope
  `hasVarRatio`, which `selectVarianceRatio` READS BACK at three guards and
  throws on. **Drilled with TSAN, the real bug as its own plant**: `d36e7277`
  reports 5 races at t=4 and 3 at t=8 on `src/SAIGE_test.cpp:1274`/`:1279` (the
  flag's own two writes); the fixed tree reports none. `bytes-vs-serial` was
  `identical` in BOTH trees, so no golden ever moved -- the byte gates were
  right and blind at the same time.
- **`8f4d5db7` -- `ScratchGuard`.** Both cell-test lambdas cleared the worker
  scratch on the line AFTER `getMarkerPval`, and the ER arm's
  `NeedsMainThread` deferral is a routine throw that the pool catches and
  continues past with the same scratch. Now cleared from a destructor. Only
  the TEST lambdas need it: the prepares' only throw is fatal.
- **`29f6f07c` -- `CarrierStore::begin()`, and this one produced WRONG
  NUMBERS.** It reset the slot map on a changed sample count but not the value
  buffers. Two `SPAGMMATtest` calls in one process (1,000 samples then 950):
  both return normally, and the second reports `GENE1 lof` at max_MAF 0.001
  with MAC 7 / p 7.860290E-02 where the truth is MAC 2 / p 3.182528E-01.
  Pinned by `tests/testthat/test-step2-region-tworun.R`, red on a pre-fix
  scratch install at exactly the second-run assertion.
- **`690ae4f9`** declares the region worker's TSAN fork edge, so the clean arm
  reads `races=0`; it had reported 1-2 main-vs-worker warnings in EVERY build,
  which meant a non-zero count there was evidence of nothing.
- **Two findings deliberately NOT changed**, reasons at the sites: the BURDEN
  memo line's unit is cells (`R/BurdenMemo.R`), and the three unused N-length
  `MarkerScratch` buffers cost ~10 MB at N = 50,000 across 8 workers.

**`tsan.sh`'s `region` workload had never been run on this work.** It existed
all along (the review's claim that it did not is FALSE and propagated into a
session summary before being caught). Run it after ANY region-threading change.

## 2. The pooled Cauchy row at alpha 0.05: ANSWERED

`docs/state/prediction-cct-dependence-2026-09-12.md`. The recorded ASSUMPTION
is confirmed and quantified: it is the COMBINATION, not the inputs. CCT is
exactly uniform at rho 0 AND rho 1 (Cauchy stability, index 1), so the
distortion is non-monotone with an interior peak; measured on SAIGE's own
`get_CCT_pvalue` it peaks at rho ~0.5-0.7 and grows with k (0.0613 at k=18,
0.0698 at k=192) while alpha 0.01 stays within +0.002-0.005. The real null
flex192 run at N = 50,000 has measured inter-set correlation **0.705** -- that
same peak band -- and goes 0.0522 (sets) -> 0.0714 (pooled), which the
simulation reproduces. **`p_FlexRV` has size ~0.07 at alpha 0.05 on this design
and is accurate at discovery thresholds; never quote its nominal size at a large
alpha.** Not measured: size AT a discovery threshold directly (needs ~1e7 reps);
the tail claim rests on the alpha 0.01 column plus Liu & Xie.

## 3. The gate is 1.97x faster, and one lever was measured at zero

`harness/bench/NOTES.md`, "The pre-push gate".

- **`b16fa83c`: one generation, both comparators.** The two compare modes never
  differed in what they GENERATED. `SAIGE_COMPARE_MODE=both` generates once and
  applies both. **420s against 826s.** The requirement is TIGHTENED: each site
  asserts byte AND table-equivalence on the same bytes. Guarded in
  `test-meta-comparators.R` by the only discriminator that exists -- a pair that
  is numerically equal and byte-different.
- **`231ea222`: `build.sh --incremental`, 167s -> 57s**, refusing if any header
  is newer than the OLDEST object.
- **`fd464fc1`: running the two modes concurrently gains EXACTLY NOTHING**
  (851.5s both arms). The predictor that was right is saturation: one mode
  already runs at 93.8% CPU. Kept behind `SAIGE_PREPUSH_CONCURRENT=1`.
- **`d59e0ef4`**: the gate reports its own wall time and whether the tree was
  DIRTY (it used to print the last commit, which on 2026-09-12 named a commit
  that did not contain what was tested).
- CLAUDE.md now says when the full gate is NOT required.

Two facts any future gate timing needs: **page-cache warmth moves it 33%**
(more than any change tried), and measurement effort should scale with the
effect -- the paired alternating design that the concurrency null required was
ceremony at 2x.

## Open, and needing the owner

- **LEDGER #129 + #130 together** -- same site (the ultra-rare pseudo-marker
  pass). #130 needs a SIGN CONVENTION chosen (a flipped pseudo-marker's burden
  is a burden of the REFERENCE allele); #129 needs a choice between two fixes.
  No golden moves for either, so a fix needs a new golden that reaches the flip.
- **#116 + #128** -- #128's fix is best decided with #116, whose mechanism is
  still undiagnosed, so that package starts with a diagnosis.
- **The binary FlexRV arm at 50k** -- still the only untested scale arm.
- **The pooled row's size AT a discovery threshold**, if anyone wants it pinned
  rather than argued.
