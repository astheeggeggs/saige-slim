# Handoff 2026-09-05 (evening): the binary single-variant scan at N = 50,000

Read after `docs/state/handoff-2026-09-05-review.md`. This session profiled the
binary single-variant scan on the `sim50k` fixture, explained why it read as
"not faster than upstream", removed the allocation churn from the Firth and SPA
kernels byte-neutrally, and benchmarked two non-neutral Firth variants behind
flags for the owner to decide on. Full numbers: `harness/bench/NOTES.md`, section
"The binary single-variant scan at N = 50,000"; the strip-vs-upstream table's
binary row is corrected in an addendum to `harness/sim/upstream_compare/RESULTS.md`.

## What was found

- **Firth is 56% of the binary scan's CPU** (17.5 ms per fit, 9 Newton
  iterations, 7,891 fits at `--pCutoffforFirth=0.1`), decode 15%, SPA 10%.
  Per Firth-eligible marker the cost is 33 ms on BOTH arms (per-chunk
  regression of the existing logs), so it was never a strip regression.
- **The 0.94x binary row was not like-for-like**: upstream's LEDGER #28 deflates
  its AF > 0.5 p-values, so it ran 569 Firth fits and 276 SPA corrections there
  against strip's correct 1,633 and 730. Like-for-like the arms were at parity.
- **The +48% user CPU at 4 threads** was large-allocation churn (mach VM
  map/unmap of >127 KB Armadillo temporaries, ~90 per Firth marker, sys 19 s ->
  61 s), on top of a uniform 1.3-1.4x turbo/hyperthreading inflation that the
  quantitative scan also shows and that is this 4-core i5, not code.

## What changed (uncommitted at handoff; three logical changes: items 1+2, item 3, items 5+6)

**Items 1+2, byte-neutral, gated:** `SAIGE::KernelWorkspace` (src/SAIGE_test.hpp)
on `Marker::MarkerScratch`, passed to `getMarkerPval` via a trailing
`KernelWorkspace * t_ws = nullptr` (the 9 region-path call sites are untouched
and build one lazily per call). `fast_logistf_fit_simple` now takes the adjusted
genotype and the workspace and keeps its n x 2 design, weighted design, Q and
five n-vectors across iterations; the six `Korg/K1/K2` kernels have `_ws`
twins writing into strict views of four shared buffers (`SPAWorkspace`,
src/SPA_binary.hpp), with `SPA_ws`/`SPA_fast_ws` in src/SPA.cpp. Every
expression is the original verbatim; the exported R-visible kernels keep their
signatures and build a local workspace. `SPA_binary.hpp` and `SPA.hpp` gained
include guards (they are now included twice). Gate: 50k binary output
`cmp`-identical at 1 and 4 threads to the pre-change run; fast tier 0 failures;
standard tier 0 failures / 2,718 passes on the items-1+2 build, and again
0 failures / 2,718 passes on the final tree with both flags at their defaults.
1-thread user CPU 240-268 s -> 205 s, sys 13-19 -> 5 s; 4-thread wall 111-118 s
-> 83.5 s, sys 49-55 -> 11 s.

**Items 5+6, DECIDED 2026-09-06: `closed` is the default, the warm start is gone.**
`--firth_solver=closed|qr` (extdata/step2_SPAtests.R -> `SPAGMMATtest(firth_solver=)`
-> `pushNullModelsToCPP` -> `setSAIGEobjInCPP`/`addSAIGEobjInCPP` -> `m_firthSolver`;
RcppExports regenerated). `fast_logistf_fit_closed` replaces the QR of the n x 2
weighted design with the 2 x 2 Fisher matrix and closed-form hat values; `qr` keeps
`fast_logistf_fit_simple`, the SAIGE 1.x path. What the decision rested on, all on
the 50k binary scan at 1 thread:

| flags | user s | output vs the QR default |
|---|---|---|
| `--firth_solver=closed` | 170.8 vs 235.5 (-27%) | **byte-identical**, 70,709 rows, 7,891/7,891 converged |
| `--firth_warm_start=TRUE` (deleted) | 172.4 | 137 rows differ, 118 fits NOT converged, BETA up to 4x off with sign flips (MAC 20-88) |

Gate for the switch: the Firth golden P31 (`plink_marker_firth.txt`) is byte-identical
under the new default (standard tier unchanged, 2,724 passes), the 50k binary output is
byte-identical to the pre-session run, and `prepush.sh` is green in both modes (byte
2,724 / numeric 2,924). `tests/testthat/test-step2-firth-solver.R` (standard tier)
pins that `--firth_solver=qr` still reproduces P31 as a PARSED comparison at 1e-10 --
never bytewise, two paths -- and that an unknown solver name is refused. The warm
start was removed rather than left switchable: the null intercept was the right
start value (the offset is centred), so the failure is the penalised Newton step from
~0.01 weights saturating at `maxstep` on rare markers, not a plumbing error, and a
flag that breaks 1.5% of fits has no use.

`harness/sim/04_run_saige.sh` gained `${S2EXTRA:-}` on the single-variant
command so flags can be appended by the benchmark driver.

**Item 3, byte-neutral, gated:** `PlinkClass::countOneMarker` tallies the four
2-bit codes (a packed-lane byte table on the whole-file-in-order path, a
per-sample loop otherwise) and `prepareOneMarker` decides QC gate 1 from them
before any genotype is expanded, through the same `summariseCounts` arithmetic
the full read uses; `GenoReader::countOneMarker` defaults to "unsupported" so
PGEN is unchanged, and the pre-screen is skipped under the male X non-PAR
rewrite. Quantitative 50k scan 60-65 s -> 45.3 s user; binary unchanged within
noise (decode was 15% of it). Outputs byte-identical for both traits at 50k,
`prepush.sh` green in both modes. The first version (four masked increments per
byte) bought nothing -- see NOTES.md -- which is why the table exists.

R lint: `R/MultiModel.R` and `extdata/step2_SPAtests.R` clean;
`R/SAIGE_Test_main.R` at 108 lints before and after (legacy file, no new ones).
lintr lives in a scratch lib this session installed under its scratchpad; it is
not in the pinned env.

## Not done

- `getGPosNeg_Binom` (14% of `SPA_fast`) and the four per-SPA-marker
  `gNA/gNB/muNA/muNB` gathers still allocate; small, and a hand loop would have
  to reproduce `accu`'s two-accumulator order to stay byte-identical.
- The case/control gather (~1.5%) is unchanged; `arma::mean` has the same
  accumulation-order constraint.
- PGEN has no counting pre-screen; pgenlib's decompressed genovec could be
  tallied the same way (`GenoarrCountFreqsUnsafe`), saving the expansion only.
