# Handoff: step 3b -- the per-cell tests of passes B and D threaded (2026-09-06)

One commit on `strip` after `3eeb554e`, gated. Plan:
`~/.claude/plans/can-you-summarise-what-shimmering-snail.md` (steps 1-3b
done, 4-5 not started). Previous handoff:
`docs/state/handoff-2026-09-06-region-threads.md`. Numbers:
`harness/bench/NOTES.md`, the last section.

## What landed

`mainRegionInCPP`'s ultra-rare pass (B) and BURDEN emitter (D) each end in
one `getMarkerPval` per (anno x max_MAF x weight) cell -- 15 on the bundled
and the `sim50k` fixtures -- and at step 3 those ran serially after the
threaded marker loop. Now each pass is three phases:

- **B1 / D1, serial, in cell order**: everything that precedes the test and
  is order-dependent or feeds it -- the `genoSumcount_noweight` and
  `genoSumMat` updates, the residualisation (SKAT-O), the running chunk index,
  the test context -- plus the test-independent writes (marker fields, MAC
  groups, case/control dosages). Each tested cell's operands go into a
  `URCell` / `BurdenCell` record (the collapsed vector after
  `imputeGenoAndFlip`, the index vectors, `altFreq`, `flip`, `w0`, the
  `TestContext`, the ER predicate).
- **B2 / D2, the tests**: `runCellTests(n, fn)` runs `fn` over the tested
  cells -- literally serial with the frame scratch at `--nThreads=1`, one cell
  per range on the pass-A `MarkerStatePool` otherwise (`RegionCellWorker`,
  `RefuseMainThreadGuard`, failure recorded per cell, lowest-index error
  rethrown, `NeedsMainThread` cells re-run on the main thread). A cell writes
  its own record and, in B2, its own `P1Mat` row / `P2Mat` column, so two
  cells never touch the same memory.
- **B3 / D3, serial, in cell order**: the pre-split output writes over the
  record values -- B3 the per-index vectors, the conditional columns and the
  G1 row; D3 is the old emitter loop verbatim, reading each cell's result
  instead of computing it, so the carry-forward cache sees the sequence it
  always did.

Every expression is the pre-split one over the same values; the output is
byte-identical whatever the thread count. `harness/tsan.sh plant` gained the
new worker as a fourth anchor (`m_fn((unsigned int)k, *state->scratch[0])`).

**The records persist across regions** (function statics, like the P1/P2
scratch): the collapsed vector is copied into an existing N-length buffer and
the index vectors into `std::vector`s that keep their capacity, handed to
`getMarkerPval` as strict non-owning `arma::uvec` views. The first version
copied, then moved, fresh N-length vectors per cell; both doubled the run's
sys time (30 x 400 KB mach VM allocations per pass per region -- see NOTES).
Persistence has one obligation: **every record's `tested` flag is reset at
the top of each region**, because D3 reads it per cell. The version without
the reset re-emitted an earlier region's BURDEN row for a cell the current
region did not test; `sim50k` caught it (the bundled BURDEN golden did not:
both its regions test the same cells) and
`test-step2-region-threads.R` now pins region-order independence with a
region whose missense cell is present but untested (drilled red on the
un-reset build, green after).

## Gate

- `test-step2-region-threads.R` 27/27 (P2 escalation, P3 conditional, P5
  BURDEN at 8 threads, P36/P37 quantitative, thread count asserted as
  printed; the region-order test at 1 and 4 threads); fast tier FAIL 0;
  lintr clean.
- `sim50k`, both traits, SKAT-O and BURDEN, 1 and 4 threads, every output
  `cmp`-identical to step 3's 1-thread outputs (`u1`/`v1`): the committed
  build 12/12 SKAT-O files (`gate2_s3b.log`) + 16/16 BURDEN files
  (`ab5_s3b.log`); the first build 24/24 (`gate_s3b_cells.log`). All under
  `harness/work/region_profile/`.
- Goldens all PASS; `prepush.sh` FAIL 0 in both modes (byte 2758, numeric
  2958 passes) on the committed build (and on the first build, 2746 / 2946).
- TSAN on the `region` workload, committed build: clean 0 / 0 races at 4 / 8
  threads with identical bytes; plant 102 / 103 races with DIFFERING bytes
  (first build 185 / 62; step 3's plant had 95 / 109 with three anchors).
- SKAT-O wall at 4 threads is unchanged (104.0 / 80.8 s against step 3's
  105.5 / 75.2 s, load 3-5): its cell tests were 2.6% / 7% of the region C++.
- Final `prepush.sh` with the two new tests: FAIL 0 both modes (2769 / 2969).
- **The cell worker's ER deferral is NOT reached by any bundled run.** Drill
  (`er_drill.py` in the session scratchpad: the worker records a
  `NeedsMainThread` as an error instead of deferring): P2 at 4 threads and P5
  at 8 stayed green; so did both with `--max_MAC_for_ER=200` (P33's value).
  `SAIGE_test.cpp:727` says why: the ER kernel runs only when the cell's
  score statistic exceeds the SPA cutoff, and resampling (the only part that
  needs R's RNG) only when the exact enumeration is infeasible, which these
  collapsed cells (MAC <= ~20 at N = 1,000) never are.
  A third drill with `--MACCutoff_to_CollapseUltraRare=100` (cells of MAC
  ~100) and `--max_MAC_for_ER=200` also stayed green: no bundled cell has a
  score past the SPA cutoff at a MAC that needs resampling. A fixture that
  reaches it is the missing piece of coverage for this path.
  The code path mirrors pass A's deferral, which P2 does exercise, and
  `test-step2-region-threads.R` now holds `--max_MAC_for_ER=200` runs of the
  SKAT-O and BURDEN cell tests byte-identical at 1 vs 4 threads so the arm
  stays covered if a fixture ever reaches it.

## Timing

Paired A/B against a scratch install of `3eeb554e`, `sim50k` BURDEN at 4
threads, arms alternated, 2 reps, both traits, load 5.7-7.7 (table and the
accounting in `harness/bench/NOTES.md`, last section):

| | step 3 | step 3b | ratio |
|---|---|---|---|
| y_binary, 4 threads | 36.3 / 35.5 s | **31.2 / 30.2 s** | 0.85-0.86 |
| y_quant, 4 threads | 19.5 / 18.9 s | **18.4 / 18.5 s** | 0.94-0.98 |
| 1 thread, either trait | equal within noise (wall, user, sys) | | |

Scaling step 3's recorded 1.57x / 1.64x by those ratios: **~1.84x / ~1.71x**
at 4 threads. The 2.2x prediction is not met, for a located reason: the
per-cell tests were 22% / 9% of the 1-thread run, so 14% / 6% of the wall was
the most their threading could remove at this inflation floor, and that is
what was removed; the rest of the prediction assumed step 3's 1.57x, which
measures 1.43x today on a box with ~1.6 cores taken by VS Code's C++ indexer
and Spotlight. User CPU is +12% at 4 threads (cell tests on workers: residual
allocation contention plus TBB spin between three parallel regions per
gene). Every intermediate build's numbers are in
`harness/work/region_profile/results_s3b_ab{,2,3,4,5}.tsv` with the loaded
library's md5 on each row.

## Not done

- **`docs/state/next.md`'s pointer line only** was updated (another session
  holds an uncommitted section in that file); the step-3b item under
  "Decided order" should be marked done by the owner.
- Steps 4 and 5 as the plan describes them.
