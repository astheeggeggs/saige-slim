# Benchmark data quality notes

## What the step-2 workloads can and cannot measure (2026-08-31)

**`s2_marker6k` and `s2_group` are ~90% fixed process overhead and cannot
detect a marker-kernel change.** Measured from `prof2.log` / a whole-genome run
at N=1000: ~2.7 s of startup (R, `library(SAIGE)`, `ReadModel`, bim/fam parse)
then ~0.33 s per 5,858-marker chunk, i.e. ~50 us/marker. `s2_marker6k` runs
exactly one chunk's worth, so only ~0.35 s of its ~5.6 s is marker testing; a
2x kernel speedup would move it ~3%, inside the run-to-run spread below. This
is why the `phase4-batch2` s2 comparison was recorded as "unresolved" -- it was
never resolvable, on any machine, however quiet.

Use **`s2_marker_all`** (all 128,868 markers, 22 chunks) for anything touching
the marker path, and read the **`work_s` / `throughput` columns**, not
`elapsed_s`. `bench.sh` derives them from the cumulative `proc.time()` SAIGE
already prints per chunk. `s2_marker6k` is retained only so the older
`phase4-*` rows stay comparable; do not quote it as a kernel ratio.

For step 1, read **`pcg_iters`** -- the summed PCG iteration count, which is the
GRM-pass budget (each iteration is one full pass over the GRM markers).

## Profiles

`sample_*.txt` are macOS `sample` call trees taken ad hoc during Phase 4, with
SYSTEM R (`/Library/Frameworks/R.framework/Versions/4.5-x86_64`) rather than the
pinned pixi env, at an unrecorded commit -- treat them as leads, not evidence.
`harness/profile.sh` now does this reproducibly against the pixi build.

The `prof*.txt` files are NOT profiles: they are step-2 result tables (ten
byte-identical 16 MB copies were committed in 64b0f677 and untracked again on
2026-08-31; `harness/bench/prof*` is now gitignored). The `prof*.log` run logs
for one step-2 and both step-1 runs are kept as the evidence behind the numbers
above.

### First reproducible profile: s2_marker_all, current build, N=1000

4,940 samples total, of which `mainMarkerInCPP` is 3,066 (62%). Within it:

| frame | samples | % of marker loop |
|---|---|---|
| `writeOutfile_single` | 682 | 22% |
| `SPA_fast` | 471 | 15% |
| `imputeGenoAndFlip` | 400 | 13% |
| `__cxa_throw` + unwinding | 269 | 8.8% |

**This corrects the lead noted in the C1 commit message.** The old
`sample_s2_marker.txt` showed 1455/4605 (~32%) in `__cxa_throw`; on the current
pixi build it is 8.8%. That profile was taken with a different R, a different
commit and a different dyld cache, so the 32% figure should not be quoted.
Exceptions on the SPA path are real but are not the headline.

**This ranking does NOT extrapolate to biobank N, and that is the point.** The
fixture is N=1000, where per-marker FIXED costs dominate. `writeOutfile_single`
formats ~19 fields per row at ~5 us/row regardless of N; the exception cost is
likewise per-marker. Both stay flat as N grows. What scales with N is the
genotype decode, `imputeGenoAndFlip`, the score test and SPA -- so at N = 500k
the top two rows here become rounding error. Optimise `writeOutfile_single`
only if small-N runs are themselves a target; its printed digits are the golden
bytes, so it is also the most hazardous thing in the table to touch.

- `phase4-baseline` (2 reps) and `phase4-batch1` (2 reps): CLEAN — machine quiet.
- `phase4-batch2` (2 reps): SUSPECT — a design workflow's agents were running
  concurrently; s1 numbers plausible, s2 numbers unresolved.
- `phase4-batch2b` (4 reps): CONTAMINATED — two editor agents ran concurrently;
  all workloads inflated with large variance. Disregard entirely.
- Rule (learned the hard way): benches are only valid with NO agents, builds,
  or test suites in flight. One bench at a time, quiet machine, and say so in
  the label if unsure.
- `phase4-batch2-quiet` (4 reps): CLEAN — no agents in flight. But the machine
  itself was intermittently noisy by this point in the day (see next entry).
- `phase4-batch1-requiet` (3 reps): CLEAN same-day A/B control — batch-1 HEAD
  rebuilt and re-benched on the same machine state as batch2-quiet. Result:
  huge variance (s1_full/nt1 72.5±35.4) — the machine, not the code. Verdict
  from the contemporaneous pair: batch 2 is equal-or-faster on dense step-1
  (1.38× nt1, far more stable) and par-within-noise on s2 workloads
  (s2_marker 0.87× on means, overlapping distributions, n=3/4 — unresolved).
  The tight morning `phase4-batch1` rows remain the best baseline estimates.
  Re-measure s2 on a quiet morning machine before the threading work.
