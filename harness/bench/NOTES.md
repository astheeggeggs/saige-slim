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

## Step-1 blocked/memoised GetTrace: measured A/B (2026-09-01)

Branch `gettrace-blocked`. Workload is the `binary_denseGRM` golden config,
which is the `s1_full` bench workload; `bash harness/env.sh` pinned, one thread.
Both arms built on the SAME DAY on the same machine by swapping only
`src/SAIGE_fitGLMM_fast.cpp` between 34fa687d and the branch tip and doing a
full clean rebuild each way -- cross-day rows are not comparable on this
machine, which is why the `phase4-c3-*` rows could not serve as the baseline.

Reported figure is "fitting the NULL model took", 3 reps:

| | user (s) | elapsed (s) |
|---|---|---|
| before | 30.238 / 30.054 / 30.698 | 31.41 / 31.27 / 32.10 |
| after  | 22.601 / 22.572 / 23.598 | 24.16 / 24.64 / 25.74 |

**1.32x user, 1.27x elapsed, non-overlapping ranges.**

**This SUPERSEDES the "1.40x" in commit 9213390e's message**, which came from a
single rep per arm. The mechanism is unchanged and still reconciles -- the
varcomp block goes 2.469/10.080/10.068 -> 0.086/7.641/7.770 s per round, round 1
collapsing to a cache hit and rounds 2-3 each saving ~2.44 s against a predicted
30 * 0.086 = 2.58 s -- but the headline ratio was optimistic. One rep per arm is
not a measurement on this machine; the run-to-run spread here is ~2% within an
arm and was ~13% for the s2 workloads.

No `results.tsv` rows are added for this: `bench.sh` records
sys_s/maxRSS/pcg_iters from its own protocol, and back-filling those columns
with NA from a different harness would put two protocols in one table. Re-run
`bench.sh` on both builds if structured rows are wanted.

## S3, blocked multi-RHS PCG (2026-09-01, commit 32ed5f1e)

Same protocol as the section above: only `src/SAIGE_fitGLMM_fast.cpp` swapped
between 34fa687d and the tip, full clean rebuild each way, same day.

**Null model fit** ("fitting the NULL model took", user s, 3 reps):

| | rep1 | rep2 | rep3 | mean |
|---|---|---|---|---|
| before | 34.843 | 33.393 | 33.229 | 33.82 |
| after  | 14.160 | 14.305 | 14.413 | 14.29 |

**2.37x user, 2.34x elapsed**, non-overlapping. Per-round varcomp block
2.469/10.080/10.068 -> 0.085/4.181/4.225 s.

**Whole step-1 script** (`/usr/bin/time`, 2 reps): user 41.60/43.09 -> 27.58/
27.22 (**1.55x**); real 43.43/46.11 -> 30.78/29.44 (**1.49x**).

### Quote the right denominator -- this cost a false "overshoot"

The 2.37x first looked like it beat the 1.81x projection. It did not: the
projection was for the WHOLE RUN and the measurement was of the NULL FIT.
Redone on the right denominator the model is fine:

- fit-phase PCG iterations 273, varratio-phase 111, **both unchanged** by this
  change (it only touches GetTrace);
- null-fit budget is 442 - 111 = 331 GRM passes, 300 of them blockable;
- after = 7 blocked passes * 14.6 + 31 single = 133 units, i.e. **2.49x
  predicted vs 2.37x measured** -- reconciles to 5%.

The whole-run figure is lower than its own 1.81x projection for the opposite
reason: the pass-count model counts ONLY GRM passes and ignores R startup,
`library(SAIGE)`, the bim/fam parse and the ~2.2 s `setgeno` ingest, which
together are ~8.5 s of the 42 s baseline and are untouchable by any of this.

**So: 2.37x on the null model fit, 1.55x on the step-1 run as a user
experiences it.** Both are real; quoting either without saying which is what
made the numbers look inconsistent.

### Cumulative, against the original plan

| | null fit |
|---|---|
| plan's projection for blocked GetTrace | 3.6x |
| measured ceiling once D/A was known | 1.81x whole-run / 2.49x null-fit |
| S1+S2 (memoise + block Au) | 1.32x |
| S1+S2+S3 (+ blocked multi-RHS PCG) | **2.37x** |



## The sparse marker path, and what a 2% change looks like (2026-09-02)

`s2_marker_all` cannot see anything in the marker-path `scoreTest` or in
`getPCG1ofSigmaAndGtilde`: at default flags every marker goes to
`scoreTestFast_noadjCov` and the PCG never runs. **`s2_marker_all_sparse`** is
the same 128,868 markers through the sparse-GRM branch (`--sparseGRMFile` plus
the `binary_sparseGRM_vr` model, the only one whose variance-ratio file carries
both `sparse` and `null_noXadj` rows). Baseline at `b87fcaba`, 3 reps:
**5.449 s work** against `s2_marker_all`'s 2.654 s, i.e. the covariate-adjusted
score test plus one PCG solve per marker costs **2.05x** the fast path, ~22 us
per marker at N = 1000.

**Three reps could not resolve the F8 allocation change; eight paired reps
could.** Unpaired 3+3 gave 5.449 -> 5.443 s (0.1%, inside a 1.2% within-arm
spread) -- a non-result that reads as "no effect". Installing the two builds
side by side (`R CMD INSTALL -l <scratch lib>` for the baseline, the pixi env
for the change) and alternating arms within each rep gives:

| rep | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | mean |
|---|---|---|---|---|---|---|---|---|---|
| before | 5.368 | 5.307 | 5.396 | 5.315 | 5.435 | 5.367 | 5.482 | 5.382 | 5.382 |
| after | 5.262 | 5.301 | 5.194 | 5.177 | 5.240 | 5.277 | 5.569 | 5.201 | 5.278 |

paired mean difference **0.104 s = 1.9%**, t = 2.91 on 7 df, **p = 0.023**,
95% CI 0.4%-3.5%. Same conclusion as the reader comparison reached at N=50,000:
**pair the arms**. Two libraries and an alternating loop cost one extra install
and turn "unmeasurable" into a number with a confidence interval. The null
control for the setup is in the script: each arm prints the `find.package`
path it actually loaded, and the two `SAIGE.so` files are checksummed
different; the outputs of the two arms are byte-identical over all 128,868
markers, which is a stronger equality check than the golden comparison (same
workload, 128,868 rows, not 70).

Load average was 3-6 throughout, i.e. NOT a quiet machine -- which is exactly
why the pairing matters and why the unpaired comparison failed.

**Do not extrapolate the 1.9% to biobank N.** The mechanism is 3 allocations
and 2 n-double copies removed per PCG iteration plus one n-vector reciprocal
per call; at N = 1000 an n-vector is 8 KB, well inside malloc's fast path and
L1-resident, so this is close to the smallest the change can ever be. Above
n = 16,384 doubles each of those allocations crosses malloc's 128 KB mmap
threshold and becomes an mmap/munmap pair, and the copies become DRAM traffic.
That is a reason to expect the effect to grow, not a measurement of how much.

One bookkeeping caveat on those rows: both the `f8-baseline` and `f8-after`
labels in `results.tsv` carry commit `b87fcaba`, because the after-arm was
measured from the working tree before the change was committed. The LABEL is
the discriminator there, not the commit column, and the authoritative numbers
for that change are the paired ones in the table above, not the 3-rep rows.


## S1: what the Rcpp::List traffic in the SPA kernels was costing (2026-09-02)

Two instruments, because the goldens are byte-identical either way.

**Allocations, exactly.** `Rprofmem(threshold = 0)` around a whole-genome
`SPAGMMATtest` call (128,868 markers, 5,572 SPA-corrected) records every large
vector allocation and every new page of the small-vector pool:

| | R-heap page allocations | large-allocation buckets |
|---|---|---|
| before | 2,958 | 292x4048, 218x8048, 205x1072, 64x4056, 63x46912, ... |
| after | 1,736 | identical, bucket for bucket |

41% of the marker loop's R-heap page traffic gone, and the untouched
large-allocation histogram is the control: the change removed small-pool
traffic and nothing else. ~1,222 pages over 5,572 SPA markers is on the order
of 30 R allocations per SPA marker, matching the threading verdict doc's
estimate of ~20.

**Note the instrument's limits, twice over.** `Rprofmem` does NOT record every
small allocation -- only the pages the small-vector pool takes to serve them --
so 1,222 is a page count, not an object count. And the first instrument tried,
counting GCs with `gcinfo(TRUE)`, returned **zero GCs in both arms**: after
`ReadModel` has grown the R heap, 20k-odd small lists never reach the trigger.
That is a real result, not a broken probe -- the same probe reports 6 GCs for
50,000 list constructions in a bare R process. GC pressure is simply not where
this cost lives at N = 1000.

**Time: 11.8%, the largest step-2 win measured so far.** `s2_marker_all` with
`--is_output_moreDetails=TRUE`, two blocks of 8 paired reps with both builds
installed side by side:

| block | before mean | after mean | paired mean diff |
|---|---|---|---|
| before-arm first | 3.433 | 3.042 | 0.391 s |
| after-arm first (ordering control) | 3.601 | 3.163 | 0.438 s |
| combined, 16 pairs | 3.517 | 3.102 | **0.414 s = 11.8%** |

t = 7.53 on 15 df, 95% CI 8.4%-15.1%; outputs byte-identical across arms over
all 128,868 markers in both blocks. **The ordering control was necessary**: in
the first block the `before` arm always ran first in each pair, so a cache
effect would have inflated the result; swapping the order made the difference
slightly LARGER, which rules that out. Per SPA-corrected marker the saving is
~74 us -- large for ~30 allocations, and it is not really the allocations: an
`Rcpp::List` named assignment reallocates the list AND its names vector, and a
named read is a linear string search plus an R conversion, so the ~30
allocations come with ~30 name lookups. Cross-check: the profile above put
`SPA_fast` at 15% of the marker loop, and this removed ~12 percentage points of
the loop, i.e. most of what was in that frame was list plumbing rather than
saddlepoint arithmetic.

**This fraction will SHRINK at biobank N, not grow** -- the opposite of the F8
change. The list overhead is per SPA marker and flat in N, while the
saddlepoint arithmetic it sits beside is O(n) per marker. At N = 1000 the
plumbing dominated the kernel; at N = 500k it will not.

## Before ANY timing run, check the throttle as well as the load (2026-09-02)

`uptime` is not sufficient. Check `pmset -g therm` and require
`CPU_Speed_Limit = 100`.

Found the hard way during S6: the laptop's battery ran flat under the test
suite (8 workers on 8 cores, ~600-800% CPU sustained, ~150 step-2 child
processes per tier), and after it came back on AC the CPU was still clamped to
**60%**. Any measurement taken in that window is worthless, and — the part
worth remembering — **a paired A/B does NOT rescue it**. Pairing cancels noise
that is shared between the arms; a throttle that lifts or drops partway through
a block is not shared, it lands on whichever arm was running, which is exactly
the correlated-error case pairing cannot see.

Same reasoning applies to anything else with a duty cycle: video playback in a
browser measured at ~23% CPU here, which is small against 8 cores but BURSTY
(ads, resolution changes), so it is worth pausing for the few minutes a timing
block takes rather than modelling.

Memory was never the constraint: 36 R processes totalled 4.3 GB on a 16 GB
machine with 71% free. The cost of this suite is cores, not RAM.
