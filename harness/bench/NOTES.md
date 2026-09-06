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

## Step-2 marker parallelism: both mechanisms measured (2026-09-02, S6/S7)

`harness/bench_threads.sh <label> <reps> <threads> [inproc|fork]`. Paired A/B in
the sense this file established: the two arms run adjacent within each rep and
the arm ORDER is swapped between the two blocks. Unlike the F8/S1 comparisons no
second installed library is needed -- both arms are the same binary with a knob
moved -- which removes a confounder rather than adding one. Byte-identity
between the arms is checked on every pair; a speedup measured on a run that
stopped producing the right answer is not a speedup.

Three workloads over the same 128,868 markers, because the ER finding means one
workload cannot represent the change: `binary` (default `max_MAC_for_ER=4`, so
ER-eligible markers are deferred to the main thread), `binary_sparse` (same,
but with the PCG and the marker-path `scoreTest` on the clock, so the
threadable fraction is much larger), and `quant` (never enters the ER branch).
16 pairs each, `CPU_Speed_Limit = 100` throughout.

**In-process (`nThreadsMarker`), 1 vs 4 threads, on `work_s`:**

| workload | 1 thr | 4 thr | speedup | paired p | maxRSS |
|---|---|---|---|---|---|
| binary        | 2.725 | 1.316 | **2.07x** | 4e-15   | 274.6 -> 275.3 MB |
| binary_sparse | 5.755 | 2.277 | **2.53x** | 4.1e-14 | 276.3 -> 276.8 MB |
| quant         | 1.973 | 0.942 | **2.09x** | 7.4e-16 | 270.9 -> 271.2 MB |

The ordering control behaves: block 1 (serial arm first) and block 2 (threaded
arm first) give the same paired difference to within 7%, so the effect is not a
warm-cache artefact of whichever arm ran second.

**Head to head against the fork path, at 4 threads.** Both mechanisms were run
against the SAME serial baseline rather than against each other, so each is
measured through a common reference:

| workload | inproc elapsed | fork elapsed | difference | maxRSS: inproc / fork |
|---|---|---|---|---|
| binary        | 6.121 | 6.358 | -0.237s, p = 0.14 | +0.7 / **+89.6 MB** |
| binary_sparse | 7.177 | 7.076 | +0.101s, p = 0.53 | +0.6 / **+89.6 MB** |
| quant         | 5.793 | 5.897 | -0.104s, p = 0.49 | +0.3 / **+47.4 MB** |

**On wall clock the two mechanisms are indistinguishable** -- all three p-values
are null, and the serial baselines of the two blocks agree (7.71/7.86,
10.74/10.69, 6.69/6.77), which is the internal check that the blocks are
comparable at all. **They differ only in memory**, and there in-process wins by
two orders of magnitude: fork's copy-on-write children cost +47 to +90 MB at
N = 1000 and 4 shards, and that grows with both N and the shard count, while
in-process costs under 1 MB.

**`work_s` IS NOT VALID FOR THE FORK ARM** and the table above deliberately does
not use it. `s2_breakdown` recovers startup and work from SAIGE's per-chunk
cumulative `proc.time()` lines; under fork those lines come from four different
processes with four independent cumulative clocks, so the decomposition is
meaningless there. `elapsed_s` is the only metric the two mechanisms share.

**What this fixture CANNOT tell you, which is the important part.** Only
**5,292 of 128,868 markers (4.1%)** have MAC <= 4 on this fixture, so the ER
serialisation that in-process pays -- and fork does not -- is a 4% tax here.
That is why in-process matched fork even on the binary workloads. It is not a
rare-variant scan, and the regime the retirement criterion warns about is
exactly the one this fixture cannot reach. With a deferred fraction f, the
in-process speedup at T threads is bounded by `1/(f + (1-f)/T)`: at f = 0.041,
T = 4 that bound is 3.6x (not binding); at f = 0.5 it is 1.6x; at f = 0.9 it is
1.08x, i.e. threading buys nothing.

The two workloads reconcile under one model, which is the check that the above
is arithmetic rather than storytelling. Quant (f = 0) measured 2.09x on work_s,
implying a non-ER serial fraction s = 0.305 -- the writer and the per-chunk R
work, neither of which threads. Feeding that s and f = 0.041 back in predicts
**2.00x** for binary; measured **2.07x**. So the ~2x ceiling at 4 threads is set
by the serial writer, not by ER, on this fixture.

## The f-sweep: the measurement that actually discriminates (2026-09-02)

The table above ends with "this fixture cannot reach the regime that decides
it". `harness/make_rare_fixture.py` builds that regime: five marker files, same
1000 samples, same 128,868 markers, same step-1 model, differing ONLY in the
fraction f of markers that are ER-eligible (binary, MAC <= 4). Every marker is a
byte-for-byte copy of a real marker's genotype block from the bundled fixture --
nothing is simulated -- drawn from a rare pool (5,292 real markers with
1 <= MAC <= 4) and a common pool (123,576 with MAC > 4), interleaved uniformly
so every `parallelFor` range sees the same f.

Two things had to be got right or the answer is an artefact. **Interleave, don't
concatenate**: contiguous ranges mean a file with its rare markers in one block
dumps the whole serial burden on one worker and measures scheduling.
**Renumber positions**: `SPAGMMATtest` sorts by CHROM/POS, so duplicated markers
carrying their original positions get reordered and the interleaving is lost.

The cheap alternative -- raise `--max_MAC_use_ER` on the existing fixture -- was
tried and rejected on inspection. ER cost climbs steeply with MAC and above ~21
carriers the exact test switches to the resampling arm, so raising the cutoff
changes the WORKLOAD rather than only f. Validated: cutoff 4/35/81/264/1000
gives f = 0.041/0.25/0.50/0.75/1.00, but the deferred markers at the top of that
range cost nothing like the ones at the bottom.

**In-process, 4 threads, `work_s`.** The model is fitted on the f = 0 point ONLY
(one parameter, s = the non-ER serial fraction) and then PREDICTS the other four:
`speedup = 1/(s + (1-s)f + (1-s)(1-f)/T)`, s = 0.285.

| f | work 1 thr | work 4 thr | speedup | predicted | ord. control b1 / b2 |
|---|---|---|---|---|---|
| 0.00 | 3.035 | 1.407 | 2.157 | 2.157 (fit) | 1.662 / 1.594 |
| 0.25 | 3.113 | 1.902 | 1.637 | 1.673 | 1.294 / 1.128 |
| 0.50 | 2.725 | 2.020 | 1.349 | 1.367 | 0.744 / 0.666 |
| 0.75 | 2.493 | 2.237 | 1.114 | 1.155 | 0.205 / 0.306 |
| 1.00 | 2.304 | 2.657 | **0.867** | 1.000 | -0.261 / -0.445 |

The one place the model misses is f = 1, and it misses in the direction and by
the mechanism predicted in advance: a deferred marker is decoded on a worker,
thrown away, then re-read and re-tested on the main thread, so at f = 1
in-process is not merely no faster but **13% SLOWER than serial**. The ordering
control agrees in both blocks at every f, including the negative differences.

**Head to head, `elapsed_s` (the only metric both mechanisms share).**

| f | in-process | fork | maxRSS in-proc / fork |
|---|---|---|---|
| 0.00 | **1.248** | 1.213 | 268 / 367 MB |
| 0.25 | 1.123 | 1.200 | 269 / 354 MB |
| 0.50 | 1.076 | 1.220 | 269 / 353 MB |
| 0.75 | 1.020 | 1.159 | 268 / 366 MB |
| 1.00 | **0.954** | 1.170 | 269 / 350 MB |

**Fork's speedup is FLAT in f** (1.16-1.22 across the whole range) -- each child
is its own process with its own R main thread, so the ER deferral costs it
nothing -- while in-process falls monotonically to a net slowdown. That is the
discrimination the bundled fixture could not provide, and it says plainly that
the two mechanisms are **not substitutable**. f = 0.041 is precisely where they
look identical.

Read the SPEEDUP RATIOS, not the raw elapsed means, and not the head-to-head
p-values: the ratios are paired with the arms adjacent, whereas a cross-mechanism
comparison is unpaired across blocks run at different times under different
ambient load, which is why those p-values scatter (0.06, 0.01, 0.68, 0.04, 2e-6)
while the ratios are clean and monotone.

## Step B: the ER deferral removed, and the cliff with it (2026-09-03)

The sweep above was re-run unchanged after `Main.cpp` stopped deferring binary
`MAC <= max_MAC_for_ER` markers to the main thread. Same fixtures, same machine,
12 pairs per cell, `work_s` speedup at 4 threads:

| f | before Step B | after Step B |
|---|---|---|
| 0.00 | 2.157 | 2.063 |
| 0.25 | 1.637 | **2.080** |
| 0.50 | 1.349 | **2.171** |
| 0.75 | 1.114 | **2.233** |
| 1.00 | **0.867** | **2.077** |

Spread across f collapses from **1.290 to 0.170**: in-process is now flat in f,
which is the property that made fork the only viable mechanism for rare-variant
binary work. At f = 1 -- every marker taking the exact test -- it goes from 13%
SLOWER than serial to **2.08x faster**, a 2.4x swing on the workload SAIGE is
mostly used for. The f = 0 cell moved -0.094, which is within the run-to-run
variation of this bench and is not read as a regression.

**Head to head, elapsed_s at 4 threads:**

| f | in-process | fork | maxRSS in-proc / fork |
|---|---|---|---|
| 0.00 | 1.243x | 1.224x | 270 / 368 MB (+98) |
| 0.25 | 1.222x | 1.214x | 270 / 355 MB (+84) |
| 0.50 | 1.222x | 1.215x | 270 / 353 MB (+84) |
| 0.75 | 1.231x | 1.178x | 270 / 367 MB (+97) |
| 1.00 | 1.173x | 1.166x | 269 / 350 MB (+81) |

In-process now matches or slightly beats fork at every f while using 80-98 MB
less. The trade-off that made the retirement hard is gone; what remains before
fork can actually be deleted is criterion item 2 (feature parity), not a
measurement.

## PGEN threads too, and it costs almost nothing here (2026-09-03, S8)

`PgenReader::clone()` landed, so `--nThreads` now threads a pgen marker scan.
Paired A/B on the 128,868-marker fixture, 6 reps x 2 blocks with the arm order
swapped between blocks, 1 vs 4 threads, `CPU_Speed_Limit = 100`:

| workload | n | work_s speedup | 95% CI | elapsed_s | maxRSS 1t -> 4t |
|---|---|---|---|---|---|
| `pgen` (binary) | 12 | **2.177x** | 2.051-2.295 | 1.256x | 245.6 -> 246.2 MB |
| `pgen_quant` | 12 | **2.040x** | 1.922-2.150 | 1.186x | 241.7 -> 242.0 MB |

That sits alongside PLINK's 2.06-2.23x at 4 threads on the same fixture, which is
the expected answer: threading the marker loop is reader-agnostic once the reader
can hand out per-worker handles. `elapsed_s` is much lower than `work_s` for the
usual reason -- ~2.1 s of the ~7 s run is process startup, which threading cannot
touch (see the top of this file).

**The memory line needs a caveat, and it is the one thing that differs from
PLINK.** PLINK clones share one mmap and cost nothing per worker. A PGEN clone is
a whole second `plink2::PgenReader` with its own arena: **~130 bytes per SAMPLE
per handle** (arithmetic from the allocation expression in
`PgenClass::initReader`, not a measurement), most of it two transpose buffers
that `Read()` never touches. At N = 1000 that is ~0.13 MB per thread, which is
exactly why the table above shows nothing. At biobank N it is tens of MB per
thread, and **nothing here can measure that** -- the bundled fixture has 1,000
samples. The obvious follow-up is to stop allocating the buffers the
single-variant read path never uses.

**A measurement note that cost two runs.** `harness/bench_threads.sh` was still
driving its threaded arm through `SAIGE_MARKER_THREADS`, which S7 deleted, so
every "threaded" arm it produced after S7 would have been serial; `harness/tsan.sh`
had the same defect. Both now pass `--nThreads`, and both keep their
DID-NOT-THREAD check, which is what makes the failure visible rather than silent.
The `fork` mechanism dimension of `bench_threads.sh` refuses with a pointer to
S7 rather than pretending.

---

## Group tests at N = 50,000: the first region-path measurement (2026-09-03)

Every number above is N = 1,000 and single-variant. This is the group path at
`sim50k` (N = 50,000, 210 genes, 5 annotations x 3 max_MAF, y_binary,
single-threaded under `harness/env.sh`), and it settles what `--is_fastTest`
does for a group run. Written up as **LEDGER #51**.

**Environment caveat, stated first.** Load average was 5-6.5 on 8 cores
throughout (Chrome + WindowServer, fluctuating), which is exactly the scale
CLAUDE.md warns has moved measurements by more than the effect being measured.
Mitigations: **user CPU seconds, not wall clock**; two reps per arm; arms
interleaved with the order swapped between blocks (A,B,C then C,B,A) so drift is
shared rather than assigned. Replicate spread came out 4-6%, which is the
resolution -- do not read differences smaller than that.

| arm | config | escalated | region passes | user_s (mean of 2) |
|---|---|---|---|---|
| A | `--is_fastTest=TRUE`, cutoff 0.05 (default) | 22 | 232 | 187.4 |
| B | `--is_fastTest=TRUE`, cutoff 1 | 156 | 366 | 286.8 |
| C | `--is_fastTest=FALSE` | 0 | 210 | 161.1 |
| D | no sparse GRM (`scoreTestFast` throughout), 1 rep | n/a | 210 | 147.7 |

**A, B and C are byte-identical** -- region file and all 9,050 `singleAssoc`
rows. That is the finding: arm A does 22 extra full region passes and arm B does
156, and neither changes a digit, because the region path picks its variance
ratio per marker from MAC and never reads the fastTest flag (LEDGER #51).

**The reconciliation, which is why this is arithmetic and not a story.** On an
equal-cost-per-pass model, per-pass cost derived independently from each arm is
0.767 (C), 0.808 (A) and 0.784 (B) seconds; predicted-vs-measured ratios are
A/C 1.105 vs 1.163, B/C 1.743 vs 1.781, B/A 1.578 vs 1.530. All within the
replicate spread. Three arms, one parameter, no leftovers.

**Two numbers to quote:**
- `--is_fastTest=TRUE` costs **+16.3%** on a group run and returns identical
  output.
- Swapping the PCG score test for `scoreTestFast` across a whole pass (arm D)
  saves only **8.3%**. That is the CEILING on any repair to fastTest here, and
  it is the more important number: **~90% of gene-test time is not the
  per-marker score test.**

**So the next group-path profile should not go looking at the score test.** The
unaccounted ~90% is some mix of the ultra-rare collapse, the P1Mat/P2Mat gene
covariance assembly, and the R-side SKAT-O post-processing
(`R/SAIGE_SPATest_Region_Func.R`). Which of those it is decides whether
threading the C++ region loop is worth anything at all -- if the time is in R,
it is not.

**Per-gene cost, for extrapolation:** ~0.77 s CPU per region pass at N = 50,000
with this annotation/max_MAF grid. A 20,000-gene exome is then ~4.3 CPU-hours
per phenotype. `[EXTRAPOLATION: linear in gene count from 210 genes, and it
assumes this fixture's gene-size distribution. Not measured at exome scale.]`

**Escalation rates** (free, from the pipeline's existing logs -- `grep -c
"Non-fast test is performed"`), the same 210 genes:

| trait | escalated | rate |
|---|---|---|
| `y_null` | 6/210 | 2.9% |
| `y_quant` | 10/210 | 4.8% |
| `y_binary` | 22/210 | 10.5% |
| `y_binary_pc` | 27/210 | 12.9% |
| `y_quant_pc` | 31/210 | 14.8% |

`y_null` is a clean negative control: escalation fires on Cauchy p < 0.05, so
under the null it should fire at ~5% (10.5 of 210 expected, 6 observed, -1.4 sd
on Binomial(210, 0.05)) -- conservative, and nowhere near fabricating
escalations. The signal arms span 4.8-14.8%, so the rate IS signal-dependent
across ~5x, which is worth knowing before quoting a per-gene cost from any one
fixture. But the aggregate effect is small: at ~10% escalation a whole run moves
~10%, NOT the ~2x that an earlier draft of `docs/state/next.md` claimed. The 2x
is per escalated gene only.

---

## Where the OTHER 90% of a gene test goes (2026-09-04)

The section above ended with a question: the per-marker score test is only 8.3%
of a group run, so what is the rest? It named three candidates -- the ultra-rare
collapse, the P1Mat/P2Mat assembly, and the R-side SKAT-O post-processing -- and
said which one it is decides whether threading the C++ region loop is worth
anything.

**It is the third, and it is not even our code.** `Rprof` over the same `sim50k`
arm A (210 genes, 5 anno x 3 max_MAF, `y_binary`, `--is_fastTest=TRUE`, 22
escalations, single-threaded under `harness/env.sh`):

| bucket | rep A | rep B |
|---|---|---|
| C++ `mainRegionInCPP` | 95.0 s, **43.9%** | 92.9 s, **43.8%** |
| R `region_postprocess` | 116.3 s, **53.8%** | 114.8 s, **54.1%** |
| .. of which `get_SKAT_pvalue` | 105.5 s, 48.8% | 104.4 s, 49.2% |
| .. of which `SKAT_davies` (`.C`) | 97.0 s, **44.9%** | 96.1 s, **45.3%** |
| .. R glue (post minus SKAT) | 10.9 s, 5.0% | 10.4 s, 4.9% |
| everything else | 4.8 s, 2.2% | 4.5 s, 2.1% |

The single largest line item in a SKAT-O group run is
`SKAT:::Met_SKAT_Get_Pvalue` -> `SKAT_META_Optimal` -> `SKAT_Optimal_PValue_Davies`
-> `integrate()` -> `SKAT_davies`, i.e. **R adaptive quadrature over the SKAT
package's Davies routine, 45% of the run**. `region_postprocess`'s own R code is
5%.

**Three things make this trustworthy rather than a profiler story.**
- *The profiled run is the same run.* 188.1 s user against NOTES's 187.4 s for
  arm A, and the same 22 escalations, so `Rprof` at a 0.02 s interval costs
  nothing measurable and this is comparable to every number above.
- *Replicates.* A and B agree within 0.3 pp on every bucket. Proportions within a
  run are far more stable than the 4-6% cross-run spread this file warns about,
  because both halves are slowed equally by load (which was 3.7-7.4 throughout).
- *A negative control that predicts a zero and gets one.* The same command with
  `--r.corr=1` (BURDEN, which never calls SKAT-O) puts `SKAT_davies` at **0.0%**
  and `region_postprocess` at 0.2%, with C++ at **93.7%**. If the attribution
  were an artifact of sampling inside `.C`/`.Call`, it would not switch off.

**So: threading the C++ region loop is worth little for SKAT-O and a great deal
for BURDEN.** Amdahl on the measured C++ fraction, threads on the region loop
only:

| test type | C++ fraction | 4 threads | 8 threads | infinite |
|---|---|---|---|---|
| SKAT-O | 43.9% | **1.49x** | 1.62x | 1.78x |
| BURDEN | 93.7% | **3.36x** | 5.55x | 15.9x |

That is the answer to the question this file asked, and it inverts the natural
plan: the C++ region loop is NOT the thing to thread first if users run SKAT-O.

**What the 45% actually is, and why it is attackable.** Each (annotation,
max_MAF) cell of each region gets its own SKAT-O p-value, and each one is an
independent 1-D integral. This grid is 5 x 3 = 15 cells, so the run computes
about 210 x 15 + 22 x 15 = 3,480 of them in ~105 s, i.e. **~30 ms each**
`[ARITHMETIC from the grid and the escalation count, not a counted call tally]`.
They share no state, so they are embarrassingly parallel and, unlike a float
reduction, computing them concurrently cannot change a digit.

**How to reproduce, since it needs no code change.** Put an `Rprof` around the
CLI script rather than instrumenting the pipeline -- the script's own `optparse`
reads the same `commandArgs`, so the arguments pass straight through:

```r
# wrap.R, run with cwd = extdata and the ordinary step2 arguments appended
Rprof(Sys.getenv("SAIGE_PROF_OUT"), interval = 0.02)
source("step2_SPAtests.R")
Rprof(NULL)
```

Then read it with `summaryRprof()$by.total`, whose rownames are QUOTED
(`"region_postprocess"`, not `region_postprocess`) -- looking them up unquoted
returns zeros for every bucket and a tidy, entirely wrong "100% in everything
else" table. Time spent inside `.Call`/`.C` is attributed to the enclosing R
frame, which is exactly what makes the C++/R split readable this way.

**Follow-up: where inside the 45% (same day, per-cell trace).** `get_SKAT_pvalue`
and `SKAT:::SKAT_davies` were shadowed in their namespaces for one run (no
package edit; the region file and all singleAssoc rows came back byte-identical
to the uninstrumented run, which is the control that the instrumentation was
inert). 2,527 SKAT-O cells, **1,123,456 `SKAT_davies` calls** -- 445 per cell,
~94 us each. So the cost is the NUMBER of integrand evaluations, not the size of
any one of them.

- **Cost is diffuse, not concentrated.** Top 1% of cells hold 8.9%, top 10% hold
  48.4%, top 25% hold 75.6%. There is no small hot set to special-case, which
  means parallelising across cells or regions is the shape that would work, and
  a targeted fix for "the big genes" is not.
- **Cost is NOT monotone in cell size.** m = 1 cells are free (0.49 ms, one
  davies call). The most expensive band is m in (2,5] at 170 ms/cell and 39.9% of
  all the cost; m in (100,250] is 100 ms/cell and 1.3%. Whatever this is, it is
  not O(m).
- **Cost is ANTI-correlated with significance.** Cells with p < 1e-6 cost 11.4 ms
  and 137 davies calls; p in (0.1, 0.5] cost 54.1 ms and 541 calls. That is the
  opposite of the obvious guess and it kills it: `SKAT_Optimal_PValue_Davies`
  passes `abs.tol = 1e-25`, which looks like the culprit until you notice the
  integral it evaluates is `1 - p`, so `rel.tol * |result|` -- the term that
  actually binds -- gets TIGHTER as p grows, not as p shrinks. `[The mechanism is
  consistent with the measured direction; it was not tested by varying the
  tolerance, and the p > 0.5 band being cheaper again (36.5 ms) does not fit it.]`

**Would a C++ rewrite help? No -- the hot code is already C.** `SKAT_davies` is a
thin R wrapper around `.C("qfc", ...)`, Davies' 1980 algorithm in compiled C.
Self times from the profile above say where the 45% actually sits:

| | self |
|---|---|
| `.C` (qfc, compiled) | 93.90 s |
| `SKAT_davies` R wrapper | 2.38 s |
| integrand `f` (R) | 1.50 s |

So the R glue is **~4% of the 45%**, i.e. **~1.8% of the run**. Rewriting it in
C++ for speed is not a lever. (C++ is still the only way to get IN-PROCESS
threading over the cells -- see below -- but that is an argument about
parallelism, not about C++ being faster.)

**Can the accuracy be loosened instead? No -- measured, and it breaks the
p-values that matter.** The integrand asks Davies for `acc = 1e-6` while the
outer `integrate()` only resolves `rel.tol` ~ 1.2e-4, a 100x mismatch that looks
like free money. Sweeping it (an identical copy of
`SKAT_Optimal_Integrate_Func_Davies` with that one constant lifted to a variable;
the `acc = 1e-6` arm reproduces the unmodified run BYTE-IDENTICALLY, which is the
control that the copy is faithful):

| `acc` | user CPU | davies calls | printed p unchanged | max rel diff among p < 1e-4 |
|---|---|---|---|---|
| 1e-6 (SKAT's) | 174.5 s | 1,123,456 | -- | -- |
| 1e-5 | 160.2 s | 1,304,193 | 1911/3338 | **8.0e-01** |
| 1e-4 | 132.9 s | 1,331,294 | 1720/3338 | **5.0e-01** |
| 1e-3 | 122.3 s | 1,330,885 | 1718/3338 | **5.9e-01** |

24-30% cheaper and up to **80% relative movement on genes with p < 1e-4**.
Rejected. Note the call count RISES as `acc` loosens -- each call is cheaper but
the outer quadrature subdivides more, so the saving is much less than the
accuracy ratio suggests. Do not re-open this without a reason to believe the
significant tail is safe.

**So the lever is parallelism, and the vehicle matters more than the language.**
The cells are independent, so any correct scheme is bit-identical. Measured on
the real workload (445 `SKAT_davies` calls per cell = the measured mean, 12 cells
per region, in a process carrying the 13 MB model plus ballast, at load average
6-8.6 on 8 cores, so these are PESSIMISTIC):

| scheme | 4 workers | 8 workers | overhead per region |
|---|---|---|---|
| `mclapply` (fork per region) | **1.12x** | 0.91x | 55-90 ms |
| `makeForkCluster` (fork once) | **2.81x** | 3.02x | 10-20 ms |
| `makeCluster` PSOCK (once) | **3.02x** | 3.18x | 8-18 ms |

**Fork-per-region is dead** -- 55-90 ms of overhead against ~350 ms of work per
region, and at 2 and 8 workers it is a net SLOWDOWN. A PERSISTENT cluster is
another matter: ~3x on the parallel part for ~10 ms a region. Against the whole
run that is 1/(0.55 + 0.45/3) = **~1.43x**, less ~1% for the cluster overhead
(210 regions x ~10 ms on a 188 s run).

In-process C++ threading over the cells would reach 1/(0.55 + 0.45/8) = 1.64x
instead of 1.43x -- 15% more -- but it means reimplementing
`Met_SKAT_Get_Pvalue`, `SKAT_META_Optimal`, the Liu fallback and the QUADPACK
integration, and vendoring Davies' `qfc`, with every group golden encoding SKAT's
exact output as the gate. `[Also unchecked: whether `qfc` is thread-safe. It is
Fortran-derived C from 1980 and may carry file-scope state; that has NOT been
looked at, and it would have to be before any in-process scheme.]`

**A silent failure found by accident while doing this.** Two runs of the
identical command disagreed: one reported GENE0002 as 16 all-`NA` rows, exit 0,
no warning. It is not nondeterministic arithmetic -- four instrumented scans
agree bit-for-bit on `sum(Score)`, `sum(Phi)` and `sum(Phi*Phi)` for all 2,527
cells -- it is an unchecked `arma::Mat::load()` of a chunk temporary file. Full
write-up in **LEDGER #57**; the tell was that the failing log has 716
`Analyzing chunks (` lines where every other run has 717.

**Two caveats before anyone quotes the 45%.** It is proportional to the size of
the annotation x max_MAF grid: a 1 x 1 run has one SKAT-O p-value per region
instead of fifteen, and the C++ fraction would dominate again. And it is
SKAT-O-specific -- the BURDEN row above is the same fixture with the same genes.

## SKAT-O in parallel: 1.30x, and the two things that nearly killed it (2026-09-04)

Built on the measurements above: the SKAT-O cells are independent, so a
persistent worker pool evaluates them while the region loop stays serial.
`--nThreads` drives it on a group run (it used to be discarded there).
`sim50k`, 210 genes, `y_binary`, two blocks with the arm order reversed between
them, load ~4.3:

| nThreads | real (A) | real (B) | mean | speedup | user | sys |
|---|---|---|---|---|---|---|
| 1 | 203.3 | 205.6 | 204.5 | -- | 170-173 | 26 |
| 4 | 157.5 | 155.9 | **156.7** | **1.30x** | 188-190 | 29-30 |
| 8 | 161.8 | 158.1 | 160.0 | 1.28x | 202-211 | 31 |

Replicate spread ~1%, so the effect is well clear of it. 8 workers buys nothing
over 4, which is what a 1.43x Amdahl ceiling minus per-task overhead looks like.
Every arm is BYTE-IDENTICAL to the one-thread run, region file and singleAssoc.

**The first version of this measured 0.97x at 4 workers and 1.06x at 8, with
user CPU UP 33 s and sys CPU up 40 s.** Two causes, and the interesting one is
not the obvious one:

- **The obvious suspect was wrong.** The design drives `region_postprocess`
  twice (collect the cells, evaluate them, serve the answers) rather than
  restructuring its 175-line triple-nested loop, so the loop's non-SKAT work
  runs twice. That looked like the problem. Measured in isolation -- collect,
  evaluate SERIALLY, serve -- the double drive costs **180.9 s user against
  180.6 s for the single drive**, i.e. nothing. Checking it before acting on it
  saved a pointless refactor.
- **The real cause was closure serialization.** The worker function was defined
  INSIDE `region_postprocess_parallel`, and R serializes a closure with its
  enclosing environment -- which held that frame's `...`, i.e. the whole
  `outList` including `VarMat`, shipped to a worker on EVERY task. ~13 ms per
  task against a ~40 ms task. Hoisting it to a top-level function (whose
  environment is the package namespace, which fork workers already have) is the
  entire fix: sys CPU 63 s -> 30 s and 0.97x -> 1.29x. **If a worker function
  ever wants a variable from its caller, pass it in the job, do not close over
  it.**

Also switched to `clusterApplyLB`. Cell cost spans three orders of magnitude
(0.5 ms to 540 ms) and the top 10% of cells hold 48% of the total, so static
round-robin would routinely hand one worker most of the work.

**What did NOT work, kept so it is not retried:** `mclapply` per region, i.e.
forking at each region rather than once. 55-90 ms of overhead against ~350 ms of
work, 1.12x at 4 workers and a net SLOWDOWN at 8.

## Phase 0 of the multi-phenotype work: how much of a marker's cost is SHARED?

**2026-09-04**, N = 50,000 (`SAIGE-sim/harness/work/sim50k`, 210,515 markers on
chr20, `--minMAC=20`, `--nThreads=1`, warm cache, 2 reps per trait). `uptime`
before the run: load 4.56 -- there WAS background load, which matters for the
absolute seconds and not for the ratio, which is what this measurement is for.

The instrument is `SAIGE_PROFILE_MARKER` (off by default; three
`std::chrono::steady_clock` buckets inside `testOneMarker`, serial arm only,
printed once per chunk). Buckets:

* **prepare** = read + decode + QC gate 1 + `imputeGenoAndFlip` + QC gate 2 +
  the two `conv_to` index copies. This is the work a multi-model scan does ONCE
  per sample-set group however many null models are tested.
* **test** = the score/SPA/ER calls and the fastTest re-test. Per model.
* **summary** = the row fields, including the binary case/control gather. Per
  model.

Whole-run totals (seconds of wall time inside `testOneMarker`):

| run | prepare `s` | test | summary | per-model `m` | `m/s` | P=5 factor | ceiling |
|---|---|---|---|---|---|---|---|
| y_binary rep1 | 63.27 | 24.39 | 6.50 | 30.89 | 0.488 | **2.16x** | 3.05x |
| y_binary rep2 | 58.43 | 21.88 | 6.00 | 27.88 | 0.477 | **2.18x** | 3.10x |
| y_quant rep1 | 63.82 | 5.88 | 0.02 | 5.89 | 0.092 | **3.74x** | 11.8x |
| y_quant rep2 | 74.32 | 5.88 | 0.02 | 5.90 | 0.079 | **3.86x** | 13.6x |

The factor is `P(s+m)/(s+Pm)` -- P serial runs against one P-model run. The
ceiling is `(s+m)/m`, what infinitely many phenotypes would approach.

Three things worth keeping:

1. **The prepare bucket is ~91% of the run's wall time on the quantitative arm
   and ~67% on the binary one.** `testOneMarker` accounts for 94.2 s of a
   103.9 s run at N=50k, so the flat/startup terms that dominate the N=1,000
   benches (`s2_marker6k`, ~90% startup) are gone at this scale. A speedup
   predicted from these buckets is therefore a prediction about the run.
2. **`summary` separates the two trait types completely**: 6.0-6.5 s for binary,
   0.02 s for quantitative. That whole difference is the two
   `GVec.elem(m_case_indices)` gathers, which are O(N) per marker per model --
   the term the "make the case/control summary O(MAC)" item in
   `docs/state/next.md` would remove, and the reason binary's ceiling is 3x
   where quantitative's is >11x.
3. The plan's a-priori arithmetic said ~2.5x (rare binary) / ~1.7x (common
   binary) and ~4-5x (quantitative) at P=5. Measured: 2.17x and 3.80x. Same
   ranking, same order of magnitude, binary at the pessimistic end of its
   predicted band -- which is expected, since `--minMAC=20` on this fixture is
   the COMMON-variant regime the 1.7x row describes, not the rare one.

Not measured here, and additive to the above: the `.bed` is read once instead of
P times (2.63 GB on this fixture) and there is one reader/mmap instead of P.
Working against it: the serial writer formats P rows per marker, so the Amdahl
serial fraction under `--nThreads` rises with P.

**Step 1.** `R/SAIGE_fitGLMM_fast.R` now prints the two timers it used to
compute and discard (the `system.time` around `setgenoFromInput`, and the one
around `glmmkin.ai_PCG_Rcpp_*`). At N=1,000 on the sparse-GRM golden the ingest
is `time: 2164 ms` of a ~6.5 s run, so sharing it across P=5 phenotypes predicts
~1.4x there; the ingest is O(M x N) against the sparse fit's O(nnz x iter), so
its share -- and the saving -- grows with the number of GRM markers.

---

## Phase 3: does multi-phenotype scanning deliver what Phase 0 predicted?

**2026-09-04**, N = 50,000 (`SAIGE-sim/harness/work/sim50k`), five traits
(`y_quant y_binary y_quant_pc y_binary_pc y_null`), `--nThreads=1`, warm cache,
paired arms run back to back in ONE batch. `/usr/bin/time -l`; **user CPU**, as
`bench_threads.sh` does.

### Step 2: 5 serial scans vs ONE manifest scan

210,515 markers on chr20, `--minMAC=20`, the fixture's own committed step-1
models, so both arms test exactly the same five.

| trait | serial user (s) | peak RSS |
|---|---|---|
| y_quant | 64.74 | 2761 MB |
| y_binary | 107.81 | 2174 MB |
| y_quant_pc | 81.73 | 2782 MB |
| y_binary_pc | 92.66 | 2779 MB |
| y_null | 65.34 | 2783 MB |
| **serial total** | **412.28** | (one at a time) |
| **one manifest run** | **137.86** | **2843 MB** |

**Speedup 2.99x on user CPU**, 3.28x on wall clock (481.9 s -> 146.9 s), and
**all five output files BYTE-IDENTICAL to their serial counterparts.**

Peak resident memory for the five-model run (2843 MB) is essentially that of a
SINGLE serial run (2174-2783 MB): the five null models are small beside the
reader's buffers, so the whole scan costs one process's memory instead of five.

**Against the Phase 0 prediction.** Phase 0 measured, per marker, shared
`s ~ 65 s` and per-model `m ~ 29 s` (binary) / `5.9 s` (quantitative) over the
whole scan. For this 2-binary / 3-quantitative mix that predicts
`(5s + 2m_bin + 3m_qnt) / (s + 2m_bin + 3m_qnt)` = `400.7 / 140.7` = **2.85x**.
Measured **2.99x**, 5% above the prediction -- the residual being the per-process
startup and model load that the serial arm pays five times and the buckets never
counted. Prediction and measurement agree; nothing here needs explaining away.

### Step 1: 5 serial fits vs ONE `--phenoCol=` run

Sparse GRM (pedigree truth `.mtx`), `--useSparseGRMtoFitNULL=TRUE`,
`--useSparseGRMforVarRatio=TRUE`, LOCO, VR markers from `data_for_vr`.

| | user CPU | peak RSS |
|---|---|---|
| 5 serial fits | 61.60 s (11.7-13.6 each) | ~390-417 MB each |
| one multi run | **37.42 s** | 440 MB |

**Speedup 1.65x**, one ingest instead of five (asserted from the log's
`setgeno mark2` count). The plan predicted ~1.4x from the N=1,000 ingest share;
1.65x here because this fixture's VR marker set is larger relative to the fit.
**All five `.rda` files are BYTE-IDENTICAL** across the arms.

**The `varianceRatio.txt` files are NOT byte-identical for phenotypes 2-5, and
that is NOT attributable to this feature.** Measured both magnitudes:

| comparison | max relative difference |
|---|---|
| SAME command, two batches, standalone (pre-existing run-to-run noise) | **1.0e-7** |
| multi vs serial, within ONE batch (the feature) | **1.1e-7** |

Same magnitude, and in BOTH cases confined entirely to **MAC category 2**
(category 1 is exactly 0.000e+00 everywhere).

**CORRECTION to the first version of this section.** It said "three standalone
`y_quant` runs in one batch reproduce each other exactly, so the run is
self-reproducing within a batch", and inferred that the noise was a
between-batch effect. **That inference was drawn from n = 3 and is wrong.**
Eleven runs of ONE identical command (same binary, same machine, back to back,
plus a rebuild with no source change, plus two under deliberate 8-way load)
produced **SIX DISTINCT `varianceRatio.txt` files** -- and the within-batch runs
disagreed among themselves (`a1 == a2`, `a1 != a3, a4, a5`, `a1 == a6`). There
is no "stable within a batch": the command is simply not reproducible. This is
recorded as **LEDGER #56**. It makes the conclusion below STRONGER, not weaker,
but the reasoning that reached it was faulty and is replaced here.

The measurement, over those 11 identical runs:

| row | distinct values | spread (relative) |
|---|---|---|
| `sparse 1`, `null 1`, `null_noXadj 1` | **1** | **0.000e+00** |
| `sparse 2` | 6 | 3.8e-12 |
| `null 2` | 6 | 1.1e-7 |
| `null_noXadj 2` | 6 | 1.5e-7 |

All ELEVEN `.rda` files are byte-identical. So the null model is fully
reproducible and only the HIGH-MAC variance-ratio category is not.

The multi-vs-serial difference (1.1e-7) sits inside that 1.5e-7 spread, in the
same rows, and the third control still stands: forcing the multi run to
re-ingest per phenotype (shared PROCESS, no shared ingest) does not remove the
difference, it merely changes WHICH phenotypes differ -- which is what a draw
from a noisy distribution does and a systematic effect does not.

LEDGER #56 records the defect and what is and is not known about its mechanism.
`[ASSUMPTION: the mechanism is one of the RcppParallel `parallelReduce` sites in
`src/SAIGE_fitGLMM_fast.cpp` (`parallelCrossProd` at :1773 for the `null` rows,
`CorssProd_usingSubMarker` at :4031 for `sparse`) or the `tbb::concurrent_vector`
insertion order in `indicesRelatedSamples` (:231, :2101). Basis: those are the
only concurrency on this path and the magnitudes match a float reduction whose
order changed. NOT established -- `setThreadOptions(1)` is called
unconditionally at `R/SAIGE_fitGLMM_fast.R:824` and `nThreads` is forced to 1
here, which SHOULD make all of them deterministic. That contradiction is the
thing to chase.]` At N = 1,000 with `--useSparseGRMforVarRatio=FALSE` -- the
test suite's fixture -- the variance ratio files ARE byte-identical, in both
phenotype orders, which is what `tests/testthat/test-step1-multipheno.R`
gates.

**What is safe to claim:** the step-1 null MODEL (`.rda`) -- everything step 2
consumes -- is byte-identical, unconditionally, at both N. The variance ratio is
byte-identical wherever the configuration is reproducible at all, and otherwise
agrees to the precision that configuration has.

### Honest bottom line

Step 2 is where the win is: **~3x CPU at P = 5**, one genotype read instead of
five, one process's memory instead of five, and byte-identical output. Step 1
gains **1.65x**, because sharing the ingest saves only the ingest.

---

## strip vs upstream at N = 50,000, whole pipeline, paired (2026-09-05)

Full write-up: `harness/sim/upstream_compare/RESULTS.md` (drivers, plots and every
run log beside it). Upstream `main` 365deeff built in `SAIGE-main` with its own pixi
env pinned to the SAME conda builds as strip's (R 4.4.3 `h46d4a84_2`, OpenBLAS 0.3.30
`openmp_h30af337_0`, same clang) -- unpinned it resolved R 4.5.2 / OpenBLAS 0.3.34,
which would have been a confounder and which also fails to link (`-lheapt_w`).
Both arms on `sim50k`, frozen step-1 models, `--nThreads=1`, 3 paired reps with the
arm order alternated per rep, `CPU_Speed_Limit = 100`; every log records the library
loaded and the `SAIGE.so` md5; every arm reproduces itself byte for byte over the reps.

**User CPU, upstream / strip** (paired t on 3 reps):

| workload | upstream | strip | ratio | 95% CI on diff |
|---|---|---|---|---|
| step 1 sparse fit + VR, binary / quant | 9.03 / 8.95 s | 6.58 / 6.68 s | **1.37x / 1.34x** | [2.4, 2.5] / [2.2, 2.4] s |
| single-variant scan, binary | 236.3 | 252.1 ± 14.2 | **0.94x** (n.s., p = 0.15) | [-45.7, +14.1] s |
| single-variant scan, quant | 80.2 | 62.5 | **1.28x** | [15.0, 20.5] s |
| gene tests (SKAT-O, 15 cells/gene), binary | 222.5 | 167.5 | **1.33x** | [30.7, 79.3] s |
| gene tests, quant | 201.7 | 147.3 | **1.37x** | [13.2, 95.6] s |

The quant scan reproduces the 09-01 reader comparison (1.37x bed then, 1.27x on
`work_s` now; different day, so not a regression claim). **The binary scan is the
new fact: strip is not faster there, and its rep spread is 4x upstream's.** Nothing
in Phase 4 measured a binary scan at this N. The binary path adds SPA, ER, Firth on
10% of markers (`--pCutoffforFirth=0.1`) and the O(N) case/control gather; at 4
threads its user CPU rises 48% against the quant scan's 27%, so the extra is in the
per-marker kernels, not the writer. Profile it before touching anything else.

Strip only, `--nThreads=4`, wall clock vs its own 1-thread mean, all outputs
byte-identical to 1 thread: single-variant **2.35x** (binary) / **2.72x** (quant);
gene tests 1.21x / 1.59x (threads go to SKAT-O post-processing; NOTES above measured
1.30x on this fixture at 4 workers). Against upstream on wall clock: 2.2x / 3.5x
single, 1.6x / 2.2x gene.

**maxRSS: strip 2.3-2.8 GB vs upstream 0.35-0.38 GB on the scan is the mmap'd 2.6 GB
`.bed`'s resident pages**, not allocation -- `SAIGE_PLINK_NO_MMAP=1` peaks at 312 MB
with byte-identical output. Group runs: 1.0 GB vs 0.62 GB.

**P-values.** With the same frozen model, strip's step 2 is bit-identical to upstream's
on every gene-test cell (3,338 x 2 traits), every single-in-group row both print, and
every single-variant row with `AF_Allele2 <= 0.5`; it differs on every AF > 0.5 row
(LEDGER #28, upstream's), tests 8 more markers (#29), and prints 965 / 305 single-in-
group rows upstream loses in its 22 / 7 escalated regions (#50). #28 is not small:
hand recomputation from the raw `.bed` reproduces strip's `var` = sum((g-gbar)^2) x VR
on both allele sides and shows upstream's AF > 0.5 `var` 2-7x too large; on `y_null`
upstream reports 0.7% of AF > 0.5 markers below p = 0.05 (5% expected; strip 4.0%).
Flip invariance checked directly (alt-first vs ref-first, 667 markers, 124 AF > 0.5:
`var` and p bit-identical, `BETA` sign-flipped). Causal-gene ranks identical between
arms on both binary traits.

With each arm's OWN step-1 model the remaining movement is all step 1 and all
reconciled row by row (a chi-square rescale by the arms' reported `var`, or a match to
the arms' variance-ratio rows): #53 moves upstream's category-1 ratio 0.14%, and, on
quantitative traits, **LEDGER #97, found here**: the #65 rescale multiplies the
`sparse` variance-ratio rows by 1/s^2 although that ratio is scale-free (both sides go
through Sigma^-1), so every HEAD quantitative model fitted with
`--useSparseGRMforVarRatio=TRUE` has `sparse` rows 11-12% high. It reaches the region
path and the escalated sparse re-test: quantitative gene-test `Pvalue` moves by up to
0.6 in -log10 on this fixture, `SE_Burden^2` strip/upstream = s^2 exactly on the cells
that use the row. No committed quantitative golden has `sparse` rows, which is why
the tier stayed green. One-line fix, test named in the row.

## The binary single-variant scan at N = 50,000: profile, and what it bought (2026-09-05)

Why: the strip-vs-upstream table above has the binary scan at 0.94x while the
quantitative one is 1.28x, and the binary scan's user CPU grew 48% under 4 threads
against 27% for the quantitative one. Same fixture (`sim50k`, frozen `step1_y_binary`,
`--minMAC=20 --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_fastTest=TRUE`),
macOS `sample` at 1 ms over the whole R process (`pgrep -f step2_SPAtests.R`, comm
== R), 245,658 samples at 1 thread; a quantitative run as the control. Load was
3.7-4.5 throughout (Chrome, WindowServer), so shares and paired ratios are what is
trusted here, not absolute seconds.

**Where the CPU goes, 1 thread, binary** (inclusive; children summed by symbol):

| frame | share | per marker |
|---|---|---|
| `fast_logistf_fit_simple` (Firth) | **56.1%** | 17.5 ms x 7,891 fits, 9 Newton iterations each |
| `PlinkClass::getOneMarker` | 15.1% | 0.18 ms x 210,515 decoded |
| SPA: dense `SPA()` 6.4% + `SPA_fast()` 3.8% | 10.2% | 15.6 ms / 3.2 ms per corrected marker |
| `imputeGenoAndFlip` | 5.6% | |
| score tests + fastTest re-test + `getadjGFast` | 7.6% | |
| case/control gather | ~1.5% | |
| `writeOutfile_single` | 0.09% | |

The quantitative control is decode 60%, impute 22%, `getMarkerPval` 13%; its decode
and impute sample counts match the binary run's to 0.3%, so the entire binary
excess is `getMarkerPval`. Inside Firth: `exp` 16%, `qr_econ` 12%, three `dgemv`
15%, Fisher `syrk` 8%, column fills and zeroing 16%, allocator 4%. A per-chunk
regression of the existing timing logs (user CPU on markers tested and Firth count,
21 chunks x 3 reps) gives **33 ms per Firth-eligible marker on BOTH arms** (upstream
R^2 = 0.96) and ~0 per other marker, so this is a shared cost, not a strip regression.

**Why the binary comparison read 0.94x: it was not like-for-like.** On AF <= 0.5 both
arms fire Firth on 6,258 markers and SPA on 3,192, identically. On AF > 0.5 strip
fires Firth on 1,633 and SPA on 730 against upstream's 569 and 276 -- upstream's
LEDGER #28 deflates its own p-values there (0.7% below 0.05 on `y_null` against 5%
expected). Strip therefore did ~1,000 more Firth fits and ~450 more SPA corrections,
roughly 22 s of work upstream skips; like-for-like the arms are at parity, and the
Phase 4 gains that show on `y_quant` are diluted because 71% of the binary CPU was
in Firth/SPA code strip had never touched.

**Four threads.** Pure-compute leaves inflate uniformly 1.3-1.4x (`exp` 1.43x, decode
1.27x, impute 1.36x): this 4-core i5's turbo and hyperthreading, and the whole of the
quantitative scan's +27%. Above that baseline sat only allocation paths:
`mach_vm_deallocate` 6.8x, `bzero` 2.9x, Armadillo temporary constructors 2.2-2.6x,
`memmove` 2.25x, sys time 19 s -> 61 s. Allocations over 127 KB go through mach VM on
macOS and serialise on the task's map lock: the Firth fitter built ~10 n-length
temporaries per iteration and each SPA kernel evaluation 3-4, i.e. ~90 map/unmap
pairs per Firth marker and ~80 per dense-SPA marker.

**Done, byte-neutral (items 1 and 2): per-worker kernel workspace.**
`SAIGE::KernelWorkspace` on `Marker::MarkerScratch` (SAIGE_test.hpp), handed to
`getMarkerPval` through a trailing pointer that defaults to null (region-path callers
untouched; they get one lazily built instance per call). The Firth fitter keeps its
n x 2 design, weighted design, Q and the five n-vectors across iterations and
markers; the `Korg/K1/K2` kernels write into strict views of four shared buffers
(`SPAWorkspace`); every expression is the original verbatim. Gate: 1- and 4-thread
50k outputs `cmp`-identical to the pre-change run (7,891 Firth + 3,922 SPA rows
exercised), fast tier 0 failures, standard tier 0 failures / 2,718 passes.

| run | real s | user s | sys s |
|---|---|---|---|
| before, 1 thread (3 reps, earlier today) | 255-285 | 240-268 | 13-15 |
| after, 1 thread | 214 | **205** | **5.0** |
| before, 4 threads (2 reps) | 111-118 | 363-382 | 49-55 |
| after, 4 threads | **83.5** | 297 | **10.9** |

Wall speedup at 4 threads is now 2.56x (was 2.35x); the 4-thread sys time is a fifth
of what it was. Microbenchmark of the fitter alone (Rcpp, n = 50,000, 6 markers x 20
reps): 11.1 -> 8.5 ms per fit, bit-identical BETA/SE.

**Benchmarked behind flags, NOT default (items 5 and 6), same fixture, 1 thread, one
batch, base first (load 6.0 at its start, 3.7-5.1 for the rest):**

| flags | user s | vs base | output vs base |
|---|---|---|---|
| default (`qr`, cold start) | 235.5 | -- | byte-identical to pre-change |
| `--firth_solver=closed` | 170.8 | **-27%** | **byte-identical**, all 70,709 rows, 7,891/7,891 converged |
| `--firth_warm_start=TRUE` | 172.4 | -27% | 137 rows differ; 118 fits NOT converged (7,773/7,891); BETA moves up to 4x, sign flips, MAC 20-88 |
| both | 138.8 | -41% | same 137 rows as warm start alone |

`closed` replaces the QR of the n x 2 weighted design with the 2 x 2 Fisher matrix
and `h_i = w_i x_i' F^-1 x_i` in closed form; the microbenchmark put its rounding
difference at <= 8e-13 relative, which never reaches the 6 printed digits here. The
warm start (intercept at logit(mean(y)) = -4.60, the offset being centred) is
harmful: from the null intercept the weights are ~0.01 instead of 0.25, the first
slope step on a MAC 20-90 marker saturates at `maxstep` and the iteration lands in
the wrong basin or oscillates to `maxit`. Recommendation: keep `closed` as an
opt-in (or promote it once the owner accepts a flag-gated path), delete the warm
start. Neither flag can move a golden: every golden runs the defaults.

**Done, byte-neutral (item 3): QC gate 1 from the code tallies, before any decode.**
`PlinkClass::countOneMarker` (src/PLINK.cpp) tallies a marker's four 2-bit codes
and returns the metadata and the gate quantities through the SAME arithmetic the
full read uses (`summariseCounts`, factored out of `getOneMarker` so the two cannot
drift); `prepareOneMarker` (src/Main.cpp) asks for it first through
`GenoReader::countOneMarker` (default: unsupported, PGEN unchanged) and only reads
the marker in full when it passes. Skipped when the male X non-PAR rewrite is on,
since that runs after the read. When the analysed samples are the file's samples
in file order the tally is one table lookup and one add per byte (four 8-bit
lanes, unpacked every 63 bytes); otherwise it is the per-sample loop without the
double store, the missing index or the genotype map. **The first version did
four masked increments per byte and bought nothing** (quant 61.4 s against
62-65 s): a scalar per-sample loop is the decode's cost whatever it stores.

| run (1 thread, 50k) | before | after |
|---|---|---|
| quantitative, user / sys | 60-65 s / 3 s | **45.3 s** / 3.3 s |
| quantitative, last chunk's `prepare` bucket (515 markers) | 0.180 s | 0.095 s |
| binary, user | 197-205 s (items 1+2 build) | 203 s (noise; decode was 15%) |

210,515 markers tallied, 70,709 expanded. Both branches are gated: the 50k
fixture's model is in fam order (table branch) and every bundled golden runs
1,000 of 10,000 samples (subset branch). Outputs `cmp`-identical for both traits;
`prepush.sh` green in both modes (byte 2,718 / numeric 2,918 passes).

**Decided 2026-09-06: `--firth_solver=closed` is the default; the warm start is
deleted.** P31 is byte-identical under the new default, so no golden moved; the 50k
binary output is byte-identical to the pre-session run; `prepush.sh` green in both
modes. `test-step2-firth-solver.R` pins `qr` against P31 at 1e-10 (parsed) and the
refusal of an unknown solver. The adoption run's timings (542 s wall for the scan,
40 min for the standard tier) are NOT benchmarks: the machine was at load 100 from
VS Code's C++ indexer, Docker and a second session; the -27% above is the number.

## The serial step-2 writer does NOT cap `--nThreads`, at any model count (2026-09-06)

The `next.md` item "format step-2 output rows on the workers" rested on an
Amdahl argument: P models format P rows per marker on the one writer thread, so
the serial fraction grows with P. Measured at HEAD `b11e6b07` (library rebuilt
from it, md5 `41844865`), `sim50k`, `--minMAC=20`, the `step2_single` flags,
`bash harness/bench_writer.sh w1 2` (driver + `writer_profile.py` +
`writer_analyse.py` beside this file; rows in `writer.tsv`, trees in
`profiles/writer_w1.*`). Decision rule fixed before the run: implement if
w(5 models, 8 threads) >= 0.15, close if < 0.10.

**w = writer / main-thread time inside `mainMarkerInCPP`** (1 ms `sample`):

| cell | rows | T=1 | T=8 | us per row |
|---|---|---|---|---|
| 1 binary | 70,709 | 0.002 | 0.005 | 2.7 |
| 1 quant | 70,709 | 0.004 | 0.014 | 2.0 |
| 3 quant (worst case: smallest per-model kernel) | 212,127 | 0.008 | 0.021 | 1.9-2.0 |
| 5 mixed (2 binary + 3 quant) | 353,545 | 0.003 | 0.009 | 2.4 |

**Verdict: CLOSED.** The best case for the item is 1/(1-0.021) = **1.02x**.
The writer's mechanism is real -- 74% of its samples are number formatting
(`__put_num`/`num_put`/`snprintf`/`__dtoa`), 15% the `filebuf` write -- its
magnitude is 2-3 us per row against 200-2,000 us of marker work, and it does
not grow with P faster than the work does.

Instrument controls, all before the verdict was read: writer samples are
invariant with T within 3% on every cell (194/191, 142/143, 403/415, 843/849)
and scale with ROWS not models (per-row cost 2.4 us at P=1 vs 2.1 us at P>1;
the binary-vs-quant spread 2.7 vs 2.0 is row width, binary rows carry 4-8
more columns). `sample` attaches after the R prologue and its cadence slows
under 8 busy threads (main-thread samples 57-90% of timed wall), which
affects absolute seconds and not the within-profile ratios that w is.

**Where the lost thread scaling actually is: CPU inflation, not a serial
section.** Wall time, mean of 2 reps, T = 1/4/8; user CPU beside it:

| cell | wall 1 | wall 4 | wall 8 | speedup 4 / 8 | user 1 -> 4 -> 8 |
|---|---|---|---|---|---|
| 1 binary | 146.6 | 69.6 | 70.4 | 2.10x / 2.08x | 141 -> 242 -> 368 |
| 1 quant | 44.4 | 17.3 | 14.5 | 2.56x / 3.07x | 41 -> 52 -> 67 |
| 3 quant | 66.0 | 28.8 | 26.9 | 2.29x / 2.46x | 60 -> 84 -> 120 |
| 5 mixed | 274.4 | 111.1 | 118.1 | 2.47x / 2.32x | 266 -> 391 -> 724 |

An Amdahl fit reads these as s = 0.19-0.41, but the profiles refuse that
reading: at T=8 the main thread spends 91-95% of its time doing marker work as
the eighth worker, 2-8% waiting on the pool and <= 1.2% in writer + post-passes.
The same work costs 1.3-2.7x more user CPU at T=8 than at T=1 -- the
allocation churn and hyperthreading (8 logical on 4 physical cores) that
`docs/state/handoff-2026-09-05-binary-profile.md` located -- and T/inflation
predicts the T=4 speedups to within ~10% (binary 2.34 vs 2.10, 5-mixed 2.72 vs
2.47), the residual being two desktop processes (VS Code's C++ indexer, a
Chrome renderer, ~30% of a core each) that also cost binary rep 2 up to 18%
against rep 1. **On this machine T=8 is no faster than T=4 for anything with a
binary model**, and the P=1 "serial fraction" of 0.16-0.23 from the upstream
comparison is this inflation, not serial code.

Side findings, same runs: after the 09-05/06 kernel work the binary scan is
146.6 s (252 s at `fa0ec961`) and the quant scan 44.4 s (62.5 s); the 5-model
manifest now buys **1.55x** over five serial runs (426 s -> 274 s; Phase 3
measured 2.99x) and the 3-quant manifest 2.02x, because the shared prepare has
less per-model work left to amortise. All outputs byte-identical across T and
against the single-model runs (driver-enforced `cmp` on every pair). Swap
(3 GB) went 1.53 -> 2.34 GB peak with >= 730 MB free and 57-64% system memory
free throughout; the scan's 2.5-2.9 GB RSS is the mmap'd `.bed`, and the swap
growth is idle desktop pages, not pressure.

### Addendum, same day: the remaining thread inflation is the hardware floor, not churn

The paragraph above inferred ~1.25x of removable allocation churn on the binary
path at 4 threads by subtracting the 1.3-1.4x turbo/hyperthreading floor from a
1.71x inflation measured under desktop load. Re-measured on a quiet machine
(`bench_writer.sh w2`, load 2.25): the 4-thread binary scan is **56.6 s wall,
196 s user against 141 s at 1 thread = 1.39x**, i.e. the floor with nothing
left over. The 4-thread profile agrees: allocator + lock frames are 3.1% of
self time (1.4% at 1 thread), pool waiting 4.3%, and every compute leaf's
share FALLS (`exp` 0.89, Firth 0.86, impute 0.77, decode 0.80) while the
memory-bound BLAS `ddot` rises 1.36x -- the signature of memory bandwidth and
shared cores, not of malloc. At 8 threads the allocator reaches 5-7.5%
(`mach_vm_deallocate` 3.7% on the quantitative path: the N-length doubles
freed per marker) and waiting 4.4-6.9%, but 8 threads is already no faster
than 4 here. **Closed: a per-worker workspace for the score-test/impute
temporaries would buy at most ~2% at 4 threads on this machine.** The
`--nThreads` help now says to use the physical core count.
