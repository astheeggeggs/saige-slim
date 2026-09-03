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
