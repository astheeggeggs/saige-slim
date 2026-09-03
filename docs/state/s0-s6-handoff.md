# Step-2 threading sequence: S0-S6 landed (handoff, 2026-09-02)

Branch `strip`, **fully pushed to `slim/strip`**; everything through
`0c42575b` is green on CI and S6 (`7814e24f`) was pushed at the end of the
session -- check its run before building on it. Working tree clean.
Commits, in order:

| commit | what |
|---|---|
| `b87fcaba` | **S0** — seven marker-path goldens (P27-P33) + step-1 golden `binary_sparseGRM_vr`; LEDGER #44; `rel_err` non-finite fix |
| `8959b5bd` | **F8** — PCG allocation hygiene (the two items C3 excluded for want of coverage) |
| `fcb3feeb` | **S2** — ER seeds the R RNG only on markers that resample (closes F4) |
| `357fb6df` | **S1** — POD returns for the SPA kernels; LEDGER #45, #46 |
| `0f828a14` | **S3** — `TestContext` + const marker path; LEDGER #32 **fixed** |
| `c10181fb` | this handoff record |
| `1c7690e1` | CI — plink-ng headers guarded separately from the archive |
| `f1d03251` | comparator — parse the non-finite spellings, not just their text |
| `456f47d6` | docs — S4 anchored on names, current line numbers |
| `1c15e3b7` | **S4** — `MarkerResult` + `MarkerScratch` + `testOneMarker`; the marker path's `vector<bool>` is gone |
| `feb19ded` | comparator — the covariate-collinear rows excluded from numeric compare (CI forced it, again) |
| `bd752aec` | **S5** — reader clone with shared ownership; `getOneMarker` const; the pread fallback tested for the first time |
| `0c42575b` | docs — the fork-path retirement criterion; a stale #32 comment corrected |
| `7814e24f` | **S6** — in-process marker threading, default OFF; F4's inventory corrected |

The authoritative plan is the status block at the top of
`docs/design/step2-threading-verdict.md`, which is current as of this commit.
`LEDGER.md`, `CLAUDE.md` and `harness/bench/NOTES.md` are all updated in the
commits that changed the facts they record; nothing is carried only in a commit
message.

## State of the gates

- Full `make_goldens.sh` + `compare_goldens.sh` cycle run on each of the five
  code commits: every step-0 and step-2 golden byte-identical (54 files at S4),
  all five step-1 goldens at `max_rel_err 0.000e+00`.
- Standard tier green in **both** compare modes at each of them.
- S4 and S5 are byte-identical to the goldens by construction, so each carries
  its own separate measurement -- see "Things worth knowing" below.
- The installed package in `.pixi` was built from this exact tree (clean
  objects). No rebuild needed to pick the work up.

## What pushing cost, and what it was worth

CI was red twice before it was green, and NEITHER failure was in the work S0-S3
did. Both are worth knowing about because both were latent and would have bitten
whoever pushed next.

1. **The build step, `pgenlib_read.h: No such file or directory`.** The workflow
   cached only `.pixi/envs/default/lib/libplink2_includes.a` while
   `harness/build.sh` guarded its whole fetch-and-compile block on that one
   file -- but `plink-ng/` is gitignored, not cached, and `src/Makevars`
   compiles everything with `-I../plink-ng/2.0/include`. So a cache HIT restored
   the archive, skipped the fetch, and lost the headers: correct on a cold
   cache, broken on every warm one. The 2026-09-01 green run is exactly what
   populated the cache; ours was the first hit. Fixed in `1c7690e1` by guarding
   the header tree and the archive independently (a warm run now re-fetches the
   source in seconds and still skips the compile) and by caching both.
2. **`P27` failed on Linux with "string column differs [SE]".** A comparator
   gap, not a numeric one: `.is_numeric_col` asked `!anyNA(as.numeric(x))`, and
   `as.numeric("nan")` returns NaN, which `anyNA()` cannot tell from a parse
   failure -- so any column holding a NaN was compared as TEXT. And text is what
   is not portable: glibc prints a sign-bit-set NaN as `-nan` where macOS prints
   `nan`. P27/P30's rs20 and rs22 carry exactly that (LEDGER #44). Fixed in
   `f1d03251`.

The lesson for the next stage is the general one: **`SAIGE_COMPARE_MODE=numeric`
passing locally is not the same check CI makes**, because on the golden-
generating platform the numeric path can be reached with text that happens to
match. Push early rather than at the end of a stage.

3. **And then the same two rows again, one layer down.** With the spelling
   parsed, the comparison reaches the VALUES on those rows, and CI went red on
   `456f47d6` -- a **docs-only** commit, which is what proves it is the goldens'
   property and not any code change: `rel_err` 1.2e-01 [BETA] and 1.6e-01
   [Tstat] on P27/P30. `rs20`/`rs22` ARE the model's covariates, so on any
   covariate-adjusting path `|Tstat|` is ~3e-15 and `var` is 0 or ~3e-14 against
   [0.021, 14.5] and [0.020, 39.9] for the other 68 markers -- thirteen orders
   of magnitude of pure cancellation noise. `BETA = Tstat/var` is then a RATIO
   of two noise quantities, O(0.1) and arbitrary, which is why an absolute
   tolerance floor would not have caught it: **the number is not small, it is
   meaningless**. Note the gate is INTERMITTENT, not deterministically red --
   the very next run, on `1c15e3b7`, passed with the same goldens and the same
   test code. That is the expected signature: low bits set by summation order,
   so a different runner microarchitecture flips them. Handled by excluding
   those two markers from the NUMERIC comparison of P27/P30 only, with controls
   in `test-meta-comparators.R` proving the exclusion is narrow (dropping a
   different row still fails; a relabelled or absent marker fails). **Fixing
   LEDGER #44 retires the exclusion** -- emit `NA` under the predicate that
   already forces `p = 1` -- and that is now the main argument for doing it.

## What S4 did

`src/MarkerLoop.hpp` is new and holds the two types. `testOneMarker` is the loop
body; `mainMarkerInCPP` is now the slot vector, the scratch, the progress line,
the `break`, the Firth tally and the writer call; `writeOutfile_single` takes a
`const std::vector<Marker::MarkerResult> &`. Both of the things S4 was supposed
to fix are fixed structurally rather than by convention:
- `std::vector<bool> isSPAConvergeVec` is gone. **No `std::vector<bool>` remains
  anywhere on the marker path** -- the two left in `Main.cpp` are the region
  path, which is not threaded -- so two markers can no longer share a word.
- The loop-hoisted buffers are a named unit a worker owns, instead of six
  variables in `mainMarkerInCPP`'s frame that happen to be reused.

The Firth counters were the loop's only cross-marker accumulators; they are now
per-marker flags on the result, summed after the scan. `is_FirthConverge` is
still only read inside `if(is_Firth)`, which is what keeps it off the
uninitialised path -- the Firth fit is the only thing that writes it.

One dead store was dropped rather than carried into `MarkerResult`: the loop
built a `CHR:POS:REF:ALT` string into an `infoVec` that the writer was never
handed.

## What S5 did

`PlinkReader` holds a `std::shared_ptr<PLINK::PlinkClass>`; clones share it and
`close()` releases only that handle, so the map dies with the LAST owner. The
thing that made it sound is that **`PlinkClass::getOneMarker` is now `const`** --
and NOT by marking the fallback's state `mutable`, which would have made const
meaningless. The shared scratch buffer `m_OneMarkerG4` is gone and the
fallback's `fseek`+`fread` is one `pread`, which takes the offset as an argument
instead of mutating the file description's. Neither path writes a member.

Two things surfaced that nobody had noticed:
- **Step-2's `closeGenoFile` is never called from R.** The map and the `FILE*`
  survived to process exit. `~PlinkClass` releases them now, and the class is
  non-copyable so the destructor cannot double-close.
- **The pread fallback had never been executed by anything** -- no test, no
  golden, no platform. It was reachable only by an mmap failure, which does not
  happen for a small local file. `SAIGE_PLINK_NO_MMAP=1` forces it; it now has a
  unit test over decoded vectors and an end-to-end run reproducing golden P33.

`isForkSafe()` still answers `isMemMapped()`, but the reason changed: pread
removed the shared-offset hazard, so the fallback is fork-safe in fact and the
`false` is conservatism about untested scale. Do not read it as a live hazard.

## What S6 did, and the thing it found

In-process threading of the marker loop, **default OFF**: `nThreadsMarker` on
`SPAGMMATtest` (no CLI flag; S7 routes `--nThreads` here). At 1 the loop is
literally serial and touches no threading machinery; above 1 it is
`parallelFor` with an explicit `numThreads`, a checkout pool of one
`MarkerScratch` + reader clone per thread, and four serial post-passes that
restore what workers cannot decide locally: rethrow the LOWEST-indexed error,
truncate at the earliest short read, re-run deferred markers, emit progress
lines in order.

**F4's inventory was incomplete, and the missing entry is the binding one.**
F4 named `SL_setseed`, on the ER resampling arms only, which S2 closed.
`Rf_lchoose` is called on EVERY ER marker via `HyperGeo::lCombinations`, and
cannot be called off the main thread at all -- measured in a 20-line program
with no SAIGE in it, aborting for every argument pair tried, while
`lgammafn`/`pnorm`/`qnorm`/`pchisq` are all fine. So the whole ER branch is
main-thread-only. **The cost is real**: at the default `max_MAC_for_ER = 4`
every binary marker with MAC <= 4 runs serially, which on a rare-variant scan
is most of them; quantitative traits thread fully. **The fork path does not
have this limit** -- each child is its own process with its own R main thread --
so the two mechanisms differ in capability, not just plumbing. That belongs in
the retirement comparison, and it is now written there.

## The two deferred S6 measurements are DONE (2026-09-02, later the same day)

Neither changed any code; both are recorded where the facts live
(`harness/bench/NOTES.md` for the numbers, the verdict doc's status block for
the decision, `CLAUDE.md` for the standing rules). New in the tree:
`harness/tsan.sh` and `harness/bench_threads.sh`, plus a `SAIGE_TSAN_ANNOTATIONS`
block in `src/Main.cpp` that is compiled out of every normal build.

**TSAN is clean, and the controls are what make that a statement.** The raw run
was not clean -- 13-35 reports per run -- but every one was main-thread-vs-worker
across the `parallelFor` join, none worker-vs-worker, because `libtbb` is
uninstrumented and TSAN cannot see the join. Declaring that edge takes them all
to zero across three workloads at 4 and 8 threads, output byte-identical to
serial; and with the SAME annotation still active, the shared-`MarkerScratch`
defect fires in every threaded cell (205-1130 races, worker-vs-worker, wrong
answer in all six) while its own 1-thread arm stays at zero.
`harness/tsan.sh plant` is that positive control. Four macOS-specific
requirements are baked into the script because each one fails by reporting *no
races* rather than by failing loudly -- see the header comment.

**Speed and memory are measured, at `CPU_Speed_Limit = 100`.** In-process at 4
threads: 2.07x / 2.53x / 2.09x on `work_s` for under 1 MB of maxRSS. Against
fork through a common serial baseline, the two mechanisms are **indistinguishable
on wall clock** (p = 0.14 / 0.53 / 0.49) and differ only in memory, where fork
costs +47 to +90 MB.

**The measurement did not settle the retirement, and that is the result.** Only
5,292 of 128,868 markers (4.1%) on the bundled fixture have MAC <= 4, so the ER
serialisation in-process pays and fork does not is a 4% tax here; the
rare-variant regime the criterion warns about is unreachable with what this repo
ships. The speedup bound is `1/(f + (1-f)/T)` in the deferred fraction f: 3.6x at
f = 0.041, 1.6x at f = 0.5, 1.08x at f = 0.9. Do not read the null p-values as
"the mechanisms are equivalent".

## The fixture was built, and it settles it (same session)

`harness/make_rare_fixture.py` varies the one thing the bundled fixture cannot:
the fraction f of ER-eligible markers (binary, MAC <= 4). Real genotype blocks,
uniformly interleaved, positions renumbered -- both of those matter, see the
docstring. Across f = 0 -> 1 at 4 threads, **fork's speedup is FLAT (1.16-1.22)
while in-process decays 1.248 -> 0.954**: at f = 1 in-process is 13% SLOWER than
serial, because a deferred marker is decoded on a worker, discarded, then re-read
and re-run on the main thread. A one-parameter model fitted on f = 0 alone
predicts the whole in-process curve and misses only there, in the direction
named in advance. **The two mechanisms are not substitutable**, and f = 0.041 is
exactly where they look identical.

**And the deferral may not be necessary at all.** `harness/probes/lchoose_offthread.cpp`:
`Rf_lchoose` off the main thread aborts only because of R's C-stack check, and
with `R_CStackLimit` neutralised it returns **bit-identical** values. Caveats in
the verdict doc -- it does not prove the rest of the ER machinery is thread-safe,
and `harness/tsan.sh` has never run with the deferral disabled.

## Steps A and B landed: the ER branch threads, and in-process now dominates

Option (ii) was taken and it worked, without the risky part. **Step A**
(`415a3db3`) computes `lchoose` from `Rf_lbeta` the way R does
(`HyperGeo::lchoose_threadsafe`), bit-identical over 2,009,004 pairs and on a
worker with R's stack check LEFT ON -- so `R_CStackLimit` is never touched and
genuine overflow protection is never surrendered. It changed no behaviour alone.
`7f897285` pins it, because the goldens provably cannot: they execute the branch
(a hard `Rf_error` there fails the tier twice) but a `+1e-9` perturbation leaves
the entire tier green.

**Step B** deleted the blanket deferral, leaving the `SL_setseed` throw as the
sole guard -- the principle S6 stated and this finally leans on: the throw is the
safety property, a predicate is only ever an optimisation.

Result: in-process speedup goes flat in f (spread 1.290 -> 0.170), at f = 1 from
0.867x -- slower than serial -- to **2.077x**, and head to head it matches or
beats fork at every f while using 80-98 MB less.

The safety gate was possible for the first time here, since a worker now runs
the ER machinery: `harness/tsan.sh` clean at 4 and 8 threads on
binary/binary_sparse/quant **and on `rare100`, where every marker takes the exact
test on a worker**, with the plant control still firing there (157/200 races,
output DIFFERS, 1-thread clean). Use `SAIGE_TSAN_WORKLOADS=rare100`: the bundled
binary fixture sends only 4.1% of markers to ER.

## Where to start

**The retirement, and it is now engineering rather than measurement.** Item 3 is
settled in favour of in-process: equal-or-better throughput at every f, 80-98 MB
lighter, and it supports restart where fork refuses.

What remains is **criterion item 2 alone** -- feature parity item by item (group
tests, restart, non-fork-safe readers) plus the recorded decision on
`varianceRatioBinsHaveGap` -- then route `--nThreads` to `nThreadsMarker`, retire
`SAIGE_MARKER_THREADS`, and delete `splitMarkerShards`, `appendShardBytes`,
`mergeShardOutputs` and `runMarkerScan`'s fork branch. **Keep E0 and E1** (both
recorded in `CLAUDE.md`); they are findings about this codebase, not about fork.

Worth doing with it: have the marker loop report how many markers took ER, the
way it reports its thread count. Nothing else makes that visible.
- **Criterion item 2 is still entirely open**: feature parity with fork item by
  item -- group tests, restart, non-fork-safe readers -- plus the recorded
  decision on `varianceRatioBinsHaveGap`.

One timing rule, unchanged and re-earned this session: check `pmset -g therm`
reads `CPU_Speed_Limit = 100` first, and check what is actually running. A
Spotlight re-index at ~60% CPU and a `ReportCrash` at ~97% (fallout from the
planted-race abort) each delayed a block here; a paired A/B cancels noise SHARED
by the arms, not a burst that lands on whichever arm was running.

S6 deliberately leaves TWO parallelism mechanisms in the tree at once -- the
existing fork/`mclapply` path in `R/MarkerShard.R` and the new in-process one,
default OFF. **That is a temporary state with a written expiry**: the retirement
criterion for the fork path is in the status block of
`docs/design/step2-threading-verdict.md` and lists the evidence required before
either one is deleted (byte-identity, feature parity item by item, a paired A/B
showing in-process is not slower, and maxRSS). Do not let S7 close without
resolving it -- two mechanisms with different semantics is exactly how this
codebase grows a path nothing tests, which is what S5 found in the pread
fallback and what S0 found in the marker-path arms.

While reading that area: `varianceRatioBinsHaveGap` in `R/MarkerShard.R` had a
comment claiming LEDGER #32 was "preserved (branch policy)". S3 FIXED #32 --
`selectVarianceRatio` throws now, both arms stop, and the guard survives only
because its message arrives earlier and reads better than the C++ throw. The
comment is corrected; the guard's fate is one of the parity items above.

## Things worth knowing before touching this area

- **The goldens cannot see a thread-safety change.** S2, S3 and S4 are all
  byte-identical by construction, so each carries its own measurement instead:
  S2 counts `set.seed` calls by tracing `base::set.seed` in a child process
  (`tests/testthat/test-step2-er-seeding.R` keeps that as a standing test), S3
  used the compiler as an exhaustive audit, S4 used the type system (`grep` for
  `std::vector<bool>`; every `g_model->`/`g_cfg.` inside `testOneMarker` read
  back and confirmed a read), **S6 ran the goldens at 1/2/4/8 threads and drilled
  BOTH directions -- a deliberately shared scratch buffer must fail the
  comparator (it does, 16 failures) and 1-thread-through-TBB must pass (it does),
  which is what separates "TBB is involved" from "markers ran concurrently"** --
  and **S5 planted defects and watched the new tests catch them** -- a wrong pread offset (equivalence test reports
  `identical: FALSE`, 2 of 5 markers diverge) and the pre-S5 `close()` semantics
  (clone-after-close fails with "called after the .bed file was closed").
  Plan for that; a green golden run is not evidence for this class of change.
  What the golden cycle DID gate on S4 is the thing it is good for: that a
  460-line function was transcribed correctly.
- **Paired A/B with two installed libraries** is the only timing method that
  resolved anything here: `R CMD INSTALL -l <scratch lib>` for the baseline,
  the pixi env for the change, arms alternated within each rep, and the arm
  order swapped between blocks as an ordering control. Unpaired 3+3 reps called
  F8 a non-result and would have called S1 one too. See `harness/bench/NOTES.md`.
- **`s2_marker_all_sparse`** is the bench workload that puts the marker-path
  `scoreTest` and the PCG on the clock; `s2_marker_all` cannot see them at all.
- Three ledger entries were opened by reading code in these commits and are
  untouched: **#44** (BETA/SE unguarded when the covariate-adjusted variance
  underflows -- two bundled markers ARE covariates), **#45** (a dead
  `== NA_REAL` NaN guard in both root finders), **#46** (the two root finders
  disagree on the sign-change test when an operand is exactly zero).
