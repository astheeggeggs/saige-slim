# Handoff: the two remaining region-path speedups (written 2026-09-04)

This is the standalone brief for the two items that were left on the board after
the SKAT-O parallel work landed, written so the next session can start from the
numbers rather than re-deriving them:

* **Item 3 — batch the SKAT-O cells.** Close the gap between the measured 1.30x
  and the 1.43x ceiling. Small, low risk, gate already exists.
* **Item 4 — thread the C++ region loop.** The other 43.9%. Composes with item 3.
  The big one, and the one with an unresolved safety question.

Everything below is measured on `sim50k` (N = 50,000, 210 genes, 5 annotations
x 3 max_MAF, `y_binary`) unless it says otherwise. The fixture lives OUTSIDE
this repo, in the `SAIGE-sim` worktree; it is ~3.1 GB, not in git, and ~65 min
to rebuild. Full tables in `harness/bench/NOTES.md`.

## The one table both items are decided by

`Rprof` at N = 50,000, single-threaded, two replicates agreeing within 0.3 pp:

| bucket | share |
|---|---|
| C++ `mainRegionInCPP` | **43.9%** |
| R `region_postprocess` | **53.8%** |
| .. of which `SKAT_davies` under `integrate()` | **44.9%** |
| .. R glue (post minus SKAT) | 5.0% |
| everything else | 2.2% |

The negative control is the reason to trust the attribution rather than treat it
as a profiler story: the same command with `--r.corr=1` (BURDEN, which never
calls SKAT-O) predicts `SKAT_davies` = 0 and measures **0.0%**, with C++ at
93.7%. If the 45% were an artifact of sampling inside `.C`, it would not switch
off.

**Amdahl on that fraction, threading the region loop only:**

| test type | C++ fraction | 4 threads | 8 threads | infinite |
|---|---|---|---|---|
| SKAT-O | 43.9% | **1.49x** | 1.62x | 1.78x |
| BURDEN | 93.7% | **3.36x** | 5.55x | 15.9x |

So the natural plan is inverted: for SKAT-O users the C++ region loop is NOT the
thing to thread first. It was not, and item 3 is what took its place.

---

## Item 3: batch the SKAT-O cells -- DONE 2026-09-06 (step 5, commit `1c706bd2`)

The first of the two shapes below landed: `clusterApplyLB` kept, one task per
cell, but the cells of a whole BATCH of regions (a group-file chunk, or fewer
under `SAIGE_SKATO_BATCH_MB`) dispatched in one call sorted expensive-first by
`ncol(Phi)`, the escalated regions' cells in a second call per batch. The
driver holds each region's C++ output list until its answers are back, which
is why the batch has a memory budget. Measurement and the handoff:
`docs/state/handoff-2026-09-06-region-skato-batch.md`, `harness/bench/NOTES.md`
(step 5 section, prediction written before the run). The text below is kept as
the brief it was.


### Where it stands

`R/SkatoParallel.R` (161 lines) drives `region_postprocess` three times —
**collect** the cells through an injected `skat_fun`, **evaluate** them on a
persistent fork pool, **serve** the answers back in call order — so the 175-line
triple-nested cell loop did not have to be restructured. `--nThreads` on a group
run drives it; it used to be discarded there.

Measured, two blocks with the arm order reversed, load ~4.3, ~1% replicate
spread:

| nThreads | mean real | speedup |
|---|---|---|
| 1 | 204.5 s | — |
| 4 | **156.7 s** | **1.30x** |
| 8 | 160.0 s | 1.28x |

Byte-identical to the one-thread run in both the region file and the
`singleAssoc`. **8 workers buys nothing over 4**, which is exactly what a 1.43x
ceiling minus per-task overhead looks like.

### The change to make

One task per worker per region instead of one task per cell. `clusterApplyLB`
currently dispatches every cell individually; at ~15 cells a region and 4
workers that is ~4 tasks per worker per region, each carrying its own `(Score,
Phi)` payload and its own round trip.

**Do NOT simply chunk into equal-sized groups.** Cell cost spans three orders of
magnitude (0.5 ms to 540 ms) and the top 10% of cells hold **48%** of the total,
which is why the current code is load-BALANCED rather than round-robin. A
naive static chunking re-introduces exactly the imbalance `clusterApplyLB` is
there to avoid. Two shapes that keep the balancing:

* keep `clusterApplyLB` but make each job a LIST of cells, sorted so the
  expensive ones are dispatched first (longest-processing-time-first, which
  needs a cost proxy — `ncol(Phi)` is the obvious one and is free);
* or keep one task per cell but batch the *dispatch* across the whole chunk of
  regions rather than per region, so the pool is fed once per chunk.

The second is likely better and is closer to the existing structure, but neither
has been tried and no measurement exists for either.
`[ASSUMPTION: that per-task overhead is what the 1.30x-vs-1.43x gap is. The
basis is that 8 workers buys nothing over 4, which is the signature of a
per-task cost, not a measurement of the overhead itself.]`

### How to gate it

`tests/testthat/test-step2-skato-parallel.R` already exists and asserts both
halves: byte-identity against the serial run, **and** that the pool was actually
built ("SKAT-O post-processing will use N workers", printed only once the pool
exists). The second assertion is the important one — without it a run where the
option never arrived would pass as a green parallel run, which is the LEDGER #30
shape. Keep it.

Byte-identity is the right gate here, despite the standing rule against
comparing two code paths: each cell's SKAT-O p-value is a pure function of that
cell's `(Score, Phi, r.corr)` — no cross-cell reduction, no shared accumulator,
no RNG — so this is one code path re-run over a partition.

### Two traps, both of which cost a full measurement cycle

1. **A worker function must be defined at TOP LEVEL.** R serializes a closure
   with its enclosing environment, so a worker defined inside
   `region_postprocess_parallel` shipped that frame's `...` — the whole
   `outList`, `VarMat` and all — on EVERY task: ~13 ms against a ~40 ms task,
   0.97x at 4 workers, sys CPU 63 s. Hoisted to `skatoEvalJob`, it is 1.29x and
   30 s. **If a worker needs a value, put it in the job.** Batching changes the
   job shape, so this is live again.
2. **Check the obvious suspect before refactoring around it.** The three-pass
   drive looked like the reason the first parallel version was slow. Measured in
   isolation: 180.9 s vs 180.6 s for a single drive. Free. That check saved a
   rewrite of the most intricate function in the package.

### What was measured and REJECTED, so it is not retried

* **Rewriting SKAT-O in C++** — ~1.8% of the run. `SKAT_davies` is already a
  thin wrapper around `.C("qfc")`, so the hot code is compiled already.
* **Loosening Davies' accuracy** (`acc = 1e-6`) — worth 24-30%, and moves genes
  with p < 1e-4 by up to **80% relative**. Not a trade worth making silently.
* **`mclapply` per region** — forking at each region costs 55-90 ms against
  ~350 ms of work: 1.12x at 4 workers and a net slowdown at 8.

---

## Item 4: thread the C++ region loop -- DONE 2026-09-06 for the marker loop

Landed the same day as the `qfc` audit above, after the region C++ was
bucketed and de-churned first (`docs/state/handoff-2026-09-06-region-threads.md`,
`harness/bench/NOTES.md`). What remains of this item is the per-cell tests of
the ultra-rare pass and the BURDEN emitter, which stay serial and cap BURDEN
at ~1.6x on 4 threads instead of the 3.36x projected below; that projection
assumed the whole C++ threads. The text below is kept as the design brief it
was.

### The prize, and the two things blocking it

43.9% of a SKAT-O run and 93.7% of a BURDEN run. Threading it alone gives 1.49x
/ 3.36x at 4 threads. **Composed with item 3** the serial fraction drops to
~11%, which is where the real argument for doing it lives — neither item is
worth much alone on the SKAT-O path, and together they are.

Two things are unresolved and both must be settled BEFORE any code is written:

1. **`qfc` is NOT thread-safe -- SETTLED 2026-09-06, from the binary.** No SKAT
   source ships in the pixi env (SKAT 2.2.5, `libs/SKAT.so` only), so the
   shipped Mach-O was read: `nm -m` on its `__bss` segment lists Davies' 1980
   declaration list item for item as internal-linkage (`_ZL`-mangled, i.e.
   file-scope `static`) objects -- `env` (a 148-byte `jmp_buf`), `r`, `lim`,
   `c`, `n`, `lb`, `nc`, `count`, `intl`, `ersm`, `ndtsrt`, `fail`, `th`,
   `sigsq`, `lmax`, `lmin`, `mean` -- with no thread-local storage anywhere in
   the object. The disassembly shows the four hazards: `qfc_1` opens with
   `setjmp(env)` and `counter()` `longjmp`s to it (a second caller overwrites
   the buffer mid-integration); the landing pad `free`s the global `malloc`'d
   workspace `th` and the entry re-`malloc`s it (double free / use after
   free under concurrency); `lb`/`nc`/`n` are raw pointers into the CALLER's
   `.C` copies, so a second caller makes the first compute over the wrong
   eigenvalues -- silently; and `count`/`lim`/`fail`/`ndtsrt`/`intl`/`ersm`
   are one shared budget and one shared accumulator for every caller.
   **Consequence: SKAT-O can never run in-process on threads with this
   dependency.** The fork pool (`R/SkatoParallel.R`) is correct precisely
   because each forked worker has a private copy of the `.so` data segment. A
   C++ port of `Met_SKAT_Get_Pvalue` that linked this `qfc` would be WRONG,
   not merely 1.8% faster; the only in-process route is a reentrant
   reimplementation of Davies with its state in a struct, gated on every
   group golden. Not planned.
2. **The SKAT-O side stays on the fork pool (see 1), so the shape is
   "threaded C++ region call, then pooled R post-processing", alternating per
   region.** Item
   3's design is a *fork pool driven from R between C++ calls*. That works only
   because the C++ region loop is serial and hands control back to R after each
   chunk. Threading the C++ loop means the R post-processing can no longer sit
   between iterations, so either the p-values come back to a serial R pass
   afterwards (which re-serialises 54% of the run and throws away item 3), or
   SKAT-O moves into C++ (rejected above at ~1.8%), or the two are pipelined.
   **Pipelining is the only shape that keeps both wins and nobody has designed
   it.**

### What is already known to be safe, and what is not

The marker path was threaded in S0-S12 and the machinery it needed is all
present and reusable:

* `SAIGE::TestContext` (`src/SAIGE_test.hpp`) carries the variance ratio, the
  sparse-GRM flag and the noadjCov flag, so **neither step-2 path writes to the
  model** — S3 did the marker path, **S12 did the region path** and DELETED
  `m_varRatioVal`, `m_flagSparseGRM_cur`, `m_isnoadjCov_cur`, `currentContext()`
  and the whole `assign*`/`set_*_cur` family. The region loop is therefore
  already `const` with respect to the model at its three selection sites (the
  rare-marker loop, the ultra-rare loop, the BURDEN collapse) and in
  `assign_conditionMarkers_factors`. **This is the single biggest reason item 4
  is now plausible at all**, and it is easy to miss because S12 was landed as a
  refactor rather than as a threading step.
* Both readers clone: `PlinkReader` shares one `PlinkClass` behind a
  `shared_ptr` (safe because `getOneMarker` is `const` — do NOT add a `mutable`
  member), and `PgenReader::clone()` gives each worker its own pgenlib decode
  handle over shared metadata. The pgen arena is ~130 bytes per SAMPLE per
  handle, so pgen threading costs memory where PLINK's shared mmap does not.
* `harness/tsan.sh` exists and has a **plant mode that must be run too** — a
  clean TSAN run alone gates nothing, because a TSAN build that cannot see this
  code also reports zero races. `plant` (one shared `MarkerScratch`) and
  `plantreader` (a `PgenReader::clone()` that returns a handle onto the SAME
  `PgenClass`: 0 / 192 / 54 races at 1/4/8 threads) are the positive controls.

What is NOT already solved, and is specific to the region path:

* the region path still has two `std::vector<bool>`s in `Main.cpp` — the marker
  path's were removed in S4 precisely because bit packing makes writes to two
  different markers writes to the same word;
* `P1Mat`/`P2Mat` are written to `<SAIGEOutputFile>_P1Mat_Chunk_<i>.bin` and
  read back (LEDGER #57), so chunk file naming becomes a per-thread concern;
* the escalation signal `g_regionEscalatedPass` lives in `Main.cpp` and is
  deliberately NOT a `TestContext` field — sharing it with the per-marker
  sparse-GRM choice is exactly what LEDGER #27 and #50 were. **If threading
  makes you want to move it into the struct, that is the collision coming
  back.**

### Where the numbers come from if they need re-measuring

`harness/bench/NOTES.md`, sections "Group tests at N = 50,000", "Where the OTHER
90% of a gene test goes" and "SKAT-O in parallel". The profile needs no code
change — wrap the CLI script in `Rprof` and let `optparse` read the same
`commandArgs`; the recipe is in NOTES.md. **Check `uptime` before timing
anything**: the first `--nThreads` benchmark ran at load 6-8 and read
0.97x/1.06x; the same build at load ~4.3 read 1.30x/1.28x.

## A cheaper thing to do first, if item 4 looks too big

`--r.corr=1` (BURDEN) is 93.7% C++ and would get **3.36x at 4 threads** from
item 4 alone, with none of item 3's design tension — there is no SKAT-O
post-processing on that path to pipeline against (`region_postprocess` is 0.2%
of a BURDEN run). Threading the region loop **for BURDEN only** is a much
smaller change with a much larger multiple, and it is a genuine use case rather
than a benchmark artifact. It also builds and tests the region-loop threading
machinery against the easy case before the SKAT-O pipelining question has to be
answered.
`[ASSUMPTION: that a BURDEN-only threaded path is separable from the SKAT-O one
in `mainRegionInCPP`. Not verified in the source; it is suggested by the fact
that `--r.corr=1` already takes a different writer (the C++ BURDEN emitter,
which is what golden P9 pins) rather than R post-processing.]`
