# Step-2 threading sequence: S0-S5 landed (handoff, 2026-09-02)

Branch `strip`, pushed to `slim/strip`. Everything through `feb19ded` is green
on CI (runs 33631137151 and 33632167654); **S5 is the only unpushed commit**.
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
| *this commit* | **S5** — reader clone with shared ownership; `getOneMarker` const; the pread fallback tested for the first time |

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

## Where to start

S6, the threading itself. Everything it was waiting on is in place: the marker
path writes no model state (S3), per-marker state is slot-owned (S4), the reader
clones with a const decode (S5). What S6 still owns:
- **The ER resampling deferral** (G1's second half). S2 removed the
  unconditional `SL_setseed`, but a marker that genuinely DRAWS still enters R's
  RNG, so those markers must run on the main thread. The criterion is
  `binary && MAC <= max_MAC_for_ER && nNonZero >= 21`; golden P33 is the one
  that reaches it.
- `parallelFor` with an EXPLICIT `numThreads` (F5: the explicit argument
  overrides `RCPP_PARALLEL_NUM_THREADS`, so `harness/env.sh`'s pin does not
  neuter a threaded test run), contiguous grain ~64, break-truncation, progress
  re-emit.
- The G2 drills: a deliberately-raced scratch build must FAIL the byte
  comparator, and 1-thread-through-TBB must pass.

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
  back and confirmed a read), and **S5 planted defects and watched the new tests
  catch them** -- a wrong pread offset (equivalence test reports
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
