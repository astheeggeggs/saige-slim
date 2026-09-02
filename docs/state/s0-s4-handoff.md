# Step-2 threading sequence: S0-S4 landed (handoff, 2026-09-02)

Branch `strip`, pushed to `slim/strip`. S0-S3 plus the two CI fixes were green
on CI (run 33621566011); **S4 has not been pushed yet** -- see "Where to start".
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
| *this commit* | **S4** — `MarkerResult` + `MarkerScratch` + `testOneMarker`; the marker path's `vector<bool>` is gone |

The authoritative plan is the status block at the top of
`docs/design/step2-threading-verdict.md`, which is current as of this commit.
`LEDGER.md`, `CLAUDE.md` and `harness/bench/NOTES.md` are all updated in the
commits that changed the facts they record; nothing is carried only in a commit
message.

## State of the gates

- Full `make_goldens.sh` + `compare_goldens.sh` cycle run on each of the five
  code commits: every step-0 and step-2 golden byte-identical (55 files at S4),
  all five step-1 goldens at `max_rel_err 0.000e+00`.
- Standard tier green in **both** compare modes at each of them.
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

## Where to start

**Push S4 first** -- it is the only unpushed commit, and the standing lesson
below is that local green is not CI green.

Then S5 in the verdict doc: the reader `clone()`, which still throws
("Deferred to the Phase-4 shared-metadata implementation") in both
`GenoReader` implementations. `PlinkReader`'s mmap path is member-write-free
(F6) so clones can share it behind shared ownership; `closegenofile` needs to
become idempotent so a clone cannot outlive the map. `PgenReader::isForkSafe()`
is unconditionally false and no pgenlib work has been done, so PGEN is
out of scope for threading in this pass.

## Things worth knowing before touching this area

- **The goldens cannot see a thread-safety change.** S2, S3 and S4 are all
  byte-identical by construction, so each carries its own measurement instead:
  S2 counts `set.seed` calls by tracing `base::set.seed` in a child process
  (`tests/testthat/test-step2-er-seeding.R` keeps that as a standing test), S3
  used the compiler as an exhaustive audit, S4 used the type system (`grep` for
  `std::vector<bool>`; every `g_model->`/`g_cfg.` inside `testOneMarker` read
  back and confirmed a read). Plan for that; a green golden run is not evidence
  for this class of change. What the golden cycle DID gate on S4 is the thing it
  is good for: that a 460-line function was transcribed correctly.
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
