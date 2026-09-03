# S8: PGEN gets marker parallelism (handoff, 2026-09-03)

Branch `strip`. Picks up exactly where `docs/state/s0-s7-handoff.md` left off:
its "Next" section named `PgenReader::clone()`, and this is that.

## What changed

**`PgenReader::clone()` is real, so `--nThreads` threads a pgen marker scan.**
Before this, a pgen input degraded to one thread no matter what was asked for,
and since S7 retired the forked path that was the ONLY route to parallel PGEN.

The implementation is deliberately not `PlinkReader`'s shape, and the reason is
the whole design:

| | PLINK (S5) | PGEN (S8) |
|---|---|---|
| what a read touches | `const` pointer arithmetic into an mmap | MUTATES the reader: `fp_vidx`, fread buffer, LD cache, workspaces |
| so clones... | SHARE one `PlinkClass` | each get their OWN `plink2::PgenReader` |
| what is shared | the mapping and the `FILE*` | `PgenShared`: pvar/psam metadata, sample subset, `PgenFileInfo` + its allocations |
| ownership | dies with the last owner | same -- `PgenShared` by `shared_ptr`, declared first so it is destroyed last |

`src/PGEN.{hpp,cpp}` is split into `PgenShared` (read-only after construction)
and `PgenClass` (one decode handle; `clone()` builds a second against the same
shared state). That is pgenlib's own **mode 2** contract, not an invention:
`PgrInit` with a non-null `fname` copies the `PgenFileInfo` by value into the
reader, the first reader inherits the header's already-open handle and each
later one opens its own, and the arrays that copy shares -- `var_fpos`,
`vrtypes`, `allele_idx_offsets`, `nonref_flags` -- are written only by
`PgfiInitPhase2` (grepped, not assumed).

`PgrInit` is called under a mutex on the shared state, because it MOVES
`shared_ff` out of the `PgenFileInfo` for the first reader and
`MarkerStatePool::make()` can run on a worker. Reads take no lock.

## The bug that fell out: LEDGER #48

`PgrInit`'s `fname` was the **pvar** path, not the pgen. Unreachable while
`clone()` threw -- the first reader onto a file never looks at `fname` -- and
reachable the moment it stopped. Drilled rather than reasoned about: in a
scratch build with the pvar argument restored (installed package untouched), the
clone's FIRST READ fails with `PgrGet1D() error 6`, kPglRetMalformedInput.
`PgrInit` itself succeeds -- it opens the pvar and seeks into it quite happily.
Loud on this fixture; that is a property of the pvar's bytes, not a guarantee.

## The gates, and what each one is worth

- **Goldens: 54/54 byte-identical**, five step-1 goldens at `max_rel_err
  0.000e+00`. Run twice -- once on the code build, once again after the final
  clean build so the installed package is the tree that ships. *Worth little on
  its own*: every golden runs serially, so a broken clone would leave them all
  green. It is the "did not change the serial answer" gate.
- **Standard tier: FAIL 0 | PASS 1769 in BOTH compare modes.** +21 against the
  baseline; 18 measured directly (`test-unit-genoreader.R` in isolation goes
  26 -> 44) and 3 from the replaced pgen threading test (5 assertions out, 8 in).
- **The pgen threading tests are the real gate.** `test-step2-pgen.R`'s
  degrade-to-one-thread test is INVERTED (for the second time -- it began life
  asserting a hard error) and now asserts the run REPORTED 4 in-process threads
  and reproduced the golden byte for byte, plus a second run over the
  128,868-marker fixture where the workers genuinely overlap for seconds.
  Asserting the thread-count LINE is what stops a silent fallback from passing
  as a green threaded run.
- **TSAN clean on `pgen` and `pgen_quant` at 4 and 8 threads: races=0,
  byte-identical to serial** -- and for the first time **with pgenlib itself
  instrumented**, which is what makes that zero mean anything about the decoder.
- **`plantreader`, the matching positive control**: `clone()` hands back a handle
  onto the SAME `PgenClass`. It fires exactly where the design predicts --
  worker-vs-worker inside `plink2::InitReadPtrs`, reached through
  `PgenClass::Read`. Matched pair at 2,000 markers: clean 0/0/0 races and
  identical at t=1/4/8, plantreader 0 / **192** / **54** races with output
  DIFFERS at t>1. The split is load-bearing.
  **The plant is capped at 2,000 markers deliberately** (`SAIGE_TSAN_MAX_MARKERS`,
  defaulted for this mode): a shared decoder is not just wrong, it is
  pathologically slow -- several threads thrash pgenlib's LD cache, and the
  full-fixture arm did not finish chunk 1 of 22 in eighteen minutes where the
  clean build does all 22 in about two. A control nobody can run to completion is
  not a control. The cap applies to the 1-thread arm too.
- **Speed, paired A/B at 4 threads** (6 reps x 2 blocks, arm order swapped
  between blocks, `CPU_Speed_Limit = 100` checked first): `work_s` speedup
  **2.177x** (95% CI 2.051-2.295) binary, **2.040x** (1.922-2.150) quantitative
  -- alongside PLINK's 2.06-2.23x on the same fixture. Details in
  `harness/bench/NOTES.md`.
- **A direct decode drill** for the unit test's power: reading the same rows
  through the PLINK reader matches the pgen clone exactly, and reading rows
  shifted by one does not.

## Two harness scripts were silently useless, and that is the reusable lesson

`harness/tsan.sh` and `harness/bench_threads.sh` both still drove the threaded
arm through **`SAIGE_MARKER_THREADS`, which S7 deleted**. Every "threaded" arm in
either script had been running serially since then. This is the LEDGER #30 shape
for the second time: a knob that fails to arrive looks exactly like a green run.
Both scripts carry their own DID-NOT-THREAD check, which is how it surfaced;
keep those checks in anything similar. `bench_threads.sh` also still had the
`fork` mechanism dimension -- it now refuses with a pointer to S7.

**And a new trap of the same family.** Instrumenting pgenlib for TSAN is not
just "compile it with `-fsanitize=thread`": R puts `-L<pixi env>/lib` on the link
line BEFORE `PKG_LIBS`, so `-lplink2_includes` resolves to the uninstrumented
archive no matter where an extra `-L` is prepended. The build succeeds, the run
works, and TSAN just sees nothing. `tsan.sh` links the instrumented copy by
ABSOLUTE PATH and the check is `otool -L .../SAIGE.so | grep plink2` -- it must
print a path under `harness/work/tsan`, not a bare filename.

## What is NOT established

- **PGEN is still not fork-safe**, and `isForkSafe()` returning `false` is now a
  live hazard rather than conservatism: a pgenlib reader owns a `FILE*` it
  fseeks. Nothing in the pipeline forks any more (S7), so nothing asks.
- **Memory at scale.** The per-handle pgr arena is ~130 bytes per SAMPLE --
  arithmetic from the allocation expression, not a measurement -- most of it two
  transpose buffers that single-variant reads never touch. At N = 1000 maxRSS
  moves ~1 MB between 1 and 8 threads; at biobank N this is tens of MB per
  thread, and PLINK's shared mmap pays none of it (measured: maxRSS moves +0.6 MB
  going 1 -> 4 threads at N = 1000, which is what ~130 B/sample predicts and is
  why the bundled fixture cannot see this). **The obvious follow-up is to stop
  allocating the buffers `Read()` never uses**, which needs care with the arena
  arithmetic and a fixture big enough to measure.
- **Group tests still run on one thread**, for every format.
- The one remaining reader that refuses to clone is PLINK's pread fallback (S5),
  still refused for want of scale testing rather than for a hazard.

## Where to start next

**Correction to the S0-S7 handoff, which this doc first copied:** step-1's
blocked `GetTrace` is NOT "sitting on branch `gettrace-blocked` and never
merged". `gettrace-blocked` is an ANCESTOR of `strip` (`git merge-base
--is-ancestor` confirms), `bda95327` is in this history, and
`getPCG1ofSigmaAndVectorMat` is live at two call sites in
`src/SAIGE_fitGLMM_fast.cpp`, unchanged since. That item is done and shipping;
do not re-do it.

So what is actually left:

- **LEDGER #47** -- restart after a hard kill DUPLICATES rows (3,826 of them on
  the bundled fixture), because the `.index` and the output stream flush at
  independent points. User-visible, pre-existing, orthogonal to threading, and
  with a measured reproduction already written down. The obvious next fix.
- The other open ledger entries: #6 (deliberately deferred), #14 (unverified),
  #17 (uninitialised `w0`, real UB on a reachable path), #18-#26, #35-#39, #42
  (a halved log in the `log.p=TRUE` SE back-calculation), #45, #46.
- **Group tests still run on one thread, for every format.** The marker path is
  done; the region path has had none of this work.
- The PGEN arena trim above, if a large-N fixture appears.
