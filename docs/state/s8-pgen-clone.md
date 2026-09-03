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

- ~~LEDGER #47~~ -- **FIXED in S9, see below.**
- The other open ledger entries: #6 (deliberately deferred), #14 (unverified),
  #17 (uninitialised `w0`, real UB on a reachable path), #18-#26, #35-#39, #42
  (a halved log in the `log.p=TRUE` SE back-calculation), #45, #46.
- **Group tests still run on one thread, for every format.** The marker path is
  done; the region path has had none of this work.
- The PGEN arena trim above, if a large-N fixture appears.


---

# S9: LEDGER #47 fixed -- a restart no longer duplicates rows

## The bug and the fix

The `.index` was written once per chunk from R while the output stream flushed on
its own ~192 KiB boundary. Those boundaries do not align, so a killed run left
rows on disk belonging to chunks the index never recorded, and
`--is_overwrite_output=FALSE` appended them a second time (measured: 3,826
duplicated markers, 132,695 rows against a complete run's 128,869).

Three pieces, and the ORDER is the fix:

1. `flushOutfile_singleInCPP()` (new export in `src/Main.cpp`) flushes
   `g_out.single`.
2. The chunk loop flushes, THEN takes `file.size(OutputFile)`, THEN writes the
   index line -- which now reads `Have completed the analysis of chunk 7 bytes
   123456`. Without the flush the recorded count would omit that chunk's own
   tail, and a restart would discard good rows: the fix would have inverted the
   bug rather than removed it.
3. A restart truncates the output to the recorded size before reopening it.

An index written before this carries no byte count; that restart still runs, with
a warning that rows may be duplicated, rather than being refused.

## LEDGER #49, which the regression test found

The truncation target is the size recorded for the chunk **before** the one the
run resumes at, because **a restart resumes AT the last completed chunk and
re-runs it** -- contradicting `checkOutputFile`'s own message ("we restart the
analysis from the indexChunk + 1 chunk"). The first version of the fix truncated
to the resumed chunk's own offset and the test failed with a whole duplicated
chunk, which is how this surfaced. Left as OPEN: with the truncation in place the
output is correct, and it costs one chunk of recomputation.

## Gates

- **The regression test reconstructs the interruption rather than racing a
  `kill`**: take a complete run, truncate its output to a point PAST a recorded
  chunk boundary, truncate its index to that chunk. That is the same state a kill
  leaves, but deterministic -- a flaky test for a duplication bug is worse than
  none. It asserts the restart reproduces the complete run byte for byte.
- **Drill-verified**: disabling `truncateOutputToIndex` fails it.
- Standard tier **FAIL 0 | PASS 1784 in both compare modes**. +15 on 1769, all
  accounted: 7 from the new test, 8 from `test-unit-util.R`, whose whole-list
  `expect_identical` on `checkOutputFile`'s return was pinning "this list has
  exactly three elements" -- a property nothing depends on, which made an
  additive field look like a regression. Now field-by-field.
- **54/54 goldens byte-identical**: no `.txt` result changed. The marker
  `.index` goldens were refreshed (chunk lines only) -- but note
  `compare_goldens.sh` EXCLUDES `*.index` from the byte gate, so those committed
  files are a record of what the pipeline writes, not a check on it.

## What is still not covered

- **Region-level restart is impossible** (LEDGER #25), so none of this applies
  there, and the region index still records no byte count.
- The output stream is flushed per chunk, not per row: a kill still loses the
  current chunk's work, which is the intended granularity.
- Nothing here defends against a machine crash -- `flush()` reaches the OS, not
  the platter.

---

# S10 and S11: the two remaining live defects

Both were reachable-in-principle bugs that no golden could see, and they were
done in this order for one reason: **#17 fires with fixtures already in the
repo, #42 does not fire at all below n ~ 20,000**. Doing #42 first would have
meant landing it on a byte-identical tier, which says nothing.

## S10 -- LEDGER #17, the uninitialised `w0`

A K=1 CUSTOM-weighted group file plus `--condition` read an uninitialised double
and multiplied the conditional export by it. `w0 = 1.0` is the fix, and 1.0 is
the VALUE rather than a placeholder: a custom set's UR collapse is WEIGHTED
(`max(w_r * G)`), so the weight is already inside the pseudo-marker; multiplying
by `w_r` again would double-count. The K>1 path already does exactly this
(exports the G1-raw row, `region_postprocess` applies 1 for custom sets), so the
fix makes the two paths agree rather than merely making the read defined.

**The measurement, on the same two runs before and after** (K=1 custom vs the
K=2 run's CADD set -- the two group files carry byte-identical weights):

| | before | after |
|---|---|---|
| conditional columns vs K=2 CADD | disagree on 12/12 cells, up to **2.1e-01** rel (BETA_Burden_cond) | **0.000e+00** on all |
| cells with `Pvalue_cond == Pvalue` exactly | 6 (all the all-ultra-rare ones) | 0 |
| UNCONDITIONAL columns | 0.000e+00 | 0.000e+00 |

That last row is what localised the defect to the conditional export rather than
to the collapse. New golden **P34** `plink_group_cond_customweight.txt` (proven
byte-reproducible over two runs before being accepted); the custom-route
reconciliation in `test-step2-multiweight-cond.R` is a comparison the fix made
possible for the first time.

## S11 -- LEDGER #42, halving a log

`qnorm(log(p)/2, log.p = TRUE)` evaluates the quantile at **sqrt(p)**, not p/2.
Three sites, not the two the ledger named -- the marginal SPA SE, the
conditional SPA SE, and the marginal FIRTH qval in the same block. They now share
one named function, `SAIGE::qnormTwoSidedLogP`, so they cannot drift apart again;
the linear-scale arms are untouched and spelled out rather than left polymorphic,
which is why all 56 goldens stay byte-identical.

**The [ASSUMPTION] on reachability is discharged.** The defective line sits
inside `if(t_isSPAConverge)`, so it mattered whether the saddlepoint converges
once the score p has underflowed. It does: at n = 20,000 with |z| = 40 the score
log p is -804 (past `log(DBL_MIN)`) and SPA returns log p = -827 with
`Isconverge` true. Error against the correct quantile: **SE inflated 1.4188x at
the edge (log p = -710), tending to sqrt(2) deeper in** -- conservative-looking,
so a meta-analysis weighting by 1/SE^2 silently down-weights the strongest hits.

**A dead end worth recording.** Shrinking the variance ratio in the model's
`.varianceRatio.txt` DOES drive the printed p into the log format (23 of 70
markers on the bundled fixture), but it cannot reach the defect: `q` is scaled by
`1/sqrt(varRatio)`, which pushes it outside the CGF's range, so the root finder
declines and every one of those 23 came back `Is.SPA = false` -- the line never
runs. Self-consistent data is what reaches it, and that needs a sample size the
bundled fixtures do not have. Hence a kernel-level gate
(`tests/testthat/test-unit-logp-quantile.R` via the `testQnormTwoSidedLogPInCPP`
hook, drill-verified by restoring the halved log) rather than a golden. **Here a
byte-identical tier is the correct outcome, not a coverage failure** -- the
distinction the standing rule exists to make.

## Where this leaves the ledger

Still open: #6 (deferred by choice), #14 (unverified), #18-#26, #35-#39, #45,
#46, and **#49** (a restart re-runs the last completed chunk -- wasted work,
right answer, found by S9's own regression test). The live-defect queue that
S8-S11 worked from is empty.
