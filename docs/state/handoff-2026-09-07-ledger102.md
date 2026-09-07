# Handoff: LEDGER #102 fixed, the backlog pushed (2026-09-07)

Two code/harness commits on top of `9347f8e7`, both gated, and the 14-commit
backlog is now on `slim/strip`. Previous handoff:
`docs/state/handoff-2026-09-07-overnight.md`, which opened #102 and recorded the
push the pre-push gate refused.

| | |
|---|---|
| start | `9347f8e7`, `slim/strip` at `b11e6b07` |
| end | code at `96f16942`, comparison at `99541489`, plus this file's own commits; pushed from `b11e6b07` in steps |
| CI | all green: `34101270841` (the two code commits), `34102287183` (this handoff), `34104489021` (section 3). `gh run list --repo astheeggeggs/saige-slim` |

Sections 1 and 2 were written and pushed before section 3 existed, so this table
names commits later than the file's first version.

## 1. `cca6b497` -- LEDGER #102, the out-of-bounds read in `setStartEndIndex`

**Deleted the loop, the `Msub` member and its `getMsub()` accessor**, rather
than bounding the loop by `alleleFreqVec.n_elem`. The check the brief asked for
came back unambiguous: nothing reads `geno.Msub`. Its only accessor had two
callers, `CorssProd_LOCO::m_Msub` (assigned in both constructors, never read in
`operator()` or `join`) and a local in `getDiagOfSigma_LOCO` that the function
never uses -- that function divides by
`Msub_MAFge_minMAFtoConstructGRM_in_b - ..._singleChr_in_b`, a different
quantity. Both go with it. `setStartEndIndex` now records only
`startIndex`/`endIndex`/`chromIndex`, which the LOCO cross-product does read.

Deleting rather than bounding, because **the two index spaces are not the same
space**: `i` runs over file order and `alleleFreqVec` is in QC-passing order, so
a bounded loop would have kept computing a wrong number and left it available
for a future reader to trust.

**Golden prediction, made before the gate: none move.** Confirmed --
`run_tests.sh full` in byte mode, which regenerates the goldens and compares
them, came back `[ FAIL 0 | WARN 2 | SKIP 0 | PASS 2919 ]`. This is the case
CLAUDE.md's first trap warns about ("a fix that moves no golden is a statement
about coverage") reached from the other end: the value is dead, so stability is
the prediction rather than a gap.

### The characterization, and why it is Guard Malloc

A repeat count is not a test of an intermittent overrun, and the brief's two
suggestions do not both work:

- **A bounds-checked (`ARMA_NO_DEBUG`-less) build catches NOTHING here.** The
  access is `alleleFreqVec[i]`, and armadillo's `operator[]` is unchecked in
  every build; only `operator()` is guarded, and only when `ARMA_NO_DEBUG` is
  absent. Recorded so nobody tries it again.
- **Guard Malloc works and is deterministic.**
  `/usr/lib/libgmalloc.dylib` puts every allocation on its own page with an
  unmapped page after it, so the read faults on the instruction that performs
  it. New harness `harness/guardmalloc.sh`, driven in the standard tier by
  `tests/testthat/test-step1-loco-oob.R` (macOS only, ~11 s).

**Drilled both ways.** Against a snapshot of the `9347f8e7` install the clean
arm is `status=139`, traceback frame 1
`setStartEndIndex(startIndex, endIndex, j - 1)` -- the same frame the pre-push
gate crashed in -- and the test reports 2 failures. After the fix the arm exits
0 and the test passes 4.

**The plant arm is not optional, and this is the finding to carry forward.**
Guard Malloc's standard failure mode on this box is to look like it is working:

1. **`Rscript` loses `DYLD_INSERT_LIBRARIES`.** It routes through a
   SIP-protected binary, which strips `DYLD_*`. The GuardMalloc banner still
   prints -- from the pre-exec process -- while the R process that runs the fit
   allocates from ordinary malloc. Measured: with the banner in the log, an
   `arma::fvec` sized 1666 read 334 elements past its end and returned zeros,
   and `posix_memalign` handed back ordinary heap addresses. **So the banner
   proves nothing.** `$R_HOME/bin/exec/R` is invoked directly instead, the same
   way `harness/tsan.sh` already does it.
2. Bypassing `bin/R` also loses the `DYLD_FALLBACK_LIBRARY_PATH` that
   `etc/ldpaths` sets, without which `SAIGE.so` cannot resolve
   `libplink2_includes.a`.
3. The plant is one deliberate `Get_OneSNP_StdGeno(1999L)` injected at the same
   program point -- index 1999 of a 1,666-element `alleleFreqVec`, and that
   function's first act is `alleleFreqVec[SNPIdx]`. It must survive with the
   guard off and be killed by it. It passes in BOTH directions of the drill, so
   the red and the green were both taken with a tool proved live.

One trap inside the plant: **survival is signalled by creating a file, not by
printing.** R deparses a tracer into its segfault traceback, so a sentinel
string in the expression appears in the log of the run that crashed at it --
the first version of the check reported the plant as surviving when it had died.

## 2. `96f16942` -- `harness/tsan.sh` gets a per-arm timeout

The gap item 2 of the previous handoff named. `harness/timeout.sh` is a shared
bash watchdog (macOS ships no `timeout(1)` and the pinned env has no coreutils);
`SAIGE_TSAN_TIMEOUT` defaults to 3600s per arm. It kills the arm **and the
children it forked** -- a TERM to the parent alone leaves SkatoParallel's pool
alive and reparented -- and gives it `/dev/null` on stdin, which is what ends
R's interactive abort prompt. An arm that dies or times out is reported as
`died (exit N)` / `TIMED OUT after Ns` instead of as `races=0` with a byte
comparison against a truncated file; under `plant` that is the mode working,
under `clean` it now makes the script exit 1 (`tsan.sh` previously always
exited 0).

Verified on the binary workload capped at 2,000 markers, one build per arm of
the drill: with the default timeout, `races=0` at t=1/4/8, bytes-vs-serial
identical, `rc=0` -- unchanged from before the edit; with
`SAIGE_TSAN_TIMEOUT=1`, all three arms report `TIMED OUT after 1s` and the
script exits 1.

`tests/testthat/test-meta-timeout.R` pins the watchdog in the standard tier,
because a watchdog that fails to kill leaves the harness hanging exactly as
before and reports nothing. It asserts statuses pass through unchanged (0, 3,
139 for SIGSEGV), that `read` comes back 1 at EOF rather than 124, and that a
TERM-ignoring process and the child it forked are both gone afterwards --
observed through files the processes touch, not `pgrep`, whose pattern would
have to be general enough to match whatever else the box is running. Drilled
with two plants: dropping the `pkill -P` lines, and dropping every kill, each
leave the child's end-marker behind and fail the test.

## Gates

| gate | result |
|---|---|
| `run_tests.sh full`, byte (regenerates + compares goldens) | `[ FAIL 0 \| WARN 2 \| SKIP 0 \| PASS 2919 ]` |
| `prepush.sh` byte | `[ FAIL 0 \| WARN 2 \| SKIP 4 \| PASS 2901 ]` |
| `prepush.sh` numeric | `[ FAIL 0 \| WARN 2 \| SKIP 4 \| PASS 3101 ]` |
| `tsan.sh clean`, binary, 2,000 markers | races=0 at t=1/4/8, bytes identical, rc=0 |
| `guardmalloc.sh both` | clean 0, plant survives without the guard, dies under it |

`WARN 2` is the standing baseline (`test-unit-varianceratio.R`). The pre-push
gate ran twice at `96f16942`, once by hand and once from `.githooks/pre-push`,
green both times; it was NOT bypassed. The byte counts are +12 on the previous
green (2889): 4 from `test-step1-loco-oob.R`, 8 from `test-meta-timeout.R`.

## 3. `y_quant_pc` and `y_null` own models refitted -- the table is now uniformly post-#97

Item 3 of the previous handoff, done after the push. `MODELS=own run_arm.sh strip
{step1,s2single,s2group} <trait> 1` for both traits, then `compare_all.sh 1 1`;
only the strip arm needs it, since #97 is a strip-only fix. Predictions were
written down first (`sparse` quotient → 1 in both MAC categories, `null` rows
unmoved, gene agreement up ~2 orders of magnitude in max|d|, `VR-mix` resolving
INTO `VR` rather than vanishing, the single-variant tables barely moving) and
all of them held.

**The negative control on the procedure**: of the 22 blocks in
`cmp_summary.txt`, exactly the six belonging to the two re-run traits changed;
every `frozen` block and every block of `y_binary`, `y_quant` and `y_binary_pc`
reproduced byte-for-byte.

| | `y_quant_pc` | `y_null` |
|---|---|---|
| `sparse` VR quotient strip/upstream | 1.11979 -> **1.00000** | 1.11169 -> **1.00000** |
| gene `Pvalue` Pearson, max\|d\| | 0.999716 -> **0.999995**, 0.609 -> **0.180** | 0.999143 -> **0.999966**, 0.218 -> **0.116** |
| gene `Pvalue_Burden` | 0.999226 -> **0.999996**, 0.609 -> **0.078** | 0.999249 -> **0.999997**, 0.201 -> **0.014** |
| gene cells matching a VR row / between | 1626/1455 -> **2421/677** | 1625/1466 -> **2617/490** |
| single-in-group `var` quotient MAC<=20.5 | median 1.1197910 -> **1.0000000** | 1.1116919 -> **1.0000000** |

The `20:17775656:A:T` fast-test straddle on `y_quant_pc` is still the only
`UNPREDICTED` row, before and after -- it is a cutoff straddle, not a
variance-ratio effect, so #97 was never expected to move it, and none appeared.

`RESULTS.md` gains a second "After the #97 fix" subsection and its preamble now
says all three quantitative traits are post-fix; the *(pre-#97)* row labels stay
as the record of what #97 did.

One thing to know when reading the two subsections together: the `y_quant`
re-run was on the pre-#102 build and these two on the post-#102 build. #102
deleted a dead value and moved no golden, so the builds are output-equivalent,
but the `so_md5` recorded in `results.tsv` differs between them (`29e41ccd...`
against `decbe51d...`) and that is why.

## Still open

1. LEDGER #91's `WHERE` cites `src/SAIGE_fitGLMM_fast.cpp:2504`/`:2533`; the
   two floors are now at `:2508`/`:2536`. It was already stale before this
   session's deletion shifted it by one, and was left alone rather than churn a
   row that is not being fixed.
