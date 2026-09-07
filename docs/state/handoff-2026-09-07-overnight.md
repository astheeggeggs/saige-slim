# Handoff: overnight gates on the step-5 build (2026-09-07)

No code changed this session. Four gates were run sequentially on an otherwise
idle box against `strip` HEAD **`30832166`**, with one intervention (a hung
TSAN plant arm, killed -- see item 2) and one document updated
(`harness/sim/upstream_compare/RESULTS.md`, item 3). Previous handoff:
`docs/state/handoff-2026-09-06-region-skato-batch.md`.

Driver and logs: `harness/work/region_profile/overnight/`
(`run_overnight.sh`, `timeline.txt`, one `<stage>.log` and `<stage>.rc` each).
Nothing was timed and no timing claim is made below; the wall figures quoted
are there only to identify the hang.

## The HEAD the brief named

The brief named `e1214627`. HEAD is `30832166`. They are the same handoff
commit before and after an `--amend`: same parent `1c706bd2`, same message, and
`git diff e1214627 30832166` is one line of `docs/state/next.md` (the pointer).
No code differs, so the gates ran on the step-5 tree either way. The installed
library was confirmed to be that build before starting -- `SAIGE.so` mtime
21:45 against a newest source of 21:44 -- and the **installed** `DESCRIPTION`
carries `Config/testthat/parallel: true` and `edition: 3`.

## 1. Test suite, full tier, both compare modes -- PASS

| mode | result | parallel |
|---|---|---|
| `SAIGE_COMPARE_MODE=byte` | `[ FAIL 0 \| WARN 2 \| SKIP 0 \| PASS 2915 ]` | `Starting 8 test processes.` |
| `SAIGE_COMPARE_MODE=numeric` | `[ FAIL 0 \| WARN 2 \| SKIP 1 \| PASS 3109 ]` | `Starting 8 test processes.` |

`WARN 2` is the standing baseline, identical to `gate_s5.log`: two deliberate
warnings from `test-unit-varianceratio.R` exercising the no-`varianceRatioFile`
path. The numeric mode's single SKIP is `test-full-pipeline.R:18`, reason
"full tier is byte-exact and therefore same-platform only" -- the golden
regeneration, which by design runs in byte mode only, and did, with `SKIP 0`.
`ledger lint OK: 101 rows, 5 open, 4 decided` in both.

## 2. TSAN -- clean PASS; plant reported races, then one arm segfaulted

**`otool -L harness/work/tsan/lib/SAIGE/libs/SAIGE.so`** names both required
things, so the instrumentation is real:

```
/Users/dpalmer/Repositories/SAIGE/harness/work/tsan/libplink2_includes.a
@rpath/libclang_rt.tsan_osx_dynamic.dylib
```

(that is the plant build; `harness/work/tsan/build.clean.log`'s link line names
the same absolute archive, since both modes build into one `lib/`).

**clean -- 12 / 12 arms, 128,868 markers each:**

```
binary         t=1/4/8  races=0  bytes-vs-serial=identical
binary_sparse  t=1/4/8  races=0  bytes-vs-serial=identical
quant          t=1/4/8  races=0  bytes-vs-serial=identical
region         t=1/4/8  races=0  bytes-vs-serial=identical
```

**plant -- races where they are required, none where they are not:**

```
  binary         t=1  races=0    bytes-vs-serial=n/a
  binary         t=4  races=118  bytes-vs-serial=DIFFERS
  binary         t=8  races=126  bytes-vs-serial=DIFFERS
  binary_sparse  t=1  races=0    bytes-vs-serial=n/a
  binary_sparse  t=4  races=722  bytes-vs-serial=DIFFERS
  binary_sparse  t=8  races=330  bytes-vs-serial=DIFFERS
  quant          t=1  races=0    bytes-vs-serial=n/a
  quant          t=4  races=215  bytes-vs-serial=DIFFERS
  quant          t=8  races=121  bytes-vs-serial=DIFFERS
  region         t=1  races=0    bytes-vs-serial=n/a
```

Six threaded arms report races and a wrong answer; all four serial arms report
none. The control is discharged: a TSAN setup that saw none of this code would
report zero everywhere, and this one does not.

**`region t=4` under the plant crashed** at 01:03:40 and the arm then hung for
4h32m until killed. `harness/work/tsan/runs.plant/region.t4.log` ends, verbatim:

```
Analyzing chunks (1/1, 1/1)........

 *** caught segfault ***
address 0x0, cause 'unknown'
```

and nothing more. R's segfault handler printed that and blocked on its
interactive abort prompt; with `--no-echo` and no stdin it never returns, which
is the 83.8%-CPU master with four idle `SkatoParallel` fork children that the
process table showed. Crash reports `~/Library/Logs/DiagnosticReports/`
`R-2026-09-07-010245.ips` and `R-2026-09-07-010302.ips`.

**This is the plant, not a defect in the code, and the controls bracket it.**
The plant shares ONE `MarkerScratch` across every worker, which is licensed to
produce anything; the same `region` workload is clean at t=1/4/8 on the
uninstrumented-defect build (above), and `region t=1` under the plant passes,
so the crash needs plant AND threads. `region t=8` was not run.

`bash harness/tsan.sh plant` was killed (`rc=143`, recorded in `timeline.txt`);
the driver advanced. The plant was not restarted: `region t=8` would repeat the
same segfault-hang and would add nothing to six arms that already reported
races.

**Harness gap, not fixed here (wants its own commit).** `harness/tsan.sh` has
no per-arm timeout, so a plant arm that crashes hangs the script indefinitely
rather than failing. A plant is *expected* to misbehave arbitrarily, so this
will recur. Suggested: wrap each arm in a timeout and treat "arm died or timed
out" as a plant SUCCESS signal distinct from `races=0`, which is the outcome
the mode exists to rule out.

## 3. Own-model quantitative re-run after the #97 fix -- prediction confirmed

`MODELS=own bash run_arm.sh strip s2single y_quant 1` (rc=0), the same for
`s2group` (rc=0), then `bash compare_all.sh 1 1` (rc=0). Both runs recorded
loading the strip library with the same `so_md5` `29e41ccd…`, which is the
wrong-build null control `run_arm.sh` exists to catch.

`cmp_summary.txt` was snapshotted before and after (`3_cmp_summary_{before,after}.txt`).
**The diff is five hunks, all inside the three `y_quant models=own` blocks** --
every `frozen` table and every other trait reproduced its previous text exactly,
which is the negative control on `compare_all` itself.

The `sparse` variance-ratio quotient strip/upstream went **1.11126 → 1.00000**
(both categories; the other four rows unmoved). Headline movements:

| table | before | after |
|---|---|---|
| gene `Pvalue` | Pearson 0.998997, max\|d\| 0.389 | 0.999991, max\|d\| 0.027 |
| gene `Pvalue_Burden` | 0.999224, max\|d\| 0.235 | 0.999995, max\|d\| 0.051 |
| gene cells matching a VR row / between rows | 1630 / 1452 of 3115 | 2547 / 552 |
| gene tags | `VR` 1618, `VR-mix` 1330, `Cauchy=1` 7, `SKATO=1` 1 | `VR` 1987, `VR-mix` 495, `tiny` 471, `Cauchy=1` 0, `SKATO=1` 0 |
| single-in-group `var` quotient, MAC ≤ 20.5 | median 1.11126 | median 1.0000000 |
| single `var` quotient, MAC ≤ 20.5 | range [1.0021371, 1.1112647] | [1.0000000, 1.0021427] |

`next.md` predicted the own `gene y_quant` tags would "lose their `VR`/`VR-mix`
mass". Directionally yes, but the precise statement is that the mass
**resolved** rather than vanished: `VR-mix` fell 63% and `VR` ROSE, because
`VR-mix` means a cell's quotient lay between the arms' two variance-ratio rows,
and with the `sparse` row no longer displaced by 1/s² those cells now match one
row cleanly. The reconciliation counts say it directly (1630 → 2547 matching).
The residual `VR` mass is the `null` row quotient 1.00530 -- #65/#53, a
different discrepancy #97 never touched, which is why the single-variant table
barely moves (0.973418 → 0.973446) while the gene tables move two orders of
magnitude in max|d|. No row became `UNPREDICTED`; the one `UNPREDICTED` in the
file is the documented `y_quant_pc` fast-test straddle, in a trait not re-run.

`harness/sim/upstream_compare/RESULTS.md` updated: four table rows marked
*(pre-#97)*, a note that the table is the 2026-09-05 state, and a new subsection
"After the #97 fix: `y_quant` own models re-run (2026-09-07)" with the above and
its reasoning. **`y_quant_pc` and `y_null` own models were never refitted after
#97** (their `.rda` files are the 2026-09-05 16:50 fits), so their rows still
describe pre-fix behaviour and are labelled as such -- if a fully post-#97
own-model table is wanted, those two traits need step 1 re-run first.

## 4. Push -- REFUSED by the pre-push gate. NOT pushed.

`git push slim strip` ran `harness/prepush.sh`, which passed byte mode
(`FAIL 0 | WARN 2 | SKIP 4 | PASS 2889`) and then FAILED numeric mode, so
git rejected the push. Verbatim:

```
FAILURE: 'test-step1-lowmem-loco.R:91:3' ----------
CLI exited with status 139; log tail:
  leave chromosome  1  out

 *** caught segfault ***
address 0x7fe596000000, cause 'memory not mapped'

Traceback:
 1: setStartEndIndex(startIndex, endIndex, j - 1)
...
Error in `load(path, e)`: empty (zero-byte) input file
[ FAIL 2 | WARN 2 | SKIP 4 | PASS 3074 ]
error: failed to push some refs to 'https://github.com/astheeggeggs/saige-slim.git'
```

The two FAILs are one event: the step-1 CLI segfaulted, so its `.rda` was
zero-byte and the next `load()` failed on it.

`slim/strip` is still at `b11e6b07`, a clean ancestor of local HEAD -- the
rejection is the hook, not a non-fast-forward. **There is no CI run id to
report.** Local commit `441e4417` stands, unpushed. The gate was NOT bypassed
(`SAIGE_SKIP_PREPUSH=1` exists and was not used); whether to push a docs-only
commit over a pre-existing intermittent crash is the owner's call.

### What the crash is: LEDGER #102, opened here

Not a flake and not caused by this session's commit (docs only). It is an
out-of-bounds read that has been in `setStartEndIndex` all along:

- `src/SAIGE_fitGLMM_fast.cpp:3316` loops `i` over `geno.M` -- ALL markers in
  the plink file, set once at `:1014` from `t_src.numMarkers()` and never
  reduced -- and reads `geno.alleleFreqVec[i]`.
- `alleleFreqVec` is sized at `:1206` to
  `numberofMarkerswithMAFge_minMAFtoConstructGRM`: the markers that PASSED QC
  (`maf >= minMAFtoConstructGRM && missingRate <= maxMissingRate`, `:175`).
- So the read runs `M - nPass` elements past the end on any LOCO fit with a
  marker under `--minMAFforGRM` -- i.e. on the default (0.01) with real data.

Quantified independently of SAIGE, by decoding the fixture's `.bed` in R:
`sm_two_chr_fixture` has **2,000 markers, 1,666 passing `--minMAFforGRM=0.01`**,
so the loop reads **334 floats = 1,336 bytes past a 6,664-byte allocation**.
That size of overrun normally lands in mapped heap and returns garbage
silently; it faults only when the allocation ends near the end of a mapped
region.

Observed rate, all on this build the same night: **1 failure in 8 standard-tier
runs** -- the two full tiers, the byte and the numeric half of the gate, and
four deliberate numeric-mode repeats afterwards, which passed 4/4 at
`FAIL 0 | WARN 2 | SKIP 4 | PASS 3089` -- and **0 failures in 12** runs of
`test-step1-lowmem-loco.R` on its own. Eight parallel test processes vary the
heap enough to hit it; one process did not, in 12 tries.

**The impact is bounded to the crash, because the value is dead.** `geno.Msub`
is read only through `getMsub()`, whose two live callers are `:1737`
(`m_Msub`, never read again) and `:2520` in `getDiagOfSigma_LOCO`, an unused
local -- that function divides by
`Msub_MAFge_minMAFtoConstructGRM_in_b - ..._singleChr_in_b`, a different
quantity. No output depends on it. That is also why every golden is stable, and
it makes the prediction for a fix unambiguous: **a fix must move no golden**,
which here is the expected result rather than the coverage gap CLAUDE.md warns
about.

Proposed fix, for a separate commit under the one-row-one-commit rule: bound
the loop by `alleleFreqVec.n_elem` (or delete the loop and the dead `Msub`
member). A deterministic characterization is the hard part -- `ARMA_NO_DEBUG`
is what hides the read -- so the honest gate is an ASAN or bounds-checked build
over `test-step1-lowmem-loco.R` rather than a repeat count.

## Open items for the owner

1. Decide whether `441e4417` (docs only) may be pushed over the #102 crash, or
   whether #102 is fixed first and the two go together.
2. `harness/tsan.sh` needs a per-arm timeout (item 2).
3. `y_quant_pc` and `y_null` own models are still pre-#97 fits; a fully
   post-#97 own-model table needs their step 1 re-run first (item 3).
