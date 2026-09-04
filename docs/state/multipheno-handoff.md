# Multi-phenotype step 1 / step 2: running handoff

Plan: `~/.claude/plans/can-we-work-to-purring-platypus.md` (written 2026-09-03,
approved 2026-09-04). This file is the resumable record: after every gated
commit it says what is done, which gates passed with their numbers, what is
next, and anything skipped or failed and why. A fresh session should be able to
pick up from the last "DONE" row without re-deriving anything.

Branch: `strip`. Starting commit: `5a33f785`.

## Scope decisions (accepted with the user; do not re-open without new evidence)

1. Step 2 marker path only. `--groupFile` with more than one model stops.
2. Models grouped by `identical()` sampleID vector; one reader, one decode and
   one QC pass per group. No union sample set, no gather.
3. CLI follows upstream's manifest:
   `--GMMATmodel_varianceRatio_multiTraits_File`, 3 columns, no header.
4. Step 1 = shared ingest + sequential fits. Batched multi-RHS PCG deferred.

## Standing rules for this work

* **Never regenerate or re-baseline a golden.** The design predicts every
  existing golden stays byte-identical; a changed golden is a finding to write
  up, not a file to refresh.
* A failed gate, or a drill that does not fail when it should, ends that phase.
  Record it here and move to the next independent phase.
* No scope widening: no region support, no union sample set, no batched PCG.
* Commit locally only, never push.

## Status

| phase | state |
|---|---|
| Phase 0 (instrument + measure) | **DONE**, gated |
| Phase 1 M1-M4 (step 2 multi-model) | not started (M1 C++ WIP saved to scratch) |
| Phase 2 N1-N4 (step 1 shared ingest) | not started |
| Phase 3 (end-to-end benchmark) | not started |

## Log

### Phase 0 -- instrument and measure (byte-neutral)

What: `SAIGE_PROFILE_MARKER`, an off-by-default three-bucket stopwatch inside
`testOneMarker` (`src/Main.cpp`), serial arm only, printed once per chunk; and
`R/SAIGE_fitGLMM_fast.R` now PRINTS the two step-1 timers it already computed
and discarded (the ingest `system.time` and the one around
`glmmkin.ai_PCG_Rcpp_*`).

Measured at N=50,000 (`SAIGE-sim/harness/work/sim50k`, 210,515 markers,
`--minMAC=20`, 2 reps per trait) -- full table in `harness/bench/NOTES.md`:

| arm | shared `s` | per-model `m` | predicted P=5 | ceiling |
|---|---|---|---|---|
| y_binary | 58-63 s | 28-31 s | **2.16-2.18x** | ~3.1x |
| y_quant | 64-74 s | 5.9 s | **3.74-3.86x** | ~12-14x |

`testOneMarker` is 94.2 s of a 103.9 s run at this N, so these buckets are a
prediction about the whole run rather than about a fraction of it. The entire
binary-vs-quantitative gap is the two `GVec.elem(m_case_indices)` gathers in the
summary bucket (6.5 s vs 0.02 s).

**This is the number Phase 3 is checked against.**

Gate: **PASSED.** `make_goldens.sh harness/work/cand_p0` +
`compare_goldens.sh`: 58 step-0/step-2 files byte-identical, every step-1
comparison `max_rel_err=0.000e+00`, "all comparisons PASS".
`run_tests.sh standard` in BOTH compare modes:
`[ FAIL 0 | WARN 0 | SKIP 4 | PASS 1868 ]` each. Nothing moved, as designed --
the instrumentation is behind a `getenv` and the R change only prints.

Commit: see `git log` for the Phase 0 commit.

### M1 -- the C++ generalisation, with P still fixed at 1

What: `src/Main.cpp` now carries `g_models` (manifest order), `g_readers` (one
per sample-set group) and `g_modelGroup`; `g_model`/`g_reader` became
NON-OWNING aliases of element 0, so the ~80 `g_model->` sites outside the
marker loop compile and mean exactly what they meant. `testOneMarker` split
into `prepareOneMarker` (model-independent: read, both QC gates, impute/flip,
index copies) and `testOneModel` (the S4 body with `g_model->` -> `model.`,
calling `model.getMarkerPval` directly rather than the `Unified_` wrapper,
which dispatches to model 0). Results are `vector<vector<MarkerResult>>`,
model-major. `g_out.singles` is a vector of `unique_ptr<SingleWriter>`;
`openOutfile_single` takes a MODEL INDEX and derives the header from that
model; close/flush loop. `MarkerWorkerState` holds one scratch and one reader
clone PER GROUP. New exports `addSAIGEobjInCPP`, `addGenoReaderInCPP`,
`setSingleOutputPrefixesInCPP` (replacing the caller-less fork vestige
`setSingleOutputPrefixInCPP`); `assign_conditionMarkers_factors` takes a model
index and reads through that model's group.

R touched only where signatures forced it: `openOutfile_single(0L, ...)`,
`mainMarkerInCPP` without `traitType`, `assign_conditionMarkers_factors(...,
0L)`.

`harness/tsan.sh`'s plant anchors were updated for the new worker shape (a
worker state now holds one scratch per group, so the defect is "every worker
uses the first state's group-0 scratch"). Without that the plant would have
died on its own drift assert -- which is exactly what it is for.

Gate: **PASSED.**
* goldens: 58/58 step-0/step-2 byte-identical, all step-1 comparisons PASS.
* standard tier: `[ FAIL 0 | WARN 0 | SKIP 4 | PASS 1868 ]` in BOTH
  `SAIGE_COMPARE_MODE=byte` and `numeric`.
* `harness/tsan.sh clean`: races=0 at t=1/4/8 on all three workloads,
  bytes-vs-serial identical.
* `harness/tsan.sh plant`: races 271/279 (binary), 118/523 (binary_sparse),
  411/546 (quant) and bytes-vs-serial DIFFERS everywhere -- so the clean run
  above is gating something.

### M2 -- the manifest, the R plumbing, the CLI and the tests

What: `--GMMATmodel_varianceRatio_multiTraits_File` (3 columns, no header) on
`SPAGMMATtest` and on `extdata/step2_SPAtests.R`. New `R/MultiModel.R`:
`readModelManifest`, `resolveModelSpecs` (mutual exclusivity + the suffix rule
-- a ONE-ROW manifest writes the unsuffixed file, so golden P1 is a negative
control for the whole feature), `groupModelsBySampleSet` (`identical()` on the
sampleID vector), `addGenoReaderForGroup`, `restartChunkForModels`,
`prepareNullModel` (the per-model half of what `SPAGMMATtest` did inline) and
`pushNullModelsToCPP`. `SAIGE.Marker` takes vectors of output files and index
files and still makes exactly ONE `mainMarkerInCPP` call per chunk.
`R/Geno.R`'s `setGenoInput` now returns the RESOLVED file paths so a second
group's reader is built from the same inputs.

Three refusals for P > 1, all raised BEFORE any C++ call: `--groupFile`,
`--is_rewrite_XnonPAR_forMales`, `--subSampleFile`.

Gate: **PASSED.**
* goldens: 58/58 byte-identical, "all comparisons PASS".
* standard tier, both compare modes: `[ FAIL 6 | PASS 1953 ]` -- and all six
  failures are the three `test-unit-step1-reingest.R` cases, which are
  **Phase 2's** test running against a binary built before Phase 2's C++ fix.
  That is drill **D1** observed for free: without the `setGenoObjFromSource`
  reset, the re-ingest probe does not even complete. Every other test passes,
  including the ten new multi-model ones.
* Two models over one sample set, run by hand at 100 markers: `_binary`
  byte-identical to golden P1 and `_quant` to golden P32, from ONE genotype
  read.

Known-not-done, deliberately: `man/*.Rd` was NOT regenerated. It is already
stale on this branch (`SPAGMMATtest.Rd` still documents `bgenFile`, stripped
long ago), so regenerating would produce a large unrelated diff; the roxygen
`@param` block IS written, so a future regeneration picks it up.

Verification detail worth keeping: the multi-model conditional test first
"failed" against golden P12 because it used `--condition=rs30,rs79`. P12 is
`--condition=rs13,rs79`. Not a bug -- but the shape of that mistake (a golden
compared against a run it does not describe) is exactly what the byte gate is
for, so it is recorded rather than quietly fixed.

## Phase 2 -- step 1 shared ingest

### N1 -- the re-ingest reset (LEDGER #55, FIXED)

What: `setGenoObjFromSource` reset only four members; everything else APPENDS.
A second ingest in one process doubled
`numberofMarkerswithMAFge_minMAFtoConstructGRM`, returned 2M QC flags (wrong
LOCO chromosome ranges), and made variance-ratio sampling index the FIRST
ingest's never-freed `genoVecofPointers_forVarRatio` entries.
`closeGenoFile_plink` cleared none of it and leaked that whole array. Fixed by
zeroing/clearing every accumulator at the TOP of the ingest plus a new
`freeGenoArrays()` that the close also uses. No shipped path re-ingests, so the
fix is byte-neutral BY CONSTRUCTION.

Gate: **PASSED.**
* goldens: 58/58 step-0/step-2 byte-identical; and all TWENTY step-1 files
  (`.rda` + `.varianceRatio.txt`, PLINK and PGEN) byte-identical by `cmp`, which
  is strictly stronger than `compare_goldens.sh`'s toleranced step-1 arm.
* standard tier: `[ FAIL 0 | WARN 0 | SKIP 4 | PASS 1971 ]`.
* **Drill D1 (positive control): PASSED, emphatically.** With the reset
  `#if 0`-ed out and rebuilt, the second ingest reports **17,166** markers
  passing the GRM gate against **8,494** fresh, and then SEGFAULTS in
  `get_GRMdiagVec()` -- `[ FAIL 6 | PASS 0 ]`. The bug is real and the test sees
  it. Fix restored and rebuilt before the gate above.

One test bug worth recording, because it looked exactly like the leak:
`indicatorGenoSamplesWithPheno` is indexed by FILE sample, not by kept sample,
so a length-900 vector for a 1000-row `.fam` is read out of bounds in C++. Being
UB, it read different memory in a fresh process than in one that had already
ingested -- i.e. it produced a difference between the two arms that was NOT the
bug under test. The test now builds the full-length indicator.

### N2 -- the shared-ingest state machine, P = 1 only

What: new `R/SAIGE_fitGLMM_multi.R` holds `.step1Geno` (loaded / key /
nIngests), `step1IngestKey`, `ensureGenoLoaded` and `releaseGeno`. The four
step-1 ingest call sites in `R/SAIGE_fitGLMM_fast.R` (both fitters, both
`skipModelFitting` reload paths) now go through `ensureGenoLoaded`, which
reduces to exactly the old `setgenoFromInput` call whenever nothing is loaded --
i.e. on every single-phenotype run. `fitNULLGLMM` gained
`releaseGenoAtEnd = TRUE` and now RETURNS, invisibly, the paths it actually
wrote (`modelOut` is reassigned under `--isLowMemLOCO`, and `outputPrefix` is
mutated in place by `--FemaleOnly/--MaleOnly`), which is what the driver needs
for the manifest. `ingestGates` (minMAF, maxMissingRate, and the VR MAC bounds)
is part of the ingest's identity, because those setters are pushed BEFORE
`setgeno` and decide which markers it keeps.

The reuse test is `identical()` on the whole key: anything not provably the same
ingest gets a fresh one, so reuse can only ever be a saving and never a
difference.

Gate: **PASSED.** goldens 58/58 byte-identical; **all 20 step-1 files
byte-identical by `cmp`**; exactly ONE `setgeno mark2` per golden step-1 run
(the ingest is not being repeated); standard tier
`[ FAIL 0 | WARN 0 | SKIP 4 | PASS 1971 ]` in BOTH compare modes; lint clean.

### N3 -- the driver, the CLI, and the step-1 multi-phenotype tests

What: `fitNULLGLMM_multi(phenoCol, traitType, invNormalize, ...)` in
`R/SAIGE_fitGLMM_multi.R`. One phenotype IS `fitNULLGLMM` -- legacy names, no
manifest. Several are grouped by complete-case sample set (`step1PhenoGroups`,
ORDERING ONLY: correctness rests on `ensureGenoLoaded`'s key), fitted
sequentially with `set.seed(seed)` before every fit, named
`<outputPrefix>_<phenoCol>`, and finished with a manifest
`<outputPrefix>.multiTraits.txt` that step 2 consumes directly.
`extdata/step1_fitNULLGLMM.R` accepts comma-separated `--phenoCol`,
`--traitType` and `--invNormalize`.

**The option_list metadata diff is EXACTLY what the plan predicted**, verified
by `identical()` over the parsed old and new `option_list` objects, all 57
options: three help strings, plus `--invNormalize` moving from
`type="logical", default=FALSE` to `type="character", default="FALSE"`. Nothing
else, including every `action`, which is `store` for all 57 -- so
`--invNormalize=TRUE` was always the required form and no existing invocation
changes meaning.

Gate: **PASSED.**
* goldens 58/58 byte-identical; all 20 step-1 files byte-identical by `cmp`.
* standard tier: `[ FAIL 0 | SKIP 4 | PASS 2017 ]` byte,
  `[ FAIL 0 | SKIP 4 | PASS 2192 ]` numeric.
* Drill **D2** (drop the per-fit `set.seed`): the TWO-GROUP test's
  `varianceRatio.txt` stops matching its standalone run
  (`[ FAIL 1 | PASS 14 ]`) while every single-group test stays green
  (`[ FAIL 0 | PASS 31 ]`) -- which is the point of the drill: the single-group
  tests could not have caught it.
* Drill **D4** (force `ensureGenoLoaded` to always re-ingest): every byte
  comparison still passes and only the two "exactly one `setgeno mark2`"
  assertions fail (`[ FAIL 2 | PASS 29 ]`), separating "shared ingest" from
  "correct bytes" -- the LEDGER #30 shape.

One comparator finding: `expect_output_matches_golden` is a step-2 TEXT-table
comparator, so in numeric mode it fed a binary `.rda` to `read.table`
("embedded nulls", "invalid UTF-8"). New `tests/testthat/helper-step1multi.R`
compares step-1 outputs properly -- byte-identical in byte mode; in numeric
mode, `.rda` field by field at 1e-5 relative over the ten sample-indexed numeric
fields (names and sampleID exactly), and `.varianceRatio.txt` with the ratio
numeric at 1e-5 and the type/category LABELS exact, because a run that put the
right number in the wrong bin is not an equivalent run.

## Phase 3 -- end-to-end benchmark on the N=50,000 fixture

Five traits, `--nThreads=1`, paired arms back to back in one batch, user CPU.
Full tables in `harness/bench/NOTES.md`.

| arm | serial | multi | speedup | outputs |
|---|---|---|---|---|
| **step 2**, 210,515 markers | 412.28 s | **137.86 s** | **2.99x** | **all five BYTE-IDENTICAL** |
| **step 1**, sparse GRM + LOCO | 61.60 s | **37.42 s** | **1.65x** | all five `.rda` BYTE-IDENTICAL |

Wall clock for step 2: 481.9 s -> 146.9 s (3.28x). Peak RSS for the five-model
scan is 2843 MB -- essentially ONE serial run's (2174-2783 MB), not five.

**Phase 0 predicted 2.85x for this 2-binary/3-quantitative mix; measured 2.99x**
(5% high, the residual being per-process startup the buckets never counted).
Prediction and measurement agree.

**The one thing that is NOT byte-identical, and why it is not this feature's
doing.** Step-1 `varianceRatio.txt` differs for phenotypes 2-5 at up to 1.1e-7
relative, confined to MAC category 2. Measured against its own control: the SAME
standalone command run in two different batches differs from ITSELF by up to
1.0e-7, in the same category. Three controls establish it:
1. three standalone `y_quant` runs in one batch reproduce each other exactly;
2. the same standalone command in a different batch does NOT reproduce itself
   (`y_binary` 5.5e-8, `y_null` 1.0e-7);
3. forcing a fresh ingest per phenotype (shared PROCESS, no shared ingest) does
   not remove the difference -- it changes WHICH phenotypes differ, which is
   what noise does and a systematic effect does not.

**CORRECTION (same day).** Control 1 above was drawn from n = 3 and its
conclusion was wrong. Eleven runs of ONE identical command -- same binary, back
to back, plus a rebuild with no source change, plus two under 8-way load --
produced **six distinct `varianceRatio.txt` files**, disagreeing WITHIN the
batch as well as across batches. There is no "stable within a batch"; the
command is simply not reproducible on this configuration. Recorded as
**LEDGER #56**. Per row over the 11 runs: category 1 has ONE distinct value
(spread 0.000e+00), category 2 has six (`null_noXadj 2` spread 1.5e-7), and all
eleven `.rda` files are byte-identical.

This makes the Phase 3 conclusion stronger -- the multi-vs-serial difference of
1.1e-7 sits inside the identical command's own 1.5e-7 spread, in the same rows
-- but the reasoning that first reached it was faulty and is replaced.
Control 3 (a fresh ingest per phenotype does not remove the difference, it
changes which phenotypes show it) still stands and is now the load-bearing one.

At N=1,000 with `--useSparseGRMforVarRatio=FALSE` (the test suite's fixture) the
files ARE byte-identical in both phenotype orders, which is what the gated tests
assert. **No golden was regenerated and none changed.**

## Status: all four phases complete

| phase | state |
|---|---|
| Phase 0 (instrument + measure) | **DONE**, gated, `705d469e` |
| Phase 1 M1-M4 (step 2 multi-model) | **DONE**, gated, `9962dc60` `e2810821` `19ab5097` |
| Phase 2 N1-N4 (step 1 shared ingest) | **DONE**, gated, `75ce9067` `79dc672d` `999143d3` `c9647f3c` |
| Phase 3 (end-to-end benchmark) | **DONE** |

### Left undone, deliberately

* `man/*.Rd` not regenerated (already stale on this branch, pre-strip).
* Optional golden P35 (`pgen_marker_quantitative.txt`) from M3's list: the
  multi-model PGEN test reconciles against an in-test single-model run instead,
  so the gap it would fill is covered without adding a golden.
* Everything in the plan's "deliberately out of scope" list: region tests with
  P > 1, union sample sets, batched multi-RHS PCG, an O(MAC) case/control
  summary, male-X and `--subSampleFile` with P > 1, restart resume-at-min.
