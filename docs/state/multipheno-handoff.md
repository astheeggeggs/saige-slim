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
