# Handoff 2026-09-14 (night): Phase G -- the god functions, as prescribed

Phase G of the 2026-09-14 audit (`docs/state/handoff-2026-09-14-codebase-audit.md`
section 6 and the "Phase G" block), taken exactly as written there: G1 in
five gated slices, G2 as the mechanical type hoist and nothing more. Every
gate was `prepush.sh` at `both` modes; every prediction before every gate
was "no golden moves", and none did. Line numbers below are at the commits
named, not the audit's (which were at `9deab7c9`).

## G1 -- `SAIGE.Region`, five slices (R/RegionState.R is new, lint-clean)

| slice | commit | what left the body | prediction argument |
|---|---|---|---|
| 1 `regionStateFor` | `6db4fd1e` | the ~100-line per-region state record: WEIGHT canonicalisation by set count, the negative-weight guard, the "Analyzing Region" line, SCORE/ISLOF, the two degenerate-score warnings | moved verbatim; stdout/stop/warnings unchanged |
| 2 `regionTestPlan` | `ac623215` | r.corr -> method, r.corr, regionTestType, is_single_in_groupTest | RETURNS r.corr, which `SKAT_Check_Method` promotes to the 7-element rho grid every SKAT-O consumer reads; dropping it would move every SKAT-O golden |
| 3 `regionDoneMessage` | `ac623215` | the return string; the dead first `.markerInfo` message deleted | a return value the CLI prints; rode on slice 2's gate |
| 4 `composeWeightSets` | `48833f1a` | custom -> transform -> beta -> No_Weights composition, the two refusals, the three stdout lines | scoped to composition only: the restart schema guard, the `openOutfile` loop, `nline_per_gene`, `region_weight_names` stay (I/O and C++ state) |
| 5 `scan <- new.env()` | `2162c88a` | the seven superassigned accumulators + `cth_chunk_to_output`, `iswriteRegion`, `iswriteMarkerList`; `OutList.all` (written twice, never read) deleted | a binding rename, 67 lines, whole-word outside string literals; an environment so no stage forces a rebind (the #81 shape) |

Not done, by the audit's instruction: the chunk writer is NOT unified (three
distinct `(append, col.names)` rules, #59's asymmetry deliberate, and no
test combines the marker list with a multi-chunk run); the file-opening
sites are NOT reordered (the SKAT-O pool forks before the first C++ call;
the BURDEN `openOutfile` needs `nWeightSetsUsed`); `regionRestartPlan` /
`writeRegionIndex` untouched. `SAIGE.Region` went from 1,225 lines (at `05e567a0`) to
1,041; what is left is orchestration and I/O, which is what the
audit said it would be.

## G2 -- `mainRegionInCPP`, the type hoist only

Twelve function-local types (the audit's ten plus `WeightSet` and
`CellResult`, which Phase E added after the count) moved VERBATIM, each with
its own introductory comment block, into `src/RegionLoop.hpp` (namespace
`Region`); `mainRegionInCPP` opens with `using namespace Region;` so every
spelling in the body is unchanged. Main.cpp includes the header right
before the function, AFTER its file-scope helpers (`MarkerStatePool`,
`MarkerWorkerState`, `RefuseMainThreadGuard`, `RegionChunkStore`, the
TSAN macros), which the types use as complete types: the header says on its
first lines that it is a section of Main.cpp lifted out and is not
includable on its own. Pass-level `// ----` section headers stayed in
Main.cpp. Main.cpp 5,800 -> 5,567 lines. Refused, as the audit refused
them: a `Work` struct, a phase lift, a new translation unit (armadillo's
partial-accumulator unrolling is a 1-ULP hazard for the printed goldens).

Because two of the types are `RcppParallel::Worker`s, the region TSAN
workload was run in both arms after the gate -- see the G2 commit for the
numbers. Warning ratchet unchanged at 92.

## What Phase G did not change

No behaviour. No golden. No test text. No CLI output. The one visible
difference is the shape of the two files, which is the point: the next
feature has somewhere to land beside the frame instead of inside it.
