# Handoff: step 4 -- several null models per region (2026-09-06)

Two commits on `strip` after `e7b346ba` (step 3b), each gated. Plan:
`~/.claude/plans/can-you-summarise-what-shimmering-snail.md` (steps 1-4 done,
5 not started). Previous handoff: `docs/state/handoff-2026-09-06-region-cells.md`.
The marker-path design this mirrors: `docs/state/multimodel-step2.md`.

## What landed

**4a, C++.** `mainRegionInCPP` takes `t_modelIdx` (0-based indices into
`g_models`) and returns one output list per model (a list of lists; R unwraps
a single model's `[[1]]`). Everything that depends on the MODEL is in a
`RegionModelWork` per model: the per-marker tests (`RegionPrep.m[k]`), the
P1/P2 scratch pair and chunk store, VarMat, the residualised `genoSumMat`
(a copy of the shared accumulation, residualised in place per model at the
ultra-rare pass), the per-index outputs, the case/control counts, the
ultra-rare and BURDEN cell tests (`URCellModel` / `BurdenCellModel` beside the
shared `URCell` / `BurdenCell`), the BURDEN emission, the output files and the
exported list. Everything that depends only on the GENOTYPES is computed once:
read, impute, QC, weights, the ultra-rare collapse, the MAC/count
accumulators, the chunk bookkeeping. Pass A prepares each marker once and
tests it for every model on the same worker; the cell tests run one (cell,
model) per task. Models in one call must share a sample set and agree on the
conditioning set-up and the sparse-GRM flag (checked at entry).
`Step2Output.regions` is one `RegionWriter` per model; `openOutfile`,
`openOutfile_singleinGroup`, `copy_singleInGroup`,
`RegionSetUpConditional_binary_InCPP` and
`assign_conditionMarkers_factors_binary_region` take the model index. Every
region call site of `getMarkerPval` now calls the model directly with the
worker's `KernelWorkspace`. A single model is index 0 at the run prefix:
byte-identical files and names.

**4b, R.** A model manifest (`--GMMATmodel_varianceRatio_multiTraits_File`)
is accepted on a region run when every model was fitted on the SAME sample
set (`identical()` step-1 `sampleID` vectors, the marker path's grouping
rule); otherwise it stops before any state reaches C++. `SAIGE.Region` takes
a list of fitted-value vectors and vectors of output files and trait types;
one C++ call per region tests every model, each model is post-processed with
its own `mu`, escalation is per model (SKAT-O: the models whose Cauchy p
crosses the cutoff are re-run together; BURDEN: the models C++ did not write
on the fast pass), non-escalated models copy their single-in-group rows, and
each model writes its own result file, `.markerList.txt` and `.index`
(restart requires the index files to agree, as on the marker path). The
conditioning scale factors are computed per binary model. A one-row manifest
writes the unsuffixed file.

## Gate

4a, commit `1e02a154` (`harness/work/region_profile/gate_s4a.log`): the
region test files 38 / 15 / 13 / 98 green (goldens at 4 and 8 threads,
escalation, chunk spill, conditional multi-weight); `sim50k` SKAT-O and
BURDEN at 1 and 4 threads on both traits, 24/24 files `cmp`-identical to
step 3's outputs; goldens all PASS; `prepush.sh` FAIL 0 both modes (2769 /
2969); `tsan.sh` region workload clean 0 / 0 at 4 / 8 threads, plant 166 / 96
races with differing bytes.

4b: `test-step2-region-multimodel.R` 93/93 -- two models (binary + quant, one
sample set) on SKAT-O with fastTest escalation, quantitative-first as well,
and on BURDEN at 4 threads, each model's region file, single-in-group file and
marker list byte-identical to its own single-model run; a one-row manifest
reproduces P2 and P5 unsuffixed; models on different sample sets (a
manufactured 990-sample model) refused before C++; and the deterministic
escalation-subset test (below). `test-step2-multimodel.R` and
`test-unit-multimodel.R` inverted where they asserted the old refusal.
`sim50k`, five-trait manifest vs five serial runs, 1 thread: 15/15 + 15/15
files identical; SKAT-O 788 s vs 934 s (1.19x), BURDEN 97 s vs 178 s (1.83x)
(`harness/bench/NOTES.md`, last section). Goldens all PASS. `prepush.sh` FAIL 0 both modes (2863 / 3063).

**Two defects `sim50k` caught that the bundled fixture cannot** (both need an
escalation subset that is a strict, non-leading part of the manifest):
`pvalRegionK[[k]] <- NULL` deleted the list element and shifted the models
after it; the escalation call handed C++ the run's first trait type while the
subset's first model was quantitative. The test file now forces such a subset
in seconds -- `--pval_cutoff_for_fastTest` between the two models' fast-pass
Cauchy p for GENE1, both manifest orders, the "Non-fast test is performed"
count asserting which models escalated -- and is red against a scratch
install with the two defects re-introduced. The lesson for the next
multi-model change: run THAT before any hour-long comparison; the `chr20`
simulation config (10,000 samples, 550 genes, five traits; not built on this
machine) is the middle tier if one is wanted.

## State for the next session (written for a context clear)

- Commits on `strip`: `e7b346ba` (step 3b), `1e02a154` (4a), then 4b (this
  handoff's commit). `docs/state/next.md`'s pointer line names this file; the
  rest of that file has another session's uncommitted section -- do not
  commit it wholesale (`harness/work/region_profile/drivers/stage_next_pointer.sh`
  stages the pointer line alone).
- Nothing is running. The installed library is the 4b build (R-only
  reinstalls after 4a's `build.sh`); `git status` should show only `next.md`.
- Drivers, logs and results of the day: `harness/work/region_profile/`
  (`drivers/README.md` lists them). The sim50k fixture is in the `SAIGE-sim`
  worktree (`~/.claude` memory `saige-bench50k-fixture`).
- Open in `next.md`: step 5 (batch the SKAT-O cells across regions); thread
  across REGIONS for BURDEN at scale; a fixture that reaches the cell worker's
  ER deferral; the `chr20` middle-tier fixture.

## Not done / caveats

- Models on DIFFERENT sample sets in one region run are refused (the marker
  path groups them and reads per group; the region path would need one
  genotype pass per group -- a later step if it is ever wanted).
- Mixed trait types in one manifest are fine: `SAIGE_Test_main.R` zeroes
  `max_MAC_use_ER` when model 1 is quantitative, but only AFTER the cutoff has
  reached C++ (`setAssocTest_GlobalVarsInCPP`) and nothing reads the zeroed
  value, and C++ gates the exact test on each model's own trait type. Checked
  while writing this; noted so nobody "fixes" the dead assignment into a live
  one.
- The SKAT-failure sidecar (LEDGER #57) counts over all models and is written
  beside the FIRST model's file (the per-cell record has no model column).
- Step 5 (batch the SKAT-O cells across regions) as the plan describes it.
