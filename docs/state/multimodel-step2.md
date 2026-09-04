# Multi-model single-variant scanning (Phase 1, commits M1-M4)

**What it is.** `SPAGMMATtest` can test several null models against one
genotype read. One process, one pass over the `.bed`/`.pgen`, one decode per
marker per SAMPLE SET, and one output file per model. It is the step-2 half of
the multi-phenotype work; the step-1 half is `docs/state/step1-multipheno.md`.

**The claim that makes it checkable.** Each model's output file is
BYTE-IDENTICAL to what that model's own single-model run produces. Not
"equivalent at a tolerance" -- identical. That is not a comparison between two
code paths (the standing rule against which still holds); it is one code path,
`prepare once then test each model`, checked against the goldens that same path
produced when there was only ever one model.

## The CLI

```
--GMMATmodel_varianceRatio_multiTraits_File=manifest.txt
```

Three whitespace-separated columns, no header, one row per model:

```
y_binary   /path/step1_y_binary.rda    /path/step1_y_binary.varianceRatio.txt
y_quant    /path/step1_y_quant.rda     /path/step1_y_quant.varianceRatio.txt
```

Mutually exclusive with `--GMMATmodelFile`/`--varianceRatioFile`. Output goes to
`<SAIGEOutputFile>_<name>` per model, each with its own `.index`. **A ONE-ROW
manifest writes the unsuffixed `<SAIGEOutputFile>`** and is byte- and
name-identical to the equivalent `--GMMATmodelFile` run -- upstream's rule,
kept because it makes golden P1 a negative control for the entire feature.

Step 1 writes this file itself when given `--phenoCol=y1,y2`; see
`docs/state/step1-multipheno.md`.

## Sample sets: grouping, not unioning

Models whose step-1 `sampleID` vectors are `identical()` form a GROUP. Each
group gets one `GenoReader`, and per marker each group does one read, one
missing-index scan, one `imputeGenoAndFlip` and one pair of index copies; each
MODEL in the group then runs its own score test over that one prepared buffer.

`identical()` rather than `setequal()` on purpose: the order of `sampleID` is
the reader's output order, so a permutation is a different group and gets its
own reader. Being conservative here is what makes the byte-identity claim true
by construction rather than by argument.

**Upstream's `multitrait_1.4.8.1` branch unions the sample sets instead**, and
that is the decision this design deliberately does not copy: a union means
per-marker QC, missing rate and N are computed over samples that are not in the
model, and it needs a per-model recount that would have to replay the reader's
accumulation order (PGEN sums dosages serially, `src/PGEN.cpp:471-483`). See
the plan's appendix for the full list of that branch's defects.

## What is refused, and why

| refused for P > 1 | why |
|---|---|
| `--groupFile` (region tests) | escalation is a whole-region C++ re-run keyed to `g_regionEscalatedPass`; per-model escalation needs the region loop restructured, and ~90% of region time is per-model work anyway (LEDGER #51) |
| `--is_rewrite_XnonPAR_forMales` | `g_cfg.indexInModel_male` is ONE vector of MODEL-order indices in run config, applied inside every read |
| `--subSampleFile` | applies one subsample to one model |

All three stop in R **before any state is pushed into C++**, which is why
`tests/testthat/test-unit-multimodel.R` can assert them in-process.

## Where it lives

**C++ (`src/Main.cpp`).** `g_models` owns the models in manifest order,
`g_readers` owns one reader per group, `g_modelGroup` maps model to group.
`g_model` and `g_reader` are NON-OWNING aliases of element 0, so the region
path, the S5 test hooks and `set_flagSparseGRM_cur_SAIGE_org` keep compiling
and keep meaning what they meant -- which is exactly why P > 1 is refused for
those paths rather than half-supported.

`testOneMarker` split in two:

* `prepareOneMarker(reader, ..., n, scratch, rows, MarkerPrep&)` -- the read,
  both QC gates, the impute/flip and the index copies. Returns -1 short read /
  0 QC-dropped / 1 ready. Marker fields go into EVERY row of the group; they are
  duplicated per row rather than shared because two models in different groups
  legitimately disagree about AC/AF/missingRate/N.
* `testOneModel(model, ..., prep, scratch, res)` -- the S4 body with `g_model->`
  become `model.`. It calls `model.getMarkerPval` directly rather than the
  `Unified_getMarkerPval` wrapper, which dispatches to `g_models[0]`.
  `getMarkerPval` is `const` and never writes `t_GVec` (S3/S12), which is what
  lets every model in a group test the same buffer with no copy.

`MarkerWorkerState` holds one `MarkerScratch` and one reader clone PER GROUP --
per group, not per model, because a group shares the decode that lands in the
buffer. Post-pass 1 rethrows in serial order (marker, group, model); post-pass 2
truncates EVERY model from the earliest short read (one loop, however many
models); post-pass 3 re-runs only the deferred models of a group on the main
thread, which is byte-neutral because prepare is deterministic.

Output is `std::vector<std::unique_ptr<SingleWriter>> g_out.singles` --
unique_ptr because an open `std::ofstream` must never be moved and a vector
reallocating on push_back would do exactly that.

**R.** `R/MultiModel.R` holds the manifest reader, `resolveModelSpecs` (the
mutual-exclusivity stop and the suffix rule), `groupModelsBySampleSet`,
`addGenoReaderForGroup`, `restartChunkForModels`, `prepareNullModel` (the
per-model half of what `SPAGMMATtest` used to do inline) and
`pushNullModelsToCPP`. `R/SAIGE_SPATest_Marker.R` takes vectors of output files
and indexes; there is still exactly ONE `mainMarkerInCPP` call per chunk,
because the models are the inner loop.

## Restart

Each model has its own `.index` carrying its own byte count, because two models
write different numbers of bytes per chunk. A restart truncates each file to ITS
OWN recorded size (LEDGER #47's fix, per file). If the index files DISAGREE
about the chunk to resume at, the run stops and names them: a disagreement means
they were not all written by the same interrupted run, and guessing which is
authoritative is how rows get silently duplicated.

## What is NOT gated

Only stdout ordering moved: `prepareNullModel` prints the sparse-GRM and
variance-ratio lines before the male-X block rather than after it, because the
models are now prepared as a unit. No output file is affected, and the golden
gate is on files.
