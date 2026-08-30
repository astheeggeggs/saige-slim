# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**This is the `strip` branch**: SAIGE reduced to step 0 (createSparseGRM), step 1
(fitNULLGLMM), and step 2 (SPAGMMATtest: single-variant + region/gene tests),
with **PLINK bed/bim/fam as the only genotype input**. The BGEN/VCF/SAV/PGEN
readers, step 3 (LD matrices), and RareEffect were stripped in commits
"Strip 1"–"Strip 7", each verified byte-exact against `main`-HEAD goldens.
`main` holds the full upstream fork (docs: https://saigegit.github.io/SAIGE-doc/).
Binary, quantitative, and survival traits are all still supported (survival's
fate is a refactor-phase decision).

**The working plan lives at `~/.claude/plans/rosy-munching-pudding.md`**:
Phase 0 (baseline/harness) and Phase 1 (strip) are done; next are Phase 2
(tests), Phase 3 (refactor; survival decision; reentrancy), Phase 4 (speedup).

## Rules of this branch

- **Bug-for-bug preservation.** Known bugs are recorded in `LEDGER.md`, NOT
  fixed. Fixes happen in a dedicated post-refactor phase, one commit per ledger
  entry, with regenerated goldens. Read `LEDGER.md` before "fixing" anything.
- **Every change is verified against the goldens** (see Harness below). Step-2
  and step-0 outputs must stay byte-identical; step-1 within tolerances.
- After adding/removing any `[[Rcpp::export]]` or `.cpp` file: update the
  explicit `OBJECTS` list in `src/Makevars` and re-run
  `pixi run Rscript -e 'Rcpp::compileAttributes(".")'`.

## Build

```sh
pixi install
# lintools is a CRAN LinkingTo dep (linear systems, not a linter), not in conda:
pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")'
pixi run R CMD INSTALL .
```

Works natively on macOS and linux (no plink-ng step anymore — it died with the
PGEN reader). `src/Makevars` still hardcodes `-I../.pixi/envs/default/include`,
so the in-tree `.pixi/` env is required. `make` at the repo root runs the
Docker build.

## Harness (the test suite until Phase 2 formalizes it)

- `source harness/env.sh` before ANY run — pins
  `RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.
  Without it, full-GRM step-1 fits are nondeterministic (TBB float
  `parallel_reduce`; empirically confirmed — see LEDGER.md #7). Sparse-GRM
  fits are byte-stable either way.
- `bash harness/make_goldens.sh <dir>` — regenerates the 12 runs (step0 sparse
  GRM; step1 binary ×2 / quantitative / survival; step2 PLINK marker, group
  fastTest, group conditional, marker conditional). Step-2 runs consume frozen
  model inputs: our committed step-1 goldens plus
  `extdata/output/example_binary_sparseGRM.{rda,varianceRatio.txt}` (the only
  legacy committed artifact still HEAD-compatible — the 2022-era ones crash,
  LEDGER.md #10).
- `bash harness/compare_goldens.sh harness/goldens <candidate-dir>` — byte
  compare for step0/step2; toleranced step1 (θ/varianceRatio 1e-2 relative —
  AI-REML `tol=0.02` makes tighter meaningless; vectors 1e-3).
- `bash harness/make_controls.sh <dir>` — permuted-phenotype negative control
  through the full pipeline + KS uniformity check (`harness/check_null.R`).
- Baseline goldens (committed): `harness/goldens/`. Retired bgen/vcf goldens:
  `harness/goldens_retired/`. Full verification cycle:
  `compileAttributes && R CMD INSTALL . && make_goldens.sh harness/work/candidate && compare_goldens.sh`.
- PLINK group files key on **bim column 2 (rs IDs)** — use
  `extdata/input/group_new_snpid.txt`; `--condition` takes rsIDs too.

## Architecture

### Entry points → R functions

| Script (extdata/) | Function | Role |
|---|---|---|
| `createSparseGRM.R` | `createSparseGRM()` in `R/SAIGE_createSparseGRM.R` | Step 0: sparse GRM (`.mtx` + `.sampleIDs.txt`) |
| `step1_fitNULLGLMM.R` | `fitNULLGLMM()` in `R/SAIGE_fitGLMM_fast.R` | Step 1: null GLMM (AI-REML + PCG) + variance ratios |
| `step2_SPAtests.R` | `SPAGMMATtest()` in `R/SAIGE_Test_main.R` | Step 2: single-variant or region tests (branches on `groupFile`) |

### R ↔ C++: configure-then-run over global state

All C++ via Rcpp (`RcppExports` generated; `NAMESPACE` uses `exportPattern`).
R pushes state into C++ **globals** via `set*InCPP()` calls, then invokes one
long-running `main*InCPP()` loop that writes output files from C++
(`std::ofstream`s in `Main.cpp`; R only writes the `.index` progress file).

- `src/Main.cpp` — step-2 bridge: static `ptr_gPLINKobj`/`ptr_gSAIGEobj`
  singletons, ~25 `g_*` settings, `mainMarkerInCPP`, `mainRegionInCPP`,
  `Unified_getOneMarker`. Forgetting a `set*_GlobalVarsInCPP` call silently
  leaves stale values. `Main.hpp` declarations are stale (LEDGER.md #5) —
  don't trust them as a contract.
- `src/SAIGE_fitGLMM_fast.cpp` (~7400 lines) — step-1 engine, separate from
  Main.cpp. Single **global** `genoClass geno` holds the bit-packed GRM
  genotypes; there is a second, independent PLINK reader (`src/PLINK.cpp` /
  `PlinkClass`) for step 2. Algorithms exist as 4–6 near-duplicates: base,
  `_LOCO`, `_q` (quantitative), `_Surv` (survival) — survival removal is
  ~110 edit sites across 8 files (deferred to refactor phase by design).
- `src/SAIGE_test.cpp` — `SAIGEClass`: `getMarkerPval` cascade (score → SPA if
  |Z| > cutoff → ER exact test if binary & MAC ≤ `max_MAC_for_ER` → Firth),
  per-MAC variance-ratio selection. `src/Binary_*`/`src/ER_binary_func.cpp`
  are the live ER machinery (imported from SKAT) — keep.
- Step 0 cannot be separated from step 1: both drive the same global `geno`,
  and `createSparseKinParallel` is defined in `R/SAIGE_fitGLMM_fast.R`.
- One dataset per R process; no reentrancy. Parallelism is process-level
  (`step2_SPAtests.R` splits markers and `mclapply`s) plus RcppParallel/TBB
  in step-1 GRM kernels (source of LEDGER.md #7).

### The step1 → step2 contract (implicit, name-based)

- `<prefix>.rda` must contain a variable literally named `modglmm`;
  `R/readInGLMM.R::ReadModel` `load()`s it blind and mutates it (no schema
  or version field). Non-LOCO fits go to `<prefix>_noLOCO.rda`;
  `--isLowMemLOCO` writes per-chromosome `<prefix>_chr<j>.rda`.
- Variance-ratio file: headerless `value type category` rows, `type` ∈
  `{sparse, null, null_noXadj}` (`Get_Variance_Ratio` in `R/readInGLMM.R`).
  Step-1 sparse-GRM fits emit only `null`/`null_noXadj` rows → step 2 must
  NOT be given `--sparseGRMFile` for those models; group tests need a model
  whose VR file has `sparse` rows. Multiple rows per type = categorical VRs
  matched to `cateVarRatioMinMACVecExclude`; C++ picks per marker by MAC.
- LOCO is three mechanisms: step-1 global-state mutation per chromosome, the
  save-time file layout, and step-2 `ReadModel` overwriting fitted values from
  `LOCOResult[[chrom]]` (`--chrom` mandatory when LOCO=TRUE).
