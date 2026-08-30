# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

SAIGE (v1.5.2): an R package with a large C++/RcppArmadillo core for genome-wide
association tests in biobank-scale data. Covers single-variant SPA tests under a GLMM
(SAIGE), gene/region-based rare-variant tests (SAIGE-GENE+), LD matrices for
meta-analysis (step 3), and per-gene effect-size estimation (RareEffect). Supports
PLINK bed, PGEN, BGEN, VCF/BCF/SAV inputs. User documentation:
https://saigegit.github.io/SAIGE-doc/ — the optparse help strings in `extdata/*.R`
are the de-facto CLI reference (`man/` has only 4 Rd files).

## Build

Dependencies are managed with **pixi** (conda-forge/bioconda: superlu, savvy, boost,
openblas, zstd, R packages). `docker/Dockerfile` is the canonical, working build
recipe; replicate its steps for a local build:

```sh
pixi install
pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")'
# plink-ng 2.0 headers must be compiled into a static lib the linker can find:
curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx
mv plink-ng-2.0.0-a.6.16 plink-ng
pixi run <cc> -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib \
  -o libplink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate
mv libplink2_includes.a .pixi/envs/default/lib
pixi run R CMD INSTALL .
```

- `src/Makevars` hardcodes `-I../.pixi/envs/default/include -I../plink-ng/2.0/include`
  and `-L.. -lplink2_includes` — the in-tree `.pixi/` env and `plink-ng/` checkout are
  required for `R CMD INSTALL` to succeed (`libplink2_includes.a` works at the repo
  root or in `.pixi/envs/default/lib`, where the Dockerfile puts it).
- `src/Makevars` has an explicit `OBJECTS` list: **adding a new .cpp requires editing it**.
- `make` at the repo root just runs the Docker build (`docker build -f docker/Dockerfile .`).
- `lintools` is a CRAN LinkingTo dependency (linear systems, not a linter) installed
  outside pixi.
- `src/SAIGE/` is a checked-in *installed package tree* (build artifact, incl. a 13 MB
  .so) — not source. `src/script_for_vcf/` is dead legacy code not in `OBJECTS`.

## Testing

**There is no automated test suite** — no tests/, testthat, or CI. Validation is
hand-run example invocations in `extdata/`, executed with cwd = `extdata/` against
`extdata/input/` fixtures, with committed reference outputs in `extdata/output/`:

- `extdata/cmd.sh` — ~10 step1/step2 invocations (bgen/vcf/plink, group tests,
  conditional analysis, survival)
- `extdata/cmd_LDmat.sh`, `extdata/cmd_extractNeff.sh`, `extdata/cmd_test_sparseGRM.sh`

These are comment-interleaved command collections, not scripts to run wholesale; to
run "a single test", copy one `Rscript step1_fitNULLGLMM.R ...` / `step2_SPAtests.R ...`
block from `cmd.sh` (prefix with `pixi run` outside Docker). The Dockerfile's smoke
test is `createSparseGRM.R --help && step1_fitNULLGLMM.R --help &&
step2_SPAtests.R --help && step3_LDmat.R --help`. After a change, diff fresh outputs
against the committed reference files in `extdata/output/`.

## Architecture

### Entry points → R functions

Thin optparse wrappers in `extdata/`, each around one exported function:

| Script | Function | Role |
|---|---|---|
| `createSparseGRM.R` | `createSparseGRM()` in `R/SAIGE_createSparseGRM.R` | Step 0: sparse GRM (`.mtx` + `.sampleIDs.txt`) |
| `step1_fitNULLGLMM.R` | `fitNULLGLMM()` in `R/SAIGE_fitGLMM_fast.R` | Step 1: null GLMM (AI-REML + PCG) + variance ratios |
| `step2_SPAtests.R` | `SPAGMMATtest()` in `R/SAIGE_Test_main.R` | Step 2: single-variant or region tests (branches on `groupFile`) |
| `step3_LDmat.R` | `generate_LDMat_forMataRegion()` in `R/SAIGE_SPATest_Region_LDMat.R` | LD matrices for meta-analysis |
| `run_RareEffect.R` | `run_RareEffect()` in `R/RareEffect.R` | Post-step-1 per-gene effect sizes |

### R ↔ C++: configure-then-run over global state

All C++ is reached via Rcpp (`RcppExports` are generated; `NAMESPACE` uses
`exportPattern`, so ~everything is exported). The protocol is not functional: R pushes
state into C++ **globals** via `set*InCPP()` calls, then invokes one long-running
`main*InCPP()` loop that writes output files itself from C++ (`std::ofstream`s in
`Main.cpp`; R only writes the `.index` progress file).

- `src/Main.cpp` — step-2 bridge. File-static singleton pointers (one reader per
  format + `ptr_gSAIGEobj`), ~25 `g_*` analysis settings, `mainMarkerInCPP`,
  `mainRegionInCPP`, and `Unified_getOneMarker` which dispatches on the string
  `"plink"|"pgen"|"bgen"|"vcf"` (no base class; each reader's `getOneMarker` has a
  slightly different signature). Forgetting a `set*_GlobalVarsInCPP` call silently
  leaves stale values.
- `src/SAIGE_fitGLMM_fast.cpp` (~7400 lines) — step-1 engine, separate from Main.cpp.
  A single **global** `genoClass geno` object holds the bit-packed GRM genotypes; all
  step-1 exports read/write it implicitly. Most algorithms exist as 4–6 near-duplicate
  variants: base, `_LOCO`, `_q` (quantitative), `_Surv` (survival), and combinations.
- `src/SAIGE_test.cpp` — `SAIGEClass`: score tests, `getMarkerPval` cascade
  (score → SPA if |Z| > cutoff → ER exact test if binary & MAC ≤ `max_MAC_for_ER`
  → Firth if p < cutoff), per-MAC variance-ratio selection.
- Readers: `src/PLINK.cpp`, `src/PGEN.cpp` (plink-ng), `src/BGEN.cpp`, `src/VCF.cpp`
  (VCF/BCF/SAV via savvy; also `--vcfFilters`).
- SPA: `src/SPA.cpp` dispatches to `src/SPA_binary.cpp` / `src/SPA_survival.cpp`.
- `src/Binary_*`, `src/ER_binary_func.cpp` — SKAT/exact-resampling tests for
  ultra-rare binary variants (imported from SKAT/MetaSKAT).

Consequences: **one dataset per R process**, no reentrancy. Parallelism is
process-level (`step2_SPAtests.R` splits marker lists and `mclapply`s separate
`SPAGMMATtest` calls) plus RcppParallel/TBB inside step-1 GRM kernels only.
`src/LDmat.cpp` declares its own static reader pointers that are never assigned —
dead shadows; it actually uses Main.cpp's globals via `Unified_getOneMarker`.

### Format handling (R side)

`R/Geno.R`: `checkGenoInput` detects the format; `setGenoInput` builds a `markerInfo`
data.table and calls the right `set*objInCPP`. `genoIndex` semantics differ per format
(PLINK/PGEN row number; BGEN byte offset + `genoIndex_prev`; VCF iterator-driven) and
indices are passed to C++ **as strings** to survive R's lack of int64. `--AlleleOrder`
defaults differ: PLINK `alt-first`, BGEN/PGEN `ref-first` (enforced).
**`R/Geno_LDMat.R` is a drifted near-copy of this logic for step 3** — format-handling
changes must be made in both.

### The step1 → step2 contract (implicit, name-based)

- `<prefix>.rda` must contain a variable literally named `modglmm`;
  `R/readInGLMM.R::ReadModel` `load()`s it blind and mutates it (no schema/version
  field, only `is.null()` back-compat probes). Non-LOCO fits go to
  `<prefix>_noLOCO.rda`; `--isLowMemLOCO` writes per-chromosome `<prefix>_chr<j>.rda`.
- Variance-ratio file: headerless `value type category` rows with `type` ∈
  `{sparse, null, null_noXadj}` (`Get_Variance_Ratio` in `R/readInGLMM.R`; absent
  types become `-1` sentinels, which e.g. make `--is_fastTest` a hard error).
  Multiple rows per type = categorical variance ratios whose count must match
  `cateVarRatioMinMACVecExclude`; C++ picks the row per marker by MAC.
- LOCO means three different mechanisms: step-1 global-state mutation of `geno` per
  chromosome; the save-time file layout; and step-2 `ReadModel` overwriting
  fitted values from `LOCOResult[[chrom]]` (`--chrom` is mandatory when LOCO=TRUE).
