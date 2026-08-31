# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**This is the `strip` branch**: SAIGE reduced to step 0 (createSparseGRM), step 1
(fitNULLGLMM), and step 2 (SPAGMMATtest: single-variant + region/gene tests),
with **PLINK bed/bim/fam as the only genotype input**. The BGEN/VCF/SAV/PGEN
readers, step 3 (LD matrices), and RareEffect were stripped in commits
"Strip 1"–"Strip 7", each verified byte-exact against `main`-HEAD goldens.
`main` holds the full upstream fork (docs: https://saigegit.github.io/SAIGE-doc/).
Binary and quantitative traits are supported; survival-trait support was
removed in Phase 3a of the refactor.

**The working plan lives at `~/.claude/plans/based-on-this-repository-silly-shore.md`**
(it supersedes `rosy-munching-pudding.md`, which holds the Phase 0-5 history):
Phases 0-3 are done (baseline/harness; strip; test suite; refactor — survival
removed in 3a, GenoReader interface + de-globalized Step2Config/Step2Output/
owned model in 3b-d, ten ledger fixes in 3e); Phase 5 (multi-weights) done
through C11; **Phase 6 done — LEDGER #28 and #29 FIXED** (326abe8f, 86a08f47:
the AF>0.5 marker path is now exactly flip-invariant, full tier 1225/1225).
Next: Phase 4 (speedup), re-scoped to step 2 — make step 2 measurable (both s2
bench workloads are ~90% process startup and can measure nothing), make
`--nThreads` actually work (it is a silent no-op without `--idstoIncludeFile`),
then per-marker allocation hygiene. Step-1 blocked-GetTrace is parked behind
that measurement. LEDGER.md still-OPEN: #6 (PCG tolerance, deliberately
deferred), #14 (unverified suspicion), #17-#27, #32 (marker-shard
variance-ratio carry), and #33-#39 from the type-polymorphic-builtin audit
(`diag`, `sample`, `1:n`, missing `drop = FALSE`) -- #33 and #34 are the two
live ones, and #33 is baked into goldens P3 and P10.
Step 1 now takes `--seed` (default 1).

## Rules of this branch

- **Bug-for-bug preservation.** Known bugs are recorded in `LEDGER.md`, NOT
  fixed. Fixes happen in a dedicated post-refactor phase, one commit per ledger
  entry, with regenerated goldens. Read `LEDGER.md` before "fixing" anything.
- **Every change is verified against the goldens** (see Harness below). Step-2
  and step-0 outputs must stay byte-identical; step-1 within tolerances.
- After adding/removing any `[[Rcpp::export]]` or `.cpp` file: update the
  explicit `OBJECTS` list in `src/Makevars` and re-run
  `pixi run Rscript -e 'Rcpp::compileAttributes(".")'`.
- **R code must pass the lint check.** `.lintr` at the repo root is the
  config (lintr defaults, 100-char lines, camelCase and the codebase's mixed
  `name_with_parts` style accepted). Run it with lintr from a scratch lib
  (not in the pinned env):
  `R_LIBS_USER=<scratch-lib> Rscript -e 'lintr::lint("<file>")'` — new or
  restyled R files must come back clean. Restyling is behavior-neutral by
  definition: gate it like any other change (goldens byte-identical), and for
  CLI option blocks verify old-vs-new `option_list` metadata value-identical
  (see the generator pattern used for the extdata scripts). The three extdata
  CLI scripts are lint-clean; the `R/` sources are legacy-style and are NOT
  retroactively restyled — the rule applies to files you add or rewrite.
- **After ANY header (`.hpp`) change, delete `src/*.o` before `R CMD INSTALL`.**
  R's make has no header dependency tracking: a class-layout change in a
  header with stale dependents' objects produces intermittent memory
  corruption (learned the hard way in Phase 4: two members added to
  `PlinkClass`, `GenoReaderRegistry.o` stale → segfaults in `readFamFile`).
  `harness/run_tests.sh` does not build; gate scripts should always clean.

## Build

```sh
pixi install
pixi run R CMD INSTALL .
```

(The former CRAN `lintools` step left with survival support in Phase 3a.)

Works natively on macOS and linux (no plink-ng step anymore — it died with the
PGEN reader). `src/Makevars` still hardcodes `-I../.pixi/envs/default/include`,
so the in-tree `.pixi/` env is required. `make` at the repo root runs the
Docker build.

## Test suite (Phase 2)

`bash harness/run_tests.sh [fast|standard|full]` — testthat suite under
`tests/testthat/` (sources the pinned env itself):
- **fast** (~15s): unit tests over the R seams (Get_Variance_Ratio,
  ReadModel, Geno.R, group-file parsing, Util.R, checkArgs validators),
  R-vs-C++ CCT reconciliation, pure-kernel tests pinned to closed forms
  (the binomial CGF and its derivatives in `test-unit-spa.R`; the exported
  `fast_logistf_fit` against `glm` in `test-unit-firth.R`), the degenerate
  `region_postprocess` corner, data-invariant checks over the goldens
  (score-test algebra, allele accounting, ACAT recomputation of the Cauchy
  rows), and comparator self-tests (the suite's own null + positive controls
  in `test-meta-comparators.R`).
- **standard** (~2.5 min, default): fast + the four PLINK step-2 CLI runs
  compared against `harness/goldens/step2/`, the multi-weight runs, and the
  live Firth path (`test-step2-firth.R`).
- **full** (~20+ min, same-platform only): standard + full golden regeneration,
  step-1 determinism double-run, permuted-phenotype negative control.

`SAIGE_COMPARE_MODE=numeric` switches step-2 comparison from byte-exact to
parsed-table tolerance (1e-5 rel) — required off the golden-generating
platform; `.github/workflows/tests.yml` uses it. The suite's detection power
is drill-verified: a planted 0.1% variance-ratio perturbation makes the
standard tier fail, and reverting restores green.
Characterization tests pin ledger bugs (#13-#16, #19-#26) — if one "fails"
after a fix, update the ledger and the test together in the dedicated fix
phase. They assert what the code does TODAY, including the malformed error
messages (#21) and the unreachable guards (#24, #25); read the LEDGER entry
before "correcting" one.

The suite runs test files **in parallel** (8 workers here; `TESTTHAT_CPUS`
overrides). This needs BOTH `Config/testthat/parallel: true` and
`Config/testthat/edition: 3` in `DESCRIPTION`, and `run_tests.sh` must pass
`package = "SAIGE", load_package = "installed"` to `test_dir` -- the default
`load_package = "none"` hands each callr worker an empty package name and dies
with "attempt to use zero-length variable name" before a single test runs.
**The trap: `find_parallel()` reads the INSTALLED DESCRIPTION, not the source
one**, so editing `DESCRIPTION` without re-running `R CMD INSTALL` leaves
parallelism silently OFF -- same wall time, no warning, no "Starting N test
processes" line. Check for that line before believing a timing. It bought
4:09 -> 2:29 on standard; wall time is bounded below by the slowest single
file, so splitting the long poles is the lever: `test-step2-multiweight.R`
was 103s and is now split three ways (schema / `-cond` / `-behaviour`, ~31s /
35s / 37s) over the shared `helper-multiweight.R`. Because `mw_cache` is
per-process, a run two of those files both need executes twice -- 3 such runs,
which is exactly the +3 `expect_cli_ok` assertions (1471 -> 1474) and the
+15s of extra serial work. `test-step2-threads.R` (~78s) is now the pole. Total CPU goes UP (223s -> 319s): every worker loads SAIGE afresh,
which is the 1.08s floor each of the ~150 step-2 child processes also pays.

Adopting edition 3 tightened `expect_equal` to waldo semantics and surfaced
two assertions in `test-unit-cct.R` whose stated 1e-12 tolerance was tighter
than the arithmetic allows: CCT round-trips through
`0.5 - atan(tan((0.5 - p) * pi)) / pi`, and for small p that argument sits
beside pi/2 where tan is ill-conditioned, costing 1.55e-12 R-vs-C++ and
6.55e-12 on the scalar identity, both entirely from p = 1e-5 (every vector
case is <= 5.3e-15). Edition 2 had been accepting them because its numeric
comparison was looser than the number written in the test. Relaxed to 1e-10.

Three standing traps when adding tests here:
- **Never assert bit-identity between two different code paths**, and gate new
  cross-run comparisons in BOTH compare modes (`SAIGE_COMPARE_MODE=numeric
  bash harness/run_tests.sh standard`) before pushing. Comparing a K>1 run's
  No_Weights column against a K=1 run with `expect_identical` on the parsed
  strings passed locally and went red on CI: algebraically identical, 1 ULP
  apart, and whether that reaches the printed digits is platform-dependent
  (macOS no, Linux yes — and one CI run hit it while the previous one did
  not). Compare counts and labels exactly, statistics at ~1e-12 relative.
  Byte-exactness is a claim about ONE code path re-run, which is what the
  goldens test; it is not a claim about two paths agreeing.
- **Do not force `--pCutoffforFirth` upward** to make the Firth arm fire.
  Firth is gated behind it because the penalised fit is expensive, and above
  the intended cutoff the fit stops converging and drifts to
  `maxstep*maxit` = 750 (LEDGER #26). Use 0.05, which fires on 5 of the 70
  bundled markers.
- **No bundled fixture has a marker with AF_Allele2 > 0.5** (max 0.498, and
  0.5 exactly in the chunked set), so the whole AF>0.5 half of the marker path
  is unexercised by every golden. Reach it by reading the SAME bed with
  `--AlleleOrder=ref-first`, which mirrors every marker across 0.5 — that is
  how LEDGER #28 (variance mis-centred by `scoreTestFast_noadjCov`) and #29
  (markers at exactly `--minMAC` dropped by the gate-1 float round trip) were
  found, and golden P11 `plink_marker_reffirst.txt` pins the path. It is a
  proxy for "ALT is the major allele", not a native fixture — a native AF>0.5
  fixture needs an edited `.bed`, not just a `.bim`, and is still worth adding.
  **Both bugs are FIXED (Phase 6, commits 326abe8f/86a08f47), so the assertions
  in `tests/testthat/test-step2-allele-order.R` are now inverted**: the default
  path is EXACTLY flip-invariant on all 70 markers, matching what was
  previously only true of the `--is_noadjCov=FALSE` control arm. That control
  arm is still the thing to reach for — it was the reference that proved the
  ref-first construction sound while the default path was broken, and it is
  what any future change to this path should be measured against.
- **A statistic over forced/degenerate output can invent a bug.** The #26
  entry was first written as "Firth effect sizes are wrong" on the strength
  of a Pearson correlation that two ±750 outliers had dragged to zero;
  Spearman over the same data is 0.85 and the estimator is fine. Inspect the
  handful of in-regime rows before trusting an aggregate.

### Measuring coverage (not part of any gate)

`covr` is not in the pinned env; install it to a scratch lib. Three
non-obvious things, each of which fails silently rather than loudly:
- Start from `getOption("covr.flags")` and only blank `FLIBS`. The defaults
  carry `CXX11/14/17/20FLAGS`, and `src/Makevars` sets `CXX_STD = CXX14`, so a
  hand-rolled vector with just `CXXFLAGS` instruments **nothing** — the link
  gets `--coverage`, no compile does, and no `.gcda` is ever written. `FLIBS`
  must go because conda clang ships no `libgcov` (`--coverage` at link pulls
  `libclang_rt.profile` instead).
- Use `type = "none"` + `code = 'setwd("<repo>/tests/testthat"); test_dir(".")'`.
  With `type = "tests"` covr runs from the installed tree, where
  `helper-saige.R`'s `getwd()/../..` repo-root guess lands on the temp lib
  (which has its own DESCRIPTION, so the `SAIGE_REPO` fallback never fires)
  and every `extdata` fixture path breaks.
- Run `gcov` with cwd = `src/`; from anywhere else it emits header-only
  `.gcov` files with no line data. Point `covr.gcov` at a wrapper calling
  `llvm-cov gcov` (version-matched to the conda clang, not Apple's).

covr builds in-place, so afterwards delete `src/*.{gcno,gcda,o,so}` and the
stray `covr/` dir, then re-run `harness/build.sh` before any gate.
Child processes ARE captured (they inherit `R_LIBS`), so CLI-driven code
counts. Baseline at C8: R 55.2% full-tier, live step-2 C++ 72.2% line /
42.9% branch (dead code and `RcppExports.cpp` excluded — see LEDGER for the
dead files).

## Harness (golden generation; the test suite drives comparison)

- `source harness/env.sh` before ANY run — pins
  `RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.
  Without it, full-GRM step-1 fits are nondeterministic (TBB float
  `parallel_reduce`; empirically confirmed — see LEDGER.md #7). Sparse-GRM
  fits are byte-stable either way.
- `bash harness/make_goldens.sh <dir>` — regenerates the 8 runs (step0 sparse
  GRM; step1 binary ×2 / quantitative; step2 PLINK marker, group
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

## Multi-weight region tests (Phase 5, through C8)

Group tests run K >= 1 weight sets in one pass; results are per
(anno × max_MAF × weight set), weight-innermost, with a `Weight` label column
right after `max_MAF` plus `Pvalue_ACATV`/`Pvalue_ACATO` (SKAT-O path) — the
schema appears only when K > 1 or `--is_output_weight_column=TRUE`; K = 1
without the flag is byte-identical legacy output. One pooled Cauchy row per
region, unchanged in count. Set composition (order is normative, identical in
R and C++): **custom** sets from group-file weight lines (`weight`,
`weight:NAME`, or `weight,NAME` as line type; same names/order in every
region) → **beta** sets from `--weights.beta="a1,b1;a2,b2"` (each pair a set
named `Beta_a_b`; composed ONLY when >1 pair — a single pair keeps the legacy
implicit path, incl. the LEDGER #19 R-vs-C++ asymmetry) → **No_Weights**
(all-ones; appended by `--is_no_weight_in_groupTest=TRUE` when other sets
exist, replacing the old hard error for weighted files). Per-set semantics
(design D2): custom w0=w_r with weighted UR collapse max(w_r·G); beta
w0=dbeta(MAF,a,b) with unweighted UR collapse; No_Weights all 1s. The
normative flat index is `jmr = j*q_maf*q_weight + m*q_weight + r` (weight
fastest) on both sides. K multiplies UR-collapse columns and output rows, not
marker reads — real markers are read and score-tested once. Fixtures:
`extdata/input/group_new_snpid_withWeights.txt` (K=1),
`group_new_snpid_multiWeights.txt` (K=2: CADD, REVEL); tests in
`tests/testthat/test-step2-multiweight{,-cond,-behaviour}.R` (standard tier,
sharing `helper-multiweight.R`) and the fast-tier
unit files. Conditional analysis (`--condition`) composes with K > 1
and the forced Weight column (C11, design D7): the conditioning-weight vector
stays a SINGLE vector across sets; at K > 1 C++ exports the G1-raw
`G1tilde_P_G2tilde` matrix (plus `VarInvMat_G2_cond_w0`/`Tstat_G2_cond_w0`)
and skips its baked `TstatAdjCond`/`VarMatAdjCond`, with R applying set r's
G1-side weight per cell; at q_weight == 1 the baked C++ path runs verbatim
(P3 byte-identical). Golden P10 `plink_group_cond_multiweight.txt` pins the
K=2 conditional run. Reconciliation tests validate the R-side adjustment
against the K=1 baked path through the BETA route only — a K=1
custom-weighted conditional run reads uninitialized `w0` on ultra-rare rows
(LEDGER #17, K=1-only after C11) and cannot serve as the comparison
baseline.

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
  `PlinkClass`) for step 2. Algorithms exist as near-duplicate families —
  base, `_LOCO`, `_q` (quantitative) — after the `_Surv` variants were
  removed with survival support in Phase 3a.
- `src/SAIGE_test.cpp` — `SAIGEClass`: `getMarkerPval` cascade (score → SPA if
  |Z| > cutoff → ER exact test if binary & MAC ≤ `max_MAC_for_ER` → Firth),
  per-MAC variance-ratio selection. `src/Binary_*`/`src/ER_binary_func.cpp`
  are the live ER machinery (imported from SKAT) — keep.
- Step 0 cannot be separated from step 1: both drive the same global `geno`,
  and `createSparseKinParallel` is defined in `R/SAIGE_fitGLMM_fast.R`.
- One dataset per R process; no reentrancy. Parallelism is process-level
  (step-2 `--nThreads`: `R/MarkerShard.R` partitions the canonical marker
  order and `mclapply`s — see below) plus RcppParallel/TBB in step-1 GRM
  kernels (source of LEDGER.md #7).

### Step-2 `--nThreads` (Phase 4 C2)

`R/MarkerShard.R` owns it, not `extdata/step2_SPAtests.R` — the CLI cannot,
because the order that must be reproduced is the CHROM/POS-sorted `genoIndex`
built inside `SPAGMMATtest`, and because the fork has to happen **after**
`ReadModel` and after the model reaches the C++ singletons (children then share
one copy-on-write model instead of each re-loading it). `runMarkerScan`
dispatches: `nThreads == 1` is the old `SAIGE.Marker` call verbatim; above
that, `splitMarkerShards` cuts `min(nThreads, M)` contiguous near-equal blocks
of the **input** index, each shard runs the same serial code, and
`mergeShardOutputs` concatenates raw bytes. Single-variant only — group tests
warn and drop to one thread; restart (`--is_overwrite_output=FALSE` with a
partial index) is refused.

Three things this rests on, each verified rather than assumed:
- **Chunk boundaries cannot reach the output bytes** (experiment E0): golden P7
  is byte-identical at `--markers_per_chunk` ∈ {1000, 1500, 2000, 3000, 10000}.
  `markers_per_chunk` reaches C++ only as a stdout progress counter.
- **A forked child MUST close `g_out.single` explicitly** (experiment E1):
  `parallel:::mcexit` ends a child with `_exit()`, no static destructor runs,
  and the child silently loses everything past the stream's 192 KiB buffer
  boundary *while still returning success*. Hence `closeOutfile_single()` at
  the end of `SAIGE.Marker` — byte-neutral serially, mandatory in a shard.
- **The reader must be fork-safe.** The mmap path is (pointer arithmetic into
  a `MAP_PRIVATE` mapping); the `fread` fallback is not (children would share
  one file offset). `isGenoReaderForkSafeInCPP()` errors at `nThreads > 1`
  rather than silently decoding each other's records.

Byte-exactness IS the right claim here, despite the cross-path trap below:
each shard is ONE code path re-run over a disjoint subset, no marker's
arithmetic depends on another's, and the only cross-marker aggregates go to
stdout. `make_goldens.sh` is deliberately NOT modified — the gate's whole
claim is that a threaded run reproduces the existing *serial* goldens.

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
