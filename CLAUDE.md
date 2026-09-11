# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**This is the `strip` branch**: SAIGE reduced to step 0 (createSparseGRM), step 1
(fitNULLGLMM), and step 2 (SPAGMMATtest: single-variant + region/gene tests),
with **PLINK bed/bim/fam and PGEN pgen/pvar/psam as genotype input** for every
step. The BGEN/VCF/SAV readers, step 3 (LD matrices), RareEffect and survival
traits were stripped, each verified byte-exact against `main`-HEAD goldens;
the PGEN reader was added back as fresh code behind the `GenoReader`
interface (`docs/state/pgen-step01.md`, `docs/state/s8-pgen-clone.md`).
`main` holds the full upstream fork (docs: https://saigegit.github.io/SAIGE-doc/).
Binary and quantitative traits are supported.

**Where to start.** Read the handoff that `docs/state/next.md` names in its
first line (several handoffs can share a date, and `ls | tail -1` sorts by
name), then `docs/state/next.md` (what is NOT
done, provenance marked on each item -- including the REFUTED candidates that
must not be re-opened), then `python3 harness/tools/ledger.py open`. Handoffs
are dated and never renamed. Plan files under `~/.claude/plans` are
session-local and unversioned: anything a later session needs goes in
`docs/state/`. The pre-trim version of this file, with every fix narrative and
gate count, is `docs/state/history-claude-md.md`; a fix's story lives in its
LEDGER row's STATUS cell, which is the authoritative version.

**Do not read `LEDGER.md` whole.** It is ~90 rows of up to ~2,600 characters
each. `python3 harness/tools/ledger.py open` lists the open rows, `ledger.py 56`
prints one wrapped, `ledger.py grep ingest` searches, `ledger.py lint` checks
that every row parses (it runs before every test tier). Never write a
hand-counted "N open of M" here or in `next.md`; the lint rejects it, because
the count drifted every time it was written by hand.

**Phase status.** Phases 0-7 are done (baseline harness; strip; test suite;
refactor -- `GenoReader`, de-globalised step-2 state, `TestContext`; speedup;
multi-weight region tests; the AF > 0.5 fixes; PGEN for every step) and so is
the step-2 threading sequence S0-S12 (`docs/design/step2-threading-verdict.md`).
**Phase 3 of the ledger-fix work swept the fixable rows on 2026-09-05**
(`docs/state/handoff-2026-09-05-phase3-sweep.md`): one commit per ledger row,
characterization test inverted in the same commit, red at unpatched HEAD, golden
movement predicted before the gate runs. What stays open is a list of decisions
for the owner, in `docs/state/next.md`.

**Step 1 fits several phenotypes in ONE process** (`--phenoCol=y1,y2`, writes
`<outputPrefix>_<pheno>.*` plus a `.multiTraits.txt` manifest that step 2
consumes; each phenotype's outputs byte-identical to a standalone run's):
`docs/state/step1-multipheno.md`. **Step 2 tests several null models in ONE
genotype pass** (`--GMMATmodel_varianceRatio_multiTraits_File`, on the marker
path and -- since 2026-09-06, `2ed00165` -- on the region path too, where every
model must have been fitted on ONE sample set; each model's file byte-identical
to its own single-model run): `docs/state/multimodel-step2.md`,
`docs/state/handoff-2026-09-06-region-multimodel.md`.

## Traps

Standing facts that have each cost at least one session. One line of why each.

- **A fix that moves no golden is a statement about coverage, not a pass.**
  #40, #44, #52, #53, #74 and #63 all left every golden byte-identical while
  being real. Predict which goldens should move BEFORE the gate, from the
  goldens themselves where possible (#66: the guard's ratio was readable from
  the printed `var`/`var_c` columns).
- **Every non-BURDEN group golden already takes the chunked region path**
  (`nchunks == 2`, the ultra-rare collapse appends a chunk), so
  `--markers_per_chunk_in_groupTest` changes nothing there. Since 2026-09-06
  the chunks stay in memory (`RegionChunkStore`); the chunk-file FAILURE is
  reached only by `test-step2-region-chunkio.R` with the spill forced
  (`SAIGE_REGION_CHUNK_MEMORY_MB=0`) and a directory planted where the chunk
  file must go.
- **`rs20` and `rs22` in `genotype_100markers.bed` ARE the covariates `a9`/`a10`**
  (R^2 = 1.000000). Any covariate-adjusting path reports `BETA=NA SE=NA` on them
  (LEDGER #44); expect two degenerate rows and do not read them as a signal.
- **No bundled marker has AF_Allele2 > 0.5.** Reach that half of the marker path
  with `--AlleleOrder=ref-first` on the same bed (golden P11); the default path
  is now exactly flip-invariant (#28/#29 fixed). A native AF > 0.5 fixture is
  still absent.
- **Do not force `--pCutoffforFirth` upward** to make Firth fire: above the
  intended cutoff the penalised fit drifts to `maxstep*maxit` = 750 (#26). Use
  0.05, which fires on 5 of the 70 bundled markers.
- **A statistic over degenerate output can invent a bug.** #26 was first
  written as "Firth effect sizes are wrong" from a Pearson correlation two
  forced +-750 outliers had dragged to zero; Spearman is 0.85. Inspect the
  in-regime rows before trusting an aggregate.
- **P35 is byte-identical to P2 on purpose** (LEDGER #51: the escalation pass
  recomputes what the fast pass had). It pins the other WRITER, which #50 lost
  19 rows through while P2 stayed green. Do not delete it as a duplicate.
- **A golden that reproduces P1 gates nothing.** `--max_MAC_for_ER=25` (and 20,
  50, 100) is byte-identical to P1; the ER resampling arm needs 200 on this
  fixture (P33). Check that a new golden actually differs.
- **`createSparseGRM` is not self-reproducing without a seed**: the
  `set.seed(1)` is in `extdata/createSparseGRM.R`, not the function, so a
  step-0 comparison that bypasses the CLI measures the RNG.
- **The sparse step-1 goldens read a frozen `.mtx` from `extdata/output`** that
  the harness does not regenerate and that still carries the #53 `geno2` bug.

## Rules of this branch

- **Bug-for-bug preservation.** Known bugs are recorded in `LEDGER.md`, NOT
  fixed in passing. Fixes are one commit per ledger row, with the row's
  characterization test inverted in the same commit and any golden that moves
  re-baselined after being proved byte-reproducible twice. Read the ledger row
  before "fixing" anything.
- **Every change is verified against the goldens** (Harness below). Step-2 and
  step-0 outputs must stay byte-identical; step-1 within tolerances.
- **Ledger rows are single markdown-table lines of exactly five cells.** A `|`
  in the prose -- `|Beta|`, `a || b`, `x |= y`, inside backticks or not -- must
  be written `\|`. `ledger.py lint` enforces it; before it existed, five OPEN
  rows were invisible to `ledger.py open` for days.
- **Push through the pre-push gate.** `bash harness/prepush.sh` runs the ledger
  lint and the standard tier in BOTH compare modes; `.githooks/pre-push` runs it
  on every push to `slim` once `git config core.hooksPath .githooks` is set (see
  Build). The byte-vs-numeric trap fired twice, once leaving CI red for five
  hours while the session recorded the tier as green.
- After adding/removing any `[[Rcpp::export]]` or `.cpp` file: update the
  explicit `OBJECTS` list in `src/Makevars` and re-run
  `pixi run Rscript -e 'Rcpp::compileAttributes(".")'`.
- **R code must pass the lint check.** `.lintr` at the repo root is the config
  (lintr defaults, 100-char lines, camelCase and `name_with_parts` accepted).
  Run it from a scratch lib, not the pinned env:
  `R_LIBS_USER=<scratch-lib> Rscript -e 'lintr::lint("<file>")'`. New or
  restyled R files must come back clean; the legacy `R/` sources are NOT
  retroactively restyled. Restyling is gated like any change (goldens
  byte-identical; CLI option blocks value-identical old vs new).
- **After ANY header (`.hpp`) change, delete `src/*.o` before `R CMD INSTALL`.**
  R's make has no header dependency tracking; stale dependents' objects after a
  class-layout change give intermittent memory corruption. `harness/build.sh`
  always cleans; `harness/run_tests.sh` does not build.

## Build

```sh
pixi install
bash harness/build.sh            # clean objects, compileAttributes, R CMD INSTALL
git config core.hooksPath .githooks   # once per clone: arms the pre-push gate
```

Works natively on macOS and linux. `src/Makevars` hardcodes
`-I../.pixi/envs/default/include`, so the in-tree `.pixi/` env is required.
`make` at the repo root runs the Docker build.

**plink-ng (pgenlib) is fetched and compiled automatically** by
`harness/build.sh` into `.pixi/envs/default/lib/libplink2_includes.a` (not
vendored; `.gitignore` excludes `plink-ng/`), guarded separately on the archive
and the headers because CI caches only the archive. `docker/Dockerfile` carries
the equivalent layer. `plink2` the CLI is a pixi dependency only for generating
PGEN fixtures (`linux-64`/`osx-64` only, hence the per-platform split in
`pixi.toml`); the installed package never needs it.

## Test suite

`bash harness/run_tests.sh [fast|standard|full]` -- testthat under
`tests/testthat/`, sourcing the pinned env itself; `ledger.py lint` and
`docdrift.py` run first.
- **fast** (~20s): unit tests over the R seams, R-vs-C++ CCT reconciliation,
  pure kernels pinned to closed forms (`test-unit-spa.R`, `test-unit-firth.R`),
  data-invariant checks over the goldens, and the suite's own comparator
  self-tests (`test-meta-comparators.R`, null + positive controls).
- **standard** (~1.5-2.5 min, default): fast + the step-1 genotype-store checks
  (`test-step1-genostore.R`, `test-step1-pgen.R`), the Guard Malloc pass over a
  step-1 LOCO fit (`test-step1-loco-oob.R`, macOS only), every step-2 CLI run
  against `harness/goldens/step2/` (goldens, marker paths, multi-weight, Firth, PGEN
  mirrors, threads, escalation, file integrity), and the characterization tests.
- **full** (~20+ min, same-platform only): standard + golden regeneration,
  step-1 determinism double-run, permuted-phenotype negative control.

`SAIGE_COMPARE_MODE=numeric` switches step-2 comparison from byte-exact to
parsed-table tolerance (1e-5 rel), required off the golden-generating platform;
CI (`.github/workflows/tests.yml`, on the `slim` remote
`astheeggeggs/saige-slim`, `gh run list --repo astheeggeggs/saige-slim`) uses it.
Detection power is drill-verified: a planted 0.1% variance-ratio perturbation
fails the standard tier.

**`harness/tools/docdrift.py` is a blocklist of statements this branch has
already written down wrongly**, each with the commit that falsified it and the
evidence to re-check. It is a regression guard, not a drift detector: it cannot
notice that a doc has gone stale, only that a known-false sentence is back --
which is the observed failure mode, because the five found on 2026-09-09 had
been COPIED from each other (CLAUDE.md into the roxygen, the roxygen into
`man/`, `man/` into the CLI help) and one twin had already been corrected while
the other was not. It scans only the CURRENT-state surfaces (`CLAUDE.md`,
`docs/state/next.md`, `R/`, `man/`, `extdata/`); dated handoffs and `LEDGER.md`
are correct about their own date and are never scanned. `docdrift.py selftest`
is the control that no rule has quietly stopped matching, and every rule is
drilled against `4c2ca05a^`, where the tool reports 10 hits. Add a rule only
after a claim has actually been written wrongly; a line that must quote a false
claim on purpose carries the marker the tool names.

Characterization tests pin ledger bugs, asserting what the code does TODAY
including malformed messages and unreachable guards. When one "fails" after a
fix, invert it and update the ledger row in the same commit.

**Parallel test files.** Needs `Config/testthat/parallel: true` and
`Config/testthat/edition: 3` in the INSTALLED `DESCRIPTION` (`find_parallel()`
reads the installed one; editing without `R CMD INSTALL` leaves parallelism
silently OFF -- check for the "Starting N test processes" line before believing
a timing) and `run_tests.sh` passing `package = "SAIGE", load_package =
"installed"`. Wall time is bounded by the slowest file, so long poles are split
over shared `helper-*.R` files (`mw_cache` is per-process, so a split can add a
run). `Config/testthat/start-first` measured no benefit and names FILES, so it
goes stale on rename. Check `uptime` before timing anything.

Three traps when adding tests:
- **Never assert bit-identity between two different code paths**, and gate
  every new cross-run comparison in BOTH compare modes before pushing. Two
  algebraically identical paths differ by 1 ULP and whether that reaches the
  printed digits is platform-dependent. Compare counts and labels exactly,
  statistics at ~1e-12 relative -- or, for step-1 quantities, at the harness's
  tolerances (`harness/compare_step1.R`: 1e-2 for theta and variance ratios,
  1e-3 otherwise). A test that demands more precision than the estimator has
  measures the platform's libm. Byte-exactness is a claim about ONE path re-run.
- Fixture facts above (Traps): `--pCutoffforFirth`, AF > 0.5, `rs20`/`rs22`.
- Edition 3 tightened `expect_equal`; the CCT round-trip through `tan` near
  pi/2 costs ~6.6e-12 at p = 1e-5, so `test-unit-cct.R` uses 1e-10.

### Measuring coverage (not part of any gate)

`covr` is not in the pinned env; install it to a scratch lib. Three
non-obvious things, each of which fails silently rather than loudly:
- Start from `getOption("covr.flags")` and only blank `FLIBS`. The defaults
  carry `CXX11/14/17/20FLAGS`, and `src/Makevars` sets `CXX_STD = CXX14`, so a
  hand-rolled vector with just `CXXFLAGS` instruments **nothing**. `FLIBS`
  must go because conda clang ships no `libgcov`.
- Use `type = "none"` + `code = 'setwd("<repo>/tests/testthat"); test_dir(".")'`.
  With `type = "tests"` covr runs from the installed tree, where
  `helper-saige.R`'s repo-root guess lands on the temp lib and every `extdata`
  path breaks.
- Run `gcov` with cwd = `src/`; point `covr.gcov` at a wrapper calling
  `llvm-cov gcov` version-matched to the conda clang, not Apple's.

covr builds in-place: afterwards delete `src/*.{gcno,gcda,o,so}` and the stray
`covr/` dir, then re-run `harness/build.sh` before any gate. Child processes
ARE captured. Baseline at C8: R 55.2%, live step-2 C++ 72.2% line / 42.9% branch.

## Harness (golden generation; the test suite drives comparison)

- `source harness/env.sh` before ANY run -- pins `RCPP_PARALLEL_NUM_THREADS=1
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. Without it full-GRM step-1 fits are
  nondeterministic (TBB float `parallel_reduce`, LEDGER #7).
- `bash harness/make_goldens.sh <dir>` regenerates step 0, the five step-1
  configurations (binary sparse x2, quantitative sparse, binary DENSE, binary
  sparse-VR), the step-0/step-1 PGEN mirrors of all of them (committed even
  though byte-identical to their PLINK twins -- `compare_goldens.sh`'s step-0
  loop flags "EXTRA in candidate", and `make_goldens.sh` `cmp`s each pgen output
  against its PLINK twin), and the step-2 runs P1-P37 plus 13 PGEN mirrors. The
  command for each run, and for P27-P33 the per-arm counts from an instrumented
  build, are recorded beside the run in the script. Step-2 runs consume frozen
  model inputs: the committed step-1 goldens plus
  `extdata/output/example_binary_sparseGRM.*`.
- `bash harness/compare_goldens.sh harness/goldens <candidate-dir>` -- byte
  compare for step 0/2; toleranced step 1 (theta/VR 1e-2, vectors 1e-3). It
  excludes `*.index` from the byte gate.
- **Before accepting any new golden, prove it byte-reproducible at unpatched
  HEAD** (two `make_goldens.sh` runs into separate trees, `cmp`); until that
  holds it gates nothing. `binary_denseGRM` is the ONLY step-1 golden on the
  dense-GRM path; the other four never leave the sparse branch.
- **P27-P33 each reach one distinct arm of the `getMarkerPval` cascade** (sparse
  `scoreTest` + PCG, `scoreTestFast`, the fastTest re-test block, Firth, ER
  resampling, unconditional quantitative), each P1's command with ONE flag
  changed. The sparse marker path needs a VR file with both `sparse` and
  `null_noXadj` rows, which only `--useSparseGRMforVarRatio=TRUE` produces
  (step-1 golden `binary_sparseGRM_vr`).
- `bash harness/tsan.sh [clean|plant|plantreader]` -- ThreadSanitizer over the
  in-process marker scan. **Always run a plant mode too**: `clean` must report
  zero races and `plant` (one `MarkerScratch` shared by every worker) must
  report races, because every macOS setup mistake fails by reporting NO races.
  Check `otool -L .../SAIGE.so | grep plink2` names the instrumented archive
  under `harness/work/tsan`; `libtbb` is not instrumented, so main-vs-worker
  races are expected and annotated under `SAIGE_TSAN_ANNOTATIONS`. Each arm runs
  under `SAIGE_TSAN_TIMEOUT` (default 3600s, `harness/timeout.sh`) and an arm
  that dies or times out is reported as that, not as `races=0` -- licensed under
  `plant`, and under `clean` it makes the script exit non-zero. Without it a
  crashed arm hangs on R's interactive abort prompt indefinitely (once for
  4h32m).
- `bash harness/guardmalloc.sh [clean|plant|both]` -- Guard Malloc
  (`/usr/lib/libgmalloc.dylib`, macOS only) over a two-chromosome step-1 LOCO
  fit: every allocation on its own page with an unmapped page after it, so a
  read past the end faults on the instruction that does it. This is what
  characterizes LEDGER #102; a bounds-checked (`ARMA_NO_DEBUG`-less) build does
  NOT, because armadillo's `operator[]` is unchecked in every build. **Always
  run the plant arm** (a deliberate out-of-range read at the same program
  point): `Rscript` routes through a SIP-protected binary that strips `DYLD_*`,
  so the GuardMalloc banner can print from the pre-exec process while the R
  process that runs the fit allocates from ordinary malloc -- the banner proves
  nothing, only the plant does. `test-step1-loco-oob.R` drives both arms.
- `bash harness/make_controls.sh <dir>` -- permuted-phenotype negative control
  through the full pipeline + KS uniformity check, through both readers.
- `bash harness/prepush.sh` -- the pre-push gate (Rules).
- Committed goldens: `harness/goldens/`; retired bgen/vcf: `harness/goldens_retired/`.
- PLINK group files key on **bim column 2 (rs IDs)**; `--condition` takes rsIDs.

## Multi-weight region tests (Phase 5)

Group tests run K >= 1 weight sets in one pass; results are per
(anno x max_MAF x weight set), weight-innermost, with a `Weight` label column
after `max_MAF` plus `Pvalue_ACATV`/`Pvalue_ACATO` -- schema only when K > 1 or
`--is_output_weight_column=TRUE`; K = 1 without the flag is byte-identical
legacy output. Set order (normative, identical in R and C++): **custom** sets
from group-file `weight[:NAME]` lines, then **beta** sets from
`--weights.beta="a1,b1;a2,b2"` (composed only when > 1 pair), then
**No_Weights** (`--is_no_weight_in_groupTest=TRUE`). Flat index
`jmr = j*q_maf*q_weight + m*q_weight + r` on both sides. Custom sets collapse
ultra-rare markers WEIGHTED (so `w0 = 1.0` on the conditional path, #17);
beta sets unweighted. Fixtures: `group_new_snpid_withWeights.txt` (K=1),
`group_new_snpid_multiWeights.txt` (K=2, whose CADD line is byte-identical to
the K=1 weight line -- one statistic computed two ways). Tests:
`test-step2-multiweight{,-cond,-behaviour}.R` over `helper-multiweight.R`.

## FlexRV (2026-09-09; F0-F4 landed, F5 measurement open)

FlexRV (Schwartzentruber et al., bioRxiv 2025.11.04.686320) is 16 score
transforms x 12 MAF transforms = 192 fixed weight vectors per gene, a burden
test each, Cauchy-combined into `p_FlexRV`. On this branch it is a fourth
weight-set kind of the Phase 5 engine, `SET_TRANSFORM`, expanded inside C++
from one `GENE score:NAME v1 v2 ...` group-file line (after `anno`, before any
weight line) and the marker's IN-SAMPLE MAF/MAC after QC. Kernels
`src/FlexRV.cpp`, spec parser and set names `R/FlexRV.R`, the pinned table
with every closed form and every resolved ambiguity
`docs/design/flexrv-design.md` (the MAF beta transform is DECREASING,
`(1 - MAF/M)^(b-1)`, zero at and above M; read A1 before touching it).
Custom-set semantics throughout: weighted ultra-rare collapse, `w0 = 1` on
the pseudo-marker, ER only for `unw_unw`; `--weights.beta` is inert on them.

The run: `--flexRV_score=NAME --flexRV_maxMAF=M --annotation_in_groupTest=<ONE
group> --r.corr=1`, optionally `--flexRV_scoreTransforms` /
`--flexRV_MAFTransforms` (spec strings, defaults are the paper's),
`--flexRV_lofAnno=lof` (exact labels, each must occur in the file). ONE
annotation group and ONE max MAF per run, or the pooled Cauchy row is not
`p_FlexRV`; `--maxMAF_in_groupTest` is replaced when absent and rejected when
different; `--condition` needs `--weights_for_condition`; M != 0.001 warns
and names the beta sets `mafbeta<b>@M`. Every set's row is written
(`--flexRV_output_sets` was decided but is NOT implemented). Three traps:
- **The bundled fixture's in-sample N is 1,000** (the models', not the fam's
  10,000), so max MAF 0.01 is MAC <= 20 and is the only cell with real markers
  at the default collapse cutoff; the goldens P58-P61 pin it, and there
  `MAC25`..`MAC500` are duplicates of `unw` (`helper-step2-flexrv.R`).
- **`MAC<=1,3,5,10` sit at or below the default collapse cutoff**, so those
  sets test the `max()` pseudo-marker alone; `--MACCutoff_to_CollapseUltraRare=0`
  is the paper's unbinned quantitative run (P59).
- **A weight of exactly 0 inside a binary SKAT-O cell returns NA** (LEDGER
  #117, pre-existing, custom weights too), which the `lof` mask and every MAC
  threshold produce; on a binary trait FlexRV is a `--r.corr=1` run.
- **A set's burden is tested with the CELL's unweighted MAC** (LEDGER #118).
  This is NOT upstream's semantics -- upstream has no per-set dimension at all
  (`jmr` appears zero times in its `Main.cpp`; `MAC_GroupVec`, `NumRare_GroupVec`
  and `NumUltraRare_GroupVec` are all `q_anno_maf`, per CELL, and mutually
  consistent). Phase 5 widened all three to `q_anno_maf_weight` but left the
  accumulations unconditional, so the non-collapsed part of `MAC` and both
  counters still describe the cell while the collapsed part describes the set.
  Variance-ratio band, ER eligibility, the min-group-MAC gate and the printed
  `MAC` column all read the cell quantity, so a `MAC<=3` set differs from the
  same markers run as a `--maxMAC_in_groupTest=3` cell by exactly the VR ratio
  when the two MACs fall in different bands. `test-step2-flexrv-reconcile.R`
  characterizes it.
- **`--minGroupMAC_in_BurdenTest` applies on BOTH region paths and defaults to
  1** (LEDGER #126). The SKAT-O path reports `Pvalue_Burden`, so gating burden
  testing only on `--r.corr=1` was incoherent; the SKAT-O half is enforced in
  `region_postprocess` on `MAC_GroupVec`, the C++ half in the BURDEN emitter on
  `genoSumcount_noweight`, and those are the same quantity computed two ways.
  The default is 1 because SAIGE-GENE+ describes no minimum-MAC filter at all
  and upstream's 5 has no published rationale. **Gate 1 is inert by
  construction** (a surviving cell has at least one allele), so a test that
  wants the gate must name a value.
- **The pooled (Cauchy) row now discloses how it was made** (LEDGER #127):
  `<SAIGEOutputFile>.pooledTests.txt` records per region the rows combined, the
  cells the gate refused, and WHICH CCT rule fired (`Cauchy`, or the
  `min(1, n*min(p))` Bonferroni fallback of #86, or `none`), and a region in
  which nothing was tested prints `NA` rather than 0.5, which is the CCT of no
  inputs. Sidecar only, written only when a region has something to disclose --
  no result file moves, and the arithmetic is unchanged because it is the
  paper's own.
`opt$flexRV_score` partial-matches `--flexRV_scoreTransforms` when the flag is
absent: the CLI reads the five flags with `[[ ]]`. The memory line a FlexRV run
prints once is the formula in the design doc's "Implementation notes", not a
measurement; the N = 1,000 fixture is startup-bound and cannot see the 192x
work, so RAM/CPU at scale is F5's open measurement on the 50k fixture.

## Architecture

### Entry points -> R functions

| Script (extdata/) | Function | Role |
|---|---|---|
| `createSparseGRM.R` | `createSparseGRM()` in `R/SAIGE_createSparseGRM.R` | Step 0: sparse GRM (`.mtx` + `.sampleIDs.txt`) |
| `step1_fitNULLGLMM.R` | `fitNULLGLMM()` in `R/SAIGE_fitGLMM_fast.R` | Step 1: null GLMM (AI-REML + PCG) + variance ratios |
| `step2_SPAtests.R` | `SPAGMMATtest()` in `R/SAIGE_Test_main.R` | Step 2: single-variant or region tests (branches on `groupFile`) |

### R <-> C++: configure-then-run over global state

All C++ via Rcpp (`RcppExports` generated; `NAMESPACE` uses `exportPattern`).
R pushes state into C++ **globals** via `set*InCPP()` calls, then invokes one
long-running `main*InCPP()` loop that writes output files from C++; R writes
only the `.index` progress file. Forgetting a `set*_GlobalVarsInCPP` call
silently leaves stale values. `Main.hpp` declarations are stale (LEDGER #5).

- `src/Main.cpp` -- step-2 bridge. `g_models`/`g_readers` own the null models
  and readers (`g_model`/`g_reader` alias element 0); `mainMarkerInCPP`,
  `mainRegionInCPP`, `Unified_getOneMarker`. Format-specific code is only
  `setPLINKobjInCPP`/`setPGENobjInCPP` and four `t_genoType ==` guards -- grep
  for them before assuming the list is exhaustive.
- `src/GenoReaderRegistry.cpp` -- `PlinkReader` (wraps `PLINK::PlinkClass`) and
  `PgenReader` (wraps `PGEN::PgenClass`). **Both clone, differently, and the
  difference is forced.** `PlinkReader` clones SHARE one `PlinkClass` by
  `shared_ptr`, which is safe only because `getOneMarker` is `const` -- do NOT
  add a `mutable` member to work around that. A pgenlib read MUTATES its
  reader, so `PgenClass` is one decode handle over a shared `PgenShared`
  (metadata, sample subset, `PgenFileInfo`), and `clone()` builds another
  handle (~130 bytes per sample each). `PgenReader::isForkSafe()` is `false`
  as a live hazard (a `FILE*` it fseeks). The `.bed` `pread` fallback runs only
  when mmap fails; force it with `SAIGE_PLINK_NO_MMAP=1` (it refuses to clone).
- `R/SAIGE_fitGLMM_fast.R::fitNULLGLMM` is ONE body for both traits (the
  glm family and the fitter are chosen by trait; the four lines that differ
  in arithmetic are marked "by trait"). **A binary fit's inner (alpha, b)
  pass converges to `tolCoef` (1e-5, `--tolCoef`), decoupled from `tol`
  since LEDGER #125**: tied to `tol = 0.02` it left the outer REML score
  four times noisier than the exact score, which is what the #122 guard
  was built to survive; `--tolCoef=0.02` reproduces the pre-#125 fitter.
  Its phenotype preparation is
  `preparePhenotypeFrame` in `R/PhenoPrep.R`, shared with `getNeff`
  (`R/SAIGE_extractNeff.R`, LEDGER #100). Region weight sets are described by
  `weightSetDescriptor` in `R/WeightSets.R`, the only constructor
  `region_postprocess` accepts (no NULL descriptor).
- `src/SAIGE_fitGLMM_fast.cpp` (~7400 lines) -- step-1 engine over a single
  global `genoClass geno`. PGEN enters by RE-ENCODING: `MarkerByteSource` has
  `BedByteSource` and `PgenByteSource` (pgenlib decode, then re-encode as a
  PLINK-1 record) and `setGenoObjFromSource` is the one body both drive, which
  is why reader byte-identity is structural. A second `setgeno` in one process
  REPLACES the accumulators (#55). Algorithm families: base, `_LOCO`, `_q`.
- `src/SAIGE_test.cpp` -- `SAIGEClass::getMarkerPval` cascade (score -> SPA if
  |Z| > cutoff -> ER exact test if binary & MAC <= `max_MAC_for_ER` -> Firth),
  per-MAC variance-ratio selection. `src/Binary_*`/`src/ER_binary_func.cpp` are
  the live ER machinery.
- Step 0 cannot be separated from step 1: both drive the same global `geno`.
- **Neither step-2 path writes to the model.** `SAIGE::TestContext`
  (`src/SAIGE_test.hpp`) carries the variance ratio, the sparse-GRM flag and
  the noadjCov flag from the loop into a `const` test; the `assign*`/`set_*_cur`
  family was DELETED (S3, S12). Build a context where the choice is made.
  **The fastTest ESCALATION signal is deliberately NOT a `TestContext` field**:
  it is `g_regionEscalatedPass` in `Main.cpp`, read once at the top of
  `mainRegionInCPP`. Sharing it with the per-marker sparse-GRM choice was
  LEDGER #27 and #50.
- **Marker-path per-marker state is split by ownership** (`src/MarkerLoop.hpp`):
  a `Marker::MarkerResult` per marker (the writer streams them) and a
  `Marker::MarkerScratch` per worker. Add an output as a `MarkerResult` field,
  a buffer as a `MarkerScratch` field; never a q-length vector, and never a
  `std::vector<bool>` (bit packing makes two markers' writes one word).
- One dataset per R process; no reentrancy. TBB also drives the step-1 GRM
  kernels (source of LEDGER #7).
- **A SKAT kernel failure is counted, not hidden** (LEDGER #57): the cell's
  rows stay NA, the run prints `SKAT failures: N ...` before `Analysis done!`,
  and `<SAIGEOutputFile>.skatFailures.txt` names each failed cell with the
  kernel's message when N > 0 (`R/SkatFailures.R`). No bundled input makes SKAT
  fail; `test-step2-skat-failures.R` plants one by patching `SKAT` in a child.

### Step-2 marker parallelism: `--nThreads`, one mechanism

`--nThreads` threads the C++ marker loop IN PROCESS via
`RcppParallel::parallelFor` with an EXPLICIT `numThreads` (which overrides
`env.sh`'s pin -- without that every threaded test would silently run serial).
One process, one model, one output file, byte-identical to a one-thread run; at
`--nThreads=1` the loop is literally serial. `R/MarkerScan.R` validates and
dispatches. The forked path was retired (S7) after measurement
(`docs/design/step2-threading-verdict.md`, `harness/bench/NOTES.md`).

**Byte-identity is the right gate here** despite the rule against comparing
two paths: one path over a disjoint partition, no cross-marker arithmetic.
`test-step2-threads-inprocess.R` holds it and **asserts the run REPORTED its
thread count** -- a knob that fails to arrive would otherwise pass as green
(the LEDGER #30 shape). Drilled both ways.

**Group tests thread too (2026-09-06):** `--nThreads` drives the C++ region
marker loop (one marker per worker inside each block of `markers_per_chunk`,
the order-dependent accumulation serial, byte-identical at any count;
`test-step2-region-threads.R`) AND the SKAT-O post-processing pool
(`R/SkatoParallel.R`, forked BEFORE the first C++ call because forking a
process that owns a TBB pool is the hazard; the worker function MUST be
top-level, or R ships the whole `outList` per task). The two alternate per
region. SKAT-O itself can never run in-process on threads: SKAT's `qfc` is
file-scope statics and a global `jmp_buf` (`docs/state/region-perf-handoff.md`).
The ultra-rare pass and the BURDEN emitter thread their per-cell TESTS as
well (`runCellTests` in `src/Main.cpp`, one cell per worker); what stays serial
in cell order is the prepare that feeds each test and the write-back that
consumes it. That took BURDEN from 1.57x/1.64x at 4 threads to ~1.84x/~1.71x
(`e7b346ba`, `harness/bench/NOTES.md`). **Identical cells are evaluated
once (2026-09-09):** a cell whose (Score, Phi, r.corr) is `identical()` to
an earlier cell's is served that answer (`skatoDedupe`/`skatoMemoKernel`,
`R/SkatoParallel.R`), the rows still emitted; the run prints
`SKAT-O cells: N asked for, M evaluated, K served`. Printed-equal rows are
not always bit-equal cells, so **K is a measurement and no test may assert its
value**: P6 has 4 twin rows and this box serves 3 (GENE2 at max_MAF 0.01 is one
Phi entry ~2 ULP off), while CI has served all 4 with the golden bytes unchanged
-- same commit, CI red for a docs-only push (run 34403657864). **What makes the
count differ is NOT identified**: ubuntu-24.04 image 20260907.300.1 served 4 in
run 34403657864 and 3 in run 34538797394, so the runner image is falsified as
the discriminator and no other cause has been measured. That is the reason the
bound, not the number, is what a test may assert. Assert the cell grid, that the line adds up,
that K >= 1 and K <= the twin count, and the golden bytes. One deliberate
degradation: a
reader that cannot clone (only the PLINK pread fallback now) drops to one
thread with a printed reason and an R `warning()`.

Three findings that outlive the fork path: **E0** chunk boundaries cannot reach
the output bytes (`markers_per_chunk` drives only a stdout counter); **E1**
`closeOutfile_single()`/`flushOutfile_singleInCPP()` exist because the output
stream and the `.index` flushed independently, which was LEDGER #47 (restart
duplicated rows; FIXED S9 by recording the file size per chunk line, `*.index`
is excluded from the byte gate) -- #49 remains (a restart re-runs the last
completed chunk); and `varianceRatioBinsHaveGap` survives in `R/MarkerScan.R`
only because its message beats the C++ throw's.

### The step1 -> step2 contract (implicit, name-based)

- `<prefix>.rda` must contain a variable literally named `modglmm`;
  `R/readInGLMM.R::ReadModel` `load()`s it blind and mutates it (no schema or
  version field). Non-LOCO fits go to `<prefix>_noLOCO.rda`; `--isLowMemLOCO`
  writes per-chromosome `<prefix>_chr<j>.rda`.
- **Quantitative step-1 models are fitted standardised and SAVED in phenotype
  units** (LEDGER #65, `R/PhenoScale.R`): `fitNULLGLMM` divides the phenotype
  by the sd of its covariate-adjusted values, fits and estimates variance
  ratios on that, then scales the `.rda` (theta by s^2, linear terms by s) and
  the `null`/`null_noXadj` variance ratios (by 1/s^2) back; the `sparse` ratio
  is scale-free and is NOT scaled (LEDGER #97). Anything that reads the model
  or the ratios sees phenotype units; anything inserted between the fit and
  the final rescale sees standardised ones.
- Variance-ratio file: headerless `value type category` rows, `type` in
  `{sparse, null, null_noXadj}` (`Get_Variance_Ratio` in `R/readInGLMM.R`).
  Step-1 sparse-GRM fits emit only `null`/`null_noXadj` rows, so step 2 must
  NOT be given `--sparseGRMFile` for those models; a group test needs `sparse`
  rows ONLY IF it is given `--sparseGRMFile` (measured 2026-09-08: a region run
  on a `null`-only LOCO model with no `--sparseGRMFile` produces rows normally --
  an earlier version of this line said group tests need them unconditionally).
  Multiple rows per type = categorical VRs
  matched to `cateVarRatioMinMACVecExclude`; C++ picks per marker by MAC.
  **A LOCO model that can also drive a group test needs BOTH halves and only one
  combination gives them**: `--useSparseGRMtoFitNULL=FALSE` (a dense fit, so
  LOCO survives `R/SAIGE_fitGLMM_fast.R:836`) with
  `--useSparseGRMforVarRatio=TRUE` (which supplies the `sparse` rows). Golden
  `binary_loco_vr_2chr`.
- LOCO is three mechanisms: step-1 global-state mutation per chromosome, the
  save-time file layout, and step-2 `ReadModel` overwriting fitted values from
  `LOCOResult[[chrom]]` (`--chrom` mandatory when LOCO=TRUE).
