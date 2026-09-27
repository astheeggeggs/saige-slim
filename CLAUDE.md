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
- **No HARD-CALL bundled marker has AF_Allele2 > 0.5, and the one dosage marker
  that does is useless for the question.** Reach that half of the marker path
  with `--AlleleOrder=ref-first` on the same bed (golden P11); the default path
  is now exactly flip-invariant (#28/#29 fixed). The earlier wording of this
  trap said "no bundled marker", which is false of the repo: `rs26` in the PGEN
  dosage fixture is AF 0.500034 (AC 1000.07 on N = 1,000) and appears in all
  three dosage goldens -- `pgen_marker_dosage.txt`, `_mac100`, `_nozerod` -- so
  it IS natively ALT-major and IS folded by `imputeGenoAndFlip`. But it sits
  3.4e-5 above 0.5, so it discriminates nothing about weighting, MAC bands or
  annotation binning, which is what #131's revisit trigger needs. A native
  AF > 0.5 fixture with real separation is still absent; the cheap way to make
  one is `harness/tools/make_dosage_fixture.py`, which already generates that file.
- **Every binary step-1 AND step-2 run on the bundled fixture is REFUSED unless
  the gates are off** (2026-09-23, `R/FitGates.R`). `--minCaseCount` defaults to
  100 and the fixture has 98 cases; `--minCovariateCellCount` defaults to 10 and
  `a10` level 2 has 7 cases. Both numbers are properties of the fixture, not
  accidents: it is a mechanical fixture, not a valid analysis. `harness/env.sh`
  sets `SAIGE_FIT_GATES=0`, which every harness script and `run_tests.sh`
  inherit, so nothing in the repo carries an opt-out flag -- and a run with the
  gates off PRINTS a WARNING line saying so, which is how you tell the two
  states apart in a log. `tests/testthat/test-step1-fit-gates.R` is the only
  place the shipped defaults are exercised; it sets `SAIGE_FIT_GATES=1` per
  test. If a step-1 or step-2 run refuses with "below --minCaseCount", the
  question is whether `env.sh` was sourced, not whether the gate is wrong.
  The thresholds are the branch's OWN measurement (`harness/work/nullrep-20260918/`,
  LEDGER #159): binary SKAT-O null tail 2.7x nominal at 1e-4 with 100 cases,
  1.5x at 1e-3 with 500; the quantitative arm is calibrated, so
  `--minSampleCount` defaults to 0 and is a power knob, not a correctness one.
  Separation is judged as a FRACTION, not as any pinned sample:
  `--separationMaxFraction` (default 0.02) -- a pinned sample contributes ~0 to
  the score variance instead of `mu*(1-mu)`, so it costs one sample's worth of
  information, and the rule this replaced refused a real model over fewer than 20
  samples in 2,737. Above zero but at or below the fraction is a warning, not a refusal.
  `--convergenceCheck` (default error) refuses a fit whose `converged` flag is
  FALSE -- which covers BOTH the loop running out at `--maxiter` AND the
  bail-out at `max(tau) > tol^(-2)`; only the step-1 log tells them apart, and
  `--tol` sets both the criterion and that bail-out (0.02 aborts above tau
  2500, 0.2 above tau 25). All 13 committed step-1 goldens are `converged=TRUE`,
  so that gate is inert here and a test asserts it stays so.
  Two halves that are NOT redundant: step 1 scans covariate cells and the
  covariate GLM; step 2 reads the loaded model's own `y` and `fitted.values`,
  which is the only check a model fitted by an older SAIGE ever meets. And the
  step-1 cell scan sees `model.matrix` DUMMIES, where the reference level has
  no column -- so it reconstructs the factor (a sample is in the reference
  level when every dummy is 0); measured, a per-dummy scan misses a zero-case
  reference level that a factor scan catches.
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
- **When the full gate is NOT required.** The gate is **~7 min per compare
  mode, 14 min total** (measured 2026-09-12, `harness/bench/NOTES.md`), and one
  full both-modes run was spent that day on a diff that could not possibly move
  a golden. A golden moves only if the INSTALLED PACKAGE's behaviour moves, so:
  - **docs-only** (`docs/`, `*.md`, comments): no tier. Run `ledger.py lint`
    and `docdrift.py` (and `optiondrift.py` if an `extdata/*.R` wrapper
    changed) -- they are what a docs change can break, and all are seconds.
  - **TSAN-annotation-only** (`SAIGE_TSAN_HAPPENS_*`): cannot move a golden,
    the macros are `((void)0)` unless `harness/tsan.sh` defines
    `SAIGE_TSAN_ANNOTATIONS`. Run the thing it actually affects --
    `SAIGE_TSAN_WORKLOADS=<arm> bash harness/tsan.sh clean`.
  - **test-only**: run the changed file (`test_file`) at `standard`. Escalate
    to BOTH compare modes only if it adds a CROSS-RUN comparison -- that is the
    byte-vs-numeric trap's shape, and nothing else in a test file has it.
  - **anything under `src/` or `R/`**: the full gate, both modes, no
    exceptions. That includes a change you are sure is inert; "inert" is what
    the goldens are for.
  The gate is also ~33% faster on a warm page cache than a cold one, so a
  timing taken right after a clean build is not comparable to one taken after a
  previous run -- and a timing without its conditions is not a measurement.
- **Push through the pre-push gate.** `bash harness/prepush.sh` runs the ledger
  lint and ONE standard tier at `SAIGE_COMPARE_MODE=both`, which generates every
  output once and applies BOTH comparators to it -- byte-exact and parsed-table
  tolerance -- so the both-modes requirement is tightened rather than relaxed
  (each comparison site now asserts both on the SAME bytes; before, each of two
  runs asserted one). **7m00s instead of 13m46s, measured 1.97x**, because the
  two modes never differed in what they GENERATED, only in how they compared.
  `SAIGE_PREPUSH_MODES="byte numeric"` restores two separate tier runs, which is
  the thing to use if a generation is ever suspected of being non-deterministic
  -- that is the one question two runs can answer and one cannot. CI stays on
  `numeric` with its own generation, because byte-exactness is only a valid
  claim on the platform that made the goldens. A bad mode value is now rejected;
  it used to be silently treated as numeric, so `SAIGE_COMPARE_MODE=bytes` ran
  the weaker comparator and said nothing; `.githooks/pre-push` runs it
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
- **Compiler warnings are a ratchet, not noise (2026-09-14).** `src/Makevars`
  carried `-w` for the branch's whole life, which also cancelled the `-Wall` R
  passes by default; `harness/build.sh` now counts the src-anchored warnings
  of a CLEAN build (generated `RcppExports.cpp` excluded) against the first
  line of `harness/warnings-baseline.txt` and refuses a build with more. Fewer
  is reported so the baseline can be lowered in the same commit. Local clang
  only: CI (`SAIGE_WARNINGS_RATCHET=0`) promotes error CLASSES through
  `SAIGE_EXTRA_CXXFLAGS` instead, because a count is not portable across
  compilers. Do not "fix" `SAIGE_test.cpp`'s `m_is_Firth_beta & (pval <= ...)`
  to `&&`: the parentheses are the existing parse made explicit, the operator
  is deliberate. The four `-Wsometimes-uninitialized` sites are LEDGER #152.
- **After ANY header (`.hpp`) change, delete `src/*.o` before `R CMD INSTALL`.**
  R's make has no header dependency tracking; stale dependents' objects after a
  class-layout change give intermittent memory corruption. `harness/build.sh`
  always cleans; `harness/run_tests.sh` does not build.

## Build

```sh
pixi install
bash harness/build.sh            # clean objects, compileAttributes, R CMD INSTALL
bash harness/build.sh --incremental   # keeps src/*.o (167s -> 57s); REFUSES if a
                                      # header is newer than the oldest object
git config core.hooksPath .githooks   # once per clone: arms the pre-push gate
```

**`pixi.lock` is committed (2026-09-15)**: `pixi install` reproduces the
pinned environment; changing `pixi.toml` regenerates the lock and the diff
lands in the same commit; moving a version is a deliberate `pixi update`
with its own gate. **`make docker`** builds `docker/Dockerfile`, a two-stage
image whose runtime stage is the `runtime` pixi environment
(`[feature.runtime]`, `no-default-feature`, 109 packages) with the
compilers pruned, the installed package and the four scripts: 595 MB,
reproduces golden P1 byte-identical on Linux. A `vX.Y.Z` tag equal to
`DESCRIPTION`'s Version runs `.github/workflows/release.yml`, which builds,
smoke-tests against P1 and pushes `ghcr.io/astheeggeggs/saige-slim:X.Y.Z`
(`docker/readme.md`).

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
`docdrift.py` and `optiondrift.py` (LEDGER #151: the CLI wrappers' defaults against
the exported formals, with an allowed-divergence list keyed to ledger rows) run first.
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
`both` applies the two comparators to one generation and is what the local gate
uses (validated: a value that is none of the three is an error, not a silent
downgrade to numeric);
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
drilled against `4c2ca05a^`, where the tool reports 45 hits (10 before the six
#136 rules of 2026-09-13, 42 before the #139 rule; the #142 rule's sentence
postdates that commit and adds none; the `arma-no-debug-noop` rule of
2026-09-19 adds 2, measured; the `170-lofanno-must-occur` rule of 2026-09-27
adds 0, its sentence postdates that commit). Add a rule only
after a claim has actually been written wrongly; a line that must quote a false
claim on purpose carries the marker the tool names.

**Documentation is generated, and `man/` is gated (2026-09-15).** `bash
harness/docs.sh` (from the pixi env: `pixi run bash harness/docs.sh`)
regenerates `man/` from the roxygen with `roxygenise(roclets = "rd",
load_code = load_installed)` -- the rd roclet only, because `NAMESPACE` is
hand-written; `load_installed` because pkgload would recompile `src/` in-tree
-- and builds the pkgdown site into `docs/site/` (gitignored; `_pkgdown.yml`
sets `destination`, since pkgdown's default is `docs/`, which holds the
hand-written pages). `.github/workflows/docs.yml` does the same on every push
to `slim`, runs `git diff --exit-code man/` (a roxygen edit that was not
regenerated fails the push), and deploys to GitHub Pages. The CLI-options
article `vignettes/articles/cli_options.Rmd` is rendered at build time from
`optiondrift.py --json`, so the published flag table IS the map the test
tiers gate on. Three pkgdown facts that each cost a build: article topics
are named `articles/<file>`, a name with a hyphen is parsed as subtraction
(use underscores), and plain `.md` files are not articles (the two
`docs/*.md` pages are wrapped as `.Rmd` copies by the script). A roxygen
edit is docs-only for the gate, but run `harness/docs.sh` before pushing or
CI is red on the `man/` step.

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
  **Run the `region` workload after ANY change to region threading.** It was
  in the default list from the start but not in the levers' gate list, and the
  first run after the 2026-09-12 cell-pool work found a real worker-vs-worker
  race on `selectVarianceRatio`'s out-flag (5 at t=4, 3 at t=8; `e72e4186`)
  that every byte gate had passed -- `bytes-vs-serial` was `identical` in the
  racing build too. The clean arm reads `races=0` only since the FORK edge is
  declared as well as the join edge; before that it reported 1-2 main-vs-worker
  warnings in every build, so a non-zero count there was not evidence of
  anything (`harness/bench/NOTES.md`).
- `bash harness/guardmalloc.sh [clean|plant|both]` -- Guard Malloc
  (`/usr/lib/libgmalloc.dylib`, macOS only) over a two-chromosome step-1 LOCO
  fit: every allocation on its own page with an unmapped page after it, so a
  read past the end faults on the instruction that does it. This is what
  characterizes LEDGER #102; armadillo's own bounds checks do NOT, because the
  site uses `operator[]`, which is unchecked in every build by construction
  (`Mat_meat.hpp`: "no bounds check"), as are `.at()` and `.at(r,c)`.
  **`-DARMA_NO_DEBUG` in `src/Makevars` is a NO-OP and always has been**
  (measured 2026-09-19): on armadillo 15.x it undefines only `ARMA_DEBUG`,
  which is off by default, while bounds checking is gated on
  `ARMA_CHECK_CONFORMANCE`, defined by default and removed only by
  `-DARMA_DONT_CHECK_CONFORMANCE`, which nothing sets. So the shipped build
  HAS armadillo's checks compiled in -- 250 `arma_stop_bounds_error` sites, and
  `strings` on the installed `SAIGE.so` shows `Mat::operator(): index out of
  bounds`. They cost ~1% and they are load-bearing: **do not add
  `-DARMA_DONT_CHECK_CONFORMANCE`**, which would turn #102's throw into silent
  UB no golden can see. **Always
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
transforms x 12 MAF transforms = up to 192 fixed weight vectors per gene
(fewer at small N, LEDGER #142 below), a burden
test each, Cauchy-combined into `p_FlexRV`. On this branch it is a fourth
weight-set kind of the Phase 5 engine, `SET_TRANSFORM`, expanded inside C++
from one `GENE score:NAME v1 v2 ...` group-file line (after `anno`, before any
weight line) and the marker's IN-SAMPLE MAF/MAC after QC. Kernels
`src/FlexRV.cpp`, spec parser and set names `R/FlexRV.R`, the pinned table
with every closed form and every resolved ambiguity
`docs/design/flexrv-design.md` (the MAF beta transform is DECREASING,
`(1 - MAF/M)^(b-1)`, zero at and above M; read A1 before touching it).
Custom-set semantics throughout: weighted ultra-rare collapse, `w0 = 1` on
the pseudo-marker; EVERY set is eligible for the exact test at MAC <= `--max_MAC_for_ER` since 2026-09-26 (LEDGER #166; it used to be `unw_unw` alone, on a premise #166 falsified), the vector is scaled to a largest entry of 1 before that test (LEDGER #167: its tie epsilon is absolute), and a set whose saddlepoint fails above the cutoff falls back to the exact test at up to 200 carriers, counted in the `SPA fallbacks:` line and named with the reason in `<SAIGEOutputFile>.spaFallbacks.txt`; `--weights.beta` is inert on them.

The run: `--flexRV_score=NAME --flexRV_maxMAF=M --annotation_in_groupTest=<ONE
group> --r.corr=1`, optionally `--flexRV_scoreTransforms` /
`--flexRV_MAFTransforms` (spec strings, defaults are the paper's),
`--flexRV_lofAnno=lof` (exact labels; one absent from the file WARNS, and if
none is present the `lof` sets keep no marker and are refused in every region
-- LEDGER #170, 2026-09-27: the old error killed a whole chromosome of a
filtered AoU rerun whose only region had no pLoF variant). ONE
annotation group and ONE max MAF per run, or the pooled Cauchy row is not
`p_FlexRV`; `--maxMAF_in_groupTest` is replaced when absent and rejected when
different; `--condition` needs `--weights_for_condition`; M != 0.001 warns
and names the beta sets `mafbeta<b>@M`. Every set's row is written
(`--flexRV_output_sets` was decided but is NOT implemented). Three traps:
- **The bundled fixture's in-sample N is 1,000** (the models', not the fam's
  10,000), so max MAF 0.01 is MAC <= 20 and is the only cell with real markers
  at the default collapse cutoff; the goldens P58-P61 pin it. **A MAC rung at
  or above `2*N*maxMAF` is dropped before the sets are composed** (LEDGER
  #142, decided 2026-09-13: it would weight every admitted marker 1 and
  duplicate `unw`, evaluated and pooled again under another name), so on this
  fixture `MAC25`..`MAC500` are never emitted and a default run tests 112
  sets (16 x 7), not 192; at N = 50,000 and the paper's max MAF the rungs
  100, 200 and 500 go (144 sets). The run says so in a warning naming the
  bound and the rungs. The ladder's values stay absolute; only what cannot
  differ at this N is dropped, and at the paper's N (~424k) nothing is.
- **`MAC<=1,3,5,10` sit at or below the default collapse cutoff**, so those
  sets test the `max()` pseudo-marker alone; `--MACCutoff_to_CollapseUltraRare=0`
  is the paper's unbinned quantitative run (P59).
- **LEAVE THE COLLAPSE CUTOFF AT 10 FOR FlexRV; do NOT recommend 0, and
  especially not on a binary trait** (decision, 2026-09-13). The paper's unbinned
  run IS available at 0 and P59 pins it, but SAIGE-GENE+ introduced the collapse
  FOR CALIBRATION, this repo has never measured the regime where that matters,
  and both the bundled fixture and the 50k sim study run with the collapse ON --
  so every calibration result this branch has, null and power, describes cutoff
  10 (`docs/state/handoff-2026-09-13-flexrv-sweep-rerun.md`). Turning it off is
  not a neutral "more faithful to the paper" choice: it removes the mechanism the
  method's own calibration was built on, in a regime nobody here has checked.
  **The gate on changing this is a replicate NULL sweep at cutoff 0 on a binary
  trait showing the type-I error holds.** Until that exists, 0 is a
  quantitative-only, paper-reproduction setting.
- **A weight of exactly 0 inside a binary SKAT-O cell returns NA** (LEDGER
  #117, pre-existing, custom weights too), which the `lof` mask and every MAC
  threshold produce; on a binary trait FlexRV is a `--r.corr=1` run.
- **A set's MAC is the set's own, and it is an allele count** (LEDGER #118,
  FIXED 2026-09-11). `MAC`, `MAC_case`, `MAC_control`, `Number_rare` and
  `Number_ultra_rare` count only the markers the set keeps (`w_r != 0`), the
  ultra-rare pseudo-marker contributes its UNWEIGHTED collapse
  (`genoURMat_noweights`, widened to `jmr`), and `genoSumcount_noweight` --
  which the variance-ratio band, ER eligibility and the min-group-MAC gate read
  -- is per set. So a `MAC<=3` set now equals the same markers run as a
  `--maxMAC_in_groupTest=3` cell, which is what
  `test-step2-flexrv-reconcile.R` asserts. Two consequences to hold on to: a
  CUSTOM-weighted cell's printed MAC stopped being the weighted pseudo-marker
  sum (the fractional MACs -- 3.34, 82.517 -- are gone, and every weight set of
  a cell now prints the SAME MAC as the plain unweighted run of that cell,
  which is the check to make); and on the user-weights path that is a knowing
  divergence from upstream, whose pseudo contribution stays weighted.
- **`--minGroupMAC_in_BurdenTest` applies on BOTH region paths and defaults to
  1** (LEDGER #126). The SKAT-O path reports `Pvalue_Burden`, so gating burden
  testing only on `--r.corr=1` was incoherent; the SKAT-O half is enforced in
  `region_postprocess` on `MAC_GroupVec`, the C++ half in the BURDEN emitter on
  `genoSumcount_noweight`, and those are the same quantity computed two ways.
  The default is 1 because SAIGE-GENE+ describes no minimum-MAC filter at all
  and upstream's 5 has no published rationale. **Gate 1 is inert by
  construction** (a surviving cell has at least one allele), so a test that
  wants the gate must name a value. **A raised gate can kill a SKAT-O run**
  (LEDGER #128, open): a region with two or more cells all refused dies in
  `regionFinish`'s `rbind`; one cell or `--r.corr=1` escapes. The MAC 1-4
  cells gate 1 admits are calibrated on sim50k, both traits
  (`docs/state/prediction-126-null-calibration.md`).
- **The pooled (Cauchy) row now discloses how it was made** (LEDGER #127):
  `<SAIGEOutputFile>.pooledTests.txt` records per region the rows combined, the
  cells the gate refused, and WHICH CCT rule fired (`Cauchy`, or the
  `min(1, n*min(p))` Bonferroni fallback of #86, or `none`), and a region in
  which nothing was tested prints `NA` rather than 0.5, which is the CCT of no
  inputs. Sidecar only, written only when a region has something to disclose --
  no result file moves, and the arithmetic is unchanged because it is the
  paper's own.
`opt$flexRV_score` partial-matches `--flexRV_scoreTransforms` when the flag is
absent: the CLI reads the five flags with `[[ ]]`. **The region path's three
per-(sample x cell) matrices live in carrier space since 2026-09-12**
(`CarrierStore`, `src/Main.cpp`: one column per carrier of the region; the
per-cell N-vector a test needs is scattered into the worker's scratch and
cleared after; the reductions over the carriers reproduce armadillo's
two-accumulator sums bit for bit; the dense `urFill` path it shipped with is
now DEAD, because LEDGER #130 stopped folding the pseudo-marker and a folded
one was its only writer), the per-set prepares run on the cell pool, and identical BURDEN
cells are evaluated once (`BURDEN cells: ... served`, a measurement like the
SKAT-O line). The memory line a FlexRV run prints once is an estimate of the
SKAT-O terms only (0 N-dependent bytes per cell on BURDEN), not a measurement.
Measured on sim50k, 192 sets, BURDEN: peak RSS 874 -> ~520 MB, wall 57 -> 19 s
at 1 thread (`docs/state/handoff-2026-09-12-flexrv-ram-levers.md`).

## Allele convention (DECIDED 2026-09-12, LEDGER #131)

Which allele a number describes was answered ad hoc in five places before this
was written down. The rule, and the one deliberate exception:

1. **Every REAL variant is folded to its minor allele before it is tested.**
   `imputeGenoAndFlip` (`src/UTIL.cpp`) does `t_GVec = 2 - t_GVec` when
   `t_altFreq > 0.5`, on every marker of both the marker path
   (`src/Main.cpp:836`) and the region path (`:2746`). This is upstream's
   behaviour (`git show main:src/Main.cpp`, `:1269`) and is KEPT.
2. **Selection and weighting read the FOLDED MAF.** `--maxMAF_in_groupTest`,
   `--minMAF`/`--maxMAF`, `dbeta(MAF, a, b)` and the FlexRV MAF transforms all
   take `MAF = std::min(altFreq, 1 - altFreq)`. Upstream's too.
3. **A MAC is an allele count of the rarer allele** (#118), and the collapse
   criterion is that count: `MAC = min(altCounts, 2N - altCounts)` at `:2747`,
   then `MAC > --MACCutoff_to_CollapseUltraRare` decides whether a marker is
   collapsed. Three places broke this and are now FIXED: the pseudo-marker's
   band/ER decision read the WEIGHTED collapse sum (#134), the BURDEN
   emitter folded a cell's AGGREGATE count as if it were a frequency so its
   gating MAC fell as the cell grew past N (#135, fixed 2026-09-13: the gating
   MAC is the count itself, golden P66 at `--minGroupMAC_in_BurdenTest=1500`),
   and `MAC_case`/`MAC_control` were un-folded halves that did not sum to the
   folded `MAC` (#133).
4. **Reported frequencies and effects are per ALT allele.** `imputeGenoAndFlip`
   undoes its own fold for `altFreq`/`altCounts`, and the single-variant path
   negates the effect with `Beta * (1 - 2*flip)`. So an `AF_Allele2` above 0.5
   in output means "this marker was folded for testing", not an error.
5. **THE EXCEPTION: a collapsed ultra-rare pseudo-marker is NOT a variant, is
   never folded, and its reported MAC is its RAW carrier-allele count** (#130
   fixed 2026-09-12; the count half is #133, fixed 2026-09-13 -- `MAC` used to be
   derived from a FOLDED frequency while the vector it describes was not, so it
   disagreed with its own `MAC_case + MAC_control` above alt frequency 0.5).
   **Its MAF stays FOLDED and that is deliberate**: `MAFVec` is never written to
   any output, and its only consumers are the two `dbeta()` calls in
   `R/SAIGE_SPATest_Region.R` (`:324`, `:430`) -- i.e. the BETA weight. The
   pseudo-marker's only reported frequency is `AF_Allele2` in the single-variant
   file, which is `altFreq` and already raw. Making MAF raw would change the
   WEIGHTING, not the reporting: measured, it left every p-value byte-identical
   and blew `BETA_Burden` up to `-1.01393e+13`, because `dbeta(0.778,1,25)` =
   5.1e-15. Whether `dbeta` should see a collapse's CARRIER frequency at all is
   parked on #131 with that number attached. Do not "restore consistency" here --
   three measured reasons: `max_m(w_m g_m)` means "carries at least one of the
   set's ultra-rare variants" and folding inverts it to "carries NONE"; for a
   weighted set its entries are `w*g` (7.7 to 20 on the bundled multiweight
   fixture), so `2 - G` is not an allele complement and the old `altFreq > 0.5`
   predicate compared a weighted sum against a frequency; and in a cell holding
   both rare and ultra-rare markers the pseudo-marker carries 28-53 % of the
   MAC (golden P65), so folding it would set that share against the rest and
   make the burden a contrast rather than a sum.
6. **An ALT-major annotated site stays in its annotated bin**, weighted as rare
   on its folded MAF, with no warning -- upstream's behaviour, which has no
   guard for it anywhere (grepped). This is the DEFERRED half of #131 and the
   one part of this convention expected to change: at such a site the burden
   counts the REFERENCE allele while the annotation describes the ALT one, so a
   cell labelled "rare LoF" can be driven by a rare allele carrying no
   annotation. **The revisit trigger is one number** -- the fraction of
   annotated sites in a real callset whose ALT allele is major -- and it cannot
   be measured on the bundled fixture, whose 70 HARD-CALL markers are all below
   AF 0.5 (`--AlleleOrder=ref-first` inverts all 70 rather than giving a
   mixture). The one native exception in the repo does not help: `rs26` in the
   PGEN dosage goldens is AF 0.500034, i.e. ALT-major by 3.4e-5, which exercises
   the fold branch but separates nothing. The
   three candidate answers, and why the chosen one is (a) for now, are on the
   row.

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
silently leaves stale values. `Main.hpp` was DELETED on 2026-09-14 (audit F):
it was a third, unchecked copy of the step-2 signatures that re-drifted after
LEDGER #5 had aligned it (a differing parameter list is a new overload, so only
the declarations nothing depended on could drift); the compiler-checked
declarations are `RcppExports.cpp`'s, and the four internal forward
declarations sit at the top of `src/Main.cpp`.

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
  **Since 2026-09-16 the marker path asks the reader for a record in CARRIER
  SPACE first** (`GenoReader::getOneMarkerSparse`): `PgenReader` returns the
  difference list pgenlib stored (`PgrGetDifflistOrGenovec`, cap
  `2 * (N/8)`), `prepareSparseMarker` prepares it over the listed samples
  with the dense path's arithmetic (`imputeSparseAndFlip`, `src/UTIL.cpp`,
  states each identity), and the carriers are scattered into the worker's
  all-zero `GVec` for the test and cleared after. `PlinkReader` has no such
  form and reads densely as before. The REGION path does the same since the
  same evening (`prepareAndTest`: scatter only for a singly-tested marker;
  case/control tallies by rank), BURDEN 1.56x / SKAT-O 1.13x on sim50k. The run prints `Sparse reads: K of M
  ...` per chunk; `SAIGE_PGEN_SPARSE=0` is the kill switch. Measured on
  sim50k: quantitative scan 1.87x, binary 1.21x, every golden byte-identical
  (`docs/state/handoff-2026-09-16-pgen-sparse-read.md`). The zero-genotype
  index is no longer materialised on the marker path for either reader: its
  size is `m_n - nnz` at the test and the ER branch builds the complement.
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
- **A saddlepoint failure is never reported as the normal p** (LEDGER #163/#164,
  2026-09-26, found on real data as three genes at p 1e-120..1e-256). The Newton
  root finder overflowed when the carriers' fitted probabilities were of order
  1e-3 and below (first step `(q - K1(0))/K2(0)` past the exp underflow; fixed
  by a bracketed rescue that runs only after the plain Newton has failed, so
  no golden moved), and SAIGE then printed the normal-approximation p as
  `p.value` with only the sidecar's `Is.SPA = FALSE` to tell (which also reads
  TRUE after the exact test, LEDGER #80). Now: on failure the exact test runs
  at any weight set up to 13 carriers, the exact test's enumeration regime
  (`--max_MAC_for_ER` itself is capped at 200, LEDGER #165), every fallback is counted in `SPA fallbacks: N ...` at the end of
  every step-2 run and named in `<SAIGEOutputFile>.spaFallbacks.txt`,
  `SAIGE_SPA_PLANT_FAIL=1` is the test plant that forces it, and
  `SAIGE_SPA_DUMP_DIR=<dir>` writes the vector every failed saddlepoint saw
  (individual-level: analyse it where the data lives). A synthetic SPA
  vector must be centred the way `getadjGFast` centres: with W = mu(1-mu), so
  for tiny-mu carriers the negative mass is ~1e-2 and the MIRRORED q is
  unattainable -- a tail of probability 0, not a failure (LEDGER #168).
  Centring without the weights hides that and made the #163 rescue look
  complete when it converged almost nothing on real data.
  **The CGF itself is overflow-safe at a large root (LEDGER #169, 2026-09-27)**:
  `exp(g t)` in `Korg_*` overflowed past `g t` ~709 while K1/K2, written in
  `exp(-g t)`, stayed finite, so the finder converged where K was Inf and the
  saddle was refused -- the 20 BenIntNeo dumps, every one. The root is that
  large when the positive mass is dominated by carriers with tiny fitted
  probabilities and S is within a hair of gpos (`root * max|g|` ~ 1/mu of the
  dominant carrier). A refused q-side saddle is now labelled `saddle(q)
  invalid`, not `root(mirrored q) not found`.
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

**The SKAT-O kernel's outer integral is screened (LEDGER #158, 2026-09-18,
`R/SkatoScreen.R`).** SKAT's SKAT-O has two layers: eight per-rho Davies
p-values (4 ms per cell) and Lee's outer rho integral with Davies at every
quadrature node (80 ms, 95% of the kernel and ~70% of a SKAT-O run).
`--SKATO_screen_threshold` (default 0.05) evaluates SKAT's own closed-form
(modified Liu) version of the outer integral first and reports it as `Pvalue`
when it is above the threshold; below it SKAT's Davies path runs verbatim.
Measured: zero screened cells with a Davies p below threshold/10 over ~250k
null cells of three trait/prevalence arms and ~4.5k cells with signal; a
screened `Pvalue` is within 10% of Davies' (median 1.00 above 0.05); ~96% of
the kernel's Davies time gone, ~2.5x on a SKAT-O run. **The per-rho layer is
never screened, and must not be**: `Pvalue_SKAT`, `Pvalue_Burden`,
`BETA_Burden` and `SE_Burden` (the signed burden Z a Stouffer meta-analysis
reads) come from it, and Liu's per-rho burden value is 0.13x-4.6x Davies'.
`--SKATO_screen_threshold=0` is the upstream kernel byte for byte
(`test-step2-skato-screen.R` pins it, and pins that the run REPORTS its mode).
Provenance is the sidecar `<SAIGEOutputFile>.skatoMethod.txt`, one row per
tested cell (`liu-screen` / `davies` / `single-rho`), committed with the
goldens like `.pooledTests.txt`. **The tail below 1e-5 is the saddlepoint (step 2, same day)**: there SKAT as
shipped is 0.05x-7.5x off an exact reference, because its per-rho Davies at
an ABSOLUTE tolerance of 1e-6 fails and modified Liu is substituted;
`--SKATO_tail_method=saddlepoint` (default) recomputes both layers with
Kuonen's saddlepoint and takes the outer integral in its direct form,
piecewise between the kinks of the rho envelope (QUADPACK reports "extremely
bad integrand behaviour" otherwise), 5-12% conservative and stable on both
traits; `--SKATO_tail_threshold` (default 1e-5, REGENIE's) says where it takes
over, and 1 is how `test-step2-skato-tail.R` forces every bundled cell through
it against Davies, plus an exact synthetic control. Those cells' `Pvalue_SKAT`,
`Pvalue_Burden` and `SE_Burden` change too, as they should. No bundled cell is
below 1e-5, so the default tail moves no golden.
Ultra-rare few-case inflation (100 cases: 2.7x at 1e-4) is upstream of every
kernel option and is LEDGER #159, decided: left as is.

Three findings that outlive the fork path: **E0** chunk boundaries cannot reach
the output bytes (`markers_per_chunk` drives only a stdout counter); **E1**
`closeOutfile_single()`/`flushOutfile_singleInCPP()` exist because the output
stream and the `.index` flushed independently, which was LEDGER #47 (restart
duplicated rows; FIXED S9 by recording the file size per chunk line, `*.index`
is excluded from the byte gate) -- #49 is FIXED too (2026-09-05: `checkOutputFile`
returns the chunk to START AT, so a restart no longer re-runs the last completed
chunk; `R/Util.R:224` carries the comment). **#140 carried S9 to the REGION path
(2026-09-13)**: `flushOutfile_regionInCPP()` flushes the BURDEN and single-in-group
streams at each chunk boundary, the region index line records BOTH files' sizes
(`... chunk N single_bytes S bytes B`, writer `indexChunkLine()` beside its reader
`parseIndexChunkLine()` in `R/Util.R`), and a restart truncates the results file
and `.singleAssoc.txt` to them. What is still open there is only the two
end-of-run sidecars, `.pooledTests.txt` and `.skatFailures.txt`, which a restart
rewrites over the resumed regions alone (recorded on the row). **E2**
`varianceRatioBinsHaveGap` survives in
`R/MarkerScan.R` only because its message beats the C++ throw's.

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
- **Sample IDs are IID alone, and a repeated IID is refused at every step**
  (LEDGER #141, fixed 2026-09-13; #67 is the same rule for the phenotype
  file). FID is discarded everywhere, so two fam/psam rows sharing an IID are
  two individuals SAIGE cannot tell apart: steps 0/1/getNeff refuse in
  `readGenoSampleIDs` (`R/Geno.R`, both readers) and step 2 in each reader.
