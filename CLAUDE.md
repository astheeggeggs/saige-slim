# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**This is the `strip` branch**: SAIGE reduced to step 0 (createSparseGRM), step 1
(fitNULLGLMM), and step 2 (SPAGMMATtest: single-variant + region/gene tests),
with **PLINK bed/bim/fam and PGEN (plink2 native format) pgen/pvar/psam as
genotype input**. The BGEN/VCF/SAV readers, step 3 (LD matrices), and
RareEffect were stripped in commits "Strip 1"–"Strip 7", each verified
byte-exact against `main`-HEAD goldens; the PGEN reader (also stripped there,
in Strip 4) was added back in Phase 7 as fresh code implementing the
`GenoReader` interface (single-threaded until S8, which made
`PgenReader::clone()` real so `--nThreads` threads a pgen scan too), not a
restore of `main`'s pre-refactor `PgenClass`. **PGEN is no longer step-2 only:
step 0, step 1 and `extractNglmm` take `--plink2File` / `--pgenFile` /
`--pvarFile` / `--psamFile` too, and every output is BYTE-IDENTICAL to the
PLINK run on the mirrored fixture — see `docs/state/pgen-step01.md`.** `main` holds the full upstream fork (docs:
https://saigegit.github.io/SAIGE-doc/). Binary and quantitative traits are
supported; survival-trait support was removed in Phase 3a of the refactor.

**Do not read `LEDGER.md` whole.** It is 56 rows in two tables (#1-#39, then
#40 on, split by an audit section) whose rows are single lines of up to ~2,600
characters -- reading it costs a large slice of context before any work starts.
Query it instead: `python3 harness/tools/ledger.py open` lists the open entries
one line each, `ledger.py 56` prints one entry wrapped, `ledger.py grep ingest`
searches. Read `docs/state/next.md` (120 lines) to CHOOSE work, and pull the one
ledger entry you are actually working on.

**`LEDGER.md` and `docs/state/` are the authoritative record of where this branch
is** -- read those first, starting with `docs/state/handoff-2026-09-04.md` (the
most recent session handoff: what landed, the gate numbers, and where CI lives --
`origin` is the UPSTREAM `saigegit/SAIGE`, while CI runs on the `slim` remote,
`astheeggeggs/saige-slim`, so push there and query it with
`gh run list --repo astheeggeggs/saige-slim`); `docs/state/next.md` is the short list of what is NOT
done, with provenance marked on each item. There is a plan at
`~/.claude/plans/based-on-this-repository-silly-shore.md` (superseding
`rosy-munching-pudding.md`, which holds the Phase 0-5 history), but it was last
written 2026-09-01 and so predates S8-S12 entirely; treat it as history, not as
a statement of what is left. Plan files in that directory are written by
whichever session is running and there may be several -- do not assume the one
named here is the only or newest one.
Phases 0-3 are done (baseline/harness; strip; test suite; refactor — survival
removed in 3a, GenoReader interface + de-globalized Step2Config/Step2Output/
owned model in 3b-d, ten ledger fixes in 3e); Phase 5 (multi-weights) done
through C11; **Phase 6 done — LEDGER #28 and #29 FIXED** (326abe8f, 86a08f47:
the AF>0.5 marker path is now exactly flip-invariant, full tier 1225/1225).
**Phase 4 (speedup) is done, and so is the step-2 threading sequence S0-S12**
(see `docs/state/s8-pgen-clone.md`, which carries S8-S11, `docs/state/
s12-region-testcontext.md`, and the status block at the top of
`docs/design/step2-threading-verdict.md`): step 2 is measurable,
`--nThreads` threads the marker loop IN PROCESS for both readers, and the forked
path is retired. Step-1 blocked `GetTrace` also landed (`bda95327`) -- do not
believe an older doc that says it is parked on a branch.
**LEDGER.md still-OPEN**: #6 (PCG tolerance, deliberately deferred), #14
(unverified suspicion), #18-#26, #35-#39 from the type-polymorphic-builtin audit
(`sample`, `1:n`, `rbind` recycling) -- all latent behind a guard -- #45, #46,
#49 (a restart re-runs the last completed chunk: wasted work, right answer),
and #54 from the PGEN-for-steps-0/1 work (the unvalidated `.bed` header --
recorded, not reached by the new path), .
**#58 and #59 are FIXED** (found by the 2026-09-04 silent-failure audit,
reproduced by hand, then fixed): #58 the mmap `.bed` read had no bounds check
where the pread fallback did -- a `.bim` claiming one marker more than the `.bed`
holds gave BYTE-IDENTICAL output at exit 0, three more killed R with SIGSEGV --
now one size check in `setPlinkobj` naming both files and the shortfall; #59 a
region results file could be written with data rows and NO header when the first
counted chunk emitted none, now an `iswriteRegion` flag, which is what the
sibling `markerList` writer twenty lines below always had. Both are gated by
`tests/testthat/test-step2-file-integrity.R`, whose inputs are DERIVED at run
time -- a phantom `.bim` row, and a group file whose first region is
monomorphic -- because no committed fixture is internally inconsistent and none
has a barren first chunk. **#52 is FIXED**: `readPsamFile` keyed
the sample map on psam column 0 and called it IID, which is FID in a
`#FID IID SEX` psam, so step 2 died on any psam with FID != IID and silently
undercounted on repeated FIDs; it now finds `IID` by name from the header and
refuses a duplicate. All 13 PGEN goldens stayed byte-identical -- every
committed PGEN fixture has FID == IID, so no golden could ever have gated it,
and `tests/testthat/test-step2-pgen-psam.R` derives the psams it needs instead.
**#53 and #56 are FIXED, and they were the same bug.** `unsigned char geno2;`
in `Get_OneSNP_Geno_atBeginning` was read uninitialised and `setGenotype` is an
OR, so the first packed byte of EVERY marker -- samples 1-4 -- was built on an
indeterminate value. It was NOT masked, as #53 had assumed: at unpatched HEAD
the loaded store disagrees with an independent `readBin` decode of the same
`.bed` on 23 of the first 50 markers, and at N = 50,000 the byte varies BETWEEN
runs, which is #56 (six distinct `varianceRatio.txt` over eleven identical runs,
high-MAC category only). `getAlleleFreqVec()` is computed BEFORE the repack, so
the frequencies stayed right while the genotypes went wrong -- an internal
inconsistency, and the cheapest detector; that is what
`tests/testthat/test-step1-genostore.R` now asserts, alongside the independent
decode. **The fix moved four committed goldens, which WERE re-baselined**
(both step-0 sparse GRMs and both `binary_denseGRM` step-1 fits), each proved
byte-reproducible first by two independent `make_goldens.sh` runs plus `cmp`.
Every step-2 golden and every sparse-GRM step-1 golden was byte-identical
throughout -- the sparse fits read a frozen `.mtx` from `extdata/output`, which
the harness does not regenerate and which still carries the bug.
**#57 is FIXED** (2026-09-04): the chunked region path writes `P1Mat`/`P2Mat` to
`<SAIGEOutputFile>_P1Mat_Chunk_<i>.bin` and reads them back, and every
`arma::Mat::save`/`.load` return value was discarded -- Armadillo warns at
`arma_warn(3, ...)` and `ARMA_WARN_LEVEL` defaults to 2, so a failure was silent
by construction. A failed load skipped a `VarMat` block, SKAT threw for every
cell, and `get_SKAT_pvalue`'s `try(silent = TRUE)` turned the whole gene into
`NA` rows at exit 0. It stops now, naming the file and the region.
**The trap worth remembering: every non-BURDEN group golden ALREADY takes this
path.** Each bundled region has three or four non-ultra-rare markers, far under
`markers_per_chunk_in_groupTest`, so it looks like one chunk -- but the
ultra-rare collapse always appends another, so `nchunks == 2`. Do not reach for
`--markers_per_chunk_in_groupTest` to "force" chunking here; it changes nothing.
The FAILURE is what no golden reaches, and
`tests/testthat/test-step2-region-chunkio.R` injects it by planting a DIRECTORY
where the first chunk file must go.
Recently FIXED: #47 (S9, restart duplicated rows), #17 (S10, the uninitialised
`w0` on the custom-weighted conditional path), #42 (S11, three `log.p=TRUE`
arms halved a log, i.e. took a square root, inflating SE ~sqrt(2) wherever the
score p had underflowed) and **#50 (S12a)**.
**#50 is the twin of #27, in the branch #27 did not touch**: the non-BURDEN
write decision still read the mutated `m_flagSparseGRM_cur`, so an ESCALATED
region test parked its single-variant-in-group rows in the `.singleAssoc.txt_temp`
bin that the next region truncates -- **all 19 rows lost on the bundled fixture,
exit 0, no warning**. No golden could see it: `pval_cutoff_for_fastTest` is not a
CLI option here or upstream, and at its default neither bundled region escalates,
so P2 covered only the `copy_singleInGroup` route and P6 only the unconditional
write -- **P2 and P6 are byte-identical to each other in both files**, a pair that
looks like coverage and is not. Gated by
`tests/testthat/test-step2-region-escalation.R`, which forces escalation through
a direct `SPAGMMATtest` call and asserts the escalation HAPPENED before checking
the rows.
**#41 is FIXED**: the conditional SPA p-value was PRINTED halved
(`src/SAIGE_test.cpp:950` printed `SPApval_c/2` while `:929` built the SE from
the un-halved value), making every `p.value_c` from the SPA path -- and, via
`std::stod` at `src/Main.cpp:2019`, the `--r.corr=1` BURDEN Cauchy row -- 2x
anti-conservative; upstream since v1.1.5 (2022-08-17) and baked into every
conditional golden. Validated at 128,868 markers (10,085 halved) before the
one-line fix; goldens P3/P10/P12 + PGEN mirrors re-baselined (8 rows, `p.value_c`
only, each exactly doubled) and P25/P26 added for the `--r.corr=1 --condition`
path that had no coverage. The marginal ancestor of the same geometry was
upstream's SE bug in v1.1.6-1.2.0 (fixed in 1.3.0) -- see the historical note
in LEDGER.md; there the p was right and the SE inflated, the mirror image.
**#27, #33 and #34 are FIXED**, the three live ones, done together: #33 (the
binary conditional arm subtracted a diagonal instead of the matrix;
re-baselined goldens P3/P10), #34 (`diag()` handed a scalar variance -- an Inf
beta below 1, a terabyte allocation above it; fixed at the landmine so no
caller can re-arm it), and #27 (BURDEN + fastTest silently dropped every
escalated region; `m_flagSparseGRM_cur` was both the VR selector and the
escalation signal). **#40 too**: uninitialised reads in the conditional SPA
block, reached via the ER path by 2 of 70 markers on the default binary fixture
under `--condition`, one shifting a conditional p-value 30x anti-conservatively.
Golden **P12** `plink_marker_cond_binary.txt` was added with that fix -- the
conditional SPA block is gated on `traitType != "quantitative"`, so P4
(quantitative) could never reach it and NO golden covered the path, which is why
it survived. When a fix leaves every golden byte-identical, that is a statement
about coverage, not a pass.
Step 1 now takes `--seed` (default 1).

**Step 1 fits several phenotypes in ONE process** (Phase 2, commits N1-N4;
`docs/state/step1-multipheno.md`). `--phenoCol=y1,y2` (with `--traitType` and
`--invNormalize` per phenotype, or one value for all) writes
`<outputPrefix>_<pheno>.*` plus a manifest `<outputPrefix>.multiTraits.txt` that
step 2 consumes directly; ONE phenotype keeps the unsuffixed names and writes no
manifest. Phenotypes are grouped by complete-case sample set and each group does
ONE `setgeno`; the fits stay sequential. **Each phenotype's `.rda` and
`.varianceRatio.txt` are byte-identical to a standalone run's**, which rests on
`set.seed(seed)` before every fit -- the only RNG-dependent ingest input is the
`randi(1000)` variance-ratio pool. This needed **LEDGER #55** fixed first: a
second `setgeno` in one process used to APPEND to a dozen accumulators rather
than replace them (doubled GRM marker count, 2M QC flags, variance-ratio
sampling into the first ingest's freed arrays). Saving is the ingest only: ~1.4x
at P=5 on the N=1,000 sparse fixture, approaching (P-1) ingests at biobank M,
and ~1% on the DENSE-GRM path where the fit is 300-450 GRM passes.

**Step 2 can test several null models in ONE genotype pass** (Phase 1, commits
M1-M4; `docs/state/multimodel-step2.md`).
`--GMMATmodel_varianceRatio_multiTraits_File` takes a 3-column manifest (name,
`.rda`, variance-ratio file; no header) and writes `<SAIGEOutputFile>_<name>`
per model, each with its own `.index`. **Each model's file is BYTE-IDENTICAL to
that model's own single-model run** -- models whose step-1 `sampleID` vectors are
`identical()` share one reader, one decode and one QC pass per marker, and
models with different sample sets get their own reader; nothing is unioned and
nothing is gathered, so the identity is structural rather than a coincidence.
A **one-row manifest writes the unsuffixed file**, which makes golden P1 the
negative control for the whole feature. Marker path only: `--groupFile`,
`--is_rewrite_XnonPAR_forMales` and `--subSampleFile` STOP for more than one
model, each before any C++ call. In C++, `g_models`/`g_readers`/`g_modelGroup`
own the state and `g_model`/`g_reader` are non-owning aliases of element 0, so
the region path is untouched; `testOneMarker` is split into `prepareOneMarker`
(model-independent) and `testOneModel`. Measured at N=50,000: the shared prepare
is 67% (binary) to 91% (quantitative) of a marker's cost, so P=5 predicts
2.2x / 3.8x -- see `harness/bench/NOTES.md`.

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

Works natively on macOS and linux. `src/Makevars` still hardcodes
`-I../.pixi/envs/default/include`, so the in-tree `.pixi/` env is required.
`make` at the repo root runs the Docker build.

**plink-ng (pgenlib) is fetched and compiled automatically.** PGEN support
(Phase 7) needs `libplink2_includes.a`, built from a plink-ng source tarball;
it is not vendored in git (`.gitignore` excludes `plink-ng/`). `harness/build.sh`
fetches+compiles it into `.pixi/envs/default/lib` on first run (guarded on the
`.a` not already existing there, so repeat builds are a no-op), using the
pixi-env compiler (`$CXX`) so it works on both macOS and Linux; `docker/Dockerfile`
carries the equivalent layer for the Docker image. `plink2` (the CLI, bioconda,
`linux-64`/`osx-64` only -- no `osx-arm64` build, hence the per-platform
`[target.*.dependencies]` split in `pixi.toml`) is a pixi dependency too, but
only for generating PGEN test fixtures in the harness scripts -- never needed
at runtime by the installed R package.

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
- **standard** (~2.5 min, default): fast + `test-step1-genostore.R`, which
  checks the loaded GRM genotype store against an INDEPENDENT `readBin` decode
  of the same `.bed` and against its own `getAlleleFreqVec()` (the gate LEDGER
  #53/#56 slipped through: every golden was byte-reproducible while encoding
  wrong genotypes for samples 1-4, because every golden came from the same
  packer) + the PGEN-for-steps-0/1 checks in
  `test-step1-pgen.R` (the loaded GRM genotype store compared
  `expect_identical` between the two readers — `M`, `N`, `alleleFreqVec`,
  `MACVec`, `getQCdMarkerIndex`, 50 decoded markers — plus the step-0 golden
  and the dosage hard-call warning; its five-configuration step-1 byte
  comparison is **full** tier) + the four PLINK step-2 CLI runs
  compared against `harness/goldens/step2/`, the multi-weight runs, the
  marker-path goldens P27-P33 (`test-step2-marker-paths.R`, kept as its own
  file so its seven CLI runs go in parallel with `test-step2-goldens.R`'s),
  and the live Firth path (`test-step2-firth.R`, mirrored through the pgen
  reader in `test-step2-pgen.R`). `test-step2-region-escalation.R` (S12a) is the
  only test that reaches the fastTest ESCALATION pass of a region run: it calls
  `SPAGMMATtest` directly in a child process because `pval_cutoff_for_fastTest`
  is not a CLI option, and at its default neither bundled region escalates.
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
+15s of extra serial work. `test-step2-threads.R` (82s) is split in two the
same way over `helper-threads.R` -- units + P1/P4 + the refusal guards, vs
P7/P11 + the repeat-identity check (~38s / ~43s); `th_run` has no cache, so
that split is exactly assertion-neutral. The pole is now
`test-step2-multiweight-behaviour.R` at ~47s, against a 1:33 tier.

`Config/testthat/start-first` is set but **measured as no benefit both times it
was tested** (2:29.26 -> 2:26.17, and later 1:33.50 -> 1:42/1:41/1:50 -- the
second comparison confounded by ~2.5 cores of unrelated load). It is kept only
as insurance for when the file count grows; do not assume it is buying
anything, and note it names FILES, so it goes stale whenever one is renamed.
When timing this suite, check `uptime` first: a background Docker daemon at
~196% CPU moved the tier by ~10s, which is larger than most changes worth
measuring. Total CPU goes UP (223s -> 319s): every worker loads SAIGE afresh,
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
  **The same trap fired a second time, from the other end**: the multi-phenotype
  step-1 comparators in `helper-step1multi.R` shipped with a 1e-5 numeric
  tolerance, passed on the golden-generating machine, and went red on CI at rel
  err 4.2e-5 in `Y`. **Use the harness's tolerances** -- `harness/compare_step1.R`
  has always allowed 1e-2 for theta and the variance ratios (AI-REML runs at
  `tol=0.02`) and 1e-3 for other numeric fields (float PCG at `tolPCG=1e-5`). A
  test that demands more precision than the estimator has is measuring the
  platform's libm, not the code. Note CI was red for five hours on this and the
  session that pushed it recorded the tier as green -- locally, in byte mode.
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
- **Two of the 100 bundled markers ARE covariates.** `rs20` and `rs22` in
  `genotype_100markers.bed` are, genotype for level, the categorical covariates
  `a9` and `a10` (R² = 1.000000 against the model matrix; next highest marker
  0.0223). Any path that regresses the covariates out of the genotype therefore
  tests a vector that is exactly zero: both covariate-adjusting paths now report `BETA=NA SE=NA` (LEDGER #44 is
  FIXED, in two parts: an absolute `var1 <= DBL_MIN` guard for the sparse
  `scoreTest`, whose var is 0 exactly, and a RELATIVE cancellation test for the
  two `scoreTestFast` paths, whose var is ~2.7e-14 and so ~294 orders above
  `DBL_MIN`). The non-adjusting default path never degenerates and still reports
  them as ordinary null markers
  (LEDGER #44 — the `var1 <= min` guard covers the p-value only). The default
  `is_noadjCov=TRUE` path does not adjust, so P1 shows them as ordinary null
  markers and nothing noticed until S0. Expect two degenerate rows in any new
  covariate-adjusting marker golden, and do not read them as a signal.
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
- `bash harness/make_goldens.sh <dir>` — regenerates step0 (sparse GRM), step1
  (binary ×2 sparse / quantitative sparse / **binary dense** / **binary sparse-VR**),
  **the step-0 and step-1 PGEN mirrors of all six of those** (see below),
  and step2: 22
  PLINK runs (P1-P12 + P25 + P27-P34, with P8b as a
  separate composed-set golden — single-variant marker/marker-conditional
  ×2 (quantitative P4, binary P12), group fastTest/conditional/BURDEN/SKAT-O
  ×2 (K=1 P2/P3/P5/P6, multi-weight P8/P9/P8b/P10), chunked marker P7, and the
  ref-first characterization golden P11, the BURDEN-conditional golden P25
  added with the LEDGER #41 fix, the seven marker-path goldens P27-P33
  added by S0 — see the marker-path bullet below — and P34, the K=1
  CUSTOM-weighted conditional run added with the LEDGER #17 fix) plus 13 PGEN mirrors added across
  Phase 7 and its coverage-extension follow-up: P13-P15 (marker/region/
  conditional parity, mirroring P1/P2/P3) and P16-P24 (mirroring, in order,
  P4/P12/P5/P6/P7/P8/P9/P8b/P10 — quantitative conditional, the LEDGER #40
  binary-conditional path, BURDEN, non-fastTest SKAT-O, chunking, and the full
  multi-weight schema including composed sets and conditional multi-weight),
  plus P26 mirroring P25.
  Two PGEN fixtures were added with the steps-0/1 work:
  `extdata/input/plinkforGRM_1000samples_10kMarkers.{pgen,pvar,psam}` (the
  step-0 fixture's mirror, `--make-pgen` from its bed/bim/fam; verified at
  generation time that pvar ALT == bim col 5 and pvar REF == bim col 6 for all
  10,000 markers and that psam IID order == fam order — that correspondence is
  what makes the readers comparable, and is worth re-checking for any new
  fixture), and `extdata/input/dosage_20markers_100samples.{pgen,pvar,psam}`,
  20 markers × 100 samples whose ALT dosage cycles 0.1 … 1.0 so that 1,800
  genotypes are fractional and 200 exactly integral — it pins the dosage
  hard-call warning AND that an integral dosage is not miscounted as rounded.
  P11 has no PGEN mirror: `--AlleleOrder=ref-first` is a PLINK bim-specific
  toggle, and pvar's fixed REF/ALT columns give it no ref-first/alt-first
  ambiguity to mirror. Step-2 runs consume frozen model inputs: our committed
  step-1 goldens plus
  `extdata/output/example_binary_sparseGRM.{rda,varianceRatio.txt}` (the only
  legacy committed artifact still HEAD-compatible — the 2022-era ones crash,
  LEDGER.md #10). The PGEN fixtures
  (`extdata/input/genotype_100markers.{pgen,pvar,psam}` for P13-P19/P21-P24,
  `extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.{pgen,pvar,psam}`
  for the chunked P20) were each generated once from their bed/bim/fam
  counterpart via `plink2 --bfile <prefix> --make-pgen --out <prefix>` and
  committed; P13-P15 were byte-identical to P1/P2/P3 on this platform
  (verified at generation time), which is a stronger result than the
  numeric-tolerance reconciliation the test suite asserts — see the standing
  rule below about not asserting bit-identity between two different code
  paths as a general claim. P16-P24 are reconciled against their PLINK
  mirrors at float tolerance in `tests/testthat/test-step2-pgen.R` (same
  reasoning); they are not asserted to be byte-identical to the PLINK
  goldens, only to their own frozen PGEN goldens and to the PLINK mirror at
  tolerance.
- `bash harness/compare_goldens.sh harness/goldens <candidate-dir>` — byte
  compare for step0/step2; toleranced step1 (θ/varianceRatio 1e-2 relative —
  AI-REML `tol=0.02` makes tighter meaningless; vectors 1e-3). It globs
  `step1/*.rda`, so a new step-1 golden is picked up with no change; no test
  enumerates that directory, so adding one is purely additive.
- **The step-0/step-1 PGEN runs ARE committed as goldens**, even though each is
  byte-identical to its PLINK twin. Leaving them out was tried first and is
  wrong: `compare_goldens.sh`'s step0 loop has an **"EXTRA in candidate"** arm
  that flags any candidate file the goldens do not know about, so two
  uncommitted step-0 pgen outputs failed the full tier. (Step1 has no such arm —
  it globs `$GOLD/step1/*.rda` — so its five would have passed silently, which
  is the worse outcome of the two.) Suppressing the EXTRA check would have been
  disabling a safety net to hide a real asymmetry. `make_goldens.sh`
  additionally `cmp`s each pgen output against its PLINK twin **inside the
  candidate directory**: that is byte-exact, where `compare_goldens.sh`'s step-1
  arm is only toleranced, and `tests/testthat/test-step1-pgen.R` asserts the
  same thing from the suite.
  **`createSparseGRM` is not self-reproducing without a seed**: it draws 5,000
  of 9,588 MAF-passing markers with `sample()` and the `set.seed(1)` lives in
  `extdata/createSparseGRM.R:99`, NOT in the function — so a step-0 comparison
  that bypasses the CLI measures the RNG, not the reader (this actually
  happened; the control that caught it was that the PLINK arm did not reproduce
  itself either).
- **`binary_denseGRM` is the only step-1 golden on the DENSE GRM path.** The
  other three pass `--useSparseGRMtoFitNULL=TRUE`, so `getCrossprodMatAndKin`
  never leaves its sparse branch and the dense TBB kernel (`parallelCrossProd`,
  `GetTrace`'s per-probe `Au`) was covered by nothing at all until 2026-09-01 —
  a change to it would have left every golden byte-identical while being
  completely untested. Its config mirrors the `s1_full` bench workload so bench
  rows and golden runs describe the same code path. **Before accepting any new
  golden, prove it byte-reproducible at unpatched HEAD** (two runs, `cmp` on
  both the `.rda` and the `.varianceRatio.txt`); until that holds it gates
  nothing. Step-1 goldens ARE byte-reproducible on one machine — the 1e-2
  tolerance above is for cross-platform regeneration, so a same-machine change
  claiming bit-identity should be checked with `cmp`, which is strictly
  stronger than what `compare_goldens.sh` asserts.
- **P27-P33 close the step-2 MARKER-path coverage gap (S0).** Before them every
  single-variant golden (P1/P4/P7/P11/P12) took the same one arm of the
  `getMarkerPval` cascade — `scoreTestFast_noadjCov` — so the sparse `scoreTest`
  + PCG, the `scoreTestFast` (`--is_noadjCov=FALSE`) arm, the fastTest re-test
  block at `Main.cpp:436-475`, the Firth arm, the ER **resampling** arm and the
  unconditional quantitative path were all gated by nothing. Each is now one
  golden that is P1's command with ONE flag changed, with the per-arm counts
  measured on an instrumented build and recorded beside the runs in
  `make_goldens.sh`; `tests/testthat/test-step2-marker-paths.R` drives them.
  Two things worth keeping: the sparse marker path needs a variance-ratio file
  with **both** `sparse` and `null_noXadj` rows, which only
  `--useSparseGRMforVarRatio=TRUE` in step 1 produces (hence the fifth step-1
  golden `binary_sparseGRM_vr`, whose `.rda` is byte-identical to
  `binary_sparseGRM_a`'s); and `--max_MAC_for_ER=25`, the value the threading
  verdict doc proposed for ER coverage, is **byte-identical to P1** here (so are
  20/50/100) — reaching the resampling arm needs a marker with ≥21 carriers,
  i.e. `--max_MAC_for_ER=200` on this fixture. A golden that reproduces P1 gates
  nothing; check that a new one actually differs. **P33 is the gate for the ER
  seeding change**: `ComputeExact::Run` used to call `SL_setseed(1)` -- which
  enters the R API -- on every ER marker; it is now conditional on a resampling
  arm actually running, and P1 (all-exact) plus P33 (one resampling marker) are
  what make that checkable.
- `bash harness/tsan.sh [clean|plant|plantreader]` -- ThreadSanitizer over the
  in-process marker scan, the selected workloads at 1/4/8 threads
  (`SAIGE_TSAN_WORKLOADS`; `binary binary_sparse quant` by default, plus
  `rare100` and, since S8, `pgen`/`pgen_quant`). **Always run a plant mode too**:
  `clean` must report zero races, and `plant` (one `MarkerScratch` shared by
  every worker) must report races -- a TSAN build that sees none of this code
  also reports zero, so the clean run alone gates nothing. Four macOS-specific
  requirements are baked in and each of them fails by reporting NO RACES rather
  than by failing loudly (preload rather than link the runtime; `DYLD_*` cannot
  be set from outside, since SIP strips it through R's `/bin/sh` wrappers;
  bypassing those wrappers loses `DYLD_FALLBACK_LIBRARY_PATH`; and
  `R CMD INSTALL` test-loads and then DELETES the install, hence
  `--no-test-load`). The build goes to a copy of the tree, never `src/`.
  Expect the raw run to report main-thread-vs-worker races anyway: `libtbb` is
  not instrumented, so TSAN cannot see `parallelFor`'s join. `Main.cpp` declares
  that edge under `SAIGE_TSAN_ANNOTATIONS`, which only `tsan.sh` defines.
  **Two S8 repairs to this script, both of which had been failing silently.**
  (a) It drove the arm through `SAIGE_MARKER_THREADS`, which S7 deleted, so every
  "threaded" arm had been running serially; the knob is `--nThreads` now, and the
  script's own DID NOT THREAD check is what makes that visible. (b) pgenlib comes
  from `libplink2_includes.a`, which `harness/build.sh` compiles WITHOUT
  instrumentation, so a clean run said nothing about races inside the decoder --
  the same blind spot uninstrumented `libtbb` has. `tsan.sh` now builds an
  instrumented copy and, crucially, PREPENDS its `-L` in the package copy's
  `src/Makevars`: `PKG_LIBS` puts `$(LAPACK_LIBS)` first, so the env's
  uninstrumented archive would otherwise win the link and nothing would say so.
  Check with `otool -L .../SAIGE.so | grep plink2` -- it must name the absolute
  path under `harness/work/tsan`, not a bare filename. `plantreader` is the
  matching positive control for S8: `PgenReader::clone()` returns a handle onto
  the SAME `PgenClass`, which is the arrangement the `PgenShared`/`PgenClass`
  split exists to prevent (0 / 192 / 54 races at 1/4/8 threads with output
  DIFFERS, against 0/0/0 and identical for the same fixture clean). It is capped
  at 2,000 markers by `SAIGE_TSAN_MAX_MARKERS` because a shared decoder thrashes
  pgenlib's LD cache and the full-fixture arm does not finish -- the cap applies
  to the 1-thread arm too, so the byte comparison stays honest.
- `bash harness/make_controls.sh <dir>` — permuted-phenotype negative control
  through the full pipeline + KS uniformity check (`harness/check_null.R`),
  run through both the PLINK and PGEN readers (Phase 7) over the same
  permuted-phenotype model.
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
K=2 conditional run, and P34 `plink_group_cond_customweight.txt` the K=1
custom-weighted one. Reconciliation tests validate the R-side adjustment
against the K=1 baked path through TWO routes: the beta route, and — since
**LEDGER #17 is FIXED** — the custom route, which is sharper because
`group_new_snpid_withWeights.txt`'s weight line is byte-identical to the CADD
line of `group_new_snpid_multiWeights.txt`, so it is one statistic computed two
ways with no beta algebra in between. (Before the fix a K=1 custom-weighted
conditional run read an uninitialised `w0` on every collapsed ultra-rare row and
could not be a baseline for anything: measured, its conditional columns
disagreed with the K=2 CADD set on all 12 cells by up to 2.1e-01 relative, with
the adjustment vanishing outright on the 6 all-ultra-rare ones, while the
UNCONDITIONAL columns agreed exactly — which is what localised it. The fix is
`w0 = 1.0`, because a custom set's UR collapse is WEIGHTED, so the weight is
already inside the pseudo-marker; afterwards all 12 cells agree to 0.000e+00.)

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

- `src/Main.cpp` — step-2 bridge: a single static `std::unique_ptr<Geno::GenoReader>
  g_reader` (format-agnostic; see `src/GenoReader.hpp`/`src/GenoReaderRegistry.cpp`),
  ~25 `g_*` settings, `mainMarkerInCPP`, `mainRegionInCPP`, `Unified_getOneMarker`.
  The only per-format code left is `setPLINKobjInCPP`/`setPGENobjInCPP` (each
  builds a `spec` list and calls `Geno::makeReader("plink"/"pgen", ...)`) and
  four `if (t_genoType == "plink" || t_genoType == "pgen")` guards
  (`Unified_getOneMarker`, `Unified_getSampleSizeinGeno`,
  `Unified_getSampleSizeinAnalysis`, `closeGenoFile`) — grep for
  `t_genoType ==` before assuming that list is exhaustive if you add a format.
  Forgetting a `set*_GlobalVarsInCPP` call silently leaves stale values.
  `Main.hpp` declarations are stale (LEDGER.md #5) — don't trust them as a
  contract.
- `src/GenoReaderRegistry.cpp` — the two `GenoReader` implementations:
  `PlinkReader` (wraps `PLINK::PlinkClass`, `src/PLINK.{hpp,cpp}`) and
  `PgenReader` (wraps `PGEN::PgenClass`, `src/PGEN.{hpp,cpp}`, added Phase 7 —
  fresh code implementing the interface, not a restore of `main`'s
  pre-refactor `PgenClass`, which predates `GenoReader` and used out-param
  overloads instead of `MarkerMeta`). **`PgenReader` clones (S8), but not the
  way `PlinkReader` does** — and the difference is forced, not stylistic. A
  pgenlib read MUTATES its reader (fread buffer, LD cache, `fp_vidx`,
  workspaces), so clones cannot share one decoder the way PLINK's share one
  mmap. `PGEN::PgenClass` is therefore split in two: `PgenShared` holds
  everything read-only after construction (the pvar/psam metadata, the sample
  subset, the `plink2::PgenFileInfo` and its allocations) and is shared by
  `std::shared_ptr`; `PgenClass` is ONE decode handle, and `clone()` builds a
  second one against the same shared state. That is pgenlib's own mode-2
  contract — `PgrInit()` with a non-null `fname` copies the `PgenFileInfo` into
  the reader, the first reader inherits the header's open handle and each later
  one opens its own — and no read path writes the arrays that copy shares
  (`var_fpos`/`vrtypes`/`allele_idx_offsets`/`nonref_flags` are written only by
  `PgfiInitPhase2`). Two consequences worth knowing: the per-handle arena is
  ~130 bytes per SAMPLE (arithmetic from the allocation expression, not measured
  at scale — most of it two transpose buffers that single-variant reads never
  touch), so PGEN threading costs memory where PLINK's shared mmap does not; and
  `PgenReader::isForkSafe()` is still `false`, now as a live hazard rather than
  conservatism — a pgenlib reader owns a `FILE*` it fseeks, so a fork would share
  one offset.
  **`PlinkReader` clones (S5).** It holds a `std::shared_ptr<PLINK::PlinkClass>`
  and clones SHARE that one instance; `close()` releases only that handle, so
  the map and the `FILE*` die with the LAST owner and closing the primary under
  a live clone is safe by construction. This rests on
  `PlinkClass::getOneMarker` being `const` — it writes no member on either read
  path — so do NOT add a `mutable` member to work around it; that would silently
  remove the only thing making concurrent reads safe. `clone()` refuses when the
  .bed fell back to pread (conservatism about untested scale, not a hazard:
  pread carries its own offset). **The two .bed read paths.** mmap is the
  default; the `pread` fallback runs only if the map fails, which for a local
  file it does not — so it was executed by NOTHING until S5. Force it with
  `SAIGE_PLINK_NO_MMAP=1`; `tests/testthat/test-unit-genoreader.R` compares the
  two paths' decoded vectors (plus the LEDGER #16 mid-stream row-0 pattern) and
  `test-step2-marker-paths.R` reproduces golden P33 end-to-end through the
  fallback. Note step-2's `closeGenoFile` is **never called from R**, so before
  S5 the map and FILE* survived to process exit; `~PlinkClass` releases them now
  and the class is non-copyable so the destructor cannot double-close.
- `src/SAIGE_fitGLMM_fast.cpp` (~7400 lines) — step-1 engine, separate from
  Main.cpp. Single **global** `genoClass geno` holds the bit-packed GRM
  genotypes. **It now takes PGEN as well as PLINK, by RE-ENCODING rather than
  by a second pipeline.** `MarkerByteSource` has two implementations —
  `BedByteSource` (the original `seekg(3 + nbyte*i)` + `read`) and
  `PgenByteSource` (pgenlib decode via `PGEN::PgenClass::Read`, then re-encode
  as a PLINK-1 `.bed` record) — and `setGenoObjFromSource` is the single body
  both drive, so `Get_OneSNP_Geno_atBeginning`'s unpack, MAF/missingness QC,
  `round(2*altFreq)` imputation and 2-bit repack are ONE code path with two
  front-ends. That is why byte-identity between the readers is structural
  rather than a coincidence to be re-checked, and why the pack/unpack round
  trip is worth its cost. `setgeno_pgen` is the Rcpp entry beside `setgeno`.
  Steps 0/1 use the `PgenShared(pgen, psam, pvar)` overload (all file samples,
  identity mapping) and never touch `posSampleInPgen` — deliberately, which is
  what kept them clear of LEDGER #52 (fixed 2026-09-04) entirely.
  Algorithms exist as near-duplicate families — base, `_LOCO`, `_q`
  (quantitative) — after the `_Surv` variants were removed with survival
  support in Phase 3a.
- `src/SAIGE_test.cpp` — `SAIGEClass`: `getMarkerPval` cascade (score → SPA if
  |Z| > cutoff → ER exact test if binary & MAC ≤ `max_MAC_for_ER` → Firth),
  per-MAC variance-ratio selection. `src/Binary_*`/`src/ER_binary_func.cpp`
  are the live ER machinery (imported from SKAT) — keep.
- Step 0 cannot be separated from step 1: both drive the same global `geno`,
  and `createSparseKinParallel` is defined in `R/SAIGE_fitGLMM_fast.R`.
- **NEITHER step-2 path writes to the model** (S3 for markers, **S12** for
  regions). `SAIGE::TestContext` (`src/SAIGE_test.hpp`) carries the three values
  that used to be mutable members between a loop and the test -- the variance
  ratio, the sparse-GRM flag and the noadjCov flag -- and `getMarkerPval`, the
  score tests, `getadjGFast`, the PCG and the Firth fit are all `const`. S12
  gave the region path the same treatment at its three selection sites (the
  rare-marker loop, the ultra-rare loop, the BURDEN collapse) and in
  `assign_conditionMarkers_factors`, then **DELETED** `m_varRatioVal`,
  `m_flagSparseGRM_cur`, `m_isnoadjCov_cur`, `currentContext()` and the whole
  `assign*`/`set_*_cur` family. Build a context where the choice is made; the
  setters no longer exist to resurrect.
  **The fastTest ESCALATION signal is deliberately NOT a TestContext field.** R
  raises it between region calls to mean "this is the accurate re-run, DO write";
  sharing it with the per-marker sparse-GRM choice is exactly what LEDGER #27
  and **#50** were. It lives in `Main.cpp` as `g_regionEscalatedPass`, written
  only by `set_flagSparseGRM_cur_SAIGE`/`_org` and read once at the top of
  `mainRegionInCPP`. If you find yourself wanting it in the struct, that is the
  collision coming back.
- **The step-2 MARKER path's per-marker state is split by ownership** (S4,
  `src/MarkerLoop.hpp`). `testOneMarker` analyses one marker; everything it
  produces goes into a `Marker::MarkerResult` (one slot per marker, which the
  writer streams) and everything it reuses comes from a `Marker::MarkerScratch`
  (the loop-hoisted buffers -- shared across markers on purpose, so one instance
  per worker when this is threaded). The ~30 parallel q-length vectors are gone,
  and with them `std::vector<bool> isSPAConvergeVec`, whose bit packing made
  writes to two different markers writes to the same word. The two `std::vector
  <bool>`s that remain in `Main.cpp` are the REGION path, which is not threaded.
  Add a per-marker output by adding a field to `MarkerResult`; add a per-marker
  buffer by adding one to `MarkerScratch`. Do not add a q-length vector.
- One dataset per R process; no reentrancy — which is why step-2 `--nThreads`
  threads the marker loop IN PROCESS rather than forking (`R/MarkerScan.R`
  dispatches; see below). RcppParallel/TBB also drives the step-1 GRM kernels
  (source of LEDGER.md #7).

### Step-2 marker parallelism: `--nThreads`, one mechanism (S7)

**There is exactly one now.** `--nThreads` threads the C++ marker loop IN PLACE:
one process, one copy of the null model, one output file, and a result
byte-identical to a one-thread run. At `--nThreads=1` the loop is literally
serial and touches no threading machinery at all, so the default configuration
cannot be affected by any of it. `R/MarkerScan.R` validates and dispatches;
`mainMarkerInCPP` runs `RcppParallel::parallelFor` with an EXPLICIT `numThreads`
(which overrides `harness/env.sh`'s `RCPP_PARALLEL_NUM_THREADS=1` — without that
every threaded test would silently run serial and prove nothing).

The forked `mclapply` path that `--nThreads` used to drive (C2) was **retired in
S7** after both were measured: equal throughput at every ER-eligible fraction,
80–98 MB less resident memory at 4 threads, and restart works where fork refused
it. `splitMarkerShards`, `appendShardBytes` and `mergeShardOutputs` are gone. The
full evidence is in `docs/design/step2-threading-verdict.md` and
`harness/bench/NOTES.md`.

**Byte-identity is the right gate here**, despite the standing rule against
comparing two code paths: this is ONE path re-run over a disjoint partition, no
marker's arithmetic depends on another's, and the only cross-marker aggregates go
to stdout. `tests/testthat/test-step2-threads-inprocess.R` and
`test-step2-threads.R` hold it, and the former **asserts the run REPORTED its
thread count** — without that, a knob that failed to arrive would pass as a green
threaded run, which is the LEDGER #30 shape. Detection power is drill-verified in
both directions: sharing one `MarkerScratch` across workers fails it, and forcing
1 thread THROUGH TBB rather than down the serial shortcut passes, which separates
"TBB is involved" from "markers ran concurrently".

**Two behaviours are deliberate degradations, not oversights.** Group tests keep
the C++ marker loop on one thread -- but since 2026-09-04 `--nThreads` is not
DISCARDED there: it drives the SKAT-O post-processing instead (`R/SkatoParallel.R`),
which `Rprof` puts at 45% of a group run. Measured **1.30x at 4 workers**,
byte-identical, and 8 buys nothing over 4. Two traps are recorded in
`harness/bench/NOTES.md` and both cost a full measurement cycle: **a worker
function must be defined at TOP LEVEL**, because R serializes a closure with its
environment and one defined inside the driver shipped the whole `outList` per
task (0.97x, sys CPU 63s -> 30s when fixed); and `mclapply` per region is a net
slowdown, so the pool is forked ONCE for the scan. A reader that cannot be CLONED also drops to one thread,
with the reason both printed and raised as an R `warning()`. Degrading rather
than erroring is the agreed S7 decision: erroring would break anyone scripting
one `--nThreads` across mixed inputs. **Since S8 there is exactly one such
reader**: PLINK when the `.bed` could not be mmapped and the pread fallback is in
use (S5). PGEN threads — `PgenReader::clone()` hands each worker its own pgenlib
reader over shared metadata — and the pgen half of the byte-identity gate is in
`test-step2-pgen.R`, at 100 markers and again over the 128,868-marker fixture.

**Three findings that must outlive the retired fork path**, all about THIS
codebase rather than about fork:
- **E0** — chunk boundaries cannot reach the output bytes. `markers_per_chunk`
  only drives a stdout progress counter; golden P7 is byte-identical at 1000,
  1500, 2000, 3000 and 10000.
- **E1** — a forked child had to close `g_out.single` explicitly, because
  `mcexit` calls `_exit()`, no static destructor runs, and the child silently
  lost everything past the stream's 192 KiB buffer boundary *while still
  returning success*. `closeOutfile_single()` is byte-neutral without the fork
  and is kept, because that same buffer WAS the mechanism behind **LEDGER #47**:
  the `.index` and the output stream flushed independently, so restarting after a
  hard kill DUPLICATED rows (measured: 3,826 duplicates, identical serially and
  threaded — pre-existing, not a threading defect). **#47 is FIXED (S9)**: the
  chunk loop calls `flushOutfile_singleInCPP()` and records the output file's
  size in the index line, and a restart truncates to the size recorded for the
  chunk BEFORE the one it resumes at. The marker `.index` goldens were refreshed
  with it (chunk lines only; no `.txt` result changed, 54/54 byte-identical) —
  but `compare_goldens.sh` excludes `*.index` from the byte gate, so those files
  are a record of what the pipeline writes, not a check on it. The "before" is **LEDGER
  #49**: a restart resumes AT the last completed chunk and re-runs it, which
  contradicts `checkOutputFile`'s own message and costs one chunk of
  recomputation.
- **`varianceRatioBinsHaveGap`** survives in `R/MarkerScan.R`, gated on
  `nThreads > 1`, purely because its message is better than the C++ throw's. Its
  original rationale — that gapped bins made a threaded run unreproducible — was
  **retired by S3's fix to LEDGER #32**, and the message no longer claims it.


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
