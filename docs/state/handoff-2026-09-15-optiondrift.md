# Handoff 2026-09-15: LEDGER #151's guard, `harness/tools/optiondrift.py`

The ledger's five open rows are all decisions or an unreproduced crash
(#146b, #150, #151, #153, #154). The one piece of prescribed work among them
was #151's STATUS: the tool that stops the NEXT CLI-vs-API default
divergence. It exists now and runs before every tier.

## What it does

`harness/tools/optiondrift.py` (pure Python, no R, seconds) paren-matches
every `make_option(` block of the four `extdata/*.R` wrappers, the exported
function's formals (`fitNULLGLMM_multi` + `fitNULLGLMM` for step 1, the
multi's own winning on overlap), and the wrapper's call site to map flag to
formal (`opt$x`, `opt[["x"]]`, or a local whose multi-line assignment
mentions exactly one `opt$x`). Values are compared, not strings: `1L`/`1`,
`1e-05`/`1e-5`, a comma string against `c(...)`, and three coercions the
readers perform (`""`/NULL/`character(0)`; `"1"`/1 for the sex codes;
`"FALSE"`/FALSE for `--invNormalize`), each documented in `equivalent()`
with the reader line it rests on. An allowed divergence carries its ledger
row and BOTH values; a moved value fails the run as STALE, a divergence that
stops existing fails it too.

Controls, all drilled: `selftest` plants one divergence and three
representation-only agreements in a synthetic wrapper (caught / passed);
deleting an ALLOWED entry fails the run with `DRIFT`; changing one's value
fails it with `STALE`. `run_tests.sh fast` green with the wiring.

## What its first run found

Twelve divergences, not four. The eight new ones are on #151's STATUS and
the ALLOWED list. Two are LIVE in the arithmetic: step 1 `--relatednessCutoff`
0 vs formal 0.125 (the `getsubGRM` cutoff on a loaded sparse GRM), and step 0
`--numRandomMarkerforSparseKin` 2000 vs 1000 (the row had called that split
inert, which is true only of `getNeff`; the step-0 golden pins 2000 because
the CLI made it). The rest: step-0 `--nThreads` 16 vs 1,
`--sampleIDColinphenoFile` `IID` vs `""` (steps 1 and getNeff: an API caller
on the default indexes the phenotype file by an empty column name), and
`--outputPrefix` `~/` vs `""` (steps 0, 1, getNeff).

Not compared: the two FlexRV spec defaults (a named constant on both sides,
reported as `skip` under `--all`), and the `[default=...]` prose in help
strings and `@param` lines -- the `--X_PARregion` third value stays a
recorded, unguarded fact on the row.

## Gate

No `src/` or `R/` change; no tier owed. `ledger.py lint`, `docdrift.py`,
`optiondrift.py selftest`, `run_tests.sh fast` all green.

## What is still open

Every open row is the owner's: #146(b) all-zero weight matrix (recommend the
stop stands), #150 (code stays, surfaces corrected), #151 (twelve
per-divergence calls, now with a guard so the count cannot grow unnoticed),
#153 (`getNeff`'s `tau = (0, 1)`), #154 (one teardown segfault under load,
0/8 in isolation).

## Later the same day: #151 DECIDED, #153 FIXED

**#151.** All twelve divergences settled by the owner, one route per value,
three commits (`4cdf670b`, `81e3202f` and the tool's `d1d7237c`); the
ALLOWED list is empty. The live ones: step-1 `relatednessCutoff` formal 0
(in step 1 it is only `drop0(tol)` on a loaded `.mtx`, so it can prune step
0's pairs and never add one); `numRandomMarkerforSparseKin` formals 2000
(SE of an unrelated pair's kinship ~ 1/sqrt(M): the 0.125 cutoff is 3.9 SE
out at 1000, 5.6 at 2000; the owner runs 5000); `maxMAF_in_groupTest`
three cells; `isCovariateOffset` TRUE (help now says what it does);
`sampleIDColinphenoFile` `IID`; `outputPrefix` `""` (no stop: in a container
the only legitimate location is a mounted volume the user names, and `~/`
is the container's own home); step-0 `nThreads` 1; `X_PARregion` `""` on
both routes = GRCh38 PAR `10001-2781479,155701383-156030895` (GRC page,
checked) filled in with a banner and `warning()` only when the male rewrite
is on. One test moved: `test-step1-ai-projection.R` fitted through the API
without naming `isCovariateOffset`, and the new default shrank its
quantitative right-vs-wrong margins under the 1e-4 bar (kernel checks still
passed); it names FALSE now. That is the shape to expect from the default
changes: an API caller that relied on a default now gets the CLI's value.

**#153.** Not the formula. `tau = (0, 1)` makes Nglmm `1' K^-1 1`, the
effective number of unrelated samples, 402.4 for the fixture's 100 families
of 10. The 1e16 was a NON-CONVERGED PCG's partial sum printed as the answer:
the sparse GRM file was passed but not `--useSparseGRMtoFitNULL=TRUE`, so it
was never loaded and the solve ran on the dense marker GRM with no residual
term, where the inverse is set by near-zero eigenvalues (rank <= min(N, M)).
Upstream's `cmd_extractNeff.sh` always passed the flag. Fix: `getNeff`
requires the sparse GRM and always uses it; the dense path is refused, not
guarded, because no setting makes it well posed without changing the
question. Two-way check in the test: printed `Nglmm` vs `sum(solve(K, 1))`
from the same `.mtx`, 402.4175 both. Gate: see the commit.

**Open now:** #146(b) (recommend the stop stands), #150 (code stays by
decision), #154 (unreproduced teardown segfault).

## Evening: #146(b) FIXED, #154 drilled

**#146(b).** The owner's framing settled it: a discrete weighting scheme
legitimately scores some genes 0 everywhere, and that is the same as the
annotation not being present. Measured first: with ONE of GENE1's two weight
lines zeroed, that set's six GENE1 cells already vanished silently -- the
BURDEN emitter's rule that a cell with no non-zero weighted genotype is
neither refused nor tested, the same rule a bin with no markers gets. So the
all-zero matrix now composes its K sets like any other (`hasCustomWeights`
is `n_rows == q0` alone), every GENE1 cell is empty under that rule, the
pooled row is NA with no inputs and the sidecar names the gene; GENE2 is
byte-identical to the untouched fixture's run. Test inverted; gate FAIL 0,
no golden moved, as predicted. Still undisclosed, upstream's convention on
the whole region path: empty cells of any cause print nothing.

**#154.** Load recreated (24 concurrent 4-thread P3 runs beside a TSAN
build, load 12-23): 24 x exit 0, all 48 output files byte-identical. TSAN
region clean arm: races=0 at t=1/4/8, bytes identical to serial. 0/32 and a
clean detector is not a mechanism; the row stays as the re-run rule, with
the gate's own shape (eight testthat workers forking SKAT-O pools) as the
one untested difference.

**Open:** #150 (owner: delete the four dead chrX assignments or keep them)
and #154.

## Night: #150 FIXED -- the ledger is down to #154

Owner's decision: delete. `minInfo <- 0`, `minMAF <- 1/(2N)` (male-X
block) and `max_MAC_use_ER <- 0` (non-binary) were dead locals, pushed to
C++ before they were assigned and read by nothing; each is now a comment
naming the row, the corrected stdout sentence stays. Gate FAIL 0, no golden
moved, `test-step2-xnonpar-missing.R` unchanged. `ledger.py open` = #154.

## Later: #154 chased to a documented candidate; #155 found and FIXED on the way

**#154.** The `cleanup(kill = SIGKILL, detach = TRUE, shutdown = TRUE)`
frame is `parallel`'s exit finalizer `clean_pids`, registered at load and
run at EVERY R exit; `mc_cleanup` (R 4.4.3 `fork.c`) has no null
dereference, and the SKAT-O pool's nodes are forked estranged, so its
child list is empty and the finalizer has nothing to do. The main thread
was simply there when another thread faulted. Live candidate: oneTBB
2021.18 through RcppParallel 5.1.11.2, which installs no finalizer for it,
and oneTBB's spec says a process does not wait for its worker threads to
clean up at exit -- a documented window for a fault at address 0 after
every file is written. The remedy oneTBB provides (`task_scheduler_handle`
+ `finalize()`) has its own reported abort (oneTBB #1040), and a fix
nothing can reproduce cannot be verified, so nothing is changed. Evidence:
24 + 200 concurrent 4-thread P3 runs at load up to 60, 0 x 139; TSAN region
races=0 at 1/4/8. Stays OPEN as the re-run rule.

**#155.** The stress drill's one non-zero exit was a different defect:
`makeForkCluster` listens on one port fixed at package load from the RNG
and the clock (11000-11999), so two SAIGE processes on one box can collide
(~2.8 % per batch of eight, the gate's shape) and the second dies in
`startSkatoCluster` before any test. `startSkatoCluster` now retries a
failed bind up to ten times on a port derived from pid and attempt (no RNG
draw). `test-characterization-ledger155.R` holds a port in a helper process
and forces the run onto it with `R_PARALLEL_PORT`; the control is that the
first bind fails in the log. No golden moves (the first attempt is
upstream's call, verbatim).

## Night: the release image

`docker/Dockerfile` is two stages: `build` (the `default` pixi environment,
pgenlib, `R CMD INSTALL` into the `runtime` environment's library via
`SAIGE_R_LIBRARY`, then the toolchain conda's r-base drags in is pruned
THERE, because a layer never shrinks) and `runtime` (`ubuntu:24.04` + that
tree at the SAME path, since a conda environment is not relocatable, +
pgenlib's shared object at the path SAIGE.so records, + the four scripts
with their `pixi run` shebang rewritten to `Rscript`). `.dockerignore` keeps
`.pixi`, `harness/` (except `build.sh`), tests, docs and the example data
out. `pixi.lock` is committed so the build does not solve. Measured: 2.87 GB
with the default feature leaking into `runtime` (pixi includes it unless
`no-default-feature = true`), 1.94 GB with the prune in the wrong stage,
595 MB now; the image's step-2 run on the 100-marker fixture is
BYTE-IDENTICAL to golden P1 (macOS bytes, all 23 columns).
`.github/workflows/release.yml` builds on a `vX.Y.Z` tag equal to
`DESCRIPTION`'s Version, checks the P1 run at 1e-5 relative, and pushes
`ghcr.io/astheeggeggs/saige-slim:X.Y.Z` and `:latest`. Not yet exercised on a
runner: the first tag is the test. Also that day: the host disk hit 926 MB
free and corrupted Docker's store; the VM disk was deleted and rebuilt
(64 GB back), which is why the image was built from an empty daemon.
