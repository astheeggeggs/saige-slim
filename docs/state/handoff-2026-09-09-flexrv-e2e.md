# Handoff, 2026-09-09 (late): the FlexRV end-to-end suite, LEDGER #118

Previous handoff: `docs/state/handoff-2026-09-09-flexrv-f2-f4.md`. Design and
decisions: `docs/design/flexrv-design.md` ("What enters the test" has the new
paragraph, "Implementation notes" item 5 and the "End-to-end suite"
paragraph). This file says what a later session needs.

## What was asked, what was done

Asked: end-to-end tests for FlexRV, and whether the docs are up to date.

Done, in the working tree and NOT committed (the owner did not ask for a
commit; see "The checkout" for how to split it):

1. **Four new test files** over `helper-step2-flexrv.R`, which gained the run
   builders (they live in the helper because lintr's object-usage check
   cannot see another file's objects from inside a function body):
   - `test-step2-flexrv-reconcile.R`: all 192 sets against 192 external
     `weight:NAME` lines, binary with the collapse on (MAC from P1's golden)
     and quantitative with it off (MAC from the run's own `.singleAssoc.txt`,
     and the two MAC sources checked against each other); the MAC ladder
     against `--maxMAC_in_groupTest=1,3,5,10` cells of a plain run; the #118
     characterization, with a K = 1 custom 0/1 mask as the control.
   - `test-step2-flexrv-cond.R`: `--condition=rs30,rs79
     --weights_for_condition=1,1`, internal against external on the
     conditional columns, the pooled row's conditional CCT, the unconditional
     columns identical to the unconditional run's.
   - `test-step2-flexrv-paths.R`: `--flexRV_lofAnno=pLoF` on a relabelled
     fixture (bytes equal after the label is mapped back); `lof,missense`
     makes the `lof_*` rows the `unw_*` rows; `--markers_per_chunk_in_groupTest=10`
     and SKAT-O `--nThreads=4` byte-inert, each asserted to have arrived.
   - `test-step2-flexrv-multimodel.R`: a two-model manifest
     (`binary_sparseGRM_a` + `quant_sparseGRM`), 192-set BURDEN with
     single-in-group and the reduced-spec SKAT-O; each model's files
     byte-identical to its own run; the expansion and memory lines printed
     once per pass.
   The existing external-route test in `test-step2-flexrv.R` now uses the
   helper's builder. All six FlexRV files: 0 failures, ~1 min wall with 2
   workers on this (loaded) box. Every file lints clean.
2. **LEDGER #118** (below), the design doc, CLAUDE.md's FlexRV traps,
   README.md, next.md, the plan header.

## LEDGER #118, the one finding

A burden cell's variance-ratio band, ER eligibility, `--minGroupMAC_in_BurdenTest`
gate and printed `MAC` all read `genoSumcount_noweight`, the cell's UNWEIGHTED
MAC, shared by every weight set of the cell. So FlexRV's `MAC<=3` set on GENE1
(10 markers, real MAC 19, band 1) is tested with the cell's 121 (band 2) and
differs from the same 10 markers run as a `--maxMAC_in_groupTest=3` cell by
exactly the ratio of the two `null` variance ratios (1.02369) in BETA and SE^2;
where both MACs share a band (GENE2 at 3, every set at 5 and 10) the rows agree
to the printed digit. The `MAC<=1` set (4 singletons) is written although a
plain cell of MAC 4 is refused by the default gate of 5.

It is the engine's and upstream's semantics, not FlexRV's: a K = 1 custom 0/1
mask over the same markers reproduces the FlexRV row, not the plain cell, and
`main` reads the same count. Recorded, not fixed. A per-set count would move
P58/P59/P61 and needs a ruling on the weighted ultra-rare pseudo-marker.

## LEDGER #119: the owner's beta(1,25) question

Asked mid-session: can the engine's beta(1,25) be applied on top of a FlexRV
weight? Traced every evaluation site (design doc, "Implementation notes"
item 6): no on every K > 1 path, already bounded by the byte-inert
`--weights.beta=1,1` test. The one path that could have reached it -- a spec
selecting ONE set, K = 1, SKAT-O -- STOPPED instead, on LEDGER #19's
`stopifnot` guard in `region_postprocess`'s legacy fall-through. Fixed in the
working tree (`R/SAIGE_SPATest_Region.R`, both K = 1 branches treat a
transform set as a custom set), pinned by `test-step2-flexrv-k1.R`, no golden
moves (none has a single transform set); `test-step2-flexrv-cond.R` gained a
`--weights.beta` byte-inert check on the conditional run. This is an R-only
change: `R CMD INSTALL` (or `harness/build.sh`) before the tier sees it.

## Docs audit (the second question)

Up to date and mutually consistent before this session: CLAUDE.md's FlexRV
section, the roxygen in `R/SAIGE_Test_main.R`, `man/SPAGMMATtest.Rd`, the five
CLI help strings, the design doc (its two paper quotations were checked
verbatim against the PDF at the repo root), the plan header, next.md.

Stale and fixed: **README.md** had no FlexRV at all (its "Several weight sets"
section listed three set kinds and no `score:` line), and its multi-model
section still said region tests were refused with a manifest, which
`2ed00165` (2026-09-06) made false. Both corrected. README.md is not on
`docdrift.py`'s surface list; the multi-model line is the shape it catches.

## The checkout

Another session's uncommitted edits were present throughout and left alone:
`CLAUDE.md` (the SKAT-O served-count paragraph), `harness/sim/README.md`,
`tests/testthat/test-step2-skato-parallel.R`, plus the untracked `flexrv.pdf`.
CLAUDE.md now carries BOTH that hunk and this session's #118 trap bullet; a
commit of this work should stage only the FlexRV hunk (`git add -p`, or apply
the filtered diff with `git apply --cached`), not the other session's.

Suggested split: (1) the four end-to-end test files + helper + LEDGER #118 +
the design doc's #118 paragraphs; (2) LEDGER #119: `R/SAIGE_SPATest_Region.R`
(two K = 1 branches), `test-step2-flexrv-k1.R`, the `--weights.beta` check
in `test-step2-flexrv-cond.R`, the design doc's item 6 -- red at the
unpatched install (the K = 1 SKAT-O probe exited 1 with `identical(ws_kind,
"beta") is not TRUE`), green after `R CMD INSTALL`, no golden moved; (3)
README, CLAUDE.md (the FlexRV hunk only), next.md, plan header, this handoff.

Gates run: `harness/prepush.sh` green in both compare modes before the R
fix (4089 passes per mode, 2 warnings and 4 skips that predate this
session), and again after it -- see the end of this file for the second
result.

## NEXT UP

Unchanged from the previous handoff: **F5 on the 50k fixture** (wall time and
peak RSS of the 192-set BURDEN run against a one-set run, on a quiet box,
against a scratch install), then its controls. Two owner's decisions are now
queued on FlexRV: #117 (zero weight in a binary SKAT-O cell) and #118 (the
cell's MAC drives a set's test).

## Gate after the #119 fix (2026-09-10, 00:05)

`harness/prepush.sh` at 898e0b62 + this working tree, package reinstalled:
byte mode FAIL 0 / WARN 2 / SKIP 4 / PASS 4155, numeric mode FAIL 0 / WARN 2 /
SKIP 4 / PASS 4355, exit 0. The same 2 warnings and 4 skips appeared in the
gate run BEFORE the R change, so they predate this session. Every committed
golden is byte-identical, as predicted (no golden has a single transform set).

## F5 timing on sim50k (2026-09-10, after the gate above)

Asked next: can twin FlexRV sets skip their work, like the SKAT-O memo?
Measured rather than guessed; the full record is `harness/bench/NOTES.md`
("FlexRV F5 at N = 50,000"), the driver and analysis are
`harness/work/flexrv50k/{run.sh,analyse.R}` (gitignored). Findings: SKAT-O
twins are already served by the existing memo; on BURDEN the 192 sets cost
35 s / 805 MB vs 8 s / 506 MB for one set, the K-dependent C++ is 39%
ultra-rare pass + 18% emitter bookkeeping + 32% tests, twins are 42% of rows
at N = 50k and shrink with N, so a test memo alone is ~12% of a 1-thread run
and was decided against; the memory levers are the item. On the way, LEDGER
#120 (the FlexRV memory line left `std::cout` at one digit; fixed in
`src/Main.cpp`, stdout only, gate re-run). NOTE: this C++ fix was installed
into the SHARED `.pixi` library (`SAIGE.so` md5 `57285a58`), so other
sessions now run it; it changes no output file.
Gate after #120: `harness/prepush.sh` green in both compare modes (byte
FAIL 0 / PASS 4155, numeric FAIL 0 / PASS 4355, the same 2 warnings and 4
skips as before this session).

## FlexRV simulation study, S0 (2026-09-10, later)

Owner asked for a simulation study with realistic MAF and annotations at
large N, then proposed calibrating the annotation x MAF structure from real
summary statistics rather than the latent model. Design:
`docs/design/flexrv-sim-study.md` (decisions recorded: N = 50k, no scale
arm; half the causal genes directional; annotation from a real table). S0
landed and smoke-tested (see the design's "S0 status"): `--anno-table`,
`--frac-directional`, latent export in `02_build_dataset.py`; `run_all.sh`
knobs; `harness/sim/flexrv_study/{make_scores.py,run.sh}`. The chr20
correctness fixture's `.bed`/`.pgen` (5.5 GB) were deleted with the owner's
go-ahead to make room (`harness/work/sim/DELETED-2026-09-10.txt`; trees kept,
`FROM=02` restores it in ~30 min). The paper's own calibration checks were
read from `flexrv.pdf` and recorded in the design (no simulation study
there; permutation + synonymous + BHR; a per-set 5-carrier gate that
bears on LEDGER #118). BLOCKED on the owner's annotation + BMI summary
statistics files; the 50k rebuild (~70 min) runs once they are in.

## FlexRV simulation study: the table is built, the rebuild is running (2026-09-10 ~13:30)

Owner supplied All of Us v8 EUR BMI results (`harness/work/flexrv_study/aou/`,
gitignored): SAIGE gene results without MAC columns, the single-in-group
variant file (ultra-rare rows pruned, N = 171,473), and a chr20 per-variant
BRaVa annotation table with REVEL/CADD. `make_anno_table.py --anno ... --single
...` built `harness/work/flexrv_study/anno_table_aou_chr20.tsv` (the table and
its trend check are in the design doc, "The table"). The 50k rebuild with it
and `FRAC_DIRECTIONAL=0.5` started at ~13:30 (`rebuild50k.log`; ~70 min).

**Later the same afternoon** the owner's effect-model question led to three
more generator changes (design doc, "Effect model reassessed"): the
within-class score-effect link restored (`--within-class-sd`, mandatory),
30% causal genes (`FRAC_CAUSAL=0.3`) and per-gene link shapes
(`LINK_MIX=step,lof,graded`), all smoke-tested. The rebuild that was already
running (stage 01 done 13:49, stage 02 streaming at 14:05) picked up the
within-class fix but NOT the last two knobs, so its build stage must be
re-run. EXACT continuation, in order (every script under
`harness/sim/flexrv_study/`; `D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k`;
the rebuild's log is `harness/work/flexrv_study/rebuild50k.log`, it ends with
`exit=0` and a `real` time when done):

1. Checksum gate on the finished rebuild: every INVARIANT line of
   `harness/work/flexrv_study/sim50k-invariants-2026-09-01.txt` must match
   `md5 -q` of the rebuilt files (bed, bim, fam, pedigree, GRM, age/sex
   columns, `var` lines). A mismatch halts everything: the generator change
   consumed a random draw it must not have. (`y_null` is now expected to
   change: FRAC_CAUSAL redraws it.)
2. The final build stage:
   `cd $REPO; CONFIG=bench50k SIMDIR=$D NEED_GB=12 FRAC_DIRECTIONAL=0.5
   FRAC_CAUSAL=0.3 LINK_MIX=step,lof,graded
   ANNO_TABLE=$REPO/harness/work/flexrv_study/anno_table_aou_chr20.tsv
   FROM=02 TO=03b bash harness/sim/run_all.sh` (~35 min). Then the checksum
   gate again (same invariants), and read `02_build.log` for the printed
   class x bin table, the within-class correlations, the link counts and
   `63 causal genes`.
3. `bash harness/sim/flexrv_study/step1.sh` (four fits, ~10 min each, NT=4),
   then `python3 harness/sim/flexrv_study/make_scores.py --dir $D --out
   $REPO/harness/work/flexrv_study/scores` (miniconda python:
   `/Users/dpalmer/miniconda3/bin/python3`; the pixi env has no Python).
4. `make_scores`'s printed cor(score, true class) per arm is the first
   number to record; then `carriers2.py --dir $D --groupfile
   $REPO/harness/work/flexrv_study/scores/groupFile_informative.txt` (the
   collapse bound, BEFORE any p-value).
5. `bash harness/sim/flexrv_study/study.sh 3` (~136 step-2 runs, three at a
   time, resumable; refuses at load > 2.5 per run -- lower
   `MINLOAD_REFUSE` if the box is busy and say so in the write-up), then
   `Rscript harness/sim/flexrv_study/analyse.R` (controls first; an OUTSIDE
   flag on a null line halts interpretation).
6. Results: a section in `harness/bench/NOTES.md`, the design doc's status,
   next.md, memory. Nothing of the study is committed yet; the scripts under
   `harness/sim/flexrv_study/` and the generator changes are the commit.

## FlexRV simulation study: rebuilt, scored, blocked at the S1 gate on y_binary (2026-09-10 ~16:00)

Done this session, in the handoff's order:
1. The interim rebuild (no `FRAC_CAUSAL`/`LINK_MIX`) finished 14:21; every
   INVARIANT md5 matched. The FINAL rebuild (`FRAC_CAUSAL=0.3
   LINK_MIX=step,lof,graded`, `rebuild50k-final.log`, 29m43s) finished 14:57:
   bed/bim/fam/pedigree/GRM/var-lines identical, `y_null` changed as predicted,
   and **age changed too** -- located, not interpreted: `age` is drawn from the
   shared rng at `02_build_dataset.py:663`, after the causal-gene draws
   (452-598) whose count moved 21 -> 63; `sex` (line 311) reconciles exactly to
   the fam (pheno `sex` is a 0/1 male indicator, fam is 1/2). The invariants
   file carries the correction. Build log: 63 causal genes (32 directional;
   links step 18 / lof 17 / graded 28), class x bin table within 0.01 of the
   AoU table in every cell, cor(anno, true class) 0.696 latent / 0.477 cut.
2. `step1.sh`: four fits in **37 s total** (the "~10 min each" above was a
   guess; the `time: 1816.57` line in the logs is an ingest figure, not
   seconds). Gate "tau2 > 0": `y_null` 0.377, `y_null_bin` 0.331, `y_quant`
   0.434 pass; **`y_binary` ends at tau2 = 0** (trace 2.05, 1.74, 0.26, 0).
3. `make_scores.py` (`harness/work/flexrv_study/scores/`): cor(score, true
   class) on medium/low = informative 0.523, realistic 0.357 (target 0.3-0.4),
   uninformative 0.005, synonymous -0.002; on all coding variants 0.567 /
   0.426 / 0.125 (the 0.125 is the LoF = 1 anchor every arm shares).
4. `carriers2.py` (`carriers2.log`): 5,269 ultra-rare markers in the
   high:medium x 0.001 cell, 12,205 carriers, **34 carry >= 2 (0.28%)**;
   per-gene fraction quantiles 0/50/90/100 = 0 / 0 / 0 / 0.072. The 0-vs-10
   collapse comparison can only move a handful of genes.
5. `run.sh`/`study.sh` gained the `grid` method (owner's question: "will we
   get the standard thresholded encodings?"): SKAT-O over `high,high:medium` x
   maxMAF 1e-4,1e-3,1e-2, its Cauchy row is the standard-practice comparator;
   `analyse.R` reads it unchanged (Cauchy row = p, best-cell Bonferroni =
   "oracle"). Study is now 144 runs. NOT started.

**The block, characterized with controls (all in the scratchpad `s1/`, ~10 s
each; the four study models in `$D` are untouched):**
- `y_binary` from `--tauInit=1,0.5` and `1,0.1`: both collapse to 0.
- NEW `y_binary_pc` (h2 = 0.2 positive control; 27 case-concordant
  first-degree pairs vs 0.94 expected): trace 3.7, 3.2, 0.6, 1.7, 1.8, 1.3,
  0.5, 0.6, 2.2, 2.0, 0.8, **0**. Not information-starved: `y_binary` has 16
  concordant pairs, `y_null_bin` 5, and `y_null_bin` converged.
- OLD `y_binary_pc` phenotype (rebuilt from the 2026-09-01 `.rda`'s `y`/`X`)
  under the CURRENT code: **1.848389, the 2026-09-01 value to 7 digits**.
- NEW `y_binary` under `harness/work/lib-old` (ddedfa6d, its own CLI): the
  identical trace to 7 digits, 0.
  => not a code regression; the new binary phenotypes send the AI-REML
  through a step that lands below `tol` = 0.02. Mechanism, two lines:
  `src/SAIGE_fitGLMM_fast.cpp:3733-3737` floors any `tau < tol` to 0 BEFORE
  the step-halving loop at 3739-3745 (`while (tau < 0) step *= 0.5`), so a
  negative tau is already 0 when the loop tests `< 0.0` and **the halving is
  unreachable**; then `R/SAIGE_fitGLMM_fast.R:213` `if (tau[2] == 0) break`.
  One AI step that overshoots to below 0.02 ends the fit (upstream code,
  `float` arithmetic; candidate LEDGER row, not written). Every binary trace
  here oscillates by 3-5x between iterations (also the converged ones), i.e.
  the AI step is ill-conditioned at 1% prevalence with 15% of samples in
  families; the 2026-09-01 draw converged, this one does not.

**Owner's decision (not taken):** (a) run `y_binary` with tau2 = 0 (a GLM,
no relatedness adjustment, VR = 1; the method comparison stays fair since
every method shares the model, but the binary power arm is no longer the
mixed-model setting and `y_null_bin` IS); (b) a LEDGER row + fix for the
absorbing boundary (move the floor after the halving loop so the halving
can run, as it was evidently meant to), gated as any fitter
change (`binary_denseGRM` etc. goldens must not move -- none of them hit the
boundary, so the fix moves no golden, i.e. coverage, not a pass); (c) a
fixture with tau identifiable on a binary trait (more families or higher
prevalence; a rebuild). The three fitting traits can run the sweep now
(`TRAITS="y_null y_null_bin y_quant" bash study.sh 3`, 108 runs), and the
first thing `analyse.R` prints is the null ECDF controls.

## LEDGER #121 fixed; the S1 gate passes for all four traits, with a caveat (2026-09-10 ~18:00)

Owner's decision: fix it. Two commits, the branch's pattern: `5281de8d`
extracts the binary AI-REML tau update into the Rcpp-exported
`tauStepBinary(tau0, Dtau, tol)` and pins today's rule in
`test-unit-tau-step.R` beside an R transcription of GMMAT's; the fix commit
gates the first floor on `tau0 < tol` (GMMAT's rule, and what
`fitglmmaiRPCG_q`'s `zeroVec` already did), inverts the two #121 tests
(red at `5281de8d`, drilled through a scratch install), and every gate
held: 30 step-0/step-1 goldens byte-identical at both commits (regenerated
into scratch and `cmp`'d), fast tier green, prepush (see the ledger row).

**The study's models were refitted with the fix** (`step1.sh`,
`harness/work/flexrv_study/step1-fix121.log`): `y_null`, `y_null_bin`,
`y_quant` byte-identical to the pre-fix files (`$D/prefix121/`, kept);
`y_binary` now tau2 = 1.31, `converged = TRUE`, VRs 0.98. The collapsed
tau2 = 0 model is kept in `$D/collapsed/` for the inflation measurement.

**Caveat that outlives the fix:** every binary trace on this fixture swings
3-5x between consecutive iterations (`y_binary`: ... 0.38, 0.62, 1.34,
1.31; the h2 = 0.2 trait runs to `maxiter` = 20 unconverged at 0.19),
identically at every `--traceCVcutoff` (the CLI does not expose it; the
sparse path is deterministic), so it is the (alpha, tau) fixed-point
iteration itself, not trace noise. `y_binary`'s "convergence" is two
neighbouring iterates within 2%, on an oscillating sequence. The value is
usable for the method comparison (every method shares it) but is not a
reliable tau2. Open question in `next.md`.

Continuation: the sweep, `bash harness/sim/flexrv_study/study.sh 3` (144
runs, resumable), then `analyse.R`; the collapsed-vs-fixed inflation on
`y_binary`'s non-causal genes is one extra `run.sh` pair (point `D` at a
directory holding the collapsed model).

## #121 null series: "drop to logistic" vs the GMMAT rule (2026-09-10 evening)

Owner's challenge: a rescued tau might be worse than a silent logistic
fallback; check null fits before settling the fix. Done, in
`harness/work/flexrv_study/nullseries/` (`pheno_series.txt`, the fits under
both rules, `s2_*.txt`, `analyse_series.R`):
- Null traits (y_null at 1/2/5/10%, permuted 1%): the two rules are
  byte-identical on every one. The boundary is only reached with a strong
  signal (the family-shared big-effect variants of `y_binary`).
- Permuted 10%: tau2 = 0 in one step. Permuted 1%: 0.13, converged. The
  estimator finds a true zero when the trait identifies it; at 1% it parks
  at a floor, so an old-rule zero cannot be read as "no relatedness".
- Same `y_binary`, collapsed (0) vs fixed (1.31) model, study cell: median
  per-gene log10 p shift +0.005 null / +0.027 causal; fixed model inside the
  binomial CI on all six null lines, collapsed outside on two; 9/63 causal
  detected under both. The fix stands; at this relatedness fraction the
  tests barely notice tau, and the fixed model errs conservative.
- Also seen: a mild excess at alpha = 0.05 for SKAT-O on every 1% model
  including the permuted control (0.076, CI [0.044, 0.121]), n = 210 genes;
  not resolved here.

## Convergence of the binary AI-REML: what helps (2026-09-10, late)

Owner's question after #121: can anything raise the likelihood of
convergence? Measured on `y_binary` (all under the fixed rule; logs and
the prototype driver `step1_damped_v3.R` in
`harness/work/flexrv_study/nullseries/`):
- **CLI knobs do nothing.** `--tauInit` 0.1/0.33/0.5 all enter the same
  period-4 cycle (~1.35 high, 0.2-0.4 low); `--maxiter=100` runs the cycle
  longer; `--tol=0.005` needs 28 iterations to catch two neighbours; the
  h2 = 0.2 trait with `--maxiter=40` cycles 24 times and then dies at 0
  anyway (a halved overshoot under `tol`, allowed by GMMAT's rule too). So
  #121 removes the single-overshoot death, not the cycle.
- **Damping the AI step collapses the cycle.** Prototype patched into the
  loaded namespace (no install): `tau2 <- tau0 + d * (proposal - tau0)`,
  `d = 0.5` fixed, or adaptive (halve `d` when the step changes sign), with
  the convergence criterion evaluated on the UNDAMPED proposal (on the
  damped step it is trivially met early: 0.90 / 1.00 / 1.31 by `d`).
- **Two self-consistent states, resolved by tolerance.** A one-AI-step
  "score sign" from cold starts brackets a root near 1.45, and a damped
  run started at 1.45 "converges" in one iteration at 1.45 at tol 0.001;
  warm paths from 0.3 / 1.0 / 2.5 go to 0.97. At tol 0.0002 every start
  including 1.45 ends at 0.9696-0.9705 in 4-20 iterations: the joint
  (alpha, b, tau) fixed point is **0.970**, and the cold-state root is the
  inner pass not yet converged. SAIGE's undamped fitter reports 1.31, 0 or
  0.19 for this trait depending on the budget; none is the fixed point.
- **Same on the h2 = 0.2 trait:** damped 0.5 / 0.3 / adaptive at tol 0.005
  give 1.32 / 1.32 / 1.34 where undamped gives 0 or 0.19.

What would change in the code (owner's decision, a ledger row, gated on
the step-1 goldens whose binary fits converge in 4 undamped iterations and
would take a different path): damp the step (adaptive sign-flip halving
needs no constant), test convergence on the undamped proposal, tighten
`tol` for tau (0.02 is the criterion that let a cycle pass) and require it
twice running. The variance ratio of the 0.97 model vs the 1.31 model
differs by well under 1%, so for the study the ranking is unaffected; the
right value matters where relatedness is dense.


## Estimator tests for step 1: design and the measurements behind it (2026-09-10, night)

Owner's questions: can "the algorithm works" be part of the suite, and is
the quantitative fitter affected too? Everything below is in the
scratchpad `est/` (planted phenotypes, references `pql_ref.R`, SAIGE
logs) and is NOT yet a test; the design is what to build.

**Fixture.** The bundled step-1 fixture is 100 families of 10, everyone
related (2,700 pairs, mean off-diagonal 0.38), N = 1,000, with the sparse
GRM in `extdata/output`. Phenotypes are planted from that GRM's Cholesky:
quantitative `y = 0.5 x1 + 0.3 x2 + b + e`, `b ~ N(0, tau2 G)`, tau2 in
{0, 0.03, 0.1, 0.3, 1}; binary the same on the logit scale with intercept
-1.5, tau2 in {0, 0.5, 1, 2} (197-298 cases). Seeds 31/32.

**References, computed in dense R at N = 1,000 (milliseconds):** exact
REML for the quantitative traits (profile over tau2/tau1 with Brent);
for the binary traits the same PQL AI-REML iteration as GMMAT's
`glmmkin.ai`, inner loop to 1e-8 and a bracket-and-bisect safeguard so
the reference cannot cycle. **The reference must use SAIGE's covariate
convention**: the CLI default `--isCovariateOffset=TRUE` fits the
covariates in the preliminary GLM and carries them as an offset, so the
mixed model has an intercept only (owner: a speed feature); with the flag
FALSE the model is the full one. Both were replicated.

**What the measurements say (SAIGE tight: tol 1e-4..1e-8, PCG 1e-9, trace
CV 5e-5 / 500 probes):**
- Quantitative: within 0.2-3% of exact REML on all five traits (absolute
  0.0006-0.004). Binary: within 1-3% of the PQL fixed point under the
  offset convention. At SAIGE's converged state the reference reproduces
  its linear predictor to 1e-7 and its `YPAPY` to four decimals; the
  whole residual is the TRACE: SAIGE's Hutchinson estimate (Rademacher
  probes, `set_seed(200)` on every call) sits 0.20 above the exact
  tr(PG) with 500 probes, where the exact per-probe sd (from tr((PG)^2))
  gives an se of 0.28. The default 30 probes give se(tau2) = 0.0245 on
  a binary tau2 of 0.13 (19%) and 0.0050 on a quantitative 0.14 (3.5%),
  and the seed makes the error the same sign in every fit and iteration,
  so it reads as a deterministic offset (+0.005 on every quantitative
  trait at defaults, closing to +0.0003 at 300 probes, as predicted).
- So: both fitters reach the exact fixed points up to the trace noise.
  The PQL bias itself is visible (planted 0.5/1/2 -> 0.14/0.29/0.63) and
  is the estimator's, which is why the test compares against the
  algorithm's fixed point, never the planted value for binary traits.
- Quantitative cycling: none. The 50k weak-signal series (y_null mixed
  with its permutation at a = 0.5..0.05) converges from cold and warm
  starts to the same values (0.056; 0 for a <= 0.3); a first "expected
  a^2 x 0.377" table was an assumption that ignored invNormalize and
  standardisation and was retracted the same hour.
- Found on the way: `getAIScore`/`getAIScore_q` apply the covariate
  projection to the solved vector (`Sigma_iXt * PAPY_1`) instead of to
  `APY`, so the AI curvature is off by 0.1-0.8% here (LEDGER #124);
  step lengths only, no effect on any fixed point.

**Test design (standard tier, ~10 SAIGE fits at N = 1,000, ~1 min):**
1. `helper-step1-planted.R`: plant the phenotypes from the bundled GRM
   (fixed seeds), write the pheno file, and hold the dense references
   (REML exact; PQL under both covariate conventions) plus, per trait,
   the exact `se(tau2) = sqrt(2 (tr((PG)^2) - sum diag(PG)^2) / nrun) /
   AI` at the reference fixed point.
2. `test-step1-estimator.R`: run step 1 with `--nrun=300
   --traceCVcutoff=1e-4 --tol=1e-4` and assert `|tau2_SAIGE - tau2_ref|
   <= 3 se + 1e-3` for every planted trait, both trait types, both
   covariate conventions; assert the planted-zero traits come back at 0
   or within 3 se of 0; assert monotonicity in the planted value.
3. Controls in the same file: a negative control that the comparator
   FAILS when the reference trace is replaced by tr(Sigma^-1 G) (a
   planted wrong formula moves the binary fixed point by 12%, far
   outside 3 se); a positive control that the reference reproduces
   SAIGE's YPAPY at SAIGE's own state (the reconciliation above).
4. Kernel-level (fast tier), once #122 lands: the safeguarded Newton as a
   pure function on synthetic score functions where plain Newton is
   proven to cycle.
The 50k cycle (#122) is not reachable at N = 1,000 with everyone related;
the test above pins the fixed point, #122's own test pins the iteration.

## LEDGER #122 implemented: adaptive damping, not a bracket (2026-09-10, night)

Owner's go-ahead with the constraint "nothing worse than what we have". The
bracket design's first drill refuted it (stale score signs from an
unconverged inner pass; details in the row), so the guard is adaptive
damping armed on the second direction reversal, zero as an iterate,
convergence on the raw proposal, never on iteration 1, twice once a step
has been damped. Every healthy fit is byte-identical (30/30 goldens, traces
and iteration counts identical, 50k nulls identical); the 50k trait
converges from four starts to 0.968-0.973 (fixed point 0.970) where it
reported 1.31 / 0 / 0.19. Two design mistakes were caught by the tests, not
by me: the boundary exit crashed the post-loop refit (state not handed
over; my scratch drill had only grepped the printed tau), and keying the
twice rule on "armed" cost the dense golden an iteration. Files:
`R/TauStepGuard.R`, the loop in `R/SAIGE_fitGLMM_fast.R`, `tauProposalBinary`
in `src/SAIGE_fitGLMM_fast.cpp`, `test-unit-tau-guard.R`,
`test-step1-guard.R`. #122 is commit 74c237d3; #123 (the `converged` flag
from the exit taken, not the iteration count) follows it as its own commit;
#124 (AI projection slip) still open.
