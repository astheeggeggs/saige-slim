# What is worth doing next (2026-09-03, updated 2026-09-12; read docs/state/handoff-2026-09-12-review-fixes-and-gate.md first)

**Current handoff: `docs/state/handoff-2026-09-12-ur-pseudomarker-129-130.md`**
(LEDGER #129 and #130 FIXED, one commit each -- the two decisions and the
reasoning that settled them, the exact three-way reconciliation #130 bought, the
18-file blast radius predicted before the gate, and LEDGER #131 opened by the
same question)
(previous: `docs/state/handoff-2026-09-12-review-fixes-and-gate.md`,
the adversarial review's three fixes -- one of them live and producing wrong
numbers -- the CCT answer, and the gate at 1.97x)
(previous: `docs/state/handoff-2026-09-12-flexrv-ram-levers.md`, the
five RAM/CPU levers on the FlexRV region path, one gated commit each; before
it `docs/state/handoff-2026-09-12-flexrv-rerun-129.md`)
(previous: `handoff-2026-09-12-126-null-calibration.md` -- the #126 null
calibration and #128; then `handoff-2026-09-11-evening-117-115-114.md` -- what #118 cost and
bought, and the disk clear-out; then `handoff-2026-09-11-percell-mac-implemented.md`, which is what
#118/#126/#127 ARE and how they were gated -- still worth reading; then
`handoff-2026-09-11-percell-mac.md` -- the DECISIONS handoff for
#118/#126/#127, still the place to read WHY each choice was made, but two of its
predictions are corrected by the one above; then
`handoff-2026-09-10-flexrv-sweep-124.md`, then
`handoff-2026-09-10-fitter-guard.md`, then
`handoff-2026-09-09-flexrv-e2e.md`, then `handoff-2026-09-09-flexrv-f2-f4.md`, then
`handoff-2026-09-09-flexrv-f0.md`, then
`handoff-2026-09-09-skato-identical-cells.md`)
-- this line is the pointer and is updated with every new handoff, because
several handoffs can share a date and `ls | tail -1` sorts by NAME, not by time.
The three 2026-09-08 handoffs are superseded; of those,
`handoff-2026-09-08-coverage-loco-flexrv.md` still reads correctly as history but
its "NEXT UP: eight open ledger rows" is out of date (five are now closed), and
`handoff-2026-09-08-ledger-zero-and-coverage.md`'s claim that the 2026-09-01
`step2_group_y_null.txt` is "a genuine BEFORE baseline" is WRONG -- ~80 commits
separate it from HEAD. `docs/state/handoff-2026-09-04.md` is where CI lives.

**2026-09-12 (morning): the FlexRV sweep was RE-RUN at HEAD** (prediction `docs/state/prediction-flexrv-rerun-2026-09-12.md`, outcome appended; NOTES.md "the sweep re-run at HEAD"). The beta-weight SKAT-O lines reproduce 2026-09-10 to the digit (the pipeline control); every FlexRV detection count fell by 1-2 through #118 (flex192 now 11 = grid 11 > unw 8 on `y_binary`). It found **LEDGER #129 (FIXED 2026-09-12)**: the dosage-zeroing filter was applied to the WEIGHTED ultra-rare pseudo-marker in weight units, so transform sets with weights <= 0.2 vanished at collapse 10 and partially-zeroed sets were tested on the surviving carriers (`GENE0181 sbeta10`: 1.7e-3 at collapse 0, 2.9e-9 at collapse 10, 1.7e-3 again with `--dosage_zerod_cutoff=0`). The filter no longer touches the pseudo-marker (fix (b) of the two on the row). **The sweep's collapse-10 columns are therefore stale and the collapse-0 caveat is now MOOT** -- the numbers in `prediction-flexrv-rerun-2026-09-12.md` and the study memo were measured with the bug in, and a re-run is the only way to know what they are; the per-trait presence counts (20/45/44/21 sets present at collapse 0 and absent at 10) should now be ZERO, which is the cheapest check that the fix reached the sweep and not only the bundled fixture. The bundled blast radius was predicted exactly (18 goldens moved, 1 deleted) and is on the row. The three study [DECISION]s stand; the memo is `docs/state/flexrv-study-decisions-2026-09-12.md`.

**2026-09-12: the #126 null-calibration check is DONE and it passed** (`docs/state/prediction-126-null-calibration.md`, prediction before the runs, outcome appended; summary on the #126 row). The MAC 1-4 cells the default now admits are uniform on the quantitative null in every stratum and conservative on the binary null by discreteness; gate 5 removes exactly those rows; the 63 causal genes (not 21 -- the 2026-09-10 rebuild) are detected equally at both gates on burden and FlexRV and 16 -> 15 on the SKAT-O grid by CCT dilution within the bound. It found **LEDGER #128 (OPEN)**: a SKAT-O run dies when a region with two or more cells has every cell refused by a RAISED `--minGroupMAC_in_BurdenTest` (unreachable at the default 1); characterization test in place, fix candidates on the row, best decided together with #116, which it re-observed (an empty region is an `NA` Cauchy row on the burden path and nothing on the SKAT-O path). **RESOLVED 2026-09-12**, the assumption CONFIRMED and quantified (`docs/state/prediction-cct-dependence-2026-09-12.md`): the pooled Cauchy row read 0.06-0.08 at alpha 0.05 on the quantitative null in every arm while alpha 0.01 was on target, and it is the COMBINATION, not the inputs. CCT is EXACTLY uniform at rho 0 and at rho 1 (Cauchy stability, index 1), so the distortion is non-monotone with an interior peak; measured on SAIGE's own `get_CCT_pvalue` it peaks at rho ~0.5-0.7 and grows with k -- 0.0613 at k=18, 0.0698 at k=192 -- while alpha 0.01 stays within +0.002-0.005. The real null flex192 run at N=50,000 has measured inter-set correlation 0.705, that same peak band, and goes 0.0522 (sets) -> 0.0714 (pooled), which the simulation reproduces; the pooled row reconciles with the CCT of its own set rows to the file's printed precision. With n=210 genes the real arm's 95% CI is [0.0405, 0.1151] and cannot reject nominal by itself -- the precision is the simulation's. **So `p_FlexRV` has size ~0.07 at alpha 0.05 on this design and is accurate at discovery thresholds; never quote its nominal size at a large alpha.** An earlier reading of this session, that the cell-level 0.055-0.059 above showed the inputs already liberal and CCT possibly exonerated, is REFUTED on that arm (inputs at nominal, pooled not); where both are present (`wgridr`/`wplain`, k=9) they add, CCT's share about +0.009.

**Read the handoff, then this file, then
`python3 harness/tools/ledger.py open`.** Never write a hand-counted "N open of
M" here; the lint rejects it, because the count drifted every time.

**2026-09-12 (later): an adversarial review of the UNPUSHED lever diff found
two concurrency regressions, both now fixed and one DRILLED BOTH WAYS.**
`e72e4186` gives the two cell-pool prepares their own `hasVarRatio` (lever 2
had threaded them while they still shared the frame-scope out-flag that
`selectVarianceRatio` reads back as control flow -- the same bug TSAN forced
out of the marker loop on 2026-09-06); `8f4d5db7` adds `ScratchGuard`, so the
worker scratch is restored when a cell test leaves by THROWING, which the ER
deferral does by design and which the pool catches and continues past.
**TSAN is the evidence for the first**: `d36e7277` reports 5 races at t=4 and 3
at t=8 on `src/SAIGE_test.cpp:1274`/`:1279`, the flag's own two writes, and the
fixed tree reports none -- the real bug as its own plant. `bytes-vs-serial` was
`identical` in both, so no golden ever moved. A third commit declares the
region worker's FORK edge, which had left the clean arm at 1-2 warnings in
every build and is why a non-zero count there meant nothing. Of the review's three
remaining survivors, ONE was real and turned out to be **much worse than the
review rated it (medium/low): it produced WRONG NUMBERS on a reachable path**.
`CarrierStore::begin()` reset the slot map on a changed sample count but not
the VALUE buffers -- the changed-n branch clears `carriers`, so the
`used = carriers.size() * q` the else arm computes came back 0 and nothing was
zeroed, while `ur`/`urNw`/`sum` still held the previous call's values in
exactly the slots the next carriers were about to be handed, and both
accumulations assume a zero start (`u = max(u, ...)` for the collapse, `+=`
for the burden). MEASURED, not inferred: two `SPAGMMATtest` calls in ONE R
process, the full 1,000 samples then a 950 subset, on the pre-fix build --
both return normally and the SECOND run's output is wrong, `GENE1 lof` at
max_MAF 0.001 reporting MAC 7 and p 7.860290E-02 where the truth is MAC 2 and
p 3.182528E-01, and at 0.01 MAC 69 against 65. Inflated MACs are the
mechanism's own signature. CLAUDE.md's "one dataset per R process" put this at
the edge of the contract, but NOTHING REFUSES IT -- which is the worst of the
three possibilities, because a user looping over sample subsets gets silently
wrong answers. Fixed, and pinned by `tests/testthat/test-step2-region-tworun.R`
(drilled both ways: red on the pre-fix build at exactly the second-run
assertion, with the first-run control green). The other two are
DECIDED-NOT-TO-CHANGE, with the reasons recorded where a reader will hit them:
the BURDEN memo line's unit is cells and deliberately so (`R/BurdenMemo.R`;
re-basing it on tests would move a count `test-step2-flexrv.R` asserts), and
the three N-length `MarkerScratch` buffers the single-variant marker scan
never reads cost ~10 MB at N = 50,000 across 8 workers, which does not buy a
change to that path's hot scratch. Also: the claim that
`tsan.sh` does not cover the region cell pool, which the review made and which
propagated into my own summary, is FALSE -- the `region` workload was there all
along; it had merely never been run on this work.

**2026-09-12 also opened LEDGER #132 and #133 (both OPEN, both pre-existing),
from the adversarial review of the #129/#130 diff.** #132: ACAT-V's combination
weights are `c(w, rep(1, numofUR))`, so a uniform rescaling of a group file's
weights -- a no-op for every other p-value, measured byte-identical -- shifts
`Pvalue_ACATV` by up to 10.5 %, because the pseudo-marker's fixed 1 does not
scale and CCT normalises by `sum(w)`; the control is that with the collapse OFF
every column is identical. #133: a cell's `MAC` is folded while
`MAC_case`/`MAC_control` are not, so they do not sum to it above AF 0.5 (444
against 151 + 1405 in P62/P63), and the suite's existing identity assertion runs
only on a golden where no cell exceeds 0.5. The review also corrected four false
statements in that diff -- read section 8 of the handoff before writing anything
that cites it.

**2026-09-12 opened LEDGER #131 (OPEN, a DECISION)**, out of the same question
that settled #130: the MAF filter and the beta weight still FOLD to the minor
allele while the annotation they select on describes the ALT allele, so a site
whose ALT allele is common (AF 0.778) is admitted to a max-MAF cell as MAF 0.222
and handed `dbeta(0.222, 1, 25)` = 0.0606, the weight of a 22 %-frequency
variant, with `MAC` 444 printed rather than 1556. After #130 the collapse and the burden
describe one vector, still minor-coded per constituent, while the weight and
the admission use the folded frequency. It
is upstream's behaviour and the published SKAT convention, so it is a modelling
question and the row carries three candidate answers; it deliberately has NO
characterization test, because what one would pin is the choice. **WIDENED the
same day:** the GENOTYPE fold belongs to this row too and is its larger half --
a group test folds every marker above AF 0.5 to its minor allele before the
burden or the collapse is built, so #130 could not and did not make the burden
ALT-coded. Reachable on
bundled data only through P62/P63's raised collapse cutoff or
`--AlleleOrder=ref-first`.

**2026-09-12 opened AND CLOSED LEDGER #130.** Opened writing the flip drill's
characterization test (`tests/testthat/test-step2-region-flip.R`): on a
`--r.corr=1` run, an ultra-rare pseudo-marker whose collapsed alt frequency
exceeds 0.5 entered the BURDEN sum as `2 - G` over the UN-folded carriers with
every non-carrier's 2 dropped, so one collapsed marker was tested as two
different vectors and the emitted column was neither. **The owner's decision
was that the pseudo-marker is never folded at all**, and the reason is the one
the row had not considered: the group file's annotation is a property of the
ALT allele, so a folded burden is a burden of the allele the annotation does
NOT describe (and for a collapse, a burden of ABSENCE). The three-way
reconciliation it bought, the goldens that moved and the dead `urFill` it left
are on the row. Two things learned that outlive it: the SKAT-O arm DOES read
the burden column on a binary trait (`g.sum` in
`get_newPhi_scaleFactor_traitType`), which the row had denied; and the same
ALT-vs-minor question one level up -- the MAF filter and the beta weights still
fold, so a site with AF(ALT) 0.778 is admitted to a max-MAF cell as MAF 0.222
and handed a rare-variant weight -- is a statement about RARITY, is NOT fixed,
and is written up as its own row.

**2026-09-11 (evening) also closed #117 and #115**, and removed ~40% of the CPU
#118 added (`harness/bench/NOTES.md`). #117's fix -- a weight of exactly 0 drops
the marker -- means **binary FlexRV is no longer restricted to `--r.corr=1`**.
**#114 was DECIDED the same evening: RETAIN the preconditioner option**, with a
note that we may want to remove it. The unreachable arms stay, commented at the
three places a reader would tidy them away. It is a REMOVAL CANDIDATE, not a
resolved design: revisit if a dense-GRM fit at biobank N ever needs a better
preconditioner, or at the next dead-surface audit if it is still uncalled. One
row remains OPEN, #116 (a region analysed but contributing no row leaves a
zero-byte file; mechanism NOT diagnosed -- a NEIGHBOURING route was diagnosed
and fixed with #126 and is recorded on that row as explicitly NOT #116's).

**2026-09-11 closed #118, #126 and #127** (the per-set MAC, the min-group-MAC
gate on both region paths at a default of 1, and the pooled row's
`.pooledTests.txt` sidecar). Twelve group goldens were re-baselined, predicted
per file and per value beforehand in
`docs/state/prediction-118-126-127.md`.

**2026-09-09 closed #107-#113 and opened #115/#116.** What remains is one
owner's decision (#114), two one-line help fixes (#115), one row that needs
its mechanism diagnosed before it can be decided (#116), and #117 (a weight of
exactly 0 inside a binary SKAT-O cell returns NA), which is now the only place
where a zero weight is still not honest. The 50k
strip-vs-upstream benchmark was re-run with a Firth-off arm:
`harness/sim/upstream_compare/RESULTS-2026-09-09.md`.

Short list, kept deliberately short. `LEDGER.md` says what is broken;
`docs/state/*.md` says what was done and why; this file says what is NOT done
and is worth keeping in mind. It is not a plan and has no phases.

**Provenance is marked on every item, because it varies.** Two came from a
read-only assessment written by another session on 2026-09-03 and are recorded
here because they are expensive to rediscover -- neither is backed by a
measurement taken for it. Two came out of the S12a / LEDGER #50 work in the same
session. Nothing here has been profiled at realistic N, which is itself the
point: see "the gate" below.

## Resolved: the four known-wrong goldens WERE re-baselined (2026-09-04)

[Grounded; LEDGER #53/#56] Both step-0 sparse GRMs and both `binary_denseGRM`
step-1 fits were regenerated after the `geno2` fix and committed. Each was
proved byte-reproducible first, as CLAUDE.md requires: two independent
`make_goldens.sh` runs into separate candidate trees, `cmp` on all six files.
`compare_goldens.sh` is back to "all comparisons PASS" and the standard tier is
green again. Nothing else moved -- every step-2 golden and every sparse-GRM
step-1 golden was byte-identical throughout.

## The gate on all of it: profile at realistic N first -- NOW DONE

Every step-2 profile so far is N = 1,000, where `harness/bench/NOTES.md:352-359`
shows the serial writer -- not the work -- sets the ~2x ceiling at 4 threads.
That ranking is not expected to hold at biobank N, so choosing among the items
below before measuring is guessing. The sim fixtures are the instrument
(`harness/sim/`; the N=50k `bench50k` build lives outside the repo).

**MEASURED 2026-09-03, so this section is now evidence rather than a worry.**
See `harness/bench/NOTES.md` (the group-test section) and LEDGER #51. At
N = 50,000 over 210 genes: `--is_fastTest=TRUE` costs **+16.3%** and returns
BYTE-IDENTICAL output, because the region path never reads the flag. Escalation
rates run 2.9% (`y_null`, against a 5% nominal -- a clean negative control) to
14.8% on signal arms.

**An earlier draft of this file said an underpowered fixture would understate
group cost "by up to ~2x". That was wrong as stated** and is corrected here: the
2x applies per ESCALATED gene, but at a ~10% escalation rate the aggregate bias
is ~10%. Still worth reporting the rate alongside any per-gene cost -- it is
signal-dependent across ~5x -- but it is not the large effect claimed.

**The measurement moved the target, and then the profile moved it again.**
Swapping the PCG score test for `scoreTestFast` across a whole pass saves only
**8.3%**, so ~90% of gene-test time is NOT the per-marker score test. That was
written as an open question -- ultra-rare collapse, P1Mat/P2Mat assembly, or
R-side SKAT-O post-processing? **Answered 2026-09-04 by `Rprof` at N = 50,000**
(`harness/bench/NOTES.md`, "Where the OTHER 90% of a gene test goes"):

| bucket | share |
|---|---|
| C++ `mainRegionInCPP` | **43.9%** |
| R `region_postprocess` | **53.8%** |
| .. of which `SKAT_davies` under `integrate()` | **44.9%** |
| .. R glue | 5.0% |

So the largest single item in a SKAT-O group run is the SKAT package's Davies
routine under R adaptive quadrature, and **threading the C++ region loop caps at
1.49x on 4 threads** for SKAT-O. For BURDEN (`--r.corr=1`) the same fixture is
93.7% C++ and the same threading would give 3.36x. Replicates agree to 0.3 pp and
the BURDEN arm is the negative control -- it predicts `SKAT_davies` = 0 and
measures 0.0%.

**DONE (2026-09-04): the cells are evaluated in parallel.** `R/SkatoParallel.R`
drives `region_postprocess` three times -- collect the cells, evaluate them on a
persistent fork pool, serve the answers -- so its 175-line triple-nested loop did
not have to be restructured. `--nThreads` on a group run drives it (it used to be
discarded there). **1.30x at N = 50,000, byte-identical**, 8 workers no better
than 4. The remaining ceiling is 1.43x, so most of the available win is taken.

Still open on this path, in rough order of value. **Both now have a standalone
brief: `docs/state/region-perf-handoff.md`**, which carries the profile, the
Amdahl table, the two traps that each cost a measurement cycle, what was
measured and REJECTED so it is not retried, and the one design question that
must be answered before item 2 below is coded.
* per-task overhead still costs the gap between 1.30x and 1.43x; batching cells
  into fewer, larger tasks would narrow it and has not been tried. Small, low
  risk, and the gate already exists.
* [DECIDED 2026-09-05, LEDGER #51] `--is_fastTest` on group runs is pure cost
  (+16.3% at N = 50,000, byte-identical output) and stays so by decision: the
  region path never reads the flag, the R side still escalates, and a warning
  says so. The potential speed-up is bounded at the **8.3%** measured for
  `scoreTestFast` throughout, and taking it would move every fastTest region
  golden (P2, P3, P8, P8b, P10, P14, P15, P21, P23, P24, P34, P35) and retire
  P35's reason to exist. Only worth revisiting once the two items above are
  taken and the region loop is the bottleneck again.
* ~~the C++ region loop is the other 43.9%; threading it would compose with the
  above~~ **DONE 2026-09-06** (`3eeb554e`): `--nThreads` drives the C++ region
  marker loop, byte-identical at any count, and CLAUDE.md documents it. This
  bullet was left stale here until 2026-09-09; the remaining region-path item is
  the task-batching one above.

## From the 2026-09-04 upstream-issue-tracker mine (LEDGER #63-#85)

The 186-issue upstream tracker was mined for defects live on this branch; the
plan and the raw-evidence pointers are in
`docs/state/plan-tracker-mine.md` (**do not re-run the sweep** --
it cost ~160 agents). Phase 1 is DONE: **23 entries, LEDGER #63-#85**, one row
each, no fixes. Every line number was re-resolved against HEAD and every
arithmetic claim re-derived before the row was written, which is what the plan's
Tier 0 retraction was about.

**#63 is the worked example and is the only one with a test**:
`tests/testthat/test-characterization-ledger63.R` (standard tier), detection
power drilled BOTH ways -- 20/20 at HEAD, 8 failures with the `- 1` removed and
the package rebuilt into a scratch library.

**Reproduced by hand this session, not merely traced**: #63 (the dropped GRM
pair, with the committed `.mtx` reproduced exactly as the positive control),
#68 (`x1^2` silently dropped from a formula), #72 (the `std::stod` boundary is
DBL_MIN, i.e. every subnormal -- not the reported 1e-324), #73 (conditioning on
a monomorphic marker: exit 0, and the `_c` columns identical to the marginal
ones on 70/70 rows, against 67/70 changed for a real conditioning marker), #74
(the data.table join-type message reproduced verbatim), #76 (a region run at
`--chrom=2` writes a ZERO-BYTE results file at exit 0 and says "Analysis
done!"), #80/#81/#82/#84 (verified against the committed goldens and `git
blame`). The rest are traced-at-HEAD with the consequence explicitly marked as
not reproduced -- read each row's PROVENANCE, it says which.

**Four candidates were REFUTED during this pass and deliberately have no row.**
Recorded here so nobody re-opens them:
- **empty ultra-rare collapse** (upstream #12/#30/#60/#67) -- the C++ handles it
  coherently, checked at four points: the pseudo-marker slot keeps its NA/NaN
  initialisation, `annoMAFIndicatorMat.row(i)` stays all-zero so the cell is
  excluded from the region test, `i1`/`i1InChunk` are NOT incremented so VarMat
  stays square and aligned with `mPassCVVec`, and `writeOutfile_singleInGroup`
  skips any row whose p is `"NA"`. If a defect lives here it is in the
  `noNAIndices`-vs-`VarMat` alignment in `region_postprocess`, not in the
  collapse.
- **sparse-GRM nnz overflow** (upstream #188) -- the proposed half-to-full
  storage doubling does not happen: in the pinned env (Matrix 1.7.6) a
  `dsCMatrix` survives `*`, `diag<-` and `[i, i]` with half storage intact
  (checked). What remains is only that `Matrix:::readMM` is called unguarded at
  `R/SAIGE_fitGLMM_fast.R:2109` and `:2160`, with no mechanism located.
- **SKAT-O compact-support p = 1** (upstream #209) -- the outer integral is
  inside `SKAT:::Met_SKAT_Get_Pvalue`, a dependency, not in this tree.
- **`null_noXadj` cardinality check in the wrong arm** -- this is already
  **LEDGER #39**, which names the unreachable `stop()` explicitly. #39's line
  numbers are stale by ~17 lines (`:405-414` -> `:388-399`).

**Phase 2 is COMPLETE; Phase 3 (fixes) has STARTED** -- #74 and #66 are done
(`22be2a69`, `c9f7e9ae`), each with its characterization test INVERTED in the
same commit. The plan orders Phase 2 by gate
gap rather than severity. Item 1 is DONE: goldens **P36/P37**
(`plink_group_burden_quant.txt`, `plink_group_nofast_quant.txt`) close the
QUANTITATIVE region-path gap -- every previous group golden is binary -- and
`tests/testthat/test-step2-group-quant.R` drives them. Both were proved
byte-reproducible over three independent runs before being accepted.

**That golden immediately paid for itself: LEDGER #86.** P37 reports SKAT-O
`Pvalue` **exactly 1** on 5 of its 12 cells while their own `Pvalue_Burden` and
`Pvalue_SKAT` are 0.82-0.99, which the suite's existing invariant
(`test-invariants-goldens.R:106-110`, SKAT-O lies between its component tests)
forbids -- it was asserted on one BINARY golden only. Measured on the kernel
directly: `SKAT:::Met_SKAT_Get_Pvalue` returns exactly 1 for 55 of 400 random
inputs, all with `min(p.val.each) > 0.777`, and never below 0.763. And it is not
inert -- both CCT implementations fall back to `min(1, n*min(p))` when any input
is 1, so one saturated cell silently turns the region's `Cauchy` row into a
Bonferroni (confirmed to the last digit on GENE2). **This corrects the
refutation recorded above**: the sweep's SKAT-O finding was dismissed as living
in the `SKAT` dependency and therefore unlocatable. The integral does live
there, but the consequence is here and is now pinned.

**Item 2 is DONE**: `tests/testthat/test-step2-cond-highld.R` builds a
conditioning fixture AT RUN TIME (seeded, no golden -- the genotypes come from
R's RNG, so a golden would pin the RNG rather than the code) over the same
10,000 fam rows as `genotype_100markers`, so the committed
`binary_sparseGRM_a` model is reused unchanged. It carries an LD ladder against
one conditioning marker (r^2 = 1.000 / 0.90 / 0.49 / 0.000, MEASURED on the
1,000 analysed samples inside the test) plus a marker built to correlate 0.35
with the null model's own residual, so a large conditional statistic survives
conditioning.

Both target regimes are reached, and both controls behave: at r^2 ~ 0
conditioning is a no-op (var/var_c = 1.001), at r^2 ~ 0.5 the variance really
drops (1.89, against 1/(1-r^2) = 2.0).
- **#84**: the associated marker gives `stat_c = 91` against a gate of 4 with
  `var/var_c = 1.51`, and `p.value_c` (2.8e-18) differs from `p.value.NA_c`
  (1.3e-21) by three decades -- the conditional SPA block executes at real
  conditioning for the first time on this branch. DRILL: on golden P12, **zero**
  of 70 rows satisfy both `stat_c > 4` and `var/var_c > 1.3` (its three
  high-ratio rows are the conditioning markers themselves, `stat_c` ~ 0.003).
- **#66 -- FIXED in Phase 3 (2026-09-05, `c9f7e9ae`); what follows is the
  measurement that located it.** Under `--is_noadjCov=FALSE` a perfectly
  collinear marker gives
  `var_c = 2.5e-14` -- cancelled, and yet ~294 orders ABOVE the `DBL_MIN` guard,
  which is exactly the geometry LEDGER #44 had to fix with a relative test.
  `BETA_c` and `SE_c` are then ratios of two cancelled quantities (SE_c = 6.3e6)
  and `p.value_c` is forced to 1. DRILL: adding a relative cancellation guard to
  `src/SAIGE_test.cpp` and rebuilding into a scratch library fails 6 assertions;
  source reverted and reinstalled.

**SETTLED 2026-09-05 (recorded in LEDGER #82).** At the CLI default
(`--is_noadjCov=TRUE`) the same perfectly collinear marker reports a non-zero
`var_c` because the tested marker's variance is the UNADJUSTED, allele-frequency-
centred kernel while the conditioning cross-term is built from the covariate-
ADJUSTED marker. The residual is exactly `vr_noXadj * [centring-mismatch term +
W-weighted covariate term]`: predicted from the model to every printed digit
(0.362754 on the run-time fixture) and, from the two runs' printed columns alone,
`var(def) - (vr_noXadj/vr_null) * var(adj)` to printing precision. The earlier
`var * R^2_OLS` prediction was low because it dropped the centring term (30% of
the residual) and used an unweighted R^2. Pinned in `test-step2-cond-highld.R`.

**Item 4 is DONE**: `tests/testthat/test-step1-varratio-cursor.R`. The
variance-ratio marker loop advances `indexInMarkerList` once unconditionally and
a SECOND time on each skip branch, so a skipped marker consumes the marker after
it, which is never read at all -- and because the exhaustion guard tests
EQUALITY, a double step from the last position jumps over it and the next read is
an NA.

Reproduced by hand, and the arithmetic reconciles two ways. With a run-time
fixture of M markers at MAC 40 (all in the single VR pool) over the step-1
fixture's own 1,000 fam rows:
- NEGATIVE CONTROL, all autosomal: exit 0, exactly 30 markers read (the
  `--numRandomMarkerforVarianceRatio` default), a valid ratio file.
- **All non-autosomal, EVEN pool (M = 40)**: exactly **20 of 40** markers read --
  the cursor visits 1, 3, 5, ..., which is the arithmetic signature; a
  single-increment loop reads all 40. Lands exactly on the guard, exits the loop
  cleanly, then dies downstream in `calCV(NULL)` with no marker accepted.
- **All non-autosomal, ODD pool (M = 41)**: 22 reads, and the failure MODE
  changes -- `Error in if (genoInd == 0) { : missing value where TRUE/FALSE
  needed`, which is #75's predicted abort and upstream #32/#151 verbatim. **One
  marker of difference flips the outcome**, which no single-increment loop can do.
- Both leave a 139 KB `.rda` beside a **0-byte** `.varianceRatio.txt` -- the
  unusable model pair #75 describes.
- The SILENT case is the partial one: 20 of 40 non-autosomal gives **exit 0**
  with the ratio estimated from 12 markers instead of 30, and 14 pool markers
  never examined. Reconciled: 12 accepted + 14 skipped = 26 reads, cursor
  = 1 + 12 + 2*14 = 41, exactly the guard value for a pool of 40.

DRILL: removing the second increment and rebuilding into a scratch library fails
5 assertions (read counts become 40 and 41, and the NA abort disappears); the
negative control correctly stays green. Source reverted and reinstalled.

**Item 5 is DONE**, and **#74 is now FIXED** (Phase 3, `22be2a69`; the test
below was inverted with it, and a pgen arm added because the pvar half of the
fix was otherwise ungated). What follows is the characterization as first
written. The bim WAS
read with no `colClasses`, so data.table INFERS the ID column's type; every
bundled fixture uses `rs<N>` (character), and an all-numeric ID column infers
INTEGER. The fixture is derived at run time and reuses the committed `.bed` and
`.fam` UNCHANGED -- a `.bed` is indexed by row, so rewriting bim column 2 changes
no genotype, which is what makes the control exact.

- NEGATIVE CONTROL: a plain marker scan is UNAFFECTED -- same row count, and
  `p.value` identical string-for-string between the two bims. The ID column is
  not joined on that path.
- `--condition=13` -> exit 1, `Incompatible join types: x.ID (integer) and
  i.SNP (character)`, upstream #137's message verbatim.
- a group file with numeric IDs -> exit 1, the same `bmerge` failure.
- POSITIVE CONTROL: both commands with the committed rs-spelled bim exit 0.

DRILL: adding `colClasses = list(character = 2)` to the `fread` and rebuilding
into a scratch library makes both numeric-ID runs SUCCEED, failing 5 assertions;
the plain-scan control stays green. That one argument is also the Phase 3 fix.
Source reverted and reinstalled. **Phase 3 then took exactly that argument**,
plus `list(character = 3)` on the pvar, and strengthened the positive control:
the numeric-ID and rs-spelled runs must now agree on EVERY column but
`MarkerID`, which they do.

**Item 6 is DONE**: `tests/testthat/test-step1-duplicate-pheno-id.R`. Nothing in
the step-1 phenotype path checks for duplicate sample IDs -- the only
`duplicated()` call in `R/SAIGE_fitGLMM_fast.R` is on `SampleIDIncludeFile` -- so
a repeated ID multiplies that sample's rows through the merge at `:1012-1014`.

A phenotype file of 1,001 rows with 1,000 distinct IIDs runs to **exit 0** and
reports "1001 samples will be used for analysis"; the unmodified file reports
1000 (negative control).

The consequence was DERIVED FROM THE CODE BEFORE IT WAS MEASURED, which is what
makes it a test rather than a fishing expedition. In `setGenoObjFromSource`,
`Nnomissing = subSampleInGeno.size()` counts the DUPLICATED list, while the
allele-count accumulation is gated on `indicatorGenoSamplesWithPheno_in[ind]`
with `ind` running over the FILE's samples, so each distinct sample is counted at
most once. `altFreq = alleleCount / float(Nnomissing * 2)`
(`src/SAIGE_fitGLMM_fast.cpp:677`) therefore puts the duplicated sample in the
DENOMINATOR but not the NUMERATOR. Predicted gap: exactly that sample's genotype.
Measured, on a fixture where the duplicated sample carries the ALT allele in half
the markers: store holds 41 alleles and the frequency implies 40 where it is a
carrier, 39 vs 39 where it is not -- `store_AC - reported_AC == g_dup` on all 40
markers to float32 precision. Both arms in one run, so the probe cannot be
fabricating the signal.

DRILL, two patches in one scratch build because the two tests pin different code:
deduplicating the merge in R breaks the sample-count assertion, and changing the
`altFreq` divisor to the indicator count breaks the three frequency-gap
assertions. 4 failures, exactly the intended ones. Both reverted and reinstalled.

**Phase 2 is COMPLETE.** See the handoff for the gate numbers.

## From the 2026-09-04 silent-failure audit

A 26-agent fan-out hunted the bug class this branch keeps finding by accident
(silent failure that every golden reproduces). 15 raw findings -> 11 unique ->
10 adversarially verified by two skeptics each, both instructed to default to
REFUTED. Five survived; five were killed, including a plausible `raw_variant_ct`
claim and a `dosage_zerod` one.

**Two were reproduced by hand the same day and are FIXED (LEDGER #58, #59)** -- the
unbounded mmap `.bed` read (a `.bim` claiming one more marker than the `.bed`
holds gave byte-identical output at exit 0, while the pread path stopped), and a
region results file written with data rows and no header. Both have positive and
negative controls in the ledger entries, and both are gated by
`tests/testthat/test-step2-file-integrity.R`.

**The other three were reproduced on 2026-09-04 and are now FIXED
(LEDGER #60, #61, #62).** None was dropped; all three were real, and two were
worse than the finding had claimed. Each carries a positive and a negative
control in its ledger entry:
* **#60**, `src/SAIGE_fitGLMM_fast.cpp` (`BedByteSource::fill`) -- the step-1
  twin of #58. `.bed` truncated to 60 of 100 markers, `.bim` left at 100: exit
  0, no warning, GRM from 18 QC-passing markers instead of 33. Truncating
  MID-record is worse than the finding predicted -- half a real marker plus
  fabricated homozygotes PASSED QC and entered the GRM, a marker the intact file
  rejects. Fixed with #58's size check plus a `gcount` check in `fill()`.
* **#61**, `R/SAIGE_Test_main.R:602` -- and the effect is not "all its markers
  discarded" as the finding read it, it is that the bin comes out a
  byte-for-byte DUPLICATE of the bin below, which looks like an answer.
  `--maxMAC_in_groupTest=100` on the 1,000-sample fixture (implied MAF 0.05
  against a largest cutoff of 0.01) gives `Number_rare` 3/4/4/1/3/3 where the
  same bin asked as a MAF gives 4/7/7/2/7/7, and moves a p-value 3x. Fixed by
  re-issuing `setRegion_GlobalVarsInCPP` with the extended vector.
* **#62**, `R/Util.R` (`updateChrStartEndIndexVec`) -- confirmed exactly as
  described: one chr1 marker inside chr2 passes both arms of the guard and
  chr1's window then spans 16 markers for its 11, dropping five chr2 markers
  from both LOCO GRMs. Fixed by asserting contiguity, which is the invariant the
  windows actually need.

Detection power drill-verified in both directions for all three: the new tests
are red at unpatched HEAD (1 / 5 / 4 failures) and green after.

## From the strip-vs-upstream comparison at N = 50,000 (2026-09-05)

[Measured; `harness/sim/upstream_compare/RESULTS.md`] Upstream `main` (365deeff)
built beside `strip` in its own pixi env pinned to the same R 4.4.3 / OpenBLAS
0.3.30 builds, both arms on `sim50k`, paired and alternated.

- **LEDGER #97 -- FIXED 2026-09-05** (`R/SAIGE_fitGLMM_fast.R`, the `sparse` rows
  of `varRatioTable` no longer multiplied by `varRatioScale`). The #65 rescale had
  multiplied the scale-free `sparse` variance-ratio rows by 1/s^2, so every
  quantitative model fitted with `--useSparseGRMforVarRatio=TRUE` between
  `66a65bb9` and the fix carried `sparse` rows 11-12% high on the N = 50,000
  fixture (1/s^2 ~ 1.5 on the bundled one); refit those models. Gate:
  `tests/testthat/test-step1-vr-sparse-scale.R` (standard tier) -- a quantitative
  step-1 with `--useSparseGRMforVarRatio=TRUE` whose `sparse` rows are invariant
  to y * 100 while its `null` rows scale by 1e-4, and whose `sparse` rows are ~1.
  No committed golden moved (none has quantitative `sparse` rows).
- **RESOLVED 2026-09-05 (evening): the binary single-variant scan's 0.94x was not
  like-for-like, and its kernels were profiled and de-allocated.**
  (`docs/state/handoff-2026-09-05-binary-profile.md`; `harness/bench/NOTES.md`.)
  Firth is 56% of the binary scan's CPU on both arms; upstream's #28 made it skip
  ~1,000 Firth fits and ~450 SPA corrections on AF > 0.5. A per-worker
  `KernelWorkspace` (byte-identical at 1 and 4 threads, standard tier green) took
  the 1-thread user CPU to 205 s (upstream 236) and the 4-thread sys time from
  55 s to 11 s. `--firth_solver=closed` (-27% more, output byte-identical, P31 unmoved) is
  now the DEFAULT with `qr` selectable; the warm start (118 non-converged fits)
  was benchmarked and deleted. The PLINK reader now
  decides QC gate 1 from 2-bit code tallies before expanding a marker
  (quantitative 50k scan 62 -> 45 s, byte-identical).
- Strip's step-2 peak RSS at N = 50,000 reads 2.3-2.8 GB against upstream's
  0.35-0.38 GB. Measured to be resident pages of the mmap'd 2.6 GB `.bed`, not
  allocations: `SAIGE_PLINK_NO_MMAP=1` gives 312 MB with byte-identical output.
  Worth a sentence in the user docs, not a fix.

## DECISIONS for the owner (2026-09-05): all twelve answered

The eleven rows the Phase 3 sweep left open, plus #18's second half, were
each decided and implemented the same day
(`docs/state/handoff-2026-09-05-decisions.md`). `python3 harness/tools/ledger.py
open` lists what is still open (documentation rot, #92/#93/#95/#96, and the
Sigma floors #91); `ledger.py decided` the four left-as-is rows.

## Time the gate at the next run (2026-09-05)

The standard tier took ~4 min 20 s per compare mode on the last pre-push
(`harness/work/prepush-{byte,numeric}.log` mtimes 21:28:34 and 21:32:54, on a
machine still at load 10-14 from the N = 50,000 gate), so `prepush.sh` is ~9
min against the ~5 CLAUDE.md implies. Not yet measured idle. What this session
added to the standard tier, all step-1 fits (wall time is bounded by the
slowest FILE): `test-step1-nglmm-prep.R` (two extractNglmm runs, ~20 s each),
`test-step1-lowmem-loco.R` (four 3 s fits), `test-step1-stale-model-file.R`
(three fits), `test-step1-vr-sparse-scale.R` (two ~5 s fits), and the
morning's `test-step1-pheno-scale.R` (three fits + three step-2 runs).
At the next gate: check `uptime` first, run the tier once with
`reporter = "progress"` and collect the per-file `[Ns]` markers (the parallel
result object carries no timings), then split or demote the long poles
(`nglmm-prep` is the likely one; its second arm could share the first's run).
Also worth asking whether numeric mode needs to re-run every step-1 fit, since
the two modes differ only in how step-2 tables are compared.

## From the 2026-09-05 code-quality review of the Phase 3 sweep

[`docs/state/handoff-2026-09-05-review.md`; measured on the two-chromosome
dense-GRM fixture that `tests/testthat/helper-pheno-scale.R` derives at run time]
The review's refactors landed (LEDGER #97-#100 fixed on the way; the fold of
`fitNULLGLMM`'s binary/quantitative halves proved byte-identical by `cmp` of
every model file against the pre-fold build). Three things it made visible and
deliberately did NOT change:

- **Binary LOCO refits are start-dependent, and the two layouts start
  differently.** `glmmkin.ai_PCG_Rcpp_Binary`'s LOCO loop passes `alpha`/`eta`
  that it overwrites with each chromosome's result, so chromosome j warm-starts
  from chromosome j-1's refit; `fitNULLGLMM`'s `--isLowMemLOCO` loop starts every
  chromosome from the full model's `alpha0`/`eta0`. At `tol.coef = 0.02` they stop
  at different points: chr1 agrees to 1e-10, chr2 to ~1e-5 (coefficients 8.5e-6,
  fitted values 2.4e-5, cov 7.7e-6) and the LOCO entry's `obj.noK$S_a`, the score
  at the refit, not at all. The quantitative refit is one linear step and agrees
  to 1e-10 on both. Pinned as found in `test-step1-lowmem-loco.R`; which start
  upstream intends is not known, and either choice moves every binary LOCO model.
  [ASSUMPTION: the warm start is the whole cause -- derived from the code, and the
  chr1-exact/chr2-loose pattern matches it; not isolated by forcing the start.]
  **RESOLVED by LEDGER #125 (2026-09-10, late night):** with the inner pass at
  `tolCoef = 1e-5` the two layouts reach the same point from their different
  starts (chr2 coefficients identical, `S_a` within 5e-4); the test now asserts
  agreement. The start-dependence was the loose inner tolerance, not the start.
- **Two by-trait asymmetries in the low-memory refit, preserved.** The binary
  refit passes the glm's offset (none unless `isCovariateOffset`) and
  `tol.coef = tolCoef` (was `tol`, LEDGER #125); the quantitative refit passes the model's covariate offset
  (`X[, -1] %*% coef[-1]`, alongside the same X and alpha) and the default
  `tol.coef = 0.1`. Both mirror what each trait's fitter does in its own LOCO
  loop, so the layouts agree per trait; whether the quantitative offset is
  double-counted is a derivation nobody has done. Marked "by trait" in the code.
- **The N = 50,000 own-model quantitative tables are a mixture and want one
  re-run.** The review's gates refitted the strip own models (sparse rows now
  1.00, the #97 fix) but did not re-run own-model step 2, so
  `harness/sim/upstream_compare/RESULTS.md`'s own-model `y_quant` section
  describes outputs from the pre-#97 models. From
  `harness/sim/upstream_compare/`: `MODELS=own bash run_arm.sh strip s2single
  y_quant 1`, the same for `s2group`, then `bash compare_all.sh 1 1`; ~10 min.
  Expected: the own `gene y_quant` tags lose their `VR`/`VR-mix` mass.
- **`#65`'s save-then-rescale is a design the owner may still change** (review
  finding 1, second half): the run saves the standardised quantitative model,
  reloads it for the ratio estimation, and rescales the files it wrote at the
  end (`savedModelFiles`, #98). The alternative -- save in phenotype units at the
  three save sites, standardise on the reload -- deletes `rescaleSavedQuantModels`
  and the CLAUDE.md sentence about standardised files, at the cost that the ratio
  estimation reads a `x s x 1/s` round-tripped model (ULP level), so P4, P16, P32,
  P36, P37 (+`singleAssoc`) and the quantitative step-1 goldens would need a
  second rounding-level re-baseline. Not done; owner's call.

## Resolved DECISION: SKAT kernel failures are counted and summarised

[LEDGER #57, the open half; decided 2026-09-05] `get_SKAT_pvalue`'s
`try(..., silent = TRUE)` still turns a kernel error into NA rows, by choice: the
user picked count-and-summarise over warn-and-NA and abort, so that what drives
the NAs can be determined. The run now prints `SKAT failures: N ...` with a
by-message tally before `Analysis done!` and writes
`<SAIGEOutputFile>.skatFailures.txt` naming each failed cell when N > 0. Result
files are byte-identical. Tests: `test-unit-skat-failures.R`,
`test-step2-skat-failures.R` (planted kernel failure in a child process).

## Items

**Resolved: formatting step-2 output rows on the workers is CLOSED by
measurement (2026-09-06).** [Measured; `harness/bench/NOTES.md`, the writer
section; `docs/state/handoff-2026-09-06-writer-bench.md`] The item predicted
that P rows per marker on the one writer thread cap `--nThreads` as P grows.
At N = 50,000, HEAD `b11e6b07`, the writer is 2-3 us per row and takes
w = 0.2-0.9% of the main thread's `mainMarkerInCPP` time at 1 thread and
0.5-2.1% at 8, on 1, 3 and 5 models alike (worst case 3 quantitative models,
w = 0.021, ceiling 1.02x). The thread-scaling loss the Amdahl fit called a
serial fraction is CPU inflation (1.3-2.7x more user CPU for the same work at
8 threads: allocation churn and hyperthreading on 4 physical cores); the main
thread is 91-95% busy as a worker at T=8. Do not re-open without a fixture
where marker work is under ~20 us (N well below 1,000), which is not the
target regime. Driver: `harness/bench_writer.sh`; reader:
`harness/bench/writer_profile.py`; table: `harness/bench/writer_analyse.py`.

**If phenotype batching is ever attempted, one trap is already known.**
[From the same assessment; recorded because it is cheap to write down and
expensive to hit.] Computing all P traits' score statistics as a single
`G' x [res_1..res_P]` matrix product CHANGES THE REDUCTION ORDER and therefore
breaks the byte gate; the per-trait dot products have to keep their serial
order. This is the same tension as the step-1 GetTrace GEMM decision, which was
resolved the same way. Nothing else about batching is designed, and the natural
unit for it -- a per-trait `TestContext` -- only became possible on the region
path with S12.

**Resolved: `--is_fastTest` warns on group runs, and
`--pval_cutoff_for_fastTest` is a CLI option.** [LEDGER #51 / S12a; done
2026-09-04] The warning sits beside the existing Firth/r.corr one, where
`is_fastTest` is final -- after the conditional path has forced it FALSE. #51
itself stays OPEN: the asymmetry is upstream's and unchanged, only the silence
is fixed. Exposing the cutoff made golden **P35** possible (P2's command with
the cutoff at 1, so both regions take the ESCALATION pass) -- **its bytes are
identical to P2's and that is #51 restated**, not a defect: the second pass
recomputes what the first already had. It pins the other WRITER, and under
LEDGER #50 it lost all 19 single-variant rows while P2 stayed green.

**Resolved: the psam sample map keyed on FID, not IID.** [LEDGER #52, FIXED
2026-09-04] `readPsamFile` now finds the `IID` column by name from the header,
strips `\r`, and refuses a duplicate IID instead of overwriting it. All 13 PGEN
goldens stayed byte-identical, exactly as predicted -- which is why the gate had
to be a fixture: `tests/testthat/test-step2-pgen-psam.R` derives psams at run
time (one FID shared by all 10,000 samples; `#IID SEX` with no FID at all; a
duplicated IID) and drives golden P13 through each. 6 failures at unpatched
HEAD, 0 after.

**PCG stops on an absolute residual, and that is a recorded hazard, not a
fix.** [DECIDED 2026-09-05, LEDGER #6] All eight PCG loops in
`src/SAIGE_fitGLMM_fast.cpp` stop on `sum(r^2) > tolPCG` with no right-hand-side
normalisation, so `--tolPCG=1e-5` is a different relative tolerance at every N
and for every phenotype scale. On the bundled fixtures it only constrains the
harness's step-1 tolerances (theta/VR 1e-2). At biobank N with a DENSE GRM the
stop is relatively looser than at N = 1,000, which is where it could bite; the
sparse-GRM fits never enter a PCG loop at the default flags (direct `spsolve`),
so a normalised rule would move only `binary_denseGRM` and its pgen mirror
[ASSUMPTION: traced, not instrumented]. If revisited: perturb `tolPCG` first to
see which goldens actually move, then normalise by `dot(b, b)` and re-baseline.

**Dosage pgens are hard-called for the GRM, and that is a real limitation, not
just a warning.** [Grounded, 2026-09-03] `genoClass`'s store is 2 bits per
sample, so `PgenByteSource` rounds any dosage to the nearest of {0,1,2} and
warns once with the count. For imputed data that discards information the GRM
and the variance ratios could have used. Carrying dosages properly means a
second storage path through `genoClass` and touches the GRM kernels, LOCO and
every step-1 golden -- deliberately out of scope, and worth measuring (does the
GRM actually move?) before being attempted.

## REFUTED: port the repo to a faster language (2026-09-06)

[Grounded; `harness/bench/NOTES.md` profiles of 2026-09-04/05 and the
2026-09-06 writer-bench handoff. Owner decision: no port, recorded here so it
is not re-opened.]

The owner asked whether the clean branch should now be rewritten in a faster
language. The measurements already in NOTES say the CPU is in compiled code:
- Binary single-variant scan, N = 50k, 1 thread: Firth 56.1%, `getOneMarker`
  15.1%, SPA 10.2%, impute 5.6%, score tests 7.6%, R writer 0.09%. All C++.
- SKAT-O gene tests: the 45% under `region_postprocess` is `.C("qfc")` (93.9 s
  self against 3.9 s of R glue), so R is ~1.8% of the run; BURDEN is 93.7% C++.
- The N = 1,000 benches are ~90% process startup, the one term a compiled
  binary would remove, and it is gone at 50k.
- Thread-scaling loss at T = 8 is the turbo/HT floor, not a serial R section.
- `src/` is ~20k lines of C++; `R/` is ~10k lines of CLI, orchestration,
  step-1 glue and region post-processing. The split a port would produce
  already exists.

What a port would cost: the byte-exact goldens and every ledger row of
bug-for-bug behaviour would have to be re-argued under a different BLAS, libm
and reduction order, and the SKAT-package parity on the region goldens would
have no reference. That is the branch's whole verification method.

The two levers that remain if speed is wanted, neither scheduled: the Firth
fit (33 ms per eligible marker, 9 Newton iterations each -- warm start, QR
reuse or a cheaper convergence rule inside the existing C++), and threading
SKAT-O cells in-process, which needs a C++ port of ONE routine
(`SKAT_Optimal_PValue_Davies` + `qfc`, gated by bit-identity against the SKAT
package), not of the repo. [Qualified the same day by the `qfc` audit in
`docs/state/region-perf-handoff.md`: SKAT's `qfc` cannot be LINKED on threads
-- file-scope statics, one global `jmp_buf` -- so that port would be a
reentrant reimplementation of Davies, not a copy. Not planned.]

## The region path (2026-09-06): measured, de-churned, split, threaded

[Measured; `docs/state/handoff-2026-09-06-region-threads.md`, the last four
sections of `harness/bench/NOTES.md`; plan
`~/.claude/plans/can-you-summarise-what-shimmering-snail.md`.] Four commits
after `7a9efcb8`: `SAIGE_PROFILE_REGION`; P1/P2 chunks in memory (42% of the
region C++ had been churn and a disk round trip, the VarMat dgemm is 1.6%);
pass A split into prepare-and-test / serial accumulate; `--nThreads` threads
the region marker loop, byte-identical at any count, TSAN clean + plant. At 4
threads on `sim50k`: SKAT-O 1.69x / 1.93x wall, BURDEN 1.57x / 1.64x. The
`qfc` audit closed in-process SKAT-O for good (see the port section below).

Decided order for what remained. **Steps 3b, 4 and 5 all landed on
2026-09-06**; this list was left stale here until 2026-09-09.
* ~~**Step 3b, next:** thread the per-cell tests of the ultra-rare pass and the
  BURDEN emitter (15 independent cells each; they are ~45% of a BURDEN run
  with accumulate and R, which is exactly why BURDEN reads 1.57x and not the
  3.36x once projected). Predicted ~2.2x at 4 threads on this machine.~~
  **DONE** (`e7b346ba`, `docs/state/handoff-2026-09-06-region-cells.md`):
  each pass is prepare / test / write-back, only the tests on the pool, so
  the prepare and write-back stay serial in cell order. Measured **~1.84x /
  ~1.71x** at 4 threads by the paired ratios, not the 2.2x predicted
  (`harness/bench/NOTES.md`, the step-3b section).
* ~~**Step 4, then:** several null models per region (`--groupFile` with a
  manifest). Step 1's buckets size it: the shared prepare is 18-23% of the
  SKAT-O C++ (~10% of a run per extra model) and 35-59% of BURDEN's.~~
  **DONE** (`1e02a154` C++ + `2ed00165` R,
  `docs/state/handoff-2026-09-06-region-multimodel.md`): every model in a
  region manifest must share one sample set, checked before any state reaches
  C++.
* ~~Step 5 (batch the SKAT-O cells) is independent and small.~~ **DONE**
  (`1c706bd2`, `docs/state/handoff-2026-09-06-region-skato-batch.md`): the
  cells of a whole batch of regions go to the pool in one call. 2026-09-09
  added identical-cell dedupe on top (`1b5bf4c1`).
* **Still open, unscheduled:** if BURDEN at scale matters, thread across
  REGIONS rather than markers, which approaches the hardware ceiling
  (~2.9x at 4 threads here).

## FlexRV: F0-F4 LANDED, F5 (the RAM/CPU measurement) is what is left (upd. 2026-09-09 late)

[Grounded; `docs/state/handoff-2026-09-09-flexrv-e2e.md`, before it
`handoff-2026-09-09-flexrv-f2-f4.md`.] The feature runs end to end:
`--flexRV_score` expands the 192 sets in C++, four goldens (P58-P61) gate it,
and every pre-existing golden stayed byte-identical in both compare modes.
**The end-to-end suite landed 2026-09-09 (late)**: all 192 sets reconciled
against 192 external weight lines (binary collapse-on and quantitative
collapse-off, two MAC sources), the MAC ladder against `--maxMAC_in_groupTest`
cells, `--condition` internal-vs-external on the conditional columns, the
`--flexRV_lofAnno` positive path, chunking and SKAT-O threads byte-inert, and
two null models per pass (`test-step2-flexrv-{reconcile,cond,paths,multimodel}.R`).
It opened **LEDGER #118**: a set's burden is tested with the CELL's unweighted
MAC (VR band, ER, the min-group-MAC gate, the printed `MAC`), upstream's and
the engine's semantics, the owner's decision. It also found and fixed
**LEDGER #119** (a one-set FlexRV spec stopped every SKAT-O run at K = 1;
the beta(1,25) set never reaches a transform set, `test-step2-flexrv-k1.R`).
What is NOT done, in order of value:
- **The five RAM/CPU levers LANDED 2026-09-12** (`bb4003a0`, `a1f44cdd`,
  `d73a9981`, `6ffa219a`, `4535ee99`; `docs/state/handoff-2026-09-12-flexrv-ram-levers.md`,
  predictions and outcomes in `docs/state/prediction-flexrv-ram-levers-2026-09-12.md`,
  NOTES.md "FlexRV RAM/CPU levers"): carrier-space storage for the three
  per-(sample x cell) matrices, per-set prepares on the cell pool, no
  N-length cell records, the reductions over the carriers, and identical
  BURDEN cells evaluated once. sim50k k192 BURDEN: RSS 874 -> ~520 MB, wall
  57 -> 19 s at 1 thread, 41 -> 13 s at 4 (loaded box), byte-identical at
  both thread counts and every golden unchanged. Of the four open items,
  TWO ARE DONE (2026-09-12): `harness/prepush.sh` is green at `2f1c1960` in
  BOTH compare modes (4428 / 4628) and again at `17676759` with the new test
  (4497 / 4697); and the flip drill is now
  `tests/testthat/test-step2-region-flip.R`, which found LEDGER #130 (above).
  **The 4-thread timing is DONE and the loaded-box numbers are CONFIRMED**
  (NOTES.md, "RE-TIMED on an idle box"): 20.3-21.3 s at 1 thread and
  13.2-14.7 s at 4 over three alternating reps, all six outputs
  byte-identical, RSS 519-522 / 538-542 MB -- within 5% of the 19.4 / 12.6 s
  already recorded, so **1.55x at 4 threads is the number to quote** and
  re-timing will not improve it. Two things to carry forward: **on this box
  load average is NOT a quiet-box criterion** (it sits at ~3 with `top`
  reporting 91% idle and 0.31 of 8 cores busy, because it counts blocked
  threads -- use CPU idle instead); and the scaling loss is CPU INFLATION,
  1.75x the user time for 1.55x the wall, so the next lever has to remove
  work rather than contention. Still open from the handoff: the binary arm
  at 50k. The F5 bullet below is history.
- **F5 timing DONE 2026-09-10** (`harness/bench/NOTES.md`, "FlexRV F5 at
  N = 50,000"): 192 sets cost 35 s user / 805 MB against 8 s / 506 MB for one
  set on sim50k `y_null` BURDEN at one cell; the K-dependent C++ is 39%
  ultra-rare pass, 32% per-cell tests, 18% emitter bookkeeping; 42% of rows
  are twins (the MAC ladder above the cell's largest MAC, shrinking with N).
  **Decided: no BURDEN test memo on its own** (~12% of a 1-thread run, less
  threaded). The two memory levers (hoist the per-set N-copy of
  `genoURMat_noweights`; derive the zero-index vectors at test time) are
  the item, and a compute-once-and-copy memo keyed on the bitwise
  `genoSumVec` column rides along if that pass is restructured. Still open
  from F5: the >= 2-carriers-per-gene fraction, the null KS on `p_FlexRV`,
  a binary arm.
- **LEDGER #117** (a zero weight NAs a binary SKAT-O cell; pre-existing) is
  the owner's decision: drop weight-0 markers before the kernel, or leave it
  and keep FlexRV-on-binary a `--r.corr=1` run. **Measured 2026-09-11**: the
  harm is not one NA row -- the gene's Cauchy row silently becomes the
  surviving set's row verbatim (0.6787 where the control gives 0.5319). The
  fix target is now exact and needs no new golden: weighting a marker 0 must
  equal deleting it, which the BURDEN path already satisfies to every printed
  digit.
- **The FlexRV simulation study** (`docs/design/flexrv-sim-study.md`): the
  144-run S2-S4 grid ran 2026-09-10 (pre-#118) and was RE-RUN at HEAD
  2026-09-12 (`harness/sim/flexrv_study/results-2026-09-12.txt`;
  `docs/state/prediction-flexrv-rerun-2026-09-12.md`): `skato`/`grid`
  reproduce exactly, every FlexRV count fell by 1-2 (flex192 11 = grid 11 >
  `unw` 8 on `y_binary`; the collapse-0 column was the #129-free one, and
  since #129 is FIXED the collapse-10 numbers need re-running). Blocked on
  the owner's decisions -- memo `docs/state/flexrv-study-decisions-2026-09-12.md`.
- **LEDGER #118 is DECIDED (2026-09-11): the per-set MAC, everywhere** --
  variance-ratio band, ER eligibility, the min-group-MAC gate and the printed
  column. Measured over P58's 384 sets: 45.3% of printed MACs count markers
  the set zeroed (median 2.83x), 20.6% change VR band (but only 0.29% of SE on
  this model), and 8.3% would be REFUSED on the per-set MAC. It travels with
  two new rows and should be implemented as ONE piece of work: **#126** (the
  min-group-MAC gate has a single read in `src/`, inside the BURDEN branch, so
  it does not apply to SKAT-O -- which still reports a `Pvalue_Burden` column;
  decided: extend it, default 1, with the global-vs-SKAT-O-only sub-question
  open) and **#127** (the pooled Cauchy row records nothing about what entered
  it, so #117, a fully-refused region and #118's own refusals all degrade it
  silently). The implementation gate is in #126: two runs, at gate 1 and gate
  5, and the ONLY cross-path assertion permitted is that the SET of refused
  cells agrees -- never the p-values, which differ 2.76x at MAC 2 because the
  two paths do not share a burden estimator. The characterization in
  `test-step2-flexrv-reconcile.R` is inverted with it.
- `--flexRV_output_sets` (decided, not implemented); a native `sum_geno`
  collapse if F5 says the operator matters. The conditional reconciliation
  test is DONE (`test-step2-flexrv-cond.R`).
Everything below this line in this section is the pre-F2 history.


[Planned; `docs/design/flexrv-plan.md`, owner's decisions recorded there.]
The repo's "flexRV" label on Phase 5 names the multi-weight ENGINE, not the
method. FlexRV (Schwartzentruber et al., bioRxiv 2025.11.04.686320) is 192
fixed weight vectors per gene (16 pathogenicity-score transforms x 12 MAF/MAC
transforms), a burden test per vector, Cauchy-combined; nothing on this branch
builds those vectors. The plan adds a `score:NAME` group-file line and a
`transform` set kind that C++ expands from the in-sample MAF/MAC it already
has per marker (`src/Main.cpp:2012-2044`), custom-set semantics throughout,
one computation exported to R. Two things verified in code before planning,
recorded so they are not re-derived: the default `beta(1,25)` never multiplies
a custom set on the unconditional path (implicit beta set composed only when no
other set exists, `R/SAIGE_SPATest_Region.R:917-927`); it DOES reach the
conditioning markers' weights when `--weights_for_condition` is absent
(`src/Main.cpp:1852-1860`), so FlexRV + `--condition` without that flag is a
hard error in the plan. Memory scales as N x K x cells x ~4 arrays (9 GB per
array at N = 500k with 192 sets and the default 12 cells), so the plan's F5
stage measures at N = 50,000 on `sim50k` before anything is called done; the
sim's latent score `a_lat` is not written by `02_build_dataset.py`, which F5
has to fix or work around.

**Updated 2026-09-08.** **F1 is LANDED** (`97db5a88`): the group file parses an
optional `GENE score:NAME v1 v2 ...` line per region, validates it (numeric,
non-NA, inside [0, 1], one value per marker, the same name in every region) and
carries it into `RegionData` as `SAIGESC_1` and each region entry as `$SCORE`.
Nothing consumes it. Proved inert three ways: every golden byte-identical in
BOTH compare modes; an end-to-end region run with a score line byte-identical to
one without; and a planted out-of-range value exiting 1 through the real CLI --
that last is what shows the parser runs in the CLI path and not only in unit
tests. Two decisions NARROWED the plan and are recorded in it: the score line's
position is FIXED (immediately after `anno`), and `scoreName` is appended to the
`detail = TRUE` return ONLY when a score line exists -- because FOUR existing
tests pin the historical returns and messages, where the plan said two.

**Updated 2026-09-09. The paper blocker is GONE: the owner supplied the PDF,
F0 is LANDED (`722c3610`) and F2 is unblocked.** The 16x12 table is pinned
verbatim in `docs/design/flexrv-design.md` with each closed form derived and
reconciled two ways, and that file SUPERSEDES the plan's lossy-fetch table.
It resolves five ambiguities that a literal transcription would have got wrong,
of which one changes the program: the MAF beta transform is DECREASING in MAF,
`(1 - MAF/0.001)^(b-1)`, not the `(MAF/0.001)^(b-1)` a literal substitution into
the paper's own score-beta formula gives -- settled against Extended Data
Fig. 1b, whose axis is allele count 1-1000, not by intent.

A second doc commit (`ab047839`) records the F3/F4/F5 decisions taken while
reading it. The load-bearing ones: **one MAF cell and one annotation group per
FlexRV run** (`--flexRV_maxMAF`), because otherwise the region's pooled Cauchy
row is a super-omnibus over anno x maxMAF x set and is not `p_FlexRV`; the
**ultra-rare collapse collides with four of the twelve MAF transforms**
(`MAC<=1,3,5,10` sit at or below its default cutoff, so those tests are carried
entirely by the `max()` pseudo-marker and are not the paper's additive count),
with `--MACCutoff_to_CollapseUltraRare=0` as the exact fix for quantitative
traits and dropping the four transforms as the fix for binary; and **the
paper's AC x decile pseudo-variant binning is NOT ported** -- it is a STAAR
speed workaround approximating a test we can compute exactly. **F4 has a
fixture constraint**: at N = 10,000 and the goldens' default max MAFs most of
the MAF axis is degenerate, so a FlexRV golden on the defaults would gate
nothing.

`a_lat` still gates F5 only, and the plan carries a cheap workaround (anno
class plus seeded jitter, no 65-minute rebuild). **Three decisions are the
owner's and none blocks F2**: whether the quantitative/binary collapse split is
advice or an enforced default, whether the constant-score warning earns its
surface, and whether per-set rows are written by default (recommended: yes,
measured at 70 B/row). F3 is the commit that wants review: it changes `mainRegionInCPP`'s signature, which pulls in
`RcppExports`, the stale `Main.hpp`, `src/Makevars` and the delete-`src/*.o`
rule. **The plan's line references were re-resolved on 2026-09-08 and are now
SYMBOL-anchored** with the line number as a dated hint; they had drifted by up
to ~145 lines within a day.

## LEDGER #104/#105/#106: the SKAT-O burden variance (2026-09-07) -- DONE 2026-09-08

**All three are FIXED and committed** (`4ebc7c38` #106, `e338bdbd` #104/#105,
`690894f1` #103 before them), together with #44, which #106 re-opened and the
same commit re-fixed. `docs/state/handoff-2026-09-08-burden-variance-fixes.md`
is the account; the rows' STATUS cells are authoritative. Nothing below is a
decision any more -- it is kept as the record of what was decided and why.

[Diagnosed, reconciled two ways, fix verified;
`docs/state/diag-2026-09-07-region-burden-variance.md`, scripts and
reproduction in `harness/diag/region_burden_variance/`; the rows hold the
authoritative story. Companion measurement by the other 2026-09-07 session,
`docs/state/diag-2026-09-07-region-vs-burden-exact.md`: both burden
constructions against the exact model variance, the region path's excess on
multi-marker cells being entirely #106's truncation.] #104 and #105 are ONE mechanism and not the "which
variance ratio" question the rows first posed: on binary traits the R-side
`Phi` adjustment computes its burden p-value on the naive variance and the
clamp `r = min(VarQ/VarQ_2, 1)` then forbids the burden variance to go below
naive, so wherever the GRM and the ratio reduce it (every cell of P6) the
SKAT-O path reports the no-GRM, no-ratio burden test and inflates all of `Phi`
with it. That is a p-value defect (all three columns), not only an effect-size
one. The published SAIGE-GENE method has no such step. #106 is new and
separate: step 2's sparse `Sigma` is one triangle of the symmetric .mtx and
its PCG stops after one or two iterations on an absolute tolerance.

Decisions for the owner, all taken and executed:
* **Apply fix A** (four lines in `get_newPhi_scaleFactor_traitType`, reference
  `fixes.R` `phi_a`): standardise the burden statistic by `1'Phi 1` before
  the `|z| < 2` test and the saddlepoint, the marker SPA's own device. One
  commit for #104 + #105 with the characterization the row describes; the
  predicted golden movement is in the row (every binary SKAT-O group golden,
  every cell, plus the conditional group goldens; nothing else). Fix B
  (delete the `r` step, the published method) is identical under the null on
  this fixture but its `SE_Burden` understates the estimator's spread by 9 %
  under a planted effect, where A's matches; that is why A.
* **Then #106 as its own row and commit** (symmetrise the triplets in
  `setSparseSigma_new`, relative PCG tolerance), characterized on a P27 marker
  first; moves every sparse-GRM step-2 golden by up to ~0.5 % in variance.
* **Do not re-open** the "which ratio should scale a burden statistic" framing:
  with both modes on the sparse arm they still differ by the discarded GRM,
  4.4-4.9 %; the real per-row category mix is 1.5-2.3 % on the quantitative
  goldens and is the selector working as defined.

Two observations outside the ledger, recorded so the planted-effect table is
not misread: the step-1 refit recovers tau2 = 0.71 of a true 4.0 (0.19 of
0.33 at the model's own value) over hundreds of simulated replicates, the
known downward bias of the penalised quasi-likelihood fit for binary traits,
which bounds what any step-2 variance fix can achieve under strong
relatedness; and the score-test `BETA` overshoots by 16 % at a per-allele
log-odds near 1 (maximum likelihood with the true random effects recovers
0.97), on the marker path as much as in the region path. Neither is a code
defect on this branch.

## The N = 50,000 null region control (2026-09-08) -- DONE

`docs/state/null50k-2026-09-08.md`. The named NEXT UP job of the previous
handoff, and it is closed. Headline: on the `sim50k` `y_null` arm (quantitative,
sparse GRM, real relatedness, 210 genes, 3,338 region rows and 9,049
single-variant rows) the null is **calibrated at alpha = 0.05 and CONSERVATIVE at
alpha = 0.01**, with no anti-conservatism on any column or MAC stratum, and the
step-1 fit recovers the planted variance component to 0.7 % as the paired
positive control.

Three things in there that outlive the run:

* **`lambda_GC` is invalid on the SKAT-O columns.** They carry a point mass at
  `p == 1` exactly (26.7 % of the Cauchy rows; LEDGER #86, DECIDED). Dropping the
  atom reports lambda 1.47 on a perfectly acceptable null; keeping it reports
  0.73. Use the ECDF at alpha with a binomial interval. Same shape as the binary
  ER atom `harness/check_null.R` warns about, in a column nobody had checked.
* **The 2026-09-01 `step2_group_y_null.txt` is a CONFOUNDED baseline** -- ~80
  commits, not four fixes. Its `.singleAssoc.txt` is missing 255 rows in six
  whole gene blocks, and those six genes are exactly the six that escalate under
  `--is_fastTest`: that is LEDGER #50, fixed 2026-09-03. Attribution needed a
  scratch install at `cca6b497` (#103's parent).
* **#103's correction reaches -79 % at 50k**, against 1.0013-1.0367 on the
  bundled goldens: `kappa_w` has median 0.99928 (rare variants are slightly
  NEGATIVELY correlated, so `1 + (m-1) rbar` dips below 1) and max 4.826. It
  moves ZERO p-values, as designed. #106 lowered the burden variance by up to
  6.1 % and made 91-93 % of moved p-values SMALLER -- the pre-registered sign --
  while leaving `lambda_GC` to move only 0.721 -> 0.726, also pre-registered.

The BINARY SKAT-O path, where #104/#105 live, is covered by the follow-on
`docs/state/null50k-binary-2026-09-08.md`. Four arms verified at source level,
every falsifiable prediction held (#103 moves zero p-values; #106 with #104 still
present moves `SE_Burden` but not `BETA_Burden`; #104/#105 leaves
`.singleAssoc.txt` byte-identical), positive control intact (16/20 causal in the
top 20, ladder Spearman 0.682, both arms). Two things worth carrying forward:

* **#106 is ~20x smaller on binary than on quantitative** at the same N, GRM and
  genotypes -- burden-variance factor min 0.99685 against 0.9387. That is exactly
  what #106's own WHAT cell argues from `1/W ~ 11` dominating `tau2 x kinship`,
  now measured rather than argued.
* **#104/#105's `r` is 0.9878 median at 50k against 0.945-0.953 on the bundled
  fixture**, so the p-value movement is a median of 0.4 %, not the 12.6 % the
  goldens showed -- the golden headline is the maximum (the tail here still
  reaches 12.9 %), not the typical case. Predicted before the run from where `r`
  comes from, and held.

**Open, and NOT caused by the fixes: the binary tail at K = 0.01.** The 210
independent Cauchy rows of the binary null give 6 at p <= 0.01 against 2.1
expected (one-sided binomial p = 0.020) -- and 6 at `cca6b497` as well, so none of
it is attributable to #103/#104/#105/#106. The quantitative null gives 0/210 and
`y_binary`'s non-causal genes 1/189 on the same genotypes. Three realisations
reading 6, 0, 1 against ~2 is suggestive and no more. Settling it needs REPLICATE
phenotypes at K = 0.01, not another single arm; that is the next calibration job
if anyone wants one, and it is a question about SAIGE-GENE+ on rare binary traits
rather than about this branch's fixes.

## Dead-surface audit of the whole repo (2026-09-08)

Method, so it can be re-run rather than re-derived: every function definition in
`R/*.R` (133, excluding `RcppExports.R`) and every R-visible Rcpp export (156),
cross-referenced against `R/`, `extdata/`, `tests/`, `harness/` and `src/`,
matching on word boundaries so string dispatch (`do.call("f", ...)`) counts as a
use. **Its own control is `gen_spsolve_v4`**, which the audit did NOT flag --
correct, LEDGER #91 records it as the live sparse solve -- while its superseded
siblings `_v2`, `_v3` and `_inR` were all flagged.

| | count | size |
|---|---|---|
| R functions defined, never referenced anywhere | 16 | ~470 lines |
| Rcpp exports never called from R | 52 | -- |
| .. of those, used inside C++ (only the R WRAPPER is dead) | 35 | -- |
| .. of those, called from nowhere at all | 17 | ~1,438 lines |

Each of the 17 was verified by listing every occurrence of its name in `src/`:
all show exactly one, their own definition. All 17 are in
`SAIGE_fitGLMM_fast.cpp` -- superseded variants (`getPCG1ofSigmaAndVector_old`,
`gen_spsolve_v2/v3/inR`, `getCrossprodMatAndKin_usingSubMarker`), the mailman
multiplication experiment (`muliplyMailman_NbyM`,
`getCrossprodMatAndKin_mailman`, `mmGetPb_*` on the wrapper-only side) and
`testTime`.

**The behavioural part is LEDGER #114 and is the only piece written up as a
row**: two of the 17 are the ONLY setters of two globals, so `isUsePrecondM` and
`isUseSparseSigmaforInitTau` are frozen at `false`, the preconditioner-M `else`
branch in the PCG is unreachable in any build, and one disjunct is inert.

The rest is size, not behaviour, and is left as a tidiness item because deleting
it moves no golden and answers no question. If it is ever done, two things are
worth knowing. The 35 wrapper-only exports are the cheaper win: dropping their
`// [[Rcpp::export]]` shrinks `RcppExports.*` and stops exposing internals,
without touching any live code path -- but it requires the
`Rcpp::compileAttributes` + `src/Makevars` dance in CLAUDE.md. And
`region_postprocess_parallel` (`R/SkatoParallel.R:204`) is among the dead R
functions while its own roxygen justifies keeping it "for the unit tests of the
three pieces" -- measured FALSE: the tests reference `skatoCollect`,
`skatoEvaluate` and `skatoServe` 3, 4 and 5 times and the composition zero. A
comment is not a caller.

## DECISION (2026-09-08): incompatible flag combinations should be REFUSED

The owner's call, recorded here because it is a policy rather than a single fix:
**where two options cannot both be honoured, the CLI should refuse the
combination rather than silently reconcile it.** The motivating case is LOCO
with a sparse GRM -- LOCO leaves one chromosome out of a full GRM and a sparse
GRM has no chromosomes to leave out, so the two do not make sense together.
Today `R/SAIGE_fitGLMM_fast.R:836` quietly assigns `LOCO <- FALSE` and says
nothing. LEDGER #111 and #112 carry the detail.

Two things to settle before implementing, both measured:

* **A hard error breaks the DEFAULT invocation.** `--LOCO` defaults to TRUE in
  BOTH CLIs, so `--useSparseGRMtoFitNULL=TRUE` with no `--LOCO` at all is the
  recommended workflow's default command; refusing it means every existing
  sparse-GRM script must add `--LOCO=FALSE`. The alternatives are to make
  `--LOCO`'s default depend on the GRM choice, or to warn loudly and continue.
  This is the choice of who breaks, not whether to fix.
* **LOCO is the ONLY silent override in `fitNULLGLMM`**, so this is one site and
  not a pattern. Every other in-body reassignment of a user-supplied option
  announces itself: `skipVarianceRatioEstimation` (:802),
  `useSparseGRMforVarRatio` (:830), `isCovariateOffset` (:960) and `minInfo`
  (`R/SAIGE_Test_main.R:451`). The neighbour `nThreads <- 1` is announced only
  `if (bedFile != "")`, so a `--plink2File` sparse fit silently overrides both
  LOCO and nThreads with no output at all -- measured on the golden logs, PLINK
  prints the notice once and the PGEN twin zero times.

Whatever is chosen, step 2's remedy sentence ("please run Step 1 using LOCO")
needs the sparse-GRM case added, because on that path it cannot be carried out.

## Golden coverage audit (2026-09-08)

[Measured, not estimated: the 38 step-2 golden commands parsed out of
`harness/make_goldens.sh` and crossed against trait type and test type, then all
133 CLI options swept for ones no golden touches, then checked against the test
suite. Reproduce with the same parse of `make_goldens.sh`.]

**Closed by this audit:** the quantitative sparse-GRM step-2 path had NO golden
(23 of 33 binary runs passed `--sparseGRMFile`, 0 of 5 quantitative). Fixed the
same day with `quant_sparseGRM_vr` and goldens P38-P41; the pre-#106 build moves
that path's `BETA_Burden` by up to 64%, so it was worth having.

What the audit found and is NOT done:

* **Missing genotypes were never exercised -- CLOSED 2026-09-08 for hard calls.**
  Re-verified by decoding every bundled .bed, which is worse than the audit said:
  all SEVEN files have ZERO missing calls (~150 million), so step 0 and step 1
  were as blind as step 2. Closed with
  `extdata/input/genotype_100markers_missing.{bed,bim,fam,pgen,pvar,psam}`
  (generated by `harness/tools/make_missing_fixture.py`, byte-reproducible from
  that script; only the `.bed` differs from `genotype_100markers`), goldens
  P42-P48 and `tests/testthat/test-step2-missing.R`. Three things it found, now
  LEDGER #107/#108/#109: `--minInfo` is a switch at 1 rather than a filter
  (both readers hardcode `imputeInfo = 1`); `--is_imputed_data` REPLACES
  `MissingRate` with that constant rather than adding a column; and
  `--dosage_zerod_cutoff` at its 0.2 default silently turns
  `--impute_method=mean` into `minor` on every `MAC <= 10` marker, which is the
  interaction the documentation promises is gated on `is_imputed_data` and is
  not. DOSAGES are closed too, with
  `extdata/input/genotype_100markers_dosage.{pgen,pvar,psam}` (228 KB;
  `harness/tools/make_dosage_fixture.py` + one `plink2 --import-dosage`, the
  2 MB text intermediate reproduced rather than committed) and goldens P49-P51.
  That fixture is the only place the pgen DOSAGE DECODE path
  (`PgrGet1D`/`Dosage16ToDoubles`) runs in a golden at all, and the only place
  `--dosage_zerod_MAC_cutoff` is moved. Both contrasts were measured before the
  goldens were accepted and each differing set is IDENTICAL to its MAC band, not
  merely inside it: cutoff 0.2 vs 0 moves exactly the 59 markers with MAC <= 10,
  and MAC cutoff 10 vs 100 moves exactly the 17 with 10 < MAC <= 100.
  `--hard-call-threshold 0.499` on the import is deliberate -- every dosage also
  gets a hard call, so the fixture adds no missing genotypes and varies one
  thing. **The guard that matters there** is that every `AC_Allele2` must be
  FRACTIONAL: if the decode were bypassed, or plink2 had hard-called on import,
  every other assertion would still pass against what is effectively a hard-call
  fixture. Still not covered: dosages AND missingness together, i.e.
  `--impute_method` on fractional data.
  **Trap paid for once:** `missingRate = numMissing / N` is computed over the
  MODEL's sample subset (fam rows 1..1000), not the .fam's 10,000 rows. A first
  fixture sized its ladder against 10,000 and SAIGE reported 0.032/0.096/0.161
  where 0.05/0.10/0.15 were meant, putting the boundary marker on the wrong side
  of `--maxMissing`. Size any genotype-QC fixture against the analysed subset.
* **`LOCO=TRUE` was in zero golden runs -- CLOSED 2026-09-08.** Step-1 golden
  `binary_loco_2chr` plus step-2 goldens P52/P53 (the same model and command with
  only `--chrom` changed, so the pair isolates `LOCOResult[[chrom]]`), and
  `tests/testthat/test-step2-loco.R`. Only a DENSE fit can carry a LOCOResult at
  all: `useSparseGRMtoFitNULL=TRUE` forces `LOCO <- FALSE`
  (`R/SAIGE_fitGLMM_fast.R:836`), so four of the five step-1 goldens never could.
  **The harness had to be fixed first, and this is the reusable lesson:**
  `compare_step1.R` descended into a list only when the field was named
  `obj.noK`, so `LOCOResult` was skipped SILENTLY -- the golden would have gated
  theta and the base fields while comparing not one value of the thing it exists
  for. Drilled: perturbing `LOCOResult[[2]]$residuals[1]` by 1% is now a FAIL and
  was invisible before. **Fixture chosen by measurement**, not by reusing the
  other dense golden's: the 22-chromosome fixture gives a 4.7 MB `.rda` (against
  ~268 KB for every other step-1 golden) and chr1-vs-chr2 residuals differing by
  1.5e-3, while `plinkforGRM_1000samples_10kMarkers` gives 677 KB, a 7 s fit and
  a 0.19 difference -- 130x the signal, because it is 9,988 markers on chr1
  against 12 on chr2. The mechanism demonstrably reaches the output: same model,
  same chromosome, only `--LOCO` changed moves every marker (chr1 median 6.5%,
  max 52%; chr2 median 1.3%, max 6.0%), and the asymmetry's direction was
  predicted from the marker counts first. The PGEN mirror of a LOCO fit is
  byte-identical to its PLINK twin, which no earlier identity check covered.
  Found LEDGER #110 on the way. Still uncovered: `--isLowMemLOCO` has no GOLDEN
  (only `test-step1-lowmem-loco.R`), and no REGION test runs under LOCO.
* **`--isLowMemLOCO` was assessed 2026-09-08 and a golden was DECIDED AGAINST.**
  It is a memory-layout switch: with LOCO it writes `<prefix>_noLOCO.rda` plus
  one `<prefix>_chr<j>.rda` per chromosome instead of one file holding all 22
  refits, so only one refit is in memory at a time. It is dead on the RECOMMENDED
  path -- a sparse-GRM fit forces LOCO off -- and on the dense path it changes
  binary answers for chromosomes 2+ through the known warm-start divergence
  above, which this session measured reaching **BETA 3.1% and p-value 5.6%** in
  step-2 OUTPUT (the existing record measured model fields at ~1e-5; the larger
  figure is [ASSUMPTION: not tested on a balanced fixture] because
  `plinkforGRM_1000samples_10kMarkers` is 9,988 markers on chr1 against 12 on
  chr2, so chr1's refit -- the warm start -- sits far from the base). Measured
  layout on that fixture: 11 s and three files totalling 789 KB, against 7 s and
  one 677 KB file, so on two chromosomes it is SLOWER and BIGGER; the win is peak
  memory and scales with chromosome count (~204 KB per refit at N = 1,000,
  O(N) content, so [ASSUMPTION: linear extrapolation] ~2 GB at N = 500,000).
  Step 2 has no logic to derive `_chr<j>.rda` -- the name is built in exactly one
  place, the step-1 writer -- so the user must pair the file with `--chrom` by
  hand; a mismatch dies with `Not compatible with requested type: [type=NULL;
  target=double]`, which names nothing relevant (the writer stubs other entries
  with `list(NULL)`, and `is.null(list(NULL))` is FALSE, so the guard at
  `R/readInGLMM.R:248` passes and reads NULL fitted values out of a stub). Exit
  1, so no wrong answer escapes. Worth a row if anyone decides to support the
  flag properly; not written, because the flag is dead on the recommended path.
* **`r.corr` only ever takes 0 or 1 -- and that is a CLOSED DOOR, not a coverage
  gap** (measured 2026-09-08). `--r.corr=0.1/0.5/0.9` all exit 1 with "r.corr
  needs to be either 1 (BURDEN test) or 0 (SKAT-O test)"
  (`R/SAIGE_SPATest_Region.R:815`). SAIGE reuses SKAT's parameter NAME but not
  its semantics: `r.corr = 0` does not mean "SKAT", it means "run SKAT-O", which
  then sweeps the whole rho grid internally (`SKAT_Check_Method(method =
  "optimal.adj")` returns the grid and it is assigned straight back into
  `r.corr` at `:802`). So the option is a two-valued MODE SELECTOR and there is
  nothing between the values to gate. The roxygen at
  `R/SAIGE_Test_main.R:97` still says "numeric. bewteen 0 and 1", which implies a
  continuum that does not exist (and misspells "between") -- #111's family, not
  written up as its own row.
* **`binary_denseGRM` is a step-1 golden that is never a step-2 model.**
* **Quantitative had no `fastTest`, multi-weight or conditional REGION golden
  (binary has 14, 9 and 9) -- two of the three CLOSED 2026-09-08.** Goldens P56
  (`plink_group_multiweight_quant`, K=2) and P57 (`plink_group_cond_quant`), both
  on `quant_sparseGRM_vr`, the only quantitative model with `sparse`
  variance-ratio rows, so they run the same sparse-GRM region path #106 was
  first-order on. P56 brings the K>1 schema (`Weight`, `Pvalue_ACATV`,
  `Pvalue_ACATO`; 26 rows against P38's 14) and P57 the five `_cond` columns --
  neither existed on any quantitative golden.
  **The third was deliberately NOT added, and the measurement is the reason:**
  a quantitative `--is_fastTest=TRUE` region run is BYTE-IDENTICAL to P38, just
  as the binary `plink_group_fastTest`/`plink_group_nofast` pair is, because the
  region path never reads the flag (LEDGER #51). A third file carrying P38's
  bytes would gate nothing, so that invariant is asserted in
  `tests/testthat/test-step2-group-quant.R` instead, where it states for the
  quantitative path what #51 records for the binary one. Still missing on the
  quantitative side: PGEN mirrors of P56/P57 (reader parity for the quantitative
  region path is covered by P41 alone), and the BURDEN-mode twins.
* Step-1 options with neither golden nor test, mostly numerical knobs:
  `tolPCG`, `maxiterPCG`, `maxiter`, `tauInit`, `nrun`, `traceCVcutoff`,
  `ratioCVcutoff`, `isCovariateTransform`, `isCovariateOffset`,
  `isDiagofKinSetAsOne`, `skipModelFitting`, `memoryChunk`,
  `usePCGwithSparseGRM`. Step-2: `weights_for_condition`.
* Recorded covr baseline, not re-measured: R 55.2%, live step-2 C++ 72.2% line
  and 42.9% branch.

## Deliberately not listed

The consolidation work (README is still the upstream stub, `pixi.lock` is
untracked, bench logs are in git) is real but needs a decision about who this
branch is for, which is not a decision to make from inside the code. The
remaining speedup candidates named in the older plan -- step-1 variance-ratio
estimation blocking, `setgeno` ingest -- are left where they are until the
profile above exists.

- **Binary AI-REML oscillates on the sim50k fixture (found 2026-09-10 while
  fixing LEDGER #121; provenance: measured, `harness/work/flexrv_study/`
  `step1-fix121.log` and the scratchpad `s1fix/` traces).** After #121 the
  fit no longer dies at the boundary, but every binary trace at 1%
  prevalence with the pedigree sparse GRM swings 3-5x per iteration and the
  h2 = 0.2 trait hits `maxiter` unconverged; identical at every
  `traceCVcutoff`, so it is the (alpha, tau) fixed-point iteration, not
  trace noise. Not a bug that a golden can see. Candidates, none tried:
  damp the AI step (GMMAT does not), tighten `tol.coef` / `maxiter` inside
  `Get_Coef`, or check whether the PQL working variable at tau near 0 is the
  driver. **Measured later the same day** (handoff, "Convergence of the
  binary AI-REML"): CLI knobs do nothing; a damped step (0.5 or adaptive
  sign-flip halving) with the criterion on the undamped proposal and
  tol 0.0002 converges from every start to 0.970 on `y_binary` (SAIGE
  reports 1.31 / 0 / 0.19 by budget) and to 1.32 on the h2 = 0.2 trait
  (SAIGE: 0). **LEDGER #122 and #123 are FIXED (2026-09-10, late):** the
  bracket design was refuted by its first drill (stale score signs), what
  landed is adaptive damping armed on the second direction reversal
  (`R/TauStepGuard.R`), zero as an iterate, convergence on the raw proposal
  never on iteration 1, twice once damped; every golden byte-identical, the
  50k trait 0.968-0.973 from four starts. What stays open here: the 5%
  planted trait's start-0.3 run stops on a single pair at 0.671 against
  0.626/0.634 from other starts (today's 2% criterion on a monotone tail,
  left as is on purpose), and the trace estimator's own noise (19% of a
  binary tau2 at N = 1,000 with the default 30 probes), which is the
  estimator-test item below.
- **Step-1 estimator tests (designed 2026-09-10, not built; handoff
  "Estimator tests for step 1").** Planted traits from the bundled GRM,
  dense exact references, tolerance = 3 x the trace estimator's own se
  (computed exactly). Measured: both fitters reach the exact fixed points
  up to that noise; the default 30 probes make it 19% of a binary tau2 of
  0.13 at N = 1,000 and the `set_seed(200)` makes it repeatable, so a
  same-sign offset in tau2 across traits is the trace, not a bug (this
  cost an hour). LEDGER #124 (AI projection slip) fell out of it.
- **LEDGER #124 FIXED (2026-09-10, night, owner's go-ahead; provenance:
  measured, the row and the addendum in
  `handoff-2026-09-10-flexrv-sweep-124.md`).** Every golden re-baselined
  (one regeneration, byte-reproducible twice, the moved-file prediction held
  file for file). What this leaves: `make_goldens.sh` produces three step-1
  pgen twins and `make_controls.sh` two control pgen files that are not in
  the committed set, and `plink_marker_missing_imputed.txt.index` is a stray
  committed index with no run -- tidy or leave, owner's call.
- **LEDGER #125 OPEN (2026-09-10, night; provenance: measured, the row).**
  The #122 guard is inert on a cycle whose slow phase passes the 2%
  single-pair criterion before a third reversal: the planted 3% trait
  reports 1.28 / 1.24 / 0.89 from three starts under the corrected AI, and
  #122's "three starts agree" was an accident of the wrong curvature. The
  MECHANISM FOUND the same night: the inner (alpha, b) pass converges only
  to the shared `tol.coef = tol = 0.02`, which makes the outer REML score
  four times noisier than the exact score; with `tol.coef = 1e-6` the three
  starts agree (1.085-1.094), the bundled fit converges in 2 iterations, and
  the 50k trait in 7 instead of 16 to the same answer. The 50k trait under
  #124 alone still converges to 0.967-0.973 from four starts. [DECISION]:
  decouple `tol.coef` from `tol`; every step-1 golden moves again, the
  bundled fit by 1.6% (beyond the 1e-2 theta tolerance), and #122's guard
  becomes second order. Cost sweep done (handoff table): total inner solves
  rise by two or three while outer iterations halve, wall time never rises;
  recommended default `tol.coef = 1e-5`. **LANDED the same night as
  `--tolCoef` (default 1e-5), every golden re-baselined again**; #122's guard
  is now second order and could be retired or kept as belt-and-braces --
  owner's call, no rush. Details in the row and the handoff.
- **FlexRV study S2-S4 RUN (2026-09-10, night; provenance: measured,
  `harness/sim/flexrv_study/results-2026-09-10.txt`, NOTES.md).** Null
  calibration fine; the positive control FAILED on `y_quant` (no method
  detects the planted quantitative signal at 50k); two design expectations
  mis-specified (the uninformative arm keeps the `lof` set and the MAF
  ladder; the `low` class carries causal variants). [DECISION] x3 before the
  power question is answerable: raise the quantitative NCP, define the
  negative control for "score informativeness" against the best score-free
  set rather than `unw`, and decide whether `low` should be effect-free.
  Before any further sweep: REFIT the four `sim50k` null models, which were
  fitted before #124/#125 (the `y_binary` model at tau2 0.9728 will move to
  ~0.9705 and every step-2 p-value with it), and re-run `study.sh`.
- **`compare_step1.R` now floors its elementwise ratio at the field's scale
  (2026-09-10, night)** for the model's vector/matrix fields; theta and the
  variance ratios are unchanged. Nothing at HEAD changes; it exists because
  #124's gate showed the floorless ratio calling a 1e-9 move of a 1e-8 entry
  15%. Reverse it if the owner prefers the old reading.
