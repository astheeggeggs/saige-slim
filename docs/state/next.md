# What is worth doing next (2026-09-03, updated 2026-09-29; read docs/state/handoff-2026-09-29-step-back.md first: it is short and says what to drop)

**Live (2026-09-29, night): STEPPED BACK.** LEDGER #172 is the real problem and it is fixed and verified on sas. The pinned-sample guard (#174) is built, gated and frozen; its hazard has never been observed on All of Us outside the cells #172 fixes, so the work on its threshold is DROPPED. "The seven" are the models with samples at or below 1e-8: afr FemInf and MatHem, amr AMD and FemInf, eas FemInf and MatHem, eur MatHem. NEXT: the owner's answer on what to do with them (caveat, exclude, or refit on a restricted cohort); the eas stage check, control and comparison; then amr, afr, eur on the #172 image. Nothing is pushed (CI minutes). The paragraphs below are earlier states.

**Live (2026-09-29, evening): LEDGER #174's guard LANDED (`08c7e7a5`) and its THRESHOLD is an owner's decision again; #173 DECIDED with a safety net (`2c6cd345`); #175 OPEN (the fit gate's separation check has #174's blind spot).** The fixture shows a separated group is fitted wherever glm stopped -- 1.7e-9..1.2e-7 for 60 samples of 1,000, 13 below 1e-8 and 47 above -- so the agreed 1e-8 drops part of the group and the p stays wrong (0.030/0.054 against a reference of 0.0196), while 1e-6 gives 0.0210 and the agreed budget then refuses. On the planted model (group at 1e-12) the guard is complete: 7e-8 and 2e-13 become 0.0209. NEXT ACTION, the owner's: run `harness/tools/aou_mu_decades.sh` over the real models (`aou_spa_rerun.md` section 11) to learn where their separated groups sit, then choose the threshold and the budget's meaning. NOT DONE: no image built from these commits; the splice control does not read `.spaPinned.txt`; the guard has no replicate null run. Nothing is pushed (CI minutes). The paragraph below is this morning's state.

**Owner's decisions open from 2026-09-29: #173 (marker path scores the centred genotype, saddlepoints the projected one) and #174 (pinned samples stay in the saddlepoint vector).** **Live (2026-09-29): LEDGER #172 FIXED (one-atom weighted sets take the exact test; goldens P69-P71; the AoU binary FlexRV rerun is the next action, `aou_spa_rerun.md` section 9). The paragraph below is the day-before state.** **(2026-09-28): LEDGER #172 OPEN** -- three AoU FlexRV cells (amr/eas) with a CONVERGED saddlepoint p of 1e-53..1e-93 where ~1e-3 is right; never rerun because the scan flags failures only, and every kernel change since the v4 image is on a failure path, so the current image is PREDICTED to reproduce them. Diagnostic landed: `SAIGE_SPA_DUMP_BELOW` + the dump's kernel line and `nz`/`g` columns + the replaying `aou_spa_dumpcheck.sh` (`harness/tools/aou_spa_rerun.md` section 8 is the recipe: rerun the three regions with the dump on, read `RUN:` / `FAST replay` / `MONTE CARLO`, then fix). Do NOT rerun every SPA cell before that. Found on the way, **#173 OPEN** (marker path: the score is the centred genotype's, the CGF the projected one's; owner's decision, upstream-faithful vs consistent). #171 still OPEN.

**Current handoff: `docs/state/handoff-2026-09-26-flexrv-signature-repro-negative.md`** -- the FlexRV signature RESOLVED as LEDGER #163-#168 (saddlepoint root finder, the normal-p fallback, ER eligibility of every weight set, the exact test's tie epsilon, the mirrored q), five images, the AoU scan/filtered-rerun/splice tooling under `harness/tools/aou_spa_*`, SAS five of six traits spliced; open: 24 BenIntNeo cells whose mirrored root the finder still misses (dump run launched, analysis recipe section 5), then the other four ancestries. The previous handoff follows.

**Current handoff: `docs/state/handoff-2026-09-26-oldnew-comparison-and-flexrv-signature.md`** -- the v4 results compared against the 2024-03-26 BRaVa run at BOTH levels over the 127 trait x ancestry pairs whose N and case count match to 0.5%. **The two runs agree except for one trait.** Variant level: every binary MAC band at `rho = 1.0000`, median `log10(p_new/p_old)` 0.0000, `frac_10x` 0; the join is exact (`MarkerID` is `chr:pos:ref:alt` both sides, zero allele swaps, AF to 2.4e-5 = ~3 allele counts at `1/(2N)`). Gene level: the only movement is what was shipped on purpose -- the Liu screen above 0.05 and the saddlepoint below 1e-5 (#158), with the like-for-like band 1e-5..0.05 holding ~400 discordant cells of 5.2 million and `eur/Height` the visible tail case at p_new/p_old ~ 1.21. **The exception is `amr/BMI`**, which owns 92.8/99.9/99.9% of the gene-level discordance; calibration says the OLD run is broken (obs/exp 13.8x at alpha 0.05 rising to 3383x at 1e-4 -- 69% of cells below 0.05 -- against 1.25x-2.76x for the new one), and **the mechanism is NOT the 6.91 variance ratio**, which the 2026-09-25 handoff established was a truncated-GRM artefact (refitted AMR null ratios 0.970-1.966, median 0.997). **OPEN and the reason to read the handoff: a FlexRV signature that needs a repo repro.** On `afr/BenIntNeo` three genes of 18,432 sit at p_FlexRV 7.8e-120 to 2.3e-256 with a ~114-order gap to the next, while the trait's own calibration is clean and its single-variant p-values agree with the 2024 run at rho > 0.99. Every extreme set row has `Number_rare = 0` -- the cell holds nothing but the collapsed ultra-rare pseudo-marker, which the MAC1/MAC5/MAC25/MAC50 rungs give by construction -- with `BETA_Burden` 260588 against `SE_Burden` 7619 on MAC 86, an impossible |Z| = 34; it is NOT separation, and scale invariance holds where it should. Since `p_FlexRV` is a Cauchy combination it is ~the minimum over the sets and inherits these, so no FlexRV hit list is usable until each hit's driving set is checked. **It reproduces or falsifies with NO AoU data**: bundled fixture, `--flexRV_score`, a cell of only ultra-rare markers, a steep transform. Also there: the 2024-03-26 files double every UR pseudo-marker row (27.7% of rows) and carry a rebuilt `MarkerID`, variant MAC starts at 11 not 1, and a third egress alert on 2026-09-26 landed on a bulk bucket -> notebook-VM read that renders nothing, which is evidence for that channel over "displaying output is the egress". The previous handoff follows.

**Current handoff: `docs/state/handoff-2026-09-25-aou-production-step1-step2.md`** -- the AoU production run: 167 step-1 models refitted and step 2 run over both arms (FlexRV and a six-mask SKAT-O) for all 22 autosomes of every model. **AMR's two problems separated and both fixed, independently confirmed.** A TRUNCATED step-0 GRM (52-61% of samples; the correct files are under `run/make_sparse_grm/`, not `step0/`) caused the step-2 hang, and the signature everyone was reading -- AMR's null variance ratio of 6.91 -- was that artefact: refitted on the full GRM its null ratios are 0.970-1.966, median 0.997 over all 42 models, and 2.004 max at the OLD cutoff too, so the truncation and not the cutoff produced the 6.91. **So the recorded verdict "AMR CANNOT use a sparse GRM in step 2" is VOID** and every ancestry now passes `--sparseGRMFile`; `sparse` ratios are 0.99977-1.0014 across all five. Separately the 0.05 relatedness cutoff caused the step-1 non-termination: AMR's 5,000-marker GRM has 46,204,255 entries at 0.05 (about 1.8-2.1 sigma of the off-diagonal noise, so mostly noise floor) and 349,445 at 0.125, where HTN fits in 58 s having never finished, with tau2 0.0741 and a real random effect in 33 of 42 models. **The one thing with no backstop anywhere: `--relatednessCutoff` must MATCH between the steps.** Step 2 applies its own on load (`R/SAIGE_SPATest_Region.R:61`) and BOTH wrappers default it to 0, so a 0.125 model tested at the default would use 46.2M entries against a variance ratio describing 349,445, cleanly and wrongly; the other four ancestries are safe because step 0 wrote their `.mtx` AT 0.05, so both 0.05 and 0 are no-ops on it. Three LEDGER outcomes: **#161 CONFIRMED at production scale** on a build already carrying its fix (two AMR fits print `tau2 = 0 by the floor` with Newton steps 0.01295 and 0.01959, both inside the predicted band, both leaving zero at `--tol 0.01`; controls clean; owner keeps tol at 0.02); **#158 VALIDATED at production scale** (200,652 `liu-screen` cells, every one above 0.05 and none below -- the screen's contract on real data, previously only checked against synthetic nulls); and **#162 RAISED**, the `skatoMethod` sidecar covering only the 231,509 cells ASKED of the kernel against 329,393 emitted, with the 30% shortfall structured entirely by max MAF. Calibration evidence for EUR T2Diab is on the handoff and is the shape you want: `synonymous` inside its binomial interval at every alpha while pLoF runs 1.12/1.46/2.87 -- signal, not inflation -- and lambda_GC is unusable because burden p-values have atoms. Two file-format traps recorded there: SAIGE writes combined masks with `;` although the command line takes `:`, and the pooled `Group == "Cauchy"` row has every column but `Pvalue_Burden` blank. Open from it: #162 unexplained (the first step is empirical and on the row), zone attribution for the two new preemption batches not collected, the QQ rollout built but run only for EUR T2Diab, and the `tauInit` roxygen fix from the #161 diagnosis still outstanding. The previous handoff follows.

**Current handoff: `docs/state/handoff-2026-09-23-fit-gates-and-aou-refit.md`** -- the step-1/step-2 eligibility gates landed (`db5bb558`, `3e61864a`; `R/FitGates.R`, five options -- `minCaseCount` 100, `minCovariateCellCount` 10, `separationCheck` + `separationMaxFraction` 0.02, `convergenceCheck`, `minSampleCount` 0 -- with `SAIGE_FIT_GATES=0` set in `harness/env.sh` because the bundled fixture fails two of them by construction: 98 cases, and `a10` level 2 has 7). Gate FAIL 0 / PASS 6009, no golden moved. **Refusing instead of silently saving exposed LEDGER #161, a strip-branch REGRESSION, FIXED the same evening** (`python3 harness/tools/ledger.py 161` has the whole story): every default binary fit enters the AI-REML loop at tau2 = 0, and a positive score there whose Newton step is below `tol` was floored to 0 by `tauStepBinary`'s GMMAT rule while the #122 boundary exit needed `score <= 0` -- identical inputs every iteration, the full `--maxiter` at (1, 0), `converged = FALSE`, 21 of 164 real AoU binary fits refused. The `AI`-sign hypothesis was REFUTED (AI is a PSD quadratic form; measured 94-115). The exit is now keyed on the floored step being 0, the log says which case fired, and the damping guard is untouched. Upstream reports the same theta (0) on those traits, so the AoU refit of the 21 needs only the new image. Loose end from the diagnosis: the `tauInit` roxygen in `R/SAIGE_fitGLMM_fast.R` says `0,0` becomes 0.5 for binary; the code (upstream's too) sets 0 -- a docs-only fix plus `harness/docs.sh`, and a docdrift rule once corrected. Two other things from the same day: the BRaVa wrapper's `--tol` was `0.2` for binary (mislabelled "SAIGE DEFAULT") and `1e-5` for quantitative where upstream and the public All-by-All pipeline (`github.com/atgu/aou_gwas`) both use 0.02 for both and All-by-All never sets it; and a sparse-GRM step-1 fit is SINGLE-THREADED whatever `--nThreads` says (`R/SAIGE_fitGLMM_fast.R:918`), so size those machines for RAM. The previous handoff follows.

**Current handoff: `docs/state/handoff-2026-09-17-skato-pvalue-methods.md`** -- LEDGER #158 FIXED 2026-09-18 in two pushed commits (`11d09a7b` the outer-integral Liu screen, `7054100a` the saddlepoint tail; `R/SkatoScreen.R`, options `--SKATO_screen_threshold` 0.05 / `--SKATO_tail_method` saddlepoint / `--SKATO_tail_threshold` 1e-5; (0, davies) = upstream byte for byte, golden P67). **Post-#158 profile MEASURED 2026-09-19 (`docs/state/handoff-2026-09-19-post158-profile.md`): kernel 40.4 s = 42% of user CPU, C++ region buckets 52.4 s = 54%, overhead 3.8 s, on a quiet box with a BURDEN negative control pinning the overhead term. #158's whole-run effect is 2.14x (207 s -> 96.7 s), under the composed ~2.5x because the shipped threshold is 0.05 not the 0.01 the tables used. So making the ENTIRE kernel free is worth 1.72x at N=50k and 1.09x at N=400k, where the O(N) buckets scale ~8.5x and the kernel does not scale at all. (This paragraph originally ranked carrier-space Phi first on that reasoning; it was REFUTED the same day -- see below -- because the bucket it touches is 3.5% of the run, not 54%.) NOT settled: the `region_postprocess` double drive, which needs a THREADED arm (`canParallelSkato(1)` is FALSE, so `--nThreads=1` forks nothing and `profile_region.sh` hardcodes 1).** Open from it: a compiled saddlepoint if biobank-scale cell counts below 1e-5 make the R one visible; ~~a null sweep per trait on the finished build as the regression guard~~ **DONE 2026-09-19, PASSES** (`docs/state/handoff-2026-09-19-158-finished-build-guard.md`): 96 arms, 244,063 paired cells, ZERO misses by this repo's definition (`Davies p < threshold/10`). The old `nullrep/dump_arm.R` patch point is STALE -- it hooks `SKAT:::Met_SKAT_Get_Pvalue`, which a screened cell no longer reaches, so it would have scored the unscreened tail only and read as a pass; the guard hooks `SAIGE:::get_SKAT_pvalue`. Two numbers to carry: the worst skipped Davies p is 0.0475 at 244k cells (0.049/0.050 was the ~1.5k-cell figure), still nothing below 0.01; and the null yields 3 unique cells below 1e-5 in 244,063, both saddlepoint ones at q = 9 and 10 -- so a null sweep CANNOT validate the tail, and the q = 4-10 exact-reference control is the instrument that can -- **BUILT AND RUN 2026-09-19, `docs/state/handoff-2026-09-19-saddlepoint-tail-validation.md`, harness `harness/skato_pvalue/tailref/`: the saddlepoint tail IS accurate at small n.** Ratio to a cancellation-free Ruben reference is 1.036 at p = 1e-2 rising to 1.11-1.12 at 1e-12, CONSERVATIVE throughout, and FLAT in the number of terms (median 1.066-1.077 across n_lam 3..10), so Kuonen's `o(n^-3/2)` worry does not materialise; the bias is a function of depth, not of n. Positive control: Davies agrees with the reference to 6-7 significant figures at 1e-2..1e-4, is 1.58x out at 1e-6, 4x median / 18x max at 1e-8, and returns no positive value at all below. Carry this one: over 840 points `SKAT_davies` returned a NON-POSITIVE probability on 277 (33%) while reporting `ifault = 0` every time -- check the SIGN, its fault flag is not a guard. ~~the OUTER rho integral is not validated by this~~ **VALIDATED 2026-09-19 (`docs/state/handoff-2026-09-19-outer-integral-and-escalation.md`): it does NOT amplify the inner bias -- <=0.5% at p ~ 1e-4 and <=0.01% below 1e-8, because the survival-free `pchisq(xstar, df=1)` term carries 92% of the p at 1.5e-4 rising to 99.8% at 1.1e-14. The `pmin * n.r` cap was checked and does NOT bind, so the agreement is not vacuous. Open below ~1e-30, where the pchisq term underflows and the integral carries everything yet the arms still agree to 5 significant figures -- observed, unexplained, and far outside anything reported.** The same handoff costs the fastTest ESCALATION lever, which looked like 19.4% of C++ (39 of 249 calls, 10.86 s): only the 4.99 s of genotype prep is reusable, the 5.87 s of tests is legitimately re-run fast -> non-fast, so it is **1.054x at N=50k / 1.095x at 400k**, not worth touching the escalation control flow for. The tail control is a desk instrument; its fast-tier guard is `tests/testthat/test-unit-skato-tail.R`; LEDGER #159 (few-case SKAT inflation) decided, left as is. The pgen-sparse-read handoff below is the previous one. **Research note 2026-09-18: `docs/state/research-2026-09-18-skato-speedups.md`** -- a 5-angle literature sweep on what accuracy-preserving speed-up is LEFT on the SKAT-O path. Headline: the landed pair IS the literature's design (no faster Davies exists; moment matching is anticonservative below ~1e-5 so it is only ever legal as an above-threshold screen; Kuonen is the only method with a proved bounded relative error), and after #158 the kernel's share FALLS with N because Davies never sees N -- so ALL further p-value work is capped at <=1.6x at N=50k and ~1.1x at N=400k. fastSKAT/bigQF is a measured LOSS at these q and implements SKAT not SKAT-O; SMMAT-E deletes the rho grid but changes the test. Two things it found that are NOT speed: the per-rho Davies trust floor is ~1e-6 (SKAT hardcodes `acc = 10^(-6)`), NOT machine epsilon; and the landed saddlepoint tail is UNVALIDATED at q=4-10, the band carrying 56% of Davies time, where Kuonen's o(n^-3/2) bound runs the wrong way and no published validation exists below n=305 -- a missing control, not a known defect. **REFUTED 2026-09-19 (`docs/state/handoff-2026-09-19-carrier-phi-refuted.md`): the layer speedup is REAL and huge (56-357x at realistic q, benchmarked in C++ against OpenBLAS with a correctness control at 1e-15 and a passing negative control) but `varmat` -- the only bucket it touches -- is 3.365 s of 96.69 s, so the RUN prize is 1.034x at N=50k and 1.062x at 400k. The note attributed the whole 54% C++ share to a lever that touches 3.5%. Two code facts it also got wrong: there is no per-cell Phi construction in C++ (Phi is an R subset of a region-level VarMat, `R/SAIGE_SPATest_Region.R:511-516`), and Sigma^-1 is NOT diagonal for 36% of rows (3282 of 9058 take `getPCG1ofSigmaAndGtilde`, and those are the RAREST markers -- exactly where carrier space would pay). DO NOT BUILD IT.** Its top candidate (carrier-space Phi by the four-term expansion) is a FLOP COUNT, and LEDGER #156 is the located precedent that a flop count in this exact layer was already wrong by 20x (24x predicted, 1.15-1.57x measured) because it ignored BLAS level -- so bench it synthetically before building it.

**Current handoff: `docs/state/handoff-2026-09-16-pgen-sparse-read.md`,
its "State at the end of the session" section first** -- step 2 reads PGEN
records in carrier space on BOTH paths, and both are MERGED on `strip`
with CI green: the marker path (`baa9ab24`: `GenoReader::getOneMarkerSparse`,
pgenlib's difference list, cap `2 * (N/8)`, prepared over the listed samples
with the dense path's arithmetic; the zero-genotype index is no longer
materialised for either reader; sim50k 1 thread quantitative **1.87x**,
binary **1.21x**, and on a pgen FILTERED to MAF < 1% **10.4x** quant /
**6.1x** binary Firth off / **2.8x** Firth on) and the region path
(`c05f23fb`: `prepareAndTest`, BURDEN **1.56x** / SKAT-O **1.13x**,
`test-step2-region-sparse.R`; the PGEN dense loop's discarded `push_back`
went with it). Every golden byte-identical at each step. What was ranked
after them, and where each stands (decided 2026-09-17):
(a) region path -- DONE (above). (c) PGEN dense loop -- DONE (above).
(d) **carrier-space saddlepoint inputs on the marker path -- BUILT, GATED,
DEFERRED, NOT merged.** Commit `37a4d5b6` on branch `sparse-read` in the
`SAIGE-sparse` worktree (its gate: FAIL 0 / PASS 5799 at `both`, no golden
moved, `test-step2-spa-carrier.R` agrees with the dense arm to 1e-10; kill
switch `SAIGE_SPA_DENSE=1`). Not byte-neutral by construction (reassociated
sums), adds a p x N precompute per model (~70 MB at N = 424k, p = 20), and
measured **1.11x** on the rare-only pgen, noise on the full spectrum -- the
common markers take the dense saddlepoint it cannot touch. The one open
question is the gain at real N, where the O(N) terms it removes are 8x
larger; **the decision rule is a PINNED All of Us chr20 pair, carrier vs
`SAIGE_SPA_DENSE=1`, after the `submit()` pinning bug is fixed
(`handoff-2026-09-16-aou-first-results.md`, second correction): under ~10%
there, drop the branch; above, merge and check whether any row moved.**
(b) **the bed reader byte-scan -- MEASURED AND DROPPED 2026-09-17**
(`handoff-2026-09-17-bed-sparse-read.md`; the implementation, 9-test drill
and green gate are archived on branch `bed-sparse`, `a476e478`, NOT on
`strip`). Built because collaborators will not all convert to PGEN;
measured on sim50k: bed became LEVEL with the PGEN read (36 vs 35.5 user-s,
quant full spectrum, byte-identical) but only 1.2x over dense bed, 1.1x
rare-only quant, 1.0x binary either way -- because bed's lane-tally
pre-screen (Item 3, 2026-09-05) already skipped the decode of every marker
the MAC gate drops, which dense PGEN never had, so PGEN's 10x was mostly
recovering what bed already had. The owner dropped it as not worth a second
code path. **Do not re-open without a measurement that contradicts that
table.**
Not levers: dosage PGEN (unsupported by design, dense), step 1. The
remaining O(N) terms on the region path, in order of size: the R-side Davies
integration (~45% of a SKAT-O run), then the per-marker `gtilde` for P1/P2.

**SKAT-O p-value methods, MEASURED 2026-09-17 on both traits (`docs/state/handoff-2026-09-17-skato-pvalue-methods.md`, harness `harness/skato_pvalue/`).** The sparse read already feeds SKAT-O; what SKAT-O adds over BURDEN in C++ is ~13% of the run (P1/P2 copies, the per-cell residualised vector, the VarMat product) and 414 MB of RSS at N = 50k; the R-side Davies is ~70%. Modified Liu as a SCREEN at 0.01 saves 96% of Davies time with zero misses on either trait (~2.5x); as a swap it is 3-7x anticonservative below 1e-6. SKAT as shipped is off 0.05x-7.5x in (1e-8, 1e-6] on the binary trait (per-rho Davies at acc 1e-6); the saddlepoint in both layers is 5-12% conservative and stable everywhere. Recommendation: Liu screen + saddlepoint confirmation below 1e-5 (REGENIE's tiering), numeric gate, sub-1e-6 goldens should move. NOT implemented; the owner confirmed 2026-09-18 that the SKAT-O p IS consumed (not with FlexRV), and chose option C. **BOTH STEPS LANDED 2026-09-18 (`R/SkatoScreen.R`): the screen (`--SKATO_screen_threshold` default 0.05; 17 region goldens re-baselined, `Pvalue` only, 0.985-1.001; 23 sidecar goldens; P67 pins the upstream path) and the saddlepoint tail (`--SKATO_tail_method=saddlepoint` below `--SKATO_tail_threshold` 1e-5; no golden moves). #158 is FIXED. Implementation notes are handoff section 10** -- the screen must be OUTER-INTEGRAL ONLY (the per-rho layer feeds `Pvalue_SKAT`, `Pvalue_Burden` and the meta-analysis burden Z, and Liu's per-rho burden value is 0.13x-4.6x Davies'); the outer integral is 95% of the kernel so the saving is unchanged (96-98% at 0.01, zero misses on both traits); blast radius is `Pvalue` alone, by under 10%, on screened cells, so the gate is a stated comparison rule, not the byte comparator. Owner's additions 2026-09-18: threshold 0.05 by default (multi-biobank meta of SKAT-O p's; Liu unbiased above 0.05), two CLI options `--SKATO_screen_threshold` / `--SKATO_tail_method` whose (0, davies) setting must reproduce the goldens byte-for-byte, and a `.skatoMethod.txt` sidecar for provenance instead of a column. **Null replicate study 2026-09-18 (handoff section 9): 32 replicates x {quantitative, binary 500 cases, binary 100 cases}, ~100k cells per arm -- screen misses ZERO in every arm; SKAT/Liu/saddlepoint vs reference on null calls 0.99-1.01 / 0.84-1.06 / 0.97-1.09. Found on the way, shared by every option and by today's code: at 100 cases the SKAT component is inflated 2x at 1e-3 and the SKAT-O omnibus 2.7x at 1e-4 (27 cells vs 10), in tau2 = 0 and tau2 > 0 replicates alike; 500 cases 1.5x at 1e-3 only; quantitative calibrated. That is the kernel's Phi adjustment at few cases, not the inversion: **LEDGER #159, DECIDED 2026-09-18 (owner): leave as is, noted**; check the few-case All of Us phenotypes against its numbers before quoting SKAT-O tails there.** plink2 has no SKAT-O; upstream SAIGE calls the identical SKAT kernel, so comparing to it is a null control.

**Docs (2026-09-15): the pkgdown site and the `man/` guard landed** (`harness/docs.sh`,
`.github/workflows/docs.yml`, `_pkgdown.yml`, `vignettes/articles/cli_options.Rmd`
from `optiondrift.py --json`). Open, owner's call: (a) the ~15 `R/` files with no
roxygen at all (`R/Util.R`, `R/readInGLMM.R`, `R/MarkerScan.R`, `R/FlexRV.R`, ...)
-- the site documents only the six topics `man/` has; (b) `NAMESPACE`'s
`exportPattern` exports all 163 Rcpp wrappers undocumented; an explicit export
list would hide them but changes what `SAIGE::` resolves in the tests, so it is
a gated change, not a docs one; (c) `extractNglmm.R`'s 40 options are tabulated
but their live/inert status is not hand-traced in `docs/cli-options-and-group-file.md`.

**Current handoff: `docs/state/handoff-2026-09-15-optiondrift.md`** -- LEDGER
#151's guard exists: `harness/tools/optiondrift.py` compares every CLI
wrapper default against the exported formal by VALUE, runs before every
tier, and fails on any divergence not on its ALLOWED list (keyed to a ledger
row, with both values, so an entry goes stale loudly). Its first run found
eight divergences beyond the four the hand sweep had; **all twelve were
decided by the owner the same day (#151 DECIDED, ALLOWED list empty)**, and
**#153 is FIXED**: the 1e16 Nglmm was a non-converged PCG on the dense GRM
because the sparse file was given without the flag; `getNeff` now requires
the sparse GRM and always solves on it (402.4175 = the closed form). **#146(b) FIXED** the same day (an all-zero weight matrix composes its sets and
falls under the empty-cell rule, the owner's framing: same as the annotation
not being present). **#154 drilled**: 24 concurrent 4-thread P3 runs at load
12-23 all clean and byte-identical, TSAN region races=0 at t=1/4/8 -- stays
as the re-run rule. **#150 FIXED** (owner: delete): the four dead chrX/ER assignments are gone,
each replaced by a comment naming the row. **Open: #154 only**, the re-run
rule, now with its mechanism candidate on the row (oneTBB workers alive at
exit, documented; 0/232 reproductions). **#155 found by that drill and FIXED**:
the SKAT-O pool's cluster port could collide between concurrent processes;
it retries. The ledger's execution queue is empty.

**Previous handoff: `docs/state/handoff-2026-09-14-phase-g.md`** -- Phase G
DONE as the audit prescribed: G1's five slices of `SAIGE.Region` (new
`R/RegionState.R`: `regionStateFor`, `regionTestPlan`, `regionDoneMessage`,
`composeWeightSets`, and the accumulator environment `scan`, no `<<-` left)
and G2's verbatim hoist of `mainRegionInCPP`'s twelve function-local types
into `src/RegionLoop.hpp`. One gate per slice, no golden moved, region TSAN
clean + plant after G2. The ledger's open rows are all decisions (#146b,
#150, #151, #153, #154). Nothing from the audit's execution queue remains.

**Previous handoff: `docs/state/handoff-2026-09-14-ledger-143.md`** -- LEDGER
#143 is FIXED in two commits (`5564c421` the units, the next the kernel). The
binary conditional region path fed its scale-factor kernel an unweighted
score against weighted inputs and ran a pre-#104 twin; the owner's decision
the row asked for was settled by a control instead: conditioning cannot depend
on the weight given to a conditioning marker, it did (0.3-6 % on P3's
`_cond` p-values), and it does not now (1e-15). Both `r` predictions landed to
the digit (0.4733 -> 0.951743 -> 1.000000). Five goldens moved on their
`_cond` columns only; everything else byte-identical, as predicted. One
finding parked, not a row: `SPA_pval`'s absolute root tolerance makes the
region burden SPA scale-dependent at ~1e-4 in log p. **Later the same
evening: #147 and #148 FIXED** (`cab75ac5`, `6f333786`; one gate each, no
golden moved, predicted): `--minGroupMAC_in_BurdenTest` is validated by
name with `minVal = 0`, and quantitative `getNeff` with `--noEstFixedEff=TRUE`
builds its offset the way the binary arm does (decidable by measurement, not
the owner: `Nglmm` reads `fit0` only for its length). **LEDGER #153 OPENED**
from that drill: `getNeff` prints `Nglmm` ~1e16 for N = 1,000 on both traits
because `tauVec_ss = c(0, 1)` makes the solved system the GRM alone; upstream
verbatim, owner's decision on the intended `tau`. #150 is untouched: its code
stays by decision and the surfaces were already corrected. **Later still
(night): #145, #144, #146(a), #152 and #149 FIXED** (`0c2239d8`, `742daa41`,
`5aa494b5`, `adfa9990`, `4f2fbd07`; one gate each, no golden moved, predicted
each time), and
**LEDGER #154 OPENED** (a teardown segfault after `Analysis done!` on a
4-thread SKAT-O run under load 13-21; 0/8 in isolation, outputs byte-identical;
recorded, re-run gate counts). The ledger is at `ledger.py open` = the
decisions only: #146(b), #150, #151, #153, #154. Warning baseline 92.

**Previous handoff: `docs/state/handoff-2026-09-14-flexrv-p2.md`** -- the
FlexRV power fixture p2 (phenotypes only, on the 50k genotypes; causal signal
at resolvable MAC, the score informative by construction) and its 72-run
sweep at collapse 10, predicted first (`prediction-flexrv-p2-2026-09-14.md`).
**The 2026-09-13 open question is closed: on `y_quant_p2` the permuted score
gives the unweighted burden gene for gene (29 : 31, sign p 0.90) while the
informative score finds 47 of 60 genes against unw's 36 (47 : 13, p 1e-5),
monotone in score quality, zero false detections among 150 null genes.**
Two owner decisions recorded there: the collapse cutoff stays at 10 for both
traits (upstream's default), and this fixture replaces the 2026-09-10 one for
power questions. **Both loose ends were closed the same evening (section 8 of the handoff):
the binary arm rebuilt on the exact liability-threshold model (per-gene
scale; unw 38 of 60 where 22 had been, informative 43, permuted 40 = unw gene
for gene, Spearman(pred, obs) 0.76), and a five-replicate null series that
puts the pooled Cauchy rows at 0.056 (flex192) and 0.062 (grid) at alpha 0.05,
inside the interval, so the single `y_null`'s 0.09 was a high draw.**
Nothing is open on this fixture; the next FlexRV question is how the score
gain scales with cor(score, causal) and N (`make_power_pheno.py --k --pmax`).
The 2026-09-12 decisions memo's fixture and collapse-0 halves are superseded.
The audit handoff below is still the code-side pointer.

**Previous handoff: `docs/state/handoff-2026-09-14-codebase-audit.md`** --
a whole-codebase structural quality audit at `9deab7c9`. 96 findings raised,
each faced by two independent skeptics (one attacking the facts, one attacking
bug-for-bug preservation and golden movement); **62 survived, 34 refuted**, all
verifier corrections applied. It is NOT a bug hunt: where it found a latent
defect the deliverable is a LEDGER row. **Read it end to end before acting on
any section** -- its section 8 is a do-not-re-raise list, and five refutations
were refused because this repo had already decided them.

**The audit's five structural facts**, which is what the 62 findings actually
are:
1. **One rule, written N times, fixed once** -- the branch's most expensive
   recurring defect, the documented cause of #44->#66, #118->#134->#135, #91
   and #37, and live-unrepaired today in the binary conditional region kernel
   (`R/SAIGE_SPATest_Region_Func.R:433-448`, a textual twin of the #104-repaired
   function 18 lines below it).
2. **The R<->C++ boundary has no contract** -- `src/Main.hpp` WAS 303 lines with
   one includer and eight declarations that disagreed with their definitions;
   a differing parameter list is a new overload, so the drift was
   structurally uncatchable. **Deleted 2026-09-14 (Phase F)**: the class is
   now impossible, `RcppExports.cpp` is the one checked copy.
3. **~3,000 lines of dead surface**, most already inventoried by this repo's own
   2026-09-08 audit and parked -- plus a slice that audit's method could not see
   (it scanned functions, never parameters, never post-`stop()` tails, never a
   translation unit absent from `OBJECTS`).
4. **Two god functions larger than upstream's** after this branch extracted 675
   lines from one of them. Every new feature has landed inside those frames.
5. **Three controls that do not control** -- and this is why 1-4 persist.

**Nothing from the audit has been implemented.** No code changed, no row
written, no gate run. The execution queue is in the handoff; the first two
phases need no gate at all:
- **Phase A -- DONE 2026-09-14** (harness only, no gate): `harness/tsan.sh` no
  longer aborts after the first arm that dies under `plant` (drilled both ways
  with `SAIGE_TSAN_TIMEOUT=5`: unpatched one arm / exit 1, patched three arms /
  exit 0, and `clean` still three arms / exit 1); `harness/prepush.sh` logs
  `shape=serial|concurrent` from the branch that ran, warns when
  `SAIGE_PREPUSH_CONCURRENT=1` arrives with one mode, and the dirty-tree
  report is one function. Before this, every plant result on this branch
  ended at the first timed-out or crashed arm with no t=4/t=8 run.
- **Phase B -- DONE 2026-09-14**: LEDGER #143-#151 written (`ledger.py open`).
  #143 was INSTRUMENTED first: on golden P3 the conditional kernel's
  `r = 0.4733`, not 1, so the pre-#104 twin is a LIVE divergence on the
  `_cond` columns, not a coverage statement (numbers on the row). #147, #148
  and #150 were drilled (non-numeric gate dies unnamed in Rcpp; quantitative
  getNeff with `--noEstFixedEff=TRUE` dies on `covoffset`; `--minMAF=0.02`
  with the chrX flag tests 20 of 48, so the "replaces minMAF with 1/(2N)"
  sentence of #139 was false and its six surfaces are corrected -- the ONE
  `R/` change in that commit is that stdout sentence). Nothing fixed.

**The ledger went from 0 open to 9 open on 2026-09-14** (Phase B, its own
commit). Every open row is a decision or a fix for the owner; none was made in
passing.

**Phase C (the dead-surface sweep) LANDED 2026-09-14**, five commits (C2 and C3
share `Main.cpp` and went in together), one clean build, one both-modes gate,
no golden moved (predicted before the gate). Deleted: `Binary_ComputeExactMC.*`
(the uncompiled orphan), `Binary_ComputeExactSKATO.cpp`, `Binary_Permu_SKAT.*`,
`SKATO_Exact`/`SKAT_Permu`/`GetProb` and the `SL_*` leftovers (C1);
`Unified_getMarkerPval` and its steer-away comments (C2); the 37 function-scope
declarations of `mainRegionInCPP` incl. the two indeterminate-bool reads, the
frame-scope `hasVarRatio`, and the three unread writer parameters (C3);
`prepSPAContext`/`spaContextReady`, `gNA`/`muNA` from the three internal SPA
kernels (C4); `getMem.*` and `UTIL::imputeGeno` with the shifted anchors
repaired (C5); `R/Firth.R`, `ReadGroupFile`, `check_close`,
`Get_Variance_Ratio`'s post-`stop()` tail and two unread formals, the
`resultDF` return field, the `printComb(3)` call (C6). `OBJECTS` is 17.
**Not taken, by decision or scope:** the three `[[Rcpp::export]]`ed
`Unified_*`/`closeGenoFile` (public API under `exportPattern`, the #114 class --
owner's call); the `t_genoType` parameter removal; the `annoMAFIndicatorMat`
collapse in `region_postprocess` (it would ADD an assignment on the
`is_rerun && !is_single_in_groupTest` path, so it is not provably inert); the
"may ride along" R functions (`pcg`, `bigGRMPar*`, `writeOutputFile`).

**Phase D (turn the compiler back on) LANDED 2026-09-14**: `-w` out,
`-Wall -Wextra -Wno-unused-parameter -Wno-cast-function-type-mismatch` in,
`harness/build.sh` ratchets the src-anchored count against
`harness/warnings-baseline.txt`, CI promotes `uninitialized`/`return-type` to errors
via `SAIGE_EXTRA_CXXFLAGS`; `parentheses` was tried and dropped the same day
(run 34836628887: GCC fires it inside `boost/mpl/assert.hpp`, unfixable here,
and at `SAIGE_fitGLMM_fast.cpp:1179`, a `|`-around-comparisons site clang does
not flag -- note for a future row or a same-parse parenthesisation). The one source change is parentheses making
the existing parse of `m_is_Firth_beta & (pval <= m_pCutoffforFirth)` explicit.
**Phase E (the deduplications) LANDED 2026-09-14**, six commits, one gate
each, no golden moved in any (predicted before each gate); the region TSAN
clean and plant arms were run ONCE over the whole phase after E6 (E1, E4 and
E5 all touch the region prepares' frames). E1 `contextForMAC` +
`isSingleVarianceRatio` + `isERCandidate`, `hasVarRatio` demoted to a local,
and a `null`-less ratio file refused by name; E2 one `SingleVariantSchema`,
one header writer, one row writer over `MarkerResult`; E3 `formatPvalOrLog` /
`formatLogP`; E4 `struct WeightSet` and ONE emission loop (byte-identity on
`x * 1.0 == x` and a per-kind SELECTED count source); E5 `CellResult` by
inheritance and one `writeBackCell`; E6 `checkGenoInput` returns its resolved
paths, the covariate GLM is fitted only where read, the region `.index`
writer is `writeOutputFileIndex`. **Not taken, by decision:** `finishScoreTest`
(the three guards are not one predicate), the `ReadModel` merge (value-pinning
tests for the subset-only branches must come first), `getNeff`'s LOCO block
(LEDGER #149 is the owner's row), `Marker::MarkerResult` into `CellResult`
(different defaults the writer keys on).

**Phase F LANDED 2026-09-14**: `src/Main.hpp` deleted; the four internal
forward declarations (`processMale_XnonPAR`, the three `writeOutfile_*`) are at
the top of `src/Main.cpp`, copied from their definitions; `MarkerLoop.hpp` is
included directly; LEDGER #5's STATUS appended (not re-opened); a `docdrift.py`
rule keeps that stale-declarations sentence from being copied back (docdrift-allow). Only
Phase G (the god functions) remains, and it is explicitly deferred.

**LEDGER #152** records the four enumeration-guarded uninitialised reads
(two live, two in the dead exports) and, since CI run 34838510204, the six
more GCC's `-Wmaybe-uninitialized` names in `SAIGE_test.cpp` -- read them one
by one before deciding anything; CI demotes that sub-category to a warning. It adds a gate and deletes no
complexity. Next: Phase E (the deduplications), one push, several commits,
E1 additionally needs the region TSAN clean + plant arms.

**The "also outstanding" doc items in the handoff's preamble were already
closed** when checked on 2026-09-14: `R/SAIGE_Test_main.R:114` and
`man/SPAGMMATtest.Rd:210` carry the conditional "192 by default, fewer where a
mac rung ..." wording (`git log -S` puts the fix inside `9deab7c9` itself, so
the handoff's "verified still present at 9deab7c9" was written against an
earlier tree), and `docs/design/flexrv-design.md:607`'s "192-set run" describes
the external-weight-line test route, which deliberately keeps all 192 lines
(#142's STATUS). No docdrift rule was added: the sentence has no live copy and
nothing at `4c2ca05a^` matches it, so a rule would drill against nothing.


**Previous handoff: `docs/state/handoff-2026-09-13-night-135-139-142.md`** --
#135 and #139 fixed (goldens P66, P67/P68), the ONE open row #142 laid out
as a decision with its evidence and a recommendation, the corrected scratch
install recipe, and the pre-push sequencing trap. The evening handoff
(`handoff-2026-09-13-evening-six-rows.md`) holds the six earlier rows and the
`pool_rules.R` first run.

**2026-09-13 (evening): LEDGER #141, #140 and #116 are FIXED**, one commit each,
gated together at `both` modes. The ledger is at **6 open**. Read the rows
(`ledger.py 141`, `140`, `116`); the short version:
- **#141**: a repeated IID in the fam/psam is refused at every step -- steps
  0/1/getNeff in `readGenoSampleIDs` (`R/Geno.R`, both readers, #67's message
  shape), step 2 in `readFamFile` (`src/PLINK.cpp`; PGEN already refused). Not
  in the three-string `PgenShared` overload the row proposed: it never keys on
  IDs, and every ID-keyed path reads them through R first.
- **#140**: S9 carried to the region path. `flushOutfile_regionInCPP()`, an
  index line carrying BOTH files' sizes (`... chunk N single_bytes S bytes B`,
  writer `indexChunkLine()` beside reader `parseIndexChunkLine()` in
  `R/Util.R`), and restart truncation of the results file and
  `.singleAssoc.txt`. One pre-existing test pinned the legacy line verbatim
  (`test-step2-multiweight-behaviour.R`) and was loosened in the same commit.
  Still open on the row: the two end-of-run sidecars are rewritten over the
  resumed regions only.
- **#116, diagnosed then decided**: the mechanism is the INGEST QC gate (`rs1`
  and `rs2` are monomorphic in the analysed 1,000; ALT count 0 measured with
  `plink2 --freq --keep`), and the early return in `regionFinish` also skipped
  the index bookkeeping (a one-region run recorded `chunk 0` and no completion
  line). The decision is the BURDEN writer's rule on both paths: a header at
  start and one `Region Cauchy NA ...` row per region with nothing tested,
  disclosed in the #127 sidecar. `R/RegionSchema.R` states the SKAT-O schema
  once and `regionAssertSchema` stops the run if any chunk of real rows
  disagrees with it, so every committed golden gates the header. Moved no
  golden (no committed golden is 0 bytes or carries an all-NA Cauchy row; the
  row's "would move every golden with a header" was a mis-stated blast radius,
  #133's shape). `test-step2-region-allrefused.R`'s two SKAT-O tests were
  inverted with it.
**Later the same evening: #137, #138 and #136 FIXED too** (`49274ad1`,
`742726b3`, `51a8e7f2`), the ledger at **3 open** (#135, #139, #142). #137 is
the missing `--varianceRatioFile` refusal beside the manifest; #138 REFUSES any
set composed beside FlexRV's transform sets (weight line, extra beta pairs,
No_Weights), placed at set composition because weight lines are known only
once the group file is read; #136 corrects six help sentences in every copy
(a seventh twin of the LOCO claim was in `R/SAIGE_extractNeff.R`), DECIDES to
delete the `isCateVarianceRatio` claim rather than enforce it, and adds six
`docdrift.py` rules -- the drill at `4c2ca05a^` now reports 42 hits, not 10.
`man/` was regenerated with `roxygen2::roxygenise(roclets = "rd")`, which also
re-wrapped one unrelated paragraph (whitespace only). Next: #142 is a DECISION
(dedupe FlexRV MAC sets before the Cauchy row, or document that `p_FlexRV`
depends on N through the duplicate count); #139 needs a chrX fixture with
missing calls; #135 needs a fixture whose cell count approaches 2N.

**2026-09-13 (night): LEDGER #135 FIXED.** The BURDEN emitter gates on the
cell's aggregate count itself (`MAC = altCounts`); the fold and the clamp are
gone. The fixture needed no new group file: `plink2 --freq --keep` over the
model's 1,000 samples showed every P65 cell but GENE1 lof past N before the
collapse (plain folded sums 959, 5271, 3043, 5789; printed after the collapse
at 700: 748, 3064, 1915, 2997). Golden P66 is P65 at
`--minGroupMAC_in_BurdenTest=1500`: at unpatched HEAD the four cells gated as
748, 1000, 85 and 1000, ALL were refused and the file was two `Cauchy NA` rows
while the SKAT-O path refused only GENE1 lof; fixed, five of P65's six rows
are byte-identical, GENE1 lof is refused and its pooled row is the CCT of the
two twins (0.11969). `test-step2-burden-aggregate-mac.R` is the
characterization (red at the scratch install of HEAD, green after). No
committed golden moved. Read `ledger.py 135`.

**2026-09-13 (night): LEDGER #139 FIXED.** `applyMaleXnonPAR` recounts the
doubled males over the non-missing entries only (it now takes the reader's
missing index). The fixture is `genotype_100markers_missing` with the
bim/pvar chromosome rewritten to X (`extdata/input/..._missing_chrX.{bim,pvar}`)
plus the bundled `sampleid_males.txt`: at unpatched HEAD the PLINK run tested
6 of P42's 48 markers (the six without a missing call) and the PGEN run died
on `MAC nan`; fixed, goldens P67/P68 are byte-identical across readers with
all 48. `test-step2-xnonpar-missing.R` (red at the scratch install, 9
failures) also covers `--condition` on a chrX marker with missing calls. The
false chrX `cat()` is corrected and a `docdrift.py` rule pins it (the
`4c2ca05a^` drill is 43 hits now). No committed golden moved. Read
`ledger.py 139`.

**2026-09-13 (late night): LEDGER #142 DECIDED (owner: option (b)) AND FIXED.**
`flexrvTransformTable(..., N = )` drops every `mac` rung at or above
`2 * N * maxMAF` before the sets are composed and the run warns with N, the
bound and the rungs; `SAIGE.Region` passes `length(mu[[1]])`. Goldens P58,
P59 and P61 moved exactly as predicted from their own rows (193 -> 113 rows
per gene; on P58 GENE1's pooled row 0.139077 -> 0.0956982, GENE2's 0.907609
-> 0.93903), P60 stayed byte-identical (no `mac` rung in its spec).
`test-step2-flexrv-rungs.R` is the characterization the request specified
(80 dropped on the fixture, 0 with `--flexRV_maxMAF` raised past
500 / (2N)); `test-step2-flexrv.R`'s row-for-row duplicate test became the
absence assertion. The ledger has NO open row. Read `ledger.py 142` and the
night handoff for the evidence the decision rested on.

Three things a later session should hold on to: (1) a scratch install of
unpatched HEAD reusing the tree's `src/*.o` (all but the changed file's) takes
about a minute and is what every red drill above ran against -- assert
`system.file(package = "SAIGE")` in the drill's output, and check the new
guard's string is ABSENT from the scratch `.so` and PRESENT in the patched one
before believing either colour; (2) `fitNULLGLMM` pre-creates the model `.rda`
EMPTY before it reads anything (`file.create`, `R/SAIGE_fitGLMM_fast.R:851`),
so a refusal test must assert the file's SIZE, not its existence; (3) the
comment linter parses a comment line whose first token is a bare symbol
(`LEDGER #140 FIXED. ...`) as code -- start such lines with a phrase.


**Current handoff: `docs/state/handoff-2026-09-13-132-133-134.md`** -- LEDGER
**#134, #132 and #133 are FIXED**, one commit each, each gated at `both` modes
(`a2962701`, `58785e12`, `8bb5a510`). The ledger is at **9 open**, from 13 that
morning. Read that handoff, then this file. The three prediction documents beside
it (`prediction-134.md`, `prediction-132.md`, `prediction-133.md`) were each
written BEFORE their gate and had the outcome appended, including where the
prediction was falsified.

Two things from it that change what a later session should believe:
**#133's WHAT and BOTH its candidate fixes were mis-located** (candidate A
returns NEGATIVE values on a mixed cell; candidate B's blast radius was
overstated), and **#134's row explained its own doubling measurement with the
wrong reason** -- the band boundaries ARE crossed at x2, they just select the
same variance ratio, because only 20.5 is live. Both corrections are on the rows.
**#132's "normalise before ACAT-V" candidate is retired as INERT** -- CCT already
normalises internally.

Also landed: **folded ultra-rare membership is now a PINNED GOLDEN** -- six
committed `.markerList.txt` sidecars, which #118 and #134 both wanted and never
had. And the previous entry's earlier handoff pointer is superseded:
`handoff-2026-09-13-flexrv-sweep-rerun.md` is still correct about the sweep and
is the one to read for the FlexRV study state.



**2026-09-13 (later): LEDGER #128 is FIXED, crash-only, and #116 stays open on
purpose.** A SKAT-O run whose every cell is refused by
`--minGroupMAC_in_BurdenTest` used to kill the whole run one region later; the
degenerate branch now writes no pooled row, which is what that path already did
for a one-cell region. Moves no golden. Test inverted in the same commit and
drilled both ways (4 failures against a scratch install of the unpatched tree, at
exactly the crash assertions; 19 passes with the fix; the two control tests pass
in BOTH). **The owner's decision was to fix the crash and NOT the reporting
question**, because how an empty region is reported is an output-schema change
that deserves its own decision and its own goldens. So the live inconsistency is
now exactly #116 and nothing else: the BURDEN path writes an `NA` Cauchy row for
the region the SKAT-O path omits entirely. The third test in
`test-step2-region-allrefused.R` pins that disagreement so a later harmonisation
cannot land unnoticed.

**2026-09-13 DECISION: the FlexRV collapse cutoff stays at 10, and 0 is not to be
recommended for a BINARY trait until a replicate null sweep says it is safe.**
The paper's unbinned run is available at 0 and golden P59 pins it, but
SAIGE-GENE+ introduced the collapse FOR CALIBRATION, this repo has never measured
the regime where that matters, and the bundled fixture and the 50k sim study both
run with it ON -- so every calibration result here, null and power, describes
cutoff 10. Turning it off is not a neutral "closer to the paper" choice: it
removes the mechanism the method's calibration was built on, in a regime nobody
has checked. The gate on changing it is a replicate NULL sweep at cutoff 0 on a
binary trait showing type-I error holds. Recorded in CLAUDE.md's FlexRV section.

**This CONSTRAINS a live recommendation**: `docs/state/flexrv-study-decisions-2026-09-12.md`
recommends "stratified + class-only + flex24 now at collapse 0", and its arms
include `y_binary`. That memo is dated and is not edited; read it with this
constraint applied -- the collapse-0 half of that recommendation is
quantitative-only until the null sweep exists. `handoff-2026-09-09-flexrv-f0.md`
already scoped its own "cutoff 0" recommendation to quantitative and needs no
qualification.

**Current handoff: `docs/state/handoff-2026-09-13-flexrv-sweep-rerun.md`** -- the
FlexRV sweep is RE-RUN and no longer owed. All 72 collapse-10 runs at
`head=976e73ce dirty=0`, exit 0; the 72 collapse-0 runs were left alone because
their logs prove 0 ultra-rare markers. Outcome in one line: **12.0% of rows moved
(max |dlog10 p| 5.77), the null calibration did not move at all, and two CONTROL 3
lines went from `monotone=FALSE` to TRUE.** The movement is confined to
fractional-weight sets and is monotone in weight shrinkage, with `unw` and `lof`
exactly zero over 104,760 rows -- the mechanism is #129's dosage filter against
the 0.2 cutoff. The pre-fix outputs are kept at
`harness/work/flexrv_study/results-prefix-38eee5f5/`. **Two things that outlive
it**: the power arm is at its own oracle ceiling, and `uninf<=unw=FALSE` on every
line -- i.e. at this fixture's resolution FlexRV's advantage over the unweighted
burden survives PERMUTING the score, which is the open scientific question.
**Section 6 was REWRITTEN on 2026-09-13 and its "10-40x unexplained shortfall" is
withdrawn**: the MAF window costs 6 points of causal variance and the annotation
mask costs 64 (so a different `--flexRV_maxMAF` fixes nothing), 22 of 62 genes
cancel directionally, and 3.1 of the naive 4.9 expected detections came from 27
genes whose entire visible signal sits below MAC 20, where the chi-square
approximation is invalid. Against resolvable genes the prediction is 1.8 and the
observation is 1-2. **There is no anomaly.** The design constraint for any new
power fixture is 6e: control the causal variants' MAC, not just their h2.


**2026-09-13: the current handoff's closing recommendation is WRONG, and this is
the correction.** `handoff-2026-09-12-ur-pseudomarker-129-130.md` ends "where to
start next ... **#134** ... so fixing it before the sweep is re-run means running
the sweep once instead of twice". That reason is false. #134's decision points sit
inside `if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest)` at
`src/Main.cpp:3634`; every FlexRV arm of the sweep is `--r.corr=1` with
`is_single_in_groupTest` at its FALSE default, and the sweep's only SKAT-O arms
use beta weights, whose collapse is unweighted. **#134 cannot move one number in
the FlexRV sweep.** The handoff is dated and is not edited; the correction lives
here and in the #134 STATUS cell, which is authoritative.

Three consequences, in order:

1. **The FlexRV sweep re-run is UNBLOCKED and is the cheapest owed item.** It is
   genuinely still owed -- its outputs were generated at `38eee5f5`, before the
   #129/#130 fixes, and only the 72 collapse-10 runs are affected.
2. **Do #132 + #133 + #134 as ONE commit series over ONE gate cycle.** Five of
   #134's six measured movers are also #132's movers, and #133's option (b) moves
   `MAC` in every group golden. Landing them separately re-baselines the same
   files two or three times, and each re-baseline destroys the BEFORE state the
   next row's prediction must be read from -- which is this branch's own Trap
   ("predict which goldens should move BEFORE the gate, from the goldens
   themselves where possible").
3. **If one item must go alone, it is #128**, the only open row with its
   characterization test already written, moving no golden, on a configuration
   that currently DIES. It is blocked only on #116's decision -- and #116's
   mechanism is now diagnosed (below), so that block is liftable.

**#116's mechanism is diagnosed, and it is NEITHER candidate the row names.** On
the row's own bundled reproduction, `rs1` and `rs2` are both monomorphic in the
analysed 1,000 (AC 0 each; they are polymorphic only across the full 10,000 fam
rows), so at the default `--minMAC=0.5` both are rejected by the region ingest
gate at `src/Main.cpp:2751`, which `return`s before `pr.passQC = true`. Neither
marker ever reaches the ultra-rare collapse, so the min-group-MAC gate and the
max-MAF cell filters -- the row's two candidates -- never see anything. The tell
for anyone re-checking is the log line at `src/Main.cpp:3278` ("In chunks 0-N, i2
markers are ultra-rare and i1 markers are not ultra-rare"), guarded by
`if(i1InChunk != 0)`: its ABSENCE means nothing passed QC, its presence with
`i1 == 0, i2 > 0` would mean the collapse/gate candidates. That costs one grep of
an existing log, no instrumentation and no build. The consequence for the
decision: the same 0-byte file is now reachable by THREE unrelated causes (#76's
`--chrom` filter, #113(a)'s ID join, and all-markers-fail-QC), which argues for
the schema fix (always write a header) over a `--chrom`-shaped refusal.

**Also corrected 2026-09-13, each verified by hand at `f5e173a9` before being
written:** #134's blast radius is now MEASURED rather than predicted, and
`cond_customweight` is IN (6 of 6 UR-only cells cross a band; the row had said
"needs measuring") while NO FlexRV golden moves at all -- including
`plink_group_flexrv_skato`, which the first version of that correction wrongly
called undeterminable from the bytes: a transform set's `unw_unw` row has weights
identically 1, so its pseudo-marker `AC_Allele2` IS the cell's `MAC_nw`, and every
value in that golden (54 and 41 against 24.5-50.5) sits above the top band edge.
The direction is the part that generalises: FlexRV weights are at most 1, so they
bias a pseudo-marker's band DOWNWARD, opposite to CADD-scale custom weights; #134's ER decision
point is structurally inert, so it has TWO live sites, not three; #139's PGEN
failure writes `BETA=NA SE=NA` rather than `nan`, making it byte-shaped exactly
like a legitimate #44 declined-effect row; CLAUDE.md's claim that <!-- docdrift-allow: quotes the
falsified sentence on purpose -->#49 remains
(it was fixed 2026-09-05) and its claim that no bundled marker exceeds AF 0.5
(`rs26` in the dosage goldens is 0.500034) are fixed, each with a `docdrift.py`
rule -- and adding the second rule immediately found a SECOND copy of that claim
in the allele-convention section, which is the copied-twin shape the tool exists
for.

**LEDGER #142 opened 2026-09-13** from `docs/flag-surprises.md` F5, which had been
workflow-verified only: every FlexRV `mac:` threshold at or above `2*N*maxMAF` is
identically 1, i.e. a duplicate of `unw`, and the duplicates still enter the
pooled Cauchy row. Re-derived by hand from `plink_group_flexrv_burden.txt` before
the row was written: exactly **80 of 192** rows per gene are byte-identical to
their score transform's `_unw` row (precisely `MAC25`/`MAC50`/`MAC100`/`MAC200`/
`MAC500` x 16 score transforms), the CCT over all 192 reproduces the committed
pooled row to every printed digit, and over the value-distinct rows it is
**0.095698 against the golden's 0.139077 on GENE1 (1.45x) and 0.939638 against
0.907609 on GENE2 (0.966x)** -- so the distortion has no fixed sign. At the
paper's default `--flexRV_maxMAF=0.001` this is degenerate on any cohort below
N = 250,000. **Eight of the nine other leads in `docs/flag-surprises.md` remain
unharvested and have no row**; F12 is the next cheapest and is likewise already
confirmable from a committed golden.

**Structural gap worth budgeting for: eight of the thirteen open rows have no
characterization test** (#134, #135, #136, #137, #138, #139, #140, #141, #142).
The bug-for-bug rule requires the row's test to be INVERTED in the fix commit, so
each of those carries an unbudgeted "write it first, prove it red at unpatched
HEAD" step that none of the fix sizings on the rows includes.

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

**2026-09-12 also opened LEDGER #136 and #137**, from a fan-out inventory of all
139 CLI options that traced each from its `make_option()` to the code consuming
it (write-up: `docs/cli-options-and-group-file.md`, which also records which
options are refused, inert, trait-gated or path-gated, and the group-file format
with and without weights and with FlexRV). #136: SIX help strings and roxygen
lines state something the code does not do -- `--pval_cutoff_for_fastTest`
claims no effect on single-variant runs while `src/Main.cpp:1006` IS the (docdrift-allow: quotes the false sentence)
single-variant escalation gate; `--X_PARregion` documents `start:end` while the
parser splits on `-`; `--is_output_moreDetails` names four columns no writer
emits; `chr:pos_ref/alt` appears in six places where the branch matches bim (docdrift-allow: quotes the false sentence)
column 2; `isCateVarianceRatio`'s "If TRUE then LOCO=FALSE" is enforced nowhere; (docdrift-allow: quotes the false sentence)
and the `--IsOverwriteVarianceRatioFile` refusal is a `stop()` that opens
"WARNING" and describes work it will not do. Each was re-verified by hand before
the row was written. #137: the multi-model manifest's help promises to refuse
`--varianceRatioFile` and `resolveModelSpecs` refuses only `--GMMATmodelFile`,
so the flag is accepted and silently ignored. Both are docs-or-one-clause fixes
that move no golden; #136's fix must add a `docdrift.py` rule per sentence,
which is that tool's stated admission criterion.

**Twelve flag combinations that silently do the wrong thing are on ONE page**:
`docs/flag-surprises.md`, from a ten-class adversarial hunt (40 candidates, 14
verified). Three were hand-checked and have rows -- #139 (chrX male-doubling
plus a missing call: PLINK silently DROPS the marker, PGEN writes p = 1 beside
`nan` frequencies, from recounting alleles before imputation with two different
missing sentinels), #140 (a region restart can LOSE regions: the byte-count
truncation that fixed #47 is marker-path only and the region writers are never
flushed) and #141 (a duplicated IID in a `.fam` is accepted by all three steps
and step 2 then reads the wrong genotype row; PGEN's psam reader refuses it).
The other nine are LEADS, workflow-verified only -- confirm before acting.

**The five open region-path issues are summarised on ONE page**:
`docs/state/open-allele-and-weight-issues-2026-09-12.md` -- #131 (decided),
#132, #133, #134, #135, what each costs to ignore, which goldens each would
move, a suggested order, and the four unverified leads from the audit that was
stopped at four of six surfaces. Read that instead of reconstructing the list
from the rows. They share ONE root cause: SAIGE's single-variant idioms (the
fold, MAF, `2 - G`, the dosage filter) reused on quantities that are not one
variant's dosage vector -- which is also what #129 and #130 were.

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
* ~~**Still open, unscheduled:** if BURDEN at scale matters, thread across
  REGIONS rather than markers, which approaches the hardware ceiling
  (~2.9x at 4 threads here).~~ **MEASURED AND CLOSED 2026-09-19**
  (`docs/state/handoff-2026-09-19-region-threading-measured.md`). The ~2.9x had
  NO SOURCE -- it landed in `ddedfa6d` with no supporting number and is most
  likely `4 / 1.39`, an arithmetic bound assuming the serial residue goes to
  zero. Measured on sim50k BURDEN: serial 29.1 s, `--nThreads=4` 19.8 s
  (1.47x), and FOUR PROCESSES on quarters of the group file 15.7 s (**1.85x**),
  with the four quarters reproducing the whole-file 3,345 rows identically.
  **Those numbers were CPU-STARVED** (4 cores at load 3; the threads arm got 1.89 cores, the
  serial baseline 0.84). **Re-measured on an IDLE 10-core box, 7.6 free**
  (`docs/state/handoff-2026-09-19-region-parallelism-idle-box.md`): serial 11.06 s, threads4
  7.84 s (**1.41x**, only **1.54 cores busy**), procs4 5.90 s (**1.87x**, 3.18 cores busy),
  all three reps reconstructing the whole run exactly. The threads arm CANNOT use 4 cores --
  that is the ~45% serial fraction NOTES.md:1409-1419 models -- and four processes run four
  copies of it concurrently, which is exactly what this item originally proposed. So the
  original instinct was RIGHT. Nothing reaches 2.9x (best 1.87x). **The case for in-process
  region threading is EFFICIENCY: processes buy 1.87x at 1.85x CPU and 4x RSS (efficiency
  1.02), threads convert CPU far better (1.19) but cap at 1.41x; threading the region loop
  in process would aim at 1.87x with threads' efficiency.** Worth building if CPU-hours or
  memory bind; not for wall clock alone, since job-splitting already gives 1.87x. Note BURDEN
  never calls SKAT-O, so the `qfc` thread-safety blocker does not apply to this path.

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
`RcppExports`, the stale `Main.hpp` (deleted 2026-09-14, audit F), `src/Makevars` and the delete-`src/*.o`
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

**Amended 2026-09-14.** The method above could not see four things, and its
"all 17 are in `SAIGE_fitGLMM_fast.cpp`" conclusion was wrong when written: the
rule "exactly one occurrence of the name in `src/`" can never match an
`[[Rcpp::export]]` (its name also appears in `Main.hpp` and twice in
`RcppExports.cpp`), so `Unified_getSampleSizeinGeno`,
`Unified_getSampleSizeinAnalysis`, `Unified_getOneMarker` and `closeGenoFile`
(callerless since the fork import) belong in that bucket; and it scanned
functions only -- never parameters, never post-`stop()` tails, never
C++-internal symbols, never a translation unit absent from `OBJECTS`. The
2026-09-14 sweep (Phase C above) executed the disposition recorded here plus
that unseen slice; the wrapper-only exports and the `_LOCO` dead clones are
still parked.

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
