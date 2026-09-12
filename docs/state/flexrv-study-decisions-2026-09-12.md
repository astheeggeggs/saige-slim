<!-- Drafted 2026-09-12 by a three-perspective / merge / per-decision-refute workflow (10 agents, transcripts under the session's subagents/workflows/wf_d1d8a961-7ab), then vetted by hand: every file:line the memo leans on was re-read (make_scores.py:60-85, 02_build_dataset.py:497-503 and 596-607, flexrv-sim-study.md:308-312, SAIGE_Test_main.R:379-392, analyse.R:100-104) and says what the memo says. Claims the referees could not trace to a file are marked [ASSUMPTION]. Nothing here has been run; the options are proposals for the owner. -->

# FlexRV simulation study: decision memo (2026-09-12)

**Context.** The study (`harness/sim/flexrv_study/`, design `docs/design/flexrv-sim-study.md`) runs on the 50k chr20 fixture and is blocked on five owner decisions after three findings of the 2026-09-10 sweep: the uninformative-score arm beats the plain burden, the "synonymous" arm is not null on the signal traits, and `y_quant` detects almost nothing. The sweep was re-run at HEAD on 2026-09-12 (`results-2026-09-12.txt`; `docs/state/prediction-flexrv-rerun-2026-09-12.md`): it reproduces the shape (`y_binary` flex192 11 / 11 / 10 against `unw` 8 at collapse 0; `y_quant` 1-2 of 63; synonymous 7 / 6) and, as predicted, resolves none of the five. Bare line refs below are to `results-2026-09-10.txt`. LEDGER #129 (the weighted pseudo-marker zeroing) touches collapse 10 only: read the collapse-0 column.

## D1. The uninformative arm

The arm permutes only non-high scores (`make_scores.py:63-70`), so high stays 1 and the 14 nonlinear transforms remain soft-LoF masks; flex192(uninformative) 12 vs unw 8 (`:344,348`) is the expected sign, not a defect. It also shifts the in-cell missense marginal, so informative minus uninformative is confounded [ASSUMPTION: in-cell medium mean 0.75 vs 0.51, recomputed from `scores/`, untracked].

| Option | What changes | Expected effect | Cost | Risk |
|---|---|---|---|---|
| Stratified permutation + class-only arm | Permute within class, own RNG (`make_scores.py:69-70`), cmp-gated; a file with high=1, non-high=one constant | informative >= class-only >= stratified on paired per-gene log10 p; score's value = informative minus class-only | 2 python lines, 16 runs | Low |
| flex24 (annotation only) | `--flexRV_scoreTransforms=unw;lof` (`R/FlexRV.R:28`) | Below flex192(uninformative): the ratio ladder carries class information | 8 runs | Not a score isolator |
| lin (paper's baseline) | `sbeta:2` x `unw` MAF (`flexrv.txt:1400-1405`) | flex192/lin is the paper's ED Fig 6 statistic | 24 runs | Uninformative lin can beat unw |
| blind arm | Fully permuted score, own RNG | flex192(blind) <= flex24 plus slack: the pure floor | 16-32 runs | No real pipeline produces it |
| Oracle arm | e_lat percentile (`02_build_dataset.py:762-763`) | Ceiling on the score axis | 2 runs on y_binary | y_quant_pc needs a refit (D3) |

Refuted: flexscore_nolof (soft-lof copies, reads >= unw); annotation shuffle (changes the tested variants); within-class-sd=0 rebuild (medium still ranked among medium+low); drop lof in one arm (LoF still encoded); delete `ok2` (`analyse.R:102`, vacuous; replace it).

**Recommend** stratified + class-only + flex24 now at collapse 0, paired statistics only; lin and blind optional. Changes it: class-only within one gene of informative on the wide cell too (D5) means the fixture's score-effect link is the lever.

## D2. The low class

`low` is other_missense + synonymous (`flexrv-sim-study.md:309-311`) and carries most planted h2 [ASSUMPTION: 0.559, recomputed from `variants.tsv.gz`]; the synonymous arm masks it (`make_scores.py:76-85`, `run.sh:34`), so 7/6 on y_binary (`:349-350`) is a benign-missense power row. The fixture has no null class (the paper's check uses near-null real variants, `flexrv.txt:195-203`).

| Option | What changes | Expected effect | Cost | Risk |
|---|---|---|---|---|
| Truth-null mask arm | Low-annotated variants with beta == 0, scores = the existing draw subset (other files byte-identical, cmp); name it in `analyse.R:20`, `run.sh:34/45` | Pass; leakage NCP <= 0.07 per gene [ASSUMPTION: from `chr20_sim.bed`, untracked]; fail = k >= 8/63 at 0.05 or any Bonferroni detection | 16 runs, ~10 min (57.8 s flex192, 17.7 s unw, `study-2026-09-12.log`) | Truth defines the mask; label as calibration |
| Relabel | Docs only; report the arm as POWER; name the three null lines outside at both collapses (`:28/64`, `:105/141`, `:158/194`) | y_null_bin flex192 synonymous 9/210 at 0.01 (`:105`) read as a pooled-CCT calibration finding | 0 | Alone, no synonymous analogue |

Refuted: fourth synonymous class (the anno table is gone, md5s were in its header, the AoU AF file is not on disk); force low null (GENE0088 loses every causal variant, the forced draw at `02_build_dataset.py:497-503` shifts the main stream, y_null moves: 4 refits + 144 runs).

**Recommend** both, no rebuild; before any rebuild move the forced draw (`:502`) to its own RNG stream. Changes it: a fail with unw too is an engine LEDGER row; flex192-only excess points at the pooled row or the collapse.

## D3. y_quant

y_quant detects 0-3/63 with rho ~0 (`:297-332`) against 7-13 on y_binary (`:333-368`). The ladder scales the whole-gene oracle burden (`02_build_dataset.py:596-607`) while the test sees the unweighted high:medium MAF<=0.001 cell, and the binary NCP column is linearised (`02_build.log:5`), so y_binary detects through threshold convexity.

| Option | What changes | Expected effect | Cost | Risk |
|---|---|---|---|---|
| Re-specify the control | `analyse.R`: NCP of the tested cell, simulated expectation band instead of rho > 0.5 | y_quant unw E[det] 1.2 vs observed 1; y_binary E[det] 12.6 [11,15] vs observed 8, left OPEN [ASSUMPTION: refuter's Monte Carlo, untracked] | ~30 lines R | Independence approximation |
| y_quant_pc refit | Column in `pheno.txt` (h2 0.2); on-disk `.rda` dated 2026-09-01, stale; add the case to `step1.sh:44-51`, the trait to `analyse.R:74` | Written first: unw 8-10/63, cell oracle ~30 [ASSUMPTION: same computation] | Refit seconds (2026-09-10 `.rda` mtimes 10-12 s apart), 36 runs ~10 min | INT cap; read CONTROL 1b on the 147 null genes first |
| y_null + k g | Script; multiplier k^2/(1 + 0.022 k^2) | Fallback if unw < 5, k >= 5 | ~1 h | Shares noise with y_null |

Refuted: rebuild on the cell (implied coding h2 0.47, disk); ladder 10,600 (dominated, binary clamp); y_binary only (betas to 5.496, `02_build.log:50`; 8 outside [11,15]).

**Recommend** re-specify, then the refit; no rebuild. Changes it: unw well under 5 with a converged refit means an engine hunt (invNormalize, VR, collapse) first.

## D4. Score-noise arms

Three arms vs a sweep of s. Paired per gene, informative vs uninformative is null here (37/63, sign p 0.21 [ASSUMPTION: refuter's test on HEAD outputs]); counts differ by one gene (`:336,340`).

| Option | What changes | Expected effect | Cost | Risk |
|---|---|---|---|---|
| Keep three arms + paired stats | `analyse.R`: per-gene log10 p differences, sign test, gene bootstrap | 1-gene differences become interpretable | 0 runs | None |
| Provenance sidecar | Record realised within-cell cor (group files have no header syntax, `R/SAIGE_SPATest_Region_Func.R:14-25`) | Arms comparable across fixtures | 20 lines | None |
| Gated sweep | 5-6 s values, flex192 only, flex24/lin from rows | Power vs cor curve | 10-12 runs | Meaningless before an oracle gap |
| Replicates | Redraw residual and poly (`02_build_dataset.py:672-678`) | Real error bars | ~1 h + refits | Next study |

Refuted: oracle at 4 runs on y_quant_pc (stale model); sweep as specified (x-axis dominated by the low class); generator sd sweep (four refits per level, wrong sign at sd=0).

**Recommend** paired stats and the sidecar now; sweep only after a paired-resolvable oracle gap.

## D5. Cells

One cell vs a grid. A FlexRV run is one group, one max MAF (`R/SAIGE_Test_main.R:379-390`); the paper's cell is LoF + all missense (`flexrv.txt:146`), so high:medium is a mask-then-score design that hides the score.

| Option | What changes | Expected effect | Cost | Risk |
|---|---|---|---|---|
| Second one-cell run, high:medium:low | `run.sh:22` CELL env, cell into the tag (`:59`) and `analyse.R:20`; y_binary c0 x {flex192, flex24, unw} x 2 arms | flex192 minus flex24 +4 to +6 vs ~0 on high:medium [ASSUMPTION: independence NCP] | 12 runs; collapse is per carrier (`src/Main.cpp:2757`) | `low` includes synonymous-like variants |
| Third mask in `grid` | SAIGE's own 9-cell default (`step2_SPAtests.R:280`) | Fair standard-practice comparator | 4 runs at c10 | Not the paper's Regenie comparator |

Refuted: three max-MAFs (1e-2 near-duplicates 1e-3; 1e-4 is pseudo-marker only at c10); one Cauchy over the grid (not p_FlexRV; the CLI refuses); 9 x 24 one-cell runs (a different question).

**Recommend** both; y_binary now, y_quant_pc after its refit, never y_quant. Changes it: flex192 not beating flex24 by more than ~2 on the wide cell means the score-effect link, not cells, needs changing.