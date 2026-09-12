# FlexRV simulation study: realistic MAF, annotations and scores at large N

<!-- Written 2026-09-10. Status (2026-09-12): the S2-S4 grid was RE-RUN at HEAD after #118/#126/#127/#117 (harness/sim/flexrv_study/results-2026-09-12.txt, prediction+outcome docs/state/prediction-flexrv-rerun-2026-09-12.md): skato/grid reproduce 2026-09-10 exactly, every FlexRV detection count fell by 1-2 (flex192 now 11 = grid 11 > unw 8 on y_binary), LEDGER #129 found (weighted pseudo-marker zeroing; read the collapse-0 column); the three decisions stand, memo docs/state/flexrv-study-decisions-2026-09-12.md. Earlier status: S0-S4 RUN 2026-09-10 late (harness/bench/NOTES.md "FlexRV simulation study S2-S4", harness/sim/flexrv_study/results-2026-09-10.txt): null calibration fine, the positive control FAILED on y_quant, two control expectations mis-specified (NOTES). S6 write-up is that section. Decisions for the owner
     are marked [DECISION]. Builds on harness/sim/ (the msprime/stdpopsim
     generator, README there) and on the F5 timing of 2026-09-10
     (harness/bench/NOTES.md, "FlexRV F5 at N = 50,000"). -->

## What the study is for

FlexRV's claim (Schwartzentruber et al. 2025) is power: 192 fixed weightings
of one pathogenicity score, Cauchy-combined, find more genes than any single
burden weighting while staying calibrated. Nothing on this branch has tested
either half. The bundled fixture cannot (N = 1,000 in-sample, two genes), and
the F5 timing run used `y_null` only. This study asks, on simulated data whose
MAF spectrum, LD, relatedness and annotation structure are realistic:

1. **Calibration.** Is `p_FlexRV` uniform under the null, on a quantitative
   trait and on a binary trait at 1% prevalence, with the collapse on and off?
   The pooled Cauchy row over 192 dependent p-values is exactly the object CCT
   is claimed to control "at least in the tail", and the MAC ladder makes the
   pool's members nested, not merely dependent.
2. **Power, and what it depends on.** At matched calibration, how does
   `p_FlexRV` rank against the unweighted burden, the beta(1,25) burden,
   SKAT-O's ACAT-O, and the oracle "best single weight per gene" (which is
   what a CCT over 192 can at most approach) -- as a function of how
   informative the score is? An uninformative score is the negative control
   of the power question: FlexRV must then cost power (192 tests, one signal),
   not gain it.
3. **The collapse operator.** SAIGE collapses ultra-rare variants by
   `max(w*G)`, the paper by additive pseudo-variants. The gap is bounded by the
   fraction of individuals carrying two or more qualifying ultra-rare variants
   in a gene; measure it before any p-value (design doc, "What enters the
   test"), then `--MACCutoff_to_CollapseUltraRare` 0 vs 10 on the same traits.
4. **Scale.** Wall time and peak RSS of the 192-set run against the one-set
   run at the largest N the machine allows, so the memory levers can be
   decided on numbers (F5 at 50k: 35 s / 805 MB vs 8 s / 506 MB per cell).

## How the paper itself checked calibration (read 2026-09-10, `flexrv.pdf`)

No simulation study. Three empirical checks on UKB, all with unrelated
individuals, no GRM, and a common-variant PGS as a covariate:

1. **Permutation**: IRNT residuals of the covariate-only null model swapped
   at random, 5 permutations x 62 quantitative traits, 5.7 million gene-trait
   tests. Observed/expected chi-square 1.02-1.04 for p in 1e-2..1e-4,
   "negligible" below 1e-4, nothing significant at FDR 1%; the MEDIAN lambda_GC
   was 0.22, which they read as lambda being useless for an omnibus CCT --
   the same conclusion as this branch's `lambda-gc-invalid-with-atoms`. The
   binary arm is stated as "well calibrated in permutation tests" with STAAR's
   SPA (max_iter 10, re-run at 200 below p 0.01), not quantified.
2. **Synonymous variants on real phenotypes**, scores drawn with replacement
   from the gene's own missense PrimateAI-3D scores: at FDR 1% a handful of
   known-confounded loci (IGLL5, LPA/SOD2, CST3), read as false positives
   with an explanation, not as miscalibration.
3. **QQ on real traits plus BHR intercepts**, to attribute the early rise to
   polygenicity rather than confounding.

One methods detail that matters here: **a weight function is excluded for a
gene when fewer than 5 individuals carry a non-zero-weight variant** -- a
per-SET carrier gate. Since LEDGER #118/#126 (2026-09-11) this engine has the
nearest thing to it: `--minGroupMAC_in_BurdenTest` reads the SET's unweighted
MAC on both region paths. It is still not the same gate -- alleles, not
carriers -- and its default is 1 rather than 5.

What this study adds that the paper could not: known truth (power against
the generating model, and the informativeness dependence); relatedness
through the GRM; a binary trait at 1% prevalence under SAIGE's SPA/ER rather
than STAAR's; the `max` collapse, which the paper does not have; and a null
that is a genuine trait (`y_null`) rather than a permutation, since permuting
here drives tau2 to 0 and switches the sparse-GRM machinery off (LEDGER
#104). The paper's two checks are replicated in-sim as extra controls: a
"synonymous-style" mask (the `low` class with scores resampled from the
gene's `medium` scores) on the signal traits, expected null; and the
per-set 5-carrier exclusion applied post hoc to the rows before pooling, to
see whether #118's cell-level gate changes any pooled p.

## What "realistic" means here, and what already exists

`harness/sim/` already simulates most of it (README, "What is simulated"):

- **MAF spectrum and LD**: stdpopsim `HomSap` chr20 under `OutOfAfrica_2T12`
  with the HapMapII recombination map; 45.8% singletons at N = 10k, more at
  50k. That is the coalescent's site-frequency spectrum under a fitted
  demography, which is as realistic as simulated genotypes get without real
  haplotypes; real exome data would add recurrent mutation and selection
  against LoF, both absent here. [ASSUMPTION: the absent selection makes the
  simulated LoF spectrum LESS rare than real pLoF, so the MAC ladder's low
  rungs are busier here than in a biobank.]
- **Relatedness**: 500 (10k) / 1,875 (50k) nuclear families by gene-drop,
  a pedigree-truth sparse GRM, so the mixed-model machinery is exercised.
- **Annotation**: three classes (high/medium/low) cut from a latent
  deleteriousness `a_lat = d + noise`, where `d` rises with rarity; the true
  effect class is a second noisy readout of `d` (cor 0.52 with the annotation).
  High is the LoF analogue (5% of coding variants), medium the damaging
  missense analogue.
- **Architecture**: spike-and-slab with MAF-dependent effects (alpha = -0.3),
  10% causal genes, a background polygenic term, quantitative and liability
  traits, `y_null` with zero chr20 h2, and at 50k an NCP ladder across causal
  genes so power is spread rather than saturated.

What is missing for FlexRV, in order:

- **A pathogenicity score.** `a_lat` is computed but not written. The score
  the paper uses is LoF = 1 and, for missense, the PrimateAI-3D percentile
  among the gene's possible missense variants. The analogue: `high` -> 1;
  `medium`/`low` -> the within-gene percentile of `a_lat + N(0, s^2)`, with
  `s` the [DECISION] score-noise knob giving three arms -- informative
  (`s = 0`, cor with the true effect class as high as the fixture allows),
  realistic (`s` chosen so cor(score, effect class) ~ 0.3-0.4, the range a
  good in-silico predictor reaches), and uninformative (score permuted within
  gene, cor 0: the negative control). The noise comes from its OWN RNG stream
  (`seed + 90002`) so adding it leaves every existing draw, hence the `.bed`
  and every trait, byte-identical -- which is the gate on the generator change.
- **Directional genes.** Fixture defect 2 (memory note, `saige-bench50k-fixture`):
  causal effects are 90% same-sign per variant but a gene's ~28 causal
  variants cancel at gene level (`|sum b| / sum |b|` 0.02-0.50). A burden method
  cannot win there and the comparison against SKAT-O is uninformative. Add
  `--frac-directional` (fraction of causal genes whose effects are drawn
  same-sign) so the study spans both regimes. [DECISION: 0.5.]
- **NCP on the collapsed burden.** Fixture defect 1: the NCP ladder is
  computed on the uncollapsed per-gene burden while the test collapses 81% of
  variants. For this study the NCP is a nuisance (the ladder control is
  Spearman(NCP, -log10 p) within a method); it does not need fixing to compare
  methods on the same genes, and is left as is, stated.

## Annotation x MAF from real summary statistics (owner's proposal, 2026-09-10)

The latent model above invents how many pLoF and missense variants a gene
carries at each MAF (`GAMMA = 0.9` on `-log10 MAF`). Real exome summary
tables carry that joint distribution, so the builder gets an
`--anno-table` path that REPLACES the quantile cut of section 5 when given:

- **Input**: a tab-separated table with columns `class` (high / medium / low,
  mapped from pLoF / missense / synonymous-or-benign by the caller),
  `maf_lo`, `maf_hi`, `count`; optionally `gene` for per-gene tables. Pooled
  rows give `P(class | MAF bin)`; per-gene rows give it per gene-size stratum
  (genes binned by their variant count), never per named gene, because the
  simulated genes are not real genes.
- **Assignment**: for each simulated coding variant, its MAF bin (at the
  study's N: bins are matched by PROPORTION within the bin, not by absolute
  count, since singleton fractions depend on N) draws the class from the
  table's `P(class | bin)` on a separate RNG stream (`seed + 90004`);
  genotypes, `y_null` and every unrelated draw stay byte-identical.
- **Score**: `high` -> 1; `medium`/`low` -> a within-gene percentile of the
  class-conditional latent (the existing `a_lat`, re-centred on the assigned
  class), so the score arms (informative / realistic / uninformative) are
  unchanged in meaning. If the source also supplies a score distribution for
  observed missense variants (a PrimateAI-3D or AlphaMissense percentile
  table), that replaces the percentile draw for the "realistic" arm.
- **Effects**: the true effect class stays a noisy readout centred on the
  ASSIGNED class (cor with annotation as now), so effects track annotation
  and FlexRV has something to find.
- **Controls added**: the realised class x MAF-bin table of the fixture is
  printed and compared with the input table (chi-square by bin, expected to
  match to sampling noise); the singleton fraction per class is reported
  next to the source's, with the N difference stated.

Candidate sources, in order of fit (which one is downloadable is the
owner's call; the file goes under `harness/work/` and its provenance -- URL,
release, date, md5 -- is recorded in the results file): the FlexRV paper's
supplement if it tabulates pLoF/missense counts per gene under MAF < 0.1%;
Genebass gene-level results (per annotation x max-MAF group; [ASSUMPTION:
carries a variant count per group]); the gnomAD v4 constraint table (per-gene
pLoF/missense/synonymous observed counts, no MAF stratification -- gives the
class mix per gene size, not the MAF dependence); Backman et al. 2021
supplementary mask counts. The coalescent has no selection, so without this
calibration pLoF at a given MAF is commoner here than in a biobank -- which is
exactly the distortion `P(class | MAF)` from real data removes.

## Effect model reassessed for the annotation table (2026-09-10, owner's question)

Three changes, the first mandatory:

1. **The within-class score-effect link was severed by the table path and is
   restored.** Under the latent model the annotation and the true effect class
   were two readouts of one latent `d`, so within the missense class a higher
   score meant a higher chance of being causal. The first table path replaced
   `d` by the class centre, making score and effect independent GIVEN the
   class -- under which a pathogenicity score can carry nothing beyond the
   class mask and FlexRV's within-class ranking has nothing to find. Now
   `a_lat = mu[class] + d_within + noise`, `e_lat = mu[class] + d_within +
   noise`, `--within-class-sd` (default 1.0; 0 severs the link and is the
   negative control of the power claim, not a default). Smoke config:
   cor(a_lat, e_lat | medium) = 0.66, cor(anno, true class) = 0.50.
2. **30% causal genes [DECIDED]**: 63 causal, 147 null at 210 genes
   (`FRAC_CAUSAL=0.3`). Power comparisons across three score arms need more
   than 21 genes; the null control keeps 147. Redraws `y_null` (an existing
   knob consumes main-stream draws), so the `y_null` column is no longer an
   invariant of the rebuild; the genotypes, pedigree, GRM and `var` lines are.
3. **Per-gene link shapes [DECIDED]**: `--link-mix step,lof,graded` draws
   each causal gene's relationship between latent pathogenicity and causality
   from: `step` (P(causal | true class), the original), `lof` (causal only in
   the true high class), `graded` (P(causal) logistic in the latent, so the
   fine within-class ranking matters most). The paper's central claim is that
   the best weight function differs by gene; with one shape a single fixed
   weighting is always near-optimal. `truth_genes.tsv` carries `link`, so
   detections are reported by link type per method.

## Constraints of this machine, and what N is possible

Measured 2026-09-10: 9.1 GB free disk, 16 GB RAM, 8 cores. The 50k fixture is
3.1 GB (bed 2.6 GB) and takes 65 min to rebuild; a 200k fixture is ~10 GB of
bed plus the tree sequence and PGEN mirror (~15 GB), 500k ~35 GB. So:

- **N = 50,000 is the study's N on this box**, and it is a fixture that
  exists. It is large enough that MAF < 0.001 spans MAC 1-100 (the MAC ladder
  half live, half twins), that a binary trait at 1% prevalence has 500 cases,
  and that step 2 costs 35 s per trait per cell for 192 sets.
- **No scale arm [DECIDED 2026-09-10: 50k is sufficient].** Power is set by
  the NCP ladder, not by N; calibration is per test; the F5 measurement gives
  memory (1.6 MB per set per cell) and time as linear in N, so 500k
  extrapolates. What 50k cannot show is the top of the MAC ladder (MAF < 0.001
  is MAC <= 100 here, <= 900 at 450k), which the real group summary statistics
  supply directly. To make room for the rebuild the chr20 correctness
  fixture's `.bed`/`.pgen` were deleted (5.5 GB; its tree sequence is kept, so
  `FROM=02 bash harness/sim/run_all.sh` restores it in ~30 min;
  `harness/work/sim/DELETED-2026-09-10.txt`).
- RAM: 192 sets at one cell is 1.6 MB per set at 50k (F5), 3 GB at 500k;
  12 cells at 500k would not fit in 16 GB without the memory levers.

## Design

One fixture build per N (the score export changes no existing bytes), then
the study is step-1 fits (one per trait) and step-2 runs (one per trait x
method x collapse), all against frozen models, all recorded with commit,
`SAIGE.so` md5 and the fixture's seed.

**Traits** (per N): `y_null` (quantitative, chr20 h2 = 0), `y_null_bin`
(1[y_null > quantile], the recipe of 2026-09-08 that keeps tau2 non-zero;
permuting is ruled out, LEDGER #104), `y_quant` and `y_binary` on the NCP
ladder, half the causal genes directional.

**Score arms**: informative / realistic / uninformative, one `score:` line
each (three group files; the var/anno lines identical).

**Methods**, every one on the same cell `high:medium` at max MAF 0.001 (the
paper's filter; `high:medium:low` as a second cell only if time allows):
- FlexRV default 192 sets, BURDEN, collapse 10 and collapse 0;
- FlexRV score-only (`--flexRV_MAFTransforms=unw`, 16 sets) and MAF-only
  (`--flexRV_scoreTransforms=unw`, 12 sets): which axis carries the gain;
- unweighted burden, beta(1,25) burden, SKAT-O (ACAT-O row), the current
  recommendation's baseline;
- oracle: min over the 192 rows of the same run, Bonferroni-corrected by 192
  (an upper bound on what CCT can deliver; free, it is in the same file).

**Outputs per (N, trait, score arm, method)**: the region file; from it the
per-gene p, and for FlexRV the served/evaluated counts and the per-gene
distinct-row count.

## Controls before conclusions

Written before the run, checked first, a misbehaving control halts.

- **Generator**: after the score-export change the 50k `.bed`, `.bim`, `.fam`,
  `pheno.txt` and `groupFile.txt` var/anno lines are byte-identical to the
  2026-09-01 build (`cmp`); the score line has LoF = 1 on every `high`
  marker, values in [0, 1], and cor(score, true effect class) at the three
  intended levels, printed.
- **Null (negative)**: on `y_null` and `y_null_bin`, per method, the ECDF of
  the 189 + 21 gene p-values at alpha = 0.05 and 0.01 with a binomial 95%
  interval (lambda is invalid here: the burden rows carry point masses at
  p = 1 in cells with no real marker, memory note `lambda-gc-invalid-with-atoms`).
  Expected: within the interval for every method and both collapse settings.
  FlexRV's pooled row is the one to watch: CCT is conservative for
  positively dependent nested p-values in the body and exact only in the tail.
- **Null on signal traits**: the 189 non-causal genes of `y_quant`/`y_binary`,
  same statistic. Expected: the same.
- **Positive**: on `y_quant_pc`-style power (the ladder), Spearman(NCP,
  -log10 p) over the 21 causal genes per method, expected > 0.5 for every
  method (it was 0.68 for ACAT-O binary, 0.20 quant at h2 = 0.01 -- the quant
  arm may be underpowered at 50k; if so the positive control fails for every
  method and the power comparison is not run, only calibration).
- **Score-informativeness monotonicity**: FlexRV's detections must be
  monotone in the score arm (informative >= realistic >= uninformative), and
  on the uninformative arm FlexRV must not beat the unweighted burden -- if it
  does, the "gain" is the MAF axis, which the MAF-only arm isolates.
- **The collapse bound**: the >= 2-carriers fraction per gene from the
  genotypes (plink2 `--score` over the cell's ultra-rare markers, or a small
  Python pass over the `.bed`), reported before the 0-vs-10 comparison; if it
  is under ~1% on every gene the `sum_geno` question closes there.
- **Two ways to one number**: the pooled `p_FlexRV` recomputed in R as the
  CCT of the 192 printed rows (as `test-invariants-goldens.R` does) on every
  file.

## S0 status (2026-09-10, smoke-tested on the tiny config)

Landed in `harness/sim/`: `02_build_dataset.py` gains `--anno-table` (the
empirical class draw above; the realised class x bin table printed against
the input, matched to sampling noise on the smoke config), `--frac-directional`
(`directional`/`dir_index` columns in `truth_genes.tsv`), and the `a_lat`/
`e_lat` columns in `variants.tsv.gz`; `run_all.sh` passes `ANNO_TABLE` and
`FRAC_DIRECTIONAL`. `harness/sim/flexrv_study/make_scores.py` writes the four
score-arm group files (informative / realistic / uninformative / synonymous)
and prints each arm's cor(score, true class): 0.49 / 0.34 / -0.01 / 0.08 on
the smoke config at `--noise 1.5`. `run.sh` there is the study driver
(trait x arm x method x collapse, provenance per row, refuses at load > 2.5).
Invariants: the smoke `.bed` is byte-identical with and without the table
and the directional flag; the 50k build's checksums to reproduce are in
`harness/work/flexrv_study/sim50k-invariants-2026-09-01.txt`. One stated
change: under the table path the latents are re-centred on the assigned
class, and cor(anno, true effect class) after cutting comes out ~0.68 on the
smoke config against ~0.52 under the latent model; the study reports the
50k value and does not tune it.

**The table (2026-09-10, built).** Source: All of Us v8, EUR, N = 171,473,
BMI (BRaVa masks; `harness/work/flexrv_study/aou/`, md5s in the table's
header). The gene results carry no MAC columns and the single-in-group file
has the ultra-rare pseudo-variants pruned, so the route was: the per-variant
chr20 annotation table (BRaVa `ANNOTATION`: pLoF -> high, damaging_missense
-> medium, other_missense + synonymous -> low, non_coding dropped) joined to
the single-in-group file's AF for every variant that survived the collapse
(MAC >= 11, MAF <= 0.01, the cells); the ABSENT coding variants (395,747 of
417,199) form the ultra-rare bin under the stated assumption that
EUR-monomorphic and QC-dropped variants share its class mix (common variants
above 1% are a fraction of a percent of them). Realised
`harness/work/flexrv_study/anno_table_aou_chr20.tsv`:

| MAF bin | n | low | damaging missense | pLoF |
|---|---|---|---|---|
| [0, 2.9e-5) = MAC <= 10, absent set | 395,747 | 0.809 | 0.132 | 0.059 |
| [2.9e-5, 1e-4) | 14,524 | 0.864 | 0.110 | 0.026 |
| [1e-4, 1e-3) | 5,774 | 0.875 | 0.098 | 0.027 |
| [1e-3, 1e-2) | 1,154 | 0.905 | 0.078 | 0.017 |

Trend check: damaging missense rises monotonically as MAF falls; pLoF is flat
over the two exact rare bins and jumps in the ultra-rare bin, the known
singleton enrichment of LoF under selection, so the absent-set assumption is
consistent with the exact bins rather than contradicted by them. Against
the latent model this is the correction the study wanted: the simulated
`medium` class was 30% of coding variants; the real damaging-missense class
is 8-13%. The annotation table also carries REVEL and CADD per variant,
which the score arms do not use (the effects are synthetic) but which give
an empirical missense-score marginal to hold the "realistic" arm against.

The 50k rebuild with this table and `FRAC_DIRECTIONAL=0.5` was started
2026-09-10 (`harness/work/flexrv_study/rebuild50k.log`); the checksum gate,
step 1 for `y_null`, `y_null_bin`, `y_quant`, `y_binary` (all refitted on
the rebuilt fixture so the model set is consistent), then S2-S4 follow.

## Stages, cost, deliverables

| stage | what | cost | gate |
|---|---|---|---|
| S0 | `02_build_dataset.py`: write `a_lat` to `variants.tsv.gz`, emit the three score lines (own RNG stream), `--frac-directional`; run_all knob | 1 h code + 65 min rebuild at 50k | existing outputs byte-identical; score controls print |
| S1 | step 1 for `y_null_bin` (exists), `y_quant`, `y_binary` at 50k (exist), directional variants of the signal traits (rebuild changes them: [DECISION] whether to keep the 2026-09-01 traits or regenerate with directional genes -- the latter invalidates the existing 50k step-2 baselines, which the memory note already says are not usable) | ~10 min per fit | theta printed, tau2 > 0 |
| S2 | the null runs: 2 traits x 3 score arms x (192, 16, 12, unw, beta, SKAT-O) x 2 collapse ~ 60 step-2 runs at ~10-40 s | ~30 min | ECDF controls |
| S3 | the signal runs, same grid on `y_quant`/`y_binary` | ~30 min | positive control, monotonicity |
| S4 | the collapse bound and the 0-vs-10 comparison | ~10 min | reported with S2/S3 |
| S6 | write-up: `harness/bench/NOTES.md` section + a results file under `harness/sim/flexrv_study/`, next.md, a handoff | 1 h | |

Driver: `harness/sim/flexrv_study/run.sh` (arms x methods, records commit /
md5 / seed per row, refuses at load > 2.5 as `bench.sh` does),
`analyse.R` (controls first, then the tables). Everything under
`harness/work/` is gitignored; the scripts and the results tables are
committed.

## Decisions for the owner

1. **N.** 50k, no scale arm. [DECIDED 2026-09-10.]
2. **Score noise levels.** Three arms as above; the "realistic" target
   cor ~0.3-0.4. Alternative: sweep `s` on a grid and plot power vs cor.
3. **Directional genes.** Half the causal genes, signal traits regenerated.
   [DECIDED 2026-09-10; `--frac-directional` landed and smoke-tested.]
4. **Cells.** One cell (`high:medium` at 0.001, the paper's) or also the
   12-cell default grid (12x the cost, and the RAM question at scale).
