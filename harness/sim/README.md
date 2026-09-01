# Simulated chr20 cohort + SAIGE step 0/1/2

Everything here is generated; nothing is committed to git (outputs live in
`harness/work/sim/`, which is gitignored). Run order:

## Two configurations, one generator

`run_all.sh` takes `CONFIG`. Both configurations run the same scripts and
differ only in variables, so a generator fix lands in both and they cannot
drift. Anything set in the environment overrides a preset; `DRYRUN=1` resolves
the configuration and stops (the full run is ~80 min).

| | `chr20` (default) | `bench50k` |
|---|---|---|
| purpose | correctness / controls | step-2 benchmarking |
| region | whole chr20 (64.44 Mb) | chr20:34-59 Mb (25 Mb) |
| N | 10,000 | 50,000 |
| founders | 9,000 | 46,250 |
| genes | 550 | 210 |
| intergenic thinning | none | 0.070 (rare intergenic only) |
| prevalence | 0.10 | 0.01 |
| trait scaling | fixed `--h2` | NCP ladder 1-60, binary mult `auto` |
| GRM for steps 1-2 | step 0's estimate | pedigree truth |
| outdir | `harness/work/sim` | `harness/work/sim50k` |
| disk | ~7 GB | ~14 GB |

`bench50k` exists because **the bundled N=1,000 fixture inverts the ranking it
is meant to measure**: `writeOutfile_single` (22% of the marker loop) and
exceptions (8.8%) are flat in N, while decode/impute/score/SPA all scale with
it. Its window and thinning come from the measured pilot in `PILOT.md`, not
from estimates -- mean recombination rate 0.98x the chromosome mean, ~19 s/Mb,
72,178 kept sites/Mb at 46,250 founders, 80.8% of variants at MAF <= 1e-4.

### Step 0 is a control, not a dependency

`03b_pedigree_grm.py` **writes** the sparse GRM implied by the pedigree instead
of estimating it. The pedigree is known here, so there is no reason to infer
it -- and over one region it cannot be inferred well: 25 Mb of chr20 is ~41 cM,
so the number of effectively independent common markers is a few hundred
*however many markers the window holds*. **[ASSUMPTION]** -- that is an a priori
argument; no GRM has been built on a 25 Mb window here. What is measured is the
full-chromosome case below, at 2.6x the genetic length, which is already
marginal.

Values are on SAIGE's GRM scale (twice kinship): 1.0 diagonal, 0.5 first-degree.
Founders are independent coalescent draws and a family's two parents are
unrelated, so within-family first-degree pairs are the *complete* set of nonzero
off-diagonals -- a property of the construction, not an approximation. The file
is `coordinate real symmetric`, lower triangle, 1-indexed against its
`.sampleIDs.txt`, which is what `Matrix:::readMM` wants
(`R/SAIGE_SPATest_Region.R:26`).

So `RUN_STEP0` defaults **on** for `chr20` (which consumes step 0's output) and
**off** for `bench50k` (which consumes the truth). Setting `RUN_STEP0=1` under
`bench50k` runs step 0 purely to be scored against the truth. `GRMFILE`/`GRMIDS`
override which matrix steps 1-2 read.

The truth matrix is verified through SAIGE's own reader, not by inspection: at
smoke scale it comes back `dsTMatrix`, symmetric, diagonal exactly 1, all
off-diagonals exactly 0.5, nothing below the 0.125 cutoff, `.sampleIDs.txt`
identical to `.fam` order -- and the row sums land at 1.0 / 2.0 / 2.5 for
singletons / parents / children, which is the pedigree arithmetic derived
independently of the file.

### What the prior run actually measured, and one trap in reading it

Two GRM configurations were run on the SAME N=10,000 full-chr20 dataset, and
they are easy to mix up because only the current one is described above:

| | superseded (50-variant prune) | current (1 Mb prune) |
|---|---|---|
| markers offered / used | 14,047 / 5,000 | 3,055 / 3,000 |
| true first-degree pairs | median 0.493 (IQR 0.455-0.536) | median 0.500 (IQR 0.456-0.543) |
| diagonal | median 0.986, range **0.228-2.230** | median 0.998, range **0.353-1.314** |
| pairs kept | 1,118,079 at cutoff 0.05 | 4,431 at cutoff 0.125 |

**The trap: `value on OTHER pairs: median` is not a measure of background
relatedness.** The sparse matrix only retains pairs above the cutoff, so that
median sits just above wherever the cutoff was put -- 0.075 at cutoff 0.05,
0.168 at cutoff 0.125. Reading those two as "the current recipe is twice as
contaminated" is backwards; it is 250x less contaminated by pair count. Compare
the diagonal (not truncated) or the extra-pair count at a fixed cutoff.

### Open: why is the diagonal so wide? (rare-marker inflation, NOT settled)

The standing hypothesis is that the GRM is dominated by near-threshold rare
markers. It is checkable in closed form with no GRM at all, from the pruned
set's MAF spectrum alone. For standardised genotypes over M independent markers
under HWE, with `h_j = 2 p_j (1 - p_j)`, the fourth moment is `1/h_j` exactly
(the binomial n=2 fourth central moment is `h`), so

    Var(diagonal) = (1/M) * mean_j(1/h_j - 1)

**[ASSUMPTION, order-of-magnitude]** under a neutral `1/p` spectrum truncated
at MAF 0.01, `mean(1/h)` is ~12, giving sd ~0.062 and an expected range over
10,000 individuals of about [0.76, 1.24]. The observed upper tail (1.314) is
close to that; the observed lower tail (**0.353**) is nowhere near it, and a
sum of positive terms should not have the longer tail on the low side. So rare
markers alone look **insufficient**, and the hypothesis should not be treated
as the explanation until the real spectrum is used.

The test, in order of cost: (1) compute `mean(1/h)` over `data_for_grm`'s actual
frequencies and compare the predicted sd against the observed diagonal -- one
pass over a `.frq`, no GRM; (2) if it reconciles, re-run step 0 at
`minMAFforGRM=0.05`, where `1/h <= 10.5` bounds the inflation and the sd should
drop by ~2x. Doing (2) first would tell you the diagonal narrowed without
telling you whether the model that predicted it was right.

## Run order

```sh
PY=<venv with msprime,stdpopsim>/bin/python
$PY harness/sim/01_simulate.py --out harness/work/sim/chr20.trees      # ~17 min
$PY harness/sim/02_build_dataset.py --trees harness/work/sim/chr20.trees \
      --outdir harness/work/sim                                        # ~15 min
bash harness/sim/03_plink_prep.sh                                      # ~8 min
$PY harness/sim/04_add_pcs.py harness/work/sim 4
SAIGE_LIB=<lib> NTHREADS=6 RELCUT=0.125 NRANDMARK=3000 \
  bash harness/sim/04_run_saige.sh step0
.pixi/envs/default/bin/Rscript harness/sim/05_check_grm.R              # GRM vs pedigree
RELCUT=0.125 NRANDMARK=3000 NTHREADS=6 bash harness/sim/07_run_pipeline.sh
.pixi/envs/default/bin/Rscript harness/sim/06_analyse.R                # results vs truth
```

## What is simulated

- **Genotypes**: stdpopsim `HomSap` chr20 (GRCh38, 64.44 Mb), `HapMapII_GRCh38`
  recombination map, demography `OutOfAfrica_2T12` (Tennessen et al. 2012;
  AFR + EUR). 9,000 EUR founders + 1,000 AFR; only EUR is exported.
  2,267,192 EUR-polymorphic biallelic variants, 45.8% singletons.
- **Relatedness**: 8,000 unrelated founders + 500 nuclear families of 4
  (father, mother, 2 offspring) produced by gene-drop along the HapMap map
  (1.06 crossovers/meiosis). 2,500 true first-degree pairs.
- **Architecture**: spike-and-slab, slab sd scaling as `[2p(1-p)]^(alpha/2)`
  with **alpha = -0.3**; chr20 h2 = 0.01 (60% coding / 40% non-coding), plus a
  background polygenic term (h2 = 0.35) representing the unsimulated rest of
  the genome, transmitted through the pedigree.
- **Annotations**: 550 genes with exon structure. A latent deleteriousness
  score `d = 0.9*z(-log10 MAF) + N(0,1)` drives TWO noisy readouts — the
  annotation and the true effect class — so annotation predicts effect only
  imperfectly (realised cor 0.52) and high/medium are rarer than low.

| class  | n      | mean -log10 MAF | singleton % | MAF>1% |
|--------|--------|-----------------|-------------|--------|
| high   | 2,466  | 4.174           | 72.0        | 0.00   |
| medium | 14,795 | 4.081           | 60.2        | 0.14   |
| low    | 32,057 | 3.574           | 36.8        | 9.58   |

- **Traits**: `y_quant`/`y_binary` at h2 = 0.01 (as specified; underpowered by
  construction), `y_quant_pc`/`y_binary_pc` at h2 = 0.20 (positive controls),
  `y_null` with ZERO chr20 h2 (negative control). Binary = liability threshold,
  10% prevalence. Only 55 of 550 genes are causal; the other 495 are exactly
  null and are the built-in negative control for the gene tests.

## Deviations from the upstream recipe, each forced by measurement

1. **`--indep-pairwise 50 5 0.05` is wrong for WGS.** The window is 50
   *variants*; on 2.27M WGS variants that is 1,410 bp, so it removes almost no
   LD. Measured: 50-variant window -> 14,047 markers but M_eff 1,732 (GRM
   off-diagonal sd 0.0240); 1 Mb physical window -> 3,055 markers, M_eff 2,702
   (sd 0.0192) and 2.7x fewer spurious related pairs. `03_plink_prep.sh` uses
   `--indep-pairwise 1000kb 0.05`.
2. **`--relatednessCutoff 0.05` is unusable on one chromosome.** SAIGE's
   default assumes a genome-wide GRM (off-diagonal sd ~0.005, so 0.05 is ~10
   sd). chr20 alone gives sd 0.0192 -> 0.05 is 2.6 sd and retains ~528,000
   spurious pairs. Cutoff 0.125 gives 1,951 spurious vs 2,480/2,500 true pairs
   recovered (99.2% sensitivity, median value 0.500 against a theoretical 0.5).
3. **`--numRandomMarkerforSparseKin 5000` exceeds the available markers** (only
   3,055 survive 1 Mb pruning); 3,000 is used.
4. **`--useSparseGRMforVarRatio TRUE` is required** in step 1. Without it the
   variance-ratio file has no `sparse` rows and step 2 hard-stops when given
   `--sparseGRMFile` (`R/readInGLMM.R`, `Get_Variance_Ratio`).
5. **PGEN psam must put IID first** — see the reader limitation below.

## Generator invariants that are asserted, not assumed

- **Every causal gene has >= 1 causal variant.** Causal variants are drawn
  independently per variant, so a causal gene could come out with none: h2
  exactly 0 while still counting as a positive control, which drags the
  control's rank statistic toward the null and makes a correct pipeline look
  broken. Causal genes are now drawn only from genes with >= 2 coding variants
  (the group file skips the rest), any that still draw zero get one forced with
  the same per-class propensity, and an assert holds it. At 200 genes / 50
  causal, 8-14 genes per seed need forcing; after it, 0 of 50 have h2 == 0.
- **Thinning is MAF-stratified.** Uniform thinning would gut the common pool
  the GRM and the LD-pruned variance-ratio set are drawn from. Only intergenic
  variants below `--thin-maf` (default 0.01) are thinned. Verified at retention
  0.07: genic 721/721 kept, common 4,618/4,618 kept, rare intergenic
  728/10,609 = 0.069 against an expected 0.070 (0.4 binomial SE).
- **The NCP ladder delivers its drawn targets to 4e-16 relative error.** With
  `--ncp-ladder lo,hi` each causal gene gets a target non-centrality drawn
  log-uniformly and its effects are scaled to hit it, so `--h2` then sets only
  the non-coding component and total h2 is a reported outcome. A fixed h2 gives
  every causal gene the same power, which at N=50,000 saturates the positive
  control (p ~ 1e-25) exactly where BURDEN/SKAT/SKAT-O most need separating;
  the ladder makes the control `Spearman(NCP, -log10 p)` instead of a hit count.
- **The binary arm is scaled separately.** At prevalence K a binary trait buys
  only `phi(t)^2/(K(1-K))` of the non-centrality a quantitative one gets from
  the same liability h2 -- 0.0718 at K=0.01. `--h2-binary-mult auto` solves for
  the equalising multiplier from the current prevalence (13.94 at K=0.01,
  measured net ratio 1.000x) so it cannot go stale. It is **clamped loudly** when
  it will not fit the liability budget of 1, which it does not for a strong
  positive-control arm; the warning states the effective multiplier, because a
  silently weakened positive control is the failure this exists to prevent.
- The causal-gene listing prints **causal genes only**. Ranking all `n_genes`
  and taking the top 10 pads the table with null genes at h2 = 0 whenever
  there are fewer than 10 causal genes, which reads as a bug and is not.

## Findings about the codebase

- `src/PGEN.cpp:248` (`readPsamFile`) skips the `#` header line and then takes
  `line_elements[0]` as the sample ID. That is the FID for any psam written
  with the default column set, so a normal plink2 `--make-pgen` output is
  rejected with "At least one subject requested is not in psam file". The
  bundled fixture `extdata/input/genotype_100markers.psam` has FID == IID, so
  no golden covers this. Worked around by writing the psam with
  `--make-pgen psam-cols=sex` (gives `#IID SEX`). NOT fixed, per the
  bug-for-bug rule -- candidate LEDGER entry.
- PLINK vs PGEN step 2 over 268,003 markers: **byte-identical** output, with
  the PLINK run 6-way sharded and the PGEN run single-threaded.
