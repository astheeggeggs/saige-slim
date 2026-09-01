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
| outdir | `harness/work/sim` | `harness/work/sim50k` |
| disk | ~7 GB | ~14 GB |

`bench50k` exists because **the bundled N=1,000 fixture inverts the ranking it
is meant to measure**: `writeOutfile_single` (22% of the marker loop) and
exceptions (8.8%) are flat in N, while decode/impute/score/SPA all scale with
it. Its window and thinning come from the measured pilot in `PILOT.md`, not
from estimates -- mean recombination rate 0.98x the chromosome mean, ~19 s/Mb,
72,178 kept sites/Mb at 46,250 founders, 80.8% of variants at MAF <= 1e-4.

**`bench50k` stops before step 0 and says why.** A GRM estimated from one 25 Mb
window is not usable, for a genetic rather than a counting reason: 25 Mb of
chr20 is ~41 cM, so the number of effectively independent common markers is a
few hundred *however many markers the window holds*. The prior run measured
exactly that failure -- true first-degree pairs at median 0.493 against an
expected 0.5, but "unrelated" pairs at median 0.075 and the diagonal ranging
0.228-2.230, i.e. a sound estimator on an unusable input. (Population structure
is ruled out: `02_build_dataset.py` keeps EUR only.) The fix is the synthetic
independent scaffold, which is not built yet. `ALLOW_REGION_GRM=1` runs step 0
anyway as a diagnostic. Meanwhile `<prefix>.pedigree_truth.txt` holds the exact
relatedness matrix step 0 is trying to estimate, so step 0 can be scored as a
control rather than trusted as a dependency.

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
