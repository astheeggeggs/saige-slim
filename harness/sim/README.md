# Simulated chr20 cohort + SAIGE step 0/1/2

Everything here is generated; nothing is committed to git (outputs live in
`harness/work/sim/`, which is gitignored). Run order:

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
