# Handoff 2026-09-25: the AoU production step-1 refit and the first full step-2 run

Two days (2026-09-24 and 25) of production work in the All of Us Researcher
Workbench, on the `strip` branch's image `1.5.2-dev-f2fb24da`. Everything below
is measured unless marked otherwise. Nothing in this handoff changes repo code;
what it changes is three LEDGER rows (#161 confirmed, #158 validated, #162
raised) and several beliefs that were written down wrongly.

## What was run

Step 1 refitted: **167 null models**, all `ok`, across five ancestries.

| ancestry | models | location | relatednessCutoff | tol |
|---|---|---|---|---|
| afr / eas / eur / sas | 42 / 25 / 43 / 15 | `nullglmm_v4/` | 0.05 | 0.02 |
| amr | 42 | `nullglmm_v4_rc0.125/` | **0.125** | 0.02 |

Step 2 then ran both arms -- FlexRV (`--r.corr=1`, one annotation group, one max
MAF, the AlphaMissense score line) and conventional SKAT-O (six masks x three max
MAFs) -- over all 22 autosomes for every model, binned by chromosome. 416 jobs for
afr/eas/sas, 336 for amr, 860 for eur.

## The AMR story: two separate causes, and neither fix was doing the other's job

AMR had been the blocker for weeks, and the diagnosis had conflated two
independent problems. They separate cleanly:

**A truncated step-0 GRM caused the step-2 hang.** The GRMs under
`data/saige/step0/` hold only 52-61% of samples for four of five cohorts; the
correct ones are under `run/make_sparse_grm/`. This was found because an older
`.rda` contained more samples than the GRM in use could possibly hold. The
signature everyone had been reading -- AMR's null variance ratio of **6.91**
against AFR's 1.09, taken as "the GRM matters enormously for this admixed
cohort" -- was an artefact of the truncation. Refitted on the full GRM, AMR's
null ratios are **0.970-1.966, median 0.997**, measured over all 42 models. At
the OLD cutoff (0.05) with the CORRECT GRM the maximum is 2.004, so it is the
truncation and not the cutoff that produced the 6.91.

**The 0.05 relatedness cutoff caused the step-1 non-termination and the `tau2`
collapse.** AMR's GRM is built from 5,000 random markers and carries 46,204,255
non-zero entries at 0.05; thinning at 0.125 leaves 349,445, a 132x reduction.
Working backwards from the counts, 0.05 sits at only ~1.8-2.1 sigma of the
off-diagonal noise for a 5,000-marker estimate, so at that cutoff the matrix is
mostly noise floor; 0.125 is ~4.5-5 sigma, and the ~150,000 surviving pairs
exceed the Gaussian tail by one to two orders of magnitude, so they are real
relatedness rather than noise. [ASSUMPTION: off-diagonals approximately Gaussian
and pairs independent, which is what allows inverting the count into a sigma.]
At 0.125, AMR HTN fits in **58 s** where it had never finished, with
`tau2 = 0.0741`, and 33 of 42 models have a non-zero random effect where ~95% had
collapsed to zero.

Consequence for the code: **the recorded verdict "AMR CANNOT use a sparse GRM in
step 2" is void**, and every ancestry now passes `--sparseGRMFile`. The
`sparse` variance ratio is 0.99977-1.0014 across all five ancestries, which is
the condition for that path being well calibrated.

**PCG versus direct Cholesky is not a lever here.** Measured on amr/HTN at 0.125:
56 s versus 58 s, with every variance ratio agreeing to ~1e-7. Keep the default
(`--usePCGwithSparseGRM=FALSE`).

## `--relatednessCutoff` must match between the steps, and nothing checks it

This is the one thing in the pipeline where a wrong value produces a clean,
complete, silently wrong run.

Step 2 applies its OWN cutoff when it loads the GRM (`setSparseSigma_new`,
`R/SAIGE_SPATest_Region.R:61`, `which(sparseGRM@x < relatednessCutoff)`) and
**both wrappers default it to 0**, i.e. no thinning
(`extdata/step2_SPAtests.R:368`, `extdata/step1_fitNULLGLMM.R:352`). A model
fitted at 0.125 and tested at the default would be tested against all 46.2M
entries while its `sparse` variance ratio describes a denominator of 349,445.
Nothing in SAIGE compares the two.

The other four ancestries are safe **for a reason, not by luck**: step 0 wrote
those files AT 0.05 (`..._relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx`,
built by `createSparseKinParallel(..., relatednessCutoff)`), so the `.mtx`
already contains only pairs above 0.05 and both step 1's `0.05` and step 2's
default `0` are no-ops on it. AMR is the only ancestry where the flag does real
work.

Two smaller facts from the same reading. Step 1 thins with `drop0(tol = cutoff)`,
which zeroes `|x| <= cutoff`; step 2 uses `x < cutoff`. An entry exactly equal to
the cutoff is dropped by one and kept by the other -- float equality on a GRM
value, so measure-zero, but the two matrices agree only up to that boundary. And
step 2's variance-ratio length check reads `nrv <- length(ratioVec_null)`
(`R/readInGLMM.R:451`), i.e. the `null`-type rows only -- so a file with three
types and two MAC categories has `nrv = 2` and matches the shipped default
`cateVarRatioMinMACVecExclude="10,20.5"`. Both wrappers ship the same default, so
they agree by construction unless one run overrides it.

## Measured costs and machine sizing

| | per-region / per-cell | peak RSS | notes |
|---|---|---|---|
| step 1, amr at 0.125 | 58 s per fit | 1,834 MB | 2 iterations; single-threaded whatever `--nThreads` says |
| step 2 FlexRV, afr | 0.90 s/region | 548-565 MB | `--r.corr=1`, so `SKAT-O cells: 0` |
| step 2 SKAT-O, afr | 0.52 s/region | 1,035-1,115 MB | 2,794 cells on chr21 |
| step 2 FlexRV, amr | 0.95 s/region | 1,362 MB | the excess over afr is the GRM LOAD |
| step 2 SKAT-O, eur | **0.1198 s/cell** | 6,057-6,170 MB | chr19 and chr1 agree to 4 s.f. |

Three sizing findings worth keeping.

**AMR's step-2 memory floor is set by the unthinned GRM, not by the trait.** Step
2 reads the whole `.mtx` with `readMM` before subsetting and thinning, so AMR's
peak carries all 46.2M entries even though 349,445 survive -- at 16-24 bytes per
triplet that is 0.7-1.1 GB, which accounts for essentially the entire 2.5x excess
over AFR on the same arm and chromosome.

**EUR's peak is flat across chromosomes, not proportional to work.** chr1 has 36%
more cells than chr19 but peaks 2% higher (6,170 versus 6,057 MB). So it is a
property of the model -- sample count, GRM load, the widest single region -- and
that makes the machine decision robust: `n2-highmem-2` (2 vCPU, 16 GB) for EUR,
`n2-standard-2` for everyone else. This also explains the older non-monotonic
observation (chr1 1,188 MB, chr2 2,134, chr19 2,230 on AFR).

**EUR's per-cell cost is 3.19x AFR's**, against 2.81x predicted from the sample
ratio (224,301 / 79,956), so mildly superlinear -- and constant across
chromosomes, which makes it usable for bin sizing arithmetically.

## Chromosome binning: balance by size, and the largest chromosome is the floor

Contiguous bins are the worst possible split because region count falls
monotonically with chromosome number. Balanced by group-file size (the same proxy
the wrapper uses to pick its calibrator), with longest-processing-time greedy
packing.

The number that matters is the LONGEST bin, not the average, and granularity
floors it: chr1 alone is ~10.3% of the genome's cells, so no bin can ever be
smaller. Measured for EUR SKAT-O (7.9 h for a whole genome at 0.1198 s/cell):

| bins | largest bin | FlexRV h | jobs |
|---|---|---|---|
| 8 | 13.6% | 1.87 | 688 |
| 9 | 12.0% | 1.65 | 774 |
| 10 | 10.9% | 1.49 | 860 |
| 12 | 10.3% | 1.42 | 1032 |

10 is the knee: going to 12 buys 4 minutes off the longest job for 172 more jobs
and ~94 GB more in image pulls, because 10.3% IS chr1. Used 4 bins for afr and
amr, 1 for eas and sas, 10 for eur.

## Spot preemption: two more datapoints, and the rate is stable

Prior observation was 15 of 158, all in `us-central1-c`, read as correlated by
zone. Added:

- afr/eas/sas, 416 jobs on spot: **36 failures (8.7%)**, all preemption-coded.
- amr + eur, 1,196 jobs on spot: **171 failures (14.3%)** -- amr 40 of 336
  (11.9%), eur 131 of 860 (15.2%), consistent with a larger simultaneous
  footprint in the same zones.

So the rate is stable around 9-15% and rises with concurrency. The
spot-then-on-demand-mop-up pattern works and is what to use. **Zone attribution
was NOT re-measured on these two batches** -- the failure spread across all six
job-name prefixes is NOT evidence against zone correlation, because the launch
loop interleaves prefixes so any contiguous block of submissions spans all of
them.

## The resume check must key on delivered OUTPUT, not on submission status

The launcher's `already()` originally read the run registry's `rc`, which is the
`aou_dsub` SUBMISSION status -- every preempted job has `rc = 0` and no output,
so the mop-up pass would have submitted nothing and reported success. It now
tests whether the results file exists in GCS. Two further traps in that check,
both of which produced a confident wrong count before being fixed:

- The output prefix embeds the bin's chromosome list, and the job NAME carries a
  `-b{i}` suffix for multi-bin ancestries. Checking the un-suffixed name against
  files written under the suffixed one reported all 336 afr jobs as outstanding
  while eas and sas (one bin, so no suffix) were correct.
- It cannot distinguish "still running" from "never submitted", so submitting
  while anything for the same ancestry is in flight duplicates it. The guard is
  a `dstat --status RUNNING` check scoped to the ancestries being submitted.

The completion criterion is not a job count: a job that finished having lost a
chromosome still delivers a results file, so the check is the `chrStatus` sidecar
sweep -- every trait-arm's bins' `requested` lists must union to 1..22.

## EUR T2Diab: the calibration evidence

Burden p-values have atoms, so **lambda_GC is not usable** -- a point mass breaks
it in both directions. Calibration was assessed as the ECDF at fixed alpha with a
binomial interval on the count. 18,443 genes.

| mask at maxMAF 0.001 | a=0.05 | a=0.01 | a=1e-3 |
|---|---|---|---|
| **synonymous** | 1.03 | 1.09 | 1.09 |
| **other_missense** | 1.06 | 0.99 | 1.20 |
| damaging_missense | 1.10 | 1.33 | 1.53 |
| pLoF;damaging_missense | 1.15 | 1.43 | 2.17 |
| pLoF | 1.12 | 1.46 | **2.87** |
| SKAT-O pooled (Cauchy) | 1.31 | 1.50 | 2.17 |
| FlexRV pooled (Cauchy) | 1.38 | 1.77 | **2.82** |

`synonymous` -- the negative control -- is inside its interval at all three
thresholds, and so is `other_missense`. The enrichment is confined to the
functional categories and GROWS as alpha tightens. Miscalibration inflates
roughly uniformly; this is signal. Supporting: **0 degenerate `SE_Burden` of
329,393 cells** (LEDGER #117's NA case never fired), and the largest atom is
`p = 1` on 2,120 cells (0.6%), FlexRV's largest being 3 rows.

Two things to know when reading these files. The pooled row is
`Group == "Cauchy"` with `max_MAF` NA and every column but `Pvalue_Burden`
blank -- for FlexRV it IS the result (`p_FlexRV`), and anything keyed on `Region`
without filtering it takes the NA row instead. And **SAIGE writes combined masks
with `;` although the command line takes `:`** -- `pLoF:damaging_missense` in,
`pLoF;damaging_missense` out -- so hardcoded mask strings give EMPTY series with
no error.

## Validated against the previous run (2024-03-26), 166 trait-by-ancestry pairs

The previous results under `saige/run/step2_group/combined_files/` have the
IDENTICAL schema (same 13 columns, same `;`-joined mask labels, same three max
MAFs), so a full cell-level comparison was possible on `(Region, Group, max_MAF)`.

**On identical data the new pipeline reproduces the old to 0.07%.** Medians of the
per-trait median `log10(p_new/p_old)` across 157 traits, by p-value column and by
band of the OLD p-value:

| band | `Pvalue` | `Pvalue_Burden` | `Pvalue_SKAT` |
|---|---|---|---|
| p > 0.05 | -0.0001 | -0.0002 | -0.0001 |
| 1e-5 .. 0.05 | +0.0003 | +0.0004 | +0.0004 |
| p <= 1e-5 | **+0.0153** (3.6%) | **+0.0218** (5.1%) | **+0.0464** (11.3%) |

The band split is a DECOMPOSITION, not a convenience. #158 touches the three
columns differently: the per-rho layer is never screened, so `Pvalue_Burden` and
`Pvalue_SKAT` above 1e-5 isolate the MODEL change with no kernel confound, while
`Pvalue` above 0.05 carries the screen as well. Below 1e-5 all three carry the
saddlepoint. So the `1e-5 .. 0.05` row measures the step-1 model change and says
it moves nothing (0.07%), and the bottom row measures the saddlepoint and lands
inside its independently documented 5-12% conservatism.

Genotype-level agreement on the stable traits is exact, not approximate: **MAC
ratio new/old median 1.000 with IQR [1.000, 1.000]** and **`BETA_Burden` sign
agreement 100.0%** on EUR T2Diab's 4,749 cells with p < 0.01 in either run.

**137 of 166 pairs have the same samples to within 0.5%.** The 29 that differ all
trace to one of four causes, none of them the pipeline:

- **Cases recoded, same people** (`dN` ~ 1.000, `dcase` != 1, and MAC ratio
  exactly 1.000 because MAC counts alleles over all samples regardless of case
  status): `MatHem` 2.44/3.23/2.72x more cases (eur/afr/amr), `CAD`
  1.67/1.59/1.55x, `ColonRectCanc` 1.21/1.56/1.34x.
- **Sample grew, cases unchanged** (`dcase` exactly 1.000, more controls): `PAD`
  +1.4 to +3.4%, `PeptUlcer` +1.0 to +1.4%, `Urolith` +0.6 to +1.0%. `eur PAD`'s
  `dN = 1.0335` matches its MAC ratio 1.0327 to three decimals.
- **The old model was mis-specified by the truncated step-0 GRM.** `amr BMI` is
  the extreme -- rho 0.30-0.75 and `med_log10` up to **6.4**, with 82,542 cells
  below 1e-5 in the OLD run, i.e. a quarter of the genome. This is the trait whose
  old null variance ratio was 6.91, so it is INDEPENDENT confirmation of the
  truncated-GRM diagnosis from the results side rather than the variance-ratio
  side. Milder cases of the same thing, all quantitative and all heritable, with
  stable sample counts: `amr AST` +0.125, `amr TChol` +0.140, `eur HDLC` +0.127,
  `eur Height` +0.103, `afr Height` +0.173. The damage scales with how much the
  GRM mattered.
- **Cosmetic**, 0.5-1% on cases or N and absent from the p-value outlier list:
  `Asth`, `Psori`, `HTN`, `CRP`.

Two findings about the OLD data that the comparison turned up.

**The old combined files carry DUPLICATE rows.** Found on eur T2Diab: two genes
appear twice with identical values in the top-hits list, and a pooled-row join
returned 20,402 rows from 18,463 genes. Probably a stitching artefact in
`combined_files`. The summariser de-duplicates on `(Region, Group, max_MAF)`
before merging, because otherwise the merge fans out and gene ranks are distorted
by repeated entries. The new per-bin outputs cannot carry this: the chromosome
bins partition 1-22 and that is asserted before every launch.

**`EFRMB` exists twice per ancestry in the old run, under one label, with
different case counts AND different totals.** Resolvable after all: for every
ancestry one of the two matches the new case count almost exactly and the other
is 2.0-2.5x off, so the canonical old file is the LARGER-case one (afr 4,710 vs
new 4,725; eur 12,234 vs 12,351; amr 3,867 vs 3,882; eas 399 vs 403; sas 186 vs
191). The smaller-case files are a superseded definition.

Three things to carry into any future comparison of these files. `MatHem`, `CAD`
and `ColonRectCanc` are NOT the same analysis in the two runs and must not be
compared as though they were -- `MatHem` in particular has 2.4-3.2x more cases and
is 6-8x LESS significant with sign agreement down to 0.909, which is a
phenotype-definition question about whether the broader definition dilutes the
rare-variant signal, not a pipeline one. `markers_same` tracks the SAMPLE set
rather than the annotation, because the MAF and MAC filters are computed in
sample, so a different sample set admits different markers per cell (eur T2Diab
0.967 with an identical sample set, eur PAD 0.344 with a 3.4% larger one). And
lambda_GC is unusable throughout, for the reason in the EUR T2Diab section.

**Still outstanding: the variant-level comparison.** The old directory holds only
`SAIGE.gene.*` files; whether single-variant results exist elsewhere was not
established. Worth doing if they do, because single-variant p-values come from the
score/SPA/ER/Firth cascade and are untouched by #158 -- so they would isolate the
model change with no kernel confound at all, more cleanly even than the
`1e-5 .. 0.05` band above. Note the new run has `FIRTH=FALSE`, so effect sizes
would differ from an old run with Firth on, and `--minMAC=0.5` sets the marker
floor.

## The three LEDGER rows

- **#161 CONFIRMED at production scale** on a build that already carries the fix.
  Two of 42 AMR binary fits print `tau2 = 0 by the floor` with Newton steps
  0.01295 and 0.01959, both inside the predicted band `4.08e-4 <= Dtau < tol`,
  and both leave zero at `--tol 0.01` exactly as the arithmetic says. Controls
  clean. Owner's decision: keep `tol = 0.02` across all ancestries.
- **#158 VALIDATED at production scale.** 200,652 `liu-screen` cells, every one
  above 0.05 and none below -- the screen's contract holds on real data, which
  had only ever been checked against this branch's synthetic nulls.
- **#162 RAISED.** The `skatoMethod` sidecar covers only the cells ASKED of the
  kernel (231,509), while the results file emits 329,393, and the 30% shortfall
  is structured entirely by max MAF. Provenance only; `Pvalue` is unaffected.

## What is NOT done

- **#162 is unexplained.** The cheap first step is empirical and on the row.
- **Zone attribution for the two new preemption batches** was not collected;
  `dstat --status FAILURE --full --format json` carries it.
- **The QQ rollout across all 190 trait-ancestry pairs** is built but was run
  only for EUR T2Diab. Render figures to the bucket rather than inline: browser
  output is the AoU egress, copying from a bucket is not, and 190 figures inline
  is 60-150 MB of it. The primary artefact should be a per-trait table whose key
  column is the `synonymous` calibration ratio.
- **The `tauInit` roxygen fix** from the #161 diagnosis is still outstanding
  (says `0,0` becomes 0.5 for binary; the code sets 0), plus its docdrift rule.
