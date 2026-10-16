# SAIGE parameters as run on All of Us (v4 production, 2026-09-2x)

The reference for porting the AoU runs to a template elsewhere. Every value
below names its source. **Step 2 is complete** (the wrapper
`saige_step2_brava_main_original.sh`, the notebook's `submit()`, and the
option dump SAIGE prints at the top of every job log, quoted from a FlexRV
run). **Step 1 is complete** (the wrapper `saige_step1.sh` and three of its
job logs, read 2026-10-07). **Step 0 is partial**: its wrapper lives in the
workbench and is not in this repo; what is here is what the handoffs recorded.
Do not treat a blank in that section as "default".

## Images

| what | image |
|---|---|
| v4 production (both arms as scanned) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f2fb24da` |
| saddlepoint reruns (LEDGER #163-#170) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-3e92d89d` |
| the #172 fix (the notebook's `IMG9`) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-a08fb295` |
| the #176 stretch-gate fix (the notebook's `IMG12`) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-152ffd8c` |
| LEDGER #179 (the notebook's `IMG13`, twelve commits after IMG12: #179, #171, the sidecar/tool fixes; also `astheeggeggs/saige-slim:1.5.2-dev-8bffa0ca`) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-8bffa0ca` |
| `--is_output_log10P` (the notebook's `IMG14`, 1c3f89b5 + the pkgdown index fix; also `astheeggeggs/saige-slim:1.5.2-dev-524a299f`) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-524a299f` |
| **IMG16: LEDGER #181** (3a6d81fa; the region path's R-side saddlepoint NaN killed HTN chromosomes on IMG15 -- kernels refuse a non-finite saddle, the rescale's caller counts/names/dumps; `docs/state/handoff-2026-10-16-region-spa-na.md`; also `astheeggeggs/saige-slim:1.5.2-dev-3a6d81fa`) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-3a6d81fa` |
| **IMG15: the log10 fix** (f4a316d5; IMG14 LOST A CHROMOSOME with `LOG10P=TRUE`, `docs/state/handoff-2026-10-14-log10p.md` addendum; also `astheeggeggs/saige-slim:1.5.2-dev-f4a316d5`) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f4a316d5` |

The tags name the commits the images were built from; a hash that no longer
resolves in the current history still names its image.

## Step 0: sparse GRM (partial)

Read off the GRM file names and `docs/state/handoff-2026-09-25-aou-production-step1-step2.md`:

| flag | value | source |
|---|---|---|
| input | array genotypes, per ancestry (`allofus_array_<pop>`) | file name |
| `--relatednessCutoff` | 0.05 (all); amr also rebuilt at 0.125, the production amr models use that | file name; handoff 09-25 |
| `--numRandomMarkerforSparseKin` | 5000 | file name (`5000_randomMarkersUsed`) |
| marker selection | SNP-wise (`snp_wise` in the name) | file name |
| output | `run/make_sparse_grm/allofus_array_<pop>_snp_wise_relatednessCutoff_<c>_5000_randomMarkersUsed.sparseGRM.mtx` (+ `.sampleIDs.txt`) | notebook `GRM` |

The `step0/` copies of the GRMs are TRUNCATED for four of five cohorts and
must not be used (memory note `aou-sparse-grm-two-locations`).

## Step 1: null models (complete)

From the wrapper `saige_step1.sh` and three of its v4 job logs (2026-09-24,
each the wrapper's `=== config ===` block and SAIGE's dump of every option as
parsed): afr/ALT (quantitative), amr/IFHern (binary) and amr/FemInf (binary,
female-only). The two amr logs are the 0.05 jobs run before amr moved to 0.125.

| flag | value | source |
|---|---|---|
| model directory | `<bucket>/data/saige/run/nullglmm_v4/<anc>_<TRAIT>.rda` (+ `.varianceRatio.txt`), ONE directory for afr, eas, eur and sas; amr in `<bucket>/data/saige/run/nullglmm_v4_rc0.125/`. Not under `SAIGE_SLIM_TESTING/`, which holds step-2 output only. File counts 42 / 25 / 43 / 15 / 42 for afr / eas / eur / sas / amr, read from the bucket 2026-09-29 | notebook `STEP1` (a dict: ancestry to directory and relatedness cutoff), printed 2026-09-29 |
| genotypes (variance ratio) | `--bedFile/--bimFile/--famFile` `data/saige/run/nullglmm/allofus_array_<anc>_for_vr.{bed,bim,fam}` | logs |
| sparse GRM | `data/saige/run/make_sparse_grm/allofus_array_<anc>_snp_wise_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx` (+ `.sampleIDs.txt`); NOT the truncated `data/saige/step0/` copies | logs |
| `--useSparseGRMtoFitNULL`, `--useSparseGRMforVarRatio` | TRUE, TRUE (the variance-ratio files carry `sparse` rows; step 2 given the GRM REFUSES a file without them) | logs; the step-2 option/VR dump |
| `--relatednessCutoff` | 0.05; amr 0.125 (must match step 2's, nothing checks it) | logs (0.05); amr: handoff 09-25, and the rerun at 0.125 confirmed by the analyst 2026-10-07 |
| `--phenoFile` | `data/phenotypes/brava_pilot_v1.collapsed.with_covariates.txt.gz` | logs |
| `--sampleIDColinphenoFile` | `eid` | logs |
| `--traitType` | binary / quantitative per trait (`TRAIT_TYPE` in the notebook) | logs; notebook |
| `--invNormalize` | TRUE for quantitative, FALSE for binary (the wrapper derives it from the trait type) | logs (ALT TRUE; IFHern, FemInf FALSE); wrapper |
| `--covarColList` | `age,sex,age2,age_sex,age2_sex,PC1,...,PC20` (25; 26 design columns with the intercept); sex-specific traits `age,age2,PC1,...,PC20` | logs (ALT, IFHern; FemInf) |
| `--qCovarColList` | `sex`; none for sex-specific traits | logs |
| sex-specific traits | `--FemaleOnly=TRUE` / `--MaleOnly=TRUE`, `--sexCol=sex`, `--FemaleCode=1` / `--MaleCode=0` (the wrapper's `IN_SEXMODE`); in this file 1 is female | log (FemInf); wrapper |
| `--tol` | 0.02, both trait types | logs; handoff 09-23 |
| `--LOCO` | FALSE | logs |
| `--isCateVarianceRatio` | TRUE, with `cateVarRatioMinMACVecExclude` 10,20.5 and `cateVarRatioMaxMACVecInclude` 20.5 (defaults) | logs |
| `--isCovariateOffset`, `--IsOverwriteVarianceRatioFile` | TRUE, TRUE | logs |
| `--usePCGwithSparseGRM` | FALSE (the default; these jobs' wrapper did not pass it, the later one passes `IN_USE_PCG`, default FALSE) | logs; wrapper |
| `--minCaseCount`, `--minCovariateCellCount`, `--minSampleCount` | 100, 10, 0 | logs |
| `--separationCheck`, `--convergenceCheck`, `--separationMaxFraction` | error, error, 0.02 (defaults) | logs |
| `--nThreads` | cores - 1 (1 or 3 here); inert, a sparse-GRM fit is single-threaded | logs |
| everything else | build defaults, as the option dump shows | logs |

## Step 2: association tests (complete)

Common to both arms (the wrapper's fixed flags plus the notebook's `submit()` values):

| flag | value |
|---|---|
| genotypes | PGEN, the coding subsets `plink_coding/aou.exome_split.v8.qced.chr<c>.<anc>.af01.popmax01.cod.{pgen,pvar,psam}` |
| `--chrom` | `chr<c>`, one SAIGE process per chromosome, results concatenated |
| `--minMAF` | 0 |
| `--minMAC` | 0.5 (group tests); 4 (variant tests) |
| `--LOCO` | FALSE |
| `--is_output_moreDetails` | TRUE |
| `--is_fastTest` | FALSE (`FASTTEST="FALSE"` in the launch; the wrapper's `auto` would give TRUE for variant tests only) |
| `--is_Firth_beta` | FALSE (`FIRTH="FALSE"` in the launch) |
| `--GMMATmodelFile`, `--varianceRatioFile` | the step-1 pair above |
| `--sparseGRMFile`, `--sparseGRMSampleIDFile` | the step-0 GRM of the ancestry, every ancestry (the earlier "none for amr" is void) |
| `--relatednessCutoff` | 0.05; amr 0.125 (`RELCUTOFF`, required by the wrapper when a GRM is given) |
| `--is_single_in_groupTest` | TRUE in v4. **From 2026-10-08 the wrapper takes `SINGLE_IN_GROUP` (default TRUE); set FALSE for FlexRV runs**: the file was 80% of a FlexRV job (LEDGER #179; afr chr20 431/465/601 s -> 55/55/75 s with it off, gene-level rows identical), its real-marker rows are the SKAT-O run's and its pseudo-marker rows were read only by the splice control's twin rule. The wrapper then delivers a one-line placeholder; `aou_spa_scan.sh` and `aou_spa_splice.sh` recognise it. SKAT-O runs keep TRUE. |
| `--is_output_log10P` | absent in v4 (the flag is 2026-10-14). **From 2026-10-15 the wrapper takes `LOG10P` (default TRUE)**: `log10_Pvalue*` columns after `SE_Burden` on both region paths and the pooled row combined in log space, so a gene below the double floor (the top ALT gene of the eur chr20 pilot at 1.48e-323; GPT, the ALT enzyme gene, reaches 1e-146 genome-wide) is estimated, not 0. Set FALSE for a partial rerun that the splice control will compare against a v4 file (its header check refuses the extra columns). Needs IMG15 (IMG14 lost a chromosome with the flag on). |
| `--is_output_markerList_in_groupTest` | TRUE |
| `--subSampleFile`, `--minGroupMAC_in_BurdenTest`, `--is_no_weight_in_groupTest` | not passed (build defaults: none, 1, FALSE) |
| machine | `n2-standard-2` (eur `n2-highmem-2`), on demand, 50 GB disk; chromosomes in size-balanced bins per job |

SKAT-O arm (`skv4-<anc>-<trait>[-b<i>]`, output `<anc>/skato_v4/`):

| flag | value |
|---|---|
| `--groupFile` | per-ancestry `annotations/aou.v8.chr<c>.<anc>.variants_only.spliceai=0.20_cadd=28.1_revel=0.773.canonical.txt.gz` |
| `--annotation_in_groupTest` | `pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous` |
| `--maxMAF_in_groupTest` | `0.0001,0.001,0.01` |
| `--r.corr` | build default (0: SKAT-O) |

FlexRV arm (`fxv4-<anc>-<trait>[-b<i>]`, output `<anc>/flexrv_v4/`, fixed copies `flexrv_v4_fixed/`):

| flag | value |
|---|---|
| `--groupFile` | pan-ancestry `annotations/aou.v8.chr<c>.variants_only.spliceai=0.20_cadd=28.1_revel=0.773.canonical.flexrv_AM.txt.gz` (carries the `score:AM` line) |
| `--annotation_in_groupTest` | `pLoF:damaging_missense:other_missense` (ONE group) |
| `--maxMAF_in_groupTest` | `0.001` (ONE max MAF, equal to `--flexRV_maxMAF`) |
| `--flexRV_score` | `AM` |
| `--flexRV_maxMAF` | 0.001 |
| `--flexRV_lofAnno` | `pLoF` |
| `--flexRV_scoreTransforms`, `--flexRV_MAFTransforms` | not passed: the paper's `unw;lof;cbrt:0.2,...,0.9;sbeta:1.5,2,3,5,10,20` and `unw;mafbeta:2,5;mac:1,3,5,10,25,50,100,200,500` (rungs at or above `2 N maxMAF` dropped per LEDGER #142) |
| `--r.corr` | 1 (burden) |

Every other option at the build's default, as the option dump of a FlexRV job
shows (sas/BenIntNeo, chr5, `1.5.2-dev-4ff9576b`): `SKATO_screen_threshold`
0.05, `SKATO_tail_method` saddlepoint, `SKATO_tail_threshold` 1e-05,
`minGroupMAC_in_BurdenTest` 1, `maxMissing` 0.15, `impute_method` best_guess,
`markers_per_chunk` 10000, `groups_per_chunk` 100,
`MACCutoff_to_CollapseUltraRare` 10, `cateVarRatioMinMACVecExclude` 10,20.5,
`cateVarRatioMaxMACVecInclude` 20.5, `weights.beta` 1,25,
`markers_per_chunk_in_groupTest` 100, `SPAcutoff` 2, `dosage_zerod_cutoff`
0.2, `dosage_zerod_MAC_cutoff` 10, `minCaseCount` 100, `separationCheck` error,
`convergenceCheck` error, `separationMaxFraction` 0.02, `pCutoffforFirth` 0.01
(inert with Firth off), `firth_solver` closed, `is_noadjCov` TRUE,
`max_MAC_for_ER` 4, `nThreads` 1, `AlleleOrder` alt-first, `is_imputed_data`
FALSE, `maxMAC_in_groupTest` 0, `is_output_weight_column` FALSE.

## Where the exact record of any single run is

The first screen of every step-2 job log is SAIGE's dump of every option as
parsed (`$bedFile ... $nThreads`), so a question about one run is answered by
`fetch_log(<job>)` and reading that block, not by this page.

## universal-saige against the AoU runs

The analyst pipeline the AoU wrapper descends from (BRaVa-genetics/universal-saige).
Its branch `saige-slim-aou-parameters` (PR #7) moved it onto this build and these
parameters; `main` (8fa1aaf) is the earlier pipeline (`wzhou88/saige:1.3.6`, Firth
and fastTest on, `--tol` 1e-5 for quantitative traits, no `sparse` variance-ratio
rows, VCF input). The branch, as of 2026-10-07:

| step | flag | universal-saige `saige-slim-aou-parameters` | AoU v4 |
|---|---|---|---|
| all | image | `astheeggeggs/saige-slim:1.5.2-dev-152ffd8c` (Docker Hub) | `f2fb24da` / `3e92d89d` / `a08fb295` (gcr) |
| 0 | GRM | LD-pruned `--indep-pairwise 50 5 0.05`, `createSparseGRM.R --numRandomMarkerforSparseKin=5000`, `--relatednessCutoff` a flag (default 0.05) | file names agree; AoU pruning not recorded |
| 0 | VR plink | 2,000 markers at 10 <= MAC < 20 plus 2,000 at MAC >= 20 | not recorded |
| 1 | GRM, VR | `--useSparseGRMtoFitNULL=TRUE --useSparseGRMforVarRatio=TRUE`, `--relatednessCutoff` a flag | same |
| 1 | `--tol`, `--invNormalize` | 0.02; TRUE quantitative, FALSE binary | same |
| 1 | covariates | the analyst's; templates `age,age2,age_sex,age2_sex,sex,PC1..PC20`, categorical `sex` | same set |
| 1 | sex-specific | `--sex F/M` -> `--FemaleOnly/--MaleOnly --sexCol sex`; numeric codes default 0 female, 1 male (BRaVa's coding, the opposite of this file's), checked against genetic sex | `--FemaleOnly/--MaleOnly`, `--FemaleCode=1` |
| 1 | dense GRM | refused above 100 relatives per sample (amr at 0.05: ~644) unless `--forceDenseGRM` | amr refitted at 0.125 |
| 2 | genotypes | PGEN or PLINK 1; no VCF | PGEN |
| 2 | `--minMAC`, `--minMAF` | 4 variant, 0.5 group; 0 | same |
| 2 | Firth, fastTest, `--LOCO` | FALSE, FALSE, FALSE | same (`FIRTH=FALSE`, `FASTTEST=FALSE` in the launch) |
| 2 | `--relatednessCutoff` | passed (default 0.05) | passed, equal to step 1's |
| 2 | `--maxMAF_in_groupTest` | `0.0001,0.001,0.01`; FlexRV 0.001 | same |
| 2 | FlexRV | one mask, `--flexRV_maxMAF` 0.001, `--flexRV_lofAnno pLoF`, `--r.corr=1`; a score line in every region required | same |
| 2 | annotations | BRaVa labels; every label must be in the group file | AoU labels |
| 3 | Nglmm | `extractNglmm.R` at `--relatednessCutoff` (default 0.05); before 2026-10 its own default, 0.125 | not run |

