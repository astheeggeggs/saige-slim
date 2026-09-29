# SAIGE parameters as run on All of Us (v4 production, 2026-09-2x)

The reference for porting the AoU runs to a template elsewhere. Every value
below names its source. **Step 2 is complete** (the wrapper
`saige_step2_brava_main_original.sh`, the notebook's `submit()`, and the
option dump SAIGE prints at the top of every job log, quoted from a FlexRV
run). **Steps 0 and 1 are partial**: their wrappers live in the workbench and
are not in this repo; what is here is what the handoffs recorded. Do not
treat a blank in those sections as "default".

## Images

| what | image |
|---|---|
| v4 production (both arms as scanned) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f2fb24da` |
| saddlepoint reruns (LEDGER #163-#170) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-3e92d89d` |
| the #172 fix (the notebook's `IMG9`) | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-a08fb295` |

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

## Step 1: null models (partial)

From `docs/state/handoff-2026-09-23-fit-gates-and-aou-refit.md` and `-09-25-`:

| flag | value | source |
|---|---|---|
| model directory | `<bucket>/data/saige/run/nullglmm_v4/<anc>_<TRAIT>.rda` (+ `.varianceRatio.txt`), ONE directory for afr, eas, eur and sas; amr in `<bucket>/data/saige/run/nullglmm_v4_rc0.125/`. Not under `SAIGE_SLIM_TESTING/`, which holds step-2 output only. File counts 42 / 25 / 43 / 15 / 42 for afr / eas / eur / sas / amr, read from the bucket 2026-09-29 | notebook `STEP1` (a dict: ancestry to directory and relatedness cutoff), printed 2026-09-29 |
| `--useSparseGRMtoFitNULL` | TRUE (sparse-GRM fits; single-threaded whatever `--nThreads` says) | handoff 09-23 |
| `--useSparseGRMforVarRatio` | TRUE: the AoU variance-ratio files carry `sparse` rows (`ratioVec_sparse 1 1` in every step-2 log), and step 2 given the GRM REFUSES a file without them | the step-2 option/VR dump |
| `--relatednessCutoff` | 0.05; amr 0.125 (must match step 2's, nothing checks it) | handoff 09-25 |
| `--tol` | 0.02, both trait types (the refit on `1.5.2-dev-3e61864a`) | handoff 09-23 |
| `--traitType` | binary / quantitative per trait (`TRAIT_TYPE` in the notebook) | notebook |
| `--invNormalize` | not recorded here | -- |
| `--covarColList`, `--qCovarColList`, `--sampleIDColinphenoFile`, `--phenoCol` | not recorded here (the step-1 wrapper) | -- |
| `--isCateVarianceRatio`, `--cateVarRatioMinMACVecExclude` | the VR files carry 2 `null` rows and 1 `sparse` row, consistent with `10,20.5` / `20.5` | VR files |

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
| `--is_single_in_groupTest` | TRUE |
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

## universal-saige (BRaVa-genetics, `main` at 8fa1aaf) against the AoU runs

The analyst pipeline the AoU wrapper descends from: drivers `00_step0_VR_and_GRM.sh`,
`01_step1_fitNULLGLMM.sh`, `02_step2_SPAtests_variant_and_gene.sh`, templates under
`templates/`, image `wzhou88/saige:1.3.6` pulled from Docker Hub by
`download_resources.sh`, one chromosome per step-2 invocation, run through
`run_container.sh`. Where a value below is marked "AoU ?" the AoU step-0/1
wrapper (workbench only) has not been read; the correspondence is presumed
from the GRM file names, not established.

| step | flag | universal-saige main | AoU v4 |
|---|---|---|---|
| 0 | GRM markers | LD-pruned `--indep-pairwise 50 5 0.05`, `createSparseGRM.R --numRandomMarkerforSparseKin=5000 --relatednessCutoff=0.05`, `--nThreads=$(nproc)` | file name agrees (`_relatednessCutoff_0.05_5000_randomMarkersUsed`); prefix `allofus_array_<pop>_snp_wise`; amr also at 0.125; AoU ? for the pruning |
| 0 | VR plink | 2,000 markers at 10 <= MAC < 20 plus 2,000 at MAC >= 20 from the merged array data | AoU ? |
| 1 | GRM | `--useSparseGRMtoFitNULL=TRUE --relatednessCutoff 0.05` | same at 0.05; amr 0.125 |
| 1 | `--tol` | 0.00001 quantitative, 0.02 binary | 0.02 BOTH (handoff 09-23) |
| 1 | `--invNormalize` | TRUE for quantitative | AoU ? |
| 1 | covariates | template `age,age2,age_sex,age2_sex,sex,PC1..PC10`, categorical `sex` | AoU ? |
| 1 | `--useSparseGRMforVarRatio` | NOT passed: the VR file has no `sparse` rows, and step 2 with the GRM refuses it on the slim build (found running the drivers end to end, 2026-09-28) | TRUE (the files carry `sparse` rows) |
| 1 | other | `--skipVarianceRatioEstimation FALSE --IsOverwriteVarianceRatioFile=TRUE --isCateVarianceRatio=TRUE --isCovariateOffset TRUE --nThreads=$(nproc)-1`, `--SampleIDIncludeFile` | AoU ?; a sparse-GRM fit is single-threaded regardless |
| 2 | image | `wzhou88/saige:1.3.6` (Docker Hub) | `saige-slim:1.5.2-dev-f2fb24da` / `-3e92d89d` (gcr) |
| 2 | genotypes | PLINK bed or VCF (`--vcfField=DS`) | PGEN; the slim build has NO VCF reader |
| 2 | `--minMAC` | 0.5 for both test types | 0.5 group, 4 variant |
| 2 | `--is_Firth_beta`, `--pCutoffforFirth` | TRUE, 0.1 | FALSE |
| 2 | `--is_fastTest` | TRUE | FALSE |
| 2 | `--maxMissing`, `--impute_method` | 1, mean | build defaults 0.15, best_guess |
| 2 | `--relatednessCutoff` | not passed (moot when steps 0 and 1 share the cutoff; wrong for amr at 0.125) | passed, equal to step 1's, required by the wrapper |
| 2 | annotations | `pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:...` (BRaVa labels) | `pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous` (the AoU annotation's labels) |
| 2 | `--maxMAF_in_groupTest` | `0.0001,0.001,0.01` | same (SKAT-O arm); FlexRV arm 0.001 alone |
| 2 | FlexRV | absent | the FlexRV arm (flags above) |
| 2 | `--condition`, `--subSampleFile` | supported | not used |
| 2 | fit gates | absent from 1.3.6 | ON in the slim build (`--minCaseCount` 100 etc.): a small cohort's trait is REFUSED, the wrapper records it as `REFUSED_BY_GATE` |
| 2 | sidecars | none | pooledTests, skatoMethod, skatFailures, spaFallbacks, chrStatus, spaDump |
