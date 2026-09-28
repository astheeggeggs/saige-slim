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
| model directory | `run/nullglmm_v4/<anc>_<TRAIT>.rda` (+ `.varianceRatio.txt`); amr in `nullglmm_v4_rc0.125/` | notebook `STEP1` |
| `--useSparseGRMtoFitNULL` | TRUE (sparse-GRM fits; single-threaded whatever `--nThreads` says) | handoff 09-23 |
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
