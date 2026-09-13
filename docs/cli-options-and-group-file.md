# CLI options and group-file format on the `strip` branch

Inventoried 2026-09-12 by tracing each option from its `make_option()`
declaration to the code that consumes it, with every "not live" claim then
checked adversarially (4 statuses were corrected that way). **The status column
is the point of this document** -- the help strings cannot tell you which
options are refused, inert, trait-gated or path-gated, and 7 of them say
something false. Spot-checked by hand afterwards: `--X_PARregion`,
`--is_output_moreDetails` and `--pval_cutoff_for_fastTest` all confirmed.

# SAIGE `strip` — implemented CLI options

**Shape.** 139 distinct options across three wrappers: **step 0** `createSparseGRM.R` 16, **step 1** `step1_fitNULLGLMM.R` 58, **step 2** `step2_SPAtests.R` 65. Step 0 builds the sparse GRM (`.mtx` + `.sampleIDs.txt`, both named from `outputPrefix`/`relatednessCutoff`/`numRandomMarkerforSparseKin`). Step 1 fits the null GLMM and estimates variance ratios (one or several phenotypes per process). Step 2 tests markers, or regions when `--groupFile` is given — the single branch point at `R/SAIGE_Test_main.R:367` that decides whether ~25 of step 2's options do anything at all. PLINK bed/bim/fam and PGEN pgen/pvar/psam are the only genotype inputs at every step.

Legend: **L** live · **P** partly live · **R** refused · **I** inert · `+` this branch's option · `~` upstream's option, meaning changed here.

---

## 2. STEP 2

### Genotype input (12)

| option | default | what it changes | |
|---|---|---|---|
| `bedFile` | "" | PLINK .bed; non-empty **selects the PLINK reader** | L |
| `bimFile` | "" (from bed) | marker universe (CHROM/ID/POS) + REF/ALT cols; cannot itself select PLINK | L |
| `famFile` | "" (from bed) | .bed sample row order → model-to-column map | L |
| `pgenFile` `+` | "" | .pgen; non-empty **selects the PGEN reader** | L |
| `pvarFile` `+` | "" (from pgen) | marker universe on the PGEN path | L |
| `psamFile` `+` | "" (from pgen) | .pgen sample row order (IID read in C++) | L |
| `AlleleOrder` `~` | alt-first | which bim column is counted as ALT | **P** |
| `idstoIncludeFile` `~` | "" | drop every marker row not listed (rs IDs here) | L |
| `rangestoIncludeFile` `~` | "" | drop rows outside the CHROM/POS intervals (inclusive) | L |
| `chrom` `~` | "" | the one chromosome tested; also picks the LOCO fit | L |
| `subSampleFile` | "" | subset every per-sample vector of the null model | L |
| `is_imputed_data` `~` | FALSE | would swap MissingRate → imputationInfo | **R** |

### Variant QC (7)

| option | default | what it changes | |
|---|---|---|---|
| `minMAF` | 0 | drop marker below this in-sample MAF | L (no test/golden ≠ 0) |
| `minMAC` | 0.5 | drop marker below this MAC (now from integer counts, #29) | L |
| `minInfo` | 0 | INFO filter — no reader supplies INFO | **R** if > 0 |
| `maxMissing` | 0.15 | drop marker whose missing rate *exceeds* this (strict `>`) | L |
| `impute_method` | best_guess | value written into a missing call: round(2AF) / 2AF / 0 | L |
| `dosage_zerod_cutoff` | 0.2 | zero every genotype ≤ this, for markers ≤ MAC cutoff | L |
| `dosage_zerod_MAC_cutoff` | 10 | the MAC ceiling that zeroing applies under | L |

### Model in / out (10)

| option | default | what it changes | |
|---|---|---|---|
| `GMMATmodelFile` | "" | the step-1 `.rda` (must hold `modglmm`) | L |
| `varianceRatioFile` | "" | MAC-keyed variance ratios rescaling the score variance | L (see §4) |
| `GMMATmodel_varianceRatio_multiTraits_File` `+` | "" | 3-col manifest: N models, one genotype pass, one file each | L |
| `SAIGEOutputFile` `~` | "" | results path — and the prefix for `.index`, `.singleAssoc.txt`, `.markerList.txt`, `.skatFailures.txt`, `.pooledTests.txt`, `_<model>` | L |
| `is_overwrite_output` `~` | TRUE | TRUE truncates; FALSE resumes from the `.index`, truncating the results file (and, on a group run, `.singleAssoc.txt`) to the sizes the index recorded for the last completed chunk, on both the marker path (LEDGER #47) and the region path (LEDGER #140, fixed 2026-09-13). The end-of-run `.pooledTests.txt` / `.skatFailures.txt` sidecars are rewritten over the resumed regions only | L |
| `is_output_moreDetails` | FALSE | adds N_case/ctrl hom/het columns | **P** (binary only) |
| `LOCO` | TRUE | swap fitted values/residuals/obj.noK for the per-chrom fit | L |
| `sparseGRMFile` `~` | "" | score test against the sparse Σ (PCG) instead of VR-only | L |
| `sparseGRMSampleIDFile` | "" | GRM row order → subset/reorder into model order | L (no-op alone) |
| `relatednessCutoff` | 0 | drop GRM entries below this before Σ | L (no-op alone) |

### Region tests (11)

| option | default | what it changes | |
|---|---|---|---|
| `groupFile` `~` | "" | **switches step 2 to region testing**; 2+K lines, rs IDs | L |
| `annotation_in_groupTest` | lof,missense;lof,missense;lof;synonymous | the annotation masks = the `Group` column | L |
| `maxMAF_in_groupTest` | 0.0001,0.001,0.01 (API: 0.01) | the max_MAF cells; its max is also the ingest cutoff | L |
| `maxMAC_in_groupTest` `~` | 0 | extra cells expressed as allele counts (÷ 2N of the model) | L (was wrong upstream, #61) |
| `MACCutoff_to_CollapseUltraRare` | 10 | MAC ≤ this collapses into one per-cell pseudo-marker | L |
| `minGroupMAC_in_BurdenTest` `~` | **1** (upstream 5) | refuse a cell below this unweighted MAC, both paths | L (see §4/§6) |
| `r.corr` | 0 | 0 = SKAT-O (+SKAT, Burden), 1 = Burden only; nothing else | L |
| `is_single_in_groupTest` `~` | FALSE | write `.singleAssoc.txt` | **P** (forced TRUE on r.corr=0) |
| `is_output_markerList_in_groupTest` | FALSE | `.markerList.txt`: tested vs collapsed per cell | L (sidecar only) |
| `groups_per_chunk` | 100 | regions per flush: `.index` granularity, RAM, SKAT-O batch | L (no result bytes) |
| `markers_per_chunk_in_groupTest` `~` | 100 | block size of the var-cov assembly (RAM knob since RegionChunkStore) | L (no result bytes; see §4) |

### Weighting and conditioning (7)

| option | default | what it changes | |
|---|---|---|---|
| `weights.beta` `~` | 1,25 | beta(a,b) MAF weights; **semicolon pair list** = one beta set per extra pair | **P** (region only) |
| `is_no_weight_in_groupTest` `~` | FALSE | all-ones weights, or an appended `No_Weights` set | L (region only) |
| `is_output_weight_column` `+` | FALSE | forces the multi-weight schema at K=1 (Weight, ACATV, ACATO, pooled row) | L (region only) |
| `condition` `~` | "" | second set of `_c` columns with those markers projected out | L (both paths) |
| `weights_for_condition` | NULL | explicit weights for the conditioning markers | **P** (binary region only) |
| `cateVarRatioMinMACVecExclude` `~` | 10,20.5 | exclusive lower MAC edge of each VR bin | **P** |
| `cateVarRatioMaxMACVecInclude` `~` | 20.5 | inclusive upper edge (step 2 appends N as the last) | **P** |

### FlexRV (5, all `+`)

| option | default | what it changes | |
|---|---|---|---|
| `flexRV_score` | NULL | master switch; must equal the group file's `score:NAME` | L (region only) |
| `flexRV_scoreTransforms` | unw;lof;cbrt:…;sbeta:… (16) | score axis of the weight-set grid | L |
| `flexRV_MAFTransforms` | unw;mafbeta:2,5;mac:1…500 (12) | MAF/MAC axis (mafbeta decreasing, 0 at M; mac is 0/1) | L |
| `flexRV_lofAnno` | lof | which anno labels the `lof` transform keys on | L (only if `lof` in the score spec) |
| `flexRV_maxMAF` | 0.001 | the single max-MAF cell **and** the mafbeta boundary M | L |

### Test method (8)

| option | default | what it changes | |
|---|---|---|---|
| `SPAcutoff` | 2 | |Z| above which SPA runs instead of the plain score p | **P** (binary only) |
| `is_Firth_beta` | FALSE | re-fit BETA/SE by penalised Firth | **P** (binary only) |
| `pCutoffforFirth` | 0.01 | p at or below which Firth fires | **P** (binary + Firth only) |
| `firth_solver` `+` | closed | closed-form 2×2 hat values vs the 1.x QR | **P** (result-inert, −27% time) |
| `is_fastTest` | FALSE | cheap first pass, accurate re-test of what clears the cutoff | **P** (markers: results; regions: cost only) |
| `pval_cutoff_for_fastTest` `+` | 0.05 | the escalation threshold (per marker, or per region's Cauchy p) | L (help says otherwise — §4) |
| `is_noadjCov` | TRUE | TRUE = raw AF-centred score + null_noXadj ratio; FALSE regresses covariates out first | **P** (ignored by group tests) |
| `max_MAC_for_ER` | 4 | MAC ≤ this gets the exact resampling p | **P** (binary only) |

### Performance / chrX (5)

| option | default | what it changes | |
|---|---|---|---|
| `nThreads` `~` | 1 | in-process threads: marker loop; region marker loop + per-cell tests + SKAT-O fork pool | L (byte-identical) |
| `markers_per_chunk` | 10000 | markers per C++ call / flush / `.index` line | **P** (no result bytes; ignored on group runs) |
| `sampleFile_male` | "" | male IDs for the X doubling | **P** (dead unless the flag below) |
| `X_PARregion` | **""** (API: hg19 windows) | positions exempt from the doubling | **P** (dead unless below; see §6) |
| `is_rewrite_XnonPAR_forMales` `~` | FALSE | double male genotypes on non-PAR chrX (chrX-only since #69) | L |

---

## 3. STEP 1 (58) and STEP 0 (16)

### Step 1

| options | default | what it changes | |
|---|---|---|---|
| `plinkFile` / `bedFile` / `bimFile` / `famFile` | "" | PLINK GRM input; prefix fills blanks only; any one selects PLINK | L |
| `plink2File` `+` / `pgenFile` `+` / `pvarFile` `+` / `psamFile` `+` | "" | PGEN input (re-encoded to PLINK-1 records → byte-identical fits) | L |
| `phenoFile`, `sampleIDColinphenoFile` | "", IID | phenotype table; ID column read as character | L |
| `phenoCol` `~` | "" | the phenotype — **comma list = several fits in one process** + `.multiTraits.txt` | L |
| `traitType` `~` | binary | family + fitter; comma list per phenotype | L |
| `invNormalize` `~` | "FALSE" | rank-INT the phenotype (quantitative only; type is now character) | L |
| `covarColList`, `qCovarColList` | "" | fixed effects; which of them are factors | L |
| `SampleIDIncludeFile` | "" | restrict the fit to these IDs (after complete-casing) | L |
| `minCovariateCount` | −1 | drop a sparse binary covariate | **P** |
| `isCovariateTransform` | TRUE | QR-orthonormalise the design, back-transform coefficients | L |
| `isCovariateOffset` | TRUE | hold covariate effects fixed as an offset | L |
| `seed` `+` | 1 | RNG seed; re-seeded before every fit | L |
| `tol` / `maxiter` | 0.02 / 20 | outer τ convergence (also the τ^−2 divergence ceiling) and both loop caps | L |
| `tolCoef` `+` | 1e-5 | inner (α,b) tolerance, decoupled from `tol` (#125) | **P** (binary only) |
| `tolPCG` / `maxiterPCG` | 1e-5 / 500 | PCG residual (compared as a **sum of squares**) and cap | L |
| `nThreads` `~` | 1 | TBB pool for the GRM kernels | L |
| `SPAcutoff` | 2 | nothing in step 1 | **I** |
| `nrun` / `traceCVcutoff` | 30 / 0.0025 | starting probe count (grows by 10) and its CV target | L |
| `memoryChunk` | 2 | nothing | **I** |
| `tauInit` | "0,0" | starting variance components (binary reads only the 2nd) | L |
| `skipModelFitting` | FALSE | re-estimate VRs from an existing `.rda` only | L |
| `skipVarianceRatioEstimation` | FALSE | write no `.varianceRatio.txt` | L |
| `numRandomMarkerforVarianceRatio` / `ratioCVcutoff` | 30 / 0.001 | VR marker pool start (grows by 10) and CV target | L |
| `isCateVarianceRatio` | FALSE | one VR per MAC category | L |
| `cateVarRatioMinMACVecExclude` / `MaxMACVecInclude` | 10,20.5 / 20.5 | the category edges | **P** (both) |
| `useSparseGRMtoFitNULL` | FALSE | fit against the sparse GRM (takes over LOCO, nThreads, VR skip) | L |
| `usePCGwithSparseGRM` | FALSE | iterative instead of direct sparse solve | L (sparse fits only) |
| `useSparseGRMforVarRatio` | FALSE | **the only source of `sparse` VR rows** | L |
| `sparseGRMFile` / `sparseGRMSampleIDFile` | "" | the step-0 `.mtx` and its ID order | L |
| `relatednessCutoff` | 0 | `drop0` tolerance on the loaded GRM | L (0 = no-op; API default 0.125) |
| `minMAFforGRM` / `maxMissingRateforGRM` | 0.01 / 0.15 | marker QC for the GRM (and the GRM denominator) | L |
| `isDiagofKinSetAsOne` | FALSE | GRM diagonal = 1 instead of the self-product | L |
| `includeNonautoMarkersforVarRatio` | FALSE | admit parseable non-autosomes to the VR pool | L |
| `LOCO` `~` | **NULL** → `!useSparseGRMtoFitNULL` | per-chromosome refits | L |
| `isLowMemLOCO` | FALSE | `_chr<j>.rda` files + `_noLOCO.rda` instead of one `.rda` | L |
| `sexCol` / `FemaleCode` / `MaleCode` / `FemaleOnly` / `MaleOnly` | "" / 1 / 0 / F / F | sex-specific fit; suffixes every output name | L |
| `outputPrefix` `~` / `outputPrefix_varRatio` / `IsOverwriteVarianceRatioFile` | "~/" / "" / FALSE | output stems; permit clobbering the VR file | L |

### Step 0

| options | default | what it changes | |
|---|---|---|---|
| `plinkFile` `~` / `bedFile` `~` / `bimFile` `~` / `famFile` `~` | "" | PLINK input; prefix now fills blanks only (upstream ignored explicit paths) | L |
| `plink2File` `+` / `pgenFile` `+` / `pvarFile` `+` / `psamFile` `+` | "" | PGEN input; dosages hard-called with a counted warning | L |
| `nThreads` `~` | **16** | TBB pool size (pinned unconditionally, #7) | **P** |
| `memoryChunk` | 2 | nothing | **I** |
| `outputPrefix` | "~/" | stem of both output files | L |
| `numRandomMarkerforSparseKin` | 2000 | markers in the cheap kinship **screen** (a power knob) | L |
| `relatednessCutoff` | 0.125 | applied twice — to the screen and to the full-marker value | L |
| `isDiagofKinSetAsOne` | FALSE | diagonal values (count is always N) | L |
| `minMAFforGRM` / `maxMissingRateforGRM` | 0.01 / 0.15 | marker QC — and hence **which** random markers get drawn | L |

Step 0 reads only the **line count** from `bimFile`/`pvarFile`; IDs, chromosomes and alleles are never parsed there.

---

## 4. NOT FULLY LIVE

**Refused outright (step 2, both are LEDGER #107/#108, kept declared for a returning BGEN/VCF reader):**
- `--minInfo` > 0 → stop, "no genotype reader on this branch supplies an imputation INFO score". Both readers hardcode `imputeInfo = 1`, so the four surviving C++ gates can only ever evaluate `1 < 0`.
- `--is_imputed_data=TRUE` → stop, "that column could only ever read 1". The four writer either/ors in `Main.cpp` are now unreachable dead surface; golden P48 retired.

**Inert (3, all inert upstream too):**
- step 1 `--SPAcutoff` — reaches the `fitNULLGLMM` formal and is never forced; not even range-checked here.
- step 0 / step 1 `--memoryChunk` — the one read is commented out at `SAIGE_fitGLMM_fast.cpp:1054`, followed by `numMarkersofEachArray = 1`. One array per marker regardless. Only observable effect: it is part of the multi-phenotype ingest key, so changing it between phenotypes forces a redundant re-ingest.
- step 0 `--nThreads`'s `nblocks`/`ncore` arguments — `createSparseKinParallel` never reads either; the blocking function that would have is never called. The pool size itself is live.

**Trait-gated (do nothing, silently, on the wrong trait):**
- quantitative: `--SPAcutoff`, `--is_Firth_beta`, `--pCutoffforFirth`, `--firth_solver`, `--max_MAC_for_ER` (every C++ site is guarded on `traitType == "binary"`); `--is_output_moreDetails` produces **no columns at all** — proved by golden P32, generated with the flag, whose header still ends at `N`.
- binary: step 1 `--invNormalize` (quantitative only), `--minCovariateCount` (its call site is `traitType == "binary"` — the help does not say so, and the test is strict `<`, so −1 *and* 0 exclude nothing); step 1 `--tolCoef` (the quantitative fitter takes the argument and drops it at all four call sites, falling back to `Get_Coef`'s 0.1).

**Path-gated (live on one path, unread on another):**
- `--AlleleOrder`: PLINK only. The pgen spec carries no `AlleleOrder` key, so it is silently ignored on PGEN — where upstream *required* `ref-first`. No validation of the string exists anywhere, so a typo fires neither branch.
- `--is_noadjCov`: ignored by every group test (`Main.cpp:2800/3657/4112/4715` hardcode `false`) and by a sparse-GRM marker run unless `--is_fastTest` is on.
- `--is_fastTest`: changes marker results (P28/P29); on regions the escalated pass recomputes what the fast pass had — result-inert at +16.3% wall time. Four R conditions force it FALSE silently (no VR file, cutoff 0, `!sparse & !noadjCov`, `--condition`).
- `--is_single_in_groupTest`: `SAIGE_SPATest_Region.R:917` overwrites it with TRUE on every `--r.corr=0` run, so FALSE cannot suppress the file — a BURDEN-only knob. The raw value does still gate one stdout line (`SAIGE_Test_main.R:778`).
- `--weights.beta` / `--weights_for_condition`: region-only for results despite the help strings. A **single** beta pair is unused if the group file has weight lines; extra pairs on a marker run are discarded without a word; `weights.beta` is inert on FlexRV transform sets. `weights_for_condition` bites only on the **binary** region path (the quantitative region terms cancel algebraically), and an explicit all-zeros vector is indistinguishable from absent.
- `--cateVarRatioMinMACVecExclude` / `MaxMACVecInclude` (both steps): dead unless the model has >1 ratio per type (step 2) / `--isCateVarianceRatio=TRUE` (step 1) — still type-checked on the path that ignores them.
- all five `flexRV_*`, both weight booleans, `--is_output_weight_column`, the region cell/chunk options: a marker-path run swallows them, exit 0, no warning.
- `--sampleFile_male` and `--X_PARregion`: never opened, never validated unless `--is_rewrite_XnonPAR_forMales=TRUE`.
- step 1 `--usePCGwithSparseGRM`: its setter is inside the `useSparseGRMtoFitNULL` branch, so on a full-GRM fit the value never reaches C++.

**Live, but cannot move a result byte:**
- `--markers_per_chunk` (finding E0 — the only C++ consumer is a stdout counter; changes the `.index` and hence restart granularity, and is ignored entirely on group runs).
- `--groups_per_chunk` and `--markers_per_chunk_in_groupTest` (RAM, chunk count, threading block, `.index` line count, the printed SKAT-O dedupe measurement). Note: the byte-inertness of `markers_per_chunk_in_groupTest` is asserted only on a BURDEN run; on SKAT-O it changes GEMM sub-block shapes and **no test pins that** — treat it as untested, not proved.
- `--firth_solver=qr` vs `closed`: ≤8e-13 relative, byte-identical tables on 70,709 rows at N=50,000; genuinely live for which code runs and for −27% wall time. A third value is refused before a marker is read.
- `--minGroupMAC_in_BurdenTest=1` (the default) is inert **by construction** — a surviving cell has at least one allele.

**Silent no-ops needing a companion option:** `--sparseGRMSampleIDFile` and step 2's `--relatednessCutoff` without `--sparseGRMFile` (parsed, passed two levels down, dropped, no warning — both are *live* when the GRM is given, so this is a missing-companion fact, not an inert option). Step 1's `--relatednessCutoff=0` is a literal no-op (`drop0` tol 0) while the R API default is 0.125.

**Documentation defects found (each a false sentence a user reads):**
1. `--pval_cutoff_for_fastTest`: help *and* roxygen say "no effect on single-variant runs" — `Main.cpp:1007` is the single-variant escalation gate, and `test-step2-allele-order.R:124-140` pins the partition. docdrift.py shape (copied twin).
2. `--idstoIncludeFile` / `--condition`: help says `chr:pos_ref/alt` in three copied places (wrapper, roxygen, `man/`); on this branch both match bim col 2 / pvar ID. The `markerInfo$ID2` fallback is dead — nothing ever assigns ID2.
3. `--isCateVarianceRatio`: help claims "if TRUE, then LOCO=FALSE" — nothing enforces it; there is no LOCO interaction anywhere.
4. `--is_output_moreDetails`: help names columns (`homN_Allele2_cases`, …) the writer never emits; the real names are `N_case_hom`/`N_case_het`/`N_ctrl_hom`/`N_ctrl_het`. Also inherited upstream.
5. `--X_PARregion`: help says the format is `start:end`; the parser demands `start-end`.
6. `--IsOverwriteVarianceRatioFile`'s refusal message opens "WARNING", says the ratios "will be output to" the file, then stops the run.
7. Manifest exclusivity: help and roxygen on both flags claim mutual exclusivity with the manifest, but `MultiModel.R:101` enforces only `--GMMATmodelFile`; `--varianceRatioFile` is accepted and silently ignored (the manifest's column 3 is mandatory anyway, so nothing is lost — but the message is missing). Ledger-worthy: extend the condition.

**Dead code found while tracing (not options):** `SAIGE_Test_main.R:507-509` zeroes `max_MAC_use_ER` for a non-binary model *after* the C++ global is set, and the variable is never read again (inherited); `g_cfg.isOutputMoreDetails` is written and never read (upstream's `g_isOutputMoreDetails` too — so `setMarker_GlobalVarsInCPP`'s first argument could be deleted without moving a golden); step 1's `isCovariateTransform` and `tolCoef` are dead formals in both fitters; the four arguments into `getSparseSigma` are read only by its unreachable GRM-building arm.

**No step-0 or step-1 option is refused outright.** All refusals there are *combinations*: LOCO=TRUE + useSparseGRMtoFitNULL=TRUE, FemaleOnly + MaleOnly, PLINK + PGEN input (either side, any one path), an existing VR file without `--IsOverwriteVarianceRatioFile`, no genotype input at all (step 0). Step 2's combination refusals: manifest + `--GMMATmodelFile`, manifest + `--subSampleFile`, manifest + `--is_rewrite_XnonPAR_forMales`, manifest with models on different sample sets on a region run, disagreeing `.index` files on a multi-model restart, `r.corr` ∉ {0,1}, `--flexRV_score` with a differing `--maxMAF_in_groupTest` / any `--maxMAC_in_groupTest` / >1 annotation mask / `--condition` without `--weights_for_condition`, duplicate weight-set names, negative weights under the multi-weight schema, a monomorphic conditioning marker (#73).

---

## 5. THIS BRANCH'S OWN OPTIONS

**Step 0 (4):** `--plink2File`, `--pgenFile`, `--pvarFile`, `--psamFile`.
**Step 1 (6):** the same PGEN quartet, plus `--seed` (upstream hardcodes `set.seed(1)`, #8) and `--tolCoef` (#125).
**Step 2 (12):** `--pgenFile` / `--pvarFile` / `--psamFile` (upstream has only `--pgenPrefix`), `--GMMATmodel_varianceRatio_multiTraits_File`, `--is_output_weight_column`, `--firth_solver`, `--pval_cutoff_for_fastTest` (as a CLI flag; the R-API parameter and its 0.05 default are upstream's), and the five `--flexRV_*`.

**Upstream options whose capability is effectively new here:** `--phenoCol` / `--traitType` (comma lists = several fits per process, per-phenotype outputs + `.multiTraits.txt`); `--weights.beta` (semicolon pair list, one beta weight set per extra pair, #19); `--LOCO` (default NULL so the sparse-fit resolution is explicit and an explicit TRUE is refused, #112); `--is_no_weight_in_groupTest` (appends a `No_Weights` set where upstream stopped); `--rangestoIncludeFile` (non-functional upstream, #11); `--maxMAC_in_groupTest` (silently produced duplicate bins upstream, #61); `--nThreads` in step 2 (upstream forked `mclapply` and did nothing without `--idstoIncludeFile`, #30); `--is_rewrite_XnonPAR_forMales` (chrX-only, #69); `--minGroupMAC_in_BurdenTest` (default 5 → 1, both region paths, per-set MAC, #126/#118).

**Removed with the strip:** step 1's `--eventTimeCol`, `--eventTimeBinSize`, `--pcgforUhatforSurvAnalysis`; all bgen/vcf/sav input options and step 3 everywhere.

---

## 6. TRAPS

1. `--X_PARregion`'s CLI default is `""`, not the hg19 windows the help advertises — so a plain CLI `--is_rewrite_XnonPAR_forMales` run doubles males across the **whole** X, PAR included, and silently resets `minInfo=0` and `minMAF=1/(2N)` for the entire run.
2. `--minGroupMAC_in_BurdenTest` is either inert (1, by construction) or dangerous: any value >1 can abort a SKAT-O run in `regionFinish`'s `rbind` when a region has ≥2 refused cells (#128, open; died at region 4 of 1,614 on sim50k at gate 5).
3. `--maxMAC_in_groupTest` converts to MAF with the **model's** in-sample N — 1,000 on the bundled fixture, not the fam's 10,000, so the implied MAF is 10× a fam-based reading.
4. `--dosage_zerod_cutoff=0.2` overrides `--impute_method=mean` on every marker with MAC ≤ 10, so `mean` and `minor` agree there unless the cutoff is 0 (#109).
5. `--relatednessCutoff` exists in all three steps with different defaults (0.125 / 0 / 0) and downstream can only **raise** it — a lower value cannot restore pairs step 0 discarded, and step 0's screen+filter pair means the `.mtx` is an intersection over two applications of one number.
6. `--idstoIncludeFile` and `--rangestoIncludeFile` **union**, they do not intersect (both append to one vector, one `%in%` at the end), and a file matching nothing leaves `markerInfo` at 0 rows with no stop in the reader — the failure surfaces later, in `SAIGE_SPATest_Marker.R:99` or the #113 group-file guard.

Everything above rests on lines read in the tree at `39b97ea2` (clean) or via `git show main:`; read-only throughout, nothing built or run. Two claims are code-read only and marked as such in the working notes: the behaviour after an unrecognised `--AlleleOrder` value, and step-0 `.mtx` byte determinism at `--nThreads > 1` (pair order comes out of a `tbb::concurrent_vector`, so the printed log preview will differ even if the bytes do not).
---

# Appendix: the annotation group-file format

Read off the parser (`checkGroupFile` and `SAIGE.getRegionList_new`), not the
docs. Whitespace- or tab-separated; no header. Every line is
`<REGION> <TYPE> <value1> <value2> ...` and **every value list must have exactly
as many entries as the region's `var` line**.

## Hard constraints, all enforced with a `stop()`

- **Region 1 defines the layout and every region must match it**:
  `nline_per_gene = 2 + nScoreSets + nWeightSets` is computed once. A region
  with more, fewer or differently-named weight lines is an error naming the
  region.
- **Order is fixed**: `var`, then `anno`, then the optional `score` line, then
  the weight lines. A score line after a weight line, or a second score line,
  is rejected.
- Marker IDs match **bim column 2** (PLINK) or the pvar `ID` column (PGEN) --
  NOT `chr:pos_ref/alt`, whatever three copies of the help text say.
- `Cauchy` is a reserved weight-set name; duplicate set names are rejected.
- Score values must be numeric and within `[0, 1]` (the cbrt/sbeta transforms
  are undefined outside it).
- Weights must be numeric; negative weights are rejected **only** when there is
  more than one weight set, because that is when they reach a Cauchy
  combination that cannot take them.

## 1. Without weights -- 2 lines per region

```
GENE1 var  rs1  rs2  rs3
GENE1 anno lof  lof  missense
GENE2 var  rs51 rs52
GENE2 anno missense lof
```
Weights come from `--weights.beta` (default `1,25`), i.e. `dbeta(MAF, 1, 25)`.
Fixture: `extdata/input/group_new_snpid.txt`.

## 2. With weights -- 2 + K lines per region

The weight line's type field is `weight`, `weight:NAME` or `weight,NAME` (both
delimiters work). A bare `weight` gets the set name `WEIGHT`.

**K = 1** (`group_new_snpid_withWeights.txt`) -- byte-identical to legacy output:
```
GENE1 var        rs1  rs2  rs3
GENE1 anno       lof  lof  missense
GENE1 weight     7.8  3.9  8.9
```

**K = 2** (`group_new_snpid_multiWeights.txt`) -- adds the multi-weight schema
(a `Weight` label column, `Pvalue_ACATV`, `Pvalue_ACATO`, a pooled row):
```
GENE1 var          rs1  rs2  rs3
GENE1 anno         lof  lof  missense
GENE1 weight:CADD  7.8  3.9  8.9
GENE1 weight,REVEL 0.01 0.38 0.75
```

Set order is normative: **custom** sets in file order, then **beta** sets from
`--weights.beta="a1,b1;a2,b2"` (composed only when more than one pair), then
**No_Weights** if `--is_no_weight_in_groupTest=TRUE`. A single beta pair is
unused when the file has weight lines. A weight of exactly 0 drops that marker
from that set (#117).

## 3. With FlexRV -- 2 + 1 lines per region

One `score` line per region, immediately after `anno`, carrying a
per-marker pathogenicity score in `[0, 1]`. No weight lines: FlexRV expands the
score into 16 x 12 = 192 weight sets inside C++.

```
GENE1 var         rs1 rs2 rs3
GENE1 anno        lof lof missense
GENE1 score:REVEL 1   1   0.8
```
Fixture: `group_new_snpid_flexrv.txt`. The run needs
`--flexRV_score=REVEL` (the name must match the line), `--flexRV_maxMAF=M`,
**one** `--annotation_in_groupTest` group and `--r.corr=1`; `--weights.beta` is
inert on transform sets, and `--maxMAF_in_groupTest` is replaced when absent and
refused when it disagrees with M. Scores for LoF markers are conventionally 1,
which is what `--flexRV_lofAnno` keys the `lof` transform on.

A `score` line and `weight` lines can coexist in principle (the layout
arithmetic allows `2 + 1 + K`), but a FlexRV run is driven by the score line
alone.

## 4. What the annotation mask does to markers it excludes

`--annotation_in_groupTest` is a MEMBERSHIP filter applied before any
weighting, not a re-weighting. The indicator is built as
`ANNO %in% <the group's labels>` (`R/SAIGE_SPATest_Region.R:2180`) and every
accumulation is gated on it (`src/Main.cpp:2948` and `:3018`). So with
`lof:missense` requested, a `synonymous`-labelled marker contributes to the
burden column, the ultra-rare collapse, `MAC`/`MAC_case`/`MAC_control`,
`Number_rare`/`Number_ultra_rare` and `maxMAFperAnno` -- none of them. Its
weight and its FlexRV score are parsed and validated, then unused. Label
everything once in the file and let each run choose its mask; that is the
intended use.

Three consequences that are not obvious:

- **An excluded marker is still read, QC'd and SINGLE-tested.** `doSingleTests`
  is gated on the run type only, never on the annotation
  (`src/Main.cpp:2699`), so on a SKAT-O run -- or `--r.corr=1` with
  `--is_single_in_groupTest=TRUE` -- it appears in `.singleAssoc.txt` with its
  own p-value while contributing to no group row. The single-variant file and
  the region table therefore cover DIFFERENT marker sets, and a row in the
  former is not evidence the marker was in any cell.
- **A requested group that matches nothing is dropped with a `warning()`, not
  an error** ("No markers are found for at least one annotation, so region N is
  skipped"). If only some requested groups match, the indicator matrix is
  rebuilt with just those, so the surviving cells RENUMBER -- a region's cell
  list is data-dependent on an ordinary multi-group run. A FlexRV run cannot
  hit this, having one group by construction.
- **Setting LoF markers' scores to 1 is a CONVENTION, not a requirement**, and
  nothing checks it (only `[0, 1]` is enforced). The `lof` score transform
  ignores the score entirely -- `src/FlexRV.cpp:13` is
  `return t_isLof ? 1.0 : 0.0`, keyed on the annotation label -- but the other
  15 score transforms read the score for every marker, LoF included. Scoring a
  LoF marker 0.3 therefore DOWN-weights it against a missense marker scoring
  0.9 under `cbrt` and `sbeta`, silently, which is backwards from the method's
  intent. Similarly `--flexRV_lofAnno` is checked against the labels present
  anywhere in the FILE, not per region: a region with no marker carrying a
  named LoF label gets an all-zero LoF mask, and on a binary trait a weight of
  exactly 0 in a SKAT-O cell returns NA (LEDGER #117), which is part of why
  FlexRV on a binary trait is an `--r.corr=1` method.
