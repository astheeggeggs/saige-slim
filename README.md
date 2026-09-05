# SAIGE (`strip` branch)

A fork of [SAIGE / SAIGE-GENE+](https://github.com/saigegit/SAIGE) reduced to the
three steps that matter for an association scan, with **PLINK bed/bim/fam and
PGEN pgen/pvar/psam as the genotype input for every step**:

| Step | Script (`extdata/`) | R function | What it does |
|---|---|---|---|
| 0 | `createSparseGRM.R` | `createSparseGRM()` | Sparse genetic relationship matrix |
| 1 | `step1_fitNULLGLMM.R` | `fitNULLGLMM()` | Null GLMM (AI-REML + PCG) and variance ratios |
| 2 | `step2_SPAtests.R` | `SPAGMMATtest()` | Single-variant tests, or region/gene tests |

The method is unchanged: a generalized mixed model that accounts for sample
relatedness, a saddlepoint approximation that calibrates case-control imbalance,
Firth bias-reduced effect sizes, and BURDEN / SKAT / SKAT-O for rare-variant
set-based tests. Binary and quantitative traits are supported.

The upstream documentation at <https://saigegit.github.io/SAIGE-doc/> is good
background on the **method**, but it documents a superset of this branch: readers,
a step 3 and trait types that are not here, and none of the features below. This
file is the reference for what this branch actually does.

## Differences from upstream

| Removed | Added |
|---|---|
| BGEN, VCF/BCF and SAV readers | PGEN input for step 0, step 1 and step 2 |
| Step 3 (LD matrices) | Several phenotypes fitted in one step-1 process |
| RareEffect | Several null models tested in one step-2 genotype pass |
| Survival / Cox traits | K weight sets per region test, in one pass |
| | A working `--nThreads` for step 2 |
| | A golden-output test suite and a divergence ledger |

Everything else — the model, the CLI option names, the output columns — is
preserved bug-for-bug against upstream. Divergences that were found rather than
introduced are recorded in `LEDGER.md` rather than fixed in passing; see
[Known limitations](#known-limitations).

## Install

```sh
pixi install
bash harness/build.sh
```

`harness/build.sh` cleans `src/*.o`, regenerates `RcppExports`, and runs
`R CMD INSTALL`. It also **fetches and compiles plink-ng (pgenlib) automatically**
into `.pixi/envs/default/lib/libplink2_includes.a`; plink-ng is not vendored.

`src/Makevars` hardcodes `-I../.pixi/envs/default/include`, so the in-tree
`.pixi/` environment is required. Works natively on macOS and Linux.

`make` at the repo root builds the Docker image instead (`docker/Dockerfile`:
Ubuntu 20.04 + pixi, which installs the three entry-point scripts into
`/usr/local/bin` and smoke-tests each with `--help`).

`plink2` the CLI is a pixi dependency only for generating PGEN test fixtures
(`linux-64` and `osx-64` only — bioconda ships no `osx-arm64` build). The
installed package never calls it.

## Quick start

All commands below are run from `extdata/`, against the bundled 1,000-sample
fixtures, with `Rscript` from the pixi environment:

```sh
source harness/env.sh          # see "Reproducibility" — do this before any run
cd extdata
R2() { pixi run --manifest-path ../pixi.toml Rscript "$@"; }
```

**Step 0** — build a sparse GRM:

```sh
R2 createSparseGRM.R \
  --plinkFile=input/plinkforGRM_1000samples_10kMarkers \
  --nThreads=1 \
  --outputPrefix=./qs/sparseGRM \
  --numRandomMarkerforSparseKin=5000 \
  --relatednessCutoff=0.05
```

**Step 1** — fit the null model for a binary trait:

```sh
R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=./qs/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=./qs/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --isCateVarianceRatio=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix=./qs/binary
```

**Step 2** — single-variant tests:

```sh
R2 step2_SPAtests.R \
  --bedFile=./input/genotype_100markers.bed \
  --bimFile=./input/genotype_100markers.bim \
  --famFile=./input/genotype_100markers.fam \
  --AlleleOrder=alt-first \
  --GMMATmodelFile=./qs/binary.rda \
  --varianceRatioFile=./qs/binary.varianceRatio.txt \
  --SAIGEOutputFile=./qs/assoc.txt \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  --is_output_moreDetails=TRUE
```

The same three steps read PGEN by substituting the input flags — `--plink2File`
for steps 0 and 1, explicit `--pgenFile/--pvarFile/--psamFile` for step 2:

```sh
R2 createSparseGRM.R --plink2File=input/plinkforGRM_1000samples_10kMarkers ...
R2 step1_fitNULLGLMM.R --plink2File=./input/nfam_100_..._22chr ...
R2 step2_SPAtests.R --pgenFile=./input/genotype_100markers.pgen \
                    --pvarFile=./input/genotype_100markers.pvar \
                    --psamFile=./input/genotype_100markers.psam ...
```

Step-0 and step-1 outputs are **byte-identical** between the two readers on
mirrored fixtures. That is structural rather than coincidental: the PGEN front end
decodes each variant and re-encodes it as a PLINK-1 record, then hands it to the
same QC / imputation / repack body.

`Rscript <script> --help` prints every option with its default.

---

## Step 0 — `createSparseGRM.R`

Genotype input, one of:

- `--plinkFile=<prefix>`, or explicit `--bedFile` / `--bimFile` / `--famFile`
- `--plink2File=<prefix>`, or explicit `--pgenFile` / `--pvarFile` / `--psamFile`

PLINK and PGEN input are mutually exclusive.

Options that matter: `--numRandomMarkerforSparseKin` (default 2000),
`--relatednessCutoff` (default 0.125), `--minMAFforGRM`, `--maxMissingRateforGRM`,
`--isDiagofKinSetAsOne`, `--memoryChunk`. Note `--nThreads` **defaults to 16**
here, unlike steps 1 and 2 where it defaults to 1.

Two outputs, named from the prefix and the parameters:

```
<outputPrefix>_relatednessCutoff_<cutoff>_<nMarkers>_randomMarkersUsed.sparseGRM.mtx
<outputPrefix>_relatednessCutoff_<cutoff>_<nMarkers>_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
```

**Reproducibility:** the markers are drawn with `sample()`, and the `set.seed(1)`
lives in `extdata/createSparseGRM.R`, not in `createSparseGRM()`. Calling the R
function directly is therefore not reproducible across processes; go through the
script, or set the seed yourself.

Step 0 and step 1 share one global genotype store, so they cannot run in the same
R process.

## Step 1 — `step1_fitNULLGLMM.R`

Genotype input is as for step 0. Phenotypes come from `--phenoFile` with
`--phenoCol`, `--sampleIDColinphenoFile`, `--covarColList` and `--qCovarColList`
(categorical covariates must appear in **both** `--covarColList` and
`--qCovarColList`). `--traitType` is `binary` or `quantitative`.

**GRM choice.** With `--useSparseGRMtoFitNULL=TRUE` plus `--sparseGRMFile` /
`--sparseGRMSampleIDFile`, the fit uses the sparse GRM from step 0. Without it the
full (dense) GRM is built from the step-1 PLINK/PGEN file. `--usePCGwithSparseGRM`
and `--useSparseGRMforVarRatio` control the two places the sparse matrix can be
used independently of the fit itself.

**Variance ratios.** `--isCateVarianceRatio=TRUE` estimates one ratio per MAC bin;
the bins are `--cateVarRatioMinMACVecExclude` / `--cateVarRatioMaxMACVecInclude`.
`--skipVarianceRatioEstimation` skips them entirely.

**LOCO.** `--LOCO=TRUE` (the default) fits leave-one-chromosome-out;
`--isLowMemLOCO` writes one model file per chromosome instead of one large object.

### Outputs, and the contract with step 2

The interface between the two steps is implicit and name-based. It is worth
knowing exactly, because a mismatch surfaces as a confusing step-2 error:

- `<outputPrefix>.rda` must contain a variable literally named `modglmm`. Step 2
  `load()`s it blind — there is no schema or version field.
- With `--LOCO=FALSE` the model goes to `<outputPrefix>_noLOCO.rda`; with
  `--isLowMemLOCO` to `<outputPrefix>_chr<j>.rda`. With LOCO on, step 2's
  `--chrom` is mandatory.
- `<outputPrefix>.varianceRatio.txt` is headerless, three whitespace-separated
  columns `value type category`, with `type` in `{sparse, null, null_noXadj}` and
  one row per category:

  ```
  0.948661317135906 null 1
  0.939556503863059 null_noXadj 1
  0.942713872154952 null 2
  0.932272959819397 null_noXadj 2
  ```

**Which variance-ratio rows you get determines what step 2 will accept**, and this
is the most common source of confusion:

- A sparse-GRM fit emits `null` and `null_noXadj` rows only. Such a model must
  **not** be given `--sparseGRMFile` in step 2.
- `sparse` rows are emitted only by `--useSparseGRMforVarRatio=TRUE`. Region/group
  tests, and the step-2 sparse-GRM marker branch, need a model whose variance-ratio
  file has them:

  ```
  0.999983810709269 sparse 1
  0.948661317135906 null 1
  0.939556503863059 null_noXadj 1
  ...
  ```

### Several phenotypes in one process

`--phenoCol` accepts a comma-separated list. `--traitType` and `--invNormalize`
then take either one value per phenotype or one value for all of them; every other
option is shared.

```sh
R2 step1_fitNULLGLMM.R ... \
  --phenoCol=y_binary,y_quantitative \
  --traitType=binary,quantitative \
  --invNormalize=FALSE,TRUE \
  --outputPrefix=./out
```

Outputs are `./out_y_binary.rda`, `./out_y_binary.varianceRatio.txt`,
`./out_y_quantitative.*`, plus a manifest `./out.multiTraits.txt` — three
whitespace-separated columns, no header, one row per phenotype:

```
y_binary        ./out_y_binary.rda        ./out_y_binary.varianceRatio.txt
y_quantitative  ./out_y_quantitative.rda  ./out_y_quantitative.varianceRatio.txt
```

Step 2 consumes that file directly as
`--GMMATmodel_varianceRatio_multiTraits_File`. A phenotype whose variance-ratio
estimation was skipped gets a literal `NA` in column 3, and step 2 refuses that row
with a message saying why.

A **single** phenotype keeps the legacy naming (`./out.rda`) and writes no manifest,
so nothing about existing pipelines changes.

What this buys, stated honestly: the fits stay sequential and only the genotype
ingest is shared. On the sparse-GRM path that is roughly 1.4x at five phenotypes
(N = 1,000); **on the dense-GRM path it is about 1%**, because the fit is hundreds
of GRM passes against one ingest. Each phenotype's outputs are byte-identical to
what a standalone run would have produced, which is the gate the feature is held to.

### Quantitative traits and scale

A quantitative null model is fitted on the phenotype **standardised** by the
standard deviation of its covariate-adjusted values, then scaled back before being
saved: `theta` by `s^2`, linear terms by `s`, and the variance-ratio file by
`1/s^2`. Everything that reads the model or the ratios therefore sees phenotype
units, and `--tauInit` is given in phenotype units. A trait measured in any unit
fits the same way, so `--invNormalize` is not needed merely to fix scale.

## Step 2 — `step2_SPAtests.R`

Genotype input is `--bedFile`/`--bimFile`/`--famFile` **or**
`--pgenFile`/`--pvarFile`/`--psamFile` (mutually exclusive; step 2 has no
`--plink2File` prefix shorthand). `--AlleleOrder` is `alt-first` (default) or
`ref-first`. Restrict the scan with `--chrom`, `--idstoIncludeFile`,
`--rangestoIncludeFile`, `--minMAF`, `--minMAC`, `--maxMissing`.

The null model comes from `--GMMATmodelFile` + `--varianceRatioFile`, or from
`--GMMATmodel_varianceRatio_multiTraits_File` (below). Supplying `--groupFile`
switches the run from single-variant tests to region/gene tests.

### Single-variant tests

Output is one row per marker. With `--is_output_moreDetails=TRUE` on a binary
trait:

```
CHR  POS  MarkerID  Allele1  Allele2  AC_Allele2  AF_Allele2  MissingRate  BETA  SE
Tstat  var  p.value  p.value.NA  Is.SPA  AF_case  AF_ctrl  N_case  N_ctrl
N_case_hom  N_case_het  N_ctrl_hom  N_ctrl_het
```

(one tab-separated header line; the case/control columns are what
`--is_output_moreDetails` adds).

Each marker takes one arm of a cascade: score test, then the saddlepoint
approximation if `|Z|` exceeds `--SPAcutoff`, then an exact test if the trait is
binary and MAC ≤ `--max_MAC_for_ER`, then Firth if `--is_Firth_beta=TRUE` and
p < `--pCutoffforFirth`. The variance ratio is selected per marker from its MAC.

`Is.SPA` is worth reading carefully: it means **the SE was back-derived from the
p-value**, which happens for the saddlepoint approximation *and* for the exact
test — not "SPA was used".

**Conditional analysis:** `--condition` takes a comma-separated list of markers to
condition on, given as **rsIDs** for PLINK input (PLINK group and condition
handling keys on bim column 2), with optional `--weights_for_condition`.

```sh
R2 step2_SPAtests.R ... --condition=rs13,rs79
```

### Region / gene tests

`--groupFile` gives each region two or more whitespace-separated lines, keyed on
the same marker IDs the genotype file uses:

```
GENE1 var  rs1 rs2 rs3 ...
GENE1 anno lof lof missense ...
GENE1 weight:CADD  7.8 3.9 8.9 ...       # optional, zero or more
```

`--annotation_in_groupTest` selects which annotation combinations to test, with
`:` separating annotations *within* one test and `,` separating tests —
`lof,missense:lof,missense:lof:synonymous` runs three. `--maxMAF_in_groupTest`
takes a comma-separated list of MAF ceilings. `--r.corr` picks the test
(`0` = SKAT, `1` = BURDEN, omitted = SKAT-O). Markers below
`--MACCutoff_to_CollapseUltraRare` are collapsed into a single pseudo-marker.
`--is_single_in_groupTest=TRUE` additionally writes per-marker results to
`<SAIGEOutputFile>.singleAssoc.txt`.

The legacy output schema is:

```
Region  Group  max_MAF  Pvalue  Pvalue_Burden  Pvalue_SKAT  BETA_Burden  SE_Burden
MAC  MAC_case  MAC_control  Number_rare  Number_ultra_rare
```

### Several weight sets in one pass

A region test can run K ≥ 1 weight sets in a single genotype pass. Sets come from
three places, and this order is normative:

1. **Custom** sets, from group-file `weight` lines. `weight:NAME` or `weight,NAME`
   names the set; a bare `weight` gives a set named `WEIGHT`.
2. **Beta** sets, from `--weights.beta` — semicolon-separated `a,b` pairs, comma
   within a pair. Quote the value because of the semicolon:
   `--weights.beta="1,25;1,1"`. Each becomes a `beta(MAF; a, b)` set named
   `Beta_a_b`. Composed **only when more than one pair is given**; with a single
   pair the behaviour is exactly as before (group-file weights take precedence, and
   the pair still supplies single-variant weighting and the conditioning-weight
   default).
3. **`No_Weights`**, from `--is_no_weight_in_groupTest=TRUE`. With no other sets in
   play the all-ones weights simply replace the defaults, as they always did; with
   other sets present it is appended as a final set.

Every region must carry the same weight lines, with the same names in the same
order. A mismatched name, order or count is a hard error, not a silent relabel;
duplicate names and the reserved name `Cauchy` are rejected.

Results are one row per (annotation × max MAF × weight set), **weight innermost**,
with a `Weight` label column immediately after `max_MAF` and, for SKAT-O, the
additional `Pvalue_ACATV` and `Pvalue_ACATO` columns:

```
Region  Group  max_MAF  Weight  Pvalue  Pvalue_Burden  Pvalue_SKAT  Pvalue_ACATV
Pvalue_ACATO  BETA_Burden  SE_Burden  MAC  MAC_case  MAC_control  Number_rare
Number_ultra_rare
```

There is still exactly one pooled `Cauchy` row per region; its `Pvalue_ACATO` is
the ACAT over the per-row ACAT-Os.

**The schema changes if and only if** more than one weight set is in play, or
`--is_output_weight_column=TRUE` is given. K = 1 without that flag is byte-identical
legacy output. `--is_output_weight_column` exists so a pipeline can pin the wider
schema across runs with different weight configurations. Appending to, or
restarting into, a file written under the other schema is refused with a message
naming the mismatch.

One deliberate limitation: conditioning-marker weights are a **single** vector
across all weight sets (`--weights_for_condition`, else `beta(1, 25)`). Per-set
conditioning weights are not implemented.

### Several null models in one genotype pass

```sh
R2 step2_SPAtests.R ... --GMMATmodel_varianceRatio_multiTraits_File=./out.multiTraits.txt
```

The manifest is the three-column file step 1 writes (or one written by hand: name,
`.rda`, variance-ratio file). It is mutually exclusive with `--GMMATmodelFile` /
`--varianceRatioFile`.

Results go to `<SAIGEOutputFile>_<name>` per model, each with its own `.index`. A
**one-row manifest writes the unsuffixed `<SAIGEOutputFile>`** and is identical, in
bytes and in name, to the equivalent single-model run.

Models whose step-1 sample IDs are `identical()` share one genotype reader and one
decode per marker; that is where the saving comes from. Sample sets are **grouped,
never unioned**, so per-marker `AC_Allele2` / `AF_Allele2` / `MissingRate` / N are
always computed over that model's own samples.

Restrictions, all refused in R before any work starts:

- single-variant tests only — no `--groupFile`
- no `--is_rewrite_XnonPAR_forMales`
- no `--subSampleFile`

On restart, each model's `.index` is truncated to its own recorded size. If the
index files disagree about which chunk to resume from, the run stops and names
them rather than guessing.

### `--nThreads`

One mechanism, in process. For **single-variant** tests the C++ marker loop is
threaded with `RcppParallel::parallelFor`: one process, one copy of the null model,
one output file, and output **byte-identical to a one-thread run** at 1, 2, 4 and 8
threads. Restarting an interrupted threaded run behaves exactly as a serial one.

For **region tests** the C++ loop stays on one thread and `--nThreads` is spent on
the SKAT-O post-processing instead, which profiling puts at roughly half of a group
run at N = 50,000 — measured 1.30x at 4 workers, byte-identical, with no further
gain at 8. Unix only; elsewhere a group run stays serial.

A genotype reader that cannot be cloned falls back to one thread with a printed
reason and an R `warning()` — deliberately a degradation rather than an error, so
one `--nThreads` can be scripted across mixed inputs. The only remaining such case
is PLINK when the `.bed` could not be memory mapped; **PGEN threads**.

Note: the `--nThreads` help text in `extdata/step2_SPAtests.R` is stale on both
points. The behaviour described here is what the code does.

### Restart and append

Step 2 writes a `<SAIGEOutputFile>.index` alongside the results, recording the
output file size at each chunk boundary. An interrupted run restarted with the same
arguments truncates to the last recorded size and resumes at the right chunk.
`--is_overwrite_output=TRUE` starts over instead.

## Reproducibility

**`source harness/env.sh` before any run.** It pins
`RCPP_PARALLEL_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.
Without it, full-GRM step-1 fits are not deterministic: float
`parallel_reduce` noise in the TBB kernels feeds back into `.Random.seed` through
the trace estimator and changes which markers get sampled for the variance ratios.

The honest limit: that buys reproducibility on the bundled N = 1,000 fixtures, and
not everywhere. At N = 50,000 with `--useSparseGRMforVarRatio=TRUE`, eleven
identical runs — with `env.sh` sourced and threads pinned to 1 — produced six
distinct variance-ratio files (`LEDGER.md` #56).

## Known limitations

These are documented decisions, not oversights. Each has a `LEDGER.md` row.

- **Dosages are hard-called for steps 0 and 1.** The GRM's genotype store is 2 bits
  per sample, so a PGEN dosage is rounded to the nearest of {0, 1, 2} and one
  warning naming the count is raised. For imputed data this discards information
  the GRM and the variance ratios could have used. Step 2 reads dosages normally.
- **PGEN/PLINK byte-identity assumes `pvar` ALT matches `bim` column 5**, which is
  what `plink2 --bfile X --make-pgen` produces. A pgen built with REF/ALT swapped
  relative to a bim flips every marker; the GRM is invariant to that, so this costs
  byte-identity rather than correctness.
- **`--tolPCG` is an absolute residual tolerance, not a relative one** (#6). No PCG
  loop normalises by the right-hand side, so the same `--tolPCG` means something
  different at every N and phenotype scale. This matters at biobank N with a dense
  GRM; sparse-GRM fits do not enter a PCG loop at default flags.
- **`--is_fastTest` is inert for region tests** (#51). The region path never reads
  the flag; the R side escalates regardless, and the flag costs about +16% at
  N = 50,000 for byte-identical output. A warning says so.
- **A saturated SKAT-O cell degrades a region's Cauchy row to a Bonferroni** (#86).
  The `SKAT` package's kernel can return exactly 1, and both CCT implementations
  then fall back to `min(1, n * min(p))`.
- **`Is.SPA` means the SE was back-derived from the p-value** (#80), covering the
  exact test as well as the saddlepoint approximation.

A SKAT kernel failure is **counted, not hidden**: the affected rows stay `NA`, the
run prints `SKAT failures: N ...` before `Analysis done!`, and
`<SAIGEOutputFile>.skatFailures.txt` names each failed cell with the kernel's
message.

`LEDGER.md` is the full divergence ledger. Do not read it whole — it is ~90 long
rows. Use the reader:

```sh
python3 harness/tools/ledger.py open      # rows still open
python3 harness/tools/ledger.py 56        # one row, wrapped
python3 harness/tools/ledger.py grep spa  # search
```

## Building and testing

```sh
bash harness/run_tests.sh [fast|standard|full]
```

- **fast** (~20 s) — unit tests over the R seams, kernels pinned to closed forms,
  invariant checks over the goldens, and the suite's own comparator self-tests.
- **standard** (~2 min, the default) — fast, plus every step-2 CLI configuration run
  against `harness/goldens/step2/`, the step-1 genotype-store and PGEN checks, and
  the characterization tests.
- **full** (~20 min, same-platform only) — standard, plus golden regeneration, a
  step-1 determinism double-run, and a permuted-phenotype negative control.

Goldens live in `harness/goldens/`. Step-0 and step-2 outputs are compared
**byte-exact**; step-1 outputs within tolerance (1e-2 for theta and variance
ratios, 1e-3 otherwise). Byte comparison is only valid on the platform that
generated the goldens — elsewhere, and in CI, set `SAIGE_COMPARE_MODE=numeric` to
compare parsed tables at 1e-5 relative.

`bash harness/prepush.sh` is the pre-push gate: the ledger lint plus the standard
tier in **both** compare modes. `git config core.hooksPath .githooks` arms it.

Contributors should read `CLAUDE.md`, then the newest `docs/state/handoff-*.md`,
then `docs/state/next.md`.

## Licence and citation

GPL-3 (`LICENSE`). Upstream SAIGE / SAIGE-GENE+ is by Wei Zhou, Bram Gorissen,
Zhangchen Zhao, Wenjian Bi, Seunggeun Lee and Cristen Willer; cite the SAIGE and
SAIGE-GENE+ papers as you would for upstream. Upstream repository:
<https://github.com/saigegit/SAIGE>; upstream documentation, which covers a
superset of this branch: <https://saigegit.github.io/SAIGE-doc/>.
