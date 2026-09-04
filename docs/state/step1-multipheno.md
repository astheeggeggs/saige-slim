# Multi-phenotype step 1 (Phase 2, commits N1-N4)

**What it is.** `--phenoCol=y1,y2,...` fits several null models in ONE process,
sharing the genotype ingest. The fits stay sequential; what is shared is the
`setgeno` read and bit-pack, which is O(M x N) and depends on nothing about the
phenotype except which SAMPLES are in the fit.

**The claim.** Each phenotype's `.rda` and `.varianceRatio.txt` are BYTE-
IDENTICAL to a standalone run's. Same argument as `--nThreads`: one code path
re-run, over a store whose contents do not depend on which phenotype asked for
it.

## The CLI

```
--phenoCol=y_binary,y_quantitative \
--traitType=binary,quantitative \
--invNormalize=FALSE,TRUE \
--outputPrefix=out
```

writes `out_y_binary.rda`, `out_y_binary.varianceRatio.txt`,
`out_y_quantitative.*`, and

```
out.multiTraits.txt:
  y_binary         out_y_binary.rda         out_y_binary.varianceRatio.txt
  y_quantitative   out_y_quantitative.rda   out_y_quantitative.varianceRatio.txt
```

which step 2 takes verbatim as
`--GMMATmodel_varianceRatio_multiTraits_File`. **One phenotype writes the
unsuffixed `out.rda` and no manifest** -- the legacy path, unchanged.

`--traitType` and `--invNormalize` are per phenotype (one value each, or one for
all). Everything else is shared. `--isLowMemLOCO` gives
`<prefix>_<pheno>_chr<j>.rda`; a phenotype whose variance-ratio estimation was
skipped gets a literal `NA` in column 3, which step 2 refuses with a message
saying why.

## Why the bytes match

The ONLY RNG-dependent input to an ingest is the `arma::randi(1000)`
variance-ratio marker pool, drawn from R's RNG inside `setgeno`
(`src/SAIGE_fitGLMM_fast.cpp`). The driver calls `set.seed(seed)` before every
fit, so:

* the first fit of a group draws exactly the pool a standalone run would;
* later fits in the group reuse that store, so they get the same pool by
  construction.

Nothing between the seed and the ingest consumes RNG (`fread`, `glm`,
`checkPerfectSep`, `Covariate_Transform`, `getsubGRM`, and `invNormalize`'s
`rank(ties="average")` are deterministic); `GetTrace` re-seeds 200 on entry;
`extractVarianceRatio` re-seeds `seed+1` before its own `sample()`s.

**Drill D2 is the proof that the per-fit seed is load-bearing**: removing it
makes the two-GROUP case's `varianceRatio.txt` differ from standalone while
every single-group test stays green.

## Where it lives

* `R/SAIGE_fitGLMM_multi.R` -- `.step1Geno` (the loaded store),
  `step1IngestKey`, `ensureGenoLoaded`, `releaseGeno`, `recyclePerPheno`,
  `step1PhenoGroups`, `fitNULLGLMM_multi`.
* `R/SAIGE_fitGLMM_fast.R` -- the four ingest call sites go through
  `ensureGenoLoaded`; `fitNULLGLMM` gained `releaseGenoAtEnd` and now returns
  the paths it ACTUALLY wrote (`modelOut` is reassigned under `--isLowMemLOCO`,
  `outputPrefix` is mutated in place by `--FemaleOnly/--MaleOnly`).
* `src/SAIGE_fitGLMM_fast.cpp` -- the LEDGER #55 reset, without which the second
  ingest in a process is garbage. That is N1, and it comes first for that
  reason.

The ingest KEY is the resolved files, the sample subset, `memoryChunk`,
`isDiagofKinSetAsOne` and the QC/variance-ratio GATES -- the gates because
`setminMAFforGRM` / `setmaxMissingRateforGRM` / `setminMAC_VarianceRatio` are
pushed BEFORE `setgeno` and decide which markers it keeps. The reuse test is
`identical()` on the whole key, so reuse can only ever be a saving and never a
difference.

## What this does and does not buy

Sharing the ingest saves the INGEST. At N=1,000 on the sparse-GRM golden the
ingest is `time: 2164 ms` of a ~6.5 s run, so P=5 predicts ~1.4x. The ingest is
O(M x N) while the sparse fit is O(nnz x iterations), so its share -- and the
saving -- grows with the number of GRM markers; on the sparse path the whole
ingest exists only to sample ~30 variance-ratio markers, so at biobank M the
saving approaches (P-1) ingests. **On the DENSE-GRM path it buys ~1%**: the fit
is 300-450 GRM passes against one for the ingest. The win there is batched
multi-RHS PCG across phenotypes, which is deliberately out of scope.

## Out of scope, recorded so it is not rediscovered

* Batched multi-RHS PCG across phenotypes (the dense-GRM win).
* Anything that would make a phenotype's output depend on which other
  phenotypes were in the run. The byte-identity gate is what forbids it.
