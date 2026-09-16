# All of Us chr20: slim vs upstream 1.4.2, first results (2026-09-16)

Executes the plan in `handoff-2026-09-16-aou-upstream-compare.md`, which registered
its predictions before any output existed. **Two of those predictions were wrong,
and both were scoping errors in the prediction rather than findings** -- recorded
below as such.

Run inside the All of Us Researcher Workbench. Cohort `afr`, chr20, All of Us v8
exome (`aou.exome_split.v8.qced.chr20.afr.af01.popmax01`), **in-sample N = 70,272**,
537 genes, 6 annotation masks x 3 max-MAF bins = 18 cells per gene, group file
`aou.v8.chr20.afr...spliceai=0.20_cadd=28.1_revel=0.773.canonical.txt.gz`.
Both arms read the SAME step-1 model and variance-ratio file (fitted by the
production 1.4.2 pipeline), so step 1 is not a confound. `FASTTEST=FALSE`
on every comparison arm. Images by digest: `saige:1.4.2` `sha256:678b13d6...`,
`saige-slim:1.5.2-dev` `sha256:eb8de37e...`.

## 0. Compatibility: slim step 2 is drop-in on a 1.4.2 model

The blocking prerequisite passed. Slim reads the 1.4.2 `.rda` and the 3-column
`varianceRatio.txt` unchanged, printing `variance Ratio null is 1.093018 1.089024`
-- two categorical rows, in use. This is stronger evidence than "it ran": slim
hard-stops on a file with no `sparse` row when `--sparseGRMFile` is given, on a
file with no `null` row (audit E1), on the pre-1.0.6 legacy format (#13), and on a
categorical count that disagrees with `cateVarRatioMinMACVecExclude`. Four guards
passed, not four guards skipped.

**So step 2 can be swapped without refitting step 1.** Do not set
`--is_noadjCov=TRUE` on these models: the file has no `null_noXadj` rows and the
build is `-DARMA_NO_DEBUG`, so an empty ratio vector is read off the end silently.

## 1. Speed and memory

Same inputs, same grid, one thread, one task each.

| phenotype | slim s | upstream s | slim/upstream | slim peak RSS | upstream peak RSS |
|---|---|---|---|---|---|
| BMI (quantitative) | 1263 | 2775 | **2.20x faster** | 3867 MB | 929 MB |
| phe_153 (binary) | 1148 | 2665 | **2.32x faster** | 3825 MB | 864 MB |
| phe_250.2 (binary) | 1002 | 2150 | **2.15x faster** | 3748 MB | 857 MB |

**Slim is ~2.2x faster and uses ~4.2x more memory.** The speed direction matches
`harness/sim/upstream_compare/RESULTS-2026-09-09.md` at sim50k; the memory is the
trade (slim keeps region chunks in RAM via `RegionChunkStore` where upstream
spills to chunk files). `[ASSUMPTION: RegionChunkStore is the mechanism -- the
direction and magnitude fit, but this was not profiled.]`

**Node sizing follows from the RSS, and it settles threads-vs-instances:** ~4 GB
per instance means a 4-vCPU/16 GB node fits four single-threaded tasks, one per
core. Independent tasks scale linearly; threading a SKAT-O run is bounded at
~1.5x. So run 1 core / 6 GB per task and pack instances. The production submission
had `--min-ram 16`, a 4x over-provision.

## 2. `--is_fastTest` on a group run: LEDGER #51 replicates on real data

`FASTTEST=TRUE` vs `FALSE`, same inputs:

- **25 escalations** vs 0 (the flag took effect -- the zero-side alone proved nothing)
- **results file and `.singleAssoc.txt` byte-identical**
- **1283 s vs 1096 s = +17.1%**

LEDGER #51 measured **+16.3%** at sim50k (N = 50,000, 210 genes, 15 cells). This is
an independent replication at N = 70,272, 537 genes, 18 cells, on real data.
**Drop `--is_fastTest` from group runs.** It remains a real two-stage screen on the
marker path, so keep it for `TESTTYPE=variant`.

## 3. The cost model, and the q correction that overturned an analysis

**A cost analysis built on the wrong kernel dimension was refuted by its own
prediction.** The first pass computed per-cell `q = Number_rare + Number_ultra_rare`,
concluded that the all-variants union mask carried 82.6% of `sum(q^3)`, and
predicted that dropping two masks would take the run to 300-400 s.

Measured: **912 s, a ratio of 0.83** -- not 0.27-0.36, and worse than
linear-in-cells (0.67). The cell-size story is refuted for this data.

**The q definition was wrong**: ultra-rare markers are collapsed into ONE
pseudo-marker column (`src/Main.cpp`, the `max()` accumulation into
`urAt(sl, jmr)`; `Number_ultra_rare` counts what went IN), so the kernel dimension
is `Number_rare + 1`. In an exome most variants are ultra-rare, so the correct q
is far smaller and the cubic term never dominates. **Any share computed as
`Number_rare + Number_ultra_rare` must be recomputed.**

Two cell-count points, fitted exactly (no residual, so unvalidated):

```
1096 s = A + 18B      A ~ 544 s   independent of cell count
 912 s = A + 12B      B ~ 30.7 s  per cell, summed over all genes
```

So **~50% of the run does not scale with cells at all** -- genotype read, decode,
region setup, single-variant tests, SAIGE's own startup -- and the cell-dependent
half is roughly uniform across cells rather than concentrated in the wide ones.
`STARTUP_AND_DELIVERY_SECONDS = 1` measures only WRAPPER overhead; SAIGE's own
startup is inside `SAIGE_SECONDS`, so the truncation series remains the only route
to the real intercept.

**The grid question is closed by architecture, not preference.**
`--maxMAF_in_groupTest` is a single global list applied to every annotation (the
grid is the full cross product, `jmr = j*q_maf*q_weight + m*q_weight + r`,
`R/SAIGE_SPATest_Region.R:256`). So "the all mask at 0.001 only, the others at all
three bins" cannot be expressed in one run; it needs two, each paying the ~544 s
fixed cost: 1004 + 575 = 1579 s against 1096 s for the single 18-cell run --
**44% worse**, before counting a second localization of the bed, GRM and model.
**Keep all 18 cells in one run.**

## 4. p-values: slim vs upstream 1.4.2

Cells keyed `Region|Group|max_MAF`, Cauchy rows excluded from the statistics.

| phenotype | cells both arms | printed-identical | Pearson on -log10 p | max abs dlog10 p |
|---|---|---|---|---|
| BMI | 9,515 | 633 | 0.999887 | 0.1578 |
| phe_250.2 | 9,495 | 594 | 0.999059 | 0.6160 |
| phe_153 | 9,501 | 663 | 0.999995 | 0.0845 |

(`Pvalue`; `Pvalue_Burden` and `Pvalue_SKAT` are the same story, tabulated in the
notebook.) **Cell sets are IDENTICAL in every phenotype** -- 0 slim-only, 0
upstream-only.

### The two predictions that were wrong

- **#126 (predicted slim-only cells > 0).** Upstream gates burden testing only on
  `--r.corr=1`; slim's contribution was applying it on BOTH region paths, with a
  default of **1, inert by construction**. On a SKAT-O run upstream's gate never
  fires and slim's never refuses. Identical cell sets was the correct prediction.
- **#118/#133 (predicted MAC columns differ).** #118's divergence from upstream is
  on the USER-WEIGHTS path; these runs use default beta weights with K = 1, so
  there is nothing to diverge. `MAC_case`/`MAC_control` were not among the compared
  columns. Observed: `MAC`, `Number_rare`, `Number_ultra_rare` differ in **0** rows.

### What did move, and why -- reconciled

**The burden p is exactly a Wald test on the reported columns in BOTH arms**:
`2 * pnorm(-|BETA_Burden / SE_Burden|)` reproduces `Pvalue_Burden` with max
`|dlog10|` of 1.95e-14 (upstream) and 2.26e-14 (slim) over ~9,478 cells, 0 rows
disagreeing. So the entire p movement reduces to `BETA_Burden` and `SE_Burden`.

Aggregates first, and they MISLEAD: median BETA ratio 0.992427 and median SE ratio
`1/sqrt(1.015046)` = 0.992553 agree to 1.3e-4, which looks like #103's pure
rescale (both scale by `kappa_w = 1'Phi 1 / tr(Phi)`, leaving Z untouched). **Per
cell it does not hold**: `|BETA ratio - SE ratio|` has median 4.3e-3 and max 0.155.
The aggregate agreement was coincidence over a heterogeneous set -- the trap
CLAUDE.md records at #26.

The per-cell chi-square ratio is **median 0.991571, IQR [0.987228, 0.996013]**,
i.e. **Z^2 moved ~0.8%** with ~0.9% spread.

**Attributed to LEDGER #106, and it confirms that row's own caveat.** #106 (sparse
`Sigma` is one triangle; the PCG stops on an absolute residual) quotes 0.05-0.11%
between two scalings and -0.08% to +0.43% against an exact solve -- all on the
bundled N = 1,000 fixture -- and says why those are small: *"Small here only
because `1/W ~ 11` dominates `tau2 x kinship ~ 0.16`; on a quantitative trait
(diagonal `tau1 ~ 1`) or a strongly related cohort the missing triangle is
first-order."* BMI is quantitative and All of Us `afr` is a related cohort, and the
effect is ~0.8%. **This is the first measurement of #106 on real related data**
and should be added to its row.

## 5. Degenerate burden cells: an exact-equality guard at a boundary

In every run, **`SE_Burden` non-finite if and only if `Pvalue_Burden == 1`**, with
no exceptions:

| phenotype | slim | upstream | shared | upstream-only | slim-only |
|---|---|---|---|---|---|
| BMI | 37 | 36 | 36 | 0 | 1 |
| phe_250.2 | 42 | 46 | 42 | 4 | 0 |
| phe_153 | 48 | 48 | 48 | 0 | 0 |

Mechanism: `R/SAIGE_SPATest_Region_Func.R:311-318` zeroes a marker's row and column
of `Phi` when its SPA-implied variance `VarS` equals `Inf` (`zscore^2 / qchisq(p,
log.p = TRUE)` with the denominator underflowing to 0). With #103's
`SE_Burden = 1/sqrt(sum(Phi))` a zeroed `Phi` gives `Inf`; upstream's pre-#103
`sqrt(B)/A` gives `NaN`. Different non-finite value, same p = 1.

**Not a divergence introduced by #103** -- both builds produce it at nearly the
same rate. What differs is WHICH cells: 5 of ~28,500 flip, in both directions.
That is the signature of a guard resting on an EXACT equality (`VarS == Inf`)
sitting near a boundary that #106's ~0.8% variance perturbation can push a marginal
cell across. Neither build is wrong; the predicate is simply not robust to a
perturbation of that size. Example cell, BMI:
`ENSG00000062598 damaging_missense 0.001` -- upstream `SE_Burden 0.0138,
Pvalue_Burden 0.99590`; slim `SE_Burden Inf, Pvalue_Burden 1`.

**A clean exclusion predicate for downstream QC: `SE_Burden` non-finite.**

## 6. Output fingerprint, reconciled three ways

537 genes (not 538), 10,052 rows, rows-per-gene `{19: 486, 16: 48, 18: 2, 14: 1}`
(19 = 18 cells + 1 Cauchy). Missing cells: 151, accounted for exactly -- pLoF
absent in 44+44+46 = 134 cells, damaging_missense in 5+5+6 = 16,
pLoF;damaging_missense in 1. The 48 genes missing exactly 3 cells are those with no
pLoF variants at all (one whole mask, three MAF bins). One gene, `ENSG00000256222`,
is listed in the group file with no markers surviving, and SAIGE names it in the
log (`R/SAIGE_SPATest_Region.R:1765`). Totals agree by two independent routes.

**This fingerprint is the cheapest cross-arm check available** and it was identical
between builds.

## 7. p = 1 saturation and the QC design

**824 of 9,515 cells (8.7%) have `Pvalue` exactly 1** on chr20 BMI. Degenerate
burden cells (section 5) account for only ~4.5% of that, so the large majority is
genuine SKAT-O saturation on low-information genes.

Relevant because the All by All QC (Lu et al. 2026, Supplementary Figure 21) is
built on lambda_GC: computed per phenotype across genes, and per gene across
well-calibrated phenotypes, with hard filters (mean coverage > 10, N_variants >= 5),
recomputation filters (N_variants > 10, CAF > 1e-5; 1e-4 for MID-like and SAS-like)
and an expected-allele-count flag (CAF x n_cases >= 5).

**A point mass at p = 1 breaks lambda_GC**, because it maps to chi2 = 0 and drags
the median down regardless of the tail, and it does so by a DIFFERENT amount per
annotation mask (saturation rate varies with cell size) -- which is exactly the
comparison the QC makes. Two mitigations, and the first is the paper's own:

- **The recomputation filters target it.** `N_variants > 10` and `CAF > 1e-5`
  remove low-information genes, which are the saturating ones. Whether they remove
  the atom AT N = 70,272 rather than at 392,030 is an open measurement.
- **Report lambda at upper quantiles and the ECDF at alpha alongside.** lambda_95 /
  lambda_99 are unaffected by an atom at the bottom of the distribution; the
  fraction of p below alpha with an exact binomial interval is well defined with
  atoms present and is directly interpretable as type-I error. Caveats: at ~9,500
  cells per chromosome there is power at alpha = 0.05 and 0.01 only, and cells are
  NOT independent (nested MAF bins, overlapping masks), so run it within one
  (mask, bin) across genes.

## 8. LEDGER #156 raised and then argued against, in one day

`get_SKAT_pvalue` never passes `isFast` to `SKAT:::Met_SKAT_Get_Pvalue`, so every
SKAT-O cell takes the exact O(q^3) eigendecomposition although SKAT ships an
O(q^2 k) alternative. Raised on the flop ratio (~24x), then MEASURED on the pinned
env: the full call gains only **1.15-1.57x**, and the eigensolve ALONE is SLOWER
below q ~ 3,000 (0.64x / 0.37x / 0.83x at q = 500/1000/2000, reaching 1.77x only at
q = 4000). The flop count ignored BLAS level: dense symmetric `eigen` is blocked
BLAS3, `RSpectra::eigs_sym` is Lanczos over repeated BLAS2 matvecs. With
`k = min(floor(q/10), 100)` the "fast" path is Theta(q^2 k), i.e. still cubic below
q = 1000. Given the q correction in section 3, exome cells sit far below the
crossover: **isFast would most likely make these runs slower.** Row records the
measurement and argues against pursuing it.

## What is NOT done

- **The `q_eff` distribution has never been computed.** One line; it sizes #156's
  scope on real data and is the empirical note that row is waiting for.
- **The non-Cauchy control on the mask experiment.** The first attempt included the
  pooled Cauchy rows, which MUST differ (they combine over a changed cell set), and
  reported 536 mismatches against a prediction of 0. Re-run excluding them.
- **Truncation series** (group file cut to K genes) -- the only route to splitting
  the ~544 s into per-gene work and fixed startup, which are different problems.
- **A third cell-count point.** Two observations fitting two parameters has no
  residual and no validation.
- **PGEN conversion and the reader arm**; **chr1 as the extrapolation control**;
  **#131's ALT-major fraction**; **the formal `compare_pvalues.R` pass with its
  UNPREDICTED = 0 criterion**; **a Firth-off arm** (see below); **the Spot
  preemption rate**, which comes free with the production batch.

## Decisions taken

- **Grid: 18 cells, one run.** Settled by architecture (section 3), not preference.
- **`--is_fastTest`: off for group runs, on for variant runs.**
- **Node: 1 core, 6 GB, pack independent instances; do not thread.**
- **Firth is now irrelevant to the meta-analysis.** The owner is moving from
  inverse-variance weighting to Stouffer with `sqrt(Neff x CAF(1-CAF))` weights.
  Stouffer consumes only the signed Z and the weight, and Firth preserves Z exactly
  (both branches pin `|BETA|/SE = |qnorm(p/2)|`), so `--is_Firth_beta=FALSE` is a
  candidate free saving -- LEDGER's sim50k figure is 1.80x on a binary
  single-variant scan. It also retires the consistency hazard: under IVW a mixture
  of Firth-on and Firth-off cohorts would have distorted the weights; under
  Stouffer it cannot. **Still to check: that BETA signs never disagree between the
  arms**, since Stouffer needs the direction. Note Firth never touches
  `BETA_Burden`/`SE_Burden` under `r.corr == 0` in any case.

---

# CORRECTIONS AND ADDITIONS, later the same day

Everything above stands except where noted here. **Read this section before
quoting any timing from the one above.**

## C1. Every timing above is unreplicated AND unpinned

`submit()` used `--min-cores 1 --min-ram 6`. Those are FLOORS, and dsub resolves
them by building the **smallest matching CUSTOM machine type** -- so each task
got a dsub-chosen **1-vCPU** custom machine, on whatever host generation Batch
had. Between-generation single-thread differences of 20-40% are ordinary, and
GCP network bandwidth scales at ~2 Gbps per vCPU, so a 1-vCPU machine is also
capped at ~2 Gbps.

There was also **no replication**: every arm is n = 1.

So the section-1 speed table (slim 2.15-2.32x faster than upstream) and the
section-3 cost model are both single-replicate, cross-machine comparisons. The
2.2x is large enough to survive that; the smaller differences are not.

**`submit()` now pins `--machine-type n2-standard-2`** (2 vCPU, 8 GB) and records
it in the registry. `--machine-type` and `--min-cores`/`--min-ram` are MUTUALLY
EXCLUSIVE in dsub -- passing both is rejected. **That pinning is a dividing
line: nothing measured before it is comparable with anything measured after**,
because both the CPU and the NIC ceiling changed. Runs predating it lack the
`machine_type` key in the registry, which is the marker.

## C2. The cost model is NOT established

Section 3 fits `1096 = A + 18B`, `912 = A + 12B` to `A ~ 544 s`, `B ~ 30.7 s`,
and concludes per-cell cost is uniform. **A third point refutes the fit.**

`masktest-5nounion` (5 masks x 3 bins = 15 cells, the all-variants union mask
dropped) came in at **1072 s** -- a 24 s saving where the model predicted 92 s.
Taken literally that makes union cells 14.9 ms and `pLoF;damaging_missense`
cells 99.4 ms, i.e. the mask with the LARGEST kernels is the cheap one. Not
credible as a mechanism.

**24 s is 2.2% of 1096 s, and this project's replicate spread is 4-6%.** So the
saving is indistinguishable from zero and the 99.4 ms is an artefact of
attributing all 184 s to the other mask. The two-parameter fit had two
observations, hence no residual and no validation; my claim that it was
"consistent with which cells were removed" was circular.

**What survives**: 12 cells is ~17% faster than 18 cells (one paired
observation), and therefore **the union mask is at most 16.8% of runtime** --
against the ~80% its `sum(q^3)` share implied. That bound holds regardless.

**Outstanding**: two noise-floor replicates of the identical 18-cell config, on
pinned machines. Until those exist, no timing difference below ~5% is readable.

## C3. q_eff measured, and LEDGER #156 is closed by it

`q_eff = Number_rare + 1`: **median 3, p90 21, p99 62, max 482, none above
1,000** over 9,515 cells. The `sum(q^3)` shares barely moved after the
correction (union 81.3% vs 82.6%), so the RANKING survived while the MAGNITUDES
collapsed -- which is precisely why the cubic term never drove the cost. A
median 3x3 eigendecomposition is free.

This is on LEDGER #156, which is now argued against on measurement: the fast
path's crossover is q ~ 2,000-4,000 and at q = 500 its eigensolve is 0.64x, so
`isFast` would slow every cell of this workload. **Do not spend a machine on
that simulation study.**

## C4. The group file is 68.6% non-coding, and those markers are never read

Annotation tally over 388,543 `anno` entries (var/anno line lengths agree for
every gene, so the alignment is sound):

| annotation | entries | share |
|---|---|---|
| `non_coding` | 266,647 | 68.63% |
| `other_missense` | 61,920 | 15.94% |
| `synonymous` | 39,943 | 10.28% |
| `damaging_missense` | 14,749 | 3.80% |
| `pLoF` | 5,284 | 1.36% |

`non_coding` is not in `--annotation_in_groupTest`, and
`R/SAIGE_SPATest_Region.R:2210-2231` removes markers matching no requested
annotation (`rowSums(annoIndicatorMat) == 0`) from `SNP` **and from
`genoIndex`** -- the vector the C++ reader fetches by. So they are dropped in R
**before any decode**.

Two consequences:
- **`.singleAssoc.txt` is restricted to the requested annotation classes.**
  Non-coding variants never appear in it. (Answers a question asked directly.)
- They cost nothing in decode or testing, **but they are still in the `.bed`,
  which dsub localizes in full.**

## C5. The coding subset: 39% of markers, and it validates NUMERICALLY

119,863 unique coding variant IDs of 307,558 in the bim (39.0%). At 79,956
samples a bed record is 19,989 bytes, so: full **5.73 GiB**, coding subset
**2.3 GiB**, coding PGEN **13 MB** (+4.1 MB pvar, 859 KB psam). PGEN's ~180x
against the subset bed is the sparse difflist encoding meeting a spectrum where
~90% of variants are ultra-rare -- a `.bed` spends 19,989 bytes per variant
saying "almost nobody carries this".

**Two traps, both hit:**

1. **plink2 rewrites chromosome codes.** `--make-bed`/`--make-pgen` emitted `20`
   where the input had `chr20`. The wrapper's chromosome check caught it and the
   run died in 108 s. **Without that check `--chrom=chr20` would have matched
   nothing and SAIGE would have completed successfully having tested ZERO
   markers** -- a header-only output, a fast "speedup", and no error.
   `--output-chr chrM` preserves the prefix (the flag is named for what
   chromosome 26 is called, and that choice implies the whole convention).
   **When subsetting chrX, match the original's spelling**: the PAR handling
   recognises `X`, `x`, `23`, `chrX`, `chr23`.
2. **Slim's PGEN CLI is `--pgenFile` / `--pvarFile` / `--psamFile`**, three
   separate flags mirroring bed/bim/fam -- NOT upstream's `--pgenPrefix`, which
   fails with `long flag "pgenPrefix" is invalid`. Slim also ignores
   `--AlleleOrder` for PGEN (a pvar is always REF/ALT-ordered), where upstream
   1.4.2 REFUSES pgen without `--AlleleOrder=ref-first`.

**Validation: byte-identity was the WRONG GATE and I should not have registered
it.** CLAUDE.md says never to assert bit-identity between two different code
paths and to compare statistics at ~1e-12 relative; a different bim with
different marker indices IS a different path. Measured over the 9,515 shared
cells: **18,193 numeric differences, 18,189 below 1e-10, worst 2.2e-7**. No
label, count or `MAC` column moved, which rules out an allele flip or a changed
marker set (either would give order-one differences).

**The 4 outliers are LEDGER #106.** One is `SE_Burden`, a closed form
(`1/sqrt(sum(Phi))`), so a 1.8e-7 move in it means `sum(Phi)` itself moved by
~3.6e-7 -- a million times a rounding difference. #106's PCG stops on an
ABSOLUTE residual (`sumr2 > tolPCG`, 0.02, no `||b||` normalisation), so the
iteration count depends on the input's SCALE and a marginal cell can take one
iteration in one run and two in the other. A perturbation that changes nothing
mathematically flips it.

**So: `SAIGE`'s region output is not a continuous function of its inputs at the
1e-7 level.** Worth knowing before someone spends a day on "why doesn't this
rerun match".

## C6. Overhead outside the script is 13-26% of wall clock

`(end - start) - TOTAL_SECONDS` from `dstat --full` = VM boot + image pull +
localization + delocalization: **183 s for a solo task, 340-411 s when six ran
concurrently** (upstream's consistently 30-60 s higher, consistent with a larger
image). Every one of those was a 1-vCPU custom machine.

**I asserted this was bandwidth contention without establishing it.** Three
candidates remain unseparated: VM provisioning (six at once), image pull
contention, and the network path. The solo overhead implies only ~31 MB/s
overall, against 91.8 MiB/s measured from the notebook -- which suggests boot
and image pull, not transfer, are the bulk of it.

**`cod-bed` vs `cod-pgen` separates them**: same image, same machine, same boot,
2.3 GB vs 18 MB localized. Equal overhead means provisioning dominates and input
size is nearly irrelevant -- in which case the lever is fewer/longer tasks, not
smaller inputs.

A correction on method: **only the OVERHEAD comparison needs matched
concurrency.** `SAIGE_SECONDS` is measured inside the container after
localization, on a VM that shares nothing, so CPU-bound arms can run in
parallel. Serialising them was over-cautious.

## C7. First pinned result

`cod-bed` (coding subset, `n2-standard-2`): **SAIGE_SECONDS 981, PEAK_RSS_MB
3607**, against the 1,096 s / 3,946 MB baseline. **Not attributable** -- the
machine and the input changed together. The noise-floor replicates are what
separate them.

## C8. Decisions and open items changed by the above

- **Grid stays 18 cells in one run** (architecture, unchanged).
- **The all-variants union mask is retained** -- the owner needs it for the All
  by All QC (Lu et al. 2026, Supplementary Fig. 21), whose gene-level lambda_GC
  filters key on `N_variants` and CAF from one (mask, bin) cell. The cost
  argument for trimming it had already collapsed (<= 16.8%, and ~6% under the
  refuted model).
- **Machine pinning, replication, and cost-per-unit** (not wall time) are now
  prerequisites for the threads / multi-trait / batching work, not part of it.
- **Multi-trait step 2** (`--GMMATmodel_varianceRatio_multiTraits_File`) is the
  highest-value untested lever: it shares the genotype read and decode across
  models, which is the half nothing else touches, and each model's output stays
  byte-identical to its own single-model run. The region path requires identical
  step-1 sample ID vectors across models -- **the owner confirms the upcoming
  trait set satisfies this.**

---

# SECOND CORRECTION, same evening: C1 was wrong about pinning, and the noise floor is ~20%

**C1 above says `submit()` "now pins `--machine-type n2-standard-2`" and calls
that a dividing line. It was not pinning anything.** The signature gained the
parameter but the conditional splice never reached `args`, which still carried
`--min-cores`/`--min-ram`. So **every run in this document -- including the ones
C1 calls pinned (`cod-bed`, `cod-pgen`, `rep-18a/b`, `no1pct`, `firth-off`,
`firth-on`, `mask4-pinned`) -- ran on a dsub-chosen 1-vCPU CUSTOM machine.**
There is no dividing line yet. `mask4-pinned` is misnamed and its name should
not be trusted; the `machine_type` column in the registry records the PARAMETER,
not what was passed, so it is unreliable for everything up to this evening.

## The noise floor is ~20%, not 3.18%

`cod-bed` **981 s** and `cod-bed2` **1172 s** are the IDENTICAL configuration on
the IDENTICAL files: **19.5% apart.** Full-bed 18-cell runs across the day span
1096 / 1238 / 1263 / 1278 -- a 16.6% range on the same work.

**The 3.18% from `rep-18a`/`rep-18b` was luck.** Two adjacent replicates do not
estimate a spread; they estimate one difference. Quoting it as a floor was the
same error as the two-point cost fit -- a model with no residual and no
validation.

**So only differences larger than ~20% are readable from single runs.**

| finding | size | survives |
|---|---|---|
| slim vs upstream 1.4.2 | 2.2x | YES |
| PGEN memory, 1634 vs ~3800 MB | 2.3x | YES |
| 12 cells vs 18 cells, 925 vs 1219 mean | 24% | marginally |
| `--is_fastTest` on a group run, +17% | 17% | NO on timing -- but the argument for dropping it is the BYTE-IDENTICAL output, which is exact and unaffected |
| PGEN CPU penalty, ~12% | 12% | NO |
| coding subset faster, ~12% | 12% | NO |
| Firth cost, 1163 vs 1137 | 2.3% | NO -- predicted in advance as likely unmeasurable, and it was |

## The one clean result: per-cell cost is UNIFORM

Three 12-cell runs, two entirely different compositions:

```
masktest-4simple   912 s    4 masks x 3 bins
mask4-pinned       927 s    4 masks x 3 bins  (repeat)
no1pct             937 s    6 masks x 2 bins  -- drops the 0.01 bin, the widest cells
```

A **2.7% spread** across grids that share no cells with each other except by
coincidence, in a population whose noise floor is ~20%. So **cell COUNT drives
cost and cell SIZE does not** -- which retires the whole q / q^2 / q^3 line of
enquiry, and is consistent with per-cell cost being fixed R-side SKAT-O
machinery (one `integrate()` over the rho grid per cell) rather than arithmetic
that scales with the kernel.

## Also measured

- **PGEN is far more reproducible than bed**: `cod-pgen` 1207 and `cod-pgen2`
  1197 (0.8% apart) with peak RSS identical to the megabyte, against bed's
  19.5%. Consistent with bed's variance coming from mmap page residency, which
  is machine- and memory-pressure-sensitive in a way sparse PGEN reads are not.
  Reproducibility is an argument for PGEN independent of speed or size.
- **Firth is fully characterised on real data and is safe to disable** under the
  owner's Stouffer scheme: region file **identical** (`0.00e+00` over every
  column), `.singleAssoc.txt` p-values **identical** (0 of 14,066 differ), BETA
  moved on 817 markers (5.8%), and **zero sign disagreements**. Since Stouffer
  reads only the signed Z, a mixture of Firth-on and Firth-off contributors
  across BRaVa is harmless -- which retires the consistency hazard raised
  earlier.
- **`no1pct` verified as a genuine 12-cell run** (`MAXMAF = 0.0001,0.001` in the
  echo, two distinct `max_MAF` values, 6,879 rows), so the cell-count comparison
  is real and not a mislabelled duplicate.

## What follows

Chasing 10-20% effects with single runs is futile at this variance. Either
accept only large effects, or move to a paired design with several replicates
per arm and interleaved ordering -- which is what `harness/bench_threads.sh`
already does in this repo, for exactly this reason.

The decisions that matter are nonetheless supported, and none of them rests on a
timing difference below the floor: **slim over 1.4.2** (2.2x, and the p-value
comparison attributed), **`--is_fastTest` off for group runs** (byte-identity,
exact), **Firth off under Stouffer** (exact p and sign invariance), **the grid
at 18 cells** (architecture), and **PGEN** (identical output, 2.3x less memory,
far more reproducible; its ~12% CPU cost is below the floor and so is not
established either way).

**In flight at the time of writing**: `mt-dup3` (multi-trait with the same model
three times -- an UPPER bound on multi-trait's benefit, because three identical
models also collapse into the SKAT-O memo, which real distinct phenotypes would
not) and `par3-bed` (three processes on one node -- a LOWER bound, sharing only
the localization). Real multi-trait lies between. The `SKAT-O cells: N asked
for, M evaluated, K served` line on `mt-dup3` quantifies how loose the upper
bound is.

---

# MULTI-TRAIT AND PROCESS PACKING (the evening's runs)

Two arms, designed to bracket real multi-trait without needing three phenotypes
on a common sample set:

- **`par3-bed`** -- three copies of the SAME single-model run as three processes
  in one container. Shares only the task overhead (one localization), so a
  **LOWER** bound on multi-trait's benefit.
- **`mt-dup3`** -- one process, `--GMMATmodel_varianceRatio_multiTraits_File`
  listing the SAME model three times under names `a,b,c`. Shares the genotype
  read and decode, so an **UPPER** bound -- and, in principle, also collapses
  into the SKAT-O memo, which real distinct phenotypes would not.

Both on `n2-standard-4` (16 GB) -- and these are the FIRST runs where
`--machine-type` actually reached `args`, so the first that are what they say.

## The memo was nearly silent, so the upper bound is TIGHT

`mt-dup3` printed `SKAT-O cells: 21915 asked for, 21420 evaluated, 495 served`
-- **2.3% served, not the ~67% predicted** for three identical models. So the
memo is NOT doing the work, and `mt-dup3`'s time is close to what real
multi-trait with distinct phenotypes would give rather than an optimistic
ceiling. That was the run's most useful output: it converts a loose bound into
a near-estimate.

Why the memo barely fires on identical models is not explained here and was not
chased.

## Results

| | `par3-bed` (3 processes) | `mt-dup3` (3 models, 1 process) |
|---|---|---|
| time | 1340 s | 2723 s |
| vs 3 separate runs (~3657 s CPU) | 4020 CPU-s (slightly MORE) | **2723 CPU-s, 25% LESS** |
| peak RSS, concurrent | 9730 MB (3243 per process) | 4975 MB |
| peak HWM, sum per process | 11887 MB (3962 each) | 5625 MB |
| output | three files identical, `01a4b431...` | three files identical, same hash |

**Two memory facts that matter more than the times:**

- **The concurrent peak is 18% below the sum of per-process peaks** (9730 vs
  11887). The copies do not peak together, so packing is governed by 3.24 GB
  per process, not 3.96 GB. This distinction exists only because the wrapper
  was changed this afternoon to sample summed `VmRSS` alongside per-process
  `VmHWM`.
- **Each extra model costs ~554 MB**: 3867 MB for one model, 4975 MB for three.
  So base ~3313 MB on bed, and the per-model slope is what sets how wide a
  multi-trait pass can go.

## Parallel beats multi-trait on cost -- but they COMPOSE, and that is the answer

Multi-trait saves work but uses ONE core, so it loses on cost per trait:

```
multi-trait  n2-standard-2, 2 cores paid, 1 used, 2723 s -> 1815 core-s / trait
3 processes  n2-standard-4, 4 cores paid, 4 used, ~1340 s -> 1340 core-s / trait
             (as run it used only 3 of 4 cores, which was worse than separate
              tasks: fill every core or the idle one eats the saving)
```

Combined -- four processes, each carrying a 3-trait pass, twelve traits on one
node:

```
4 x 2723 s on 4 cores = 10,892 core-s for 12 traits = 908 core-s / trait
vs twelve separate tasks incl. overhead                ~1519 core-s / trait
```

**~40% cheaper, and ONE task instead of twelve.**

## And this is where PGEN finally pays

Four bed processes at 4975 MB is 19.9 GB -- does NOT fit on `n2-standard-4`.
With PGEN's lower base (1634 MB for one model, so ~1080 MB base at the same
554 MB per-model slope) four 3-trait processes is **11 GB**, comfortable in 16.

So PGEN's memory advantage, which bought NOTHING in the single-process regime
(CPU binds first at 4 GB/vCPU), is exactly what makes the combined shape
feasible. Its ~12% CPU penalty is below the noise floor and is the price.

`[ASSUMPTION: the 554 MB per-model slope measured on bed carries over to PGEN.
Plausible -- the increment should be per-sample model vectors, not genotype
machinery -- but UNVERIFIED. Measure it before committing.]`

## Reproducibility: four fragile cells, and only four

Established by accident while checking the above:

- `slim-BMI`, `rep-18a` and `rep-18b` are **byte-identical** (`f3315f02...`) --
  three runs, three machines. The region path IS reproducible run to run.
- `par3-bed` and `mt-dup3` are byte-identical to each other (`01a4b431...`) and
  differ from that reference in **exactly 4 cells at 2.20e-07**.
- The SAME 4-cell / 2.20e-07 signature appeared in three unrelated comparisons:
  full-bed vs coding-subset bed, coding bed vs coding PGEN, and now across
  machine populations.

**So this dataset has exactly four cells sitting on LEDGER #106's iteration-count
boundary, and any perturbation tips them -- a different bed layout, a different
reader, a different CPU. The other 9,511 of 9,515 are bit-stable across all of
it.** The likely cross-machine mechanism is OpenBLAS runtime kernel dispatch
differing by microarchitecture, which #106 amplifies into printed digits; that
is a hypothesis, not a measurement.

**This does NOT undermine today's identity results.** The fastTest pair and the
Firth pair both returned `0.00e+00` and were within-machine comparisons. It is
the same rule CLAUDE.md already states for the goldens -- byte-exactness is a
claim about one platform -- now with a measured blast radius of 4 cells in
9,515.

## Recommended production shape

PGEN input, coding subset, one task per chromosome, several single-threaded
processes per node, each carrying several traits in one multi-trait pass. Every
fixed cost -- VM boot, image pull, localization, R startup, GRM load, genotype
decode -- paid once per node rather than once per trait-chromosome.

Threads stay unused: bounded at ~1.5x, they cost memory (the SKAT-O pool forks),
and filling cores with independent work scales linearly. Reach for them only if
memory ever forces one process per node, in which case the idle cores are
already paid for.

**Two measurements before committing**: the per-model memory slope on PGEN, and
whether the 25% holds with DISTINCT phenotypes rather than three copies. The
near-silent memo says it should, but that is an inference.
