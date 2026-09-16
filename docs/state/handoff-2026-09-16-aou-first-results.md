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
