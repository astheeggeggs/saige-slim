# Prediction: null calibration of the cells LEDGER #126 newly admits

Written 2026-09-11 23:50, **before any step-2 run**, against HEAD `72924699`
(clean tree) installed into a session scratch library (`R_LIBS=<scratch>/lib`,
built from `git archive HEAD` with the conda compilers activated by hand, not
`pixi run`, because another session's test tier was on the pinned env).
Fixture: `SAIGE-sim/harness/work/sim50k`, the 2026-09-10 rebuild, seed 7,
N = 50,000 (1,875 families of 4 + 42,500 singletons), 210 genes.
**The fixture has 63 causal genes, not 21** (`sim_meta.json`, `truth_genes.tsv`:
63 rows with `causal = 1`). The "21" in the job description is the 2026-09-01
fixture, which the 2026-09-10 rebuild replaced. The positive control below is
against 63.

Question: `--minGroupMAC_in_BurdenTest` moved 5 -> 1 (#126), so burden now
tests cells whose unweighted MAC (#118 units) is 1-4. Those cells take the
sparse-GRM score test (`src/Main.cpp:3221`,
`useSparseGRM = (MAC > 20.5) ? false : flagSparseGRM`) and, on a binary trait
with a unit-weight set, the ER exact test (`:3230`, `MAC <= max_MAC_for_ER = 4`
and `erAllowedVec[r]`, which is true for No_Weights and `unw_unw` only). Are
they calibrated under the null?

## Design

Two null traits (`y_null`, Gaussian by construction; `y_null_bin`, its top 1%
dichotomised, 500 cases, tau2 non-zero) and one signal trait (`y_binary`, 63
causal genes). Every run at gate 1 AND gate 5, `--nThreads=1`,
`--is_fastTest=FALSE`, `--LOCO=FALSE`, pedigree-truth sparse GRM, models
`step1_<trait>.rda` from the fixture directory (refitted 2026-09-10 under the
#122 guard), collapse cutoff at its default 10.

| arm | group file | cells | weights / path |
|---|---|---|---|
| `plain` | fixture `groupFile.txt` (210 genes) | high, high:medium, high:medium:low x 1e-4, 1e-3, 1e-2 | No_Weights, `--r.corr=1` (BURDEN emitter) |
| `flex` | `flexrv_study/scores/groupFile_informative.txt` | high:medium x 1e-3, 192 FlexRV sets, lofAnno high | transform sets, `--r.corr=1` -- the #118 arm |
| `grid` | fixture `groupFile.txt` | the same 3 x 3 | beta(1,25), SKAT-O path (`region_postprocess` gate) |
| `wplain` | windowed: each gene cut into consecutive 25-marker regions (1,614 regions), same markers, same labels | high, high:medium x 1e-4, 1e-3 | No_Weights, `--r.corr=1` |
| `wgrid` | the windowed file | the same 2 x 2 | beta(1,25), SKAT-O path |

Why the windowed arms exist: on the gene-level file only 11-16 genes per MAC
stratum have a sub-5 cell (table below), so a stratified ECDF there has a 95%
interval of roughly [0, 0.2] at alpha 0.05 and can refuse only gross
miscalibration. A null phenotype cannot see how markers are grouped, so cutting
genes into windows multiplies the sub-5 cells without changing what is being
tested. The gene-level arms are kept because they are the configuration users
run.

`wplain`/`wgrid` run on the two null traits only; `plain`/`flex`/`grid` on all
three traits. 26 runs.

## Statistic, decided before the data

ECDF of the cell p-value at alpha = 0.05 and 0.01 with a Clopper-Pearson 95%
interval, per MAC stratum 1/2/3/4 and pooled; **not lambda_GC and not KS** (the
pooled row has an atom at 1 from #86, and a binary trait at K = 0.01 is
discrete). The atom fractions (`p == 1`, `p == 0`, `NA`) are reported beside
every ECDF.

Rows are not independent: within a region the max_MAF ladder and the
annotation masks nest, and the 192 FlexRV sets share supports. **Primary test:
one cell per (region, MAC stratum)**, the cell with the largest marker support
(ties broken by file order). Every-row ECDFs are printed as descriptive only.

## Predictions

### P1 (pipeline) -- the cell counts, from first principles

Computed from `allcounts.acount` (in-sample allele counts over the 50,000
model samples, OBS_CT = 100,000, no missingness) and the group file, with the
engine's rules: marker kept iff `MAF <= maxMAF` (`find(maxMAFVec >= pr.MAF)`)
and MAC >= 1; cell MAC = sum of kept markers' MACs (#118 units). Script
`predict_cells.py` / `make_windowed.py` in the session scratchpad; the per-cell
table is `predicted_cell_mac.tsv` / `predicted_w25.tsv`.

`plain` (genes with cell MAC in stratum; 0 = empty cell, no row):

| cell | 0 | 1 | 2 | 3 | 4 | >=5 |
|---|---|---|---|---|---|---|
| high 1e-4 | 5 | 14 | 16 | 12 | 11 | 152 |
| high 1e-3 | 5 | 11 | 10 | 8 | 11 | 165 |
| high 1e-2 | 5 | 11 | 10 | 8 | 10 | 166 |
| high:medium 1e-4 | 0 | 0 | 0 | 2 | 2 | 206 |
| high:medium 1e-3 | 0 | 0 | 0 | 1 | 1 | 208 |
| high:medium 1e-2 | 0 | 0 | 0 | 1 | 1 | 208 |
| high:medium:low, any | 0 | 0 | 0 | 0 | 0 | 210 |

140 sub-5 cells in total (36 / 36 / 32 / 36 at MAC 1 / 2 / 3 / 4); after
taking one per (gene, stratum) about 14 / 16 / 12-14 / 11-13.

`flex` (high:medium, 1e-3), the sets whose support is predictable without the
score: `lof_unw` 11 / 10 / 8 / 11 (the `high 1e-3` row above, one quantity two
ways); `lof_MAC1` 41 / 32 / 32 / 14; `lof_MAC3` 18 / 21 / 20 / 18;
`lof_MAC5` 16 / 17 / 15 / 14; `lof_MAC10` 14 / 16 / 12 / 11 (= `high 1e-4`,
since MAC <= 10 is MAF <= 1e-4 at N = 50,000); `unw_MAC1` 5 / 8 / 14 / 23;
`unw_MAC3` 0 / 0 / 3 / 4; `unw_unw` 0 / 0 / 1 / 1. The score-shaped sets
(`cbrt*`, `sbeta*`) zero the markers whose score is 0 and are not predicted.
Every set of a cell must print the plain unweighted MAC of ITS support, an
integer; a fractional MAC anywhere is a #118 regression.

`wplain` (W = 25, 1,614 regions):

| cell | 0 | 1 | 2 | 3 | 4 | >=5 |
|---|---|---|---|---|---|---|
| high 1e-4 | 663 | 334 | 251 | 114 | 70 | 182 |
| high 1e-3 | 601 | 303 | 230 | 105 | 63 | 312 |
| high:medium 1e-4 | 82 | 93 | 117 | 132 | 136 | 1054 |
| high:medium 1e-3 | 59 | 64 | 67 | 86 | 84 | 1254 |

**Falsifies the pipeline:** any cell whose printed `MAC` differs from the
table's value (units, MAF-filter boundary or missingness accounting); any
non-integer `MAC`; row counts off the table.

### P2 (pipeline) -- what gate 5 removes

The gate-5 output of every arm must be the gate-1 output with EXACTLY the rows
at `MAC` in {1, 2, 3, 4} deleted: every other cell row byte-identical (the gate
is read after the test, it cannot touch a surviving cell's arithmetic), the
`Cauchy` row of a region changed only where that region lost a cell, and the
`.pooledTests.txt` sidecar at gate 5 naming, per region, a refused-cell count
equal to the number of sub-5 rows that region has at gate 1. **Falsifies the
pipeline:** a row present at gate 5 and absent at gate 1; a shared cell row
whose p-value differs between gates; a refusal count that does not reconcile.

### P3 (the hypothesis, quantitative) -- `y_null`

`y_null` is polygenic + Gaussian noise and the model's residual is close to
Gaussian, so a burden score over 1-4 alleles is (a weighted sum of) 1-4
residuals and its normal-approximation p is uniform. Prediction: **ECDF(0.05)
and ECDF(0.01) inside their 95% intervals in every stratum and pooled, with no
common direction across strata**; on `wplain high 1e-4` (n = 334 / 251 / 114 /
70) the intervals at alpha 0.05 are roughly +-0.025 / 0.03 / 0.04 / 0.05 wide.

Expected magnitude of a REAL miscalibration if the sparse-GRM variance were
mis-scaled: the two `null` ratios of this model are 0.54% apart (handoff
2026-09-11 evening), so a wrong band moves SE by ~0.3% and is invisible here;
what this check can see is the estimator switch itself
(`useSparseGRM`), which changes the variance by the within-family
covariance -- for a singleton carrier nothing, for a carrier with a genotyped
relative a term of order the kinship times tau. With 7.5% of the sample in
families, a wrong switch would show as a shift of the SAME sign at every
stratum, growing with MAC.

**Falsifies the hypothesis (cells miscalibrated, code doing what it says):**
ECDF(0.05) outside its interval in the same direction at 3 or more strata, or
pooled ECDF(0.05) outside [0.035, 0.065] on `wplain`, while P1 and P2 hold.
**Falsifies the pipeline instead:** P1 or P2 failing; p-values `NA`, exactly 0
or exactly 1 in bulk on admitted cells; or a departure so large (pooled
ECDF(0.05) > 0.15 or < 0.01) that a unit error, not an estimator property, is
the parsimonious cause.

### P4 (the hypothesis, binary) -- `y_null_bin`, K = 0.01

The p-value of a cell with m carrier alleles is DISCRETE: the score is a sum of
m terms (y_i - mu_i) with mu_i ~ 0.01, so the only outcomes are "k of the m
carrier alleles are in cases". ER (plain / No_Weights, MAC <= 4) evaluates that
distribution exactly; SPA (the 191 non-unit FlexRV sets) approximates it. With
one case allele among m, the two-sided p is ~ P(X >= 1) = 1 - 0.99^m, which is
0.010 / 0.020 / 0.030 / 0.039 for m = 1..4; with none, p is ~1. Prediction:
**ECDF(0.05) per stratum ~ 1 - 0.99^m and in every case at or BELOW 0.05 within
its interval; ECDF(0.01) between 0 and ~0.01 m** (a MAC-1 cell with a case
carrier sits exactly at p ~ mu_i, on either side of 0.01 depending on the
carrier's covariates). Families: a MAC-2 cell whose two alleles are a sib pair
has correlated case status, a second-order effect on these rates.

This is the memory's point (`lambda-gc-invalid-with-atoms`): the binary null
at low prevalence reads CONSERVATIVE by construction and that is the data. What
would be wrong is the other direction. **Falsifies the hypothesis:** ECDF(0.05)
whose lower 95% bound exceeds 0.05 in any stratum (anti-conservative), or
ECDF(0.05) in a stratum more than a factor ~2 off 1 - 0.99^m on the `wplain`
counts. **Falsifies the pipeline:** P1/P2; every sub-5 p exactly 1 or `NA`
(the case-carrier cells must exist: expect ~1% of MAC-1 cells to have p ~ 0.01).

### P5 (positive control) -- `y_binary`, 63 causal genes

Detection = pooled `Cauchy` row p <= 0.05 / 210 on the gene-level arms. Reference
(sweep of 2026-09-10 at collapse 10, pre-#118 code `5dc679c3`): `flex192`
informative 13 / 63, `unw` 8 / 63, `grid` 11 / 63. #118 moved ~9% of burden
p-values by up to 20%, so HEAD at gate 5 is predicted within +-2 of each.

Gate 1 adds sub-5 cells to a gene's Cauchy combination. Adding one cell with a
non-small p to a combination of n cells moves the pooled p by a factor of at
most about (n+1)/n, so on `flex` (n ~ 180 admitted sets) and `plain` (n <= 9)
a gene can cross the threshold only if it already sits within that factor of
it. Prediction: **detections at gate 1 >= detections at gate 5 minus the
number of gate-5 detections within a factor (n+k)/n of the threshold, k = the
cells the gene gained**, and in practice equal. **Falsifies the pipeline:** a
gene detected at gate 5 and lost at gate 1 whose pooled p moved by more than
the CCT arithmetic over its added cells allows (recomputed from the cell rows:
one quantity two ways), or a gene GAINED at gate 1 through a sub-5 cell with p
< 1e-4 (a 1-4 allele cell cannot carry that much evidence at K = 0.01 unless
every carrier is a case, whose probability is 0.01^m).

## What is NOT being tested

The variance-ratio BAND (0.54% apart on this model, invisible); a native
missing-genotype case (the fixture has none); `--r.corr=1` against `--r.corr=0`
p-values on the same cell (never compared, LEDGER #126 -- only the refused SET
may be compared across paths, and it is: at gate 5 `plain` and `grid` must
refuse the same (region, group, max_MAF) cells).

---

# Outcome (2026-09-12, appended after the runs; nothing above was edited)

26 runs planned, 30 run (the reduced windowed SKAT-O arm `wgridr` was added,
below). Library: scratch install of `72924699`, `SAIGE.so` md5
`fd02953542c492c2b17c7d4f3cd94f0b`, which is also the pinned library's. Raw
analysis: `harness/sim/calib126/results-2026-09-12.txt`; scripts
`harness/sim/calib126/`.

## P1 -- held, after one correction to the RULE, not to the engine

The table above predicts the cell MAC as the plain allele sum. The engine's
count is the sum over kept NON-ultra-rare markers plus, for the ultra-rare ones
(MAC <= 10), the sum over samples of the per-sample MAX of their genotypes --
the collapse, LEDGER #118's unweighted twin. The two differ by the number of
alleles a sample carries beyond its first ultra-rare one in the cell. On the
gene file that is 387 of 1,875 cells (up to 80 alleles on `GENE0001
high:medium:low`); on the windowed file 6 of 5,051 (by 1-2); on the sub-5
cells, none. Verified from the genotypes on 7 cells first (10, 10, 11, 33,
1201, 703, 5108: each equals the printed MAC exactly), then the whole table
recomputed under that rule (`predict_collapsed.py`): **0 mismatches on 1,875
gene cells, 5,051 windowed cells and 4,128 FlexRV set rows** (the 22 sets whose
support the score does not shape), no non-integer MAC, no predicted non-empty
cell missing, in every arm and every trait. The sub-5 counts were exactly the
predicted 36 / 36 / 32 / 36 (gene file), 2,249 in the windowed strata
(794 / 665 / 437 / 353 rows), 1,849 FlexRV set rows.

## P2 -- held

In all 22 (trait, arm) pairs with both gates: 0 rows only at gate 5; every
row only at gate 1 has MAC in 1..4; no sub-5 row survives at gate 5; **0 shared
cell rows differ**; the Cauchy rows that change are exactly the regions that
lost a cell. The gate-5 sidecar's refusal counts reconcile with the distinct
sub-5 cells of the gate-1 run in every region but one (`GENE0138`, sidecar 2
against 1 distinct: its two identical `high` cells at max MAF 1e-4 and 1e-3
were both evaluated and refused, where elsewhere a refused cell's identical
successor is not re-counted -- a #127 bookkeeping detail, not a result).

## P3 -- held: the quantitative null is calibrated in every stratum

Primary cells, one per (region, MAC stratum); ECDF at 0.05 with the 95%
interval; alpha 0.01 in brackets.

| arm | MAC 1 | MAC 2 | MAC 3 | MAC 4 | pooled 1-4 | MAC 5-20 reference |
|---|---|---|---|---|---|---|
| `wplain` (n 392/353/240/201) | 0.059 [.038,.087] | 0.057 [.035,.086] | 0.054 [.029,.091] | 0.050 [.024,.090] | 0.056 [.043,.070] (0.009) | 0.059 (0.011) |
| `wgrid` gate 1, `Pvalue_Burden` and SKAT-O `Pvalue` | identical to `wplain` in every stratum (a sub-5 cell is one pseudo-marker) | | | | | |
| `wgridr` (n 268/254/148/117) | 0.060 | 0.063 | 0.047 | 0.051 | 0.057 [.042,.076] | 0.059 |
| `flex` (n 44/46/60/50) | 0.045 | 0.065 | 0.050 | 0.060 | 0.055 [.028,.096] | 0.057 |
| `plain` = `grid` (n 14/16/14/13) | 0.000 | 0.062 | 0.143 | 0.077 | 0.070 [.019,.170] | 0.052 |

No stratum outside its interval, no common direction, no atom (0 rows at
p = 1 or 0, 0 NA). The hypothesis falsifier did not fire; neither did the
pipeline one.

## P4 -- held: the binary null is conservative by discreteness, never liberal

| arm | MAC 1 (exp .010) | MAC 2 (.020) | MAC 3 (.030) | MAC 4 (.039) | pooled | MAC 5-20 |
|---|---|---|---|---|---|---|
| `wplain` (ER path) | 0.020 [.009,.040] | 0.023 [.010,.044] | 0.025 [.009,.054] | 0.025 [.008,.057] | 0.023 [.015,.033] | 0.058 (0.004 at 0.01) |
| `wgridr` | 0.026 | 0.020 | 0.027 | 0.026 | 0.024 [.015,.037] | 0.055 |
| `flex` (SPA, non-unit sets) | 0.000 | 0.022 | 0.000 | 0.040 | 0.015 [.003,.043] | 0.061 |
| `plain` = `grid` | 0/14 | 0/16 | 0/14 | 0/13 | 0.000 [0,.063] | 0.065 |

Every upper bound is below 0.05 at MAC 1-2 and every point estimate is below
0.05 everywhere; the p-value distribution is the predicted two-point shape (the
1st percentile 0.015 / 0.0095 / 0.019 / 0.037 at MAC 1-4 -- the case-carrier
cells -- and a median of 0.93 / 0.90 / 0.87 / 0.85 for the rest). MAC 1 reads
0.020 against 0.010 (8 of 392; P(>= 8 | 0.010) ~ 0.05, and 0 of 44 and 0 of 14
on the other arms), MAC 3-4 read below expectation; neither is significant and
neither is in the liberal direction. At alpha 0.01 the MAC 5-20 reference on
this trait is 0.004 [.001,.009], the same conservatism one stratum up.

## P5 -- held, with one loss that the pre-registered arithmetic explains

63 causal genes, threshold 0.05 / 210 on the pooled row:

| arm | gate 5 | gate 1 | 2026-09-10 sweep, pre-#118 (like-for-like only for `flex`) |
|---|---|---|---|
| `plain` | 14 causal, 0 non-causal | 14, 0 | 8 (`unw`, ONE cell: high:medium x 0.001) |
| `flex` | 11, 0 | 11, 0 | 13 (informative, same cell, same 192 sets) |
| `grid` | 16, 1 (`GENE0094`) | 15, 1 | 11 (`grid`, SIX cells: high, high:medium x 3 MAFs) |

The one loss is `GENE0183` on `grid`: pooled p 2.17e-4 at gate 5 (9% inside
the threshold 2.38e-4) and 2.62e-4 at gate 1, ratio 1.21 against the (n+k)/n
bound 9/6 = 1.5 -- the three added `high` cells are one MAC-2 cell copied
across the three max MAFs at p 4.5e-4, and `(6 x 1470 + 3 x 707) / 9` gives
the printed 2.62e-4 exactly. Over every gene with gate-5 pooled p < 0.01 the
ratio never exceeds the bound by more than 1.1% (a cell at p > 0.5 pulls
harder than the bound assumes, by that much). The two sub-5 FlexRV rows at
p < 1e-4 on the signal trait (`GENE0010`, `GENE0150`, `lof_MAC3`, MAC 4 with
3 case alleles) are both causal genes and both already detected. **CORRECTION
(2026-09-12 02:50):** the `plain` and `grid` "references" are NOT like-for-like
and the +-2 prediction did not apply to them -- the sweep's `unw` is one cell
and its `grid` six, where these arms run nine (adding `high` alone and
`high:medium:low`), so 14 against 8 and 16 against 11 compare different cell
grids, not two code versions. Only `flex` (11 against 13) is the same
configuration, and it is inside +-2. The like-for-like numbers come from the
re-run of the sweep at HEAD (`prediction-flexrv-rerun-2026-09-12.md`).

## What the check found that it was not looking for

1. **LEDGER #128 (new, OPEN).** Both windowed SKAT-O runs at gate 5 died at
   region 4 of 1,614 (`GENE0001_w04`, two cells at MAC 1): a region with two or
   more cells all refused kills the run in `regionFinish`'s cross-region
   `rbind`. Reproduced on the bundled fixture 8 of 8 (`GENET var rs28 rs32`
   with two annotation groups), independent of order and collapse cutoff; one
   cell escapes; `--r.corr=1` escapes. The SKAT-O gate-5 control therefore ran
   on the 1,254 windowed regions whose largest cell has MAC >= 5 (`wgridr`),
   where P2 held. Unreachable at the default gate 1.
2. **#116 re-observed**: `GENE0151_w08` and `GENE0156_w03` are analysed and
   produce no cell (their labelled markers are monomorphic or above max MAF);
   the burden path writes an `NA` Cauchy row and the SKAT-O path writes
   nothing, so the two paths already disagree on an empty region.
3. **The pooled Cauchy row at alpha 0.05 reads 0.06-0.08 on the quantitative
   null in every arm at BOTH gates** (`wgridr` 0.064 [.051,.079] / 0.067
   [.054,.082]; `wplain` 0.058 / 0.063; gene-level `plain` 0.076, `flex`
   0.081, n = 210) while alpha 0.01 is on target (0.011-0.019). Not a gate
   effect, and [ASSUMPTION: the known behaviour of the Cauchy combination of
   nested, positively dependent cells at a large alpha, accurate in the tail
   and mildly liberal at 0.05] -- not tested here; a separate item.
4. **A driver trap**: `echo "exit=$?"` placed after a `$(date)` substitution
   reports date's status. `harness/sim/flexrv_study/run.sh` had it (every
   `exit=` it ever logged was 0) and is fixed; the calibration driver captures
   `rc` first.
