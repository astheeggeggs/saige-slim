# Vectorised weighted burden as a screen: measured, not built (2026-10-01)

**Question.** A collaborator computes each weight set's burden by multiplying
once-computed per-marker score statistics by the weights, instead of one test
per set. Do we do that, and could the BURDEN path (`--r.corr=1`, FlexRV's)?

**Answer.** The SKAT-O path already does (`R/SAIGE_SPATest_Region.R`,
`wStatVecList`, `VarMat * w w'`). The BURDEN path does not: it runs
`getMarkerPval` per (cell, set) on the weighted sum vector (`src/Main.cpp`,
`burdenTestFn`). A screen -- vectorised Z for every set, the existing full test
only above the SPA cutoff or on exact-test cells -- is feasible and would skip
~91% of the per-set tests under the null. **Nothing in `src/` or `R/` was
changed; this is a measurement.** HEAD `7b7b447a`, installed `SAIGE.so` md5
`8654f0c7`, load average 4-5 throughout (shares and counts are robust to that,
seconds are not).

Fixture: sim50k (N = 50,000, SIMULATED, so counts are written plainly), 210
genes, FlexRV 192 sets, one annotation group, max MAF 0.001, sparse GRM,
1 thread. Scripts, logs and dumps: `harness/work/flexrv50k/vecscreen/`
(untracked work directory).

## 1. How much there is to win

| arm | wall | C++ total | per-set tests (`emitTest`) | share of C++ | share of wall |
|---|---|---|---|---|---|
| quantitative (`y_null`), rep 2 | 22.2 s | 18.0 s | 11.6 s | 65% | 52% |
| binary (`y_null_bin`), rep 2 | 24.1 s | 19.3 s | 13.1 s | 68% | 54% |

`BURDEN cells: 30167 asked for, 23498 evaluated, 6669 served` in every run. The
stop rule was 30%; both arms are above it. Binary costs only 13% more than
quantitative per test (~0.5 ms each), so the cost is the O(N) score on the
scattered vector, not the saddlepoint.

## 2. Reconciliation: w'S and w'Phi w against the per-set test

Method: the SKAT-O path run with `get_SKAT_pvalue` (quantitative) or
`get_newPhi_scaleFactor_traitType` (binary, i.e. BEFORE the saddlepoint rescale)
traced to dump `sum(Score)` and `sum(Phi)` per cell; the BURDEN run's
`BETA_Burden`/`SE_Burden` give T = BETA/SE^2 and V = 1/SE^2. **The comparison is
limited to printed precision (6 significant digits, ~1e-5), not the 1e-12 the
plan named**: the BURDEN path prints no more.

**The SKAT-O path's own `Pvalue_Burden`/`SE_Burden` columns cannot be used for
this**: SKAT evaluates "rho = 1" at 0.999 and `SE_Burden` is back-derived from
that p, so on GENE0001 Z reads 0.2038 against the BURDEN path's 0.1977 while
T and V agree to 0.1%. Use the kernel's inputs, not its outputs.

Controls (K = 1, unweighted set, quantitative):

| control | T'/T - 1 (max abs) | V'/V - 1 | reading |
|---|---|---|---|
| positive: no sparse GRM, variance-ratio bands EQUALISED | 1.2e-5 | max 9.7e-6 | agrees at printed precision |
| ratio control: no sparse GRM, bands as fitted (1.0892 / 1.0951) | 1.2e-5 | max 5.36e-3 | predicted bound 1 - 1.0892/1.0951 = 5.36e-3: met exactly |
| production shape: sparse GRM | 1.2e-5 | median 1.6e-3, max 3.0e-2, both signs | a different variance ESTIMATOR, below |
| negative: sums compared with the NEXT weight set's row (k192) | median 0.49 | median 0.44 | breaks, as it must |

**Why V differs with the sparse GRM (observed in `contextForMAC`,
`src/SAIGE_test.cpp`):** `useSparseGRM = (MAC > 20.5) ? false : flag`. A marker
at MAC <= 20.5 gets its VarMat row from the sparse GRM; the CELL's burden test,
whose MAC is in the hundreds, gets the variance-ratio approximation. So
1'Phi 1 and `varT` are two estimators of one variance, not one number computed
two ways. Pseudo-marker-only cells agree at printed precision (max 9.9e-6 over
16,384 cells) because there the cell IS the one row. This is the finding that
decides the design: T is linear and exact; V is exact only if the vectorised
path uses the cell's estimator.

k192, quantitative, sparse GRM (30,153 of 30,167 cells joined to the dump):
T agrees at printed precision in every cell, including custom/transform
pseudo-markers (they carry a per-set row). V'/V - 1: pseudo-only max 9.9e-6;
real + pseudo cells median 1.5e-3, 99th percentile 2.1e-2, max 3.6e-2. Max
|Z' - Z| = 0.038.

k192, binary, sparse GRM (all 30,167 cells): sign of Z' and Z agree in every
cell; among cells below the screen, |p'/p - 1| is median 8e-8, 99.9th
percentile 2.2e-3, max 3.5e-3. Tighter than the quantitative arm because this
model's two estimators are closer (null ratio 0.994 against 1.09).

## 3. What a screen would escalate (null phenotype)

| | quantitative | binary |
|---|---|---|
| \|Z\| > 2 on the current path | 1,392 | -- (not printed on this path) |
| \|Z'\| > 2 | 1,389 | 1,703 |
| decisions that differ at exactly 2 | 3 | -- |
| \|Z\| > 2 but \|Z'\| <= 1.8 (missed with the margin) | **0** | **0** (smallest p among non-escalated 0.072; none differs from the normal p by more than 0.35%) |
| escalated at \|Z'\| > 1.8 | 2,105 (7.0%) | 2,157 (7.2%) |
| ER-eligible (MAC <= 4), not already escalated | n/a | 477 |
| cells that skip the full test | 93.0% | 27,533 (91.3%) |

Expected before the run: ~5% past the cutoff plus the ER cells. Observed 7.2%
at the 1.8 margin (P(|Z| > 1.8) = 7.2% under the null: as expected) plus 1.6%.
The 99 exact-by-weight cells of this run (#172) need the carrier vector to be
identified and are not subtracted above; they are at most 0.3% more.

**NOT measured:** the power phenotype (escalation will be higher in genes with
signal, by construction); any real-data cell; more than one thread; a trait
with pinned samples.

## 4. Verdict and what building it would take

Go, as a screen, behind a switch -- but it is an owner's decision because it
moves output:

1. Z'_r = w_r'S / sqrt(V'_r) for every set from the per-marker statistics.
2. Escalate to the UNCHANGED per-set test when |Z'_r| > cutoff - margin, the
   cell is ER-eligible or exact-by-weight routed, or Z' is NaN. Those rows stay
   byte-identical.
3. The open choice is the variance of the NON-escalated rows:
   (a) report from 1'Phi 1 -- cheapest, but p moves by up to ~2% in Z on a
       quantitative model with relatedness (every multi-marker FlexRV golden
       moves);
   (b) rebuild the cell's own estimator from carrier-space cross-products
       (ratio x non-sparse g'Pg) -- the equalised-ratio control shows this is
       exact at printed precision, so goldens would move only at the ULP
       level, and it needs a second m x m product per region.
   Recommendation: (b). It keeps "the screen changes no reported number" as a
   testable claim; (a) silently swaps the estimator.
4. The run must print `BURDEN screen: N cells, M escalated`, and a kill switch
   must restore today's path byte for byte.

Ceiling on the gain: ~91% of 52-54% of wall, i.e. roughly 1.9x at 192 sets on
one thread, more as the set count grows (the per-set cost goes from O(N) to
O(m^2)). Not a reason to build it today; a reason to build it before the
weight grid grows.

## 5. For the collaborator

"Reproduces current SAIGE" over all cells is dominated by cells below the SPA
cutoff, where SAIGE reports the normal p and agreement is by construction. The
informative check is binary cells above the cutoff on an unbalanced trait,
split by how concentrated the weight is. And if their comparison used the
SKAT-O path's `Pvalue_Burden`, that column is SKAT at rho 0.999, not the
burden test.

## 6. Added 2026-10-02: a diagonal Phi, and what the collaborator's vectorisation is

**A diagonal Phi does not rescue the variance.** The estimator mismatch of section 2
is on the diagonal too, and the off-diagonals matter in the tail. Measured on the
same fixture, K = 1 unweighted set, the 186 cells with more than one Phi row (rows
are the real markers above the collapse cutoff plus one pseudo-marker: median 3,
max 25), `sum(diag(Phi)) / sum(Phi) - 1`:

| arm | median | min | max | worst Z inflation |
|---|---|---|---|---|
| quantitative (no sparse GRM, equalised ratios) | +9.2e-4 | -0.248 | +6.1e-3 | +15% |
| binary (sparse GRM, before the saddlepoint rescale) | +9.8e-4 | -0.297 | +1.1e-2 | +19% |

So the typical gene loses nothing and the worst understate the variance by a
quarter, which is the anti-conservative direction. As a screen only, a diagonal
would need the margin widened from 1.8 to about 1.65 and could not supply the
reported p. Not measured: the collapse off (many more rows), real LD.

**What the collaborator's vectorisation is** (`~/Repositories/rvas_aou` at
`ee8b236`, 2026-08-13, `ukb_weighted_burden/genebass_burden.py`; a newer SAIGE
fork, if one exists, was not read). From per-variant summary statistics:
`U_i = BETA_i / SE_i^2`, burden score `sum w_i U_i`, variance `sum w_i^2 / SE_i^2`
(diagonal), p from a 1-df chi-square, all 192 sets in one `hl.agg.array_sum`. No
burden-level saddlepoint, exact test or collapse. Against this branch: the linear
score we have (SKAT-O path) and could bring to the BURDEN path as the screen above;
the diagonal variance and the per-variant-adjusted tail we should not take. The one
thing worth borrowing is the shape -- all sets' Z' from one matrix product
(`W'S`, `diag(W' Phi W)`) rather than a loop over sets. The summary-statistic side
of this discussion (count-based null, concentration guard, rarest-variant check)
is recorded in `~/Repositories/genebass-burden/README.md`, where it applies.
