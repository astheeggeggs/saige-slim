# Prediction: FlexRV power fixture p2 (2026-09-14, written BEFORE the sweep)

The 2026-09-13 sweep re-run left one scientific question the old fixture
could not answer (`handoff-2026-09-13-flexrv-sweep-rerun.md`, section 5):
FlexRV beat the unweighted burden by the same margin with the score PERMUTED
within gene as with it informative, so either the gain is the 144/192-set
Cauchy combination itself, or the fixture was too weak to separate the two.
Section 6e's design constraint for the next fixture: **put the causal signal
where a burden test can resolve it**. This file records the fixture built to
that constraint and what it is expected to show, so the sweep is compared
against a prediction and not read into.

## What was built (`harness/sim/flexrv_study/make_power_pheno.py`)

Phenotypes only. Every genotype, the pedigree-truth GRM, the PCs, the group
file and the three score group files are the 2026-09-10 fixture's, untouched.

| knob | value | why |
|---|---|---|
| causal genes | 60 of 210, drawn among the 179 whose visible cell (`high:medium`, MAF <= 0.001) has >= 30 minor alleles | the cell every study run tests |
| P(causal) within a gene | `0.9 * sigmoid(10 * (pct(a_lat) - 0.5))`, `high` at the top | a_lat is what the informative score ranks, so the score is informative BY CONSTRUCTION; the realistic score (a_lat + N(0, 1.5²)) partly, the permuted one not |
| effect | one beta per gene, same sign (directional) | the burden is the test under study; the causal-only burden is the oracle |
| per-gene NCP | log-uniform on [10, 80] on the causal-only burden, delivered exactly on the realised burden variance | the old fixture set NCP on the whole-gene burden and 27/62 genes had it below MAC 20 |
| background | `y_null` (pedigree polygenic + covariates + noise, var 1.007) scaled by sqrt(1 - h2) | keeps the GRM structure the null models were fitted on; no new polygenic draw is possible without regenerating |
| binary | prevalence 0.10 (5,000 cases), own liability at h2 / 0.342 | at K = 0.01 a burden test on 500 cases has 0.072 of the quantitative NCP and the equal-NCP liability h2 (0.55) breaks the budget; at K = 0.10 it is 0.116 |
| binary null | `y_null` thresholded at K = 0.10 (`y_null_bin_p2`) | the 2026-09-08 recipe; permuting is ruled out (LEDGER #104) |
| seed | 20260914 | own stream; the fixture's seed 7 is not consumed |

Realised: chr20 h2 0.0397 (quant), 1,182 causal variants, causal minor
alleles per gene median 91 (q10 20, min 3), causal variants per gene median
16 of a visible median 30, mean cor(informative score, causal) 0.67.
`NCP_binary = NCP_quant` by construction (`truth_genes_p2.tsv`).

## Controls that already passed (generation)

- Ladder delivered: max relative error of `NCP_quant` against `NCP_target` is 0.
- Two ways to one number: sum of per-gene h2 = 0.0397 = var of the total
  genetic value (ratio 0.999; LD between causal genes is negligible).
- `cor(y_quant_p2, y_null)` = 0.980 = sqrt(1 - h2), as designed.
- bim order = `variants.tsv.gz` order (asserted every 997th row); fam IID
  order = phenotype row order (asserted); 0 missing calls.

## The prediction (1-df non-central chi-square, Bonferroni 0.05/210, chi2 > 13.50)

Expected detections among the 60 causal genes, quantitative = binary by
construction. `_s1` is the score itself as weight, `_s3` its cube, `_best` the
max over exponents 0.2..20 -- an UPPER bound on what a Cauchy combination
over the transform sets can return, not a prediction of it.

| burden | NCP median | E[detections] |
|---|---|---|
| unweighted (`unw`) | 14.1 | **33.7** |
| causal-only oracle | 29.5 | 49.5 |
| informative s1 / s3 / best | 18.7 / 20.8 / 24.0 | 39.0 / 41.5 / **43.5** |
| realistic s1 / s3 / best | 16.1 / 17.6 / 20.2 | 37.7 / 37.8 / 40.3 |
| uninformative s1 / s3 / best | 14.5 / 12.9 / 19.4 | 32.8 / 30.4 / **36.7** |

What the sweep is expected to show, stated as inequalities so a miss is a miss:

1. **`unw` detects 30-37 genes** on `y_quant_p2` (binomial sd of the count
   is ~3.5). Outside that band the fixture or the estimator is wrong, not the
   method: halt and reconcile per gene against `pred_unw`.
2. **`flex192` informative > `unw` by 5-10 genes** (39-43 against ~34).
   The paired per-gene comparison is the statistic; the counts are noisy.
3. **`flex192` uninformative lands between `unw` and 37**: the Cauchy
   combination of 144 sets (MAC100/200/500 are dropped at N = 50k, LEDGER
   #142) pays for the selection, so it should sit near `unw` (33 +- 3). If
   it matches the informative arm the gain is the set combination, not the
   score, and the 2026-09-13 finding stands at a resolution where it means
   something. If it sits at `unw` the finding was an artefact of the old
   fixture's power.
4. **Monotone**: informative >= realistic >= uninformative, each gap >= 2
   genes in expectation (39-43 / 38-40 / 33-37).
5. **The oracle row (min over sets, Bonferroni by the set count) is at most
   49** and above every Cauchy row.
6. **Null calibration unchanged**: `y_null` and `y_null_bin_p2` inside the
   binomial interval at alpha 0.05 and 0.01 on every line, and the 150
   non-causal genes of the signal traits likewise. The models are new fits
   (`step1_y_quant_p2`, `step1_y_binary_p2`, `step1_y_null_bin_p2`); the
   binary fits must report tau2 > 0 (a zero is LEDGER #121/#122's shape,
   which is data at K = 0.01 and should not recur at K = 0.10).
7. **Synonymous arm (low class only)**: no causal variant is in `low`, so it
   is a per-gene null on the signal traits: inside the interval.
8. **Binary ~ quantitative** in every count (equal NCP), within noise; SPA
   on 5,000 cases at MAC >= 20 is well inside its regime.

Collapse: 10 only (owner's decision 2026-09-13/14: the collapse stays at 10
for both traits, which is upstream's default). At MAC <= 10 of 100,000
alleles the collapse is a no-op on this fixture (6d of the sweep handoff),
so a c0 column would add nothing.

Cost: 4 traits x (3 arms x 4 FlexRV methods + 4 plain + 2 synonymous) = 72
runs; the 2026-09-13 timings (8-64 s each at 1 thread, load ~1) put it at
~25 min at 3 in parallel. The box carries an external load of ~6 today, so
wall times are not comparable to the 2026-09-13 ones and are not measured.
