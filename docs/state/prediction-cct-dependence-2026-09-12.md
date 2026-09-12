# The pooled Cauchy row's size at alpha 0.05 (2026-09-12)

Closes the open observation from `docs/state/prediction-126-null-calibration.md`
item 3: "the pooled Cauchy row at alpha 0.05 reads 0.06-0.08 on the
quantitative null in every arm at BOTH gates while alpha 0.01 is on target",
carried there as `[ASSUMPTION: the known behaviour of the Cauchy combination of
nested, positively dependent cells at a large alpha, accurate in the tail and
mildly liberal at 0.05] -- not tested here`.

**It is now tested. The assumption is CONFIRMED and quantified.**

## The analytic anchors, fixed before any simulation

`tan((0.5 - p) * pi)` is standard Cauchy for uniform `p`, and the Cauchy is
stable with index 1, so a weighted mean of INDEPENDENT standard Cauchys is
again exactly standard Cauchy. Therefore:

- at rho = 0 the pooled p is EXACTLY uniform, not asymptotically;
- at rho = 1 every `p_j` is identical, the statistic is one standard Cauchy,
  so it is exactly uniform again.

Both endpoints exact => any dependence distortion is NON-MONOTONE in rho, with
an interior peak. That is a falsifiable prediction about shape, and it is what
licenses reading the middle of the grid as real rather than as drift.

## Control 1 (negative control): SAIGE's own CCT on independent uniforms

`SAIGE:::get_CCT_pvalue`, 200,000 reps per cell, MC SE 0.00049 at alpha 0.05.
Equicorrelated two-sided z, `p = 2*pnorm(-|z|)`.

| k | rho 0 | 0.3 | 0.5 | 0.7 | 0.9 | rho 1 |
|---|---|---|---|---|---|---|
| 3  | 0.0500 | 0.0517 | 0.0528 | 0.0541 | 0.0517 | 0.0500 |
| 6  | 0.0503 | 0.0528 | 0.0568 | 0.0570 | 0.0527 | 0.0501 |
| 9  | 0.0494 | 0.0553 | 0.0588 | 0.0593 | 0.0538 | 0.0494 |
| 18 | 0.0494 | 0.0576 | 0.0613 | 0.0607 | 0.0537 | 0.0506 |

Both anchors hold within MC error, so the implementation is exact where theory
says it must be -- the negative control passes and the grid can be believed.
The peak is at rho ~ 0.5-0.7 and grows with k, as predicted. At alpha 0.01 the
whole grid stays within +0.0022 of nominal (max 0.0122 at k=18, rho=0.5): the
tail is accurate, the shoulder is not.

At k = 192 (the flex192 weight sets), 40,000 reps, inputs exact by construction:

| rho | cells at .05 | POOLED at .05 | cells at .01 | POOLED at .01 |
|---|---|---|---|---|
| 0.00 | 0.0500 | 0.0517 | 0.0100 | 0.0097 |
| 0.30 | 0.0501 | 0.0698 | 0.0101 | 0.0138 |
| 0.50 | 0.0498 | 0.0694 | 0.0099 | 0.0148 |
| 0.70 | 0.0501 | 0.0657 | 0.0101 | 0.0132 |
| 0.90 | 0.0508 | 0.0553 | 0.0101 | 0.0112 |
| 0.99 | 0.0505 | 0.0507 | 0.0108 | 0.0110 |

## Control 2 (real data): the reconciliation the calibration never did

A NULL flex192 run at N = 50,000 -- `y_null` on sim50k, 210 genes x 192 sets,
40,196 set rows -- left on disk by the same day's timing reps, so this cost no
compute. Rates averaged WITHIN gene and the SE taken BETWEEN genes, because
192 sets of one gene are not 192 observations.

- **The pooled row IS the CCT of its own set rows**: max relative error 4.6e-06
  over all 210 genes, which is the printed precision of the file (7 significant
  digits) and not a discrepancy. Computable two ways, and it reconciles.
- **Measured correlation among the sets**: mean 0.705, median 0.722, quartiles
  0.542 / 0.900, range 0.240-1.000. That is squarely the simulation's
  peak-distortion band, which is why this arm is the worst case rather than a
  typical one.
- **Size**: set level 0.0522 (SE 0.0103) at alpha 0.05 and 0.0132 (SE 0.0061)
  at 0.01; POOLED 0.0714 (SE 0.0178) and 0.0190 (SE 0.0095), n = 210 genes.

Simulation at k = 192, rho = 0.7 gives 0.0501 -> 0.0657 and at rho = 0.5
0.0498 -> 0.0694, against the measured 0.0522 -> 0.0714. They agree.

## What this settles, and what it does not

**Settled: it is the combination, not the inputs.** An earlier reading of this
session -- that the calibration doc's cell-level ECDFs (0.055-0.059) showed the
inputs already liberal, so CCT might be exonerated -- is REFUTED on this arm:
the inputs here are at 0.0522, nominal within noise, and the pooled row is at
0.0714. The input-level offset belongs to a different fixture and different
(sub-5 MAC) strata and does not transfer. Where BOTH are present -- `wgridr`
and `wplain`, k = 9 cells, cells at 0.056-0.057 -- the two contributions add,
and the k = 9 row of the grid above puts CCT's share at about +0.009.

**Not settled by the real arm alone.** With n = 210 genes the pooled 95% CI is
[0.0405, 0.1151], which contains 0.05: one null run cannot reject nominal. The
precision is the simulation's (MC SE 0.0005) and the real run's job is to
supply the correlation and confirm the direction and size. Neither arm carries
the claim by itself.

## The practical consequence

`p_FlexRV` IS the pooled row, so at alpha 0.05 its size is about 0.07 on this
design -- roughly 40% more rejections than nominal. **At discovery thresholds
it is fine**: the alpha 0.01 column is already within +0.002-0.005, the
distortion collapses toward the tail, and that tail is exactly where Liu & Xie's
dependence-robustness result applies. So: do not quote the pooled row's
nominal size at a large alpha; a gene-level threshold (0.05/210 = 2.4e-4, or
genome-wide 2.5e-6) is unaffected. If a nominal large-alpha size is ever
needed, the fix is not a different combiner but a calibrated one (permutation
or a fitted null for the pooled statistic).

Scripts: `cct_control.R`, `cct_decomp.R`, `cct_real.R` in the session
scratchpad (pure R over the installed package plus one file already on disk;
nothing to add to the harness).
