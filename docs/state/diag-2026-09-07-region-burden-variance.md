# LEDGER #103 / #104 / #105 diagnosed: the SKAT-O path's burden variance is the naive GLM variance (2026-09-07)

Status: diagnosis only. No fix applied, no golden moved, no ledger row edited
(LEDGER.md and next.md carry another session's uncommitted #103 work; the
proposed row text is at the end of this file). Reproduction is
`bash harness/diag/region_burden_variance/grid.sh` (~4 min; every number
below comes from its output under `harness/work/diag-burden/`).

## Summary

1. **#103 stands as written** (trace of `Phi` where the burden information is
   its full sum). Confirmed on every cell: `BETA_Burden == sum(Score)/tr(Phi)`
   to the last digit. Independent of everything below.

2. **#104 and #105 are one mechanism, and it is not "which variance ratio
   should scale a burden statistic".** On binary traits the region path
   rescales `Phi` by `1/min(r, 1)` in `SPA_ER_kernel_related_Phiadj_fast_new`
   (`R/SAIGE_SPATest_Region_Func.R:241-321`), where `r = VarQ / VarQ_2`,
   `VarQ = 1'Phi 1` (the C++ covariance: GRM- and ratio-corrected) and
   `VarQ_2 = S^2 / qchisq(p_burden)`. That `p_burden` is computed by
   `get_newPhi_scaleFactor_traitType` (`:324-357`) from the NAIVE burden
   variance `var1 = sum(mu(1-mu) g.sum^2)` -- no GRM, no variance ratio --
   and below `|z| < 2` it is the plain chi-square p on that variance, so
   `VarQ_2 == var1` exactly (`g.sum` is the residualised weighted sum, so
   `(q.sum - m1)^2 == S^2`). The clamp `min(r, 1)` then means the adjusted
   `Phi` can never have a burden variance BELOW the naive one. On this model
   the GRM and the ratio both reduce the variance (`r` = 0.945-0.953 on all
   12 cells of P6), so on every cell where no per-marker SPA fired
   (11 of 12) the region path's burden variance is `sum(mu(1-mu) g.sum^2)`
   to seven digits, and the whole `Phi` is inflated by `1/r` with it.
   The `--r.corr=1` path is SAIGE's marker machinery on the collapsed vector:
   `ratio x variance` with the ratio and the GRM chosen by the cell's MAC.
   So the cross-mode ratio is the null ratio (0.942146) where the cell MAC
   exceeds 20.5 [#104] and the sparse GRM's own variance reduction
   (0.946-0.954) where it does not [#105]. The 5 % in #105 is NOT the
   SPA-implied variance the row hypothesised; that accounts for 0.36 % on the
   one cell whose pseudo-marker was SPA-corrected.

3. **This is a p-value defect, not only an effect-size one.** `Pvalue_Burden`
   on the SKAT-O path is the no-GRM, no-ratio GLM burden p (0.29692 against
   0.28240 from `--r.corr=1` on GENE1 lof 0.01; with every ratio set to 1.0
   the two agree to 7 digits, 0.2967670). Because `r` rescales all of `Phi`,
   `Pvalue_SKAT` and the SKAT-O `Pvalue` carry the same `1/r` inflation.
   Conservative by 5-6 % in variance on this fixture; by `1/VR x V_null/V_GRM`
   in general, i.e. whatever the mixed model was supposed to buy.

4. **NEW, separate (candidate ledger row): step 2's sparse `Sigma` is one
   triangle, and its PCG stops after 1-2 iterations.** `setSparseSigma_new`
   (`R/readInGLMM.R`) takes `@i/@j/@x` of the `dsTMatrix` that `readMM`
   returns for the `symmetric` .mtx -- 2,700 entries below the diagonal, 0
   above -- and `SAIGEClass` builds `arma::sp_mat` from exactly those
   triplets, so `getPCG1ofSigmaAndGtilde` iterates on a lower-triangular
   matrix. Its stopping rule is ABSOLUTE (`||r||^2 < 0.02`), so an unweighted
   rare variant (`||b||^2 ~ MAC`) stops after ONE iteration and the same
   vector times its beta weight (~23) after two. An R replica of exactly that
   (`pcg_check.R`) reproduces the C++ variances to six digits; the symmetric
   replica does not. Against an exact symmetric solve the C++ sparse variance
   is off by -0.08 % to +0.43 % here, and the same vector gets a variance 0.05-
   0.11 % different depending on its scale. Currently hidden on the burden
   direction by (2); reachable on the sparse marker path (golden P27) and on
   SKAT's off-diagonals (not measured). Small on this fixture only because
   `1/W ~ 11` dominates `tau2 x kinship ~ 0.16`; on a quantitative trait or
   a strongly related cohort the missing triangle is a first-order error.

## The mechanism, line by line (binary traits, `regionTestType != "BURDEN"`)

`region_postprocess` (`R/SAIGE_SPATest_Region.R:350-377`) hands each cell's
`Score` (weighted per-marker scores), `Phi` (`wadjVarSMat[tempPos, tempPos]`,
i.e. the C++ `VarMat` with each row scaled by `sqrt(vr_i)` for THAT marker's
MAC category, times the weights), `g.sum = genoSumMat[, jmr]` (the
residualised weighted genotype sum) and `q.sum = sum(gyVec * w)` to
`get_newPhi_scaleFactor_traitType`, which computes

    m1   = sum(mu * g.sum)
    var1 = sum(mu * (1 - mu) * g.sum^2)          # naive: no GRM, no ratio
    p_burden = pchisq((q.sum - m1)^2 / var1)      # if |q.sum - m1|/sqrt(var1) < 2
             = SPA over (mu, g.sum)               # otherwise, same CGF

and passes it to `SPA_ER_kernel_related_Phiadj_fast_new`, which forms

    D_i     = (Score_i^2 / qchisq(p_i)) / Phi_ii  # per-marker SPA/ER inflation, 1 without SPA
    VarQ    = 1' D^1/2 Phi D^1/2 1
    VarQ_2  = (1'Score)^2 / qchisq(p_burden)      # == var1 when no burden SPA
    r       = min(VarQ / VarQ_2, 1)
    Phi_adj = D^1/2 Phi D^1/2 / r

so `1'Phi_adj 1 = max(VarQ, VarQ_2)`. The C++ marker SPA rescales its
statistic by `sqrt(var1/var2)` before the saddlepoint precisely so that the
ratio is not lost there (`src/SAIGE_test.cpp`, `q = t_Tstat/sqrt(t_var1/t_var2) + m1`);
the R burden SPA has no such rescale, and its no-SPA branch is the naive
chi-square outright. Then `get_SKAT_pvalue` reports `BETA_Burden =
1'Score / tr(Phi_adj)` (#103) and SKAT's rho = 1 p-value is
`pchisq(S^2 / 1'Phi_adj 1)` up to the rho = 0.999 device.

One more thing the dumps show: the collapsed ultra-rare pseudo-marker of the
max_MAF = 0.01 cells has a collapsed MAC of 41-54, so inside the region loop
it takes the COMMON arm (`useSparseGRM = false`, null ratio, category 2)
while the real markers beside it take the sparse arm. `Phi` mixes the two
approximations in one matrix; the rescale then discards both on the burden
direction. (On the quantitative goldens, where nothing is rescaled, this mix
is the whole 1.5-2.3 % cross-mode difference on multi-marker cells -- the
honest residual of the old #104 question, an order of magnitude smaller than
the row claimed.)

## Evidence 1: golden P6 (SKAT-O) against P5 (BURDEN), committed ratio file

`naive = sum(mu(1-mu) g.sum^2)`; `B_raw = 1'Phi 1` as C++ delivered it;
`B_adj = 1'Phi_adj 1`; `var_B` = the `--r.corr=1` variance implied by its
`BETA` (`S / BETA`, same `S`); `BETA_fix(R) = S / B_adj` is #103 corrected
offline. The instrumented runs are byte-identical to goldens P5 and P6.

| cell | m | MAC | z_naive | r | B_adj/naive | B_raw/naive | var_B/naive | BETA_fix(R)/BETA(B) | p_Burden R | p_Burden B | p naive |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GENE1 lof 0.001 | 1 | 5 | 2.214 | 0.99999 | 0.9963953 | 0.950190 | 0.951144 | 0.9545848 | 0.026528 | 0.026527 | 0.026802 |
| GENE1 lof 0.01 | 4 | 98 | -1.043 | 0.95287 | 1.0000000 | 0.952874 | 0.942148 | 0.9421483 | 0.296920 | 0.282397 | 0.296767 |
| GENE1 missense;lof 0.001 | 1 | 10 | 1.901 | 0.94533 | 1.0000000 | 0.945330 | 0.946330 | 0.9463297 | 0.057352 | 0.050728 | 0.057352 |
| GENE1 missense;lof 0.01 | 5 | 119 | -0.913 | 0.94722 | 1.0000000 | 0.947217 | 0.942144 | 0.9421445 | 0.361574 | 0.347023 | 0.361363 |
| GENE1 missense;lof;synonymous 0.001 | 1 | 10 | 1.901 | 0.94533 | 1.0000000 | 0.945330 | 0.946330 | 0.9463297 | 0.057352 | 0.050728 | 0.057352 |
| GENE1 missense;lof;synonymous 0.01 | 5 | 119 | -0.913 | 0.94722 | 1.0000000 | 0.947217 | 0.942144 | 0.9421445 | 0.361574 | 0.347023 | 0.361363 |
| GENE2 lof 0.001 | 1 | 5 | -0.870 | 0.95348 | 1.0000000 | 0.953484 | 0.954174 | 0.9541738 | 0.384408 | 0.373226 | 0.384408 |
| GENE2 lof 0.01 | 2 | 22 | -0.089 | 0.95331 | 1.0000000 | 0.953309 | 0.942147 | 0.9421466 | 0.927633 | 0.926719 | 0.928865 |
| GENE2 missense;lof 0.001 | 1 | 15 | -0.023 | 0.95174 | 1.0000000 | 0.951731 | 0.952224 | 0.9522241 | 0.981413 | 0.980953 | 0.981413 |
| GENE2 missense;lof 0.01 | 4 | 77 | 0.757 | 0.94923 | 1.0000000 | 0.949226 | 0.942147 | 0.9421469 | 0.449325 | 0.435598 | 0.449194 |
| GENE2 missense;lof;synonymous 0.001 | 1 | 15 | -0.023 | 0.95174 | 1.0000000 | 0.951731 | 0.952224 | 0.9522241 | 0.981413 | 0.980953 | 0.981413 |
| GENE2 missense;lof;synonymous 0.01 | 4 | 77 | 0.757 | 0.94923 | 1.0000000 | 0.949226 | 0.942147 | 0.9421469 | 0.449325 | 0.435598 | 0.449194 |

Read across: `r < 1` everywhere; `B_adj/naive = 1.0000000` on the 11 cells
with `|z_naive| < 2` (and `B_adj == VarQ_2` on all 12); the cross-mode ratio
is `var_B/naive` on every row, which is the null ratio on the four
multi-marker cells with MAC > 20.5 and the sparse-GRM reduction (times a
0.05-0.11 % PCG scale effect, item 4) on the 1 x 1 cells. `p_Burden R` is the
naive p on 11 cells up to SKAT's rho = 0.999 device (exactly equal on the
1 x 1 cells, where SKAT sees a scalar). The remaining cell (GENE1 lof 0.001)
is the only one whose marker was SPA/ER-corrected (`Is.SPA` true in the
`.singleAssoc` file): there `D = 1.0486`, the R burden SPA and the C++ ER
imply the same variance to 1e-5 (`r = 0.99999`), and `B_adj = 0.9964 x naive`.

With every ratio set to 1.0 (pair R1/B1) the multi-marker ratios are
1.0000007, 1.0000030, 1.0000006, 1.0000030 and the 1 x 1 ratios 0.9545971,
0.9463442, 0.9541888, 0.9522392 -- the numbers in the #104 and #105 rows,
reproduced by the offline #103 correction. That reconciliation "worked" for
#103 because BOTH sides were the naive variance, not because either is right.

## Evidence 2: the MAC-cutoff and ratio-file grid (predictions written before the runs)

Seven (SKAT-O, BURDEN) pairs of P6's/P5's command; `grid.sh` runs them all.

| pair | ratio file (sparse / null) | category boundary | predicted | measured (multi-marker cells; 1 x 1 cells) |
|---|---|---|---|---|
| R0/B0 | committed 0.99998 / 0.9421 | 20.5 | ratio = VR_null; sparse-GRM reduction | 0.942146-0.942148; 0.946-0.955 |
| R1/B1 | 1.0 / 1.0 | 20.5 | 1.0 exactly; sparse reduction | 1.0000007-1.0000030; 0.946-0.955 |
| R4/B4 | committed | **200** (pseudo-variant on the sparse arm, sparse ratio) | NOT 1: region discarded the GRM | 0.9512-0.9561; 0.946-0.955 |
| R3/B3 | **0.5** / 1.0 | 20.5 | region output unchanged (`r` = 0.48-0.65, rescaled away) | `Pvalue_Burden` within 2e-5 of R0; BURDEN 1 x 1 p-values move 5x |
| R6/B6 | 1.0 / **0.5** | 20.5 | region unchanged; BURDEN multi-marker halves | `Pvalue_Burden` within 6e-6 of R0; ratio 0.5000000-0.5000003 |
| R2/B2 | **2.0** / 1.0 | 20.5 | `r > 1`, clamp defused, region `Phi` = 1.9 x naive | `r = 1` on all 12; `Pvalue_Burden` 0.297 -> 0.411 |
| R5/B5 | 2.0 / 1.0 | 200 (both arms: same GRM, same ratio, no clamp) | agreement to the PCG residual | 0.99990-1.00215 on all 12 cells |

R4/B4 is the experiment the owner asked for, and it settles the question the
old #104 row posed: with the pseudo-variant moved onto the SAME arm as the
region markers (sparse GRM, sparse ratio) the modes still differ by 4.4-4.9 %,
which is `V_sparse/V_null`, the GRM's variance reduction that the rescale
threw away. R3/R6 show the region output is INVARIANT to a ratio of 0.5 on
either row type (a scale of `Phi` cancels in `VarQ/VarQ_2`), i.e. the
variance-ratio file does nothing to the SKAT-O burden variance unless it
pushes `r` above 1. R2 shows the clamp is one-sided: a ratio can only ever
inflate. R5/B5 is the reconciliation: defuse the clamp and put both arms on
one GRM, and two entirely separate code paths agree to 0.1-0.2 % -- the
residual being item 4's scale-dependent PCG truncation.

## Evidence 3: quantitative negative control (goldens P37 / P36)

`get_newPhi_scaleFactor_traitType` runs only `if (traitType == "binary")`,
and the quantitative model has no sparse GRM. Prediction: the 1 x 1 cells
agree exactly and the multi-marker cells differ only by the per-row ratio
category mix (rare rows category 1 = 1.259, the MAC 41-54 pseudo-marker row
and the whole `--r.corr=1` pseudo-variant category 2 = 1.230). Measured
`BETA_fix(R)/BETA(B)`: 0.9999978-0.9999993 on the 1 x 1 cells with identical
p-values; 0.9769-0.9848 on the multi-marker cells. The discrepancy is absent
where the rescale is absent.

## Evidence 4: the one-triangle Sigma and the absolute PCG tolerance

`pcg_check.R`, on the four distinct 1 x 1 cells of R0 (unweighted residualised
pseudo-marker `gu`; `w` its beta weight; variances before the ratio):

| cell (MAC) | C++ region (PCG on `gu`) | C++ BURDEN (PCG on `w gu`) | exact, symmetric | exact, one-triangle | replica: one-triangle, 1 it / 2 it | replica: symmetric, 1 it / 2 it |
|---|---|---|---|---|---|---|
| GENE1 lof (5) | 0.385283 | 0.385670 | 0.385605 | 0.385281 | 0.385283 / 0.385669 | 0.385225 / 0.385604 |
| GENE1 mis;lof (10) | 0.811823 | 0.812682 | 0.808515 | 0.811814 | 0.811823 / 0.812682 | 0.807607 / 0.808515 |
| GENE2 lof (5) | 0.433966 | 0.434280 | 0.432791 | 0.433968 | 0.433966 / 0.434280 | 0.432386 / 0.432790 |
| GENE2 mis;lof (15) | 1.462261 | 1.463018 | 1.456021 | 1.462329 | 1.462261 / 1.463018 | 1.454664 / 1.456020 |

The one-triangle replica reproduces both C++ columns to six digits, with the
unweighted vector stopping after one iteration and the weighted one after
two; the symmetric replica reproduces neither. `Matrix` 1.7.6 keeps a
`dsTMatrix` symmetric through `[idx, idx]`, so its `@i/@j` are one triangle.
[ASSUMPTION: an older `Matrix` returned a general `dgTMatrix` from that
subset, which would have carried both triangles and made this upstream code
correct when written; not checked.] The PCG stopping rule
`sumr2 > tolPCG` with `tolPCG = 0.02` on `||r||^2` is absolute, so the
number of iterations is a function of `||b||^2`, i.e. of the MAC and of
whether the caller weighted the vector first.

## What this means for the ledger rows (proposed text; NOT applied)

- **#103**: unchanged. Its route (b) reconciliation should be re-worded: with
  every ratio 1.0 both modes reduce to the naive variance, so the agreement
  at 1.000000x is a statement about #103 only, not about either path's
  variance.
- **#104**, WHERE -> `R/SAIGE_SPATest_Region_Func.R:324-357`
  (`get_newPhi_scaleFactor_traitType`, `var1` and `pval_noadj`) and `:296-305`
  (`r = min(VarQ/VarQ_2, 1)`, `Phi_ccadj = G2_adj_n / r`). WHAT -> on binary
  traits the SKAT-O path's burden variance is `max(1'Phi 1, S^2/qchisq(p))`
  with `p` computed on the naive variance, so wherever the GRM and the ratio
  reduce the variance (every cell of P6) the region path reports the no-GRM,
  no-ratio burden test: `Pvalue_Burden` conservative by `1/r`, `Pvalue_SKAT`
  and `Pvalue` inflated by the same rescale, `BETA_Burden` (after #103)
  deflated by `r`. The `--r.corr=1` path is unaffected. STATUS -> OPEN,
  diagnosed and reconciled two ways (Evidence 1-3). A characterization
  test can pin `1'Phi_adj 1 == sum(mu(1-mu) g.sum^2)` on P6's 11 no-SPA
  cells through the same three-function wrap `drv.R` uses; the quantitative
  pair is its negative control.
- **#105**: same mechanism as #104, one row or a cross-reference. The
  SPA-implied-variance hypothesis is refuted as the cause (0.36 % on one
  cell); the 0.05-0.11 % residual is the new row below.
- **NEW row**: one-triangle sparse `Sigma` + absolute PCG tolerance
  (Evidence 4). WHERE `R/readInGLMM.R` `setSparseSigma_new` (`rbind(sparseGRM@i,
  sparseGRM@j)`) and `src/SAIGE_test.cpp` `getPCG1ofSigmaAndGtilde`
  (`while (sumr2 > tolPCG ...)`). Reaches every sparse-GRM step-2 golden
  (P3/P6/P8/P10/P19/P27 and their PGEN mirrors, plus the multi-weight and
  conditional group goldens); step 1's sparse handling is R-side
  (`gen_sp_v2`) and was not examined.

## Fix sketch (not applied; predicted golden movement)

- For #104/#105: make the R burden SPA see the C++ variance, the way the
  marker SPA does -- scale `(q.sum - m1)` by `sqrt(var1 / VarQ_D)` (with
  `VarQ_D = 1'D^1/2 Phi D^1/2 1`) before the `|z| < 2` test and the
  saddlepoint, so that below the cutoff `VarQ_2 == VarQ` and `r == 1`, and
  above it the rescale carries only the imbalance correction. Moves EVERY
  binary SKAT-O group golden on every cell (all three p-value columns and
  both burden columns; `r < 1` everywhere here), and no BURDEN, marker,
  quantitative, step-0 or step-1 output. Prediction from Evidence 1: after the
  fix `Pvalue_Burden` on GENE1 lof 0.01 moves from 0.29692 toward P5's 0.28240
  (not onto it: the pseudo-marker row keeps its category-2 null ratio inside
  `Phi` while P5 applies category 2 to the whole vector; the quantitative
  pair shows that residual at 1.5-2.3 %).
- For the new row: `as(sparseGRM, "generalMatrix")` (or the transposed
  triplets) before `@i/@j` in `setSparseSigma_new`, and a relative PCG
  tolerance. Moves every sparse-GRM step-2 golden by up to ~0.5 % in variance
  on this fixture. Gate it as always: prove the symmetric `Sigma` against
  `Matrix::solve` on the P27 markers first (the reconciliation in
  `pcg_check.R` is the template), then re-baseline.

## Verification of the proposed fix (2026-09-07, later the same day)

Both candidate fixes were run without touching the package: `fixes.R` builds
them as replacements of the three namespace functions (fix A: standardise the
burden statistic by `1'Phi 1` before the `|z| < 2` test and the saddlepoint,
the C++ marker SPA's device; fix B: `r <- 1`, i.e. the published method; #103:
`sum(Phi)`), `drv.R` installs them per `SAIGE_FIX`, `fixruns.sh` runs them.
Summary: fix A is the one to take. It reconciles with BURDEN mode, it is
identical to B under the null on this fixture, and under a planted effect its
`SE_Burden` matches the estimator's spread where B's does not.

**1. Reconciliation** (`fix_check.R`). With the category boundary at 200, so
both modes test on the sparse arm with the sparse ratio, SKAT-O under A or B
against BURDEN mode on all 12 cells: `BETA_Burden` ratio 0.99990-1.00215,
`Pvalue_Burden` ratio 0.9976-1.0003 -- the PCG residual of Evidence 4. At the
default boundary the 1 x 1 cells sit at 1.0005-1.0011 and the multi-marker
cells at 0.9887-0.9946, the per-row ratio-category mix (the quantitative pair
shows 0.977-0.985 for the same reason). Today's code with #103 alone: 0.942-
0.955, p ratios up to 1.13. A and B agree to 1.6e-7 on every cell: on this
fixture the burden SPA never implies more variance than the per-marker
adjustments already do, so once standardised the `r` step is inert.

**2. Null calibration** (`sim.sh null`, `sim_analyze.R`). A permuted phenotype
cannot reach this defect: permutation destroys the family structure, the refit
puts tau2 near 0 and `r` near 1. So the null is parametric from the fitted
GLMM: `y ~ Bernoulli(plogis(offset + intercept + b))`, `b ~ N(0, tau2 K_sparse)`
with the model's own tau2 = 0.329, step 1 refitted per replicate in the
`binary_sparseGRM_vr` configuration, the instrumented SKAT-O run, and the
three variants applied offline to each replicate's dumps (the offline path is
asserted to reproduce the dumped output bit for bit). 100 replicates x 8
distinct cells; the rescale fired on 95.8 % of cells. Rejection rate at 0.05 /
mean of the chi-square-inverted p, today+#103 / A / B:

| column | multi-marker cells (n = 400) | single-pseudo-variant cells (n = 400) |
|---|---|---|
| `Pvalue_Burden` | 0.025 / 0.025 / 0.025; 0.857 / 0.880 / 0.882 | 0.0375 / 0.040 / 0.040; 1.009 / 1.032 / 1.032 |
| `Pvalue_SKAT` | 0.065 / 0.070 / 0.070; 0.998 / 1.041 / 1.044 | same as Burden (a 1 x 1 cell) |
| `Pvalue` (SKAT-O) | 0.0325 / 0.040 / 0.040; 0.934 / 0.972 / 0.976 | same |

Binomial SE of a rate 0.011, of a mean chi-square 0.071, both optimistic
(cells within a replicate are correlated; a KS test is not usable, the MAC-5
exact-test p-values are discrete). Nothing is anti-conservative beyond noise;
the fixes are nearer nominal on the burden and SKAT-O columns; the difference
between variants is inside the noise, as the power analysis said it must be
for `r ~ 0.95`.

**The variance identity, where the power is** (`rep_stats.R`, pure R on each
refit: `S = g~'(y - mu)` for the planted cell's burden variable, against its
candidate null variances). At the model's tau2 (n = 100): `Var_emp(S) /
mean(V)` is 0.853 for the naive variance, 0.885 for the sparse-GRM one, SE
0.14 -- inconclusive. At tau2 = 4 (`sim.sh nullx 400 0 s1 3 4.0`, step 1
only, n = 400, SE 0.071): naive 0.691, sparse 0.814, BURDEN-mode (null ratio x
naive) 0.821. The naive variance is 4.4 SE too large, the sparse one 2.6 SE.
Both are too large because the refit recovers tau2 = 0.71 of 4.0 on average
(0.19 of 0.33 at the model's tau2): the downward bias of the penalised
quasi-likelihood fit for binary traits, outside these rows and recorded here
because it bounds what any step-2 fix can achieve. The mixed-model variance
is itself conservative under strong relatedness; the rescale to naive makes it
more so.

**3. Planted effect** (`sim.sh planted 100 0.05 full`, `planted_mle.R`). A
burden coefficient of 0.05 on SAIGE's own weighted burden variable for GENE1
missense;lof 0.01 (per-allele log-odds ~ 1), 100 replicates, the burden SPA
firing in 85 % of them. Mean `BETA_Burden` / 0.05, sd, mean `SE_Burden`,
sd / mean SE, 95 % coverage:

| variant | mean/beta_b | sd | mean SE | sd/SE | coverage |
|---|---|---|---|---|---|
| today + #103 | 1.039 | 0.0160 | 0.0167 | 0.957 | 0.97 |
| fix A | 1.066 | 0.0165 | 0.0169 | 0.977 | 0.97 |
| fix B | 1.108 | 0.0187 | 0.0172 | 1.089 | 0.94 |
| BURDEN mode, `S/V` with the score SE | 1.157 | 0.0206 | 0.0175 | 1.182 | 0.88 |
| exact sparse solve, `S/V` | 1.160 | 0.0207 | 0.0175 | 1.185 | 0.88 |
| MLE, true `b` as offset | 0.970 (SE 0.029) | 0.0144 | | | |
| MLE, `b` ignored (marginal) | 0.917 | 0.0145 | | | |

The score estimator overshoots the coefficient by 16 % at this effect size --
`E[S] = V beta` holds to first order only, and maximum likelihood with the
true random effects recovers 0.97 -- so the arm does not test which variance
gives an unbiased `BETA_Burden`; the region variants sit lower because their
variance carries the SPA inflation, two errors partly cancelling. What the arm
does test is `SE_Burden`: A's SE matches the estimator's spread (0.98), B's
understates it by 9 %, the plain score SE by 18 %, with coverage 0.97 / 0.94 /
0.88. This is the one place A and B differ, and it favours A.

**4. No-op checks.** Under A and under B the quantitative golden P37 is
byte-identical (the adjustment is binary-only). On the sparse-ratio-2.0 file,
where the clamp was already inert, A leaves every p-value within 5e-7 of the
unpatched run.

## State of the tree at handoff (2026-09-07, end of session)

No commit and no push this session; the branch's end commit and CI state are
those of `docs/state/handoff-2026-09-07-ledger102.md` (`96f16942` code,
`99541489` comparison, all four CI runs green). Untracked, from this session:
this file and `harness/diag/region_burden_variance/` (scripts, all lintr-clean
against `.lintr`; outputs under the gitignored `harness/work/diag-burden/`,
about 3 h of runs, reproducible with `grid.sh`, `fixruns.sh` and the three
`sim.sh` arms). Uncommitted and NOT this session's: LEDGER rows #103-#105 and
`tests/testthat/test-step2-burden-trace.R` (the #103 characterization), the
next.md sections on the port question, the region path and FlexRV, and
`docs/design/flexrv-plan.md`, and the companion measurement
`docs/state/diag-2026-09-07-region-vs-burden-exact.md` (which reads this
session's R0 dumps and reproduces Evidence 4 independently); rows #104 and #105 were rewritten and #106 added
on top of that working copy at the owner's request, and the pointer in
next.md was moved to this file. Nothing in `R/`, `src/` or the goldens moved.

The decision left for the owner is the one at the end of "Verification":
apply fix A (four lines in `get_newPhi_scaleFactor_traitType`, reference
implementation `fixes.R` `phi_a`) as one commit for #104 + #105 with the
characterization the ledger row describes, then #106 as its own row and
commit; both predicted to move goldens as the rows say.

## Reproduction

`harness/diag/region_burden_variance/`:
- `drv.R` -- wraps `get_newPhi_scaleFactor_traitType`,
  `SPA_ER_kernel_related_Phiadj_fast_new` and `get_SKAT_pvalue` in the
  installed namespace (inputs and outputs to `.rds`), then sources the
  unmodified CLI script; outputs stay byte-identical to the goldens.
- `run.sh` / `run_quant.sh` -- one instrumented run of P6/P5 (binary) or
  P37/P36 (quantitative); `grid.sh` -- all sixteen runs, the golden `cmp`s
  and the four analyses.
- `analyze.R` -- per-cell identities and cross-mode ratios, with the
  independent `Matrix::solve` sparse variance; `summary.R`, `numdiff.R`,
  `pcg_check.R`, `quant_check.R`.
- `fixes.R` -- the candidate fixes over the installed namespace; `drv.R`
  installs them per `SAIGE_FIX` (`103`, `A`, `B`, combined); `fixruns.sh` +
  `fix_check.R` -- verification 1 and 4.
- `sim.sh <label> <n> <beta_b> <full|s1> [parallel] [tau2]` over
  `sim_prep.R`, `sim_pheno.R`, `sim_rep.sh`, `rep_stats.R`; `sim_analyze.R`
  and `planted_mle.R` -- verification 2 and 3 (`null 100 0 full`, `planted
  100 0.05 full`, `nullx 400 0 s1 3 4.0`; about 45 min per full arm at 4
  parallel on a loaded box, 15 min for the step-1-only arm). Seeds are per
  replicate id, so an interrupted arm resumes by re-running the missing ids
  with `sim_rep.sh` and `SIM_PREP` set. macOS has no `setsid`.
All R files are lintr-clean against `.lintr`.
