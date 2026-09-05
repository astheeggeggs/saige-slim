# The conditional saddlepoint uses the residualised genotype (LEDGER #84)

Decided 2026-09-05 (owner: residualise); implemented the same day in
`src/SAIGE_test.cpp`, the conditional block of `SAIGEClass::getMarkerPval`.

## What was wrong

For a tested marker with covariate-adjusted genotype `g~` (n x 1) and
conditioning markers `G2~` (n x k), SAIGE's conditional score and variance are

    T_c   = T  - (g~' P G2~) (G2~' P G2~)^-1 T2          T  = g~' res,  T2 = G2~' res
    var_c = var - (g~' P G2~) (G2~' P G2~)^-1 (G2~' P g~)

with `P` the GRM-adjusted projection and the variance-ratio scalings folded
into the stored factors (`m_P2Mat_cond`, `m_VarInvMat_cond`, `m_Tstat_cond`).
When `T_c^2 / var_c` clears the SPA cutoff, the p-value is saddlepointed. The
saddlepoint needs the cumulant generating function of the statistic being
tested. The block built `q_c = T_c / sqrt(var_c / var2) + m1` and handed
`SPA_fast` the MARGINAL marker's CGF pieces, `g~` itself with `m1 = sum(mu g~)`
and `var2 = sum(mu(1-mu) g~^2)`: the case-control-imbalance correction of a
different statistic. The error grows with the LD between the tested and the
conditioning markers (LEDGER #84's provenance: -0.007 decades at r^2 = 0.378,
-0.445 at 0.854, -1.643 at 0.982, at z_c = 5, anti-conservative).

## The identity

Write `c = (G2~' P G2~)^-1 (G2~' P g~)` (k x 1). Then

    T_c = g~' res - c' G2~' res = (g~ - G2~ c)' res = g~_c' res,     g~_c := g~ - G2~ c

so the conditional score IS the score of the residualised genotype `g~_c`, and
the statistic whose CGF the saddlepoint needs is `g~_c' y`. For a binary trait
that CGF is

    K_c(t) = sum_i [ log(1 - mu_i + mu_i exp(t g~_c,i)) ] - t sum_i mu_i g~_c,i

with mean `m1_c = sum(mu g~_c)` and variance `var2_c = sum(mu(1-mu) g~_c^2)`.
The block now forms `g~_c`, and calls the saddlepoint with

    q_c    = T_c / sqrt(var_c / var2_c) + m1_c
    qinv_c = 2 m1_c - q_c
    SPA(mu, g~_c, q_c, qinv_c, ...)

exactly parallel to the marginal block with `g~` replaced by `g~_c`. The
`sqrt(var_c / var2_c)` factor is the same variance-ratio rescaling the
marginal block applies (`var1 / var2`): the saddlepoint is taken on the
statistic rescaled to the CGF's own variance.

`c` is `(t_G1tilde_P_G2tilde * m_VarInvMat_cond)'`, already computed for
`T_c`; the only new storage is `G2~` itself (`m_G2tilde_cond`, n x k, one
column per conditioning marker, filled in `assign_conditionMarkers_factors`
beside the scaled `P1Mat` row), and the only new arithmetic is one n x k
product per marker that clears the cutoff.

## Full SPA, not SPA_fast

`SPA_fast` treats the samples with genotype 0 by a normal approximation of
their (small, covariate-adjustment-only) contribution and the carriers
exactly. `g~_c` is dense: `G2~ c` touches every sample, and for a marker in LD
with the conditioning set the carriers of the conditioning markers carry most
of the residual. The zero-genotype split of the TESTED marker is therefore not
the right partition, and the block uses the full `SPA` over all n samples. It
runs only for markers whose conditional statistic clears the cutoff.

## The `--is_noadjCov=TRUE` caveat

At the CLI default the marginal `T` is the UNADJUSTED score `g' res` with the
`null_noXadj` variance ratio, while the conditioning terms are built from the
ADJUSTED `g~` (LEDGER #82, settled: the conditional columns are a mixture). Then
`T_c` is not exactly `g~_c' res`, and `q_c` is that mixture statistic
saddlepointed against `g~_c`'s CGF. Under `--is_noadjCov=FALSE` the identity is
exact, which is why the gate runs there.

## Gate

`tests/testthat/test-step2-cond-highld.R`, block "LEDGER #84 FIXED": one
statistic computed two ways. With ONE conditioning marker, `c` is a scalar and
is read off the output as `(T - T_c) / T2`; the test rebuilds `g~` and `g~2`
from the fixture's own genotypes with the model's `obj.noK` (the same
`g - XXVX_inv (XV g)` as C++), checks `g~_c' res = T_c` to printing precision
(`Tstat` carries 6-7 significant digits, so 1e-5), and requires
`p.value_c` to equal `SPA_pval(mu, g~_c, q_c, qinv_c, ...)` to 1e-6 on the
marker that clears the gate under real conditioning (`sig` against `sigpart`,
var/var_c > 1.3). Red at unpatched HEAD by the marginal-vs-residualised gap,
green after. The Phase 2 characterization block (the conditional SPA runs and
moves the p-value by more than a decade) is kept.

## What moved, and the cost (measured 2026-09-05)

- Wall time of golden P12's command (70 markers, `--condition=rs13,rs79`):
  1.90 / 2.00 / 2.01 s before, 1.91 / 1.87 / 1.92 s after. The extra work is
  one n x k product per marker that clears the conditional cutoff.
- P12 (and its pgen twin P17): `p.value_c` moved on the 4 of 70 rows that clear
  the gate and nowhere else; `Tstat_c`, `var_c`, `BETA_c` and every marginal
  column are byte-identical. rs23 +0.0075 decades (var/var_c 1.0073), rs38
  +0.0007 (1.0015), rs56 +0.064, rs44 -0.53. The last two are the MAC 1 and 2
  exact-test markers at var/var_c 1.0000: their movement is the fast-vs-full
  saddlepoint difference on a near-singleton, not LD.
- The single-variant files of P3, P10, P15, P24, P34 moved on their 1-3 rows
  that clear the gate. Every SKAT-O region file is unchanged (its conditional
  burden p is saddlepointed on the R side). NOT predicted, and caught by the
  standard tier: the BURDEN conditional goldens P25 and P26, whose C++ burden
  score goes through this same block -- one cell (GENE1, lof, 0.001) moved from
  `Pvalue_Burden_c` 2.639e-02 to 2.398e-02, with its `seBeta_Burden_c` and the
  `Cauchy` row; plink and pgen identical.
- The high-LD fixture's `sig` marker (var/var_c 1.51, `stat_c` 91): 2.8e-18
  before, 1.6e-19 after, -1.25 decades; `p.value.NA_c` is 1.3e-21. So on this
  fixture the MARGINAL CGF was the more conservative of the two. LEDGER #84's
  provenance called the marginal CGF anti-conservative from a different toy;
  the sign depends on the skewness of `g~` against `g~_c`, and neither sign
  is general. What is general is which CGF belongs to the statistic.
