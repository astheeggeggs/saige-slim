# Exact-reference control for the SKAT-O saddlepoint tail (2026-09-19)

Answers the question `docs/state/research-2026-09-18-skato-speedups.md` left
open: is the landed saddlepoint tail (`--SKATO_tail_method=saddlepoint`)
accurate at SMALL numbers of terms, where Kuonen's `o(n^-3/2)` relative-error
bound runs the wrong way and no published validation exists below n = 305?

**A null sweep cannot answer it** (measured: 3 unique cells below 1e-5 in
244,063 -- `docs/state/handoff-2026-09-19-158-finished-build-guard.md`), and
Davies cannot be the reference down there, which is the premise of #158 step 2.

## The reference: Ruben's series summed in the UPPER tail

`ref_ruben.py`. Ruben/Farebrother gives `F(x) = sum_k a_k P(chi2_{n+2k} <= x/b)`
and `sum_k a_k = 1`, so the SURVIVAL is `sum_k a_k P(chi2_{n+2k} > x/b)`.
Taking `b <= min(lambda)` puts every `gamma_j = 1 - b/lambda_j` in `[0,1)`,
hence every `a_k >= 0`: **a sum of nonnegative terms, with no 1-F subtraction
anywhere.** That is the one property Davies, Imhof and Farebrother all lack,
and it is why this can be trusted where they cannot.

Truncation at K omits `R_K = 1 - sum_{k<=K} a_k`, which bounds the omitted
contribution because every chi-square factor is <= 1. A value is returned as
`ok` only when `R_K < rtol * S`. **The reference refuses rather than guesses.**

## Run

```sh
Rscript spectra.R  <nullrep work dir> <outdir>   # real per-rho spectra, q = 4..10
Rscript evaluate.R <outdir>                      # shipped saddlepoint + SKAT Davies at target p
python3 score.py   <outdir>                      # Ruben reference -> scored.csv
```

Spectra are REAL, taken from the guard sweep's own cells and rebuilt exactly as
`skatoScreenKernel` builds them (`R_rho = (1-rho)I + rho 11'`, pivoted chol,
`Get_Lambda`). 20,888 collected; median condition 18.7, 75th pct 50, heavy tail
to 4.6e5 at rho -> 1. The reference is tight on the bulk and refuses on the
ill-conditioned decile -- which is the regime where Kuonen's proved tail rate
applies most directly anyway.

## Controls (both must pass before any number here is believed)

- `exact_equal`: `lambda_j == c` gives `c * chi2_n` exactly. Reference ratio
  1.0000000000 at p = 1e-2, 1e-6, 1e-12. Validates `a_0` and the plumbing, but
  the series TERMINATES (gamma = 0), so it does not exercise the recursion.
- `exact_two_pairs`: `lambda = (a,a,b,b)` is `2a E1 + 2b E2`, whose survival has
  a closed form. Distinct eigenvalues, so the gamma recursion IS exercised.
- The strongest control is in the result below: **Davies and this reference are
  two completely independent algorithms, and they agree to 6-7 significant
  figures wherever Davies is trustworthy.**
