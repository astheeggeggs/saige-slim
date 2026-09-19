# The saddlepoint tail IS accurate at small q (2026-09-19)

`docs/state/research-2026-09-18-skato-speedups.md` open question 2: the landed
saddlepoint tail is unvalidated at small numbers of terms -- the band carrying
56% of this pipeline's Davies time -- because Kuonen's relative error is
`o(n^-3/2)`, i.e. it IMPROVES with n, and the published validation (Chen &
Lumley 2019) spans only n = 305..11,259. **Answered: the worry does not
materialise. The error is flat in n and bounded by 12%.**

Harness `harness/skato_pvalue/tailref/` (README has the method and the run).

## Why a new instrument was needed

Neither available check can reach this regime:
- A null sweep cannot manufacture the cells: 3 unique below 1e-5 in 244,063
  (`handoff-2026-09-19-158-finished-build-guard.md`).
- Davies cannot be the reference below ~1e-6 -- that is the premise of #158
  step 2, and it is re-measured below.

The reference is Ruben's series summed in the UPPER tail with `b <= min(lambda)`,
which makes every coefficient nonnegative, so the survival is a sum of positive
terms with **no 1-F cancellation** -- the property Davies/Imhof/Farebrother lack.
It refuses (returns `ok = FALSE`) whenever the truncation bound is not provably
below `rtol * S`.

## Result: the shipped saddlepoint, on 468 converged points over real spectra

Real per-rho spectra from the guard sweep's own cells, `n_lam` (terms in the
quadratic form) 3-10, rebuilt exactly as `skatoScreenKernel` builds them.
Ratio is shipped / reference; **above 1 is CONSERVATIVE** (p reported too
large: costs power, cannot inflate type-I error).

| target p | n | median ratio | min | max |
|---|---|---|---|---|
| 1e-2 | 88 | 1.0357 | 0.9665 | 1.0503 |
| 1e-3 | 88 | 1.0620 | 0.9795 | 1.0663 |
| 1e-4 | 88 | 1.0756 | 0.9966 | 1.0774 |
| 1e-6 | 88 | 1.0916 | 1.0063 | 1.0921 |
| 1e-8 | 84 | 1.1014 | 1.0134 | 1.1018 |
| 1e-10 | 16 | 1.1075 | 1.0252 | 1.1087 |
| 1e-12 | 16 | 1.1132 | 1.0391 | 1.1203 |

**Flat in the number of terms** -- which is the whole question, since Kuonen's
bound is stated in n:

| n_lam | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| median ratio | 1.0686 | 1.0717 | 1.0655 | 1.0765 | 1.0773 | 1.0751 | 1.0749 | 1.0714 |
| max abs log10 | 0.0467 | 0.0469 | 0.0421 | 0.0469 | 0.0493 | 0.0473 | 0.0465 | 0.0468 |

No degradation at n = 3-5 relative to n = 10. The bias is a function of DEPTH
(3.6% at 1e-2 rising to 11-12% at 1e-12), not of n. This CONFIRMS and extends
the handoff-2026-09-17 figure of "5-12% conservative" -- to 1e-12, and to small
n where it had never been checked.

## The positive control, which is also the strongest evidence

Davies and this reference are independent algorithms. Where Davies is
trustworthy they agree to 6-7 significant figures; where it is not, it fails
exactly as #158 said:

| target p | n | median Davies/ref | max abs log10 |
|---|---|---|---|
| 1e-2 | 88 | 1.00000 | 2.6e-6 |
| 1e-3 | 88 | 1.00000 | 6.4e-6 |
| 1e-4 | 88 | 0.99998 | 1.2e-3 |
| 1e-6 | 88 | 1.00540 | 0.198 (factor 1.58) |
| 1e-8 | 48 | **4.0354** | 1.255 (**factor 18**) |
| 1e-10, 1e-12 | 0 | -- | Davies returns no positive value at all |

The agreement at 1e-2..1e-4 is mutual validation of both methods. The collapse
after is the measurement #158 step 2 was built on, reproduced here from a
different direction: broadly consistent with the handoff's "0.05x-7.5x off in
(1e-8, 1e-6]", with a larger maximum (18x) because these spectra are
worse-conditioned.

**And Davies never says so.** Over the full 840 evaluation points it returned a
**non-positive "probability" on 277 (33%) while reporting `ifault = 0` in every
single case.** Its own fault indicator is not a usable guard. Anything that
consumes `SKAT_davies` output below ~1e-6 must check the SIGN, not the ifault.

## Caveats

- The reference converged on 468 of 840 points (56%). Refusals are the
  ill-conditioned spectra (cond > 5000 skipped by design) and the deepest
  targets -- by construction, not by failure.
- **The 1e-10 and 1e-12 rows rest on 16 points each**, because the reference
  refuses more often the deeper it goes. The 11-12% figure there is the
  weakest number in the table; 1e-2..1e-8 carry 84-88 points each.
- The conservative bias is systematic and grows with depth. It costs power. No
  action proposed here -- it is the price of a method with a PROVED bounded
  relative error, and the alternative in this regime returns negative numbers.

## What this closes and what it does not

CLOSES: research-note open question 2. The tail method is sound at the q that
matter, and `--SKATO_tail_method=saddlepoint` needs no change.

DOES NOT CLOSE: this validates the inner survival function
(`skatoSaddlepointSurvival`). The outer rho integral in its direct piecewise
form (`skatoSaddlepointTail`) is a separate question -- the errors above are
per-survival-call, and how they compose through the outer integral is not
measured here.

SUGGESTED, NOT DONE: this is a desk instrument, not a test. If it should gate,
the cheap form is a handful of pinned (lambda, q, expected) triples from
`scored.csv` in `test-unit-*.R` at the fast tier, plus the two exact closed
forms as the controls -- no Python, no sweep data.
