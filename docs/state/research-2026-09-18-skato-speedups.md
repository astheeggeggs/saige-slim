# Research: further accuracy-preserving speed-ups on the SKAT-O path (2026-09-18)

Deep-research run (5 search angles, 18 primary sources fetched, 90 claims
extracted, 25 adversarially verified 3 votes each, 12 confirmed / 13 killed).
Question: beyond the landed #158 pair (modified-Liu outer-integral screen,
Kuonen saddlepoint tail), what accuracy-preserving speed-ups remain on the
SKAT-O path, for N from ~3,000 to ~400,000?

**Bottom line.** Almost all the remaining prize is OUTSIDE the p-value kernel.
The literature confirms the landed design is the best-known point in the
accuracy/cost space -- there is no faster Davies and no drop-in replacement.
What is left inside the kernel is capped by Amdahl at <=1.6x at N=50k and
~1.1x at N=400k. The ranked prize is carrier-space algebra for Phi and the
ultra-rare pass, whose payoff is the only one that GROWS with N.

Epistemic key used throughout: **[MEASURED]** in this repo, **[MEASURED-EXT]**
by a cited paper on its own data, **[DERIVED]** flop count or algebra with no
measurement, **[ASSUMPTION]** gap-fill with its basis named.

---

## 1. The Amdahl ceiling that ranks everything else

After #158 the p-value kernel is roughly **30 s of a ~70-80 s SKAT-O run at
N=50k**, and its share FALLS as N grows, because **Davies never sees N** -- its
cost is set by q, the variants per cell (4-134 on sim50k, 482 on All of Us
chr20), not by sample size.

- C++ buckets, sim50k binary SKAT-O: test 14.1 s, urPass-minus-urTest 13.1 s,
  rareOther 6.8 s, flush 4.4 s, varmat 3.8 s; C++ total 56.7 s; whole-run user
  CPU 207 s; peak RSS 808 MB. [MEASURED,
  `docs/state/handoff-2026-09-17-skato-pvalue-methods.md` section 1]
- Kernel 166 s shipped -> 27-30 s with the screen. [MEASURED, section 8]
- The SKAT-O-specific O(N) C++ work is ~27 s = 13% of the pre-#158 run, and
  "at biobank N the O(N) share grows ~8.5x while Davies does not". [section 1]

So: **making the entire kernel free is worth <=1.6x at N=50k and ~1.1x at
N=400k**, where the O(N) terms are ~75% of the run and the kernel ~10%.

[ASSUMPTION: the ~70-80 s post-#158 whole-run figure is COMPOSED from two
separate PRE-#158 profiles taken on a box under load 4-13. Basis: sections 1
and 8 of the same handoff, which overlap but were not taken in one run. This is
why "re-profile" is ranked second below, not tenth.]

---

## 2. The landed design is the literature's design

No better-known replacement for Davies exists. The three error regimes are now
pinned by primary sources:

- **Davies (1980) and Farebrother (1984) compute the upper tail by subtraction
  from 1** and "break down completely if the extreme right tail probabilities
  are near or beyond machine epsilon" (Chen & Lumley 2019, CSDA 139:75-81, S2;
  same sentence in `survey::pchisqsum`).
  **CORRECTION that matters here, caught by a verifier and confirmed in-repo:
  machine epsilon is the METHOD's floor, not this pipeline's.**
  `SKAT:::Get_PValue.Lambda` calls `davies(Q[i], lambda, acc = 10^(-6))`
  [MEASURED: deparsed from the installed SKAT, 2026-09-18], so the eight
  per-rho Davies calls are untrustworthy from ~1e-6 and return negative or
  exactly-zero values by 1e-8..1e-12 -- ten orders of magnitude above eps, and
  consistent with this branch's own 0.05x-7.5x measurement below 1e-5.
  **Do not quote machine epsilon as the trust bound.**
- **Moment matching is anticonservative in the right tail** "for upper tail
  probabilities less than 10^-5" and "more variable depending on a than the
  other methods" (`survey::pchisqsum`); fastSKAT adds the mechanism --
  chi-square-shaped moment matches have tails too light versus a weighted sum
  of chi-square_1, "often misstating p-values in the vicinity of 10^-6 by an
  order of magnitude". Legitimate ONLY as an above-threshold screen, which is
  exactly the screen's design.
- **Kuonen's saddlepoint is the only method with a proved bounded relative
  error**: Chen & Lumley Appendix A Theorem 1 gives the correct exponential
  tail rate 1/(2*lambda_1) for any A and Sigma; empirically max log10 error
  ratio vs Davies < 0.07 (factor ~1.17) from p=1e-1 to 1e-13 on six simulated
  sequence quadratic forms, n=305..11,259. [MEASURED-EXT] This branch's
  measured 5-12% conservative bias sits inside that.
- Davies' cost is not fixed -- at acc=1e-13 wall time rose 0.01 -> 34.77 s
  (~3,500x) purely by scaling the largest eigenvalue, because the term count
  explodes when a few large eigenvalues dominate. **Negative control matters:
  at SKAT's acc=1e-6 the same experiment is FLAT (0.006..0.050 s), so this is
  not a cost this pipeline pays.**

---

## 3. Ranked shortlist

| # | Lever | Class | Speed-up | N-dependence |
|---|---|---|---|---|
| 1 | Carrier-space Phi + ultra-rare pass (four-term expansion) | exact to ~1e-15, numeric gate | ~1.15x @50k, potentially 2-4x @400k [DERIVED] | **grows with N**; ~nil at 3k |
| 2 | Re-profile post-#158; re-measure the `region_postprocess` double drive | measurement | -- | -- |
| 3 | Small-q validation control for the landed saddlepoint tail | accuracy, not speed | -- | q-dependent |
| 4 | Vendored thread-safe `qfc` + compiled saddlepoint | **bit-exact by construction** | retires the fork pool's double drive | N-independent; value shrinks |
| 5 | Closed-form R_rho^(1/2) rank-2 assembly of the 8 spectra | ~1e-15, moves goldens | 1.45x @q=134, 2.04x @q=482, eigen layer only | q-dependent, not N |
| 6 | SMMAT-E (Fisher 4 df, no rho grid) | **different test** | <=~1.2x on top of #158 | -- |
| 7 | fastSKAT / bigQF | approximation, anticonservative in a named regime | **measured LOSS at these q** | -- |

### 1. Carrier-space Phi and ultra-rare pass -- the highest-value lever

bigQF's *design* points at this even though bigQF itself does not apply: never
form the kernel, require only mult/tmult plus the trace, worked example a
sparse genotype matrix times a low-rank projection.

Transposed to this region path, the exact identity is

    g~_j' Sigma^-1 g~_l  =  g_j'Sigma^-1 g_l
                            - (X'Sigma^-1 g_j)'(X'Sigma^-1 X)^-1 (X'Sigma^-1 g_l)

With Sigma diagonal (no-GRM / dense-fit branch) every term is carrier-space:
`g_j'Sigma^-1 g_l` is a merge over the INTERSECTION of two carrier lists,
O(nnz_j + nnz_l); `X'Sigma^-1 g_j` is O(p*nnz_j). Total O(q^2*nnz + q*p*nnz)
against today's O(q^2 N + q p N). At MAF < 1%, nnz <= 0.02N, i.e. **>=50x fewer
operations** on the varmat / rareOther / flush buckets. [DERIVED -- flop count,
not a measurement; ignores merge overhead and dense common-variant cells.]

Same change the #158 handoff already scoped ("a carrier-space VarMat ... could
remove it and the P1/P2 memory"), with its blockers already known: not
byte-neutral (reassociated sums), cannot reach the sparse-GRM branch where
Sigma^-1 is a per-marker solve (3282 of 9058 tested rows on that recipe), caps
at ~1.15x at N=50k. **It ranks first anyway because of section 1**: at N=400k
the O(N) buckets are ~75% of the run. It also removes the largest memory term
(peak RSS 808 MB on sim50k "is the P1/P2 store"). Precedent that this class
lands byte-identically here: the 2026-09-16 carrier-space read gave BURDEN
1.56x / SKAT-O 1.13x with every golden byte-identical.

### 3. The one genuine accuracy gap found

**The landed saddlepoint tail is unvalidated in exactly the band where it fires
most.** The q=4-10 band costs 138 ms per call and carries **56% of Davies time**
here; q > 100 carries 2% [MEASURED, handoff section 2]. Chen & Lumley supply
both halves of the problem: Davies is slow to reach high accuracy "if the sum
is dominated by a small number of eigenvalues" (what a 4-10 marker cell is),
while Kuonen's relative error is "of order o(n^-3/2), so the approximation
improves with increasing n", validated empirically only over n=305..11,259, in
a paper that explicitly frames itself against prior literature dealing with
n <= 10. **The tail method is fastest and least validated in the same band.**
This is a missing control, not a known defect.

Proposed control (positive and negative before interpreting): build small-q
cells with a KNOWN exact answer -- for q<=10 the weighted sum of chi-square_1
admits the Ruben/Farebrother series (AS 204, `CompQuadForm::farebrother`) or
high-precision Imhof -- and compare the landed saddlepoint against it from
p=1e-2 to 1e-12, alongside the existing adaptive-tolerance Davies reference
(valid only to ~1e-9). **Counter-constraint: Farebrother also subtracts from 1,
so the reference must be built in the LOWER tail or at extended precision.**

### 4. Thread-safe `qfc` -- justified now by a different prize

The blocker is on record: "Davies' routine is file-scope statics and a global
`jmp_buf`" (`harness/bench/NOTES.md`, `docs/state/region-perf-handoff.md`). The
workaround -- a persistent fork cluster -- was chosen on measurement
(fork-per-region 55-90 ms against ~350 ms of work, a net SLOWDOWN at 8 workers;
persistent cluster ~3x on the parallel part) and costs a structural tax:
`region_postprocess` is driven TWICE through its injected `skat_fun`.

The OLD argument for vendoring (a compiled Davies) was measured at ~1.8% of a
run and rejected. The LIVE arguments are different: (a) the compiled saddlepoint
`next.md` already lists as open -- the R saddlepoint costs 83-251 ms per cell
against Davies' 48-61 ms; (b) **with the kernel down 166 -> 30 s, the second
drive is a proportionally larger share of what is left, and has not been
re-measured since #158.** Making `qfc` reentrant is mechanical (file-scope
state -> a struct passed by pointer; `jmp_buf` -> a per-call error return) and
is **bit-exact by construction** -- same arithmetic, same order.

### 5. Per-rho-layer micro-optimisations (both ~1e-15, neither bit-exact)

(a) **Cholesky-free closed form.** With a=sqrt(1-rho), b=(sqrt(1-rho+q*rho)-a)/q,
u=Phi2%*%1, s=1'Phi2*1:

    R^(1/2) Phi2 R^(1/2) = a^2*Phi2 + a*b*(1u' + u1') + b^2*s*11'

collapsing the 8 spectra to a rank-2 assembly. [MEASURED in this session on
SYNTHETIC beta(1,25)-weighted Phi, 1 thread, loaded box: q=20 both <0.05 ms;
q=50 0.15 -> 0.05 ms; q=134 3.20 -> 2.20 ms (1.45x); q=300 15.80 -> 11.00 ms;
q=482 47.80 -> 23.40 ms (2.04x); max relative eigenvalue error 8e-16..3.8e-15.]
Bounding context: SKAT's `Get_Lambda` already uses `eigen(only.values=TRUE)`, so
there is no free win from unwanted eigenvectors, and the eight-spectra layer is
only ~20-25% of a Davies-path cell even at q=134, ~2-5% at q<=50. **This is an
All-of-Us-chr20-scale item (q=482), not a sim50k item.**

(b) **All 8 spectra are rank-2 updates of ONE base eigendecomposition**,
O(q^3 + 8q^2) instead of 8*O(q^3) -- the standard symmetric rank-one eigenvalue
update via the secular equation (Golub 1973 / Bunch-Nielsen-Sorensen 1978 /
Cuppen 1981; LAPACK `dlaed4`), applied twice. [DERIVED, not implemented; LAPACK's
constants strongly favour a plain `dsyevr` below q~300.]

Both change the last digits, so both need the numeric comparator and a golden
re-baseline. **On a branch gated byte-exact that is a real cost -- do them only
alongside a change that already moves goldens.**

### 6. SMMAT-E -- adopt for science, not for speed

The only published way to delete the rho grid entirely: Fisher's method on 4 df
over the burden p and an adjusted SKAT p, justified by approximate asymptotic
independence (Chen et al. 2019, AJHG). Table 4 (TOPMed fibrinogen, n=23,763,
46,859 windows, single thread): SMMAT-B 81, SMMAT-S 91, SMMAT-O 266, SMMAT-E
103 CPU-min -> **2.6x** versus their own rho-grid test; marginal kernel cost 22
vs 185 min, ~8.4x on that portion. [MEASURED-EXT] Implemented in GMMAT (CRAN
1.5.0).

Why it does not rank high: the branch's kernel is already 5.5-6.9x cheaper, so
the remaining headroom collides with section 1 (**<=~1.2x**); it is an
accuracy-CHANGING substitution resting on an asymptotic-independence
approximation; and the paper itself flags that "the asymptotic issue of
single-variant scores also exists for binary traits with highly unbalanced
case-control ratios" -- the exact regime SAIGE's SPA/ER machinery exists for, so
it needs its own null calibration on both traits first. **One attractive
side-effect: it returns the burden p alongside, which is the signed Z a Stouffer
meta-analysis needs and a directionless SKAT-O p cannot give.**

### 7. fastSKAT / bigQF -- do not adopt

Killed on its own authors' crossover numbers (Chen & Lumley 2019, Table 5, s):
SKAT-package example Q0 (s=2000, m=67) SVD 0.03 vs SSVD 0.12 -- **a 4x LOSS**;
Q1 (n=305) a dead heat; Q4 (n=3985) 811.56 vs 36.06, 22x gain. The paper states
it flatly: "unless n is greater than hundreds, there is no reason to use the
leading-eigenvalue approximation as it does not save any time." **Here n is the
number of non-zero eigenvalues (~ variants per cell), NOT sample size** -- their
own within-sample-size control is Q0 (s=2000, loss) vs Q3 (s=2000, n=1063, 5.2x
gain). This branch's largest q is 134 (482 at All of Us chr20), so it is on the
wrong side of the crossover **at every N from 3k to 400k**.

Three further disqualifiers: bigQF implements SKAT only -- no rho envelope, no
Lee outer integral, so it cannot produce a SKAT-O p-value; it has a documented
failure mode ("If the remainder term ... has less than 1 degree of freedom it
will be dropped ... The approximation will then be anticonservative"); and its
sparsity win is in the OPERATOR, not the matrix -- bigQF's own docs note the
matrix "has been centered and is no longer sparse".

---

## 4. Foreclosed -- do not re-open

1. **Loosening Davies' `acc`.** [MEASURED] 174.5 -> 122.3 s going 1e-6 -> 1e-3
   (24-30%), but "up to 80% relative movement on genes with p < 1e-4", and the
   call count RISES as acc loosens because the outer quadrature subdivides more.
2. **Modified Liu everywhere.** Fastest on the table (8.0x/8.2x kernel) and
   wrong where it matters: Liu/Davies ratio 0.02-2.0 in (1e-8,1e-6] and
   0.13-1.0 in (1e-6,1e-4]; for a p-value-level meta-analysis the per-cohort
   bias adds coherently.
3. **Liu for the six NON-REPORTED per-rho values.** Verified here: only two of
   the eight are reported -- `R/SAIGE_SPATest_Region_Func.R:466-474` takes
   `p.val.each` at rho==0 as `Pvalue_SKAT` and at rho==1 as `Pvalue_Burden`
   (`SE_Burden` back-derived from the latter) -- so the other six feed only
   `pmin`, and Liu would save ~6/8 of the per-rho layer on the ~85% of cells
   that screen out, ~5 s of a ~71 s run (~7%). **Not worth it: it changes every
   screened cell's reported Pvalue and destroys the invariant that makes #158
   safe** -- that the per-rho layer is byte-identical, so `Pvalue_SKAT`,
   `Pvalue_Burden`, `BETA_Burden` and `SE_Burden` never move.
4. **Single/mixed precision in the eigen or Davies layer.** Halves that layer,
   yields eigenvalues good to ~1e-7 relative, propagating into exactly the tail
   the saddlepoint work exists to fix.

---

## 5. Open questions (controls, in priority order)

1. **What does a post-#158 profile of a whole SKAT-O run actually look like?**
   Every ranking above rests on composing two separate PRE-#158 profiles. One
   paired sim50k run on a quiet box, both traits, `SAIGE_PROFILE_REGION=1`,
   settles whether the kernel is now 25% or 45% of the run -- and that decides
   whether items 4 and 5 are worth doing at all. Re-measure the
   `region_postprocess` double drive in the same pass.
2. **Is the landed saddlepoint tail accurate at q=4-10?** (Section 3.3.)
3. **Can carrier-space Phi reach the sparse-GRM branch, and what does it buy at
   N=400k?** 3282 of 9058 tested rows take that path, where Sigma^-1 is a
   per-marker solve; a sparse Cholesky solve is O(nnz(L)) not O(N), so the
   branch may be reachable at a smaller constant -- unmeasured. Decision rule
   should mirror the carrier-space saddlepoint's: a pinned All of Us chr20 pair,
   carrier vs dense, drop it if under ~10% there.
4. **Does the q DISTRIBUTION at real biobank scale change the ranking?** All of
   Us chr20 reaches q=482 vs sim50k's 134 (3.6x). Item 5 pays 2.04x at q=482 and
   nothing at q<=50, and Davies' per-call cost is non-monotone in q. It is the
   distribution over a real exome, not its maximum, that decides whether the
   eigen layer or the small-q Davies band is the target. **Nobody has plotted it.**
5. **What does a 100x change at p=1e-29 do to the FlexRV Cauchy pooled row?**
   Still open from the #158 handoff: #127's pooled row consumes SKAT-O p-values,
   the extreme-tail values have no reference, and the saddlepoint tail changed
   them. A CCT is dominated by its smallest input.

---

## 6. Caveats on the evidence

- **Independence is weaker than the vote counts look.** Seven of the twelve
  confirmed claims trace to Thomas Lumley (CRAN docs, bigQF vignettes, the CSDA
  and Genetic Epidemiology papers). Treat the fastSKAT/bigQF cluster as **one
  voice with internal consistency, not seven**.
- Two surviving claims were 2-1 splits (Davies/eps; leading-eigenvalue <1% error).
- The 3,500x Davies blow-up is at acc=1e-13 on a deliberately adversarial
  synthetic spectrum; flat at SKAT's acc=1e-6. Not a cost paid here.
- SMMAT-E's 2.6x is author-reported at n=23,763, single thread, against their own
  rho-grid implementation -- **not against R SKAT's `qfc`** -- with no 400k figure.
- Two verifiers exhausted their WebSearch budget, so "no contradicting source"
  for three claims rests on direct fetches, not a systematic adversarial search.
- The session's own benchmarks (`scratchpad/bench2.R`, `rank2.R`) use SYNTHETIC
  beta(1,25)-weighted Phi on a loaded box. They reproduce the repo's STRUCTURE
  (outer integral dominates; eigen layer grows with q; `Get_Lambda` already skips
  eigenvectors) but not its absolute per-cell milliseconds, which are 3-10x
  larger on real cells. **Use them for ratios within a row, never as absolute costs.**
- **Time-sensitivity: low.** Davies 1980, Farebrother 1984, Kuonen 1999 and the
  Lee 2012 rho envelope are settled numerics, and SKAT 2.2.5 still ships the same
  `qfc`. The one live area is the far tail: an MCMC-CE method (bioRxiv
  2025.03.16.643492, J Comput Biol 2026) claims better accuracy than
  Davies/Imhof/Farebrother/LTZ/saddlepoint for small right-tail probabilities but
  offers no quantified counter-figure and no released speed. Re-check in 6-12
  months, not now.

## Sources

- Chen H, Lumley T (2019), *Numerical evaluation of methods approximating the
  distribution of a large quadratic form in normal variables*, CSDA 139:75-81.
  https://www.sciencedirect.com/science/article/abs/pii/S0167947319301136
- fastSKAT: https://pmc.ncbi.nlm.nih.gov/articles/PMC6129408/
- bigQF (CRAN + SKAT projection vignette):
  https://cran.r-project.org/web/packages/bigQF/index.html ,
  https://cran.r-project.org/web/packages/bigQF/vignettes/skat-proj.html
- `survey::pchisqsum`: https://search.r-project.org/CRAN/refmans/survey/html/pchisqsum.html
- Chen et al. (2019) SMMAT, AJHG: https://pmc.ncbi.nlm.nih.gov/articles/PMC6372261/
- GMMAT (CRAN): https://cran.r-project.org/web/packages/GMMAT/index.html
- SKAT `SKAT_CompQuadForm.R`: https://github.com/cran/SKAT/blob/master/R/SKAT_CompQuadForm.R
- REGENIE `src/SKAT.cpp`: https://github.com/rgcgithub/regenie/blob/master/src/SKAT.cpp
- PGEN spec: https://github.com/chrchang/plink-ng/blob/master/pgen_spec/pgen_spec.pdf
- MCMC-CE far-tail preprint: https://www.biorxiv.org/content/10.1101/2025.03.16.643492v1.full
- In-repo: `docs/state/handoff-2026-09-17-skato-pvalue-methods.md`,
  `docs/state/handoff-2026-09-16-pgen-sparse-read.md`,
  `docs/state/region-perf-handoff.md`, `harness/bench/NOTES.md`,
  `R/SkatoScreen.R`, `R/SkatoParallel.R`, `R/SAIGE_SPATest_Region_Func.R`
