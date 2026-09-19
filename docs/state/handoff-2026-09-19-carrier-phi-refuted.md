# Carrier-space Phi: layer win huge, run win 1.03x. DO NOT BUILD IT (2026-09-19)

`docs/state/research-2026-09-18-skato-speedups.md` ranked "carrier-space Phi by
the four-term expansion" FIRST, on the reasoning that it attacks the O(N) C++
work that is ~54% of the run and grows with N. **Benchmarked: the layer
speedup is real and enormous, and the run-level prize is 1.034x at N = 50,000
and 1.06x at N = 400,000. The ranking was wrong.**

Benchmark `harness/skato_pvalue/carrierbench/` (standalone C++, OpenBLAS,
single core, package compile flags). Raw sweep in
`harness/work/carrierbench/` (gitignored).

## Why the ranking was wrong -- a located error, not bad luck

The note attributed the whole 54% C++ share to this lever. It does not touch
54%. **`varmat` -- the ONLY bucket this changes -- is 3.365 s of the profiled
96.69 s: 6.4% of C++ and 3.5% of the run.** The other O(N) buckets
(`rareArm` 19.6 s, `urPass` 21.2 s, `test` 14.0 s) are the per-marker
`gtilde`/`P2` construction, the PCG solves and the tests -- different work,
untouched by a carrier-space Gram.

Two further corrections from the code census, each `file:line` verified:

- **There is no per-cell Phi construction in C++ at all.** Phi is built in R by
  subsetting ONE region-level marker x marker `VarMat` and applying the cell's
  weights (`R/SAIGE_SPATest_Region.R:511-516`), which is O(cell^2) with no N
  term. The O(N) object is `VarMat` (`P1 * P2`, `src/Main.cpp:4193`/`:4207`,
  armadillo dgemm at 13-19 GFLOP/s), built ONCE PER REGION. "Carrier-space
  Phi" is the wrong name for the lever; the flop model per region was right.
- **Sigma^-1 is NOT diagonal for 36% of the rows on this very recipe.**
  `contextForMAC` picks per marker by MAC (`src/SAIGE_test.cpp:1213`); MAC <=
  20.5 takes `getPCG1ofSigmaAndGtilde` (`:305`), a Jacobi-preconditioned CG
  with a sparse N x N SpMV per iteration. Measured from the profile run's own
  output: **3282 of 9058 rows (36.2%)**, reproducing the note's figure
  exactly. Those are the RAREST markers -- i.e. precisely the ones a
  carrier-space Gram would help most are the ones it cannot reach, because
  the PCG solve is dense by construction and common to both arms.

## What the benchmark measured

Correctness control first: four arms agree on the Gram to **~1e-15** relative.
It earned its keep -- the sample-major arm first disagreed at 1.1e-2 (a
triangle-mirroring bug) and a timing-only bench would have crowned it fastest.

Negative control (nnz/N = 0.5, common variants): carrier space **loses in all
12 cells** against a fair dense arm, by 7-26x at q >= 30. As it must.

Layer speedup at the realistic operating points (N = 50,000, rare columns at
nnz/N ~ 5e-4), pairwise merge vs the SHIPPED dense arm: **357x at q = 11**
(75th pct), **184x at q = 30** (90th pct), **56x at q = 134** (max). Against a
dense arm given one afternoon's work (`dsyrk` + contiguous writes) it is 98x /
47.6x / 19.5x.

**And LEDGER #156's mechanism is present and quantified**: a merge step costs
4.0-5.3 ns (~10 cycles, unpredictable branch, no vectorisation) against
0.028-0.09 ns per dgemm flop -- **the flop count overstates by 45-190x**. The
model reconciles once the per-element penalty is in it (predicted 27x on the
Gram at q=134, measured 22.7x). The >=50x claim survived only because it
started 2000-20000x up rather than #156's 24x.

## The number that decides it

| | N = 50,000 | N = 400,000 |
|---|---|---|
| `varmat` share of run | 3.5% | 5.8% |
| **carrier-space Gram free** | **1.034x** | **1.062x** |
| ALL SKAT-O-only C++ free (the class ceiling) | 1.30x | 1.45x |
| whole R kernel free | 1.72x | 1.09x |

Two independent routes agree on 1.034x/1.036x -- the benchmark's measured layer
ratio applied to the bucket, and the census's bound from the bucket share.

**Conclusion: do not build it.** A 357x layer speedup on 3.5% of a run is
1.03x. It is also carrier space's BEST case here -- an all-diagonal-Sigma
region, while 36% of real rows take the PCG.

## Two side findings worth keeping

1. **Pairwise merge is the wrong carrier algorithm** if this is ever revisited:
   O(q^2 * nnz) against sample-major's O(N + sum_i c_i^2), measured **10x**
   apart at q = 134, nnz/N = 5e-4 and **26x** at nnz/N = 2e-2.
2. **The shipped dense arm leaves 2.9-3.9x on the table with NO algorithm
   change** -- contiguous writes instead of the strided `P1` row, and `dsyrk`
   instead of a full `dgemm` for a symmetric result. The `dsyrk` half is
   BLOCKED: `Sigma^-1` is chosen per marker, so one `VarMat` mixes diagonal and
   PCG rows and the `sqrt(w)` factorisation does not hold; `src/Main.cpp:4195-4196`
   imposes symmetry by fiat instead. The contiguous-write half is unaffected and
   independent of any carrier-space work -- but it is 2.9-3.9x on a 3.5% bucket,
   so it is worth ~1.02x and is not a priority either.
   [ASSUMPTION: that symmetrisation-by-fiat on a mixed region is pre-existing
   rather than intended; three sites read, not traced against the goldens.]

## Standing conclusion for the SKAT-O path

With #158 landed (2.14x, measured), the profile (kernel 42%) and this
benchmark, **there is no large accuracy-preserving speed-up left on the SKAT-O
path.** The ceilings are: whole kernel free 1.72x at 50k / 1.09x at 400k; all
SKAT-O-only C++ free 1.30x / 1.45x. Everything specific that has been costed
lands at or under those, and the biggest single identified lever is worth 1.03x.
