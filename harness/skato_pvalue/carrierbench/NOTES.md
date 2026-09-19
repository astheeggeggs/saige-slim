# carrierbench -- dense vs carrier-space assembly of the SKAT-O variance matrix

Standalone, gitignored (`harness/work/`), touches nothing under `src/`, `R/`, `tests/`.

## What is benchmarked

One region's q x q variance matrix, built from q genotype columns, covariate
adjustment included, on the **diagonal-Sigma path** (the 64% of rows whose MAC is
above `cateVarRatioMinMACVecExclude.back()` = 20.5 and therefore skip
`getPCG1ofSigmaAndGtilde`). Both arms compute

    gtilde_j' W gtilde_l = g_j'W g_l - (X'W g_j)' (X'W X)^-1 (X'W g_l),   W = diag(mu2*tau)

i.e. the census's `varmat` bucket (Main.cpp:4143-4218) PLUS the per-marker P1/P2
formation that feeds it (`rareArm`/`urPass`), which is the O(N) work a
carrier-space change would remove.

Naming: the census established that there is no per-cell Phi construction in C++ at
all -- Phi is an O(|cell|^2) reweighted subset of the region-level `VarMat`, taken in
R at SAIGE_SPATest_Region.R:511-516. The object that costs O(N) is VarMat. This
bench builds VarMat and calls it that.

## Faithfulness of the dense arm (so it is not a strawman)

- BLAS3 `cblas_dgemm` over (m x N)(N x m), OpenBLAS, `OPENBLAS_NUM_THREADS=1`.
- `P1` is a column-major (100 x N) matrix written by **strided rows**; `P2` is
  (N x 100) written by **contiguous columns** -- the layout at Main.cpp:2868-2869,
  :3240-3241, :4091-4092.
- Chunked at 100 markers with the nchunks x nchunks triangle of Main.cpp:4160-4214,
  including the `submat` write-backs and the `.t()` mirror.
- Chunks copied BY VALUE into a store, as `RegionChunkStore::add` does (Main.cpp:2399).
- `getadjGFast`'s **split** form (SAIGE_test.cpp:545-554): carrier-space forward
  gather for `X'Wg`, dense N x p `dgemv` back-projection, dense N subtraction.
  Not a fully dense `XV * GVec`.
- p = 7. `m_isVarPsadj` FALSE, so the extra O(N p) correction of
  SAIGE_test.cpp:308-310 is NOT included, per the census.
- Compiled with the package's own flags (`src/Makevars` + R's `CXX14FLAGS`):
  `-march=core2 -mtune=haswell -mssse3 -ftree-vectorize -O2`. This matters: the
  compiled C++ gets no AVX2/FMA while OpenBLAS dispatches AVX2 kernels at runtime.
  That asymmetry is real in the shipped package and is the LEDGER #156 hazard in
  its concrete form.

`A_syrk` is a second dense arm: one (N x q) matrix `sqrt(w).*gtilde` and a single
`cblas_dsyrk`. Not what ships (the shipped code does a full `dgemm` on the diagonal
blocks, so it pays 2x for a symmetric result), but it is the cheap dense lever any
carrier-space prize has to beat.

## Carrier arms

- `B_merge` -- as specified in the task: sorted carrier-list INTERSECTION merge per
  pair, O(nnz_j + nnz_l); `U = X'WG` by gather O(p*nnz); then `Phi -= U'(X'WX)^-1 U`.
- `B_samp` -- same maths, sample-major (SpGEMM) accumulation, cost `sum_i c_i^2`
  instead of `q^2 * nnz`. Included so a `B_merge` loss cannot be mistaken for a
  carrier-space loss: for near-disjoint rare columns `c_i` is 0 or 1 and this is
  O(nnz_total), asymptotically better than the pairwise merge. Its CSR transpose
  from marker-major to sample-major is timed as part of the arm, because the reader
  hands SAIGE marker-major lists.

## Grid

From the parameter census: N in {3000, 50000, 400000}; q in {4, 11, 30, 134}
(measured quantiles, max); nnz/N in {5e-5, 5e-4, 1.2e-3, 5e-3, 2e-2} (MAF-bin
reality: rare surviving columns 4-5e-4, pseudo-marker 1.2e-3, 1e-2 bin boundary
2e-2) plus **0.5 as the negative control** (common variants; carrier space must
lose there).

## Controls

1. Correctness before timing: all four arms agree on Phi to < 1e-10 relative
   (achieved ~1e-15). This caught a real triangle bug in `B_samp` on the first run.
2. Negative control at nnz/N = 0.5.
3. Dense `dgemm` GFLOP/s reported, as a positive control that the dense arm is
   BLAS3-efficient (census estimated 13-19 GFLOP/s for the real code).
4. Warmup + 3 timed reps, median reported; `uptime` recorded either side.

## Not modelled (and why it bounds the result)

- The 36% of rows on the sparse-GRM PCG path. There `Sigma^-1 g` is a PCG solve
  that fills in, so `P2Vec` is dense whatever the genotype is; the q^2 N Gram term
  can still be gathered against a dense P2, but the per-marker O(N) solve cannot be
  removed at all. Carrier space is strictly inapplicable to that work.
- The 18.6% inflation from fastTest escalation re-running 39 of 210 regions.
- Cell-level work in R (`get_newPhi_scaleFactor_traitType`, SPA), which is
  downstream of VarMat.
