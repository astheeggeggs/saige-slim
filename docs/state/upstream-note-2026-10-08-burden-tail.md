# BURDEN region tests with a sparse GRM do an O(N) projection per cell that nothing reads

Summary: in `mainRegionInCPP`'s BURDEN emitter, every (annotation x MAF x weight)
cell's `getMarkerPval` call passes `is_region = true` whenever a sparse GRM is
given (`src/Main.cpp:1849-1852` at `365deeff`). Inside `getMarkerPval` that flag
only requests the SKAT-O path's by-products after the test has finished -- the
dense covariate-projected vector (`getadjGFast`), `gy`, and `P2Vec`
(`src/SAIGE_test.cpp:926-931`), which feed P1/P2 for the variance matrix and the
Phi rescale. The BURDEN path keeps no P1/P2 and its writer reads `pval`, `Beta`
and `seBeta` alone, so each cell pays an O(N x covariates) projection for a
vector that is never read. The score itself is already carrier-space
(`scoreTestFast`) for cells above the sparse-GRM MAC band.

Measured (one fork, simulated N = 50,000, 210 genes, 192 custom weight sets per
gene, one thread, `--r.corr=1`, sparse GRM): the per-cell test phase was 13.0 s
of 24.0 s wall on a binary trait and 11.6 s of 22.2 s on a quantitative one.
With the flag false in that emitter: 6.5 s and 4.4 s (wall 16.7 s and 13.7 s).
What remains is the PCG solve for cells inside the sparse-GRM band and the
saddlepoint/exact-test cells, which is real work.

Fix: in the BURDEN emitter, pass `is_region = false` (the saddlepoint, exact
test, Firth and conditional branches build their own projected vector when they
fire, so no result changes). Optionally also pass an empty zero-index vector:
the exact test builds the complement itself from an empty one, in the same
order, so the per-cell O(N) `complementIndices` can go as well.

Verification on the fork: 30,167 cells including 1,703 saddlepoint cells, 489
exact-test cells and 99 cells routed to the exact test by weight, byte-identical
before and after; the fork's full golden set byte-identical; ThreadSanitizer
clean at 1/4/8 threads.

The rule dates from the 2022 commit `3a3c5506` ("marker list and burden test
only"), so every SAIGE-GENE+ BURDEN run with a sparse GRM since then has paid it.
