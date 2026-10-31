# PC-Relate natively: `harness/tools/pcrelate_to_mtx.py` (2026-10-31)

**Why.** KING-robust is biased on admixed pairs of differing ancestry: on AoU amr at 0.05 the
third-degree links joined 37k people into one indefinite component (LEDGER #183), and PMBB
(~80k, heavily admixed, ~1,000 traits) reports that its sparse GRM "does not work" and is
considering a dense fit, which at that N and trait count is the whole compute budget. The ladder
agreed with the owner: KING-robust at 0.125 first; if that fails, PC-Relate; the dense GRM last.
PC-Relate (Conomos et al. 2016) conditions kinship on the PCs and is the standard fix, but
GENESIS's `pcrelate()` returns every pair from an R-level block loop: 3.2e9 pairs at 80k, a day.

**What.** The estimator is two matrix products, `R R'` and `A A'` with `R = G - 2 mu` and
`A = sqrt(mu (1 - mu))` over the markers both samples carry with `0 < mu < 1`, `mu` from the
regression of each marker on the PCs over an unrelated set. The tool streams marker blocks against
sample tiles through numpy's threaded BLAS in float32, keeps only pairs at or above the cutoff on
SAIGE's GRM scale (2 x kinship), and writes the Matrix Market file SAIGE reads with the same
conventions as `king_to_mtx.py` (diagonal 1 by default, or 1 + the PC-Relate inbreeding). Input is
PLINK 1 (`--bed`, a numpy reader, no dependency) or PLINK 2 (`--pgen`, through the Pgenlib python
package; or `plink2 --make-bed`). No Spark, no GPU: plain numpy; a GPU would be a one-line swap of
the two products and is not needed (below).

**Controls (`selftest`, 15 s, every claim printed with its number).** A two-population
Balding-Nichols cohort (Fst 0.12) with admixed founders and planted relatives from haplotypes:
parent-offspring and full sibs 0.25, half sibs 0.125, first cousins 0.0625, cross-population
parent-offspring (the pair KING misreads), one inbred child of full sibs (f 0.25); PCs from the
unrelated founders with everyone projected (PC-AiR's shape).
- Positive: every planted class recovered (mean error within 0.005, max 0.016 at 16k markers);
  cross-population parent-offspring 0.241.
- Negative: unrelated pairs at 0 within and across populations (sd 0.0051, the closed form gives
  0.0059 at 8k and scales as 1/sqrt(M)); nothing at or above 0.125.
- The correction is load-bearing: the SAME estimator with the PCs removed reads across-population
  pairs at -0.064 and within at +0.062, 23k unrelated pairs past the cutoff.
- The null stated first: the Gaussian expectation past the 0.05 cutoff is 0.03 pairs; 3 are seen,
  every one involving an admixed founder, and halving the markers gives 4 -- the finite-marker
  PC-noise term (`sum_m R_im R_jm` picks up `4 eps_i eps_j sum beta_m^2` from the PC-score noise,
  individual-level and shrinking as 1/M; Conomos recommends many markers for this reason). Real
  runs use 100k+ pruned markers, where it is below the cutoff.
- Permutation within marker: 0 of 110,685 pairs past the cutoff (expected 0.03), mean -8e-4.
- Two ways of blocking (one block vs 1000 x 97) agree to 3e-7; 2 % missing genotypes move planted
  pairs by at most 0.004; `.bed` round trip and the PGEN arm (`plink2 --make-pgen` + pgenlib)
  reproduce the in-memory run to 2e-16; the inbred diagonal 1.23.
- **Independent reference**: `harness/tools/pcrelate_genesis_recon.R` runs GENESIS's `pcrelate` on
  the kept fixture (`selftest --keep DIR`) with the same PCs and training set: 110,685 pairs,
  max |diff| 7.2e-4, median 4.7e-5, correlation 0.99997. GENESIS builds only as conda binaries
  (the script's header has the recipe; it does not compile against the pixi R).
- `aou_sigma_pd_check.R --grmOnly` on the output: positive definite at 0.05 and 0.125.

**Speed (this 8-core laptop, load 8-40 from other work, so lower bounds).** Same input, same
machine: 471 x 16k, tool 0.8 s vs GENESIS 10 s; 4,011 x 16k, tool 9 s vs GENESIS 198 s (22x), and
GENESIS's cost grows with the pairs it returns. 20,000 x 50,000 with 10 PCs: 499 s, 7.2 GB
peak, 1.6 cores of BLAS engaged (OpenBLAS on macOS under load; the two GEMMs ran at 173 GFLOPS
in isolation). The flops are 4 N^2 M, so 80k x 100k is 32x that: a few hours here,
[ASSUMPTION] under an hour on a 32-core Linux VM whose BLAS scales. The prepare (decode,
regression, masks) cost as much as the GEMMs when written with float64 temporaries and is now
float32 in place; it runs once per sample tile, and `--memory-gb` sizes the tile (16 GB gives
25k at N = 80k, 4 passes).

**Use.** LD-pruned array markers (100k+), PCs from PC-AiR or plink2's PCA on a KING-unrelated set,
`--unrelated` that set, `--npcs` as fitted in step 1:
```
pcrelate_to_mtx.py run --bed pruned --pcs pcair.tsv --npcs 10 --unrelated unrel.txt --out cohort_pcrelate --cutoff 0.05
aou_sigma_pd_check.R --grmOnly --sparseGRMFile=cohort_pcrelate.mtx --sparseGRMSampleIDFile=cohort_pcrelate.mtx.sampleIDs.txt --relatednessCutoff=0.125 ...
```
then step 1 with `--sparseGRMFile`/`--sparseGRMSampleIDFile` as for the KING matrix.

**Open.** Not yet run on real data (no PC-AiR PCs for an AoU cohort are in the bucket; the
step-1 PCs are the All-by-All ones and would do for a first look at amr: compare its pairs with
`amr_king.mtx` through `king_to_mtx.py compare`). The universal-saige README's KING recipe
should name this as the next rung; a PR for the owner.
