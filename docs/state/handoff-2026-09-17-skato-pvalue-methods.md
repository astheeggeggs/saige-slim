# Handoff 2026-09-17: SKAT-O p-value methods, measured on both traits

Question (owner): can the SKAT-O arm be sped up the way BURDEN was, is a
Davies drop-in available, what does modified Liu buy, and is it safe for
meta-analysis. Harness: `harness/skato_pvalue/`; raw records in the session
scratchpad only (`prof/`, `liu/`, `cells/`, `cells_q/`). Fixture sim50k, pgen
reader, y_binary and y_quant models, the `04_run_saige.sh` group recipe
(sparse GRM, fastTest, single-in-group, 5 annotation groups x 3 max-MAF),
1 thread, box under load 4-13 (shares are robust, seconds are not).

## 1. Where a SKAT-O run spends its time (C++ buckets, seconds)

| bucket | SKAT-O | BURDEN | what the difference is |
|---|---|---|---|
| test (single-variant tests) | 14.1 | 14.0 | same on both arms |
| rareOther | 6.8 | 0.1 | scaled copies of gtilde/P2 into the block matrices |
| flush | 4.4 | 0.0 | copying the block into the chunk store |
| urPass minus urTest | 13.1 | 1.6 | per-cell residualised burden vector, N x p, for R |
| varmat | 3.8 | 0 | the P1.P2 product, O(m^2 N) |
| C++ total | 56.7 | 35.0 | |
| user CPU, whole run | 207 | 35 | SKAT-O's R side is ~150 s |
| peak RSS | 808 MB | 394 MB | the P1/P2 store |

The sparse read already feeds both arms (`c05f23fb`). The SKAT-O-specific
O(N) work in C++ is ~27 s = 13% of the run; a carrier-space VarMat (four-term
expansion of gtilde_j' W gtilde_l, no N term) could remove it and the
P1/P2 memory, but is not byte-neutral, cannot reach the sparse-GRM branch
(Sigma^-1 per marker; 3282 of 9058 tested rows on this recipe), and caps at
~1.15x at N = 50k. At biobank N the O(N) share grows ~8.5x while Davies does
not; unmeasured. The RAM saving is the concrete prize (All of Us 3.8 GB vs
upstream 0.9 GB fits this by scaling; [ASSUMPTION], not profiled there).

## 2. Davies never sees N; the cost is the outer rho integral

`SKAT:::Met_SKAT_Get_Pvalue(method = "optimal.adj")` has two layers: eight
per-rho p-values via Davies at acc 1e-6 / lim 1e4 (value KEPT whatever
`ifault` says; replaced by modified Liu only when outside (0, 1]), then the
Lee 2012 integral `1 - int (1 - S(delta_min(x))) f_chi2_1(x) dx` with Davies at
acc 1e-6 inside every node, falling back to the Liu-form integral if any node
fails, capped at `8 * pmin`, and replaced by `3 * min(p_rho)` when `1 - int`
underflows to 0. Per-call cost is NOT monotone in q: the q 4-10 band costs
138 ms per call and carries 56% of Davies time; q > 100 carries 2%. fastSKAT
(a q lever for q in the thousands) is therefore not on the table for gene
tests; largest q here 134, on All of Us chr20 482.

## 3. Modified Liu, every cell (1613 binary multi-rho cells, in-run timing)

Davies 142 s total (49 ms/call) vs Liu 13 s (4.5 ms/call). As a swap it is
wrong where it matters -- Liu/Davies ratio by Davies p band: (1e-8,1e-6]
0.02-2.0; (1e-6,1e-4] 0.13-1.0 median 0.86. As a SCREEN it is clean:

| skip Davies when Liu p > | cells skipped | Davies time saved | misses |
|---|---|---|---|
| 0.05 | 83% | 91% | 0 (no skipped cell had Davies p < 0.005) |
| 0.01 | 89% | 96% | 0 (worst skipped Davies p 0.0097) |
| 0.001 | 93% | 98% | 0 (worst skipped 0.0012) |

Run-level at 0.01: kernel 142 s -> ~18 s, run 204 s -> ~80 s, **~2.5x**,
byte-identical below the threshold by construction. Same zero-miss result on
the quantitative arm (1443 of 1506 cells skipped, 0 misses at reference
p < 1e-3). The test for it must assert the miss COUNT on the goldens and on a
null sweep, not assume it.

## 4. Accuracy against a reference, both traits

Reference: Davies with an ADAPTIVE absolute tolerance (acc = 1e-4 x a
saddlepoint pilot, lim 1e8), both layers strict, split-integral form
`P = P(chi2_1 > x*) + int_0^{x*} S f` (Wu et al. 2016's formula (3) shape; the
x = t^2 substitution removes the density's singularity at 0). Valid where
every inner Davies converges: p from 1e-2 down to ~1e-9. Below that Davies
returns `ifault` 1/2 at any acc/lim -- checked on synthetic spectra against
Ruben's exact positive series, which agrees with Davies to 5 digits at 1e-10
and continues to 1e-23 where Davies fails; on real cells the series never
converges (SKAT keeps eigenvalues down to mean/1e5), so the extreme cells
have NO exact reference. The reimplementation reproduces SKAT's p bit-for-bit
on 1630/1645 binary and 1497/1506 quantitative cells (rest at 1e-13 relative).

Ratio to the reference, median (min-max), binary trait, 153 referenced cells:

| reference p band | SKAT as shipped | modified Liu | saddlepoint, direct | REGENIE tier |
|---|---|---|---|---|
| (1e-12,1e-8] n=2 | 10.6 | 0.21 | 1.12 | 1.05 |
| (1e-8,1e-6] n=22 | 1.32 (0.05-7.5) | 0.37 (0.14-1.09) | 1.11 (1.04-1.14) | 1.10 (0.78-1.44) |
| (1e-6,1e-5] n=8 | 0.99 (0.92-1.17) | 0.56 (0.40-0.91) | 1.08 (1.07-1.09) | 1.12 (1.08-1.14) |
| (1e-5,1e-4] n=28 | 1.01 (0.93-2.26) | 0.83 (0.53-2.24) | 1.09 (1.04-1.12) | 1.01 (0.93-2.26) |
| (1e-4,1e-3] n=34 | 1.01 (0.97-1.13) | 0.85 (0.64-1.13) | 1.09 (1.01-1.11) | 1.01 |
| (1e-3,1e-2] n=59 | 1.00 (0.99-1.02) | 1.00 (0.94-1.06) | 1.07 (0.94-1.09) | 1.00 |

Quantitative trait, 59 referenced cells (p 1.6e-5 to 1e-2): SKAT 1.00
(0.96-1.02) except one capped cell at 2.30; Liu 0.66 median at (1e-5,1e-4]
(0.54-2.30), 0.95 at (1e-4,1e-3]; saddlepoint 1.05-1.06 (1.01-1.12). No
quantitative cell lies below 1e-5 on this fixture, so the binary arm is the
only evidence for the tail.

Readings:
- **SKAT as shipped is off by up to 20x IN EITHER DIRECTION between 1e-8 and
  1e-6 on the binary trait** (0.05 to 7.5), i.e. at discovery thresholds.
  The per-rho layer is the cause: SKAT's min-over-rho p is 0.22-3.5x the
  reference's there, because Davies at acc 1e-6 returns junk or a
  non-positive value for a per-rho p below ~1e-6 and the Liu substitute is
  what survives; the outer integral inherits it. Above 1e-6 both layers are
  exact. On the quantitative arm the per-rho layer is exact (1.000) because
  no cell is below 1e-5. Wu et al.'s Table 1 "SKAT 0.92/0.66" is this.
- **The extreme cells (SKAT 1e-16 to 1e-48) are Liu values**: SKAT's p equals
  the `liu.mod` p to every digit there (per-rho Liu substitution, then the
  `8 * pmin` cap). The saddlepoint puts them 100-10,000x higher. Neither is
  verified; nothing reads a p there anyway.
- **Modified Liu is anticonservative by 3-7x at 1e-8..1e-6, 2x at 1e-6..1e-5,
  and only converges to Davies above 1e-3** -- Wu's 4.3x at 1e-8, reproduced.
- **The saddlepoint is 5-12% conservative and STABLE across every band and
  both traits**, with no failures in the direct form (1.04-1.14 below 1e-5).
  Its 1 - int form shows outliers (0.78, 1.44) that the direct form does not.
  This is Wu's "SP 0.91" column and Chen & Lumley's bounded relative error.
- **Davies at acc 1e-9 / lim 1e5 (REGENIE's last resort) fails to converge
  on most cells above 1e-5 and on 6 of 22 in (1e-8,1e-6]**; where it converges
  it is 0.75-1.11. It is a fallback, not a method.
- Timing per deep cell (R implementations, ms): SKAT 68, Liu 10, saddlepoint
  86 (direct 234), tight Davies 31 where it converges. The saddlepoint's cost
  here is a uniroot per node in R; a compiled one would be a few ms.

## 5. Literature, and what upstream does

REGENIE: Davies by default; Kuonen saddlepoint when Davies' p < 1e-5; then
Davies lim 1e5 / acc 1e-9 if that fails; citing Wu, Guan & Pankow 2016
(Ann Hum Genet, PMC4761292), who show Liu 9.1x / modified Liu 4.3x /
saddlepoint 0.91 / SKAT-default Davies 0.66 of nominal at alpha 1e-8 and
propose Davies acc 1e-9 with saddlepoint fallback plus the P - T integral
form (code: github.com/baolinwu/SKATR). Chen & Lumley 2019 (CSDA): Davies
usable to machine epsilon in principle, saddlepoint bounded relative error
in the extreme tail. Zhang & Shen 2020 (arXiv 2005.00905): LTZ and LTZ4
inflate type I error whenever skewness and kurtosis cannot both be matched.
fastSKAT (Lumley 2018) targets q in the thousands. Meta-SAIGE and REMETA
meta-analyse per-variant scores plus an LD matrix and never see a per-cohort
gene p-value. **plink2 has no SKAT-O** (plink 1.x `--set-test` is unrelated).
**Upstream SAIGE `main` calls the identical `Met_SKAT_Get_Pvalue(method =
"optimal.adj")`** (`main:R/SAIGE_SPATest_Region_Func.R:322`), and
`harness/sim/upstream_compare/RESULTS-2026-09-09.md` found the gene region
file byte-identical between strip and upstream on this fixture -- so a
comparison to upstream is a null control: same numbers, same defects.

## 6. Meta-analysis

Score-level meta (Meta-SAIGE, REMETA) is untouched by any of this, but needs
the LD matrices this branch stripped (step 3). P-value-level meta (Stouffer,
Fisher, Cauchy) sees the reported p directly: a Liu swap would inject 3-7x
anticonservative errors at exactly the genes being combined; the screen
design leaves every p below the threshold to the confirming method. Note
#158's standing point: the owner's Stouffer needs a SIGNED Z, which the
directionless SKAT-O p cannot give -- if SKAT-O is not consumed downstream,
`--r.corr=1` removes the whole cost with a flag.

## 7. Recommendation (not implemented; nothing under R/ or src/ changed)

1. Modified Liu as a screen at 0.01 (or 0.05) in front of the existing
   kernel: ~2.5x on SKAT-O runs, byte-identical below the threshold, both
   traits, zero misses measured. Gate: the miss count on the goldens and on a
   null sweep, both traits.
2. Replace the confirmation below ~1e-5 with the saddlepoint in BOTH layers,
   direct integral form (REGENIE's tier, Wu's integral): fixes the 20x-either-
   way band that SKAT ships today, at a known 5-12% conservative bias, no
   failures. Gate is numeric, not byte: the sub-1e-6 goldens SHOULD move.
   Needs SKAT's SKAT-O integrand reimplemented in this tree (`variants.R`
   has it, controlled) since SKAT is a dependency.
3. Do not swap Liu in wholesale. Do not use Davies acc 1e-9 as anything but a
   last resort.

Open: a compiled saddlepoint (the R one costs as much as Davies); the
extreme-tail values have no reference and the FlexRV Cauchy row consumes
SKAT-O p's (#127) -- check what a 100x change at 1e-29 does to p_FlexRV.

## 8. Cost of each option (2026-09-18), and upstream

Per-cell kernel cost, ms, median over the deep cells (R implementations; the
saddlepoint is a `uniroot` per node in R, a compiled one would be a few ms):

| | binary | quantitative |
|---|---|---|
| SKAT as shipped (Davies, both layers) | 48 | 61 |
| modified Liu | 8 | 11 |
| saddlepoint, 1 - int form | 83 | 101 |
| saddlepoint, direct form | 251 | 319 |
| Davies acc 1e-9 / lim 1e5 (where it converges) | 139 | 673 |

Whole-run kernel cost, seconds, composed from the recorded per-cell times
(1645 binary / 1506 quantitative multi-rho cells; the kernel is ~70% of a
SKAT-O run, so run-level factors are ~2.5x where the kernel factor is ~6x):

| option | binary | vs shipped | quant | vs shipped |
|---|---|---|---|---|
| A SKAT as shipped | 166 | 1x | 158 | 1x |
| B Liu screen 0.01, then SKAT | 27 | 6.1x | 23 | 6.9x |
| C B, then saddlepoint below 1e-5 (REGENIE tier) | 30 | 5.5x | 23 | 6.9x |
| C with the direct-form saddlepoint | 42 | 3.9x | 24 | 6.6x |
| REGENIE tier without the screen | 169 | 1.0x | 158 | 1.0x |
| saddlepoint everywhere (R) | 141 | 1.2x | 162 | 1.0x |
| Liu everywhere (wrong in the tail) | 21 | 8.0x | 19 | 8.2x |

Stage counts, binary: 200 of 1645 cells pass the 0.01 screen, 82 are below
1e-5, 50 below 1e-8; quantitative: 63 / 4 / 4. So the screen is the speed and
the saddlepoint confirmation is the accuracy; the confirmation costs ~3 s a
run in R and nothing measurable compiled.

**Upstream does not read in carrier space.** `main:src/Main.cpp` has no
sparse/difference-list read (0 hits); every variant is expanded to a dense
N-vector on both its readers. The clean paired measurement is
`harness/sim/upstream_compare/RESULTS-2026-09-09.md` (bed reader, load 3-4,
3 reps): gene tests upstream/strip 1.38x binary Firth on, 1.35x Firth off,
1.37x quantitative -- BEFORE the carrier-space read, which then added a
measured 1.13x on SKAT-O and 1.56x on BURDEN (pgen, byte-identical), so the
composite against upstream today is ~1.5x SKAT-O and ~2.1x BURDEN
[ASSUMPTION: multiplicative; not measured as one pair]. A pgen pair run
overnight 2026-09-18 (upstream `--pgenPrefix --AlleleOrder=ref-first`, 1416
user-s; strip 1320 user-s) is NOT a measurement: Docker's backend was at
~190% CPU throughout and the strip command that took 207 user-s at 20:07
took 1320, a 6.4x inflation that need not be equal across arms. Re-run on a
quiet box before quoting; script `harness/skato_pvalue/upstream_pair.sh`.
Peak RSS on sim50k: upstream 666 MB, strip 742 MB (the RegionChunkStore
trade is small at this N). P-values, upstream vs strip on that pair, 3135
non-Cauchy cells: median |dlog10 p| 0.0017, p90 0.010, but 33 of the 69
cells below 1e-5 differ by more than 0.1 in log10 (max 1.14, 9.7e-10 vs
1.3e-8) -- the SKAT kernel is identical, so these are the branch's ledger
fixes since 2026-09-09 (#118, #130, #133-#135 and the allele-order flag),
not decomposed here.

## 9. Null replicate study (2026-09-18): screen misses and calibration, three arms

Design: 32 replicate null phenotypes from the fixture's own generator
(`harness/sim/flexrv_study/make_null_series.py -R 32`; seeds reproduce r1-r5,
so those models were reused), each also thresholded to a binary null at
K = 0.01 (500 cases) and K = 0.002 (100 cases), the `y_null_bin` recipe.
91 step-1 fits, 96 SKAT-O runs (the `04_run_saige.sh` group recipe, pgen,
1 thread, every cell dumped), offline scoring of every kernel call: modified
Liu on all, saddlepoint (both forms) and the adaptive-Davies reference on
calls with SKAT p < 1e-2. Scripts `harness/skato_pvalue/nullrep/`; data
`harness/work/nullrep-20260918/` (gitignored). ~2.5 h on 8 cores.
Caveat on the fits: tau2 went to 0 on 4 of 32 replicates at 500 cases and
17 of 32 at 100 cases (100 cases cannot inform a variance component), so
those runs are logistic-with-SPA nulls; results are given both ways.

**Screen: zero misses in every arm.** Over 81k kernel calls per arm, no
call with Liu p > 0.01 had SKAT p < 1e-3 or reference p < 1e-3; 1.7-2.3%
of calls pass the screen. Combined with the two signal runs, the screen's
false-negative bound now rests on ~250k null cells and ~4.5k signal cells.

**Inversion methods on the referenced null calls** (SKAT p < 1e-2, 5%-95%
of the ratio to the reference): SKAT as shipped 0.99-1.01 (max 1.29), Liu
0.84-1.06 (min 0.51 at 100 cases), saddlepoint direct 0.97-1.09. At 1e-4 the
counts were SKAT 10 / Liu 14 / saddlepoint 8 / reference 10 (500 cases) and
44 / 52 / 40 / 42 (100 cases): Liu over-calls, the saddlepoint slightly
under-calls, SKAT tracks the reference. No null call reached 1e-6, so the
band where SKAT as shipped fails is not reachable by a genotype-level null.

**Calibration of the statistic itself, from the output files** (100,320
cells per arm, Cauchy rows excluded; identical for every inversion method
because they share Score and Phi). Note the dumped kernel calls are NOT the
cell set: fastTest escalation re-evaluates significant regions, so a
call-level ECDF over-counts small p (6.6% at 0.05 on the quantitative arm
against 5.3% in the file). File-level, empirical rate at alpha:

| arm | 0.05 | 0.01 | 1e-3 | 1e-4 (n, expected 10) |
|---|---|---|---|---|
| quantitative, SKAT-O | 0.053 | 0.011 | 8.2e-4 | 4 |
| binary 500 cases, SKAT-O | 0.050 | 0.011 | 1.5e-3 | 8 |
| binary 500 cases, Burden / SKAT | 0.045 / 0.052 | 0.010 / 0.012 | 1.1e-3 / 1.2e-3 | 8 / 2 |
| binary 100 cases, SKAT-O | 0.050 | 0.014 | 1.8e-3 | 27 |
| binary 100 cases, Burden / SKAT | 0.037 / 0.059 | 0.010 / 0.015 | 1.1e-3 / 2.0e-3 | 18 / 19 |

Quantitative: calibrated (5.5% at 0.05 for q 2-10 and q > 30, 5.0% at
q 11-30 and q = 1). 500 cases: calibrated to 0.01, 1.5x at 1e-3, on
target at 1e-4. **100 cases: the SKAT component is inflated across the
range (5.9% at 0.05, 2x at 1e-3, 1.9x at 1e-4) and the omnibus inherits it
(2.7x at 1e-4, 3 cells below 1e-5 against 0.8 expected); the burden is
conservative at 0.05 (3.7%, SPA discreteness) and 1.8x at 1e-4.** It is
present in both the tau2 = 0 and the tau2 > 0 replicates (1.9e-4 and
3.6e-4 at 1e-4), so it is not the variance-component collapse; it is the
kernel's Phi adjustment at 100 cases, upstream of every option here and
shared by the current code. Per-replicate counts below 1e-3 range 0-38
(expected 2.5) because the 15 cells of a gene are correlated, so
gene-level clustering, not replicate-level, is the right unit for a CI.
**Consequence for the decision: no option changes calibration relative to
today; the 100-case inflation is a separate, real item** for the SPA/ER
Phi adjustment: LEDGER #159, DECIDED 2026-09-18 by the owner, leave as is,
noted. The pinned All of Us phenotypes with few cases are where it should
be checked next.
