# strip vs upstream at N = 50,000: speed and p-values (2026-09-05)

> **SUPERSEDED IN PART -- see `RESULTS-2026-09-09.md`.** The timing table below predates
> commit `b11e6b07` (2026-09-06); the single-variant rows in particular have moved
> (binary 0.94x -> 1.78x, quant 1.28x -> 1.99x). Step 1 and the gene-test ratios still
> hold. The gene-test p-value rows predate LEDGER #103-#106 (2026-09-08); the 2026-09-09
> run reproduces this file's binary figure of 0.018 by re-running the pre-fix code, so
> what is below is correct FOR THAT CODE, not for HEAD.

`strip` HEAD (`fa0ec961` + this session's uncommitted comparison scripts; `SAIGE.so`
md5 `41a20d16…`) against upstream `main` 365deeff (SAIGE 1.5.2, `SAIGE.so` md5
`2ec37f64…`), both on the `sim50k` fixture (`SAIGE-sim/harness/work/sim50k`: 50,000
samples, chr20, 210,515 markers, 210 genes, pedigree-truth sparse GRM). How the arms
were built and run: `README.md` in this directory; every command is `run_arm.sh` over
`harness/sim/04_run_saige.sh`, so the flags are identical between arms.

**Confounder controls.** Upstream has its own pixi env (it needs `savvy`, `MetaSKAT`,
`RSQLite`, `qlcMatrix`, `dplyr`, absent from strip's) pinned to the SAME conda builds
as strip's: R 4.4.3 `h46d4a84_2`, OpenBLAS 0.3.30 `openmp_h30af337_0`, clang
`x86_64-apple-darwin13.4.0-clang++`. Every run's log records the library R actually
loaded and the `SAIGE.so` md5 (`results.tsv` columns `lib`, `so_md5`); the two arms
never loaded each other's build. Both arms reproduce themselves byte for byte over the
three timing reps (all 8 step-2 outputs x 3 reps, `cmp`), and the binary step-1 models
of the two arms are bit-identical (`theta`, coefficients, fitted values: 0 difference),
so the model files are interchangeable -- the smoke test in `README.md` fed each
build's models to the other build and got the same 213 rows back.

## 1. Speed, 1 thread, paired (3 reps, arms adjacent, order alternated per rep)

Frozen step-1 models (the fixture's committed ones) for step 2, so both arms test the
identical null model. User CPU seconds from `/usr/bin/time -l`; the ratio is
upstream / strip, so > 1 means strip is faster; the paired t-test is on 3 reps.
`CPU_Speed_Limit = 100` throughout; the only other CPU consumer was Spotlight's
`mds_stores` (~30%, indexing the new output files) -- shared between the arms, which
is what the pairing is for.

| workload | upstream user s | strip user s | upstream / strip | paired 95% CI on the difference | maxRSS upstream / strip |
|---|---|---|---|---|---|
| step 1, `y_binary` (sparse GRM fit + categorical VR) | 9.03 ± 0.06 | 6.58 ± 0.03 | **1.37x** | [2.38, 2.53] s, p = 5e-5 | 430 / 389 MB |
| step 1, `y_quant` | 8.95 ± 0.04 | 6.68 ± 0.05 | **1.34x** | [2.15, 2.38] s, p = 1e-4 | 441 / 391 MB |
| step 2 single-variant, `y_binary` (70,701 markers tested, `--minMAC=20`, Firth 0.1) | 236.3 ± 3.5 | 252.1 ± 14.2 | **0.94x** | [-45.7, +14.1] s, p = 0.15 | 377 / 2342 MB |
| step 2 single-variant, `y_quant` | 80.2 ± 3.4 | 62.5 ± 2.4 | **1.28x** | [15.0, 20.5] s, p = 1e-3 | 344 / 2789 MB |
| step 2 gene tests, `y_binary` (210 genes x 5 anno x 3 maxMAF, SKAT-O, single-in-group) | 222.5 ± 8.5 | 167.5 ± 1.5 | **1.33x** | [30.7, 79.3] s, p = 0.01 | 627 / 1072 MB |
| step 2 gene tests, `y_quant` | 201.7 ± 13.4 | 147.3 ± 4.5 | **1.37x** | [13.2, 95.6] s, p = 0.03 | 622 / 1006 MB |

Single-variant `work_s` (startup removed by regressing cumulative `proc.time()` on the
chunk ordinal, R² ≥ 0.994): upstream / strip = 0.935x binary, 1.273x quant; per
marker 1244 vs 1333 µs (binary), 383 vs 300 µs (quant). The quantitative figure
reproduces the 09-01 comparison (1.37x bed then; 1.27x now, different day, different
HEAD, so not a regression claim without a same-day A/B).

**Two things this table says that were not known:**

- **The binary single-variant scan is not faster than upstream.** 0.94x on the
  mean, the CI straddles zero, and the strip arm's spread (±14 s) is 4x upstream's.
  Every Phase 4 marker-path measurement was on `y_quant` or on the N=1,000 fixture;
  the binary path adds SPA, the ER exact test, Firth on 10% of markers and the O(N)
  case/control summary gather, and at N = 50,000 whatever strip changed there did
  not net out positive. Profile before optimising anything else (`next.md`).
- **Strip's resident set is 6-8x upstream's on the marker scan.** It is not
  allocation: it is resident pages of the mmap'd 2.6 GB `.bed`. With
  `SAIGE_PLINK_NO_MMAP=1` the same scan peaks at 312 MB with byte-identical output.
  The OS can drop those pages under pressure; it is a documentation item.

### Strip only, `--nThreads=4` (upstream's `--nThreads` was a silent no-op, LEDGER #30)

Frozen models, 2 reps, same fixture, wall clock (`real`), against the 3-rep 1-thread mean
above. Every 4-thread output is byte-identical to its 1-thread twin (`cmp`, all 8 runs,
region file and `.singleAssoc.txt` included), and each log carries the "Marker scan
will use 4 in-process threads" / "SKAT-O post-processing will use 4 workers" line, so
the knob arrived (the LEDGER #30 shape is excluded).

| workload | 1 thread real s | 4 threads real s | wall speedup | user CPU 1 -> 4 thr | maxRSS 1 -> 4 thr |
|---|---|---|---|---|---|
| single-variant, `y_binary` | 269.6 | 114.5 | **2.35x** | 252 -> 372 s | 2342 -> 2779 MB |
| single-variant, `y_quant` | 66.3 | 24.4 | **2.72x** | 62 -> 79 s | 2789 -> 2793 MB |
| gene tests, `y_binary` (threads go to the SKAT-O post-processing) | 198.1 | 163.2 | 1.21x | 168 -> 198 s | 1072 -> 970 MB |
| gene tests, `y_quant` | 172.5 | 108.5 | 1.59x | 147 -> 160 s | 1006 -> 986 MB |

Against upstream on wall clock that makes the 4-thread strip scan 2.2x (binary) and
3.5x (quant) faster than upstream's single-threaded one, and the gene tests 1.6x /
2.2x. The binary scan's user CPU rises 48% at 4 threads (252 -> 372 s) where the
quantitative one rises 27%: the extra is contention inside the binary per-marker
kernels, another pointer at the binary path.

## 2. P-values: -log10 p, upstream (x) against strip (y)

Plots in `plots/` (one panel per trait; `*_own.png` = each arm reads the step-1 model
its own step 1 fitted, i.e. the pipeline as a user runs it; `*_frozen.png` = both arms
read the fixture's committed models, isolating step 2). Per-row tables in
`harness/work/upstream_compare/cmp/` (the gene-level ones gzipped in `cmp/`; the 70k-row
single-variant ones are regenerable and not committed), the printed summary in
`cmp_summary.txt`. Every departing row carries a tag naming its cause; the tagger's
reconciliations are exact (a chi-square rescale by the two arms' reported variances,
or a match to the two arms' variance-ratio rows), not bands.

**Prediction written before the run** (plan file): no differences except LEDGER #29
(markers at exactly `--minMAC` that upstream drops), #28 (upstream's mis-centred
variance on `AF_Allele2 > 0.5`), #65 rounding on quantitative models, the #50/#59
single-in-group writer fixes, and the step-1 decode fix #53 moving variance ratios.
**Observed:** all of those, plus one thing the prediction called rounding-level and
was not -- LEDGER #97, found here.

### What moves and why (own models, all five traits)

**This table is the 2026-09-05 state, taken before LEDGER #97 was fixed.** All three
quantitative traits have since been re-run against models refitted after the fix --
`y_quant` on 2026-09-07 morning, `y_quant_pc` and `y_null` later the same day -- so the
four rows marked *(pre-#97)* below no longer describe ANY trait's own-model comparison.
They are kept as the record of what #97 did. The post-fix numbers are in the two
"After the #97 fix" subsections at the end of this section. The `frozen` tables are
unaffected either way: both arms read the fixture's committed models there.

| where | rows | cause | how it was pinned |
|---|---|---|---|
| single, 8 markers only in strip, every trait | 8 / 70,709 | **#29** upstream drops MAC exactly 20 | all 8 have MAC = 20 = `--minMAC` |
| single, `AF_Allele2 > 0.5` | 11,671 / 70,701 per trait, every one differs | **#28** upstream's variance on that path is 2-7x too large, marker by marker | hand recomputation from the raw `.bed` + null model reproduces strip's `var` = Σ(g-ḡ)²·VR exactly on BOTH allele sides; upstream's does not. Negative control: on `y_null` upstream's AF>0.5 p-values have median 0.76 and 0.7% below 0.05 (expect 0.50 / 5%); strip's are 0.54 / 4.0%, matching its own AF≤0.5 side band by band |
| single, `AF ≤ 0.5`, binary traits | ~59,000 per trait, all `|Δlog10 p| < 0.05` | **#53** (step-1 decode): upstream's category-1 variance ratio is 0.14% off, category-2 0.005% | the ratio of the arms' reported `var` equals the ratio of their `null_noXadj` rows to 7 digits (0.99864 / 1.00022); a chi-square rescale by it predicts strip's p to 1e-3 in log10 on every row |
| single, `AF ≤ 0.5`, quantitative traits  *(pre-#97)* | ~59,000 per trait | **#65** standardisation moved the `null` rows 0.2-0.6% (rounding on theta 1e-5, coefficients 1e-4, plus the absolute PCG stop #6 now acting in standardised units `[ASSUMPTION: consistent, not isolated]`), and **#97** on the ~35 MAC-20 markers that escalate to the sparse re-test (`var` quotient exactly 1.11126 = 1/s²) | same reconciliation as above, every row |
| gene tests, binary traits | max `|Δlog10 P|` 0.018 | #53 variance-ratio shift | `SE_Burden²` quotient = variance-ratio quotient |
| gene tests, quantitative traits  *(pre-#97)* | up to 0.61 in -log10; 15% of cells | **#97**: cells whose markers use the `sparse` row have `SE_Burden²` strip/upstream = 0.8999 = s² exactly (459 cells for `y_quant`); the rest 1/1.0053 (the `null` row) or in between (mixed MAC categories) | bimodal quotient table; all cells inside [min, max] of the row quotients |
| gene tests, cells with a SKAT-O p of exactly 1  *(pre-#97)* | 1-4 per quantitative trait (e.g. GENE0119, one rare marker, 0.86 → 1.0) | **#86**: SKAT's own saturation at 1, reached once #97 shifts the cell | flagged `SKATO=1` |
| gene tests, `Cauchy` rows  *(pre-#97)* | 7-10 per trait at exactly 1.0 on strip | **#86** (documented decision): SAIGE's CCT falls back to Bonferroni when any cell's p is exactly 1, and #97 pushed cells to 1 | flagged `Cauchy=1` |
| single-in-group (`.singleAssoc.txt`), binary | 965 rows upstream lacks, in exactly 22 regions | **#50**: upstream loses an escalated region's single-variant rows through a temp file the next region truncates | upstream's log has exactly 22 "Non-fast test is performed" lines; the 22 regions coincide |
| one row | `20:17775656:A:T`, `y_quant_pc`, p 0.037 → 0.050 | fast-test cutoff straddle: strip's fast p is 0.0500066 (not re-tested), upstream's was below 0.05 (re-tested on the sparse path) | the only row on 5 traits x 3 tables the tagger could not reconcile; d = 0.13 at p ≈ 0.05 |

**Positive control**: the 21 causal genes rank identically in both arms on both binary
traits (16 in the top 21) and within two places on the quantitative ones (5 and 11 in
the top 21; `y_null` 2-3, chance). **Negative control** (`y_null`): both arms are
near-uniform on AF ≤ 0.5 (KS D 0.014-0.015, 4.4% below 0.05); the AF > 0.5 split is
above.

### After the #97 fix: `y_quant` own models re-run (2026-09-07)

[Measured. `MODELS=own run_arm.sh strip {s2single,s2group} y_quant 1` against the
models refitted after the #97 fix, then `compare_all.sh 1 1`; strip HEAD `30832166`,
upstream arm untouched. Logs `harness/work/region_profile/overnight/3*`; the summary
before and after is snapshotted there beside them.]

The `sparse` variance-ratio quotient strip/upstream is now **1.00000** where it was
**1.11126**; the other four rows are unmoved (`1.00195 1.00199 … 1.00530 1.00553`).
Both `sparse` categories move, and nothing else in `cmp_summary.txt` changes -- the
re-run's diff is five hunks, all inside the three `y_quant models=own` blocks, so every
`frozen` table and every other trait reproduced its previous text exactly.

| table | before | after |
|---|---|---|
| gene, `Pvalue` | Pearson 0.998997, Spearman 0.997843, max\|d\| 0.389 | **0.999991 / 0.999980, max\|d\| 0.027** |
| gene, `Pvalue_Burden` | 0.999224, max\|d\| 0.235 | **0.999995, max\|d\| 0.051** |
| gene, `Pvalue_SKAT` | 0.999116 | **0.999933** |
| gene, cells matching a VR row / lying between rows | 1630 / 1452 of 3115 | **2547 / 552** |
| gene tags | `VR` 1618, `VR-mix` 1330, `Cauchy` 151, `Cauchy=1` 7, `SKATO=1` 1, `small` 23 | **`VR` 1987, `VR-mix` 495, `tiny` 471, `Cauchy` 158, `small` 19** |
| single-in-group, `p.value` | Pearson 0.998771, identical 0/8669, max\|d\| 0.309 | **0.999996, identical 244/8669, max\|d\| 0.022** |
| single-in-group, `var` quotient at MAC ≤ 20.5 | median 1.11126, range [1.1112519, 1.1112676] | **median 1.0000000, range [1.0000000, 1.0000089]** |
| single-in-group tags | `VR` 8630, `AF>0.5` 39 | **`VR` 5908, `tiny` 2479, `none` 244, `AF>0.5` 38** |
| single, `var` quotient at MAC ≤ 20.5 | range [1.0021371, **1.1112647**] | **[1.0000000, 1.0021427]** |

What this says, and what it does not. The prediction in `next.md` was that the own
`gene y_quant` tags "lose their `VR`/`VR-mix` mass". Directionally that is what
happened -- `VR-mix` fell 1330 → 495 (-63%) and 471 cells became `tiny`
(rounding-level) -- but the mass did not go to zero: it largely **resolved** from
`VR-mix` into clean `VR` (1618 → 1987), which is the sharper statement. `VR-mix` means
a cell's quotient lay *between* the two arms' variance-ratio rows because its markers
straddle MAC categories; with the `sparse` row no longer displaced by 1/s², cells that
had been a mixture of a wrong row and a right one now match a single row. The count
that says this directly is the reconciliation line: cells matching a VR row 1630 →
2547, cells between rows 1452 → 552.

The residual `VR` mass is a **different** discrepancy and is untouched by #97: it is
the `null` / `null_noXadj` row quotient 1.00530 / 1.00553, i.e. #65/#53 rounding on the
step-1 fit, which is why the single-variant table barely moves (Pearson 0.973418 →
0.973446) while the gene and single-in-group tables move by two orders of magnitude in
max\|d\|. #97 lived in the `sparse` row, and the `sparse` row is what the low-MAC
gene-test path uses.

Two tags disappeared entirely: `Cauchy=1` (7) and `SKATO=1` (1). Those were LEDGER #86
-- SKAT's saturation at p = 1 and SAIGE's CCT Bonferroni fallback -- reached *because*
#97 had displaced the cells. With the displacement gone, no cell reaches p = 1, so #86
is no longer triggered on `y_quant`. #86 itself is unchanged; only the input that
provoked it is.

No row became `UNPREDICTED`. The one `UNPREDICTED` in `cmp_summary.txt` is the
`20:17775656:A:T` fast-test straddle on `y_quant_pc` documented above, in a trait this
re-run did not touch.

### After the #97 fix: `y_quant_pc` and `y_null` own models re-run (2026-09-07)

[Measured. These two were the last own models still predating the fix. Full chain per
trait -- `MODELS=own run_arm.sh strip {step1,s2single,s2group} <trait> 1`, then
`compare_all.sh 1 1`; strip HEAD `96f16942`, upstream arm untouched, so the `main`
models and outputs are the same 2026-09-05 files the earlier tables used. The strip
build carries the LEDGER #102 fix that the `y_quant` re-run did not; #102 deleted a
dead value and moved no golden, so the two re-runs are output-equivalent builds.]

**Negative control on the procedure.** The comparison summary was snapshotted before
and after. Of its 22 blocks, exactly the **six** belonging to the two re-run traits
changed; every `frozen` block and every block of the three traits not re-run
(`y_binary`, `y_quant`, `y_binary_pc`) reproduced its previous text byte-for-byte.

**Prediction, recorded before the runs, and the outcome.** The `sparse` variance-ratio
quotient strip/upstream was predicted to go to 1 from the 1/s² factors LEDGER #97
measured independently, and it did, in both MAC categories, with the `null` rows
unmoved:

| trait | `sparse` before | `sparse` after | `null` / `null_noXadj` (unmoved) |
|---|---|---|---|
| `y_quant_pc` | 1.11979 | **1.00000** | 1.00030 / 1.00032, 1.00590 / 1.00614 |
| `y_null` | 1.11169 | **1.00000** | 1.00198 / 1.00202, 1.00515 / 1.00538 |

| table | `y_quant_pc` before → after | `y_null` before → after |
|---|---|---|
| gene, `Pvalue` | 0.999716 → **0.999995**, max\|d\| 0.609 → **0.180** | 0.999143 → **0.999966**, max\|d\| 0.218 → **0.116** |
| gene, `Pvalue_Burden` | 0.999226 → **0.999996**, max\|d\| 0.609 → **0.078** | 0.999249 → **0.999997**, max\|d\| 0.201 → **0.014** |
| gene, `Pvalue_SKAT` | 0.999749 → **0.999998** | 0.999116 → **0.999990** |
| gene, cells matching a VR row / between rows | 1626 / 1455 → **2421 / 677** (of 3111) | 1625 / 1466 → **2617 / 490** (of 3120) |
| gene tags | `VR` 1615 → **1895**, `VR-mix` 1360 → **619**, `tiny` 0 → **468**, `Cauchy=1` 1 → **0**, `SKATO=1` 4 → **0** | `VR` 1610 → **2096**, `VR-mix` 1351 → **441**, `tiny` 0 → **431**, `SKATO=1` 1 → **0**, `Cauchy=1` 4 → **2** |
| single-in-group, `p.value` | 0.998439 → **0.999996**, max\|d\| 0.923 → **0.019** | 0.998715 → **0.999997**, max\|d\| 0.437 → **0.015** |
| single-in-group, `var` quotient at MAC ≤ 20.5 | median 1.1197910 → **1.0000000** | median 1.1116919 → **1.0000000** |
| single, `p.value` (barely moves) | 0.997830 → 0.997831 | 0.953439 → 0.953528 |
| single, `var` quotient at MAC ≤ 20.5 | range [1.0003262, **1.1197954**] → [**1.0000000**, 1.0003266] | [1.0020140, **1.1116958**] → [**1.0000000**, 1.0020157] |

Same shape as `y_quant`, including the part that is easy to state wrongly: the
`VR`/`VR-mix` mass **resolves** rather than vanishes. `VR-mix` falls 54% and 67% while
`VR` RISES on both traits, because a cell tagged `VR-mix` had a quotient lying between
the two arms' variance-ratio rows; with the `sparse` row no longer displaced by 1/s²
those cells match a single row cleanly. The reconciliation counts say it directly.
The single-variant tables barely move, for the same reason as `y_quant`: their residual
is the `null` row quotient (#65/#53), which #97 never touched.

`SKATO=1` goes to zero on both and `Cauchy=1` to zero on `y_quant_pc` and 4 → 2 on
`y_null` -- LEDGER #86 was reached *because* #97 had displaced those cells, and #86
itself is unchanged.

**The one `UNPREDICTED` row survives.** `20:17775656:A:T` on `y_quant_pc` is still the
only cell the tagger cannot reconcile, before and after: it is a fast-test cutoff
straddle (strip's fast p is 0.0500066 and is not re-tested), not a variance-ratio
effect, so #97 was never expected to move it. No new `UNPREDICTED` row appeared on
either trait.

### Frozen models (step 2 isolated; `y_binary`, `y_quant`)

Both arms read the fixture's committed step-1 models (the strip binary models are
bit-identical to what HEAD refits; the quantitative ones predate #65). This isolates
step 2, and the result is the 09-01 picture exactly:

| table | rows | identical | differing | cause |
|---|---|---|---|---|
| single, `y_binary` | 70,701 shared + 8 strip-only | 59,310 (every AF ≤ 0.5 row) | 11,391, all `AF_Allele2 > 0.5` | #28 (+ #29 for the 8) |
| single, `y_quant` | 70,701 shared + 8 strip-only | 59,386 (every AF ≤ 0.5 row) | 11,315, all `AF_Allele2 > 0.5` | #28 (+ #29) |
| gene tests, `y_binary` | 3,338 cells | **3,338** | 0 | -- |
| gene tests, `y_quant` | 3,338 cells | **3,338** | 0 | -- |
| single-in-group, `y_binary` | 8,084 shared | **8,084** | 0; upstream lacks 965 rows | #50 |
| single-in-group, `y_quant` | 8,744 shared | **8,744** | 0; upstream lacks 305 rows | #50 |

So with the same null model, strip's step 2 reproduces upstream's gene tests bit for
bit on 6,676 cells, and its single-variant tests everywhere except where upstream is
wrong (#28) or drops rows (#29, #50). Everything in the own-model tables above that is
not #28/#29/#50 therefore comes from step 1: #53 (both trait types), #65 and #97
(quantitative).

**Flip invariance of strip's AF > 0.5 path, checked directly**: the same 667-marker
window of `y_null` run with `--AlleleOrder=alt-first` and `ref-first` gives
`AF_Allele2` summing to 1, `BETA` sign-flipped, and `var` and `p.value` bit-identical
on all 667 markers, 124 of them AF > 0.5.

### LEDGER #97, found by this comparison

Refitting the quantitative null models at HEAD against the fixture's pre-#65 models
moved the `sparse` variance-ratio rows by exactly 1/s² (1.11126, 1.11979, 1.11169 for
`y_quant`, `y_quant_pc`, `y_null`; s = sd of the covariate-adjusted phenotype, computed
independently) while the `null` rows agreed to 1e-6. Both the numerator and the
denominator of the `sparse` ratio go through Σ⁻¹ (`R/SAIGE_fitGLMM_fast.R:2578`,
`:2598`), so it is scale-free and must not be multiplied by `varRatioScale`
(`:2662`, `:2675`). Upstream's `sparse` rows are 0.99999 (correct, ≈ 1 by construction
when the same sparse GRM fits the null); strip's are 1.11. Ungated: no quantitative
golden has `sparse` rows. Recorded as LEDGER #97 (OPEN); fix and test named there and
in `next.md`. Binary traits are untouched (s = 1).

### Upstream defects this comparison re-confirms at N = 50,000

- #28 (AF > 0.5 variance) is not a rounding issue: on `y_null` upstream reports 0.7%
  of AF > 0.5 markers below p = 0.05 where 5% is expected, and on `y_quant_pc` a
  marker at AF 0.958 goes from p = 0.25 (upstream) to 4e-14 (strip). Same allele order
  in both arms; `AF_Allele2` identical everywhere.
- #29 drops exactly the MAC = `--minMAC` markers (8 here).
- #50 silently loses single-in-group rows for every escalated region (965 rows / 22
  regions here).
- #53's uninitialised byte biases upstream's variance ratios; the effect on p-values
  is `|Δlog10 p| < 0.05` on this fixture.

## Files

- `run_arm.sh`, `bench.sh`, `compare_all.sh`, `analyse.R`, `compare_pvalues.R`,
  `plot_pvalues.R`, `README.md` -- the drivers; `results.tsv` -- every timed run;
  `cmp_summary.txt` -- the comparator's printed summaries; `plots/*.png`; `cmp/*.tsv.gz`
  -- per-row comparisons; `logs/` -- every run log (outputs themselves are not
  committed: 7 MB each, regenerable).

## Addendum (2026-09-05, later): the binary single-variant row is not like-for-like

Profiled the same day (`harness/bench/NOTES.md`, "The binary single-variant scan at
N = 50,000"). Firth is 56% of the binary scan's CPU on both arms (33 ms per
Firth-eligible marker, from a per-chunk regression of the logs above, upstream R^2 =
0.96). On AF <= 0.5 the two arms fire Firth on 6,258 and SPA on 3,192 markers,
identically; on AF > 0.5 strip fires Firth on 1,633 / SPA on 730 against upstream's
569 / 276, because #28 deflates upstream's own p-values there. Strip therefore did
~1,000 more Firth fits and ~450 more SPA corrections (~22 s) than the arm it was
timed against; like-for-like the arms are at parity, and after the kernel workspace
change of the same evening strip's 1-thread binary scan is 205 s user against
upstream's 236.
