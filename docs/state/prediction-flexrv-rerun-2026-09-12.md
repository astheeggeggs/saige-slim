# Prediction: the FlexRV study sweep re-run at HEAD (2026-09-12)

Written 2026-09-12 02:40, **before any run**, against HEAD `38eee5f5` and the
pinned library (`SAIGE.so` md5 `fd02953542c492c2b17c7d4f3cd94f0b`, built from
`6ffff192`'s sources; the two commits since touch no `src/`, `R/` or `extdata/`).

## Why re-run

The 2026-09-10 sweep (144 runs, `harness/sim/flexrv_study/results-2026-09-10.txt`)
ran at `5dc679c3` (dirty), BEFORE LEDGER #118 (per-set MAC), #126 (min-group-MAC
gate on both paths, default 5 -> 1), #127 (pooled-row sidecar), #117 (zero
weight drops the marker), #115, #114 and the `49c1b60d` perf change. Its raw
outputs were deleted in the 2026-09-11 clear-out; only the `analyse.R` summary
survives. Every number quoted from it about FlexRV is therefore a number about
pre-#118 code. The step-1 models are the SAME inputs the sweep used
(`step1_<trait>.rda`, 2026-09-10 19:27, fitted under the #122 guard and before
#124 landed; kept frozen so the comparison is like-for-like).

The grid is the sweep's own: 4 traits x 2 collapse cutoffs x (3 score arms x 4
FlexRV methods + `unw`/`beta`/`skato`/`grid` + synonymous x {flex192, unw}) =
144 runs, `study.sh 4`, `MINLOAD_REFUSE=20` (the driver's default of 2.5 refuses
its own 4-wide parallelism).

## Predictions

**P1 (pipeline).** All 144 runs exit 0 (captured `rc`, not the `$(date)`
trap) with 210 regions each; #128 is unreachable at the default gate. The three
score arms' `flex192` outputs carry 40,386 data rows (210 x 192 set rows minus
the 144 sets with empty support, plus 210 pooled rows), the same at collapse 0
and 10 because support does not depend on the collapse.

**P2 (pipeline, the control).** `skato` and `grid` use beta(1,25) weights on
the SKAT-O path. On that path #118 is the identity (a beta set's ultra-rare
collapse was already unweighted; `dbeta` is never exactly 0), the #126 gate at
1 refuses nothing, #127 writes only a sidecar, #117 finds no zero weight. So
**every `skato` and `grid` line of CONTROL 1, 2 and 3 must reproduce the
2026-09-10 file exactly** (obs, rho, detections, oracle, by-link counts). `unw`
and `beta` are BURDEN-path runs with No_Weights / beta sets: #118 is the
identity there too, but the gate moved 5 -> 1 and admits the two genes whose
`high:medium x 0.001` cell has MAC 3 and 4 (predicted from the allele counts
last night: 1 gene each). So their CONTROL 1 obs may move by at most 2/210 =
0.0095, detections by at most +1, and only on the informative/plain group file
(the synonymous `low` mask has no sub-5 cell).

**P3 (what may move).** The FlexRV methods. Mechanism: the per-set MAC now
drives the variance-ratio band (worth ~0.3% of SE on these models), ER
eligibility (`unw_unw` only) and the sparse-GRM switch at MAC 20.5, which the
F5 measurement put at 2,722 -> 4,287 rows on the `y_null` informative file
(+1,565 of 40,176 set rows take the sparse estimator now); 8.8% of set
p-values moved, median 0.8%, max 20%. Consequence for the SUMMARY: CONTROL 1
lines for flex methods move by a few genes (|delta obs| <= 0.02), CONTROL 2
detections on `y_binary` by at most 2 per line and rho by at most 0.05, and
CONTROL 3's monotone flags may flip where two counts were within 1 of each
other. Direction not predicted: the sparse estimator can go either way per
set.

**P4 (what will NOT change).** The three findings that block the study are
structural, not numerical: the "uninformative" arm keeps the `lof` set and the
12 score-free MAF transforms, so it will still beat `unw`; the `low` class
carries causal variants by construction of the simulation, so the synonymous
mask will still detect causal genes on `y_binary`; `y_quant` will still be
undetectable (0-3 of 63). A re-run cannot resolve them; the owner's decisions
can.

**Falsifies the pipeline:** any `skato`/`grid` line differing at all; `unw`/
`beta` moving by more than 2 genes or a detection changing by more than 1;
any run with non-zero exit or fewer than 210 regions; any flex line moving by
more than 5 genes at alpha 0.05 or more than 3 detections.

---

# Outcome (2026-09-12 03:20, appended after the runs; nothing above was edited)

144 runs in 41 minutes wall at four wide (`harness/work/flexrv_study/study-2026-09-12.log`).
Analysis: `harness/sim/flexrv_study/results-2026-09-12.txt`; the line-by-line
comparison against 2026-09-10 is `compare_results.py` (null control: identical
inputs report nothing; positive control: three planted changes are all reported).

**P1 -- held, except one clause that a new ledger row explains.** 144 of 144
exit 0, 210 regions each, no empty output. The three score arms' `flex192`
outputs carry 40,386 (informative), 40,361 (realistic) and 40,362
(uninformative) data rows at collapse 10, and **40,406 at collapse 0 in every
arm**. "The same at collapse 0 and 10" was wrong: at collapse 10 the dosage-
zeroing filter is applied to the weighted ultra-rare pseudo-marker in weight
units, and a set whose whole support is ultra-rare with weights at or below
0.2 loses its row (LEDGER #129, with a bundled reproduction and a
characterization test). 20 / 45 / 44 / 21 rows per trait for the four arms,
never the reverse.

**P2 -- held exactly.** `skato` and `grid`: **0 of 40 lines changed**, every
obs, rho, detection, oracle and by-link count reproduced to the digit. `unw`
and `beta`: 0 CONTROL 1 lines changed, 0 detections changed; 8 CONTROL 2 lines
moved in rho by 0.01-0.02, through exactly the two genes predicted
(`GENE0039`, MAC 3, and `GENE0167`, MAC 4, whose cells the gate now admits;
both burden runs carry 210 cells where they carried 208).

**P3 -- held, inside every bound, and with one direction.** 56 of 280 FlexRV
lines changed. CONTROL 1: two lines moved by one gene (`y_binary` non-causal
at alpha 0.01, 0.0136 -> 0.0068); max |delta obs| 0.0068 against the bound
0.024. CONTROL 2 on `y_binary`, collapse 10 (collapse 0 in brackets):

| method | 2026-09-10 | HEAD | delta |
|---|---|---|---|
| flex192 informative | 13 | 11 (11) | -2 |
| flex192 realistic | 12 | 10 (11) | -2 (-1) |
| flex192 uninformative | 12 | 11 (10) | -1 (-2) |
| flexscore, all arms | 11 | 10 (10 / 10 / 9) | -1 (-1 / -1 / -2) |
| flexmaf, all arms | 9 | 8 | -1 |
| flex1 = unw = beta | 8 | 8 | 0 |
| grid | 11 (c0: 10) | 11 (10) | 0 |
| skato | 9 (c0: 7) | 9 (7) | 0 |

Max |delta detected| 2 against the bound 3; max |delta rho| 0.03 against 0.05;
the oracle (Bonferroni best set) fell further, 11 -> 7 on flex192. **Every
FlexRV count fell and none rose**: the per-set MAC sends the small-support
sets to the sparse-GRM estimator and the variance-ratio band of their own
count, and on this fixture that is uniformly more conservative. CONTROL 3's
monotone flag flipped at collapse 10 (11 / 10 / 11), counts within one as
allowed for.

**P4 -- held.** Uninformative flex192 still beats `unw` (10-11 against 8); the
synonymous mask still detects 7 of 63 causal genes on `y_binary`; `y_quant`
is still 1-2 of 63 for every method.

**What the like-for-like numbers say now.** At HEAD the full 192-set FlexRV
(11) equals the standard SKAT-O grid (11) on `y_binary` and beats the single
unweighted burden (8) by three; the 2026-09-10 margin over the grid (13
against 11) was pre-#118. The collapse-0 column is the cleaner one, because
#129 acts only through the pseudo-marker: at collapse 0 flex192 reads
11 / 11 / 10 against grid 10 and `unw` 8.

**The one large collapse discrepancy is #129, proved by switching it off.**
The COLLAPSE section's `y_binary uninformative` max |delta log10 p| went
0.477 -> 4.45; that is `GENE0181` (causal), pooled 8.7e-4 at collapse 0 and
3.1e-8 at collapse 10 with all 192 sets present at both. Its `sbeta10` sets
(12 ultra-rare markers, MAC 30, 2 case / 28 control alleles) read 1.7e-3 and
2.9e-9; re-run at collapse 10 with `--dosage_zerod_cutoff=0` they read
1.710212E-03 to the digit and the pooled row 0.000874857 is collapse 0's. The
carriers whose weighted entry was at or below 0.2 had been dropped from the
test while the MAC columns still counted them. 21 of that gene's 192 sets
differ by more than tenfold between the cutoffs with the zeroing on, none by
twofold with it off. So one of the collapse-10 FlexRV detections is a #129
artefact that happens to land on a causal gene; the collapse-0 counts carry
no such term.
