# Handoff 2026-09-26: the old-vs-new comparison, and a FlexRV signature that needs a repo repro

## What was done

The v4 production results were compared against the 2024-03-26 BRaVa run
(`all-of-us.palmer.BRaVa_v8.*.SAIGE.{variant,gene}.20240326.txt.gz`, 170 files
each under `dunc_persistent_disk/workspaces/bialleliceffectsincomplextraits
acrossbiobanks/final_summary_files{,_variant}/`), at BOTH levels, over the 127
trait x ancestry pairs whose step-1 N and case count match the old run to 0.5%.
Work was done in an AoU notebook; nothing from it is in this repo.

## Result: the two runs agree except for one trait

**Variant level.** Every binary MAC band has `rho = 1.0000`, median
`log10(p_new/p_old)` of 0.0000 and `frac_10x` of 0; quantitative sits at +0.0001
to +0.0005. The join is exact: `MarkerID` is `chr:pos:ref:alt` on both sides,
zero allele swaps, AF agreeing to 2.4e-5 (which is ~3 allele counts at
`1/(2N)`, i.e. the 20-sample difference in the fitted sets, not a mis-join).

**Gene level.** The only movement is the two changes that were shipped on
purpose. Above p = 0.05 the modified-Liu screen replaces Davies (#158), and
below 1e-5 the saddlepoint tail runs; the like-for-like band 1e-5..0.05 holds
~400 discordant cells out of 5.2 million. `eur/Height` is the visible tail case:
`rho = 0.9932`, median shift +0.085, i.e. p_new/p_old ~ 1.21, which is the
saddlepoint being conservative as specified.

**The exception is `amr/BMI`**, which owns 92.8%, 99.9% and 99.9% of the
gene-level discordance in the three bands and is the only pair of 127 with
variant-level `frac_10x` above 0.01. Calibration says the OLD run is the broken
one: gene-level obs/exp against uniform is 13.8x at alpha 0.05, 58x at 0.01,
446x at 1e-3 and 3383x at 1e-4 (pooled Cauchy rows worse), i.e. 69% of all cells
and 94% of genes below 0.05; the new run is 1.25x / 1.49x / 1.87x / 2.76x.
**The mechanism is NOT the 6.91 variance ratio** -- the 2026-09-25 handoff
established that 6.91 was a truncated-GRM artefact and the refitted AMR models
carry null ratios 0.970-1.966, median 0.997. So the shift is measured and
unexplained; relatedness under-correction in an admixed cohort is the candidate.

## Facts about the 2024-03-26 files, for anyone reading them again

- **Every ultra-rare pseudo-marker row appears TWICE**: 326,139 named
  `REGION:anno:maxMAF` plus 326,139 spelled `UR:UR:UR:UR`, equal counts and
  bit-identical AC sums, 27.7% of the file's rows redundant. Upstream has two
  emission sites (`main:src/Main.cpp:1555` sets `markerVec = "UR"`, `:1662` the
  named form); the four-field spelling is written by no SAIGE, so the combining
  step rebuilt `MarkerID` as `CHR:POS:Allele1:Allele2`, and on a UR row all four
  are the literal "UR". Drop `CHR == "UR"` on both sides before any join.
- A real marker gets one row per REGION it falls in, so overlapping genes
  duplicate it; the duplicates are the same test, so deduping is exact.
- **Variant MAC starts at 11, never 1**, on both sides: at or below
  `--MACCutoff_to_CollapseUltraRare` a marker gets no single-variant row.
- The annotation x max-MAF grids match exactly (6 masks x 0.0001/0.001/0.01,
  17.69 vs 17.695 named UR cells per region), so a gene-level merge on
  (Region, Group, mmaf) is not being silently shrunk.

## OPEN: FlexRV sets whose cell holds only the collapsed pseudo-marker

On `afr/BenIntNeo`, FlexRV, three genes of 18,432 sit at p_FlexRV 7.8e-120,
3.9e-174 and 2.3e-256 with a ~114-order gap to the next gene. The trait is fine
(ECDF obs/exp 0.71 / 1.01 / 1.18 / 1.74 at alpha 0.5 / 0.05 / 0.01 / 1e-3, and
its single-variant p-values agree with the 2024 run at rho > 0.99), so this is
not the null model.

**Signature: every extreme set row has `Number_rare = 0`** -- the cell holds
nothing but the collapsed ultra-rare pseudo-marker, which is what the
MAC1/MAC5/MAC25/MAC50 rungs give, those thresholds sitting at or below the
collapse cutoff. The same genes' sets WITH rare markers (`Number_rare` 5-18)
sit at p 0.15-0.97. `BETA_Burden` = 260588 with `SE_Burden` = 7619, ratio 34.2,
on MAC 86 with MAC_case < 20 and MAC_control 79. |Z| = 34 is impossible from
that MAC under any weighting. NOT separation (MAC_case and MAC_control are never
0), and scale invariance holds where it should: `AM_lof_unw` / `_mafbeta2` /
`_mafbeta5` give betas 792 / 807 / 850 and the SAME p to seven digits.

Two documented candidate mechanisms, NEITHER confirmed: the pseudo-marker's
weighting (CLAUDE.md's allele convention records that changing what `dbeta` sees
for the collapsed marker "blew `BETA_Burden` up to `-1.01393e+13`" -- the same
shape), and the parked `SPA_pval` absolute root tolerance from the #143 work,
measured to move log p ~5e-4 under a x7 rescale of `g.sum` where a steep score
transform rescales by ~1e7.

**Consequence: `p_FlexRV` is a Cauchy combination and so is ~the minimum over
the sets, inheriting these directly.** No FlexRV hit list is usable until each
hit's driving set is checked for `Number_rare = 0`. This also bears on the AoU
chr20 validation, which reported FlexRV results at 160 sets.

**Next step, and it needs NO AoU data**: reproduce on the bundled fixture --
`--flexRV_score`, a gene whose cell holds only ultra-rare markers, a steep
transform -- and check whether |Z| explodes and tracks the transform's
steepness. If it reproduces it is a LEDGER row against the FlexRV kernel; if it
does not, the cause is in the AoU data and is narrowed to one cell shape.

## Constraint on further AoU work

A third egress alert fired on 2026-09-26 during a notebook loop reading both
step-2 arms for every pair (~130 MB FlexRV + ~12 MB SKAT-O per pair on 6
workers). That cell renders essentially nothing -- Agg backend, figures saved and
closed -- so it is evidence for the bucket -> notebook-VM channel and against
"displaying output is the egress". Measured directory sizes: old variant 4.95 GB
and old gene 2.06 GB; new arms across all traits and bins, skato 11.17 GB and
flexrv 24.49 GB. The day's reads were ~20 GB. The FlexRV pass needed only the
`Group == "Cauchy"` rows, 0.6% of what it moved, so a dsub pre-filter removes
almost all of it. The cause of the alerts is still NOT established and the
missing measurement is what AoU meters and over what window.

## Loose ends

- The counterfactual analysis ("what would each mask alone have detected at
  0.05/20000") completed 48 of ~150 pairs before the alert; counts and hit
  identities are in `~/qq_counts.csv` and `~/qq_hits.csv` on the notebook VM,
  and the 48 completed are skewed toward the small cohorts by completion order.
  Hit p-values were not persisted, which they should have been.
- "SKAT-O pooled" in that analysis is the Cauchy combination of the cells'
  BURDEN p-values, not of their SKAT-O `Pvalue`; the pooled row carries all
  three (`R/SAIGE_SPATest_Region.R:779-793`). The BRaVa-reported number is
  `Pvalue`, and it was not recorded.
