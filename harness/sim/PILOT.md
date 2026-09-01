# Stage 0 pilot — measured 2026-09-01

Everything below is MEASURED at the final founder count (46,250 EUR), not
extrapolated from the prior N=10,000 run. Two pilots were run because the first
window was unrepresentative.

## Window choice

The first pilot used chr20:35–37 Mb, which has `mean_r = 3.42e-09` — **0.20x the
chr20 mean of 1.679e-08**. Site counts were unaffected but runtime was 5.2x
optimistic, and low recombination means longer haplotypes, fewer independent
markers and more cross-gene LD leakage. Re-piloted at chr20:14–16 Mb
(`mean_r = 1.70e-08`, 1.01x).

**Chosen region: chr20:34,000,000–59,000,000 (25 Mb), `mean_r = 1.638e-08`,
0.98x the chromosome mean, entirely on the q arm (centromere is ~26.4–30 Mb).**

## Measured, 2 Mb at 46,250 founders (92,500 haplotypes)

| | low-r window | representative window |
|---|---|---|
| `mean_r` | 3.42e-09 (0.20x) | 1.70e-08 (1.01x) |
| wall clock | 7.3 s | **38.0 s** |
| peak RSS | 331 MB | **416 MB** |
| trees | 23,496 | 90,054 |
| raw sites | 167,017 | 166,982 |

**Site count is insensitive to recombination; runtime is not.** ~19 s/Mb and
~208 MB RSS/Mb at the representative rate.

This corrects the plan's "1–2 h, the dominant cost and the top risk" estimate by
an order of magnitude: 25 Mb is **~8 min**. Simulation is not the bottleneck and
should be dropped from the risk list.

> ### CORRECTION, measured on the real run (2026-09-01)
>
> **The ~8 min figure was wrong, and the error is in the model, not the
> measurement.** Extrapolating s/Mb linearly assumes coalescent cost is linear in
> sequence length. It is not: the number of ancestral segments carried by the
> Hudson algorithm grows faster than L, so per-Mb cost rises with L.
>
> | at 25 Mb, 46,250 founders | predicted here | measured |
> |---|---|---|
> | wall clock | 7.9 min | **33.2 min** (1990.3 s) |
> | peak RSS | 5.2 GB (208 MB/Mb) | **2.07 GB** (82.6 MB/Mb) |
> | raw sites | 83,491/Mb | **83,211/Mb** |
>
> Realised runtime exponent `ln(52.4)/ln(12.5) = 1.57`, i.e. cost ~ L^1.57 over
> this range. Site counts were accurate to 0.3% and memory was 2.5x
> conservative; **only the runtime model was wrong**, so the sizing table above
> is still sound and only its `sim` column is not.
>
> Realised: 1,159,476 trees, 2,080,270 sites, 6,155,502 edges, 349 MB `.trees`.
>
> The general lesson, twice over now: this stage has been mis-estimated in both
> directions (5x optimistic from an unrepresentative window, then 4.2x optimistic
> from a linear extrapolation). Do not quote a runtime for a length that has not
> been run.

## Site-frequency spectrum (per Mb, founders)

```
raw sites in ts                                83,491/Mb
dropped multiallelic                            2,431/Mb
KEPT (biallelic, polymorphic in founders)      72,178/Mb

  MAF <= 1e-4  (== MAC <= 10, the collapse set)  58,321/Mb   80.80%
  1e-4 - 1e-3                                     9,708/Mb   13.45%
  1e-3 - 1e-2                                     1,981/Mb    2.74%
  1e-2 - 0.05                                       678/Mb    0.94%
  MAF >= 0.05  (clean GRM pool)                   1,490/Mb    2.06%

MAF >= 0.01 (createSparseGRM default floor)      2,169/Mb
VR bin MAC 11-20                                 4,326/Mb
VR bin MAC > 20                                  8,670/Mb
```

**80.8% of variants are at MAF <= 1e-4**, i.e. inside the ultra-rare collapse
set. This is why amendment A1 matters so much: the `max_MAF = 1e-4` stratum is
not a small corner of the data, it is four fifths of it.

## Derived sizing (M capped at 200,000 by the 2.5 GB .bed budget)

Retain a site iff: inside a gene (CDS), OR MAF >= 0.01, OR with probability
`p_thin`. Per Mb: coding 1,031, common-intergenic 2,138, rare-intergenic 69,009.

| L (Mb) | genes | fixed | p_thin | M | .bed | sim | GRM pool |
|---|---|---|---|---|---|---|---|
| 12 | 100 | 38,028 | 0.196 | 200,000 | 2.50 GB | 3.8 min | 17,880 |
| 16 | 134 | 50,704 | 0.135 | 200,000 | 2.50 GB | 5.1 min | 23,840 |
| **25** | **210** | **79,225** | **0.070** | **200,000** | **2.50 GB** | **7.9 min** | **37,250** |
| 40 | 336 | 126,761 | 0.027 | 200,000 | 2.50 GB | 12.7 min | 59,600 |

Pre-registered rule: `p_thin < 0.05` means shorten the region rather than thin
harder. **L = 25 Mb chosen**: 210 genes (vs 168 at 20 Mb) for the same disk and
+1.6 min, and null-gene count is the binding constraint on calibration power.

Checks at L = 25: VR bin MAC 11-20 retains `4,326 x 25 x 0.070` = 7,570 markers
(need >= 2,000); GRM pool at MAF >= 0.05 is 37,250 markers (need >= 3,000).
Both clear.

## The GRM finding is confirmed, and the scaffold is necessary

25 Mb at `mean_r = 1.638e-08` is a genetic length of **0.41 Morgans = 41 cM**, so
however many common markers the region contains, the number of *effectively
independent* ones is a few hundred. That is an upper bound set by genetics, not
by marker count, and it is why `createSparseGRM` cannot be given a clean input
from one region — confirming the synthetic-scaffold decision.

The prior run's diagnostic (`prior_run_evidence/grm_check_0.05.txt`) shows this
directly: **true pairs were estimated correctly (median 0.493, expected 0.5)**
while "unrelated" pairs sat at median 0.075 with the diagonal ranging
0.228–2.230. So the estimator is fine and the *input* is the problem.

Population structure was ruled out: `02_build_dataset.py:74` keeps EUR
individuals only and drops the AFR founders. **[ASSUMPTION, untested]** the
remaining likely cause is that the GRM is dominated by rare markers — at
`minMAFforGRM = 0.01` a carrier of a 1%-MAF variant contributes
`(2-0.02)^2/(2*0.01*0.99) ~ 198` to its diagonal, so a handful of markers
dominate. Testable by re-running the GRM at `minMAFforGRM = 0.05`; the synthetic
scaffold (founder frequencies drawn on [0.05, 0.45]) sidesteps it by construction.
