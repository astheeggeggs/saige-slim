# FlexRV: the transform table, pinned from the paper (step F0)

<!-- Written 2026-09-09. This is F0 of docs/design/flexrv-plan.md: no code.
     It exists so that F2's kernels and F2's closed-form tests are written
     against a source that was read, not against a lossy web fetch. The plan
     file's own transform table was assembled from a fetch that dropped
     characters; where the two differ, THIS file wins and the plan's table is
     superseded (only the cubic-root denominator differs, cosmetically). -->

## Source

`flexrv.pdf` in the repo root at the time of writing (NOT committed; 5.7 MB),
sha256 `200d1138fe6de734801a644da7f185ff386f98527f118104ee0a397d232c3573`.

Schwartzentruber, Fiziev, McRae, Ulirsch, Farh (Illumina), *FlexRV*, bioRxiv
`10.1101/2025.11.04.686320`, version posted 2025-11-05, CC-BY-NC-ND 4.0.
37 pages. All quotations below are from the **Methods**, section "FlexRV" (PDF
p.30-31 of the extracted text) and the **Extended Data Fig. 1** caption (PDF
p.13). Reference implementation: STAAR, not SAIGE; `github.com/Illumina/flexrv`
was still absent when the plan was written.

Text was extracted with `pdftotext -layout`; the figure was read as a rendered
page image. Both routes agree, and the numeric checks in "Reconciliation" below
were run before anything here was written down.

## The two quotations that define the weights

Verbatim, with the paper's own typography (`1/3` exponents were superscripts;
the multiplication dot is the paper's):

> To approximate the relationship between this score and a complex phenotype of
> interest, we applied 16 transformations, f(s), to this score as follows (see
> Extended Data Fig. 1):
> - Unweighted (w = 1)
> - LoF only (wj = 1 if LoF, otherwise wj = 0)
> - Cubic root with each of 8 different transition values
> - Beta distribution functions with each of 6 different shape parameters
>
> Cubic root functions, f, were defined to map input scores, s, on the interval
> [0,1] to weights, w, in [0,1], with a transition based on a "transition"
> value, t, selected from {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9} as follows:
>
>     w = (sign(s-t)*|s-t|^(1/3) + t^(1/3)) / (|1-t|^(1/3) + t^(1/3))
>
> Beta distribution functions, f, also mapped input scores, s, on the interval
> [0,1] to weights, w, in [0,1], based on the beta probability density function,
> dbeta, with fixed parameter a = 1, and parameter b selected from
> {1.5, 2, 3, 5, 10, 20} as follows:
>
>     w = dbeta(1 - |s|, 1, b) / dbeta(0, 1, b)

> We also applied 12 transformations to variant minor allele frequency (MAF):
> - Unweighted (w = 1)
> - Beta distribution functions with 2 different shape parameters
> - Hard threshold functions such that wj = 1 for variants at or below the
>   threshold and wj = 0 otherwise
>
> Beta distribution functions, f, followed the definition above, but we mapped
> input scores s to w as f(MAF / 0.001). This scales w such that the most common
> variants included in this test will have the smallest scores (~0). Parameter b
> was selected from {2, 5}.
>
> Threshold functions, f, map input scores, s, to weights w based on MAF
> thresholds, where the MAF thresholds were chosen to correspond to minor allele
> counts (MAC) in {1, 3, 5, 10, 25, 50, 100, 200, 500}.

> After obtaining a weight vector from each transformation based on the
> pathogenicity score or annotation, ws, and the transformation based on MAF,
> wMAF, we obtained a final w = wswMAF by taking the product of the two scores,
> resulting in a set of 192 weight vectors.

16 x 12 = 192 confirms both counts (1 + 1 + 8 + 6 = 16; 1 + 2 + 9 = 12).

## The pinned table

`w(i) = w_score(i) * w_MAF(i)`, evaluated per variant `i` in a gene, where `s`
is the variant's pathogenicity score in [0,1] and MAF/MAC are in-sample after QC.

### Score transforms (16), index order `c = 0..15`

| c | name suffix | closed form | notes |
|---|---|---|---|
| 0 | `unw` | `1` | |
| 1 | `lof` | `1` if the variant's annotation is a LoF label, else `0` | annotation-based, see A2 |
| 2-9 | `cbrt0.2` .. `cbrt0.9` | `(sign(s-t)*abs(s-t)^(1/3) + t^(1/3)) / (abs(1-t)^(1/3) + t^(1/3))`, `t` in {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9} | |
| 10-15 | `sbeta1.5`, `sbeta2`, `sbeta3`, `sbeta5`, `sbeta10`, `sbeta20` | `abs(s)^(b-1)`, `b` in {1.5,2,3,5,10,20} | equals the paper's dbeta ratio, see D1 |

### MAF transforms (12), index order `d = 0..11`

| d | name suffix | closed form | notes |
|---|---|---|---|
| 0 | `unw` | `1` | |
| 1-2 | `mafbeta2`, `mafbeta5` | `(1 - MAF/0.001)^(b-1)` for `MAF < 0.001`, else `0`; `b` in {2,5} | direction resolved in A1 |
| 3-11 | `MAC1`, `MAC3`, `MAC5`, `MAC10`, `MAC25`, `MAC50`, `MAC100`, `MAC200`, `MAC500` | `1` if `MAC <= t` else `0`, `t` in {1,3,5,10,25,50,100,200,500} | in-sample MAC, see A3 |

**Set ordering is score-major, MAF-minor**: the `Kt`-row transform table is
`(c, d)` in row-major order over 16 x 12, so `k = c*12 + d`, and the set name is
`<NAME>_<scoreSuffix>_<mafSuffix>` (e.g. `PAI3D_cbrt0.5_MAC10`, `PAI3D_unw_unw`).
This is normative and must be identical in R and C++, exactly as the Phase 5
`jmr = j*q_maf*q_weight + m*q_weight + r` stride is.

### Derivations (D)

- **D1, score beta.** `dbeta(x, 1, b) = b*(1-x)^(b-1)`. With `x = 1-|s|`:
  `dbeta(1-|s|,1,b) = b*|s|^(b-1)` and `dbeta(0,1,b) = b`, so the ratio is
  `|s|^(b-1)`. `s` is in [0,1] by construction, so `abs()` is defensive only.
  Endpoints: `w(0) = 0`, `w(1) = 1` for every `b > 1`. `b = 2` is the identity
  `w = s`, which is why the Extended Data Fig. 1a `Beta, b=2` panel is a
  straight line -- a free check on the whole derivation.
- **D2, MAF beta.** Same algebra applied to `x = MAF/0.001` directly:
  `dbeta(MAF/0.001, 1, b) / dbeta(0,1,b) = (1 - MAF/0.001)^(b-1)`. `b = 2` is
  linear in MAF; `b = 5` is quartic. Undefined for `MAF > 0.001` (the base is
  negative), which is why the `else 0` clamp is part of the definition and not
  an implementation detail. `MAF = 0.001` exactly gives 0.
- **D3, cubic root.** At `s = 0` the numerator is `-t^(1/3) + t^(1/3) = 0`; at
  `s = 1` it is `(1-t)^(1/3) + t^(1/3)`, which is the denominator, so `w(1) = 1`.
  At `s = t`, `w = t^(1/3) / ((1-t)^(1/3) + t^(1/3))` -- 0.3865 at t=0.2, 0.5 at
  t=0.5, 0.6753 at t=0.9. Monotone increasing on [0,1]; range exactly [0,1].
  `abs(1-t)` in the denominator is the paper's; for `t` in (0,1) it is `1-t`,
  so the plan file's `(1-t)^(1/3)` is equivalent and not a discrepancy.

## Ambiguities in the paper, and how each was resolved (A)

These are the reason F0 exists. Each is a place where transcribing the Methods
literally would produce a different program.

### A1 -- the MAF beta transform's DIRECTION (the one that matters)

The Methods say the MAF beta "followed the definition above, but we mapped input
scores s to w as f(MAF / 0.001)". Substituted literally into the score-beta
formula, `s := MAF/0.001` gives `w = (MAF/0.001)^(b-1)`, which is **increasing**
in MAF: weight ~0 for singletons and 1 at the 0.1% boundary.

That is contradicted by the paper's own next sentence -- "This scales w such
that the most common variants included in this test will have the smallest
scores (~0)" -- and by Extended Data Fig. 1b, whose `Beta, b=2` and `Beta, b=5`
panels both start at weight 1 on the left of the axis and fall to 0 at the right.

**Resolved: the weight is decreasing.** `w = (1 - MAF/0.001)^(b-1)`, i.e. the
beta pdf is evaluated at `MAF/0.001` itself, not at `1 - MAF/0.001`. Two
independent routes reach the same expression: (a) `dbeta(MAF/0.001, 1, b) /
dbeta(0,1,b)`, and (b) the score form with the argument reflected,
`s := 1 - MAF/0.001`. The plan file's table already had this; F0 confirms it
against the figure rather than against intent.

Quantitative check against the figure. Fig. 1b's x-axis is **allele count 1 to
1000 on a log scale**, not MAF, and the caption says so. `MAF/0.001 = MAC/1000`
exactly when `2N = 10^6`, i.e. at the paper's ~423,614 unrelated individuals
(`2N = 847,228`) the two axes coincide to within ~15%; the figure is the MAF
transform drawn in MAC units at their N. At MAC = 100 the pinned form gives
0.900 (b=2) and 0.656 (b=5); the literal substitution gives 0.100 and 1e-4. The
b=5 panel reads ~0.65 at MAC = 100 and is visibly above zero across the whole
decade 10-100, so the literal reading is excluded, not merely disfavoured.

**Implementation consequence:** we key on MAF, not on MAC/1000. That is
N-independent and matches the Methods text; it reproduces the figure only at the
paper's N, which is the correct behaviour and not a discrepancy to chase.

### A2 -- the LoF-only mask: annotation, or `s == 1`?

The Methods define it by annotation (`wj = 1 if LoF, otherwise wj = 0`).
Extended Data Fig. 1a draws it as a step at `s = 1`, because LoF variants are
assigned `s = 1` by construction. These differ for a missense variant whose
within-gene PrimateAI-3D percentile is exactly 1.0, which the paper's own
"taking the max" rule can produce.

**Resolved: annotation-based**, per the Methods and per `--flexRV_lofAnno`
in the plan. The figure is a rendering convenience. Documented because a reader
comparing our `lof_*` rows against the figure will otherwise see a difference
and think one of us is wrong.

### A3 -- hard thresholds: MAF thresholds, or MAC?

The paper applies **MAF thresholds "chosen to correspond to" MAC in
{1,...,500}** at their fixed cohort size. We apply the MAC comparison directly,
on the in-sample MAC after QC, which is what `MACVec.at(i)` holds and what the
plan's whole reason for building the sets inside C++ rests on. The two agree
exactly at the paper's N and differ at any other N -- and the MAC reading is the
one that keeps its meaning across cohorts, which is why the thresholds were
specified in MAC units in the first place. Recorded as a deliberate,
documented difference.

### A4 -- the CCT formula as typeset drops its normalisation

The paper gives `T = sum_i tan(0.5 - p_i)*pi` and `p_FlexRV = P[C(0,1) >= T]`.
As typeset that is (i) missing the parentheses -- it must be
`tan((0.5 - p_i)*pi)` -- and (ii) missing the `1/k` (or general weights summing
to 1) that makes `T` standard-Cauchy under the null; with a bare sum over
k = 192 tests the null is Cauchy scaled by k, and `P[C(0,1) >= T]` would be
badly anti-conservative. Their citation is Liu & Xie, whose CCT is the weighted
average.

**Resolved: use the repo's existing `CCT_cpp` / `cctVec` machinery unchanged**
(equal weights `1/k`, the Phase 5 pooled-Cauchy row). No new p-value combiner is
written for FlexRV. This is a typesetting loss in the preprint, not a
methodological difference.

### A5 -- which AC do the paper's pseudo-variant MAC thresholds read?

For genes with >= 400 variants the paper collapses AC <= 5 variants into up to
50 pseudo-variants, binned by AC (1..5) x score decile. It never says whether
the `MAC <= t` transforms then read the bin's **constituent** AC (1..5) or the
pseudo-variant's **realised** summed AC.

**Read as the constituent AC.** The 5 x 10 structure has no other purpose: if
the AC label were not used downstream, score deciles alone would give 10
pseudo-variants and the same speedup. Under the realised-AC reading, `MAC<=1`,
`<=3`, `<=5` and nearly `<=10` would select nothing in every large gene -- 30+
of the 192 sets empty -- which contradicts their claim that binned and unbinned
burden p-values were "highly similar".

We do not replicate the binning at all (see below), so this only matters if we
ever reconcile against a STAAR FlexRV run. Recorded so that reconciliation
starts from the right reading.

## Differences from the paper that are deliberate and permanent

- **Large-gene collapse.** The paper: "for genes with >= 400 variants we
  collapsed variants having total AC <= 5 into up to 50 pseudo-variants based on
  each AC from 1 to 5, and pathogenicity score in bins of 0.1 ... Each
  pseudo-variant ... was assigned the mean pathogenicity score of variants in
  the bin." This was a STAAR speed workaround ("The association p-values
  obtained from burden tests with binning were highly similar to those without
  binning"). We do NOT replicate it. SAIGE's own ultra-rare collapse
  (`--MACCutoff_to_CollapseUltraRare`, default 10, per set `max(w_r*G)`) plays
  the same role with different arithmetic. Consequence: our per-gene p-values
  for genes with >= 400 variants are not expected to match a STAAR FlexRV run
  variant-for-variant, and no test should assert that they do.
- **Null model and test machinery.** STAAR + SPA (`max_iter = 20`) there;
  SAIGE's step-1 null GLMM, its own SPA cascade and ER there here. The paper's
  own runs used unrelated individuals and no GRM.
- **Score construction is the caller's job.** The paper builds `s` as LoF = 1,
  missense = PrimateAI-3D percentile among all *possible* missense variants in
  the gene, combined by max. SAIGE reads `s` from the group file's
  `score:NAME` line (F1). We define the transforms, not the score.
- **The MAF < 0.1% cohort filter is not a weight.** The paper selects variants
  at MAF < 0.1% before testing (Methods, "Variant annotation": VEP v109.3 most
  severe consequence only, in-frame indels excluded, variants with > 5% missing
  calls excluded, missense and LoF only, "MAF < 0.1% in the full cohort", no
  lower bound -- singletons are kept, which `MAC <= 1` requires). `MAF/0.001` is
  therefore a rescaling of an already-filtered set, not a filter in disguise.
  Three mismatches when mapping it onto SAIGE:
  - **Boundary.** SAIGE's cell is inclusive (`find(maxMAFVec >= pr.MAF)`,
    `src/Main.cpp` ~:2509), so MAF exactly 0.001 is IN the cell where the
    paper's strict `<` excludes it. `mafbeta*` gives it weight 0 either way;
    `unw` and the larger `MAC*` sets do not.
  - **Which MAF.** Theirs is a fixed cohort-level property of the variant list;
    SAIGE recomputes MAF in-sample after QC. The closer analogue is to apply
    the 0.1% filter when BUILDING the group file.
  - **Missingness.** Their 5% exclusion vs SAIGE's `--maxMissing` default of
    0.15. Set it to 0.05 to match.


## What enters the test (and what the transforms cannot change)

The 192 sets are re-weightings of a variant set they do not choose. A weight of
zero can push a variant out; nothing in the transform table can pull one in.
The cascade, in order:

1. the group file's `var` line;
2. the annotation group -- exact `%in%` match on the ANNO row
   (`R/SAIGE_SPATest_Region.R` ~:1944);
3. the `--maxMAF_in_groupTest` cell, inclusive;
4. run-level QC (`--minMAF`, `--minMAC`, `--maxMissing`);
5. **the ultra-rare collapse** -- every marker with
   `MAC <= --MACCutoff_to_CollapseUltraRare` (default 10) leaves as an
   individual and is folded into ONE pseudo-marker per set, genotype
   `max(w_r * G)` (`src/Main.cpp` ~:2598), weight `w0 = 1` unscaled (LEDGER #17);
6. only then `w = w_score * w_MAF` on the surviving real markers.

**Step 5 collides with four of the twelve MAF transforms.** `MAC<=1`, `<=3`,
`<=5` and `<=10` sit at or below the default collapse cutoff, so every variant
they select has already been collapsed: every real marker gets weight 0 and the
whole test is carried by the pseudo-marker. Because the collapse is WEIGHTED and
per set, the masking still happens in the right order and each cutoff gets its
own pseudo-marker (`genoURMat` is `N x q_anno_maf_weight`) -- but the operator is
`max`, an "any qualifying carrier" indicator, not the additive weighted count the
paper's corresponding column reports. Two consequences:

- `w0 = 1` does NOT bias these four sets. When the cutoff is >= the threshold
  there are no real markers, the pseudo-marker IS the test, and a burden
  p-value is invariant to a global rescaling of `w`. The `w0 = 1` asymmetry
  only bites sets that mix a UR bucket with real markers.
- `max` and the paper's additive collapse differ only for individuals carrying
  TWO OR MORE qualifying ultra-rare variants in the same gene. That count is
  computable straight from the genotypes, before any test runs, and it bounds
  the size of the difference a priori. F5 computes it first.

**We do not port the paper's AC x decile binning.** It is a STAAR speed
workaround ("STAAR ran slowly for genes with a large number of variants") whose
own validation is that it approximates the unbinned additive burden -- and the
unbinned test is available exactly, with no new code, via
`--MACCutoff_to_CollapseUltraRare=0` (the collapse is one condition,
`if(!(MAC > g_cfg.region_minMAC_cutoff))`, `src/Main.cpp` ~:2345, and 0 is a
legal value: `R/checkArgs.R` rejects only negatives). Porting it would mean
implementing someone else's approximation to something we can compute. The
recommendation is therefore by trait:
- **quantitative** (the paper's own 62 traits): `--MACCutoff_to_CollapseUltraRare=0`;
- **binary**: keep the collapse, which is SAIGE-GENE+'s variance stabilisation
  and earns its keep on unbalanced case/control, and drop the four now-redundant
  transforms with `--flexRV_MAFTransforms=unw;mafbeta:2,5;mac:25,50,100,200,500`.

If F5 shows the operator matters, the cheap middle option is reviving upstream's
inactive `sum_geno` collapse (the commented block at `R/checkArgs.R` ~:123):
additive instead of `max`, still one pseudo-marker, no shape change. That is a
strictly better target than the paper's binning.

## When FlexRV degenerates

Burden p-values are invariant to a positive global rescaling of `w`, so two
weight sets differ only if their RATIO pattern across markers differs. The
number of distinct score values therefore bounds the number of distinct sets:

- **constant score within a cell** -- every score transform returns the same
  value, all 16 coincide, and the 192 sets collapse to 12. This is the pLoF-only
  case under the paper's own construction (LoF = 1): `cbrt_t(1) = 1` and
  `sbeta_b(1) = 1`, so the score axis is entirely degenerate and what runs is a
  12-way MAC/MAF sweep combined by CCT. Not worthless, but it is not FlexRV and
  the pooled p should not be called `p_FlexRV`.
- **score in {0, 1} exactly** -- 0 and 1 are fixed points of every transform, so
  all 14 nonlinear transforms become the same indicator mask: 16 score sets
  collapse to at most 3 (`unw`, the indicator, `lof`). A two-valued score at
  INTERIOR values does not degenerate: {0.2, 0.8} gives ratio 1:2.4 under
  `cbrt0.5` and 1:2.5e5 under `sbeta20`.
- **a LoF-only annotation group** -- `lof_*` equals `unw_*` there, 12 duplicate
  sets. Except through ER: `erAllowedVec[r]` is keyed on the set KIND, so
  `unw_unw` is ER-exempt and `lof_unw` is not, and on a binary trait with low
  MAC the two identical-weight sets can take different test paths. Keep the
  static rule (a data-dependent ER choice would be worse) and document it.

Duplicates are not invalid -- CCT is robust to arbitrary dependence in the tail
-- but they enter the pool more than once and re-weight it. The escape is the
spec string (`--flexRV_scoreTransforms=unw` collapses the score axis to 1), not
a code path. If duplicate-column memoisation is ever added it must be
COMPUTE-ONCE-AND-COPY (byte-identical output, pool unchanged), not
skip-and-reweight (which changes the row schema); and equality must be tested on
the FULL weight vector including the `w0 = 1` UR entry, since proportionality is
not sufficient once a UR bucket is present.

## Rules this fixes for F2/F3

- **ER gating.** `erAllowedVec[r]` is true only for the `unw_unw` set, whose
  weights are identically 1; every other transform set has non-unit weights and
  follows the custom-set rule. Unchanged from the plan; restated because the
  `lof_unw` set also has weights in {0,1} and is NOT exempt -- a zero weight is
  not a unit weight, and the ER path assumes an unweighted burden.
- **Zero-weight cells are expected, not errors.** Any `MAC<=t` set with `t`
  below a gene's smallest MAC gives an all-zero weight vector: BURDEN omits the
  row and the Cauchy pool skips it; SKAT-O reports NA through the SKAT failure
  path. At K = 192 on a small gene this will be many rows.
- **Memory.** `genoSumMat` is `N x (n_anno * n_maf * K)` and `P1Mat`/`P2Mat` are
  `max(m1, n_anno*n_maf*K) x N`, 8 bytes per double. At K = 192 with one anno
  group and one max-MAF cell: `192 * N * 8` = 768 MB per array at N = 500,000,
  ~3 GB over the four arrays. At the default 12 cells it is 2304 columns,
  9.2 GB per array -- which is why the recommended run is one anno group
  (`lof:missense`), one max MAF (0.001), `--r.corr=1`, and why F3 prints an
  estimate once at region start when `Kt > 0`. Whether UR storage must go sparse
  is F5's measurement to decide.

## Reconciliation (run 2026-09-09, before this file was written)

Each closed form computed two ways in R 4.x, on grids of 11 points:

| quantity | route A | route B | max abs difference |
|---|---|---|---|
| score beta, b in {1.5,2,3,5,10,20} | `dbeta(1-s,1,b)/dbeta(0,1,b)` | `s^(b-1)` | 1.11e-16 |
| MAF beta, b in {2,5} | `dbeta(m/0.001,1,b)/dbeta(0,1,b)` | `(1-m/0.001)^(b-1)` | 1.67e-16 |
| cubic root, t in {0.2,0.5,0.9} | direct | endpoints `w(0)=0`, `w(1)=1`, `w(t)` = 0.3865 / 0.5 / 0.6753 | exact |

These are the values F2's `test-unit-flexrv.R` should pin (fast tier, the
`test-unit-spa.R` shape): `s = 0`, `s = 1`, `s = t`, the `MAF = 0.001`
boundary from both sides, `MAC = t` and `MAC = t+1`, and the 16 x 12 product
ordering `k = c*12 + d`.
