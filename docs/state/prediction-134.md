# Prediction for the LEDGER #134 fix, written BEFORE the gate

Date: 2026-09-13. Written before `harness/build.sh` was run on the change, as
CLAUDE.md requires ("Predict which goldens should move BEFORE the gate, from the
goldens themselves where possible"). Every number below is read off the
COMMITTED bytes at `53a8d471`, not from a run of the patched tree.

## The change

`src/Main.cpp`, the ultra-rare pass: `ctx.useSparseGRM`, `selectVarianceRatio`
and `cm.isER` read `MAC_nw` (the unweighted collapse's allele count, already
computed ~200 lines above and already what the reported group MAC uses since
#118) instead of `MAC` (the weighted collapse sum, whose units are weight x
genotype).

## Which goldens must move, and why

The criterion is a run that (a) composes a CUSTOM or TRANSFORM weight set, so the
collapse is weighted and `MAC_nw != MAC`, AND (b) reaches the block at
`src/Main.cpp:3634`, `if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest)`.
The observable signature of (b) is a committed `.singleAssoc.txt`.

Band crossings are computed at the default `--cateVarRatioMinMACVecExclude=10,20.5`
by comparing each UR-only cell's weighted AC (its pseudo-marker's `AC_Allele2` in
`.singleAssoc.txt`) against its unweighted `MAC` in the region file.

| golden | UR-only cells | crossings of 20.5 | expected |
|---|---|---|---|
| `plink_group_multiweight` | 12 | **6** | MOVES |
| `pgen_group_multiweight` | 12 | **6** | MOVES |
| `plink_group_cond_multiweight` | 12 | **6** | MOVES |
| `pgen_group_cond_multiweight` | 12 | **6** | MOVES |
| `plink_group_multiweight_quant` | 12 | **6** | MOVES |
| `plink_group_cond_customweight` | 6 | **6** | MOVES |
| `plink_group_flexrv_skato` | 8 (all MIXED) | **0** | DOES NOT MOVE |

**Only the 20.5 boundary is live.** `selectVarianceRatio` (`src/SAIGE_test.cpp`)
maps `MAC <= 10` AND `10 < MAC <= 20.5` both to `varRatio(0)`; only `> 20.5`
selects `varRatio(1)`. An earlier draft of this table said "8 crossings" per
multiweight golden by counting crossings of EITHER boundary -- two of those
eight cross only the inert 10 and change nothing. The six that matter are the
six CADD cells.

**`plink_group_flexrv_skato` DOES NOT MOVE, and this is measured, not unknown.**
An earlier draft of this document called it "not predictable from the bytes"
because all eight of its cells are MIXED (`Number_rare` 4 or 6), and told the
reader that byte-identity there would be "a result to explain". **That was a
regression**: LEDGER #134's STATUS cell and `docs/state/next.md` already record
the correction, made earlier the same day. The way in is that a FlexRV transform
set's `unw_unw` row has weights identically 1 (`FlexRV::isUnitRow`), so ITS
pseudo-marker `AC_Allele2` IS the cell's `MAC_nw`, and every other set of that
cell is read against it: GENE1 `MAC_nw` 54 against weighted 47.85 / 50.47 /
44.4398, GENE2 41 against 37.125 / 27.08 / 24.517 -- all eight above 20.5, so
every set picks the same ratio.

### The sidecars, including the two kinds an earlier draft missed

Each mover also moves its `.singleAssoc.txt`. THREE movers additionally have a
committed **`.singleAssoc.txt_temp`**, which is inside the byte gate
(`compare_goldens.sh` excludes only `*.index`): `plink_group_multiweight`,
`pgen_group_multiweight`, `plink_group_multiweight_quant`. An earlier draft
omitted these entirely.

`plink_group_multiweight_quant.txt.pooledTests.txt` is committed and is predicted
**NOT** to move: the #127 sidecar records the rows combined and the rule that
fired, not the p-values themselves.

**Total predicted: 15 files** -- 6 result files, 6 `.singleAssoc.txt`, 3
`.singleAssoc.txt_temp`.

## Which goldens must NOT move, and why

- **Every `--r.corr=1` run WITHOUT `--is_single_in_groupTest=TRUE`**:
  `plink/pgen_group_multiweight_burden`, `plink/pgen_group_flexrv_burden`,
  `plink_group_flexrv_burden_quant`, `plink_group_ur_smallweights`. They never
  enter the `:3634` block -- the BURDEN else-arm derives no `ctx`.
  **The rule is the conjunction, not `--r.corr=1` alone.** An earlier draft said
  "every `--r.corr=1` run", which is FALSE: `plink_group_urflip_burden` and
  `plink_group_urflip_skato` pass `--is_single_in_groupTest=TRUE` alongside
  `--r.corr=1` and DO enter the block (both have a committed `.singleAssoc.txt`,
  which is the tell). They still do not move, for the OTHER reason -- no custom
  or transform weight set -- so the list was right and its stated reason was
  wrong.
- **`plink/pgen_group_composed`**: reaches the block, but composes only beta and
  No_Weights sets, whose collapse is unweighted (`src/Main.cpp:3025`), so
  `MAC_nw == MAC` identically.
- **Every non-weighted group golden** (`*_nofast*`, `*_fastTest*`, `*_burden*`,
  `*_cond_binary`, `*_cond_quant`, `*_loco_*`, `*_urflip_*`, `*_urmixed_*`):
  same reason -- no custom or transform set.
- **Every marker-path and step-0/step-1 golden**: the change is inside the region
  path's ultra-rare pass.

## The ER site is inert and no test may assert it

`cm.isER` is `&&`-gated on `erAllowedVec[r]`, set at `src/Main.cpp:2335`/`:2343`
to `arma::all(t_weight.col(c) == 1.0)` / `FlexRV::isUnitRow`. When that holds the
collapse weight `wUR` is 1.0, so the weighted and unweighted collapses are the
same vector and `MAC_nw == MAC` exactly. The site is changed for consistency
only; an ER-based assertion would be vacuous.

## The non-vacuity assertion for the new characterization test

From the row: **the SAME cell under two weight sets of different SCALE must get
the SAME band.** On `group_new_snpid_multiWeights.txt`, `GENE1 lof` at max MAF
0.001 is five alleles and today gets 38.4 (CADD, top band) against 3.34 (REVEL,
bottom band). After the fix both must select the band for MAC 5. That assertion
is independent of the goldens and cannot pass vacuously.

## The size of the effect, now MEASURED

An earlier draft of this section declined to confirm the claim that the defect
"changes a headline p-value by up to a factor of two under a pure rescale", on
the grounds that the row records a DOUBLING moving no p-value. **The claim was
right and the hedge was wrong.** A factor of 2 is precisely the factor that
cannot show it: it crosses only the INERT boundary at 10. Measured by the owner
over `x0.01 / x1 / x100 / x10,000` on the same fixture, binary trait, collapse
on:

- headline SKAT-O and burden p-values: **17 of 24 rows move at x0.01, 24 of 24 at
  x100, 17 of 24 at x10,000, worst ratio 0.53** -- i.e. very nearly the factor of
  two claimed.
- `Pvalue_ACATV` ranges over 0.54 to 1.46 and `Pvalue_ACATO` over 0.84 to 1.00
  across the four scales. **That half is LEDGER #132, not this row**, and is not
  fixed here.

The #132 mechanism, for whoever takes that row: the pseudo-marker's normalised
share of the ACAT-V combination is `1 / (sum(w) + 1)`, so as the scale falls
ACAT-V tends to the pseudo-marker's own single-variant p, and as it rises it
tends to the Cauchy of the rare markers alone. Every scale is validly calibrated
under the null; what the scale silently chooses is where between those two limits
the statistic sits, which is a POWER choice the user never made. It also explains
why FlexRV transform sets (weights in [0,1] over a handful of markers, so
`sum(w)` small, pseudo-marker dominates) and `dbeta(1,25)` weights (~25 per
marker, pseudo-marker nearly ignored) behave like different methods.

## OUTCOME: what actually moved

Measured after the fact, `compare_goldens.sh` over a full regeneration:
**exactly the 15 files predicted above, and nothing else.** Every
`--r.corr=1` run, both `composed` goldens, `plink_group_flexrv_skato`,
`multiweight_quant`'s `.pooledTests.txt`, and every marker / step-0 / step-1
golden came back byte-identical.

Two misses in the FIRST draft of this document, both caught by an adversarial
review BEFORE the gate ran rather than by the gate: the three
`.singleAssoc.txt_temp` sidecars were omitted, and `flexrv_skato` was wrongly
called unpredictable. Both are corrected above. The lesson worth keeping is that
the review caught them from the committed bytes alone -- the prediction was
checkable without running anything, and was not fully checked before being
written down.
