# Prediction for the LEDGER #132 fix, written BEFORE the gate

Date: 2026-09-13, at `a2962701` (the #134 fix). Read off the COMMITTED bytes, not
from a run of the patched tree.

## The change

The ultra-rare pseudo-marker's **ACAT-V combination weight** becomes the MEAN of
`w_r` over the markers actually folded into it in that cell (non-zero weights
only), for CUSTOM and TRANSFORM sets. Beta and No_Weights sets are untouched --
beta sets already weight the pseudo-marker by `dbeta` at its own collapsed MAF.

**Applied at the ACAT-V call site only.** `AnnoWeightsList` keeps the
pseudo-marker at 1 and continues to drive `wStatVec` and the weighted variance
sub-matrix; a second vector, `AnnoWeightsACATVList`, carries the mean and is used
only by `get_ACATV_pvalue`. The first attempt changed the shared vector and
double-counted `w_r` in the burden -- see "What the first attempt got wrong".

## Which goldens move

Criterion: a golden with ACAT columns (K > 1 or `--is_output_weight_column`),
carrying a CUSTOM or TRANSFORM set, with at least one cell holding BOTH a
collapsed pseudo-marker and a surviving rare marker (`Number_rare > 0` AND
`Number_ultra_rare > 0`). With no surviving rare marker the cell's ACAT-V is the
pseudo-marker's own p whatever weight it carries, so those cells cannot move.

| golden | cells with both | sets | expected |
|---|---|---|---|
| `plink_group_multiweight` | 12 of 26 | CADD, REVEL | MOVES |
| `pgen_group_multiweight` | 12 | CADD, REVEL | MOVES |
| `plink_group_cond_multiweight` | 12 | CADD, REVEL | MOVES |
| `pgen_group_cond_multiweight` | 12 | CADD, REVEL | MOVES |
| `plink_group_multiweight_quant` | 12 | CADD, REVEL | MOVES |
| `plink_group_flexrv_skato` | 8 of 10 | 4 transform sets | MOVES |
| `plink/pgen_group_composed` | 18 | Beta_1_25, Beta_1_1, No_Weights | **DOES NOT MOVE** |

**Six movers.** `composed` is the control: it has ACAT columns and 18 qualifying
cells, and is excluded purely by set kind -- the fix tests
`ws_kind[r] %in% c("custom", "transform")`. If `composed` moves, the set-kind
guard is wrong.

## Which COLUMNS move, and this is the sharper prediction

Only **`Pvalue_ACATV` and `Pvalue_ACATO`**. `AnnoWeightsList` is untouched, so
`Pvalue`, `Pvalue_Burden`, `Pvalue_SKAT`, `BETA_Burden`, `SE_Burden`, `MAC`,
`MAC_case`, `MAC_control`, `Number_rare` and `Number_ultra_rare` must be
BYTE-IDENTICAL in every mover. `Pvalue_ACATO` moves only because it is the Cauchy
combination of (ACAT-V, Burden, SKAT) and takes ACAT-V as an input.

**No `.singleAssoc.txt` or `.singleAssoc.txt_temp` file moves at all** -- they are
single-variant output with no ACAT columns and no dependence on the combination
weights. #134 moved six of them; #132 must move none.

That column-level split is the strongest available check that the fix is confined
to the combination and has not leaked into the statistics.

## Contrast with #134's blast radius, since they overlap but are not the same

- `plink_group_cond_customweight` MOVED under #134 and must **NOT** move here: it
  is K = 1 without `--is_output_weight_column`, so it has no ACAT columns.
- `plink_group_flexrv_skato` did NOT move under #134 (all its cells sit above the
  20.5 variance-ratio boundary) and MUST move here (it has transform sets and 8
  qualifying cells).

## What the first attempt got wrong, recorded so it is not repeated

The mean was first written into `AnnoWeightsList` itself. That vector is not only
ACAT-V's: `wStatVecList <- lapply(AnnoWeightsList, function(w) StatVec * w)`
builds the weighted score vectors, and the weighted variance sub-matrix is formed
from the same weights. The pseudo-marker's 1 is CORRECT there, because C++ has
already baked `w_r` into the collapse `max(w_r*G)` (LEDGER #17) -- weighting its
score again counts `w_r` twice.

`tests/testthat/test-step2-multiweight-behaviour.R` caught it on the first run,
failing on `Pvalue`, `Pvalue_Burden`, `Pvalue_SKAT`, `BETA_Burden` and
`SE_Burden` -- i.e. on the statistics, not on anything to do with ACAT-V. The
existing characterization did its job; the two vectors are now separate.

## The characterization inversion

`test-step2-multiweight-behaviour.R` asserted that doubling every weight in the
group file leaves `Pvalue`/`Pvalue_Burden`/`Pvalue_SKAT` byte-identical while
`Pvalue_ACATV` and `Pvalue_ACATO` MOVE (bounded at 20%). The ACAT half is
inverted: both are now invariant. Toleranced at 1e-10 rather than asserted
identical, because the mean is a fresh arithmetic operation on the weights and
ACAT-V is a CCT round trip through `tan`/`atan` on top of it.

## Addendum: the markerList sidecar, and the flag's blast radius

The fix invalidated `test-invariants-goldens.R`'s first-principles recompute of
`Pvalue_ACATV`, which hardcoded the pseudo-marker's weight at 1. The information
the test lost is already an output SAIGE can write, so it was recovered rather
than conceded.

`--is_output_markerList_in_groupTest=TRUE` writes `<out>.markerList.txt`: one row
per (region, annotation, max_MAF) with the surviving RARE rsIDs and the folded
ULTRA-RARE rsIDs (`R/SAIGE_SPATest_Region.R:1497-1556`). It is weight-independent
-- the C++ indicator marks a marker folded whatever its weight -- so the test
applies the per-set `w != 0` filter itself, which is the implementation's own
rule (`src/Main.cpp:3052`) rather than a restatement of its output.

**PREDICTED**: the flag gates only the block at `:1497`, so it adds one file per
run and changes nothing else.

**CHECKED** by a full regeneration: six new `.markerList.txt` files appeared and
`compare_goldens.sh` reported **zero `DIFFERS`** -- the only lines were
`EXTRA in candidate` for the six new sidecars. The prediction held.

The six sidecars are committed, so folded membership is now part of the pinned
record. #118 and #134 both wanted exactly this and neither had it.

### Why the interim test was not good enough to keep

Solving for the pseudo-marker's weight and bounding the root inside the set's own
weight range gives a wide window -- 0.8 to 10 on the CADD set -- so a wrong
membership, a wrong per-marker p-value or a modestly shifted ACAT-V golden still
produces a root in range. It failed only on gross errors. Drilled by planting
against the restored version, over the 12 mixed cells:

| reconstruction | cells passing |
|---|---|
| correct mean of the folded markers' non-zero weights | **12 / 12** |
| the old hardcoded 1 | 0 / 12 |
| correct mean with ONE member dropped | 0 / 12 |

### One earlier claim withdrawn

"The solved weights 6.525 over 12 and 5.247 over 17 are exact means, confirming
the implementation" was CIRCULAR -- both the solved weight and the value it was
solved from came from the same output, so it showed internal consistency only.
The non-circular confirmation is that the mean of the CADD weights of the twelve
markers the SIDECAR names (rs3, rs8, rs9, rs10, rs32, rs34, rs35, rs42, rs44,
rs45, rs46, rs49) is 6.525, computed from the group file without reference to any
ACAT-V output.
