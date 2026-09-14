# Handoff 2026-09-14 (evening): LEDGER #143 FIXED -- the binary conditional region kernel

Two commits, one ledger row, each gated at `both` modes with the golden
movement predicted before the gate. Read `python3 harness/tools/ledger.py 143`
for the authoritative numbers; this is the narrative and the method.

## What the bug was, in one paragraph

A binary-trait group test with `--condition` computes per-conditioning-marker
scale factors through a burden scale-factor kernel. That call site fed the
kernel one input of four in the wrong units (the raw per-marker score against
a weighted Phi, q.sum and g.sum) and ran it through an unsuffixed textual
twin of `get_newPhi_scaleFactor_traitType` that lacked the #104
standardisation. Only the `_cond` columns of a binary SKAT-O/SKAT/Burden
region run are affected; unconditional columns, the marker path, the
quantitative trait, and the C++ BURDEN emitter (`--r.corr=1`) never touch it.
Upstream `main` has both halves.

## How it was decided, not argued

The row's STATUS said the correct `VarQ_2` on this path was an owner's
decision because the #104 identity was not exact there. It turned out to be
decidable by a control rather than by theory:

- **Weight invariance.** Conditioning on a set of markers is a linear
  reparameterisation, so the conditional result cannot depend on the weight
  attached to a conditioning marker (one marker: any non-zero weight; several:
  a uniform rescaling). At `e508a8e4`, conditioning on rs30 alone with
  `--weights_for_condition=1` against `5` moved every `_cond` p-value on GENE1
  by 0.3-6 % and `SE_Burden_cond` by ~2e-5. That is a defect by construction.
- **Located, not observed.** `test-characterization-ledger143.R` has a
  fast-tier block showing the kernel itself IS invariant under a consistent
  rescaling of its four inputs (1e-15 below the |z| < 2 gate). So the units the
  CLI fed it were the only place the dependence could live.
- **Once the units agree, #104's identity holds mechanically** (`sum(Score)`
  = `q.sum - m1`, `VarQ_2` = `var1` before the standardisation, `sum(Phi)`
  after), so there was no second theory to choose. Both predictions were made
  from the row's instrumented P3 numbers and landed to the printed digit:
  `r` 0.4733 -> 0.951743 (= 654.2276 / 687.3994) -> 1.000000.

## What moved

Five goldens, twice each: P3, P15 (pgen mirror), P10 + pgen mirror, P34. Only
their `_cond` columns. Net on P3: conditional p-values by up to 7 % relative,
`BETA_Burden_cond` sign flips on three GENE2 rows at |beta/SE| ~ 0.02 (noise),
`SE_Burden_cond` <= 2.2 %. Byte-identical, as predicted: the BURDEN
conditional pair (P25 + mirror), every `.singleAssoc.txt` (marker-path `_c`
columns), every `.markerList.txt`, and every other golden in the suite.

The prior "`SE_Burden_cond` staying byte-identical" prediction on the row was
right for the standardisation half and wrong for the units half; the STATUS
cell says so.

## Two things learned that outlive the row

1. **`SPA_pval` is not scale-invariant at the 1e-4 level.** Its root-finder
   uses an absolute tolerance (`eps^0.25`) on a root that scales as 1/w, so
   rescaling `g.sum` by 7 moves the log p by 5.6e-4. On the marker path g is a
   genotype and the scale is fixed; on the region path `g.sum` is a weighted
   sum with dbeta weights up to 25. Not #143, not a row yet -- the test holds
   the SPA regime at 1e-3 and says why. Worth a row if anyone needs the region
   burden SPA reproducible across weight scales.
2. **The regeneration recipe** used here (`scratchpad/regen143.sh`, not
   committed): run the predicted movers AND the predicted non-movers twice into
   two trees, `cmp` the pairs, then `cmp` against the committed golden and
   print `REPRO|NONREPRO MOVED|IDENTICAL` per file. One gate run then proves
   the prediction in both directions: a predicted mover that did not move
   shows as IDENTICAL before it is copied, an unpredicted mover fails the gate.

## Next

Unchanged from `next.md`: nine open rows (#144-#152), Phase G deferred, the
FlexRV scaling question.
