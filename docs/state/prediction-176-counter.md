# Prediction: LEDGER #176's counter and the fixture that drills it (2026-10-01)

Written BEFORE the first run of `group_new_snpid_withWeights_abovecap.txt`, so
the run can falsify it. Method: decode `extdata/input/genotype_100markers.bed`
directly, restrict to the model's 1,000 samples, apply the anno and max-MAF
filters, split at the collapse cutoff of 10, build the vector each site measures,
and evaluate `atomDominance` (`src/SAIGE_test.cpp:196`) on it --
`neff = (sum|v|)^2 / sum v^2`, `material = #{ |v| >= 1e-6 * max|v| }`.

## The control that says the method is right

The same code, run on the EXISTING atom fixture (non-atom weight `1e-7`),
predicts `n_material` **1** for GENE1 and **2** for GENE2 at every max MAF.
That is exactly what the committed `plink_group_flexrv_atom.txt.exactByWeight.txt`
and the P70/P71 sidecars record. So the arithmetic is reproduced against ground
truth before it is used to design anything.

It is also the reason the fixture is needed at all: `1e-7` is BELOW
`kExactByWeightTieEps` (1e-6), so every non-atom carrier is immaterial and
P69-P71 can only ever exercise the branch that ROUTES. A counter written against
them alone would ship reading 0 everywhere, with no way to tell "correct,
nothing qualifies" from "broken, never fires".

## Predicted, at the default `--max_MAC_for_ER=4`

Cap is `kExactByWeightCarriers = 13`. A cell is ABOVE CAP when
`neff <= 4 && material > 13`.

| gene | anno | max MAF | site | n_eff | n_material | outcome |
|---|---|---|---|---|---|---|
| GENE1 | lof | 0.001 | ultra_rare | 1.08 | 5 | ROUTES (`.exactByWeight.txt`) |
| GENE1 | lof | 0.01 | ultra_rare | 2.06 | 45 | **ABOVE CAP** |
| GENE1 | lof | 0.01 | burden | 3.84 | 94 | **ABOVE CAP** |
| GENE1 | missense;lof | 0.001 | ultra_rare | 1.19 | 10 | ROUTES |
| GENE1 | missense;lof | 0.01 | ultra_rare | 2.33 | 54 | **ABOVE CAP** |
| GENE1 | missense;lof | 0.01 | burden | 4.69 | 114 | neither (`neff` > 4) |
| GENE2 | lof | 0.001 | both | 5.00 | 5 | neither (`neff` > 4) |
| GENE2 | lof | 0.01 | both | 9.00 / 22.0 | 9 / 22 | neither |
| GENE2 | missense;lof | 0.001 | both | 2.27 | 15 | **ABOVE CAP**, just over |
| GENE2 | missense;lof | 0.01 | ultra_rare | 2.85 | 41 | **ABOVE CAP** |
| GENE2 | missense;lof | 0.01 | burden | 3.77 | 74 | **ABOVE CAP** |

`GENE2 missense;lof 0.001` at `material = 15` is the sharp one: two above the
cap. `GENE1 missense;lof 0.01` on the burden path at `neff = 4.69` is the near
miss on the OTHER condition, so the fixture pins both edges.

NOT predicted and simply to be observed: the `0.0001` max-MAF tier and the
`missense;lof;synonymous` group, which this calculation did not cover. The
golden records whatever they do.

## What must NOT happen

**No existing golden moves.** The counter records and prints; it sets no
`isER`, no `ctx`, no p-value. If any committed golden moves, the change has
altered routing rather than observed it, and that is a bug in the change and not
a discovery. The `run_tests.sh standard` tier at `SAIGE_COMPARE_MODE=both` is
the gate for that claim.

Secondary: with `SAIGE_EXACT_BY_WEIGHT=0` the whole block is skipped, so the
counter must read **0** on the same fixture -- the kill switch covers the
observation as well as the routing, which is what makes the switch still mean
"the pre-#172 build".
