# Handoff 2026-10-01: the old-vs-new comparison on all four -- the ELEVEN ARE NOT THE RISK SET

Read with `handoff-2026-09-30-rollout-complete.md`, one of whose conclusions this
revises. Method and cells: section 7a of `harness/tools/aou_spa_rerun.md`.
ALL FOUR are in (sas, amr, afr, eur), on COMPLETE coverage -- see the method
warning below, which cost one verdict before it was caught.

No count of people appears here. Every number is a p-value, a gene count or a
(trait, gene) pair count.

## The aggregate, which answers nothing

| | binary genes | per-trait min | >2x up | spurious |
|---|---|---|---|---|
| eas (2026-09-30) | 275,519 | -- | 46 | 7 genes, found by inspection |
| sas | 109,589 | 18,087 | 60 | **0** |
| amr | 562,090 | 13,795 | 132 | **2** |
| afr | 562,119 | 13,797 | 80 | **0** |
| eur | 590,081 | 18,435 | 63 | **0** |

### THE METHOD WARNING, which cost a verdict before it was caught

The obvious old-side filter `source_old == "fixed"` is a BIASED SUBSET.
`flexrv_v4_fixed/` holds only the chromosomes the EARLIER filtered rerun
(#164/#166) touched, and those were selected by RECORDED SPA FAILURES -- a
criterion ORTHOGONAL to #172's damage, which is converged-but-wrong and was never
recorded as a failure. Measured: it cost eur **38 %** of its genes (366,633 of
590,081, per-trait median 11,339 against a max of 18,441) and amr almost nothing,
so it degrades one ancestry's verdict and not another's. Accept BOTH pre-#172
states -- `source_old.isin(["fixed", "v4"])` -- and CHECK the per-trait coverage
before reading any ratio; the deficit is invisible in the ratio itself. The first
pass of this comparison was run on the biased subset and read "eur is clean" from
62 % of eur's genes. Re-run complete, eur is still clean -- the concern was
legitimate and the answer did not change, which is the only reason this is a
method note rather than a correction.

The identical fractions from the biased pass, kept because they make the point
below: eas 93.8 %, sas 92.6 %, amr 97.0 %, afr 96.7 %, median ratio exactly
1.0000 in every one.

Median ratio is exactly 1.0000 in all four. **These numbers cannot answer the
question and eas proves it**: eas's row looks like the others and eas is where
seven genes moved from 3.3e-91..8.4e-9 to nominal. The signal is only in the top
movers broken down BY TRAIT. Do not read a verdict from an identical-fraction.

sas has the MOST movement per gene of the four, which was predicted to be the
least (it is the no-tail control). A plausible account -- at small N a
concentrated cell stays under the 13-carrier cap and so actually gets ROUTED,
where at large N it exceeds the cap and does not -- is POST HOC and is recorded
as a hypothesis, not a result.

## The result: sas and afr clean, amr carries the eas phenomenon

Criterion: a gene whose `p_FlexRV` was below 1e-8 before #172's fix and above
1e-3 after it was a SPURIOUS HIT being reported.

- **sas: 0.** Largest old p among its ten biggest movers is 1.1e-4 -> 3.7e-4.
- **afr: 0.** Largest is 1.0e-6 -> 1.4e-3. Nothing of the eas kind.
- **eur: 0**, on COMPLETE coverage (590,081 genes, per-trait 18,435-18,441).
  Largest among its ten is 7.5e-3 -> 2.3e-2; the rest are 0.03-0.22, i.e. null
  genes moving slightly. eur is the LARGEST cohort and the CLEANEST.
- **amr: 2 by that strict test, 3 by inspection:**

| trait | gene | v4 | v7 | best set | orders |
|---|---|---|---|---|---|
| BenCervUterNeo | ENSG00000101198 | 2.84e-51 | 5.99e-3 | `AM_cbrt0.4_MACsparse` | ~48 |
| BenIntNeo | ENSG00000196109 | 5.78e-45 | 3.79e-6 | `AM_sbeta20_mafbeta5` | ~39 |
| PAD | ENSG00000104415 | 1.39e-13 | 5.30e-3 | `AM_cbrt0.7_MACsparse` | ~10 |

That is also what amr's 3.2:1 upward skew in the aggregate was: real
spurious-hit removal, not noise.

## WHAT THIS OVERTURNS: flagged and at-risk are independent, in BOTH directions

The rollout handoff's census found ELEVEN models with mass at or below 1e-6 (AMD
in afr/amr/eur; FemInf and MatHem in afr/amr/eas/eur) and the eas comparison
found all ten of its top movers inside two of them. The natural inference --
that the eleven are the set at risk -- is now FALSIFIED both ways:

- **amr's three spurious hits are in BenCervUterNeo, BenIntNeo and PAD, none of
  which is among the eleven.** amr's flagged models contribute almost nothing:
  AMD's entry is 9.7e-4 -> 4.5e-2 and MatHem's best is 5.4e-8 -> 3.7e-5.
- **afr has FemInf and MatHem flagged and shows ZERO damage.**

**The reconciliation.** The census measured a MODEL-level property -- how many
participants sit below a fitted-probability threshold -- while the damage needs a
CELL-level one: ONE dominant carrier whose fitted probability is small enough to
send the saddlepoint to a phantom root. A single such carrier needs no tail at
all, and a large tail guarantees nothing unless a transform concentrates the
weight onto one of its members. So a model-level census cannot predict this and
must not be used as a risk filter.

The handoff's other conclusions STAND: the age tails cost nothing conditional on
#172 (QUESTION ONE), and nine of the eleven are the model being right. What does
not stand is any reading in which the eleven bound the blast radius.

**Practical consequence**: anyone who read PRE-FIX v4 FlexRV results for ANY amr
binary trait may have seen a spurious hit, not only for the flagged models. The
eas finding ("anyone who read v4 for eas FemInf or MatHem saw seven genes that
are not real") generalises past the flagged traits.

## LEDGER #176's own gene, confirmed

`amr EFRMB ENSG00000213030` moved **4.35e-18 -> 2.81e-16**: it changed, and it is
still extreme. That is exactly what #176 predicts -- the cell sits above the
13-carrier enumeration cap, so #172's routing cannot reach it -- and it confirms
that the 2.8e-16 `docs/state/next.md` calls "probably should not be there" is
live in current output. It is now at least VISIBLE, in
`.exactByWeightAboveCap.txt` (LEDGER #176's counter, landed today).

## The prediction eur refuted, for the second time

Written before eur landed: "on this evidence eur should carry spurious hits too,
and NOT confined to AMD/FemInf/MatHem". **eur is clean on complete coverage.**
That is the second eur refutation on this topic -- #176's N-scaling prediction
died on eur too (2026-09-30).

So **no covariate separates the dirty from the clean**:
- NOT cohort size. eur is the largest and cleanest; amr (N 52,192) is the
  dirtiest; eas (N ~ 4,050) is also dirty. Both ends of the range are clean and
  both are dirty.
- NOT the tail. eur AMD carries the largest tail of the eleven -- 2,827
  participants at 1e-6 -- and eur is clean.
- NOT trait identity. amr's three are in unflagged traits; eas's seven were in
  flagged ones.

**amr-specific remains the leading explanation and amr has prior form** (the
truncated step-0 GRM, the 6.91 variance-ratio artefact, amr/BMI at 13.8x obs/exp
in the pre-v4 run) -- but eas is dirty too and is not amr, so "amr-specific" does
not cover the evidence either. Nothing identifies the mechanism at ancestry
level, and the next discriminator should be sought at CELL level, not cohort
level.

## THE BLIND SPOT, which is why this whole exercise understates the problem

This comparison can only see damage that the fix REMOVED. A cell above the
13-carrier enumeration cap stays wrong AND stays extreme, so it does not move and
is invisible here BY CONSTRUCTION. `amr EFRMB ENSG00000213030` is the proof:
4.35e-18 -> 2.81e-16, still extreme, and it is known only because a human looked.

**And the cap binds MORE at large N**, because more carriers land in a cell. So
eur's cleanliness is equally consistent with eur having MORE unroutable
concentrated cells than anyone. The comparison cannot distinguish "no damage"
from "damage the fix could not reach", and it is biased toward finding less at
exactly the N where more is expected.

That is what LEDGER #176's counter (landed today) measures, and **nobody knows
the number for any ancestry**. It needs an image carrying the counter and one
run per ancestry; until then the five verdicts above are a lower bound on the
problem, not a measurement of it.

## Residual, recorded not chased

amr and afr each have at least one trait at ~75 % old-side coverage (per-trait
min 13,795 and 13,797 against medians of 18,432 and 18,433) even with both
pre-#172 states accepted. That is a gap in the ORIGINAL v4 production run, not
something the rollout introduced, and it is too small to flip a 0-or-2 verdict --
but it has never been explained and `chrStatus` for those traits would say.

## Next
- **The discordance screen** (`next.md` blocking item) should rank by the
  pre-vs-post ratio where both exist, and by FlexRV-vs-SKAT-O otherwise -- with
  `p_skato` above ~0.9 excluded (LEDGER #86 saturation).
- A gene list for the three amr hits is worth keeping: they are the only known
  real-data instances of the #172 artefact outside eas.
