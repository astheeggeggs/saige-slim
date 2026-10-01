# Handoff 2026-10-01: the old-vs-new comparison on sas, amr and afr -- the ELEVEN ARE NOT THE RISK SET

Read with `handoff-2026-09-30-rollout-complete.md`, one of whose conclusions this
revises. Method and cells: section 7a of `harness/tools/aou_spa_rerun.md`.
eur is still running (320 files against amr/afr's 122) and is NOT in this file.

No count of people appears here. Every number is a p-value, a gene count or a
(trait, gene) pair count.

## The aggregate, which answers nothing

| | binary genes | identical | >2x up | >2x down | up:down |
|---|---|---|---|---|---|
| eas (2026-09-30) | 275,519 | 93.8 % | 46 | 46 | 1.0 |
| sas | 109,589 | 92.6 % | 60 | 50 | 1.2 |
| amr | 534,244 | 97.0 % | 128 | 40 | **3.2** |
| afr | 511,461 | 96.7 % | 73 | 33 | **2.2** |

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

## Next

- **eur's comparison**, still running. Prediction: on this evidence eur should
  carry spurious hits too, and NOT confined to AMD/FemInf/MatHem. If eur comes
  back clean, cohort size is not the driver and amr-specific returns as the
  leading explanation for #176 as well.
- **The discordance screen** (`next.md` blocking item) should rank by the
  pre-vs-post ratio where both exist, and by FlexRV-vs-SKAT-O otherwise -- with
  `p_skato` above ~0.9 excluded (LEDGER #86 saturation).
- A gene list for the three amr hits is worth keeping: they are the only known
  real-data instances of the #172 artefact outside eas.
