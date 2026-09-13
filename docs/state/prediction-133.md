# Prediction for the LEDGER #133 fix, written BEFORE the gate

Date: 2026-09-13, at `58785e12` (after #134 and #132). Read off the COMMITTED
bytes, not from a run of the patched tree.

## The row's WHAT was mis-located, and that is fixed first

The row said the defect is "a cell's `MAC` is FOLDED while its `MAC_case` and
`MAC_control` are not". That framing sent both of its candidate fixes to the
wrong place:

- **Candidate A**, "fold the case/control halves the same way `MAC` is", is
  INCOHERENT for a mixed cell. Measured on `plink_group_urmixed_burden`: a
  per-stratum cell-level fold gives **negative** values (-1064, -997) on four of
  six cells, because the cell's `MAC` is a sum of PER-MARKER folded counts while
  the halves are raw cell-level sums. There is no cell-level fold that
  reconciles them.
- **Candidate B**, "report the unfolded AC in `MAC` too", was described as moving
  `MAC` in every group golden. It does not: only three goldens have any row where
  `MAC != MAC_case + MAC_control`.

**The actual defect is one marker, not the cell.** Every REAL marker is already
consistent: `imputeGenoAndFlip` flips the VECTOR and rewrites `altFreq`/
`altCounts`, and `src/Main.cpp:2746-2747` then derives `MAF`/`MAC` from the
post-flip values, so vector, counts and halves all describe the same allele. The
COLLAPSED PSEUDO-MARKER does not: since LEDGER #130 its vector is never folded,
and its halves are `cs.paritySumRanked(genoURVec_noweights, ...)` over that
unfolded vector -- but `src/Main.cpp:3444-3446` still folds its FREQUENCY through
`std::min(altFreq, 1 - altFreq)`, and `MAC`/`MAC_nw` follow from that. So the
pseudo-marker's number is folded while its vector is not.

Record A and B as mis-located so nobody re-opens them.

## The fix (CORRECTED -- the first attempt was wider than the defect)

At `src/Main.cpp:3444-3494`, `MAC` and `MAC_nw` become the RAW allele counts
directly, `altCounts * (1 - missingRate)`, instead of being derived from a folded
frequency. **`MAF` and `MAF_nw` STAY FOLDED.**

The first attempt made MAF raw as well, on the premise that it is "the reported
MAF". **That premise is false and the check is one grep**: nothing writes
`MAFVec` to any output file, and its only consumers are the two `dbeta()` calls
in `R/SAIGE_SPATest_Region.R` (`:324`, `:430`). The pseudo-marker's only reported
frequency is `AF_Allele2` in the single-variant file, which is `altFreq` and was
already raw -- which is why no `.singleAssoc.txt` moves either way.

So making MAF raw changed the WEIGHTING, not the reporting. Measured before it
was reverted, and kept as evidence for #131: every p-value stayed byte-identical
and `BETA/SE` held its ratio at -0.210239, but `BETA_Burden` went from -0.860935
to **-1.01393e+13** -- the coefficient per unit of a burden weighted to nothing,
because `dbeta(0.778,1,25)` = 5.1e-15 against `dbeta(0.222,1,25)` = 0.0605. A
golden containing -1e13 is a modelling decision, not a reporting fix, and it is
parked on #131 rather than taken here.

Nothing else changes: the halves already come from the unfolded vector, and every
real marker already folds its three counts together.

## Movers

| golden | UR rows above AF 0.5 | expected |
|---|---|---|
| `plink_group_urflip_burden` | 5 | MOVES |
| `plink_group_urflip_skato` | 5 | MOVES |
| `plink_group_urmixed_burden` | 5 | MOVES |
| every other group golden | 0 | DOES NOT MOVE |

`min()` is the identity below 0.5, and **no other golden's pseudo-marker comes
close**: the next highest anywhere is 0.1548 (the multiweight family), and most
sit at 0.0270. Measured across every committed group `.singleAssoc.txt`.

## Direction, inside the movers

`MAC` must RISE to equal `MAC_case + MAC_control`. On `plink_group_urflip_burden`
today: `GENE1 missense;lof 0.5` prints `MAC` 444 with halves 151 + 1405 = 1556,
and 444 = 2000 - 1556. After the fix `MAC` must read 1556. The three columns must
reconcile on every row of every group golden.

## Two second-order effects that must be named, not discovered

1. **#134's band and ER predicate read `MAC_nw`**, which does change value here
   (444 -> 1556). **But no p-value moves, and that was predictable**: both the
   folded and the raw value sit ABOVE the only live boundary, 20.5, so
   `selectVarianceRatio` returns the same `varRatio(1)` either way. Verified --
   every `Pvalue_Burden` is byte-identical on all six cells of both burden
   goldens. The earlier draft of this item said the p-values WOULD move; that was
   wrong, and the band edge is why.
2. **~~The beta-set pseudo-marker weight moves in the SKAT-O drill.~~
   WITHDRAWN by the corrected fix**: MAF stays folded, so the beta weight does
   not move at all and `BETA_Burden`/`SE_Burden` are unchanged (-0.860935 /
   4.09503, verified against the regenerated goldens). The paragraph below
   describes what the FIRST, wider attempt did and is kept as the #131 evidence.
   Original text:**
   `plink_group_urflip_skato` is the K = 1 legacy schema, so its one set is the
   implicit beta set and `AnnoWeights <- dbeta(MAFVec, 1, 25)` over every output
   row INCLUDING the pseudo-marker. `dbeta(x, 1, 25) = 25*(1-x)^24`, so the
   pseudo-marker's weight goes from `dbeta(0.222)` ~ 0.0606 to `dbeta(0.778)` ~
   1e-15 -- effectively zero. That is the honest reading (a pseudo-marker
   carrying 78% of the alleles is not a rare variant) but it is a large change
   and must be predicted rather than explained afterwards.

`plink_group_urflip_burden` and `plink_group_urmixed_burden` are `--r.corr=1`,
so they have no SKAT-O weighting to move; only their counts and whatever the
band change carries.

## The characterization assertion, and why it is not vacuous

For every non-Cauchy row of every group golden carrying a `MAC_case` column:

    MAC == MAC_case + MAC_control

**Red at HEAD in a known place**: it fails on exactly 15 rows, 5 in each of the
three drills, and passes on every other row of every other golden. Green after.
`test-invariants-goldens.R:105` already asserts this identity but ONLY on
`plink_group_fastTest.txt`, where no cell exceeds 0.5 -- the narrow assertion is
why the defect was never caught. Widening it to every group golden is the fix to
the test.

## The allele convention

`CLAUDE.md` rule 5 currently says the collapsed pseudo-marker is never folded.
It gains: **"and its reported MAC is its RAW carrier-allele count"** -- the count
only. The same entry records that its MAF stays FOLDED on purpose, because
`MAFVec` is never reported and feeds only `dbeta`, so that a later reader does
not "restore consistency" by making it raw and silently changing the weighting.

## Gate and commit shape

Standard both-modes run. One commit for one row: the fix, the widened
characterization, and three goldens re-baselined after the two-run
byte-reproducibility check -- the same shape as #134.


## OUTCOME with the corrected fix

`compare_goldens.sh` over a full regeneration: **exactly the three predicted
movers and nothing else**. On the drill cell `GENE1 missense;lof`:

| quantity | before | after |
|---|---|---|
| `MAC` | 444 | **1556** (= 151 + 1405, its own halves) |
| `BETA_Burden` | -0.860935 | -0.860935 (unchanged) |
| `SE_Burden` | 4.09503 | 4.09503 (unchanged) |
| `Pvalue_Burden` | 8.334812E-01 | 8.334812E-01 (unchanged) |

and on `plink_group_urmixed_burden`, `MAC` 2452 -> 3064 with `BETA_Burden`
-412.908 and `SE_Burden` 265.352 both unchanged. **The BETA/SE columns holding
still is the check that the weighting did not leak** -- they were the two
failures the first attempt produced, and their disappearance is what
distinguishes a counts-only fix from a weighting change.
