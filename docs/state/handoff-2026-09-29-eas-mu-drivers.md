# Handoff 2026-09-29 (late): what the eas tails are, and what that costs

Read `handoff-2026-09-29-step-back.md` first; it still says what to drop. This
file adds one measurement and four checks. The recipe cells are section 13 of
`harness/tools/aou_spa_rerun.md`.

No count from real data below 20 appears in this file.

## The measurement

The three `mudrv-eas-*` jobs of section 12 came back.

| model | N | cases | tail at or below 1e-6 | what carries the gap |
|---|---|---|---|---|
| eas FemInf | 4,054 | 103 | 175 | `age2` 63.3 % + `age` 36.5 % = 99.8 % |
| eas MatHem | 4,089 | 102 | 244 | `age` 50.1 % + `age2` 49.9 % = 100.0 % |
| eas BreastCanc | 4,074 | 174 | 0 | -- (control) |

No PC exceeds 0.2 % of the gap, the random effect contributes -0.1 %, no
categorical column appears, and `sex` entered none of the three fits. The gap
is -13.12 and -13.00 log-odds.

**The control is what makes this conclusive.** BreastCanc has the same twelve
covariates, the same female-only cohort and N within 20 of the other two, and
no tail at all. So the tail is not produced by the covariate design, the PC set
or the sex restriction. It is the trait's relationship with age: FemInf and
MatHem can only occur in reproductive-age women, so the model puts everyone
else at ~0 and is right to; breast-cancer risk is monotone in age, so nobody is
put at ~0.

**It reconciles two ways.** 103/4,054 = 2.54 %, logit -3.65, minus 13.12 is
`exp(-16.8)` = 5.6e-8; MatHem gives 5.8e-8. Both sit between the 1e-8 and 1e-6
rungs of section 11's decade table, where the counts are 62/175 and 77/244.

## What it changes

1. **The `[ASSUMPTION: steep age effect. Not measured.]` of the step-back
   handoff is measured and confirmed.** The competing hypothesis -- a separated
   categorical level, the fixture's shape -- is refuted for these two models.
2. **It IS separation, on a continuous covariate**, which `R/FitGates.R:249`
   says the mu check exists to catch. Only the CATEGORICAL shape is absent, so
   the covariate-cell scan is not the relevant half of the gate.
3. **LEDGER #175 is exercised by real data, not only by the fixture.** The gate
   counts at `eps = 1e-8`: 62 of 4,054 = 1.53 % and 77 of 4,089 = 1.88 %
   against a 2 % limit, so both models passed with a WARNING by a whisker. The
   tails actually sit at ~5.6e-8; counted at 1e-6 they are 4.3 % and 6.0 % and
   both would be refused.
4. **The census that produced "the seven" is at that same 1e-8.** A model with
   a shallower age effect keeps its whole tail between 1e-8 and 1e-6 and never
   appears in it. Section 11's decade scan, which counts at 1e-5 too, has been
   run for eas only.
5. **A stricter gate is the wrong fix.** `FitGates.R` refuses the whole fit:
   no per-gene, no per-sample granularity, so it would cost each affected trait
   entirely, all 22 chromosomes, and a `stop()` takes a multi-phenotype batch
   with it. The defect is who is in the cohort, not whether the trait can be
   analysed. Restricting to the plausibly-at-risk age window keeps the analysis,
   fixes the estimand, and removes the numerical hazard as a side effect.

## The gate's threshold is the wrong statistic, not the wrong number

Asked whether `eps` should scale with N. It should not, because that is not
what it scales with. Substituting deviance ~ 2 N H(p) into #175's formula, with
f the separated fraction:

    mu_rest ~ 1.6e-8 x H(p) / f

**N cancels.** Checked two ways: eas FemInf (p = 2.54 %, H = 0.118,
f = 175/4,054) predicts 4.4e-8 against the 5.6e-8 read from its driver log's
gap; the 60-sample fixture (p = 9.8 %, H = 0.321, f = 0.06) predicts 8.6e-8
against an observed 1.7e-9..1.2e-7. Both within ~1.3x, both high by the amount
null deviance overstates residual.

So the gate's percentage is f times whatever share of the group happens to rest
below 1e-8 -- an unstable slice of f. Two statistically identical situations can
read 1.5 % or 6 %. The number compared against `--separationMaxFraction` is not
the quantity that threshold was reasoned about, which is a worse problem than a
mis-set cut.

An `eps` built from `H/f` would be circular (f is what is being measured) and
would calibrate a scientific gate against glm's stopping tolerance. The gate's
own warning names the stable quantity instead: a pinned sample contributes ~0
to the score variance instead of mu(1 - mu), so the share of total
sum-of-mu(1 - mu) lost is dimensionless and needs no fitted-value threshold.
Recorded as #175's third and recommended candidate repair. **NOT BUILT**: no
damage downstream is yet demonstrated, and the last threshold hunt on this
mechanism is the session that had to be stepped back from.

## Why the SKAT-O arm is quiet, and why that is not reassurance

The hazard needs a tiny-mu carrier AND enough weight on that one carrier for
the saddlepoint root to pass `-log(mu)` (16 to 37 here). beta(1,25) mask
weights never concentrate that far -- which is why the beta-weighted `region`
TSAN workload never enters the #172 routing and `region_atom` had to be built.
FlexRV's transforms concentrate by design, and its Cauchy pooling is
approximately a minimum, so one pathological set contaminates the gene's
pooled p. Both arms share the same null model, so SKAT-O's good behaviour is
DIAGNOSTIC (it is what makes the discordance signature usable) and not evidence
that the models are sound.

## The #172 routing's dose-response, measured (2026-09-30)

eas MatHem x ENSG00000128607 -- the cell that opened the row at 1e-93 -- read
off `flexrv_v7`'s `.exactByWeight.txt.gz` on IMG9. Counts of material carriers
are withheld here as `<20` under this file's own rule; `n_eff` is Kish's
effective carrier count on the vector the test sees, a derived scale-free
statistic, and it is the finding:

| score transform | n_eff | n_material |
|---|---|---|
| cbrt0.4 | 3.309 | `<20` |
| cbrt0.5 | 2.722 | `<20` |
| cbrt0.6 | 2.421 | `<20` |
| cbrt0.7 | 2.239 | `<20` |
| cbrt0.8 | 2.128 | `<20` |
| cbrt0.9 | 2.079 | `<20` |
| sbeta3 | 2.107 | `<20` |
| sbeta5 | 1.097 | `<20` |
| sbeta10 | 1.00043 | `<20` |
| sbeta20 | 1.0000000131 | `<20` |

`n_eff` falls monotonically as the score transform steepens and reaches 1 to
eight decimal places at `sbeta20`: ONE carrier holds the whole weight. That is
the row's "one atom", measured on real data rather than on a fixture, and it is
a better argument for the fix than any p-value.

The outcome moved with it. `sbeta10` and `sbeta20` now read `Pvalue_Burden`
5.000000E-01 and 8.529162E-01 against 1e-93 in v4 -- the control of section 9
step 5(b), passed. The exact 0.5 is the MID-p of a degenerate enumeration
(nothing more or less extreme was reachable), not a failure. `BETA_Burden`
stays absurd at 1.2e+07, with `SE_Burden` tracking it to |Z| ~ 0.19, which is
the effect-size question parked on LEDGER #131 and is unaffected by this fix:
do not quote those betas.

**Three internal consistency checks pass**, which is what makes it more than
one number moving:

- `n_eff` is IDENTICAL across all six MAF transforms of each score transform.
  It must be: `Number_rare` is 0, so the cell is the pseudo-marker alone, every
  collapsed marker sits at essentially one MAF, the MAF transform is therefore a
  near-constant multiplier, and Kish's count is scale-invariant.
- The `ultra_rare` and `burden` rows carry identical `n_eff`. Also required:
  with no rare non-collapsed marker the burden vector IS the pseudo-marker.
- `unw`, `lof`, `cbrt0.2` and `cbrt0.3` are ABSENT from the listing, as the
  routing's own boundary requires -- it fires at `n_eff` at or below
  `--max_MAC_for_ER` (4), and `cbrt0.4` at 3.309 is the steepest set that gets
  under it. The gentler sets took the saddlepoint as before.

This is also direct support for the 13d prediction ahead of the per-set read:
the sets the fix rescued are exactly the concentrated ones. And it quantifies
why the SKAT-O arm was quiet throughout -- beta(1,25) mask weights never come
near `n_eff` of 4.

## The damage check: eas FemInf and MatHem are CLEAN (2026-09-30)

`aou_flexrv_summary.sh` on eas, verified to have read the post-fix files
(`flexrv: 25 v4 files, 15 fixed, 15 v7`; 459,932 (trait, gene) pairs over 25
traits, 18,397 genes each against the splice's 18,401 regions, 458,709 with both
arms).

Under a clean null the smallest of ~18,350 p-values has a median of 3.78e-5.

| trait | smallest `p_flexrv` | ratio to null median | P(min >= obs) |
|---|---|---|---|
| FemInf | 3.90e-5 | 1.03 | 0.49 |
| MatHem | 1.39e-4 | 3.68 | 0.078 |

FemInf sits on the null median; MatHem is on the CONSERVATIVE side. Nothing
resembling the 1e-53..1e-93 that opened LEDGER #172.

Three further signs the pathology is gone rather than less visible:
`set_best` is SCATTERED across the transform ladder (`cbrt0.2`, `cbrt0.3`,
`cbrt0.9`, `sbeta1.5`, `sbeta2`, `sbeta20`), which is what a null gives and the
opposite of 13d's failure mode; the worst FlexRV-vs-SKAT-O discordance is ~1,900x
(FemInf ENSG00000134285, 4.0e-4 against 0.77), an ordinary consequence of a
differently weighted test rather than the bug's 100-order gap; and `n_sets` reads
90-96, not 112, which is #142's rung dropping at N ~ 4,050 behaving correctly
(16 score transforms x 6 surviving MAF transforms, 90 where a score transform
keeps no marker).

**What this settles**: for these two models the age tails cost nothing
measurable, so the cohort question is about the ESTIMAND, not correctness. A
caveat would be honest; a refit would be an improvement, not a repair.

**What it does NOT settle, and must not be read as settling.** These are the two
SMALLEST cohorts of the eleven (N ~ 4,050) and the models with by far the largest
tails are untested: eur AMD (2,827 participants at 1e-6) and amr AMD (4,146 at
1e-6, 23.9 % of the cohort at 1e-5) have not been rerun on IMG9 at all. A 24 %
tail in 52,000 is a different regime from a 6 % tail in 4,089. **So AMD is the
priority within the amr/afr/eur reruns** -- a change from "all traits, order
immaterial", and the recount is what identified where the risk concentrates.

## THE SEVEN IS ELEVEN (recount, 2026-09-30)

`aou_mu_decades.sh` run over all five ancestries -- 115 binary models. The
counts below are of PARTICIPANTS and are all at or above 20 except where written
`<20`.

**Seven models have a sample at or below 1e-8, and that set IS "the seven"**,
which confirms the pinned scan's census was exactly the gate's `eps`. At 1e-6 it
is **ELEVEN**. Four models are invisible at 1e-8:

| model | N | cases | at 1e-8 | at 1e-6 | at 1e-5 |
|---|---|---|---|---|---|
| eur AMD | 170,540 | 2,407 | **0** | 2,827 | 16,505 |
| amr MatHem | 34,466 | 2,207 | 0 | 44 | 222 |
| eur FemInf | 102,504 | 1,853 | 0 | 38 | 483 |
| afr AMD | 60,499 | 158 | 0 | `<20` | 1,573 |

And the gate's percentage understates the group by one to two orders of
magnitude where it did fire: **amr AMD** reads 74 at 1e-8 (0.14 %, a warning far
under the 2 % limit) against 4,146 at 1e-6 and **12,469 at 1e-5, 23.9 % of the
cohort**. **eur AMD** produced no warning at all.

The full eleven: AMD in afr, amr and eur; FemInf and MatHem in afr, amr, eas and
eur. sas has none of the three (its six binary traits are Asth, BenIntNeo,
EFRMB, HTN, PAD, T2Diab).

**The 13b predictions held, including the negative half.** Only three trait
types are flagged in 115 models and nothing unpredicted appeared. FemInf and
MatHem are flagged in EVERY ancestry that fitted them. AMD is the sign-reversed
case, its tail being the young.

**Two columns say the mechanism is the same everywhere.**
`empty_decades_in_widest_gap` is **0 on all eleven**, so every case is a
continuum and the fixture's separated-cluster shape does not occur in real data
at all; and `contradicting` is 0 throughout, so no boundary sample has an
outcome that argues against its fitted value.

**EFRMB is the row that matters for the method.** Female-specific, reproductive,
fitted in both sas and eas, and flagged NOWHERE. So neither "female-specific"
nor "reproductive" is the predictor: what matters is whether the fitted age
curve actually collapses to zero over part of the cohort, which is empirical.
A biology-first list would have flagged EFRMB and been wrong. **The census is
what to trust; the a priori rule is a weak guide that happened to have no false
positives here.**

**On the `mu_rest ~ 1.6e-8 x H(p)/f` derivation of #175**: the ordering is
supported on 10 of the 11 -- eur FemInf and amr MatHem at a predicted ~3e-6 have
nothing below 1e-8, eas FemInf and MatHem at ~4e-8 have 35 % and 32 % of their
group below it. **amr AMD is a genuine exception**: predicted lowest at 4.7e-9,
yet only 1.8 % of its group sits below 1e-8.
`[ASSUMPTION: at f = 7.9 % the separated group is no longer the small
perturbation to the deviance that the derivation assumes, and the calculation
used null rather than residual deviance. Not established.]` Supported ordering,
one clear miss, NOT a confirmation.

**What this changes.** The owner's decision is about ELEVEN models, not seven,
and the two with the largest tails (eur and amr AMD) are ones the gate never
named. NEXT: the driver scan on the four newly flagged models -- AMD above all,
to test that its tail is the YOUNG -- then the same caveat/exclude/refit
question over the eleven.

## eas is VERIFIED at the file level (2026-09-30)

The rule at the top of this file is about PARTICIPANT counts; the counts below
are of tests and rows, which are not individual-level data.

Chain run on IMG9's `flexrv_v7`: tool checksums, sidecar-delivery gate, positive
control, stage check, splice control. Every step passed.

- **Sidecar gate**: all 15 binary traits delivered `.exactByWeight.txt.gz`
  (2.4-3.4 MB). The wrapper edits of section 9 were in when eas ran, so nothing
  needed resubmitting -- unlike the first sas pass.
- **Stage check**: `STAGE OK`, 15 sidecars, fallbacks 0, NA 0, traits with a
  failed chromosome 0. So no saddlepoint failed under IMG9, no cell went NA and
  vanished, and no chromosome was lost (#170).
- **Splice control** (`REGIONS=ALL`, `FAILURES=NONE`, 18,401 regions per trait):
  **15 of 15 clean**. Every suffix `all explained`, every row count identical
  old/new/spliced, zero refusals, no `ROWS` line, no `WARNING` or `NOTE`.
  `FAILURES=NONE` is legitimate because all 15 old sides are the
  `flexrv_v4_fixed/` copies of the #163-#168 round.

**The measured blast radius of #172 on eas.** ~2.41 M cells routed to the exact
test over the 15 traits (4,821,029 carried-over rows, which double-count the
`ultra_rare` and `burden` kinds); **120,923 cells actually changed, 3.2 % to
7.8 % per trait**; 400,146 changed rows in the results files, the remainder
being `MAC <= 4` updates to regions the earlier FILTERED rerun never touched.
So the fix is overwhelmingly a no-op and bites only where the saddlepoint was
wrong.

**Two predictions made before the 14 jobs ran, both confirmed.**

- Real-marker rows split on cohort size: non-zero for all NINE both-sex traits
  (6,699-6,874) and exactly zero for all FIVE female-specific ones. That is the
  `2*N*maxMAF > 10` boundary at N ~ 5,000, over 14 independent jobs.
- The unobservable-burden-cell rule -- the only one accepting at REGION rather
  than cell level -- fires only where `Number_rare > 0` cells exist, i.e. on
  both-sex traits. It fired on three (AFib 1 cell, CRF 2, T2Diab 1) out of
  2.41 M cells, three rising and one falling, so nothing systematic.

Two more quantities track cohort size the same way, which is why the picture
looks coherent rather than coincidental: routed cells carried over run
405 k-419 k on the female-specific traits against 306 k-315 k on the both-sex
ones, and pooled-test rows 54-60 against 21-24. Smaller cohort, more cells under
the ER cap, more to disclose.

**What eas still needs**: the section 7 summary and QQ jobs over
`flexrv_v7_checked/`, then the gene-level read of FemInf and MatHem for a FlexRV
pooled p that is extreme while SKAT-O is nominal. That is the damage check for
two of the seven, and it needs no further compute beyond those jobs.

## The four checks, cheapest first (section 13 has the cells)

1. **13a, the recount**: `submit_mu_decades` for sas, afr, eur, amr -- one
   model-only job each. Flag every model whose count at or below 1e-6 is not
   `<20`. This replaces guessing which traits are susceptible.
2. **13b** is a pre-registered prediction of what 13a will flag, written before
   it ran: FemInf and MatHem in every ancestry, amr AMD with the sign reversed
   (its tail should be the young), any other pregnancy-confined trait, and NO
   tail for monotone sex-specific traits of BreastCanc's shape. A flagged trait
   on none of those lines falsifies the rule.
3. **13c** is already measured except one number: the pinned table's eas/mathem
   control reads 77 of 4,089 at 2.22e-16. Prediction for the unquoted row: eas
   FemInf reads 62 of 4,054 = 1.53 %.
4. **13d**: the per-set rows of an extreme gene, from output already in the
   bucket. Prediction: the extreme rows are the concentrated sets (`sbeta20`,
   steep `cbrt`, `MAC1`/`MAC3` where `Number_rare` is 0). If `unw_unw` and the
   gentle beta sets are extreme too, weight concentration is not the whole
   story and the discordance signature should not be leaned on elsewhere.
5. **13e**: the damage check needs no special rerun. It is section 9's
   gene-level comparison on the IMG9 rerun of amr/afr/eur, which is needed
   regardless and should cover ALL traits -- a susceptible subset would cost a
   second pass and make the comparison's denominator a chosen set.

## Still the owner's, now better informed

For the seven: caveat, exclude, or refit on a restricted cohort? The evidence
points at restriction, and 13e says whether it matters at all.

## Local state

- Branch `strip` pushed to `slim/strip` at `8932505f` ; local and remote equal. The pre-push
  gate ran green on that tree at `SAIGE_COMPARE_MODE=both`. **CI says nothing
  about it**: both jobs died in 3-4 s on the account's Actions billing limit
  (runs 36615578049 `tests`, 36615577980 `docs`). `gh run rerun` when billing
  is fixed; until then every push burns a failed pair.
- The #174 guard null sweep (`harness/work/null4k-20260929/pinned100/`) stopped
  on its own during replicate 69 of 100 -- no error logged, cause not
  established. All 68 completed replicates read `off=0 / e8=0 / e6=0 / drop=0`,
  every arm agreeing. Left dead by decision; the guard is frozen and in no
  image.
