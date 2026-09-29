# Handoff 2026-09-29 (evening): the pinned-sample guard landed, its threshold did not survive the fixture

Read after `handoff-2026-09-26-flexrv-signature-repro-negative.md`, whose last
sections hold LEDGER #172 and the sas rerun. This one covers #173, #174, #175
and the null calibration of #172's routing.

No count from real data below 20 appears in this file. Every count below is
either from the bundled fixture (simulated) or a real count of 20 or more.

## What landed (four local commits, NOT pushed: CI minutes)

| commit | what |
|---|---|
| `0714dcce` | a real-data sample count below 20, quoted in four places since 2026-09-23, reworded |
| `6d0562ef` | seven `.exactByWeight.txt` sidecars that #172 produced and did not commit |
| `08c7e7a5` | LEDGER #174: the guard, the fixture, goldens P72-P76, LEDGER #175 |
| `2c6cd345` | LEDGER #173: DECIDED, with the negative-remainder safety net |

Gates: FAIL 0 / PASS 6906 (`08c7e7a5`'s tree) and FAIL 0 / PASS 6929
(`2c6cd345`'s), both at `SAIGE_COMPARE_MODE=both`. No existing golden moved in
either. Two full `make_goldens.sh` regenerations were byte-identical to each
other in all 175 step-0/step-2 files and in the four new step-1 files.

Still held locally from before: `23841f01` .. `955cee39`. The history on the
remote predates all of them.

## LEDGER #174: what the guard does

A sample is PINNED when its flip probability -- mu for a control, 1 - mu for a
case -- is at or below `SAIGE_SPA_PINNED_EPS` (1e-8). It leaves the vector the
saddlepoint's CGF is built from, in three places:

1. the marginal saddlepoint of `getMarkerPval` (pinned carriers on the fast
   kernel, every pinned sample on the dense one);
2. the conditional saddlepoint of `getMarkerPval`;
3. the saddlepoint R runs itself, `get_newPhi_scaleFactor_traitType`, which
   rescales the SKAT-O path's variance and the conditioning markers'. **The
   ledger row did not list this site; the fixture found it** when the SKAT-O
   burden p failed to move with the guard on.

A test leaves them out only while their summed flip probability is within
`SAIGE_SPA_PINNED_BUDGET` (1e-6). `SAIGE_SPA_PINNED_GUARD=0` is the kill
switch and is byte-identical to the build before the guard.

## The controls, all three together

Reference for the tested cells: **0.0196**, computed two ways that agree (an
enumeration of the ten material carriers against a normal remainder; a seeded
Monte Carlo of 1e5 draws).

| control | model | guard off | guard on, shipped (1e-8) | guard at 1e-6 |
|---|---|---|---|---|
| positive | planted, group at 1e-12, BURDEN | 7.0e-8 and 2.4e-13 | 0.0209 and 0.0209 | same |
| positive | planted, SKAT-O, all three columns | 2.2e-13 | 0.0207 | same |
| negative | planted, unweighted marker scan | golden | no p moves past the last digit | same |
| negative | a model with no pinned sample | golden | golden, `Pinned samples: 0` | same |
| **the finding** | **fitted, group at 1.7e-9..1.2e-7, BURDEN** | **0.030 and 0.055** | **0.030 and 0.054** | **0.0210 and 0.0210** |

The positive control behaves: the hazard is real in a cell that is not
one-atom, on both arms, and the guard removes it. The last row is the problem.

## Why the shipped threshold is partial

R's glm stops on the relative change of the deviance. A separated group's
share of the deviance is about 2 n mu, so the fit stops with the group near
`1e-8 x deviance / (1.26 x group size)`, whatever the data say. Measured on the
bundled 1,000 samples:

| separated group | fitted values | at or below 1e-8 |
|---|---|---|
| 8 samples | 1.7e-8 .. 4.8e-7 | 0 |
| 60 samples | 1.7e-9 .. 1.2e-7 | 13 |
| 150 samples | 1.5e-9 .. about 4e-8 for the carriers checked | 100 |
| 8 samples, no covariate offset | 1.4e-6 .. 3.9e-5 | 0 |

So a threshold on the fitted value splits a separated group wherever the fit
stopped, and the samples it leaves in the vector do the same damage as the
ones it takes out. At 1e-7 the fixture leaves 2 of 60 in and reads 0.026.

## What is NOT known, and how to find out

**The count of pinned samples measured on All of Us this morning (7 of 167
models) is a count of models with samples at or below 1e-8.** It says nothing
about how many samples of those seven models, or of the other 160, sit between
1e-8 and 1e-5. `[ASSUMPTION: the real models were fitted with the covariate
offset on, as the bundled fits are; if not, the groups sit higher still.]`

`harness/tools/aou_mu_decades.sh` reads the models and writes, per binary
model, the number of samples in each decade of flip probability and the mass
below 1e-8, 1e-7, 1e-6 and 1e-5. Every count from 1 to 19 is written `<20` and
its mass `masked`. Self-test: 16 checks, including a positive control that the
masking check trips on an unmasked table. Notebook cell: `aou_spa_rerun.md`
section 11.

## Decisions for the owner

1. **The threshold.** Candidates, with what the fixture says about each:
   - 1e-8, as agreed: partial. Complete only when the whole separated group is
     below it.
   - 1e-6: complete on the fixture. Each dropped sample leaves out at most
     1e-6 of probability.
   - not a threshold: call a sample separated when its fitted value is still
     falling geometrically under one more scoring iteration, wherever it sits.
     Needs the covariate glm re-run at step-2 set-up. Not built.
2. **The budget.** At 1e-6 the fixture's cells leave out 1.0e-6 to 1.6e-6 and
   the agreed budget of 1e-6 REFUSES them, which restores 0.030 and 0.055.
   Refusing is never better than leaving out once a pinned sample is active,
   so the budget is better a disclosure (the sidecar's `dropped_mass` column)
   than a gate.
3. **LEDGER #175**, the same blind spot in the fit gate's separation check. A
   separated group of 8 passes silently; one of 60 is reported as 13.

Recommendation: run the decade scan first. If the real separated groups sit
wholly below 1e-8, the shipped threshold is enough for these models and the
fixture is a warning about others. If they straddle it, take 1e-6 and turn the
budget into a disclosure.

## LEDGER #173

DECIDED: upstream's arithmetic stays. The safety net sends a test whose
fast-kernel remainder is negative or NaN to the dense kernel and counts it. It
has never fired; the plant `SAIGE_SPA_PLANT_NEG_REMAINDER=1` exercises it.

## LEDGER #172's routing: the null calibration at 100 replicates

`harness/work/null4k-20260929/rep100/` (gitignored), N = 4,000, replicates 11
to 110, routing on and off on the same phenotypes. All 200 step-2 runs exited 0.

| level | alpha | routing ON | routing OFF |
|---|---|---|---|
| gene (pooled row), n = 20,800 | 0.05 | 0.0420 [0.0393, 0.0448] | 0.0420 [0.0393, 0.0448] |
| | 0.01 | 0.0088 [0.0075, 0.0101] | 0.0087 [0.0075, 0.0101] |
| | 0.001 | 0.0011 [0.0007, 0.0016] | 0.0011 [0.0007, 0.0016] |
| | 0.0001 | 0.0001 [0.0000, 0.0003] | 0.0001 [0.0000, 0.0003] |
| per gene, Bonferroni minimum over routed cells, n = 11,400 | 0.05 | 0.0144 | 0.0159 |
| | 0.01 | 0.0025 | 0.0030 |
| | 0.001 | 0.0003 | 0.0004 |

Routed cells: 5.6 % of 1,947,700; of those, 5.0 % changed; no non-routed cell
changed. Among the changed, log10(p on / p off) has median 0.02 and a 5 % to
95 % range of -0.44 to 0.27. No routed cell below 1e-6 in either arm. The
routing does not move the gene-level error rate, which is at or slightly below
nominal at every level.

## All of Us state

- **sas**: rerun on the #172 image, stage check, control and gene-level
  comparison all done (earlier handoff sections).
- **eas**: the rerun is in (owner, this evening). Next are the stage check, the
  control, the summary and the comparison, each with `ANC = "eas"`. FemInf and
  MatHem ran on the #172 image, which has no guard: their rows are provisional.
- **amr, afr, eur**: not started. Hold back the models with samples at or below
  1e-8 (afr FemInf and MatHem, amr AMD and FemInf, eur MatHem).
- **No image has been built from these commits.** The default threshold is
  under question, so a candidate image would have to be rebuilt after the
  decision. The settings are environment variables, so one image can serve any
  threshold the owner picks.

## Not done

- The splice control does not yet accept a changed cell because it is listed in
  `.spaPinned.txt`. It must before the seven models are verified, and the
  wrapper and `submit()` must deliver that sidecar, as they had to for
  `.exactByWeight.txt`.
- The guard's effect on a whole null cohort is not measured. The fixture is
  twelve cells. A replicate null run on a phenotype with a separated covariate
  is the control that would show the guard does not move the error rate.
- `aou_spa_dumpcheck.sh` replays a guarded dump correctly (a dropped sample's
  entry is 0 and contributes nothing) but does not report how many were dropped.
