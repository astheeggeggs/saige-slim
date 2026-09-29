# Handoff 2026-09-29 (night): stepping back. READ THIS FIRST

The session before this one built more than the problem needed. The owner
called it ("I feel like we're overcomplicating things") and was right. This
file is the state a new session starts from. Detail is in
`handoff-2026-09-29-pinned-guard-and-threshold.md` (the guard) and
`handoff-2026-09-26-flexrv-signature-repro-negative.md` (LEDGER #172 and the
sas rerun); read them only if a question below needs them.

No count from real data below 20 appears in this file.

**HASHES.** Commit hashes cited in this file may not resolve in the current
history; the image tags are what they name, and the images are unchanged.

## The problem, in four lines

- A handful of All of Us FlexRV gene-trait pairs had implausible p-values
  (1e-53 to 1e-93) while SKAT-O was nominal.
- Cause: the saddlepoint on a weighted cell whose weight ONE ATOM carries
  (LEDGER #172, two mechanisms).
- Fix: such cells take the exact test. Image
  `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-a08fb295` (the notebook's
  `IMG9`), also on Docker Hub. Pushed, CI green.
- Verified on sas (rerun, stage check, control, gene-level comparison). The
  null calibration at 100 replicates is unmoved by the fix.

**That is the real problem and it is fixed.** Everything below is secondary.

## "The seven"

Seven null models, each ONE trait in ONE ancestry, have samples whose fitted
probability is at or below 1e-8 (read from the step-2 logs, 2026-09-29):

| ancestry | traits |
|---|---|
| afr | FemInf, MatHem |
| amr | AMD, FemInf |
| eas | FemInf, MatHem |
| eur | MatHem |

The other 160 models have none. In eas, 13 of 15 binary models have no sample
below 1e-5 at all; FemInf and MatHem have a smooth tail (62 / 103 / 175 / 310
and 77 / 141 / 244 / 394 samples at or below 1e-8 / 1e-7 / 1e-6 / 1e-5, no
gap). Both have N near 4,050, like the other female-specific traits, so the
other sex was excluded. `[ASSUMPTION: the tail is a steep age effect. Not
measured. Three jobs named mudrv-eas-{feminf,mathem,breastcanc} were running
at handoff and answer it; aou_spa_rerun.md section 12 has the reading cell.]`

## The hazard I chased, and its standing

LEDGER #174: a sample with a tiny fitted probability leaves a term in the
saddlepoint's CGF that can capture the root. **It has never been observed on
All of Us outside the one-atom cells #172 already fixes.** It was raised from
one dump and demonstrated on a fixture built to show it. For the seven, the
underlying issue is the null model (it assigns people a probability of
1e-16), which a patched test does not repair.

## The plan the owner was offered (NOT yet agreed)

1. Stop holding models back: run amr, afr and eur on `IMG9`, all traits.
2. On the seven, apply the check that found the original problem: a FlexRV
   pooled p that is extreme while SKAT-O is nominal.
3. Nothing shows: report the seven with a caveat, or leave them out.
   Something shows: look at that gene.
4. Treat the seven as a phenotype question (restrict to the at-risk group, or
   drop the trait). Owner's call.

**Open question put to the owner, unanswered at handoff:** for the seven,
caveat, exclude, or refit on a restricted cohort?

## What to DROP (do not resume without being asked)

- Tuning the guard's threshold; the budget debate; a "damage by decade" sweep.
- A new image. None is needed for the plan above.
- More All of Us scan tools. `aou_mu_decades.sh` and `aou_mu_drivers.sh`
  exist, are self-tested, and are enough.

## What exists and is frozen

- **The guard** (`08c7e7a5`): threshold 1e-8, budget 1e-6, kill switch
  `SAIGE_SPA_PINNED_GUARD=0`. Inert for a model with no sample at or below
  1e-8. Complete on the planted fixture, partial on the fitted one (P72 pins
  that). In NO image. Keeping it frozen or backing it out is the owner's call;
  the recommendation was to keep it and not develop it.
- **#173** DECIDED with a safety net (`2c6cd345`). **#175** OPEN (the fit
  gate's separation check counts only what is below 1e-8).
- **The splice control** reads `.spaPinned.txt` (`e3bf8edc`). Only needed if
  a rerun is ever made on an image with the guard; the wrapper and `submit()`
  would then have to deliver that sidecar.

## All of Us, where each ancestry stands

| ancestry | rerun on IMG9 | stage check | control | summary and comparison |
|---|---|---|---|---|
| sas | done | done | done | done, quiet |
| eas | done (all 15 binary, FemInf and MatHem included) | to do | to do | to do |
| amr, afr, eur | not started | | | |

Each step is the sas cell with `ANC` changed (`aou_spa_rerun.md` section 9).
In the eas comparison read gene ENSG00000128607 on MatHem first: it is the
positive control and should come down by many orders of magnitude.

## Local state

- Branch `strip`. **Everything up to the step-back handoff is PUSHED** (through the pre-push gate). CI failed in seconds on the Actions minutes limit, which is
  not a test result. Gates on the two code commits: FAIL 0
  / PASS 6906 and FAIL 0 / PASS 6929 at `both`; no existing golden moved. The
  rewrite left every tree identical outside the four reworded files.
- The installed package is `2c6cd345`'s tree.
- A background null run for the guard was in progress at handoff
  (`harness/work/null4k-20260929/pinned100/`, gitignored, `run.log`; 18 of 100
  replicates done, ~27 s each). It is optional under the plan above. Its
  design: 100 of 4,000 samples separated by a covariate; arms guard off,
  shipped, 1e-6, and the model refitted without the 100. No analysis script
  was written; `analyse_rep100.R` in the parent directory is the template.

## Lessons this session paid for

- A hazard that has not been seen in the data gets a CHECK ON THE DATA first,
  not a guard. The original detection signature was available all along.
- A number borrowed from a neighbouring component (the gate's 1e-8) carries
  that component's untested assumption.
- When a fix needs a threshold and the data have no natural cut, that is a
  sign the fix is in the wrong place.
