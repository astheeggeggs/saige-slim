# Handoff 2026-09-30: the rollout is complete, the tail question is CLOSED, #176 is the one thing open

Read this first. It supersedes `handoff-2026-09-29-step-back.md` on every point
they overlap. The day's detail is in `handoff-2026-09-29-eas-mu-drivers.md`
(which grew to hold all of it despite its date) and in sections 13 and 14 of
`harness/tools/aou_spa_rerun.md`.

**No count of people, alleles or carriers below 20 appears in this file; such
values are written `<20`.**

## What is done

All five ancestries are rerun on IMG9 (`1.5.2-dev-a08fb295`, the LEDGER #172
fix), stage-checked, spliced and read. 568 rerun jobs went in today for afr, amr
and eur; eas and sas were already rerun.

| ancestry | reruns | splices | summary rebuilt today | calibration read |
|---|---|---|---|---|
| eas | 15/15 | 15/15 | `15 v7` | all 15 traits calibrated |
| sas | 6/6 | 6/6 | `6 v7` | all 6 calibrated (no-tail control) |
| amr | 124/124 | **122**/124 | `122 v7` | done |
| eur | 320/320 | 320/320 | `320 v7` | done |
| afr | 124/124 | **122**/124 | `122 v7` | done |

Stage checks: `STAGE OK` everywhere, zero SPA fallbacks, no lost chromosome, no
`qsaddle=INVALID`. The four prefixes short of a splice are blocked by ONE control
gap (below), not by anything wrong with the data.

## QUESTION ONE, CLOSED: the age tails cost nothing

Every flagged model, in every ancestry that fitted one, reads calibrated or
concordant with SKAT-O:

| model | reading |
|---|---|
| eas FemInf / MatHem | smallest p at the null median, mid-pack among 15 traits |
| amr AMD / FemInf / MatHem | ratio to the null median 2.91 / 1.32 / 0.378 |
| eur AMD | extreme but CONCORDANT with SKAT-O (0.83 orders) -- a real signal |
| afr FemInf / MatHem | P(min <= obs) 0.63 and 0.031; one trait near 3 % in 31 is a uniform draw |

That spans cohorts from ~4,000 to ~170,000 and tails up to a quarter of a cohort
(amr AMD: 12,469 participants at or below 1e-5, 23.9 %).

**So no sample-exclusion method is needed.** LEDGER #174's guard stays frozen;
its hazard is now empirically absent on this data rather than merely unproven, and
#174 could reasonably be CLOSED rather than left open -- an owner's call, but it
has no open question attached.

**The one qualification, and it matters.** The tails are harmless *conditional on
#172's fix*. Before it, eas FemInf and MatHem were the ONLY models in their arm
carrying spurious hits: seven genes moved from a `p_FlexRV` between 3.3e-91 and
8.4e-9 to nominal, including the positive control ENSG00000128607 at
3.3e-91 -> 1.000000. The damage needs the COMBINATION of a tiny fitted
probability and a weight one atom carries; the tails supply the first, the
transforms the second.

## QUESTION TWO, ANSWERED AND THEN NARROWED: what the eleven actually are

"The seven" was a census taken at the fit gate's own `eps = 1e-8`. Recounted at
1e-6 across all 115 binary models it is **ELEVEN**: AMD in afr, amr and eur;
FemInf and MatHem in afr, amr, eas and eur. Four were invisible at 1e-8,
including eur AMD (2,827 participants at 1e-6, 16,505 at 1e-5, and **no warning
at all** from the gate).

**Nine of the eleven are the model being RIGHT, not failing.** Their decade
profiles hold nothing below 1e-9 or 1e-12 and rise monotonically toward 1e-5 --
the density of a population under a smooth age-risk gradient, i.e. finite
coefficients extrapolating a steep effect. Only the two SMALLEST cohorts (eas,
N ~ 4,050) carry mass at machine epsilon, fewer than 20 people per decade, with
`min_fitted` exactly `2^-52` in BOTH -- the floating-point floor, so an artefact.

`aou_mu_drivers.sh` measured the mechanism: `age` and `age2` carry 98.2 % to
99.4 % of the tail's gap on every model checked, with BreastCanc and EFRMB as
negative controls (same covariates, same cohorts, zero tail). EFRMB is the
instructive one -- female-specific and reproductive yet unflagged -- so the
predictor is whether the fitted age curve COLLAPSES, not the phenotype's biology.
**Trust the census, not the reasoning.**

## STILL OPEN: LEDGER #176

On amr EFRMB `ENSG00000213030`, `AM_sbeta10_unw` and `AM_sbeta20_unw` report
`BETA_Burden` 83.780 and 83.788 -- identical to four significant figures on the
same markers -- against `SE_Burden` 9.78 and 30.54, a factor 3.12. So `|Z|` is
bimodal across a smooth family of transforms: 0.51, 8.56, 2.74. One regime is
wrong and the 3.12 is unexplained.

Why #172's routing could not fire: it caps at 13 material carriers because the
exact test enumerates over them, and that gene's cells sit above the cap
(`n_material` 21 on the routed twins' siblings).

**Two predictions from it have FAILED.** eur (larger cohorts, 176 sets per gene)
is MILDER not worse, and its discordance is spread across the whole transform
ladder rather than concentrated on steep sets. afr adds a second failure: its
largest genuine discordance is AMD at 4.9 orders on a MILD `cbrt0.5` set.
**amr-specific is now the leading explanation**, and amr is the ancestry with
prior form (the truncated step-0 GRM, the 6.91 variance-ratio artefact).

### Next steps on #176, cheapest first

1. **A counter, not a fix.** Nothing in the code flags this: the run reported zero
   fallbacks and the splice control accepted the file. Make the run print how many
   cells were concentrated enough to route but too large to enumerate. That turns
   an invisible case into a log line, changes no result, and would have put this
   finding on screen instead of requiring someone to go looking.
2. **Reproduce it** with `SAIGE_SPA_DUMP_BELOW` and `SAIGE_SPA_DUMP_DIR` on that
   one region, replayed through `harness/tools/aou_spa_dumpcheck.sh` -- the
   instrument that found #172.
3. **A fixture**, buildable with no AoU data: a concentrated weight vector over
   more than 13 material carriers, which goldens P69-P71 do not cover.

## OTHER OPEN WORK, all recorded

- **The splice control's twin rule** (section 14b). One awk change -- read the
  pseudo-marker's own folded `AC_Allele2` instead of inferring from its 0/1 twin.
  It unblocks FOUR prefixes across two ancestries (amr Pancreat and Stroke, afr
  Asth and VTE, all currently mixed-provenance) and removes a refusal that will
  recur wherever a gene's rare variants are largely unscored by AlphaMissense.
- **LEDGER #175's repair (b)**: report the share of total sum-of-mu(1 - mu) lost.
  It changes no p-value, so its validation is a byte-identity test plus a fixture
  check. Quantified: amr AMD carries the information of ~40,000 people against a
  nominal 52,192; eur AMD ~154,000 against 170,540; the eas pair 4 % to 6 %.
  Repair (a), the guard's threshold, is parked with the guard.
- **Four old-versus-new comparisons** (`submit_summary_prefix` for sas, amr, eur,
  afr). No new compute. On eas this was the most informative read of the day.
- **afr AMD's top gene** at 7.4e-8 (P(min <= obs) 0.0014, past a genome-wide
  threshold, `p_skato` 5.8e-3, a mild set, so NOT #176's shape). Check it against
  known AMD loci before anything else.
- **The owner's decision on the eleven.** The evidence now supports a caveat for
  all of them. A restricted-cohort refit is the better ANALYSIS for the pregnancy
  traits, with the `--minCaseCount` floor of 100 to watch: eas FemInf has 103
  cases and MatHem 102, so a restriction that removes even a few cases refuses
  the trait.

## Five operational traps, now in section 14 of `aou_spa_rerun.md`

Each cost part of this session and each is invisible until it is known.

- **14. The scan stamp is PER ANCESTRY** and a wrong one reads as CLEAN. afr's is
  ten seconds and one digit from the notebook's scalar `STAMP`. Derive it from the
  bucket and `assert` the region-file count is non-zero.
- **14a. `submit_check_v7`'s old-side fallback is LOAD-BEARING off eas.** A
  hardcoded `flexrv_v4_fixed/` passes the whole eas chain and killed 130 of 444
  jobs elsewhere, because a finely batched ancestry has many batches with nothing
  to rerun.
- **14b. The splice control's twin rule has a gap** -- the rollout's first
  refusal, replicated on a second ancestry.
- **14c. The `summary_v4/` tables were STALE**, and sas's recorded "quiet"
  gene-level verdict rested on a pre-fix one. Only the log's `K v7` count and the
  table's `source` column distinguish them; the filename does not.
- **14d. A long dsub campaign stalls on accumulated Batch job records.** 5,310
  jobs in one region; the symptom is a HANG, not an error, and it silently lost a
  resubmission. Purge terminal-state jobs BEFORE a large campaign.

## A metric trap of my own making, worth not repeating

Ranking traits by `log10(p_skato / p_flexrv)` is contaminated when `p_skato`
SATURATES. afr MatHem sorts to the top at 5.77 orders and Gout at 3.95 purely
because their `p_skato` is 0.994 and exactly 1.0 -- SKAT-O returning exactly 1
once its smallest per-rho p exceeds ~0.78, which is LEDGER #86. **Flag or exclude
rows with `p_skato` above ~0.9 before ranking**, or the most null cells in the arm
present as the most suspicious.

## The CCT fallback, checked and benign

LEDGER #86's Bonferroni substitution (`min(1, n * min p)` when any set returns
exactly 1) fires on **3 of 18,401 regions** for eas MatHem, and all three have a
pooled p of exactly 1 -- every set null, so nothing was suppressed. The modified
CCT prevents the catastrophic case by construction; the conservative substitute is
rare and lands only where there is nothing to lose. Worth the same one-line check
on the other ancestries' `.pooledTests.txt`.

## Local state

- Branch `strip`, **22 commits unpushed**, all documentation. Held because the
  GitHub Actions minutes lapsed; `gh run rerun 36615578049` and `36615577980`
  re-run against what is already pushed once billing is fixed.
- The Batch registry was purged from 5,310 jobs to ~1,800 (pre-today SUCCEEDED
  only; today's kept in every state, because `dstat --full` on them is the only
  history of this session's diagnostics).
- `RUNS` in the notebook is ~2,663 entries and its `out_prefix` values were
  MUTATED by today's reruns to the v7 paths. A backup was written to
  `/home/jupyter/RUNS_backup_20260930.json`. The v4 paths are recoverable by
  substituting `fxv7` -> `fxv4` and `flexrv_v7` -> `flexrv_v4`.

## Lessons this session paid for

- **Read the error before proposing a cause.** Three diagnoses in a row came from
  reconstructing paths instead of reading a log; the log named the cause in one
  line the moment it was opened.
- **A check that can return "clean" from a wrong path is not a check.** The stamp
  trap, the stale summaries and one of my own verification cells all failed this
  way, each producing a confident and empty result.
- **Anything that deletes must print the name beside the id.** A mis-paired
  `job-id`/`job-name` list cancelled a needed job that had done 75 minutes of work.
- **State a prediction, then let it fail.** #176's N-scaling story was refuted
  within the hour by the next ancestry, and its cross-trait pattern by the one
  after. Both were withdrawn in the row rather than left to be inherited.
