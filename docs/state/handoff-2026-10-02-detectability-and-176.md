# Handoff 2026-10-02: the three detectability instruments, #175's diagnostic, and the one open number on #176

Read this with `python3 harness/tools/ledger.py 176` and `ledger.py 175`. Everything
below is committed; **13 commits are unpushed and CI is blocked on the Actions
billing**, which the owner expects resolved. Nothing here is half-applied in the
tree: `docs/state/next.md` carries ANOTHER session's uncommitted paragraph (the
vectorised BURDEN screen) and was deliberately left alone.

## What closed

**The splice control's twin-rule gap (section 14b of `aou_spa_rerun.md`).** The
fix that section PROPOSED is refuted and the refutation is written down there:
it would have accepted a changed row whenever `AC_Allele2` matched, on the
premise that `AC_Allele2` is the collapse's UNWEIGHTED count. It is the
WEIGHTED one, so the rule would have re-created LEDGER #134's defect inside the
control, in the anticonservative direction. What was built instead is
**CONTROL 1c, a statistic-identity tier** in `harness/tools/aou_spa_splice.sh`:
a row whose ONLY differing columns are `SE` and `p.value` is accepted and named
in `<pre>.spliceStatId.txt`, with the total on the manifest and a `columns that
differ old -> new` line printed on any refusal. The self-test is now 8 positive
controls, 9 planted refusals and 3 stage-checker cases; `ur` and `twin` were
flipped from refusals to positive controls, and `sidvar`/`sidbeta`/`sidac` were
added so a plant that moves `var`, `BETA` or `AC_Allele2` still refuses.
Four blocked prefixes spliced (70/70, 70/70, 84/84, 84/84 -- the 140 predicted
beforehand was 2x because `comm -3` emits both sides of a changed row).
**amr and afr are now 124/124.**

**Old-vs-new on sas, amr, afr and eur.** `handoff-2026-10-01-oldnew-amr-afr-sas.md`.
The headline is that **the eleven flagged models are not the risk set**, in both
directions: flagged models land clean and unflagged ones move. eur lands clean
on complete coverage, and the old-side filter it was first run against was a
biased subset.

**The FlexRV-vs-SKAT-O discordance screen.**
`handoff-2026-10-01-flexrv-skato-discordance-screen.md`. One residual outlier in
~1.65M binary tests; FlexRV validated against the canonical genes.

**LEDGER #140 is seven files, not two.** The restart defect was recorded as
`.pooledTests.txt` and `.skatFailures.txt`; it is every end-of-run sidecar.
Results resume byte-exactly; the sidecars are overwritten and then describe only
the resumed portion.

## The three detectability instruments

The premise (`flexrv-detectability-before-defect`): at thousands of traits the
instrumentation blocks before the bug does, so each one must be a COLUMN
somebody reads, not a log line somebody greps.

1. `harness/tools/aou_flexrv_calibration.py` -- **per-trait calibration as a
   column**, the third and last. `null_median = 1 - 0.5**(1/m)` and
   `p_cal = -expm1(m*log1p(-min_p))` per (ancestry, trait, arm). Its self-test
   has 4 arithmetic checks, 20 null traits x 2 arms, and planted FlexRV-only and
   planted-both controls. **The tolerances are RELATIVE**: the first version
   compared at 1e-12 ABSOLUTE against six-significant-figure output, and the
   tell was that the one already-relative check (`ratio`) passed while the
   others failed.
2. and 3. are the two from 2026-09-30, unchanged.

## LEDGER #175 repair (b): report the effective sample size

`effectiveSampleSize()` and `reportEffectiveN()` in `R/FitGates.R`, wired into
`R/readInGLMM.R` after `checkModelEligibility`. Binary only;
`neff = sum(mu(1-mu)) / (pbar(1-pbar))`, reported **unconditionally** -- outside
`fitGatesEnabled()`, because the number is a description of the model and not a
gate, and `env.sh` turns the gates off for every harness run.

## LEDGER #176: the counter is built, and the number is NOT yet known

`src/Main.cpp` counts, at BOTH exact-by-weight routing sites (the ultra-rare
one and the BURDEN emitter), a cell that is concentrated enough for the exact
test (`neff <= --max_MAC_for_ER`) but has more material carriers than the
enumeration cap (`km == 0`, i.e. above `kExactByWeightCarriers = 13`):
`g_exactByWeightAboveCap`, `recordExactByWeightAboveCap()`,
`getExactByWeightAboveCapRecordsInCPP()`. `R/ExactByWeight.R` prints
**`Concentrated cells above the exact test's enumeration cap: N`** and writes
`<out>.exactByWeightAboveCap.txt` **only when N > 0**. Goldens P77 (SKAT-O,
3 routed + 7 above cap) and P78 (BURDEN, 3 + 5) over the new fixture
`extdata/input/group_new_snpid_withWeights_abovecap.txt`, which is P70/P71's
atom fixture with `1e-7` raised to `0.01`.

**The symptom, reproduced exactly on amr** (image `1.5.2-dev-7b7b447a`):

| set | BETA | SE | \|Z\| | p |
|---|---|---|---|---|
| `AM_sbeta10_unw` | 83.780 | 9.782 | 8.56 | 1.086297E-17 |
| `AM_sbeta20_unw` | 83.788 | 30.544 | 2.74 | 6.085320E-03 |

Two sets, the same BETA to four digits, SE a factor of 3.1 apart. The two
symptoms are ONE: the SE is back-derived from the p-value, so a different
p-value IS a different SE.

**What is open is one line of one log**, and it is the line that decides whether
#176 is a defect at all:

```
gs://<masked-identifier>/dsub/logs/fxcap-amr-efrmb/
  <masked-identifier>/20261002/090432/
  fxcap-amr---<masked-identifier>--261002-090432-76-task-None.log
```

`gcloud storage cat <LOG> | grep -iE 'above the exact|Exact test by weight|Effective sample'`

- **`: 0`** -- neither set is concentrated-above-cap, so the two sets took
  DIFFERENT TESTS and different p-values are CORRECT. #176 is not a defect.
- **`: 1`** -- one is and one is not: different code paths, which also explains it.
- **`: 2`** -- both are, in the same regime, and the discrepancy is genuinely
  unexplained. Only then are the wiring and a dump worth doing.

**Do not read the ABSENCE of `<out>.exactByWeightAboveCap.txt` from the job's
outputs as evidence of anything** -- that was claimed in this session and
corrected in `d13fa7e1`. The file is written only when N > 0, and the wiring to
deliver it does not exist yet, so its absence is equally consistent with a
count of zero, which would itself be the answer.

## Open, in the order it is worth doing

1. **Read that log line.** Minutes, and it may close #176 outright.
2. **Wire the sidecar** -- `OUT_FILE_EXACTBYWEIGHTABOVECAP` in `submit()` plus
   the pair in the wrapper's deliver loop. Correct either way, but pointless
   before (1).
3. **The popmax MAF check** (owner's, ~2 min).
4. **Push the 13 commits** once the Actions billing clears.
5. **#176's mechanism**, only if the count is 2 or more.
6. **Instrument A, `SAIGE_EXACT_TEST_LOG`** -- specified on 2026-10-01, NOT
   built. B (the statistic-identity tier) was judged sufficient for now.

## Two traps this session paid for

- **`submit()` defaulted `geno` to `bed`** and every production job uses pgen,
  which is why the first #176 run died at 61 s (`54c8d841`). A dsub failure
  under 60 s is localization; this one was just past it.
- **Asserting ahead of the evidence** cost four corrections in an hour, all the
  same shape: an unverified claim about the sidecar, a `head -2` image check
  that always exits 0, a grep for R source inside a compressed rdb, and
  `env_extra` as a `submit()` keyword that does not exist. The repo's own rule
  covers it -- read the live state, do not infer it.
