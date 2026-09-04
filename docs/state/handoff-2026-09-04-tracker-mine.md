# Handoff: the upstream-tracker mine, Phases 1 and 2 (2026-09-04, evening)

Read `docs/state/next.md` for the detail. This file is the short version: what
landed, the gate numbers, and what is not done.

## What landed

**Phase 1 -- 24 LEDGER rows, #63 through #86, entries only, no fixes.** Mined
from the upstream issue tracker (plan:
`~/.claude/plans/can-you-take-a-dreamy-breeze.md`; **do not re-run the sweep**,
it cost ~160 agents). Every line number was re-resolved against HEAD and every
arithmetic claim re-derived before the row was written. Each row's PROVENANCE
says explicitly whether the consequence was REPRODUCED here or only traced --
read that field before trusting a row.

**Phase 2 -- five of the six planned items, one commit each.** Item 3 was dropped
because the finding behind it was refuted (below).

| item | test / golden | reaches |
|---|---|---|
| 1 | goldens P36/P37 + `test-step2-group-quant.R` | #72, and found **#86** |
| 2 | `test-step2-cond-highld.R` | #66, #84 |
| 4 | `test-step1-varratio-cursor.R` | #75 |
| 5 | `test-step2-numeric-marker-ids.R` | #74 |
| 6 | `test-step1-duplicate-pheno-id.R` | #67 |

Every test was drilled BOTH ways: green at HEAD, and failing when the pinned
behaviour is patched out (patched into a scratch library, then reverted and the
package reinstalled from clean HEAD). Every new R file is lint-clean. P36/P37
were proved byte-reproducible over three independent runs before being accepted,
and their runs are in `harness/make_goldens.sh` so a full-tier regeneration stays
consistent.

## Gate numbers

`bash harness/run_tests.sh standard`, and the same with
`SAIGE_COMPARE_MODE=numeric`. 0 failures throughout.

| commit | item | byte | numeric |
|---|---|---|---|
| `2fc00616` | Phase 1 (#63-#85) | 2128 | not run |
| `407a4041` | item 1 (+#86) | 2145 | 2320 |
| `983f39ac` | item 2 | 2182 | 2357 |
| `42eec450` | item 4 | 2199 | 2374 |
| `798f5f3c` | item 5 | 2213 | 2388 |
| `HEAD` | item 6 | 2226 | 2401 |

No timing figure is quoted from any of these runs: the machine carried heavy
unrelated load for part of the evening, which makes the tier slow but not wrong.
`harness/bench/NOTES.md` is still the place for timings.

## Refuted, and deliberately NOT written up

Recorded so nobody re-opens them. Detail in `next.md`.
- **empty ultra-rare collapse** (#12/#30/#60/#67 upstream) -- the C++ handles it
  coherently, checked at four points.
- **sparse-GRM nnz doubling** (#188) -- Matrix 1.7.6 keeps half storage through
  `*`, `diag<-` and `[i, i]`; the proposed mechanism does not exist here.
- **the `null_noXadj` unreachable `stop`** -- already LEDGER #39.
- **SKAT-O compact support** (#209) was refuted as "lives in the SKAT
  dependency, unlocatable" and then **OVERTURNED the same evening** by golden
  P37. The integral does live in the dependency; the consequence is here, and is
  now LEDGER #86. Treat the other three refutations as provisional in the same
  way -- they are "not found", not "cannot exist".

## Open, and not a ledger row

At the CLI default (`--is_noadjCov=TRUE`) a marker IDENTICAL to a conditioning
marker reports `var_c = 0.259444` out of `var = 28.0131` -- 0.93% residual
variance where there is none, and `p.value_c = 0.96` rather than 1. The obvious
mechanism (the tested marker is not covariate-adjusted while
`assign_conditionMarkers_factors` hardcodes `ctx.noadjCov = false`, so the
subtraction mixes two projections) predicts `var * R^2(G ~ covariates)` = 0.117
against 0.259 reported -- **off by 2.2x, so it is NOT established.** Worth
settling before anyone writes it up.

## What is NOT done

- **Phase 3 (fixes) has not started.** One commit per ledger entry, goldens
  regenerated, per the branch rule.
- **Most of #63-#86 have no test.** Only #63, #66, #67, #74, #75, #84 and #86 do.
  The rest are entries whose consequence is traced, not reproduced.
- Two fixes are already known one-liners, recorded with their drills: #74 wants
  `colClasses = list(character = 2)` on the bim `fread`
  (`R/Geno.R:92`), and #66 wants the relative cancellation test that LEDGER #44
  added to the marginal score tests.
- A **native AF > 0.5 fixture** is still absent (pre-existing gap, see CLAUDE.md).
