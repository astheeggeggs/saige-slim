# Handoff: Phase 3 starts -- LEDGER #74 and #66 (2026-09-05)

Two fixes, one commit each, per the branch rule. Both were the "already known
one-liner" pair the previous handoff named, and both had characterization tests
that were INVERTED in the same commit as the fix.

Read `docs/state/handoff-2026-09-04-tracker-mine.md` for the session before this
one (Phases 1 and 2 of the upstream-tracker mine, the refuted candidates, and
the one open question), and `docs/state/handoff-2026-09-04.md` for where CI
lives.

## What landed

| commit | entry | code | test |
|---|---|---|---|
| `22be2a69` | #74 | `R/Geno.R` bim + pvar `fread` `colClasses` | `test-step2-numeric-marker-ids.R` inverted, pgen arm added |
| `c9f7e9ae` | #66 | `src/SAIGE_test.cpp` relative guard + `src/Main.cpp` two writers | `test-step2-cond-highld.R` #66 block inverted |

## Gate numbers

`bash harness/run_tests.sh standard`, and the same with
`SAIGE_COMPARE_MODE=numeric`. 0 failures throughout.

| commit | byte | numeric |
|---|---|---|
| start (`8aaed1f0`) | 2226 | 2401 |
| `22be2a69` (#74) | 2237 | 2412 |
| `c9f7e9ae` (#66) | 2253 | 2428 |

**NO GOLDEN MOVED for either fix**, which was predicted in advance both times
rather than discovered afterwards -- see below. No golden was re-baselined, so
the "prove byte-reproducible twice first" rule did not have to be exercised.

## #74 -- the fix, and why no golden could gate it

`colClasses = list(character = 2)` on the bim `fread`, `list(character = 3)` on
the pvar. The number indexes the FILE's columns, not the `select`ed ones; that
was checked against data.table for both call shapes before the edit.

Every committed `.bim` and `.pvar` spells its IDs `rs<N>`, which already inferred
character, so the pin is a no-op on all 36 step-2 goldens -- the entry was right
that no golden could ever reach it.

The fixtures are derived at run time, reusing the committed `.bed`/`.pgen` and
`.fam`/`.psam` UNCHANGED, because a genotype file is indexed by row. That is what
makes the positive control EXACT rather than approximate: the numeric-ID run and
the rs-spelled run must agree on every column but `MarkerID`, and do.

**A pgen arm was added.** Without it the pvar half of the fix was ungated, and
that half is not hypothetical -- a pvar's ID column is as free to be all-numeric
as a bim's.

Drilled both ways: 6 failures at unpatched HEAD with `Incompatible join types:
x.ID (integer) and i.SNP (character)` in the log. The plain-marker-scan NEGATIVE
CONTROL passes in BOTH builds -- that path never joins on ID, so the bug was
confined to the joining paths exactly as the entry said.

## #66 -- the fix, and the prediction that made it cheap

The conditional variance is a DIFFERENCE, `var1 - G1'PG2 VarInv G2'PG1`, so a
marker collinear with the conditioning set cancels it to float noise. Its only
guard was `varT_c <= DBL_MIN` -- the ABSOLUTE test LEDGER #44 had already shown
to be the wrong question for a cancelling sum. The fix uses **#44's test, #44's
threshold and #44's helper**: `varianceCollapsed(varT_c, |var1| +
|varT_ctemp(0)|, m_n)`. No new constant is introduced.

**Measured on the constructed high-LD fixture** (n = 1000, so n*eps = 2.22e-13):

| marker | collapse ratio | vs threshold |
|---|---|---|
| cond, ld100 | 4.89e-16 | 450x below -> fires |
| sigpart | 7.19e-16 | 310x below -> fires |
| ld090 (r2 = 0.9) | 5.02e-02 | ~11 orders above |
| smallest at `--is_noadjCov=TRUE` | 4.1e-03 | ~10 orders above |

Fourteen orders of separation. The CLI default is untouched, which is #82.

**The no-golden-moves prediction was computed from the goldens themselves,
before the change.** Because `varT_ctemp = var1 - varT_c`, the collapse ratio is
`var_c / (2*var - var_c)` -- readable straight out of the printed `var`/`var_c`
columns, no instrumented build needed. Its minimum across all nine conditional
goldens is **1.655e-03**, ten orders above the threshold. The gate's zero
failures is the confirmation. Worth reusing: this trick turns "will my guard
fire anywhere it should not?" into an arithmetic question over committed files.

Drilled both ways with a clean rebuild each time: **7 failures at unpatched
HEAD**, showing the published junk verbatim -- `BETA_c = -0.0178571` with
`SE_c = 6.34e+06`, `BETA_c = 0.333333` with `SE_c = 4.84e+06`, at
`p.value_c = 9.999999E-01` (so the old guard genuinely never fired). Note the
ratio assertions PASS at HEAD: the arithmetic is unchanged by the fix, only the
guard's action is, so those lines are descriptive and the NA assertions are the
discriminating ones. That is the right shape for an inverted test.

## Two things found on the way, both recorded

- **A gap LEDGER #44 left.** `writeOutfile_singleInGroup` and
  `writeOutfile_BURDEN` printed effect sizes with a raw `<<` where
  `writeOutfile_single` used `writeNumOrNA`. So an `NA_REAL` reaching the REGION
  path would have printed as the platform-dependent `nan`/`-nan` that #44 exists
  to avoid. All four effect-size fields in each writer are routed through the
  helper now. Byte-neutral: no committed golden carries a non-finite effect size,
  which is also why nothing had ever noticed.
- **LEDGER row #66 was unparseable.** Its WHAT cell contained the literal pipes
  of `|Beta_c|`, so the markdown table split wrong and
  `harness/tools/ledger.py` had been reporting scrambled fields -- the row's
  STATUS appeared inside PROVENANCE, and `ledger.py open` never listed it. Any
  row written from now on must avoid a bare `|` in its prose. Its PROVENANCE
  also still said "not reproduced" when Phase 2 item 2 had reproduced it; both
  corrected.

## What is NOT done

- **37 of 86 ledger rows remain open.** `python3 harness/tools/ledger.py open`
  is the list; do not read `LEDGER.md` whole.
- Of the tracker mine's rows, **#63, #67, #75, #84 and #86 have tests** and the
  rest are traced, not reproduced. `docs/state/next.md` names the two coverage
  holes most of them sit behind, and the four REFUTED candidates that must not
  be re-opened.
- The **open question** from the previous handoff is still open: at
  `--is_noadjCov=TRUE` a marker identical to a conditioning marker reports
  `var_c` = 0.259444 out of `var` = 28.0131, and the obvious mechanism predicts
  0.117 -- off by 2.2x, so it is NOT established. This session's fix does not
  touch it (that path's collapse ratio is 4.1e-03, ten orders above the guard),
  and the measurement above reproduces the same regime if anyone wants to settle
  it: `test-step2-cond-highld.R`'s `def` run is exactly it.
- A **native AF > 0.5 fixture** is still absent (pre-existing gap, CLAUDE.md).
