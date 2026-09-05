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
| `d95168c2` | #75 | `R/SAIGE_fitGLMM_fast.R` two deleted increments + `>=` guards | `test-step1-varratio-cursor.R` inverted, mixed pool added |
| (this one) | #63 | `src/SAIGE_fitGLMM_fast.cpp` the `- 1` on the pair-screen bound | `test-characterization-ledger63.R` inverted |

## Gate numbers

`bash harness/run_tests.sh standard`, and the same with
`SAIGE_COMPARE_MODE=numeric`. 0 failures throughout.

| commit | byte | numeric |
|---|---|---|
| start (`8aaed1f0`) | 2226 | 2401 |
| `22be2a69` (#74) | 2237 | 2412 |
| `c9f7e9ae` (#66) | 2253 | 2428 |
| `d95168c2` (#75) | 2260 | 2435 |

| `d74057f6` (#63) | 2262 | 2437 |

**NO GOLDEN MOVED for any of the four fixes**, predicted in advance each time
rather than discovered afterwards -- see below.

**The FULL tier was run at the end of the session and is green: 2289 pass, 0
fail, 0 skip.** `SKIP 0` is the load-bearing part -- `make_goldens.sh` +
`compare_goldens.sh` is gated on `skip_if_not(compare_mode == "byte")`, so a
run in numeric mode would have skipped the regeneration silently. It did not:
everything was regenerated and compared, step-0/step-2 byte-exact and step-1
toleranced, plus the determinism double-run and the permuted-phenotype control.
**CI on `slim/strip` is green too** (ubuntu, numeric mode) at `d49aac8e`, which
is the check that matters most for #66 -- it moved two region writers from a
raw `<<` to `writeNumOrNA`, and NaN spelling is exactly the platform split that
cost two red CI runs before. No golden was re-baselined, so
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

## #75 -- the crash was the smaller half

`indexInMarkerList` advanced once unconditionally and a SECOND time on each of
the two skip branches. The fix is not a judgement call: the commented-out
`# indexInMarkerList = indexInMarkerList + 1` still sitting inside the ACCEPTED
branch records that the increment was hoisted to the top of the iteration and
the two copies were left behind. Both are deleted; the exhaustion guards become
`>=`, which with a one-per-iteration cursor is no longer load-bearing and is
kept as insurance.

Four run-time-derived pools, each measured at HEAD and after (n = 1000,
MAC 40, `--numRandomMarkerforVarianceRatio` at its default 30):

| pool | HEAD | fixed |
|---|---|---|
| 40 autosomal | 30 read, 30 tested | identical (negative control) |
| 40 non-autosomal | 20 read | 40 read, still fails |
| 41 non-autosomal | 22 read + NA crash | 41 read, clean exit |
| 30 usable + 30 not | **24 of 30 tested, exit 0** | 30 tested |

**The last row is the finding and it is bigger than the crash the entry was
written around**: a variance ratio returned at exit 0, with no warning, from a
silently thinned marker set. The all-skipped pools are what show the fix reads
the whole pool without inventing data -- with no usable marker the run still
fails, leaving the `.rda` on disk and the variance-ratio file empty.

**A trap that cost a measurement cycle, now recorded in the test file**: the
variance-ratio marker list is ordered by MAC, not by file order. A first draft
of the positive control assumed file order, predicted that HEAD would fail on
an alternating pool, and was wrong -- HEAD succeeded there. Only the TOTALS are
derivable; every assertion in the file is order-independent for that reason.

## #63 -- one character, and the goldens really do not move

The bound is `C(N,2)` now, not `C(N,2) - 1`. A live step-0 ingest prints
`totalCombination: 499500` against 499499 before.

The provenance predicted both step-0 goldens would stay byte-identical, because
the pair this fixture drops is 1-based (500,501) with a refineKin kinship of
-0.0087 against a 0.05 cutoff. That prediction was CHECKED rather than trusted:
`createSparseGRM.R` re-run through both readers, `cmp` against
`harness/goldens/step0/` -- `.mtx` and `.sampleIDs.txt`, PLINK and PGEN, all
four identical. This is the LEDGER #44 shape, and worth restating because it is
the branch's most repeated lesson: **a fix whose goldens do not move is a
statement about coverage, not a pass.**

Drilled both ways with a clean rebuild each: 22/22 with the fix, 8 failures at
unpatched HEAD naming the dropped pair as `5:6` at N = 10 and `5:11` at N = 11,
which is exactly the hand derivation.

**Deliberately NOT touched**: the debug `printComb` carries the same `- 1`, but
it is not a copy of the worker -- it folds on `j < i` where the worker folds on
`j <= i`, so it also emits a diagonal pair the worker never would. It only
prints and R calls it with N = 3; correcting it means deciding what that printer
is for. The commented-out OpenMP block carries the `- 1` too and is dead.

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

## One assertion was walked back before pushing

Worth reading if you are about to write a count-based assertion here. #75's
mixed-pool test first asserted `read == 60` and `skipped == 30`. Both are true
on this machine and both are reproducible -- but neither is DERIVED: they
depend on where the seeded `sample()` permutation happens to place the last
usable marker. An intermediate note claimed the CV check forced the loop to
exhaust the list, which is simply wrong -- the CV is 1.2e-04 against a 1e-03
cutoff, so it never triggers a marker-count increase.

They were replaced with the derived bounds, leaving `tested == 30` as the
assertion with teeth (HEAD tested 24). The file's own header claims every
assertion is order-independent, and the point is that the header and the
assertions have to agree: a test whose comment overstates what it checks is
worse than one that checks less.

## What is NOT done

- **`python3 harness/tools/ledger.py open` is the open-row list and count;**
  do not read `LEDGER.md` whole. This file first said "38 of 86", and the
  lint added on 2026-09-05 (`ledger.py lint`) then found five more OPEN rows
  that the tool had never listed because their prose carried bare pipes --
  #26, #45, #64, #65, #77 -- so the figure at this handoff was 44 of 89.
  Hand-written counts drift; quote the tool.
- Of the tracker mine's rows, **#63, #67, #75, #84 and #86 have tests** and the
  rest are traced, not reproduced. `docs/state/next.md` names the two coverage
  holes most of them sit behind, and the four REFUTED candidates that must not
  be re-opened.
- The **open question** from the previous handoff (`var_c` = 0.259 where the
  obvious mechanism predicted 0.117) was SETTLED later on 2026-09-05: see
  LEDGER #82 and `docs/state/next.md`. The mechanism was right, the arithmetic
  was missing the centring-mismatch term and the mu(1-mu) weighting.
- A **native AF > 0.5 fixture** is still absent (pre-existing gap, CLAUDE.md).
