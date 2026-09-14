# Handoff 2026-09-14: whole-codebase quality audit

**Read this file end to end before acting on any section of it.** It is long and
its sections cite each other; a section read in isolation will mislead you about
what was already refuted.

## What this is

A structural code-quality audit of the whole `strip` branch, run 2026-09-14 at
commit `9deab7c9`. It is NOT a bug hunt. Per this branch's bug-for-bug rule,
where the audit uncovered a latent defect the deliverable is **a LEDGER row, not
a fix**.

## Provenance -- what produced it, and how much to trust it

A 19-unit fan-out: 13 subsystem readers (each assigned a bounded set of files)
and 6 cross-cutting lenses (god functions, the R<->C++ contract, duplicated
knowledge, guard/error architecture, dead and thin surfaces, concurrency
ownership). Each was capped at 4 findings to force ranking rather than padding,
and each was required to record what it searched the ledger for before
proposing anything.

**96 findings raised. Each was then faced by two independent skeptics:**

- one attacking the FACTS -- are the line numbers, counts and structures real,
  recount them, is there an existing mechanism the finder missed;
- one attacking PRESERVATION and VALUE -- is this a ledgered behaviour that must
  not change, could the remedy move a golden or break an asserted string, does
  it actually delete complexity or merely relocate it, is it upstream-inherited,
  is there a recorded decision against it.

Both were instructed to default to "refuted" when uncertain. **62 survived, 34
were refuted.** A completeness critic then named 5 gaps, which were covered by
targeted finders under the same verification. 218 agents total.

**Every verifier correction is applied in the text below. Its line numbers and
counts supersede the raw findings** -- several findings were materially wrong
until a verifier recounted them, and several remedies are narrower here than
their finder proposed.

## How to read the tags

- **(2/2 verified)** -- neither skeptic refuted it.
- **(1/2, contested: X)** -- one skeptic refuted part of it. **The item shown is
  the narrowed version that survived**, and X names the claim that did not. Do
  not reinstate X.
- **Section 8 is the do-not-re-raise list.** Five of the 34 refutations were
  refused precisely because this repo had already decided them. That list is
  load-bearing, not filler.

## What is NOT done

No code was changed. No ledger row was written. No gate was run. Nothing below
has been implemented.

**The raw per-agent results are gone.** The workflow journal and per-agent
transcripts were session-local and are not preserved -- only this synthesis is.
Re-running the audit is possible but would not reproduce agent-for-agent, so
treat this document as the artefact of record.

## One known defect in this document

The synthesis agent's return was truncated to its tail; the first two-thirds was
recovered from its transcript and the two halves stitched. The join was verified
(the first block ends `...since a single` and the second opens `-model run will
now execute...`, forming one sentence in section 5.5). No other gap was found,
but a section that reads as if it is missing a premise may be worth checking
against the surrounding text rather than assumed complete.

## Tree state when this was written

Clean. `9deab7c9` ("fix: LEDGER #142 -- FlexRV drops the MAC rungs at or above
2\*N\*maxMAF") landed while the audit was running. `ledger.py lint`: **142 rows,
0 open, 7 decided**. `docdrift.py`: 18 rules over 39 files, 0 hits.

**That "0 open" is about to change.** The audit surfaced nine latent-defect
candidates (Phase B of the execution queue). Recording them flips the branch's
headline state, which CLAUDE.md and `next.md` both quote. The recommendation is
to write them in their own commit, separate from this handoff, so the state
change is separately reversible. **That decision is the next session's first
one.**

## Outstanding from the earlier `/thermo-nuclear` review of the #142 diff

Separately from the audit, the #142 diff itself was reviewed before it was
committed. Two items landed; three did not, and they are verified still present
at `9deab7c9`:

- **`R/SAIGE_Test_main.R:114`** (roxygen for `flexRV_score`) and its generated
  twin **`man/SPAGMMATtest.Rd:210`** still assert "the 16 x 12 = 192 FlexRV
  weight sets" unconditionally. The CLI-help twin at
  `extdata/step2_SPAtests.R:567` and the `flexRV_MAFTransforms` roxygen WERE
  corrected. One twin fixed, the other not -- which is verbatim the failure mode
  CLAUDE.md describes ("COPIED from each other ... one twin had already been
  corrected while the other was not") and which `docdrift.py` has no rule for.
- **`docs/design/flexrv-design.md:607`** still says "byte-inert on the 192-set
  run". (`:524` was corrected; that one landed.)
- Worth a `docdrift.py` rule for the unconditional "16 x 12 = 192" claim. The
  tool's own bar -- add a rule only after a claim has actually been written
  wrongly -- is met twice over.

Landed, no action needed: `R/BurdenMemo.R:21` now reads 224, and
`test-step2-flexrv.R:50` is now `expect_lte(nums[3], 1L)` -- a measured bound,
and better than the `2L` the review proposed.

---

## The report

62 findings survived two rounds of adversarial verification. They do not describe 62 independent problems. They describe **four** structural facts, each of which has already cost this branch ledger rows, and a fifth about the instruments used to gate changes.

1. **One rule, written N times, fixed once.** The branch's most expensive recurring defect is a decision rule duplicated across copies where a fix lands in one copy and not the other. It is the documented cause of `#44`→`#66` (two single-variant writers), `#118`→`#134`→`#135` (three MAC-gated decision sites), `#91` (two `getDiagOfSigma`), `#37` (two phenotype blocks), and it is live and unrepaired today in the binary conditional region kernel (§4.1) and the region restart path (§1.2).
2. **The R↔C++ boundary has no contract.** A 303-line header that one translation unit includes and the compiler cannot check; config pushed into globals before R has finished deriving it; a reader interface whose capability questions nobody asks; a `t_genoType` string threaded past the polymorphic object it describes.
3. **~3,000 lines of dead surface**, most of it already inventoried by the 2026-09-08 audit and parked — plus a slice that audit's method *structurally could not see* (C++-internal functions, dead *parameters*, post-`stop()` tails, an uncompiled translation unit).
4. **Two god functions** (`mainRegionInCPP` 2,501 lines; `SAIGE.Region` 1,251 lines) that grew back past upstream size *after* this branch extracted 675 lines out of one of them.
5. **Three controls that do not control**: `-w` disables the only static analysis in the build; `tsan.sh` aborts after the first arm that dies under `plant`, which is the mode where dying is the expected outcome; `prepush.sh` records a run shape it did not execute.

Nothing below is "fix this bug". Where a structural finding uncovered a latent defect, the deliverable is **a LEDGER row**, per this branch's own rule.

---

## How to read this

- **(2/2 verified)** — neither skeptic refuted it. **(1/2, contested: …)** — one skeptic refuted it; the item below is the *narrowed* version that survived, and the contested claim is named.
- Every verifier correction has been applied. Counts and line numbers here supersede the original findings.
- **Gate arithmetic.** Any change under `src/` or `R/` = full gate, both modes, `bash harness/prepush.sh`, ~7 min warm. A `.hpp` change additionally forces `rm src/*.o` → ~167 s clean build. Docs-only = `ledger.py lint` + `docdrift.py`, seconds. Test-only = the changed file at `standard`, escalating to both modes only if it adds a cross-run comparison. The gate is per **push**, so batching commits is free; batching *risk* is not.
- Where a claim is "cannot move a golden", that is a prediction to state **before** the gate, not a pass. Per CLAUDE.md's first Trap, a green gate on an inert change is a statement about coverage.

---

## §0 — The action queue

Ordered by value per unit of risk.

| # | Item | Cost | Section |
|---|---|---|---|
| 1 | `tsan.sh` early-abort under `plant` | 2 lines, **no gate** | §2.1 |
| 2 | Ledger rows for 9 latent defects the audit surfaced | `ledger.py lint` | §1 |
| 3 | `prepush.sh` logs what it ran | ~12 lines, no gate | §2.3 |
| 4 | The dead-surface sweep (~3,000 lines, 6 commits, one push) | 1 gate + clean build | §3 |
| 5 | Drop `-w` from `src/Makevars`, ratchet the warning count | 1 gate | §2.2 |
| 6 | Eight genuine deduplications in C++ and R | 1 gate (several commits) | §4 |
| 7 | `Main.hpp` deletion + doc sweep | 1 gate + clean build | §5.1 |
| 8 | `SAIGE.Region` decomposition, 4–6 sliced commits | 4–6 gates | §6.1 |
| 9 | `mainRegionInCPP` type hoist (no phase lift) | 1 gate + clean build | §6.2 |

---

## §1 — Latent defects: ledger rows, not fixes

These cost nothing to record and are the highest-value output of the audit. Every one is behaviour today; none may be repaired inside a refactor.

### 1.1 The binary conditional region path still runs the pre-#104 burden-variance kernel (2/2 verified)

`R/SAIGE_SPATest_Region_Func.R:433-448` (`get_newPhi_scaleFactor`) is a textual twin of `get_newPhi_scaleFactor_traitType` (`:451-483`) sitting **18 lines above it**, differing only in the LEDGER #104 fix (`var_b <- sum(Phi)`; the `q_adj` standardisation at `:470-471`) and the `traitType` pass-through. The repaired one carries a 14-line comment naming its ledger row; the unrepaired one carries nothing. It is live: `R/SAIGE_Test_main.R:670`, inside `if (traitType == "binary" & isGroupTest)` — i.e. **every binary group test run with `--condition`**. `ledger.py grep` finds the unsuffixed name nowhere; #104's WHERE names only the `_traitType` function.

Two verifier corrections change what the row must say:

- The #104 identity is **not exact here**. `Score = m_Tstat_cond` is the *unweighted* `TstatVec` (`src/Main.cpp:4898`) while `Phi` and `q.sum` are dbeta-weighted (`src/Main.cpp:2070`, `:4896`). So `VarQ_2 = var1 · (ΣT)²/(Σ w T)²` — `r` is not #104's measured 0.94–0.95 and **may already be clamped to 1**, in which case the twin is inert on every golden.
- `r` **cancels** out of `Phi_cond`: through `R/SAIGE_SPATest_Region.R:604-615`, `VarMatAdjCond` is r-free and only `TstatAdjCond` keeps `r^(1/4)`. So the clamp reaches `Score_cond` and the `_cond` p-values, not `SE_Burden_cond`.

**Row content:** an unrecorded textual twin of a repaired kernel, live on the binary conditional group path; *and separately*, that call site feeds an unweighted `Score` against a weighted `Phi`, so `D = 1/w²` and `Q_b`/`VarQ` are on different scales. **Before any fix:** instrument `get_newPhi_scaleFactor` in the installed namespace and dump `var1, sum(Phi), VarQ, VarQ_2, r` for rs30/rs79 on golden P3. If `r == 1` today the row is a coverage statement.

**The byte-preserving structural half is available now and worth taking:** collapse to ONE function with the standardisation behind an explicit argument, passed `FALSE` at `R/SAIGE_Test_main.R:670`, with the ledger row number in the comment. Deletes ~15 duplicated inherited lines and turns a silent textual difference into a named one, with no golden re-baseline.

### 1.2 `.markerList.txt` loses every pre-restart row on a region restart (1/2, contested: the registry the finding proposed is refuted; the defect is not)

`R/SAIGE_SPATest_Region.R:1371` seeds `iswriteMarkerList <- rep(FALSE, nModels)` where `:1372` — the adjacent line, and exactly what LEDGER #59's fix added — seeds `iswriteRegion <- rep(isappend, nModels)`. On a restart `Start` is FALSE and the first resumed chunk with rows takes the `append = FALSE, col.names = TRUE` branch at `:2015-2025`, truncating the file. #140's truncation covers the results file and `.singleAssoc.txt` and is silent about this one. Upstream verbatim (`main:…:294`, `:1172-1216`). Untested: `test-step2-region-restart.R:39-43`'s loop names `.markerList.txt` but the run never passes the flag, so the loop is vacuous.

Candidate fix is one line (`rep(isappend, nModels)`), **not** the index-size route — that needs a third count in the index line and breaks `test-step2-region-restart-midchunk.R:90`'s adjacency assertion.

### 1.3 Stale sidecars survive `--is_overwrite_output=TRUE`

`Check_OutputFile_Create` (`R/Util.R:20`) removes only the results file. `.pooledTests.txt`, `.skatFailures.txt` and `.markerList.txt` from a previous run survive when the new run produces none, and `pooledSummary` writes only when at least one row is *interesting* (`R/PooledRow.R:121-126`) — so a clean re-run can leave a sidecar describing the previous run while the log says "SKAT failures: 0". Fix is three `unlink()` calls, not a registry.

### 1.4 Two indeterminate-`bool` reads on the live region path (2/2 verified)

`src/Main.cpp:2656` declares `bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;`. `is_Firth`/`is_FirthConverge` are **never assigned** in `mainRegionInCPP` (the three inner loops at `:2840`, `:3932`, `:4286` declare shadowing copies) and are read and passed by value at `:4573`/`:4575` and `:4618`/`:4620`; `is_FirthConverge` additionally converts `bool`→`int`. Benign only because the callee never reads those three parameters. This executes on every region run with `--is_single_in_groupTest=TRUE` (26 `.singleAssoc.txt` goldens). LEDGER #130's fix note records the *same trap* firing once already: `N_case_hom` derived from the identically dead function-scope `N_case`/`N_ctrl`.

This one is repaired by deletion (§3.3) rather than by a behaviour change, so it needs a comment at the deletion site, not a row.

### 1.5 The BURDEN header and its rows decide the `Weight` column from two derivations (2/2 verified)

`src/Main.cpp:5058` writes the BURDEN header from R's `nWeightSetsUsed`; `:2353` emits the rows' `Weight` column from C++'s own `q_weight`, composed at `:2318-2391` from argument *shapes*. Nothing reconciles them, and `t_weightNames` is read with a silent `"W"+r` fallback at `:3792`/`:5588` — the code already tolerates a size disagreement.

Verifier correction: three of the four composition arms are exact pass-throughs. The **one** live divergence surface is `hasCustomWeights = (!t_weight.is_zero() && t_weight.n_rows == q0)` at `:2317`: a region whose entire WEIGHT matrix is zero makes C++ drop the custom sets and fall through to one implicit beta set. With `nWeightSets ≥ 2` that is a **header with a `Weight` column and rows without it, at exit 0**. With `== 1` it silently swaps the set's arithmetic to beta(1,25) while printing the custom label.

**Actionable now, ~3 lines, moves no golden:** assert `q_weight == t_weightNames.size()` in `mainRegionInCPP` after `:2352` and delete the two `"W"+r` fallbacks. That converts a column-shifted BURDEN file into a loud error. Making C++ honour R's layout is a *behaviour change* on a reachable input and belongs on its own row.

### 1.6 Four more, stated briefly

- **`--minGroupMAC_in_BurdenTest` is accepted unchecked** (1/2, contested: the rest of the `checkArgs.R` finding is upstream-verbatim and out of scope). `R/checkArgs.R:57` takes the formal; no line checks it; the body's only mention is a commented-out block at `:75-77`. LEDGER #126 made it load-bearing on both region paths. A non-numeric coerces to `NA`; `R/SAIGE_SPATest_Region.R:509-510` NA-hardens `cellMAC` but not the threshold, so `gateOK` is `NA` and the run dies in `if (NA)` with no argument name — the exact LEDGER #21(d) shape the leaf was hardened to prevent. One-line fix.
- **`getNeff(traitType="quantitative", noEstFixedEff=TRUE)` errors.** `covoffset` is assigned only in the binary arm (`R/SAIGE_extractNeff.R:353`); the quantitative arm replaces the phenotype instead and never sets it, yet `:446` reads `offset = covoffset`. Unexercised since LEDGER #43 first routed the flag there.
- **`getNeff`'s LOCO window block computes indices in the wrong space.** `R/SAIGE_extractNeff.R:244` uses the *raw* `readVariantChromFromInput` and runs before `setgenoFromInput`, where `fitNULLGLMM` subsets to QC-retained markers (`R/SAIGE_fitGLMM_fast.R:125-129`) because the C++ indices live in QC'd-marker space. Harmless only because nothing reads the result. Also missing the LEDGER #62 contiguity guard and hardcoding `plinkFile, ".bim"` in a message reachable from a PGEN run.
- **`writeOutputFileIndex`'s `singleBytes` parameter has no caller** (2/2 verified). LEDGER #140 correctly extended the canonical index writer at `R/Util.R:355-357` and then the region path open-coded the protocol instead (`R/SAIGE_SPATest_Region.R:1933-1936`, `:2031-2067`), re-stating `"This is a Region level analysis."` and `nEachChunk = 1` as literals against `checkOutputFile`'s independently-built versions (`R/Util.R:195-199`). The comment at `R/Util.R:351-353` is now stale. Divergence surfaces as *"probably it has been modified by user"* on the next restart, and `*.index` is excluded from the byte gate.

---

## §2 — The three controls that do not control

The golden gate is a byte oracle over one fixture. It is structurally blind to dead code, uninitialised reads, unreachable-but-exported paths, and cross-language divergences whose halves are mutually exclusive. Everything in this section buys back detection power the gate cannot provide.

### 2.1 `tsan.sh` aborts after the first arm that dies — in the mode where dying is expected (2/2 verified) — **do this first**

`harness/tsan.sh:58` sets `set -euo pipefail`. The died-arm branch ends:

```
    [ "$MODE" = clean ] && DIED=1
    return
```

Under `plant`/`plantreader` the `&&` list has status 1, the bare `return` propagates it, and `tsan_run` is called as a bare command at `:296`, so errexit kills the script. The operator sees one arm line, no summary, exit 1 — reading as "the harness broke", not "the plant worked and t=4 and t=8 were never run". `plantreader` is capped at 2000 markers precisely because it is pathologically slow and would otherwise hit the timeout (status 124, same branch).

CLAUDE.md: *"every macOS setup mistake fails by reporting NO races"* — the plant is the positive control for every `clean` result on this branch, including the `races=0` baseline the region workload's fork-edge annotation rests on.

**Fix:** `if [ "$MODE" = clean ]; then DIED=1; fi` + `return 0`, and make the other early return at `:284` explicit. Reject the `harness/lib.sh` coupling the finding proposed — it does not exist and `guardmalloc.sh`'s `verdict()` is a different shape (one expected status per call site; under `plant` *both* outcomes are licensed). Not "behaviour-preserving" — it is golden-safe and gate-free, which is the honest framing.

**Drill both ways, cheaply:** `SAIGE_TSAN_WORKLOADS=binary SAIGE_TSAN_TIMEOUT=5 SAIGE_TSAN_MAX_MARKERS=200 bash harness/tsan.sh plant`. Unpatched: one arm, exit 1. Patched: three arms, the guidance block, exit 0.

### 2.2 `src/Makevars` compiles 22k lines of C++ with every diagnostic off (2/2 verified)

`src/Makevars:1` carries `-w` beside a `-pedantic` it silences. R's own `CXX14FLAGS` on this platform is `-falign-functions=64 -Wall -g -O2`, and R places `ALL_CPPFLAGS` first — so `-w` actively *cancels* the `-Wall` R passes by default, and `-w` is order-independent in this clang, so appending `-W…` elsewhere would silently do nothing. `harness/build.sh:89` writes `harness/work/build.log` and greps it for `error` only on failure. Nothing reads a compiler diagnostic at any point.

Measured with the pinned clang, R's exact flags, `-w` removed, `-Wall -Wextra -Wno-unused-parameter`: **306 unique diagnostics, 270 src-anchored, 139 in hand-written code.** Src-anchored histogram: `unused-variable` 65, `sign-compare` 37, `deprecated-declarations` 16, `unused-but-set` 7, `sometimes-uninitialized` 6, `bitwise-instead-of-logical` 3, `uninitialized` 2, plus `cast-function-type-mismatch` 131 (all generated `RcppExports.cpp`). The `deprecated-copy*` family is header-origin, not src.

**Codegen identity verified, not asserted:** `SPA.o`, `SAIGE_test.o`, `CCT.o`, `UTIL.o` at `-O2` are `cmp`-identical built with `-w` versus with `-Wall -Wextra`.

The 8 uninitialised-use sites are the payoff: §1.4's pair, the two dead `Unified_getSampleSize*` returns (`:1727`/`:1736`), and — the ones nobody knew about — `src/SAIGE_test.cpp:377` (`var2`, `Bmu2`) and `:1151` (`q_c`, `qinv_c`, passed into `SPA_ws` at `:1161` when `m_traitType != "binary"`), both on the **live SPA / conditional-analysis path**. Those two are LEDGER material, not drive-by fixes; a repair there can move goldens.

**Corrected remedy.** Local: drop `-w`, keep `-pedantic`, add `-Wall -Wextra -Wno-unused-parameter -Wno-cast-function-type-mismatch`; have `harness/build.sh` count `warning:` lines matching `^src/` **excluding the generated `RcppExports.cpp`** and compare against a committed baseline as a ratchet. **Do not** put a count ratchet in CI — CI is Ubuntu/conda-GCC against a macOS/conda-clang baseline, `pixi.lock` is uncommitted and every dep is `"*"`, so the compiler is not pinned on either side and warning spellings differ. In CI instead promote only bug-shaped categories: `-Werror=uninitialized -Werror=sometimes-uninitialized -Werror=return-type -Werror=parentheses`. Do **not** let the first green push "fix" `SAIGE_test.cpp:1000`'s `-Wparentheses` — `m_is_Firth_beta & pval <= m_pCutoffforFirth` is correct.

This adds a gate; it does not delete complexity. Say so.

### 2.3 `prepush.sh` records a run shape it did not execute (1/2, contested: the knob is *not* dead)

`harness/prepush.sh:102` logs `concurrent=$CONCURRENT` on **both** paths. With the default `MODES=both`, `set -- $MODES` leaves `$# = 1`, the single-tier branch runs, and `SAIGE_PREPUSH_CONCURRENT=1` records `concurrent=1` in `harness/work/prepush-both.log` for a run in which nothing ran concurrently. This is the same shape as the defect this file's own header (`:45-50`) records and warns about — a gate log line that named the wrong thing and was then quoted as evidence.

**Contested and refuted:** the knob is *not* dead. `MODES="byte numeric" CONCURRENT=1` reaches `:165`; concurrency is by construction a property of running two tiers, and the retention decision is recorded in three places (`:88-90`, `harness/bench/NOTES.md:1946`, `docs/state/handoff-2026-09-12-review-fixes-and-gate.md:72`). Do **not** delete it.

**Take:** pass the fact into `run_tier` so it logs what it *did*; warn when `CONCURRENT=1` arrives with one mode; extract the byte-identical 8-line dirty-tree block at `:153-160`/`:192-199`. ~12 lines, no gate.

---

## §3 — The dead-surface sweep

The repo already ran a whole-repo dead-surface audit on 2026-09-08 (`docs/state/next.md:1333-1375`) and recorded the disposition: *"size, not behaviour … left as a tidiness item because deleting it moves no golden and answers no question."* Five of the refuted findings were refused on exactly that ground. **Frame this sweep as executing that recorded audit, not as a discovery** — and note what the audit's method structurally could not see:

- Its rule was "exactly one occurrence of the name in `src/`", which an `[[Rcpp::export]]` can never satisfy (its name also appears in `Main.hpp` and twice in `RcppExports.cpp`). So its conclusion *"all 17 are in `SAIGE_fitGLMM_fast.cpp`"* is **wrong at HEAD** and was wrong when written — at minimum `Unified_getSampleSizeinGeno` and `Unified_getSampleSizeinAnalysis` (callerless since the fork import `6867bd26`) belong in that bucket, as do `closeGenoFile` and `Unified_getOneMarker`.
- It scanned functions, never **parameters**, never **post-`stop()` tails**, never C++-internal (non-exported) symbols, and never a translation unit absent from `OBJECTS`.

Amend the audit table in `next.md` as part of this work.

Six commits, one push, one gate. Predict *no golden moves* in writing first.

### 3.1 `Binary_*`: ~1,305 lines, ~495 of them linked into the shipped `.so` (2/2 verified)

Only `SKAT_Exact` (`src/Binary_global.cpp:69-80`) is live, reached from `ER_binary_func.cpp:247`. Dead: `SKATO_Exact` (`:83-98`, sole instantiator of the 286-line `ComputeExactSKATO`), `SKAT_Permu` (`:127-139`, sole instantiator of the 209+56-line `Binary_Permu_SKAT`), `GetProb` (`:118-124`, a verbatim duplicate of the live `GetProb_new`), the commented blocks at `:100-116` and `:141-189`, `SL_GetPermu`, `Test1`, `SL_runif_int`, `SL_out`, `mat_to_std_vec`. Confirmed at the symbol table: `nm -gU` on the installed `.so` returns all of them as defined exports.

**`src/Binary_ComputeExactMC.{cpp,h}` (595 lines) is an orphan**: git-tracked, `#include`d at `Binary_global.cpp:6`, **absent from `src/Makevars:10`'s OBJECTS**, and its `Binary_ComputeExactMC.cpp:8` includes `"Binary_global.h"` — a header that does not exist. That is a latent build trap, not a size item: anyone who uncomments `SKAT_ExactMC` gets a link error, anyone who edits those 493 lines edits a file that never compiles, and any attempt to let R wildcard the sources (§3.5) fails outright.

**Constraints.** Do **not** merge `SL_Sample` with `SL_GetSample` — identical bodies, *swapped parameter order*, called with different intents (`SL_Sample(k, m_k, …)` vs `SL_GetSample(n, n, …)`), and the live ER resampling arm is byte-gated by golden P33. Make the `CalTestStat`/`CalTestStat_INV` de-virtualisation a **separate optional commit** and drop its performance justification (at `--max_MAC_for_ER=4` the recursion enumerates k ≤ 4; the indirect call is noise, and there is a fifth call site at `Binary_ComputeExact.cpp:270` the finding missed). Rewrite the S2 comment at `Binary_ComputeExact.cpp:383-386`, which cites `SKAT_ExactMC` and `SKAT_Permu` by name as part of LEDGER #2's byte-neutrality argument, and re-run `docdrift.py`. Deleting `SL_runif_int` shifts `SL_setseed` off LEDGER #2's cited `src/Binary_resampling.cpp:49-56` — either leave `Binary_resampling.cpp` alone or update the row.

### 3.2 The `Unified_*` bridge and `t_genoType` (2/2 verified; the R-facing half contested)

`Unified_getOneMarker` (`:1680-1707`), `Unified_getSampleSizeinGeno` (`:1725-1731`), `Unified_getSampleSizeinAnalysis` (`:1734-1740`), `Unified_getMarkerPval` (`:1746-1786`, 30 parameters, **not** exported, zero callers) and `closeGenoFile` (`:5032-5040`, exported, no caller — R uses `closeGenoFile_plink`) are ~90 lines with no caller anywhere.

**Correction that lowers the severity:** the fallthrough is *unreachable*. `Geno::makeReader` (`src/GenoReaderRegistry.cpp:182-196`) `Rcpp::stop`s on an unregistered type **and** re-checks `reader->formatName() != genoType` against the constructed object; `R/Geno.R` can only produce `"plink"`, `"pgen"` or `""`. So three of the four `t_genoType ==` guards sit in dead functions and one (`:372`) is live with four call sites; the two indeterminate `uint32_t` returns are reachable only from an R console.

**Take now:** delete `Unified_getMarkerPval` and its `Main.hpp:70` declaration, plus the three comments that exist only to steer readers away from it (`src/Main.cpp:921`, `:2842`, `:4861`). No export, no `compileAttributes`.

**Owner's call:** the three `[[Rcpp::export]]`ed ones are public symbols via `exportPattern`. Deleting them is a public-API removal of the same class #114 chose to RETAIN. Hand them over with the audit-method defect attached.

**Do not bundle the `t_genoType` parameter removal.** It churns the R-facing signatures of `mainMarkerInCPP`, `mainRegionInCPP` and `assign_conditionMarkers_factors` for zero behavioural gain. If it is ever done, justify it by the **latent defect** rather than line count: `mainMarkerInCPP` passes ONE `t_genoType` to `testOneMarker` for EVERY reader (`src/Main.cpp:1441-1442`) while `addGenoReaderInCPP` accepts plink and pgen independently (`R/MultiModel.R:169-182`) — the guard can be testing a string that does not describe the reader being read. Inert today only because both formats take the same branch. Docs that must move: CLAUDE.md's "four `t_genoType ==` guards", `docs/state/multimodel-step2.md:85`, `docs/flag-surprises.md:45`, `harness/bench/writer_profile.py:29`.

### 3.3 Region dead declarations and the UB parameters (2/2 verified)

`src/Main.cpp:2493-2495` and `:2653-2671` name 37 function-scope declarations, of which **35 are dead and 2 are read uninitialised** (§1.4). Four of them (`P1Vec`, `P2Vec`, `GVec`, `GZeroVec`) are `arma::vec(t_n)` allocated and freed on every region call. `numofUR0` (`:4561`) is assigned at both writer call sites and never read; `int mFirth = 0` (`:4562`) exists only to fill a slot.

Delete those, plus `t_isFirth`/`mFirth`/`mFirthConverge` from `writeOutfile_singleInGroup` (`src/Main.cpp:5648-5650`, `src/Main.hpp:263-265`) — grep of the body `:5683-5783` returns them only in the signature. Leave `writeOutfile_single`'s three alone; they are live at `:5525-5526`. Note that `numofUR` reaching R is **not** a miscount: R wants the per-cell slot count (`R/SAIGE_SPATest_Region.R:314`), which is what `:4560` supplies. Land a one-line comment at each deletion site in the style of `:3823` so the next upstream merge does not re-import the block; frame the commit as removing shadow traps in a threaded function, not as a UB fix.

Also delete `bool hasVarRatio = true;` at `src/Main.cpp:2671`: `e72e4186` gave the three inner region sites their own shadowing locals and left the frame-scope declaration standing — that is the exact variable TSAN caught (5 races at t=4), still in the frame, still armed for the next threaded prepare added there.

### 3.4 SPA: the guard with no call site and the two half-length gathers (2/2 verified)

`prepSPAContext` (`src/SAIGE_test.cpp:750-772`) and `spaContextReady` (`:749`, `:798`) have no call site. `git log -S"prepSPAContext();"` shows `8f3ded00` added the call for LEDGER #40 and `b60141f9` (the #84 residualisation) deleted it and left the lambda plus its 12-line "#40 is fixed by this" comment behind. The conditional block has computed its own CGF since (`:1141-1161`).

This is the most expensive kind of wrong documentation: the guard *reads* as armed, complete with a ledger citation and a measured reachability count. It is not armed. (Mitigating: goldens P12/P17 pin the #40 path, so a reroute would likely move bytes.)

Separately, `gNA` (`:827`) and `muNA` (`:829`) are `n/2`-length subview copies built per SPA-fast marker and read by nothing — verified through `SPA_fast_ws` → `getroot_K1_fast_Binom_core` / `Get_Saddle_Prob_fast_Binom_core`. **This is inherited** (upstream's bodies do not read them either), so **do not touch any exported signature**: drop `mu`/`gNA`/`muNA` only from the three *internal* functions (`SPA_fast_ws`, `getroot_K1_fast_Binom_core`, `Get_Saddle_Prob_fast_Binom_core`) and stop materialising at `:827`/`:829`. Keep `:721-726` — the lazy-sizing paragraph is a live measured Phase-4-C3 rationale recorded nowhere else. Keep the two `clear()` calls but fix `gNB.clear()` ×2 → `muNB.clear()`. **Append** to LEDGER #40's STATUS naming `b60141f9`; do not rewrite it.

**Reject** the `isScoreFast`/dead-PCG half: it is upstream verbatim and the `else` at `:1241-1245` is the structural twin of the live `else` at `scoreTest:263-272`. Reduce it to a comment at `:1236`.

### 3.5 `getMem`, `UTIL`, and the `OBJECTS` list (1/2, contested on scope)

`src/getMem.{cpp,hpp}` (52 lines) defines `process_mem_usage`, whose only references are four commented-out call sites. It reads `/proc/self/stat`, absent on this branch's primary platform, and on failure returns two **uninitialised** locals, not the 0.0 its header comment promises. Delete it, **all five** `#include "getMem.hpp"` lines — `CCT.cpp:5`, `SAIGE_test.cpp:16` **and** `:17` (duplicated), `SAIGE_fitGLMM_fast.cpp:17`, and `src/Main.cpp:37`, which the original finding missed — and the `getMem.o` entry.

In `UTIL`, delete `imputeGeno`: its declaration (`UTIL.hpp:19`, 3 params) **does not match its definition** (`UTIL.cpp:33-56`, 4 params), so the declared overload has no definition and the defined one no declaration — survivable only because nothing calls it. It is a second implementation of the imputation rule with an incompatible vocabulary (`"bestguess"`/`"none"` vs the live `"best_guess"`/`"minor"`). **Repair the line anchors in the same commit** — deleting it shifts `imputeGenoAndFlip` and stales LEDGER #28 (`:74-79`, `:117-122`), #109 (`:108-112`), `tests/testthat/test-step2-missing.R:109`/`:251`, `src/PLINK.cpp:508`, `docs/flag-surprises.md:30`/`:165`.

**Leave the rest of UTIL alone** — `getWeights`, `getInnerProd`, `getinvStd`, `nb`, `whichCPP`, the five unread `PlinkClass` accessors, and especially `imputeGenoAndFlip`'s signature (the most-ledgered function on the branch; `t_MAC`'s mutation is *not* dead — it is read six lines later by #109's dosage-zerod gate).

**`OBJECTS` (2/2 verified):** `src/Makevars:10` names 20 objects against 21 sources — the difference is exactly the orphan from §3.1. Order matters: **delete the MC pair first**, then optionally drop the `OBJECTS` assignment and let R derive it. And verify the second step rather than trusting it: R's `share/make/shlib.mk` expands `$(OBJECTS)` behind `@if test "z$(OBJECTS)" != "z"`, so an empty list **silently skips the link** and leaves a stale `.so` that passes the gate while proving nothing. Check the link command in the build log and `ls src/*.o | wc -l` = 20. Drop the `Makevars.win`/OpenMP commit — inherited, already dead, and it *adds* a "Windows is unsupported" claim rather than deleting one.

### 3.6 R-side dead surface

- **`R/Firth.R` — 450 lines, 9 functions, zero external references** (1/2, contested: the NAMESPACE attack is refuted). It arrived with `5ef29751` "Add RareEffect"; `376e2d01` stripped RareEffect and left this behind. Live Firth is C++ twice over (`SAIGEClass::fast_logistf_fit_simple`/`_closed`), and C++ never calls back into R except `base::set.seed`. It is a **third** implementation of penalised logistic regression with a different convergence rule, a different `maxit`, an `l2` penalty the C++ has no equivalent of, and a `# HERE IS THE PROBLEM!!!` comment at `:65`.
  **The claimed blocker does not exist.** `exportPattern` exports whatever symbols exist; deleting the file removes them with no `NAMESPACE` edit, and nothing asserts the export set. **Do not narrow `NAMESPACE`** — 19 `SAIGE::<internal>` references across three test files plus five `library(SAIGE)` child snippets calling ~22 bare names mean the real export surface is ~35+, not 6. **Do not touch `checkArgsBool`/`checkArgsNumeric`** — `test-unit-util.R:109` is LEDGER #15(a)'s inverted characterization test, and `test-unit-checkargs.R:6-9` is a written warning against exactly that consolidation.
  Also delete `ReadGroupFile` (`R/SAIGE_SPATest_Region_Func.R:3-60`; its own error branch at `:48` formats an undefined `n2`, which is proof no input ever reached it) and `check_close` (`:63-67`). `pcg`/`pcgSparse`, `bigGRMPar`/`bigGRMPar_new`/`refineKinPar`, `Util.R::writeOutputFile` may ride along; note `docs/cli-options-and-group-file.md:199` explains step-0's inert `nblocks`/`ncore` by reference to `bigGRMPar`, so that doc line moves with it.
- **`Get_Variance_Ratio`'s post-`stop()` tail** (2/2 verified). `R/readInGLMM.R:443-447` is the unconditional LEDGER #13 `stop()`; `:448-466` is 19 unreachable lines including two WARNING strings describing a format the function refuses. Two formals (`isGroupTest`, `useSparseGRMtoFitNULL`) appear nowhere in the 118-line body, and the sole production caller (`R/MultiModel.R:257-260`) passes five arguments — so `useSparseGRMtoFitNULL` is a never-forced missing promise, verbatim the shape LEDGER #85 fixed by dropping the formal. Add the mirror of #85's gate (`expect_identical(names(formals(Get_Variance_Ratio)), …)`) and fix `test-unit-varianceratio.R`'s header comment, which asserts a six-argument contract production never uses. **Keep** the `print("ratioVec_sparse")` block — it is live CLI stdout, so removing it is a behaviour change for another commit (and no test reads it, contrary to the original justification).
- **`printComb(3)`** at `R/SAIGE_fitGLMM_fast.R:1800` (1/2, contested: LEDGER #63 deliberately left the printer). Take only the single **R call line**: it prints two junk lines (`i,j 0,0` — a diagonal the real worker never emits — and `i,j 0,1`) on every step-0 and internal-GRM run. Deleting the call needs no `compileAttributes` and leaves #63's description of the function true. Fold the rest (the `[1:102]` dumps, the `numeric_limits`/`x = 1000001` block, the three dead R kinship functions) into the batched sweep, keeping `src/Binary_ComputeExact.cpp`'s `totalCombination:` line, which `test-characterization-ledger63.R:126` asserts.
- **`CarrierStore::urFill`** (1/2, contested: already adjudicated). LEDGER #130's STATUS says it and the dense clear are *"kept as a removal candidate like #114's arms"* and the handoff says removing them is *"not a passing change"*. If the owner declares this sweep to be the audit those defer to, the remedy misses a fifth `clear()` site at `src/Main.cpp:3695` and two dangling comment references (`:571`, `:3376`), plus `CLAUDE.md:460` and a LEDGER #130 STATUS amendment. No region TSAN run is required — the guard predicate is compile-time dead, not a threading change.
- **`region_postprocess`'s dead return field** (2/2 verified, narrowed hard). `ret$resultDF` (`R/SAIGE_SPATest_Region.R:870-872`) has no reader anywhere, and the unpack at `:282-284` exists only to make `exists()` true. Delete both plus the `resultDF` half of the comment at `:125-127`, and collapse `:291-293`/`:360-362` into one unconditional `annoMAFIndicatorMat <- outList$annoMAFIndicatorMat`. **Withdraw the `is_rerun` half entirely** — both branches are reached, the `:112-129` comment is accurate past-tense provenance, and `run_core <- is_rerun` is *not* redundant because the caller's predicate is computed over `S$outLists[[1]]`'s fast-pass counts while the function's would be over the model-k outList `:1546` has just replaced. **Keep** the two `%in% names(res)` guards in `takePost`: they are the same class of protection as the `[k] <- list(...)` idiom above them, whose comment records a real sim50k model-row mis-binding, and a live mock at `test-unit-skato-batch.R:44` returns only `pval.Region`.
- **`extractNglmm.R` / `getNeff`** (1/2, contested: refuse the deletion). 15 of `getNeff`'s 55 formals are never read; four are live CLI flags (`--tol`, `--maxiter`, `--tauInit`, `--numRandomMarkerforSparseKin`), parsed, echoed by `print(opt)`, documented in `man/getNeff.Rd`, and forwarded into a function that fits no variance component (`tauVec_ss <- c(0,1)` hardcoded). `--tauInit` is even validated and printed at `:347` before being discarded. **Do not delete** — `getNeff` is exported, `man/getNeff.Rd`'s 55 `\item`s are LEDGER #93's measured invariant, and `docs/cli-options-and-group-file.md:196` records the branch's chosen treatment for inert flags: keep them declared and give them a status. **Do:** a LEDGER row; correct `extdata/extractNglmm.R:7-8`, which says the script *"stops after the variance-component fit"* that never happens, with a `docdrift.py` rule so the twin cannot come back; and extend the cli-options inventory to a **fourth** wrapper (its 2026-09-12 sweep covered 139 options across three scripts and skipped this one, which is why the whole surface has no recorded status).

---

## §4 — Duplicated decision rules

These are the ones that have already produced ledger rows. Each is a real deletion; each is one commit.

### 4.1 The per-cell "how is this tested" rule, six times (2/2 verified)

The same derivation — `hasVarRatio`, `ctx.useSparseGRM`, `ctx.noadjCov`, `select{Single,}VarianceRatio` — appears at `src/Main.cpp:952-961`, `:1013-1020`, `:2828-2837`, `:3747-3754`, `:4206-4213`, `:4808-4815`. The ER predicate beside it has **five** sites in three spellings (`:992`, `:2844`, `:3756`, `:4215`, `:4860`), and `isSingleVarianceRatio` is derived in three places with two spellings that disagree at `n_elem == 0` (`== 1` at `:1373` and `:4812`, `<= 1` at `:2571`).

Drift here is the documented cause of **three** rows: LEDGER #134's STATUS says in as many words that *"#118 fixed one of the pair and missed the other"* — for two months a custom weight set's *scale* chose its variance-ratio band (5 alleles read as 38.4 under CADD, 3.34 under REVEL) — and #135 then fixed a third divergence at the same three decision points.

**Corrected remedy — four sites unify, not six.** Add `SAIGE::TestContext SAIGEClass::contextForMAC(double) const` reading only `m_cateVarRatioMinMACVecExclude`, `m_flagSparseGRM`, `m_varRatio_null`; call it at `:1013`, `:2830`, `:3748`, `:4206`. **`:950` is a different rule** (`m_isFastTest ? false : m_flagSparseGRM`, `noadjCov = m_isnoadjCov`, no MAC term). **`:4808` must keep its own derivation** — routing the conditioning-marker path through `contextForMAC` would *add* the `MAC > back()` gate and move golden P4's `_c` columns.

`isERCandidate` **cannot** be a `SAIGEClass` method: `g_cfg` is `static Step2Config g_cfg;` at `src/Main.cpp:123`, file-scope internal linkage invisible to `SAIGE_test.cpp`. Make it a free `static` beside `g_cfg`, and never re-spell `cm.isER` at `:3945`/`:4296` — it carries an extra `erAllowedVec[r]` conjunct and reads `MAC_nw`.

**Unify on `<= 1` *and* add the guard that makes the question moot.** `src/Makevars:1` sets `-DARMA_NO_DEBUG`, so an empty `m_varRatio_null` does not throw — `varRatio(0)` reads off the end silently. Reject a model whose `null` rows are absent, in `Get_Variance_Ratio` or at the top of `contextForMAC`. That is the part of this change with user-facing value.

Bundle §4.2 with it.

### 4.2 `selectVarianceRatio`'s `hasVarRatio` out-parameter (2/2 verified)

`src/SAIGE_test.cpp:1261` takes `bool & t_hasVarRatio`; the only path on which it stays false **throws** at `:1298-1307`, so on every normal return it is `true` and carries no information. It is written at 14 sites and **read at zero**. LEDGER #32's WHAT cell says so. The cost is recorded: the comments at `src/Main.cpp:2825-2827` and `:3714-3721` document that this dead bool produced a real worker-vs-worker race (TSAN, 2026-09-06), fixed by giving each worker its own dead bool.

Correction: it *is* read three times **inside** the function (`:1283`, `:1290`, `:1297`) to sequence the fallback guards — which is why a shared copy was a live race rather than a benign one. So the remedy is to **demote it to a local**, deleting the parameter, the 14 caller-side declarations and the whole "shared out-flag written by every worker" class by construction. Two comment obligations: the per-site #134 (`:3722-3746`) and #135 (`:4185-4195`) notes describe which *argument* is passed and stay at the call sites; and `src/Main.cpp:3186` justifies the declared TSAN fork edge by citing `SAIGE_test.cpp:1274`/`:1279`, which after deletion name nothing — reword it to cite `d36e7277` (racing) vs `e72e4186` (clean).

**Gate:** full both-modes **plus** `SAIGE_TSAN_WORKLOADS=region bash harness/tsan.sh clean` and a plant arm; this touches the region prepares' frames.

### 4.3 One single-variant row schema, two writers (2/2 verified)

`writeOutfile_single` (`src/Main.cpp:5424-5529`, over `Marker::MarkerResult`) and `writeOutfile_singleInGroup` (`:5645-5783`, over parallel vectors) emit the identical schema — same column order, same `if(pval != "NA")` skip, same `writeNumOrNA`, same `std::boolalpha`, same more-output order. The header emitters at `:5080-5128` and `:5135-5193` duplicate the column emission against two streams. Nothing ties them: no test compares the two files' column sets, and the goldens pin each separately, which is exactly why a fix to one leaves the other's golden green.

**The hand-sync failure is documented, not inferred:** `292b5744` (LEDGER #44) added `writeNumOrNA` to the marker writer; `c9f7e9ae` (#66) added it to the region writers a day later; `docs/state/handoff-2026-09-05-phase3-start.md:167` records it as *"a gap LEDGER #44 left"*. The helper now has **three** hand-synced sites (`writeOutfile_BURDEN:5594-5603` too).

**Corrected remedy — keep the region path's SoA.** One `writeSingleVariantRow(std::ofstream&, const SingleVariantSchema&, const Marker::MarkerResult&)`; inside the shared row writer, fill a local `MarkerResult` from the vectors at index `k`. That deletes the ~250-line duplicated emission body and the duplicate header writer, keeps the `t_isImputation` exchange in one place, and touches neither `RegionModelWork`, `accumulate`, the `OutList` contract, nor memory. Three constraints: the unified function must return **both** counts (the marker path prints `numtest` from inside the same loop at `:5533`); the schema's trait type must keep coming from **R** for the in-group writer (`R/SAIGE_SPATest_Region.R:1067`) rather than silently switching to `model.m_traitType`; and folding the four `t_isImputation` branches into two is **LEDGER #108's owner call**, not a self-granted reading.

The signature is 39 parameters (`src/Main.hpp:260-297`), and the parameter order (`N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom`) **differs from the write order** — a live hazard across two 39-argument positional call sites. Land §3.3's dead-parameter deletion first; it is independently correct and provably inert.

### 4.4 The p-value formatter, six times (2/2 verified)

`scoreTest`'s tail (`src/SAIGE_test.cpp:288-344`), `scoreTestFast`'s (`:398-453`) and `scoreTestFast_noadjCov`'s (`:512-568`) are 42 code lines each; stripped, the two fast ones are **identical** and `scoreTest` differs in exactly two tokens. The 13-line LEDGER #44 rationale is pasted verbatim three times. The formatter ladder (`log10p` → floor exponent → `fraction >= 9.95` → `sprintf("%.1fE%d")`) appears at `:315`, `:424`, `:539`, `:930`, `:1100`, `:1208`.

The printed spelling of every step-2 p-value is this branch's golden contract, and it is written six times. LEDGER #42 was one defect fixed at three sites in this function; #44 at three more. **The log-p branch is covered by nothing** — no golden contains a one-decimal-exponent p-value and no test names the 9.95 rule, so a missed edit at one of six would move no golden and fail no test.

**Take the two formatters only.** `formatPvalOrLog` at `:300`, `:411`, `:524`, `:1086`; `formatLogP` for the else-halves at `:917-931` and `:1192-1209`, leaving each site's `%.6E` line in place because `:1192` carries LEDGER #41's comment.

**Reject `finishScoreTest` as specified.** The three p-value guards are not the same predicate: `scoreTest:288` is `if (collapsed)` where `collapsed` includes `projectionCollapsed`; the two fast paths test `var1 <= DBL_MIN` alone. Passing `pvalCollapsed = false` would **delete the absolute guard** from both fast paths — unreachable on this fixture (P30 sits ~294 orders of magnitude above `DBL_MIN`), which is exactly why CLAUDE.md forbids silently removing unreachable guards. If taken anyway, the helper must apply the `DBL_MIN` test unconditionally and OR the caller's predicate on top, and must not serve the conditional block at `:1063-1120` (no variance-ratio multiply, forces `stat_c = 0` inside its guard, uses `sqrt(stat_c)` where the others use `sqrt(fabs(stat))`).

**Cheapest useful version, if nothing else:** replace two of the three 13-line #44 comment copies with `// LEDGER #44 FIXED -- rationale at scoreTest(), src/SAIGE_test.cpp:325.` That deletes ~26 duplicated lines, removes the must-be-corrected-three-times hazard, and changes zero instructions.

### 4.5 The BURDEN emitter's six loops (2/2 verified)

`src/Main.cpp:2354-2390` builds four parallel arrays over one index — `setKindVec`, `setColVec`, `betaDistVec` (indexed by `setColVec[r]`, **not** by `r`), `erAllowedVec` (whose beta and No_Weights entries rely on `vector<bool>(q_weight, true)`'s constructor default to mean "ER allowed"). The kind is re-dispatched at six sites. `:3083` is a provable identity (`pr.w0(r)` on both arms). `:3628-3684` is three kind branches × two `q_weight == 1` branches = **six copies** of one accumulation loop.

**Scope honestly:** the per-kind table *is* written down (design doc D2, and `src/Main.cpp:2306-2309`), and the "transform behaves like custom" invariant *is* gated (`test-step2-flexrv.R`, `test-step2-flexrv-reconcile.R`). Two ladders relocate rather than vanish, and a sixth site (`:2933`, the transform-only export) survives. What genuinely disappears: four parallel arrays with three index domains → one file-local `struct WeightSet { int kind; unsigned col; boost::math::beta_distribution<> dist; bool erAllowed; bool collapseWeighted; }` with `erAllowed` explicit for all four kinds; the six emission loops → one; the identity ternary; the dead `arma::vec w0_vec(q_weight)` at `:2392`.

**Declare the struct file-local in `Main.cpp`**, not in `FlexRV.hpp` (which would pull `boost::math` into every includer, and FlexRV is one *kind* of this table, not its owner). Bind the count source as `const arma::vec &` chosen once, not assumed equal to the collapse. Hoist the `q_weight == 1` predicate so the compiler unswitches — this loop runs 210 × 192 on sim50k.

**Mandatory with the commit:** carry the whole `#17`/`#118`/`#129`/`#130`/`#133`/`#134` comment block at `:3530-3684` forward at the call site — it is the authoritative provenance for five FIXED rows — and update `src/Main.cpp:3531` and LEDGER #130's WHERE cell, both of which say "six emission loops". Stale *line numbers* are tolerated here; a stale *structural* description is not. State in the commit message the two facts byte-identity rests on: `x * 1.0 == x` exactly for every finite `x`, and the count source is **selected**, not assumed.

### 4.6 Three region records that are one record (2/2 verified)

`URCellModel` (`src/Main.cpp:2518-2527`) and `BurdenCellModel` (`:2528-2536`) are field-for-field identical but for `arma::rowvec G1`; `RegionPrepModel` (`:2699-2711`) repeats the same 14 result fields. The write-back is duplicated too: `accumulate`'s at `:2949-2965` and pass B3's at `:3973-4002` assign the same 15 slots with the same `(1 - 2*flip)` corrections.

**Use inheritance, not composition** — `struct URCellModel : CellResult` leaves every existing `prm.Beta` / `cm.pval` spelling untouched, where composition forces rewriting ~114 member accesses, most of them inside six ~600-character `getMarkerPval` calls that pass 19 out-parameters each. **`G1` must not go into the base:** `URCellModel::G1` holds the raw `G1tilde_P_G2tilde_Vec` (`:3954`) while `RegionPrepModel::G1row` holds the already-weighted row (`:2855`) — one field name carrying two quantities is the shape of #118/#133/#135. Keep the shared writer file-local in `Main.cpp` unless a `RegionLoop.hpp` actually lands, and carry #130's comment explaining that the ultra-rare path's `(1 - 2*flip)` is a deliberately preserved no-op.

**Do not** unify `Marker::MarkerResult` into this: it stores post-flip values, has no `gy`, and defaults `pval` to `"NA"` where the region records default to `""` — and the writer keys on that string, so a shared default would emit rows for QC-dropped markers.

### 4.7 The ER if/else, five times (2/2 verified — take only this slice)

`src/Main.cpp:993/998`, `:2845/2847`, `:3946/3948`, `:4297/4299`, `:4866/4875` are byte-identical 30-argument call pairs differing in one boolean literal; `:4866/:4875` repeats the same five-line comment in both arms. Twelve full-argument call sites in total.

Collapse the five branches by hoisting each site's **existing** predicate — `cm.isER` verbatim at `:3945`/`:4296` (it carries `erAllowedVec[r]` and `MAC_nw`; re-spelling it sends weighted cells into ER and moves P58-P61). Leave `:1022` alone, and **do not** fold the escalation guard at `:1005` into `!isER` — the two expressions agree only because `MAC` is finite here.

**Reject** the 13-field `TestOutcome` struct and the `conditionalMarkerTest` extraction: the extraction's frame-locals inventory is wrong (`ispvallog` is written at `:1090`/`:1101` and read at `:1169`/`:1194` and passed by reference into `SPA_ws`; the lazily-constructed `ws()` deferral is defeated by taking a `KernelWorkspace&`), and none of the four records actually disappears.

### 4.8 Two more C++ dedups, smaller

- **`traceProbeBlock`** (2/2 verified, narrowed). `GetTrace`'s probe block (`src/SAIGE_fitGLMM_fast.cpp:3536-3569`) and `GetTrace_q`'s (`:3898-3931`) reduce to 30/31 stripped lines differing by exactly one (`tempVec0(i)=dot(uCol,Pu);`). Extract, with a signature that does **not** include `Sigma_iX`/`Sigma_iXt`/`cov1` (not used in the block), and delete the locals it orphans. **Keep `canBlockPCG`** — replacing it with a `pcgArm()` enum restructures the live scalar PCG loop around an arm LEDGER #114 records as uncallable, and `bda95327` recorded the dense-only gating as a verified design decision. `getDiagOfSigma` is **not** hand-inlined (`:3432` calls it). Note in the commit that the gate **cannot** prove this — step-1 goldens are toleranced; the real proof is `bda95327`'s protocol, capturing the PCG stdout log before and after on the dense golden config and requiring it identical line-for-line.
- **The `_LOCO` family** (1/2, contested: the clone-chain merge is refuted). Take only the **dead-code half**: `getCoefficients_q_LOCO` (`:3662`), `getSigma_X_LOCO` (`:3825`), `getSigma_G_LOCO` (`:3861`), `getAIScore_q_LOCO` (`:4025`), `fitglmmaiRPCG_q_LOCO` (`:4072`) have no caller anywhere except the RcppExports glue; two of them are not clones at all but stale copies of a superseded signature. Optionally collapse `Get_Coef_LOCO` into `Get_Coef(..., loco=TRUE)` in R (three call sites: `R/SAIGE_fitGLMM_fast.R:309`, `:576`, `:1321`). **Leave the three live C++ clones.** The merge would route a LOCO fit's preconditioner through `gen_spsolve_v4` — changing what LEDGER #114's owner decision retained, invisibly to every gate — and the invariant it claims to delete has fired exactly once (#91), at the leaf pair the remedy leaves duplicated. If ever done, four sites need `!t_isLOCO` guards (`:2896`, `:2899`, `:2923`, `:3053`), and the gate is the FULL tier with golden regeneration: four LOCO-derived step-2 goldens (`plink_{marker,group}_loco_chr{1,2}`) are the only thing that can falsify "behaviour-preserving".

### 4.9 R-side duplications

- **`ReadModel` / `ReadModel_subsample`** (2/2 verified) — `R/readInGLMM.R:37-177` and `:180-371` are 333 of the file's 491 lines, with **132 lines verbatim identical**, including the 8-line LOCO-missing `stop()`, the 6-line `--chrom` `stop()` and the entire 25-line `is_Firth_beta` covariate-offset block down to its `# LEDGER #68` comment. The cost is on the record twice: #110's PROVENANCE reads *"Both sites are byte-identical strings, so a fix must touch two lines"*, and #14 was a double-subsetting bug that existed **only in the copy**, undetected because no golden passes `--subSampleFile`. The copy has rotted further: `:217` nulls `obj.glm.null` (which `ReadModel` never does) and `:218` is dead code on a field NULLed at `:211`.
  **Merge on `subSampleFile`, not on a precomputed `includeIndex`** (the latter double-loads the `.rda`). Four subsample-only steps, not three. `keep()` is **not** an identity on the `res` line — `modglmm$res` is a *partial match* to `residuals` and `$<-` does not partial-match, so that line must stay behind the guard; and note that `:315` subsets the already-subset residuals a second time, LEDGER #14's exact shape, which needs its own row. Add value-pinning tests for the subset-only branches **first** — the only e2e coverage compares two runs that both traverse the changed code, so a symmetric error passes green.
  While there: add `readNullModelRda(file)` using `new.env()` (the idiom already at `R/PhenoScale.R:107`) with a guard naming the file when `modglmm` is absent. Today both functions `load()` blind into their own frame, so any variable in the `.rda` silently rebinds the formals `chrom`/`LOCO`/`is_Firth_beta`.
- **The covariate-offset GLM is fitted twice on every run** (2/2 verified). `R/SAIGE_fitGLMM_fast.R:1063` builds `mmat` and `:1065-1075` fits `modwitcov`; their only reader is `:1081`, inside `if (isCovariateOffset)`. On the FALSE branch `:1176` refits the identical GLM; on the TRUE branch `:1179` refits it as `fit0orig` and `:1180` rebuilds the same design matrix as `Xorig`. So **both** branches fit it twice, and `mmat` (N × p+1, never `rm`'d) stays live through the GRM fit on the branch that never reads it. Move the two lines inside the guard; note that the binary 0/1 phenotype check is at `:1163-1166`, so a malformed phenotype's error message moves from `glm`'s to the repo's own — state it rather than discover it. Optionally `Xorig <- mmat`, gated by an `identical()` drill against a scratch install, not by inspection.
  Extract **only** `setupSparseGRMForSamples` (the `getsubGRM`/`gen_sp_v2`/`setupSparseGRM` triple, ~13 lines × 2). **Do not** unify the covariate-offset algorithms — `fitNULLGLMM` uses `mmat[,-1] %*% coef[-1]` (intercept excluded) while `getNeff` uses `modwitcov$linear.predictors` (included) and, for quantitative, `lm` residuals. Two CLI flags carry **word-for-word identical help text** (`"Whether to estimate fixed effect coeffciets"`) with opposite defaults. That divergence is a LEDGER row plus a comment at each site, not a mode argument.
- **`getNeff`'s LOCO block** (2/2 verified) — `R/SAIGE_extractNeff.R:241-299` is a stale copy of `updateChrStartEndIndexVec` (`R/Util.R:30-89`) missing the LEDGER #62 contiguity guard, computing indices in the wrong space (§1.6), and writing C++ globals nothing on this path reads. **Do not delete the block or the `LOCO` formal** (exported signature, #114-class). Replace the inline copy with the shared helper: deletes ~60 of 65 lines, ends the drift, and fixes the `.bim`-hardcoded message for free. Collapse the two trait tails (`:402-438` / `:439-468`, 25 identical lines) **keeping the trait guard** — an unrecognised `traitType` currently falls through both arms and silently does nothing. Make the return `invisible(Nglmm)`.
- **`R/Geno.R`'s two resolvers** (2/2 verified for the "computes and discards" half; 1/2 contested for the wider claim). `checkGenoInput` (`:1-49`) derives and `Check_File_Exist`s `bimFile`/`famFile`/`pvarFile`/`psamFile` and then returns only `dosageFileType` at `:48`. `setGenoInput` re-derives all four at `:86-91`/`:131-136` and hands *those* strings to `setPLINKobjInCPP`. The branch already identified this hazard and fixed it for a different consumer — the comment at `:294-296` explains that `genoList$files` exists *"so a second sample-set group's reader can be built from exactly the same inputs without re-deriving"* — and that comment's own claim that nothing reads the field is now false (`R/MultiModel.R:165`).
  **Minimal, safe version:** have `checkGenoInput` return the resolved list as `checkGenoInputStep1` already does (`:439-443`), read it in `setGenoInput`, delete the eight `gsub` lines, flip `test-unit-geno.R:31` to `expect_identical(obj$genoType, "pgen")`. **Do not** touch `:297-304` — the `files` record carries `AlleleOrder`, which `MultiModel.R:172` reads and the resolved-path list does not have. **Do not** unify the two resolvers: their branch *predicates* differ (`bedFile != ""` vs `anyPlink`), so adopting step 1's would give step 2 a new *acceptance* — `--bimFile` alone currently stops loudly, which is LEDGER #15(c)'s fixed message pinned at `test-unit-geno.R:21-23`. The exclusivity divergence (`docs/flag-surprises.md` F10: `--pgenFile --bimFile=other.bim` silently discards the bim) is a **new hard refusal** and needs its own row.
- **The region `.index` writer** (2/2 verified) — replace `R/SAIGE_SPATest_Region.R:1933-1936` and `:2031-2067` with one `writeOutputFileIndex(...)` call, consuming the orphan `singleBytes` parameter (§1.6). Take `nEachChunk` and `AnalysisType` from one shared constant so `checkOutputFile` and the writer cannot state them differently — that is the LEDGER #25 class. Net ~26-30 lines. One stdout change to state: the helper's `cat("write to output\n")` moves after the flush and prints once per model. **The goldens do not verify this** (`*.index` is outside the byte gate) — the real net is `test-step2-region-restart.R`, `test-step2-region-restart-midchunk.R` and `test-unit-util.R`.

---

## §5 — The R↔C++ contract

### 5.1 `Main.hpp` is a third, unchecked copy of the step-2 signatures (2/2 verified)

303 lines, 28 declarations, **one includer** (`src/Main.cpp:24`; `LDmat.cpp`, the second includer #5's provenance names, was stripped). **Eight declarations disagree with their definitions today**, every one a function the multi-model or FlexRV work changed: `mainMarkerInCPP` (`:43`, a phantom `std::string & t_traitType`), `mainRegionInCPP` (`:103`, missing `t_nThreads` and `t_modelIdx`), `setSAIGEobjInCPP` (`:141`, missing `t_firthSolver`), `assign_conditionMarkers_factors_binary_region` (`:185`), `RegionSetUpConditional_binary_InCPP` (`:200`), `openOutfile` (`:204`), `openOutfile_singleinGroup` (`:207`), `copy_singleInGroup` (`:257`).

A differing parameter list is a new **overload**, not an error, so the drift is structurally uncatchable — and the only declarations that *can* drift are exactly the dead ones, since a mismatch on a real forward reference is a link error. LEDGER #5 was FIXED as "declarations aligned to definitions"; the row is accurate about what *it* fixed, but the header at large re-drifted through step-4 and FlexRV F3 and CLAUDE.md's standing trap is true again. Meanwhile `src/RcppExports.cpp` already carries the correct, compiler-checked declaration for all 24 exported ones — the header is a third copy, not a second.

**Delete it.** Only four declarations are load-bearing: `processMale_XnonPAR` (call `:1653`, def `:1711`), `writeOutfile_single` (`:1607`/`:5424`), `writeOutfile_BURDEN` (`:4461`/`:5547`), `writeOutfile_singleInGroup` (`:4570`/`:5645`). Move those into `Main.cpp` above first use. **`src/Main.hpp:7` is the only includer of `MarkerLoop.hpp` in the tree** — that include must be *added* to `Main.cpp`, not "folded". No `compileAttributes` is needed; `OBJECTS` is unaffected; headers change, so `rm src/*.o` and a clean build.

**Doc sweep, same commit** (otherwise four documents instruct future sessions to align a file that no longer exists): `CLAUDE.md`'s Architecture line, `docs/design/phase5-multiweights-design.md:73`, `docs/design/flexrv-plan.md:128`/`:217`/`:369`, `docs/state/next.md:1219`. **Append** to LEDGER #5's STATUS recording that the class is now structurally impossible; do **not** re-open it (that breaks `0 open of 142` and would require a characterization test for a bug that cannot be called). Consider a `docdrift.py` rule so *"Main.hpp declarations are stale"* cannot be copied back in.

Cheaper alternative the owner should be offered: re-align the eight in place. Same gate, no doc edits, leaves the header closer to upstream's 382-line version.

### 5.2 Config is pushed to C++ before R has finished deriving it (1/2, contested: the config-setter rewrite is refuted)

`setAssocTest_GlobalVarsInCPP` is called at `R/SAIGE_Test_main.R:345-356`, ~200 lines before `SPAGMMATtest` finishes computing what it pushes. Three later assignments therefore reach nothing: `max_MAC_use_ER <- 0` (`:508`), `minInfo <- 0` (`:551`), `minMAF <- 1/(2*length(sampleID))` (`:552`). The parameters actually pushed are `min_MAF`/`min_Info`; the locals have no reader.

**This is inherited verbatim** (`main:R/SAIGE_Test_main.R:270-272`), so it is bug-for-bug preserved code. What is *this branch's* is that six surfaces now assert the opposite:

1. `R/SAIGE_Test_main.R:544-549` — the `cat()` telling the user *"the flag replaces --minInfo and --minMAF"*, added yesterday in `c8571eea`;
2. `tests/testthat/test-step2-xnonpar-missing.R:83`, which asserts that exact string;
3. `harness/tools/docdrift.py:238`'s `139-xnonpar-minmaf-not-applied` rule, whose `why` field encodes the false mechanism — the blocklist that exists to stop false sentences recurring now contains one;
4. LEDGER #139's PROVENANCE;
5. `docs/flag-surprises.md:31`/`:149`;
6. `R/SAIGE_Test_main.R:466`, which correctly calls them *"dead locals"* — contradicting (1) 80 lines later.

**Do only the documentation defect**, as one commit: a LEDGER row recording that upstream's chrX `minInfo`/`minMAF` overrides have never been applied by either fork, plus corrections to (1)–(5). Note the observable divergence is narrow — `min_Info > 0` is a hard `stop()` (#107/#108) so only the `minMAF` half is user-visible, and only with `--is_rewrite_XnonPAR_forMales=TRUE` plus a user-supplied `--minMAF > 0`. `max_MAC_use_ER <- 0` is inert (every ER site conjoins `traitType == "binary"`).

**Reject the `setStep2ConfigInCPP` rewrite.** `setAssocTest_GlobalVarsInCPP` also constructs the output writers (`src/Main.cpp:281-289`, `g_out.regions`/`g_out.singles`) — a replacement assigning only the 19 `Step2Config` fields drops writer construction and every single-model golden dies; one that carries them and runs after `setSingleOutputPrefixesInCPP` clears the per-model writers instead. A test worth adding either way: run the chrX fixture with `--minMAF=0.02` and assert the marker count against the 1/(2N) prediction — it will match the 0.02 gate.

### 5.3 The reader layer's dead read mode (2/2 verified)

`GenoReader::getOneMarker` (`src/GenoReader.hpp:42-48`) advertises three read modes; production uses one. All live call sites pass literals (`isOutputIndexForMissing=true, isOnlyOutputNonZero=false`), and `t_isTrueGenotype` is hardcoded `true` at `GenoReaderRegistry.cpp:46`, making the guard block at `PLINK.cpp:450-455` unreachable.

The two readers implement the dead mode **differently**, and one of them pays for it per marker: `PGEN.cpp:571-573`'s `t_indexForNonZero.push_back(i)` is **not** guarded by the flag, so every PGEN marker does up to N `push_back`s into a vector the caller clears before `imputeGenoAndFlip` rebuilds it. `PLINK.cpp:505-508` carries a comment saying exactly this — *"the per-sample push_back here was pure waste"* — and the fix (`ee4a2b8d`) was applied to PLINK only. The would-be semantics also diverge: PLINK selects on the post-map value where `bufferG1` is `size_t` and missing maps to `-1` → `SIZE_MAX`, so its packed mode would emit every missing sample as non-zero; PGEN selects on the raw dosage.

**Do the one-line PGEN guard first** — it is the only part with a payoff, closes the divergence, touches no interface, and the PGEN mirror goldens are the control. The wider mode deletion (three parameters off the interface and both implementations, four call sites, two `j` counters, two `resize(j)`) is legitimate and provably unobservable, but note it replaces LEDGER #139's runtime `Rcpp::stop` in `applyMaleXnonPAR` with a structural guarantee — state that and update #139's and #69's WHERE lines. `isForkSafe` has no production consumer but is a **documented retention** (`src/Main.cpp:5287-5290`) covering a live hazard (the SKAT-O pool forks with readers alive); leave it.

### 5.4 The escalation signal travels as process state (1/2, contested: the argument conversion is refuted)

`g_regionEscalatedPass` (`src/Main.cpp:90`) is written by two exports and read once, at `:2256`. `set_flagSparseGRM_cur_SAIGE_org()` is called **only** on the path where the value is never read (`R/SAIGE_SPATest_Region.R:1876`, under `!is_fastTest`), so its write is dead by construction — and its only other effect is an unchecked `g_model->` dereference that, on a multi-model region run, reads model 0's flag on behalf of all of them. The comment at `src/Main.cpp:48-52` justifying the element-0 alias is stale: the region path has accepted P > 1 since 2026-09-06.

**Do not convert it to a parameter.** The claimed prize ("zero step-2 loop globals") is false — `g_pooledRecords`, `g_burdenAsked`/`g_burdenEvaluated`, `g_regionProfile`, `g_carrierStore` and two function-scope statics remain — and the 11-line comment at `:78-89` that the remedy counts as deleted complexity *is* the written-down protocol. **Free version:** append two sentences to that comment recording the unreachability. **If the dead setter must go:** replace `R/SAIGE_SPATest_Region.R:1875-1879` with an unconditional `set_flagSparseGRM_cur_SAIGE(FALSE)` and delete the `_org` export, updating `src/Main.cpp:50`, `src/SAIGE_test.hpp:47` and `docs/state/multimodel-step2.md:72`.

### 5.5 Output-writer allocation hidden in a config setter (1/2, contested: the hazard is loud, not silent)

`setAssocTest_GlobalVarsInCPP` allocates writers at `src/Main.cpp:282-288`; `setSingleOutputPrefixesInCPP` (`:5274-5285`) then clears and rebuilds **both** lists despite its name, and is called only `if (nModels > 1)`. On every multi-model run the region writer built by the first call is constructed and immediately discarded.

The order dependency is documented in three places and enforced by two explicit `Rcpp::stop`s (`:1330`, and the bounds-checked `regionWriter()` at `:158-163`), so drop the "silent" framing. What survives: the setter is named for half of what it does, and an unnecessary conditional. Call `setOutputPrefixesInCPP(OutputFiles)` unconditionally (for one model `OutputFiles[1] == OutputFile`), drop `t_outputFilePrefix` and the five writer lines from the config setter, and rename. Required in the same commit: update `src/Main.hpp:10-20`, and rewrite — not delete — the two comments at `src/Main.cpp:5268-5271` and `R/SAIGE_Test_main.R:358-360` that record the decision being reversed, since a single-model run will now execute a `.Call` it previously skipped. Also fix the stale comment at `R/SAIGE_Test_main.R:333-334`. The gate is **not** pre-paid — if §5.2 does not land, this is its own `src/`+`R/` change at ~7 min. Rename in a separate commit (it is an API change under `exportPattern`).

### 5.6 The weight-set layout is composed three times, in two languages (2/2 verified — see §1.5 for the actionable part)

`R/SAIGE_SPATest_Region.R:1224-1230` composes `composedNames`; `:1316-1326` composes `weight_sets_kind` in the same order, 92 lines away; `src/Main.cpp:2356-2392` composes it a third time from argument *shapes*. `weightSetDescriptor` zips the two R vectors positionally, and `stopifnot(length(...) == nWeightSetsUsed)` at `:1335` checks length, never alignment.

Verifier corrections cut this down: the R-side lengths are anchored to shared scalars so length drift is structurally impossible (the residual risk is reordering one expression and not the other), and the alignment **is** reconciled by the repo's two-routes pattern — `test-step2-multiweight-cond.R:97-119` over `helper-multiweight.R:110-127` runs the K=2 CADD rows against a K=1 CADD-only run, `test-step2-flexrv.R:133-139` does the same for FlexRV, and the composed name *order* is pinned at `test-step2-multiweight.R:51`/`:173` and `test-invariants-goldens.R:258`.

**So this is a locality finding, not a correctness hazard**, and the extraction the finding proposed relocates: the legacy no-weight set's weight *vector* is built at `:1780`, outside any plausible extraction; ~25 downstream consumers read the twelve derived locals; and the two fall-throughs go 2 → 1, not 0.

**Cheap in-place version, ~5 lines:** build `composedNames` and `weight_sets_kind` from ONE `c()` of `(kind, name)` pairs at `:1224` and read the two columns off it at `:1296` and `:1316`. Or strengthen `:1335` to assert the partition rather than the length. Leave the FlexRV block and #138's refusal alone.

### 5.7 CLI option surface: no single owner (1/2, contested: the manifest rewrite is refuted)

179 `make_option()` blocks across four wrappers restate what the exported R functions already declare; `SPAGMMATtest` repeats the text a third time in 65 roxygen `@param` lines, and `man/SPAGMMATtest.Rd` carries it a fourth (generated, per LEDGER #93 — which makes the generated file *machine-produced evidence* of the drift). Four live divergences:

- **`--X_PARregion` has three values**: CLI default `""` (`extdata/step2_SPAtests.R:350`), its own help + roxygen + `man:161` say `…155260560`, the formal (`R/SAIGE_Test_main.R:186`, `man:45`) says `…155270560`. `strsplit("", ",")` is `character(0)`, so the CLI default skips the setter entirely and `g_cfg.X_PARregion_mat` stays 0-row. LEDGER #136 edited that exact 12-line help block last week without noticing.
- **`maxMAF_in_groupTest`**: CLI `"0.0001,0.001,0.01"` (three cells) vs formal `c(0.01)` (one). A direct R caller gets a third of the group-test grid a CLI user gets. This is the #18/#82 shape, unrecorded.
- **`numRandomMarkerforSparseKin`**: 2000 (CLI) vs 1000 (formal) — inherited, and only possible because nothing reads it.
- **`--isCovariateOffset`**: TRUE (CLI) vs FALSE (formal), both help strings saying *"Whether to estimate fixed effect coeffciets"* while `:1082` prints the opposite. `docs/cli-options-and-group-file.md:144` is a fifth surface contradicting the CLI help for the same flag.

**Reject the option-spec table.** It breaks `test-step2-loco.R:153`, which regexes the literal `make_option` source for `default = TRUE` (LEDGER #111's pin); 905 of the 1,621 in-block lines are help prose that must be retyped, not deleted; 163 of 179 blocks are upstream; `man/` is already derived.

**Do this instead — `harness/tools/optiondrift.py`**, beside `ledger.py lint` and `docdrift.py`, in the same seconds-long pre-tier block. Paren-match each `make_option(` block in all four `extdata/*.R`; parse `names(formals())`/defaults and each `@param`; assert per option, **comparing values not strings** (so `1L`/`1`, `1e-05`/`1e-5`, and comma-strings vs `c(...)` after `strsplit` all reconcile): the `[default=X]` in the help equals the declared `default =`, and the CLI default equals the formal's or the pair is on an explicit allowed-divergence list carrying its ledger row. ~150 lines, no `src/`, no `R/`, no upstream churn, no gate. LEDGER #115's STATUS already records why a hand sweep is not enough (*"a sweep that is worth trusting must parse the blocks, not scan a fixed number of lines, and compare VALUES"* — its own re-sweep produced 19 false hits from a regex window running past a block's end).

File the two live divergences as rows first. Deciding whether CLI or API wins is an owner decision of the #18/#82 shape.

---

## §6 — The god functions

Last, because they are the biggest grind and the least likely to be worth a single push. Both must be sliced.

### 6.1 `SAIGE.Region` — 1,251 lines, seven closures, nine `<<-` (2/2 verified)

`R/SAIGE_SPATest_Region.R:883-2133`, behind a 43-line, 29-parameter signature. Upstream's is 1,209 — **this branch extracted the 675-line `region_postprocess` and then put ~700 lines of new work back into the same frame**. Named concerns in one body: multi-model fan-out (`:926-949`, plus nine `for (k in seq_len(nModels))` loops), restart/truncate bookkeeping (`:1008-1030`), r.corr → test-type resolution, SKAT-O pool startup, three file-opening sites 200-300 lines apart, FlexRV setup (`:1095-1167`), weight-set composition (`:1169-1335`), the step-5 batching pipeline (`:1466-1736`), the scan loop, a 145-line chunk writer (`:1932-2076`) and message assembly.

Seven closures (`:1413`, `:1429`, `:1487`, `:1510`, `:1527`, `:1571`, `:1693`) communicate with the loop through nine `<<-` over seven enclosing names (`:1660`, `:1663`, `:1669`, `:1673`, `:1674`, `:1683`, `:1715`, `:1734`, `:1735`), and `skatoBatchBytes` is additionally written with plain `<-` at `:1909` — the same binding mutated two ways.

**Corrections that shrink the claim:** the per-region pipeline already threads an explicit `S` object, so the `<<-` web is only the cross-region accumulator layer; the stage contract *is* written, at `:1466-1480`, with each `<<-` carrying its ledger citation; and `region_postprocess`, `skatoCollect`/`skatoEvaluate`/`skatoServe` and five other helpers are already extracted and unit-tested. What is left in the body is orchestration and IO — which is precisely what the 1,251 lines now are.

**Slice it, in this order, one gated commit each, after the FlexRV work settles:**

1. `regionStateFor(region, …)` — the 100 inline lines at `:1771-1870` (WEIGHT canonicalisation, the negative-weight guard, SCORE/ISLOF, the degenerate-score warnings). Clearest win, self-contained.
2. `regionTestPlan(r.corr)` (`:1032-1056`) — **must return `r.corr`**, not just the method: `SKAT_Check_Method` promotes `0` to the 7-element rho grid `c(0,0.01,0.04,0.09,0.25,0.5,1)`, consumed at `:207`, `:539` and `:1493`. Drop it and every SKAT-O golden moves or the run dies.
3. `regionDoneMessage` (`:2083-2131`), and delete the dead `.markerInfo` message at `:2087-2093` (the string occurs nowhere else; `:2108` overwrites the variable unconditionally).
4. `composeWeightSets` scoped to `:1181-1253` **only** — pure composition. Leave the restart-schema guard (`:1255-1272`), the `openOutfile` loop (`:1274-1287`), `nline_per_gene` (`:1289`) and `region_weight_names` (`:1296`) in the body: they are I/O and C++ global-state mutation, not composition. `regionFlexRVSetup` emits warnings and stdout whose exact text six test files assert — it is not pure either, and must keep `weightSetNames` for the `:1151` collision check.
5. Only then the accumulator: `scan <- new.env()` carrying **eleven** mutable fields (the seven superassigned ones plus `cth_chunk_to_output`, `iswriteRegion`, `iswriteMarkerList`, `OutList.all`), with a read-only `cfg`. Use an environment, not a returned list — a list-returning stage forces every call site to rebind, which is where a dropped field becomes a silent stale value of exactly the #81 shape.

**Do not unify the chunk writer.** The five `fwrite` sites carry **three** distinct `(append, col.names)` rules — the results file has two branches, the marker list has three including a mid-run **truncate** — and LEDGER #59's STATUS records the asymmetry as deliberate (*"The truncate-vs-append decision is UNTOUCHED — only the header moves"*). No golden and no test combines `--is_output_markerList_in_groupTest=TRUE` with a multi-chunk run, so a wrong unification would pass green. If the nesting is worth attacking, do it per *writer*, not per rule.

**Do not touch** `regionRestartPlan`/`writeRegionIndex` in the same work: `*.index` is excluded from the byte gate, so the two restart tests are the entire coverage. **Do not reorder** the two constraints inside the file-opening sites: the SKAT-O pool is forked at `:1058-1060` *before* the first `mainRegionInCPP` because forking a TBB-owning process is the hazard, and the BURDEN `openOutfile` sits at `:1282` because its header needs `nWeightSetsUsed`.

### 6.2 `mainRegionInCPP` — 2,501 lines, ten local types, five lambdas (1/2, contested: the phase lift is refuted)

`src/Main.cpp:2207-4707`, **42.7% of a 5,850-line file**. Upstream's is 1,148, so this branch added ~1,350 lines *inside* one body. Ten types are declared function-locally (`URCellModel:2518`, `BurdenCellModel:2528`, `RegionModelWork:2537`, `RegionPrepModel:2699`, `RegionPrep:2712`, `RegionBlockWorker:3120`, `CellOutcome:3167`, `RegionCellWorker:3172`, `URCell:3372`, `BurdenCell:4131`), and `urPrepFn` (`:3409-3903`) is 495 lines in one `[&]` lambda.

**Take only the mechanical type hoist:** create `src/RegionLoop.hpp` (namespace `Region`) and move those ten declarations verbatim. No `Region::Work`, no `Region::SetLayout`, no phase lift, no new translation unit.

Why the rest is refused, in the verifier's terms: the body has **89 distinct function-scope locals**, so a `Work` struct is an ~85-member relocation of the same coupling with *wider* reachability (a header type is visible to every includer where a `[&]` capture is bounded by one function); `SetLayout` is not pure (it reads the file-static `g_cfg` and makes four `Rcpp::stop` calls) and is already brace-scoped and already covered end-to-end; the marker-path analogy fails on the property that made it work (a marker is a pure function of reader/index/model, while the region body's `MAC_GroupVec`, `genoSumcount_noweight`, `i1`, `numofUR`, `mPassCVVec` are cross-phase accumulators mutated in order); and a new TU changes inlining and FP-contraction opportunities in a body whose own comment at `:3451` records that armadillo's `accu` *"unrolls it into several partial accumulators … exactly the kind of 1-ULP change that moves a printed golden"*, next to a hand-rolled `paritySumDense` written to match armadillo's parity-split sum bit for bit.

Two further unpriced costs the hoist avoids and a phase lift would not: `src/Main.cpp` has 78 commits in 30 days and is the file CLAUDE.md's allele-convention section, LEDGER #51/#118/#130/#131 and `docdrift.py`'s evidence strings all cite **by line**; and moving the two `Worker` structs is a change to region threading, which mandates `harness/tsan.sh` region clean **and** plant (hours, not seven minutes).

The one live hazard in this frame is already handled by §3.3 and §4.2: the frame-scope `hasVarRatio` at `:2671` that `e72e4186` left standing, shadowed by hand at three threaded sites.

---

## §7 — Test suite and harness

Cheap, and two of the three items buy detection power.

### 7.1 Child-process spawners lose their diagnostics (2/2 verified, narrowed hard)

Fourteen sites across thirteen test files hand-roll tempfile + `writeLines` + `system2(Rscript)` + status check. The **actionable** half is diagnosis, not duplication:

- `test-step2-region-escalation.R:60` and `test-step2-skat-failures.R:59` call `expect_cli_ok` on a bare `system2` status with no `log` attribute — a child failure prints `"<no log>"` while the log sits on disk.
- `test-step1-guard.R:137`, `test-step1-genostore.R:41`, `test-step1-pgen.R:97` and `:218` pass `stdout = FALSE, stderr = FALSE`, discarding the child's output entirely.
- `test-step2-api-defaults.R:30-33` and `test-unit-step1-reingest.R:48-53` hand-copy `expect_cli_ok`'s body.

**Add `spawn_r_script(script, args, env, log)` to `helper-saige.R`** returning a status with `attr(status, "log") <- log`, deliberately the same shape as `run_saige_cli` so `expect_cli_ok` works unchanged. Fix the seven sites. ~20 lines added, ~15 deleted, 7 files.

**Keep the string bodies.** The `deparse(quote(...))` remedy is refuted: the 243 payload lines are all *different* code, so a shared helper cannot make them lintable; `reingest_probe` (`test-unit-step1-reingest.R:13-58`) takes the child body as a parameter and is called six times with variable statement sequences — it is a code generator; and `test-step2-region-tworun.R:111-141` holds a cross-run byte comparison, so touching it escalates the change to both compare modes.

### 7.2 Forty-five identity aliases for a lint no gate runs (2/2 verified, but note the contested twin)

Forty-five identity bindings across nineteen files (`fx_cli <- run_saige_cli`, `bam_cli`, `pc_cli`, `fl_cli`, `nid_cli`, …) plus eight `get("saige_goldens", envir = parent.frame())` lookups of two constants that are lexically visible at call time. Seventeen private names for one function. No gate runs lintr: `harness/run_tests.sh:45` and `harness/prepush.sh:24` run only `ledger.py lint`, and CI's only lint-ish line installs `lintools`, a numeric package.

A parallel finding attacking the same structure was **refuted** on two of its claims, and both refutations apply here too: the aliases are *not* ungreppable (the alias line itself contains the literal, so `grep -rl run_saige_cli tests/testthat` already returns all 70 users), and `.lintr` must not be touched — `globals:` is not a lintr key, and `exclusions: tests/` or `object_usage_linter = NULL` at the root drops every linter over 127 files / 18,642 lines, or drops it for `R/` too.

**Use the mechanism already in the tree:** `# nolint start: object_usage_linter.` / `# nolint end`, as `test-step1-guard.R:30` and `test-step1-ai-projection.R:26` already do — and note that ~49 test files call these helpers inside function bodies with no alias and no suppression at all, so most aliases can simply be deleted in place. Fix the self-contradicting prose while there: `helper-step2-flexrv.R:42-44` says `helper-saige.R` *"comes after this one"*, which is false and is contradicted by `:56-58` in the same file. Nineteen test files at `standard`; no cross-run comparison added.

### 7.3 Twelve memoized CLI runners (1/2, contested: the shared runner is refuted)

Fourteen per-file `new.env(parent = emptyenv())` caches, three keyed on a hand-written **tag** rather than on the arguments (`loco_run(key, extra)`, `miss_run(key, extra, geno)` where `geno` is not in the key at all, `fl_run(extra, tag)`), and two carrying a `key <- if (length(extra)) … else "<none>"` line because `env[[""]]` errors — a footgun rediscovered twice.

**Do not build a shared runner with an args-derived key.** It is falsified by the suite as it stands: `test-step2-threads-goldens.R:34-37` runs three arms with **byte-identical args** distinguished only by the tag and then asserts `files_byte_identical(outs[1], outs[2])` — under an args-keyed cache that compares a file with itself; the restart tests have the same shape; three runners vary only by `env`; `test-step2-loco.R:284-286` reads its cache directly as data; and `urzerod_run` writes a fresh group file per call. Per-file caches are also *deliberate* — testthat runs each file in its own process, so a shared env dedupes nothing (`helper-multiweight.R:1-9` and CLAUDE.md both record the trade).

**Take only:** one `write_model_manifest(entries)` in `helper-saige.R` with `rm_manifest`'s signature (it subsumes the other three encoders of the `name rda vr` production format), one `pnum`, and de-duplicate the `rel_err_vec` shadow at `test-step2-pgen.R:325` against `:67`. **Leave `rel_err_vec` itself alone** — swapping it for `helper-saige.R`'s `rel_err` adds an absolute floor that would suppress exactly the rs20/rs22 `Tstat`/`var` cells (3.09e-15, 3.47e-29) the pgen-vs-plink reconciliation compares at full relative tolerance, and its stated justification is stale (post-#44 those goldens print `NA`, not `nan`).

---

## §8 — Leave alone

Explicitly, so nobody re-raises them:

- **`src/CCT.cpp` and `R/CCT_modified.R`.** `src/CCT.cpp` is byte-identical to upstream; the branch order is SAIGE-GENE+'s published prescription (LEDGER #86, quoted on the row); the two implementations exist because the BURDEN pooled row is written inside `mainRegionInCPP` with no R callback available, and `test-unit-cct.R` reconciles them on every branch. The four rule-namers (`CCT`, `CCT_cpp`, `pooledRule`, `pooledCCTRule`) are 2 roles × 2 languages, and `"none"` is not a CCT branch at all.
- **The `dbeta` split.** C++ computes a beta set's weight with boost; R recomputes it with `dbeta`, and both enter one cell's statistic. This is **normative** (design doc D2's fourth column is literally "R-side AnnoWeights"), and D1 makes the K=1 byte-identity structural. Unifying re-baselines ~122 group goldens and destroys the check the byte gate currently provides. Current disagreement is ~1e-15. If anything is wanted here it is a tighter reconciliation test — and one already exists at `test-step2-region-flip.R:241-260`.
- **`region_postprocess`'s two-world K=1/K>1 arms.** `docs/design/phase5-multiweights-design.md:110` records *"branch on K, not a rewrite of the K=1 arm"* as the mitigation for byte drift. The `sum(WEIGHT) > 0` discriminator is not interchangeable with `kind` dispatch (a zero-summing K=1 weight column currently stops the run; under `kind` dispatch it would silently emit nothing — LEDGER #119's residual shape).
- **`.singleAssoc.txt_temp`.** (1/2, contested — and the contest matters.) It is a filesystem round-trip of a vector that never leaves the process, never deleted, undocumented, and eight copies are byte-gated. But deleting it removes the **only** byte gate on the fastTest *first-pass* single-variant writer: `plink_group_fastTest_escalated`'s `.singleAssoc.txt` comes from the accurate pass, so a regression confined to the fast-pass writer would move the `_temp` golden and nothing else — the LEDGER #50 surface, which lost 19 rows at exit 0. The remedy also does not compile as written (`writeOutfile_singleInGroup` and `writeNumOrNA` take `std::ofstream&`, and `writeNumOrNA` is shared with the marker writer and the BURDEN emitter). **Document it; do not delete it.** The only free piece is `copy_singleInGroup`, which is dead as a call and whose `Main.hpp:257` declaration takes no arguments against a definition taking `int` — that rides §5.1.
- **`NAMESPACE`'s `exportPattern`.** It is not a blocker (deleting a file removes its exports with no NAMESPACE edit) and narrowing it needs ~35+ names, not 6, because of 19 `SAIGE::` test references and five `library(SAIGE)` child snippets calling ~22 bare names.
- **`checkArgsBool`/`checkArgsNumeric`** (`R/Util.R:1`, `:7`) — `test-unit-util.R:109` is LEDGER #15(a)'s inverted characterization test, and `test-unit-checkargs.R:6-9` is a written warning against consolidating them.
- **`checkArgs.R`'s wrapper structure, `getMarkerPval`'s body, `SPA_binary`'s `Isconverge`, `src/Makevars.win` and the commented OpenMP block** — all upstream verbatim, all inherited defects rather than local ones. `SPA_binary`'s forced `Isconverge = true` deserves a ledger row, not a fix; `SPA_pval` → a three-line wrapper over `SPA_ws` is the one genuinely free collapse in that family (same calls, same order, same expressions), and the test repoint the finding proposed is mechanically impossible (`SPA_pval` returns a `double`; `isSPAConverge` is copied by value through Rcpp) and would delete LEDGER #22's inverted test.
- **`createSparseGRM`'s `isSetGeno`/`isWritetoFiles`** (1/2, contested). The unreachable `getSparseSigma` branch returns a sample-ID vector describing the full fam beside a matrix sized by the phenotyped subset — a real latent defect, worth a **row** cross-referenced to #114, not a signature change. The stated value ("unblocks threading a `GenoReader` through step 0") is not real: the obstacle is the global `geno`, which the remedy leaves in place.

---

## Verdict

**Trending up, clearly — but the instruments are thinner than the discipline.**

The evidence for "up" is specific. This branch invented `TestContext`, `MarkerResult`/`MarkerScratch`, `CarrierStore`, `RegionChunkStore`, `preparePhenotypeFrame`, `weightSetDescriptor`, `indexChunkLine`/`parseIndexChunkLine`, `ledger.py lint`, `docdrift.py`, `tsan.sh`'s plant arm and `guardmalloc.sh`'s — every one of them a mechanism that turns a claim into a check. It closed 142 ledger rows with zero open. It runs a byte-exact gate on a 7-minute budget and re-derives the gate's own cost when it drifts. Most of the dead surface in this report was *already inventoried by the repo itself* in 2026-09-08 and consciously parked; five of the findings that were refuted were refuted precisely because the repo had already decided them. That is a codebase with a memory.

The evidence for the caveat is equally specific, and it is the same fact five times: **the only oracle is a byte comparison over one fixture, and every defect class in this report is invisible to it.** `-w` silences the one free second oracle. `tsan.sh` aborts the arm that proves the first oracle works. `prepush.sh` writes a false sentence into the artefact the branch quotes as evidence. `writeOutputFileIndex`'s `singleBytes` parameter has had no caller since the day it was added, because `*.index` is outside the byte gate. And the single most expensive structural fact — one rule written N times, a fix landing in one copy — has produced #44→#66, #118→#134→#135, #91, #37, #59's untaken half, and one live unrepaired twin (§1.1), because a fix to one copy leaves the other's golden green *by construction*.

**Three things most need doing.**

1. **Repair the controls, this week.** `tsan.sh`'s early abort (2 lines, no gate) and `prepush.sh`'s log line (~12 lines, no gate) first; then drop `-w` and ratchet the warning count locally with `-Werror=uninitialized`-class flags in CI. That last one hands the branch a diagnostic that finds uninitialised reads on live paths (`src/SAIGE_test.cpp:377`, `:1151`) which no fixture will ever reach — and it costs one gate.
2. **Write the nine ledger rows.** They cost a `ledger.py lint` run and they are the audit's real output: the unrepaired #104 twin, the markerList restart truncation, the stale sidecars, the BURDEN header/rows divergence, the unchecked `--minGroupMAC_in_BurdenTest`, `getNeff`'s broken quantitative path and its wrong-index-space LOCO block, the CLI/API default divergences, the `createSparseGRM` mismatched return. None of them may be fixed in passing.
3. **Spend one push on the dead-surface sweep and one on the deduplications.** The sweep (§3) removes ~3,000 lines including an uncompiled translation unit, a decl/def mismatch, a guard that reads as armed and is not, and three unread parameters on a live writer — one gate, one clean build, no golden movement. The deduplications (§4) are where the recurring defect actually stops: `contextForMAC` + the `hasVarRatio` demotion, the single-variant row writer, the p-value formatters, the BURDEN emitter's six loops.

**What to leave alone** is in §8, and the shape of it is worth stating plainly: *do not restructure inherited code for readability alone.* `src/CCT.cpp`, `checkArgs.R`, `getMarkerPval`'s body, the `dbeta` split, the K=1 legacy arms and `.singleAssoc.txt_temp` are all either byte-identical to upstream, protected by a written design decision, or the sole gate on something. This branch's ability to answer *"is this upstream's behaviour?"* with `git show main:` is a live, load-bearing tool — CLAUDE.md's allele-convention section and LEDGER #51's decision both rest on it — and every gratuitous rewrite of an inherited body spends it.

The one thing that would change the trajectory more than any item above: **`mainRegionInCPP` and `SAIGE.Region` are both larger than upstream's after this branch extracted 675 lines from one of them.** Every new feature — multi-weight, FlexRV, multi-model, carrier space — has landed *inside* those two frames rather than beside them. Until that reverses, the duplication tax in §4 regenerates faster than it can be paid down.
---

## The execution queue

Ordered by value per unit of risk. Phase 0 (writing this handoff) is done.
The gate arithmetic each phase assumes is in the report header.

## Phase A — repair the controls (no gate, do first after Phase 0)

The highest value-per-risk work in the audit. Each buys back detection power the
byte gate structurally cannot provide.

**A1. `harness/tsan.sh` early abort — 2 lines.** (§2.1, 2/2 verified)
`:58` sets `set -euo pipefail`; the died-arm branch ends `[ "$MODE" = clean ] &&
DIED=1` + bare `return`, so under `plant`/`plantreader` the `&&` list has status
1, `return` propagates it, and errexit kills the script at `:296`. The operator
sees one arm and exit 1 — reading as "the harness broke", not "the plant worked
and t=4/t=8 never ran". **The plant is the positive control for every `clean`
result on this branch**, including the `races=0` baseline the region workload's
fork-edge annotation rests on.
Fix: `if [ "$MODE" = clean ]; then DIED=1; fi` + `return 0`; make the other early
return at `:284` explicit. Reject the proposed `harness/lib.sh` (does not exist;
`guardmalloc.sh`'s `verdict()` is a different shape).
Drill both ways, cheap: `SAIGE_TSAN_WORKLOADS=binary SAIGE_TSAN_TIMEOUT=5
SAIGE_TSAN_MAX_MARKERS=200 bash harness/tsan.sh plant` — unpatched: one arm,
exit 1; patched: three arms, guidance block, exit 0.

**A2. `harness/prepush.sh` logs what it ran — ~12 lines.** (§2.3, 1/2 contested)
`:102` logs `concurrent=$CONCURRENT` on both paths; with the default
`MODES=both`, `set -- $MODES` leaves `$#=1`, the single-tier branch runs, and the
log records `concurrent=1` for a serial run. Same shape as the defect the file's
own header warns about.
**The knob is NOT dead** — `MODES="byte numeric" CONCURRENT=1` reaches `:165`,
and its retention is recorded in three places. Do not delete it. Pass the fact
into `run_tier` so it logs what it *did*; warn when `CONCURRENT=1` arrives with
one mode; extract the byte-identical 8-line dirty-tree block at `:153-160` /
`:192-199`.

## Phase B — write the nine ledger rows (`ledger.py lint` only)

The audit's real output. None may be fixed in passing. (§1, §5.2, §3.6)

1. §1.1 — the unrepaired pre-#104 burden-variance twin
   (`R/SAIGE_SPATest_Region_Func.R:433-448`), live on every binary group test
   with `--condition`. **Instrument before writing the row**: dump `var1,
   sum(Phi), VarQ, VarQ_2, r` for rs30/rs79 on golden P3 — if `r == 1` today the
   row is a coverage statement, not a live divergence.
2. §1.2 — `.markerList.txt` loses every pre-restart row on a region restart
   (`R/SAIGE_SPATest_Region.R:1371`). Candidate fix is one line
   (`rep(isappend, nModels)`), **not** the index-size route.
3. §1.3 — stale `.pooledTests.txt` / `.skatFailures.txt` / `.markerList.txt`
   survive `--is_overwrite_output=TRUE` (`R/Util.R:20`).
4. §1.5 — the BURDEN header and its rows decide the `Weight` column from two
   derivations; `hasCustomWeights` at `src/Main.cpp:2317` is the one live
   divergence surface (all-zero WEIGHT matrix → header with a `Weight` column,
   rows without it, exit 0).
5. §1.6 — `--minGroupMAC_in_BurdenTest` accepted unchecked (`R/checkArgs.R:57`);
   LEDGER #126 made it load-bearing. NA threshold → `if (NA)` death with no
   argument name — the #21(d) shape.
6. §1.6 — `getNeff(traitType="quantitative", noEstFixedEff=TRUE)` errors
   (`covoffset` assigned only in the binary arm,
   `R/SAIGE_extractNeff.R:353` vs `:446`).
7. §1.6 — `getNeff`'s LOCO window block computes indices in the wrong space
   (`:244`, raw vs QC'd-marker space).
8. §5.2 — upstream's chrX `minInfo`/`minMAF` overrides have **never** been
   applied by either fork, and **six surfaces now assert the opposite** —
   including `docdrift.py:238`'s own rule, whose `why` field encodes the false
   mechanism. Correct all six in this commit; the code stays.
9. §5.7 — the CLI/API default divergences: `--X_PARregion` (three values),
   `maxMAF_in_groupTest` (3 cells vs 1), `--isCovariateOffset` (opposite
   defaults, identical help text). Deciding whether CLI or API wins is an
   owner decision of the #18/#82 shape — record, do not resolve.

Also §8's parked rows if wanted: `SPA_binary`'s forced `Isconverge`, and
`createSparseGRM`'s `getSparseSigma` returning a full-fam ID vector beside a
subset-sized matrix.

## Phase C — the dead-surface sweep (one push, one gate, clean build)

~3,000 lines. §3. Frame the commits as **executing the repo's own 2026-09-08
audit**, not as discovery — and amend that audit's table in
`docs/state/next.md:1333-1375`, whose "exactly one occurrence in `src/`" rule
could never match an `[[Rcpp::export]]` and whose conclusion is wrong at HEAD.

Six commits, one push. **Predict "no golden moves" in writing first.**

- **C1** §3.1 `Binary_*`: delete `SKATO_Exact`, `SKAT_Permu`, `GetProb`,
  `SL_GetPermu`, `Test1`, `SL_runif_int`, `SL_out`, `mat_to_std_vec` and the
  commented blocks. **Delete the `Binary_ComputeExactMC.{cpp,h}` orphan first** —
  git-tracked, `#include`d, absent from `OBJECTS`, and including a header that
  does not exist. Do **not** merge `SL_Sample`/`SL_GetSample` (identical bodies,
  *swapped* parameter order, golden P33 gates the live arm). Rewrite the S2
  comment at `Binary_ComputeExact.cpp:383-386` (LEDGER #2 cites the deleted
  names) and re-run `docdrift.py`.
- **C2** §3.2 delete `Unified_getMarkerPval` + its `Main.hpp:70` declaration +
  the three steer-away comments. **Owner's call** on the three `[[Rcpp::export]]`
  ones (public API under `exportPattern`, the class #114 chose to retain). **Do
  not** bundle the `t_genoType` parameter removal.
- **C3** §3.3 the 35 dead region declarations + the two indeterminate-`bool`
  reads (§1.4) + `numofUR0`/`mFirth` + the three dead params on
  `writeOutfile_singleInGroup` + the frame-scope `hasVarRatio` at `:2671` that
  `e72e4186` left standing. One comment per deletion site in the style of
  `:3823` so the next upstream merge does not re-import it.
- **C4** §3.4 delete `prepSPAContext`/`spaContextReady` (no call site since
  `b60141f9`; it *reads* as armed, with a ledger citation). Stop materialising
  `gNA`/`muNA`; drop them only from the three **internal** functions. Keep
  `:721-726`. Fix `gNB.clear()` ×2 → `muNB.clear()`. **Append** to LEDGER #40's
  STATUS naming `b60141f9`; do not rewrite it.
- **C5** §3.5 delete `src/getMem.*` and **all five** `#include "getMem.hpp"`
  lines; delete `UTIL::imputeGeno` (declaration and definition have different
  arities) and **repair the line anchors it shifts** in LEDGER #28/#109,
  `test-step2-missing.R`, `src/PLINK.cpp:508`, `docs/flag-surprises.md`. Leave
  the rest of UTIL alone, especially `imputeGenoAndFlip`.
- **C6** §3.6 R-side: delete `R/Firth.R` (450 lines of RareEffect residue; the
  NAMESPACE blocker does not exist), `ReadGroupFile`, `check_close`,
  `Get_Variance_Ratio`'s 19-line post-`stop()` tail and two unread formals,
  `region_postprocess`'s `resultDF`, the `printComb(3)` **call line only**.
  **Do not** narrow `NAMESPACE` (~35+ real exports, not 6). **Do not** touch
  `checkArgsBool`/`checkArgsNumeric` (#15(a)'s inverted test warns against it).
  **Do not** delete `getNeff`'s unread formals (exported; #93's measured
  invariant) — that is a row, not a deletion.

Verify `OBJECTS` explicitly after C1: R's `shlib.mk` guards `$(OBJECTS)` behind
`@if test "z$(OBJECTS)" != "z"`, so an empty list **silently skips the link** and
leaves a stale `.so` that passes the gate while proving nothing. Check the link
line in the build log and `ls src/*.o | wc -l` = 20.

## Phase D — turn the compiler back on (one gate)

§2.2, 2/2 verified. `src/Makevars:1` carries `-w` beside a `-pedantic` it
silences, and `-w` *cancels* the `-Wall` R passes by default. Nothing reads a
compiler diagnostic anywhere.

Measured with the pinned clang, `-w` removed: **306 diagnostics, 270
src-anchored, 139 in hand-written code** — including 8 uninitialised-use sites,
two of them on the **live SPA / conditional path** (`src/SAIGE_test.cpp:377`,
`:1151`) that no fixture reaches. Codegen identity was verified, not asserted:
`SPA.o`, `SAIGE_test.o`, `CCT.o`, `UTIL.o` are `cmp`-identical at `-O2` either
way.

Local: drop `-w`, keep `-pedantic`, add `-Wall -Wextra -Wno-unused-parameter
-Wno-cast-function-type-mismatch`; have `harness/build.sh` count `warning:` lines
matching `^src/` **excluding generated `RcppExports.cpp`** against a committed
ratchet baseline. **Do not put a count ratchet in CI** — different OS, different
compiler, `pixi.lock` uncommitted, every dep `"*"`. In CI promote only
bug-shaped categories: `-Werror=uninitialized -Werror=sometimes-uninitialized
-Werror=return-type -Werror=parentheses`. Do **not** let the first green push
"fix" `SAIGE_test.cpp:1000` — `m_is_Firth_beta & pval <= m_pCutoffforFirth` is
correct. The two live uninitialised reads are **LEDGER material, not drive-by
fixes**; a repair there can move goldens.

This adds a gate; it deletes no complexity. Say so in the commit message.

## Phase E — the deduplications (one push, several commits)

§4. This is where the recurring defect actually stops.

- **E1** `contextForMAC` + demote `hasVarRatio` to a local (§4.1 + §4.2,
  bundled). **Four sites unify, not six** — `:950` is a different rule and
  `:4808` must keep its own derivation or golden P4's `_c` columns move.
  `isERCandidate` **cannot** be a `SAIGEClass` method (`g_cfg` is file-static).
  Unify `isSingleVarianceRatio` on `<= 1` **and** reject a model with no `null`
  rows — `-DARMA_NO_DEBUG` means `varRatio(0)` reads off the end silently. That
  guard is the part with user-facing value.
  **Gate: full both-modes PLUS `SAIGE_TSAN_WORKLOADS=region bash harness/tsan.sh
  clean` and a plant arm** — this touches the region prepares' frames.
- **E2** one `writeSingleVariantRow` over `Marker::MarkerResult`, with the region
  path's SoA filling a local at index `k` (§4.3). Land C3's dead-parameter
  deletion first. The unified function must return **both** counts; the
  in-group writer's trait type must keep coming from **R**.
- **E3** `formatPvalOrLog` + `formatLogP` only (§4.4). **Reject
  `finishScoreTest`** — the three guards are not the same predicate, and folding
  them would delete the `DBL_MIN` guard from both fast paths. If nothing else:
  replace two of the three 13-line #44 comment copies with a one-line pointer
  (~26 lines deleted, zero instructions changed).
- **E4** the BURDEN emitter's six loops → one, behind a file-local
  `struct WeightSet` (§4.5). Declare it in `Main.cpp`, **not** `FlexRV.hpp`.
  Carry the whole `#17/#118/#129/#130/#133/#134` comment block forward and
  update `src/Main.cpp:3531` + LEDGER #130's WHERE, both of which say "six
  emission loops".
- **E5** `struct URCellModel : CellResult` — **inheritance, not composition**
  (composition forces rewriting ~114 member accesses) (§4.6). `G1` must **not**
  go in the base: the two spellings hold different quantities, which is the
  #118/#133/#135 shape.
- **E6** R-side: merge `ReadModel`/`ReadModel_subsample` on `subSampleFile`
  (§4.9) — **add value-pinning tests for the subset-only branches first**, since
  the only e2e coverage compares two runs that both traverse the changed code;
  move the double-fitted covariate GLM inside its guard; replace `getNeff`'s
  stale LOCO copy with `updateChrStartEndIndexVec`; have `checkGenoInput` return
  its resolved paths; collapse the region `.index` writer onto
  `writeOutputFileIndex`.

## Phase F — `Main.hpp` (one gate, clean build)

§5.1, 2/2 verified. 303 lines, one includer, **eight declarations disagree with
their definitions today** — every one a function the multi-model or FlexRV work
changed. A differing parameter list is a new *overload*, so the drift is
structurally uncatchable, and only the dead declarations *can* drift.
`RcppExports.cpp` already carries the compiler-checked declaration for all 24
exported ones — the header is a third copy.

Delete it; move the four load-bearing declarations into `Main.cpp` above first
use. **`Main.hpp:7` is the tree's only includer of `MarkerLoop.hpp`** — that
include must be *added* to `Main.cpp`. Doc sweep in the same commit (CLAUDE.md,
`phase5-multiweights-design.md:73`, `flexrv-plan.md:128/217/369`,
`next.md:1219`), **append** to #5's STATUS rather than re-opening it (that would
break `0 open of 142`), and add a `docdrift.py` rule.

Offer the owner the cheaper alternative: re-align the eight in place. Same gate,
no doc edits, closer to upstream.

## Phase G — the god functions (deferred; 4–6 gates)

§6. Do **not** start until the FlexRV work settles.

`SAIGE.Region` (1,251 lines) and `mainRegionInCPP` (2,501 lines) are **both
larger than upstream's after this branch extracted 675 lines from one of them**.
Every new feature has landed inside those frames rather than beside them. Until
that reverses, the Phase E tax regenerates faster than it can be paid down.

- **G1** `SAIGE.Region`, sliced in this order, one gated commit each:
  `regionStateFor` (`:1771-1870`, clearest win) → `regionTestPlan`
  (**must return `r.corr`** — `SKAT_Check_Method` promotes 0 to the 7-element rho
  grid) → `regionDoneMessage` → `composeWeightSets` scoped to `:1181-1253`
  **only** → then the accumulator as `scan <- new.env()` with **eleven** fields.
  **Do not** unify the chunk writer (three distinct `(append, col.names)` rules;
  #59's asymmetry is deliberate; no test combines markerList with multi-chunk, so
  a wrong unification passes green). **Do not** reorder the file-opening sites.
- **G2** `mainRegionInCPP`: take **only** the mechanical type hoist — ten
  function-local types verbatim into `src/RegionLoop.hpp`. No `Work` struct
  (89 frame locals), no phase lift (`:3451`'s own comment records armadillo's
  partial-accumulator unrolling as "exactly the kind of 1-ULP change that moves a
  printed golden"), no new TU. Moving the two `Worker` structs would mandate a
  full region TSAN clean **and** plant run — hours, not seven minutes.

## Leave alone (§8)

Do not re-raise: `src/CCT.cpp` / `R/CCT_modified.R`; the `dbeta` R-vs-C++ split
(**normative**, design doc D2); `region_postprocess`'s K=1/K>1 two-world arms
(the recorded byte-drift mitigation); `.singleAssoc.txt_temp` (**document, do not
delete** — it is the only byte gate on the fastTest first-pass single-variant
writer, the LEDGER #50 surface); `NAMESPACE`'s `exportPattern`;
`checkArgsBool`/`checkArgsNumeric`; `checkArgs.R`'s wrappers, `getMarkerPval`'s
body, `SPA_binary`'s `Isconverge`, `Makevars.win` — all upstream verbatim.

The principle: **do not restructure inherited code for readability alone.** The
ability to answer "is this upstream's behaviour?" with `git show main:` is a
live, load-bearing tool — CLAUDE.md's allele-convention section and LEDGER #51's
decision both rest on it — and every gratuitous rewrite of an inherited body
spends it.

## Verification

- **A, B**: `python3 harness/tools/ledger.py lint` + `python3
  harness/tools/docdrift.py` (seconds). Plus A1's two-way plant drill above.
- **C, D, E, F**: `source harness/env.sh && bash harness/prepush.sh` — full
  gate, both modes, ~7 min per push. `.hpp` touched in C, F → `rm src/*.o` first.
  E1 additionally needs the region TSAN clean + plant arms.
- **Before every gate**, write the golden-movement prediction down. For C, D, F
  it is "no golden moves"; for E it is the same, and E2/E4 carry explicit
  byte-identity arguments (`x * 1.0 == x`; the count source is *selected*, not
  assumed) that belong in the commit message.
- **G**: one gate per slice. G1's step 2 can move every SKAT-O golden if
  `r.corr` is dropped — state that prediction before running it.
