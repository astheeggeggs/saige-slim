# C6 gate state (Phase 5 flexRV port — R multi-weight post-processing)

Snapshot 2026-08-30, branch `strip`, HEAD fb122646 (C5) + uncommitted C6 tree:
`R/SAIGE_SPATest_Region.R`, `R/SAIGE_SPATest_Region_Func.R`,
`tests/testthat/test-step2-multiweight.R`, `tests/testthat/test-unit-acatv.R`
(new), `tests/testthat/test-unit-groupfile.R` (extended).

## Review fix round: four confirmed defects, all fixed in place

The 3-lens adversarial review workflow (wf_88b99700-929, 11 agents) confirmed
six findings = four distinct defects. Fixes applied directly (the C6 subagent
was no longer addressable):

1. **NULL-`pval.Region` crash** (`SAIGE_SPATest_Region.R` ~:436, K≤1
   regression). `rep(NA, ncol(NULL))` errors when the pooled-Cauchy branch
   runs but no (anno, max_MAF, set) cell emitted a row (e.g. every flagged
   marker's p-value NA, or a weighted URV collapse all-zero). **Fix:**
   `is.null(pval.Region)` guard reproducing HEAD's legacy `cctVec`
   construction verbatim (every `get_CCT_pvalue(NULL)` is NA;
   `rbind(NULL, cctVec)` is a 1-row character matrix) so the run continues
   exactly as at HEAD; multi-weight reuses the same degenerate row.
   Rationale: continue-not-crash, byte-faithful in the legacy-reachable corner.
2. **False comment on the named Cauchy row** (~:440). Claimed names prevent
   positional filling — false: `rbind.data.frame` ignores names on an
   atomic-vector row and binds positionally. **Fix:** comment now states the
   real invariant — bytes are correct because `names(cctRow)` is assigned from
   `colnames(pval.Region)` IN ORDER; never reorder fields or add one via `$<-`
   (appends at the end, silently shifting later columns).
3. **Unvalidated weights crash the run in CCT** (K>1 only; one verifier rated
   it major). Negative weights → CCT `stop()`s mid-run; all-zero weight cell →
   NaN normalisation → `missing value where TRUE/FALSE needed`. **Fix:**
   parse-time non-negativity check in `checkGroupFile`, gated
   `nWeightSets > 1` with set/region in the message. K=1 files keep the
   historical no-sign-check behavior, bug-for-bug (weights enter only
   multiplicatively there; no ACAT-V at K=1).
4. **exp-underflow zeroes ACATV/ACATO** (design-inherited from D3).
   `exp(log p) == 0` for p < ~1e-308, and CCT's small-p `1/p` branch
   overflows on denormals — rows reported exactly 0. **Fix:**
   `get_ACATV_pvalue` drops zero-weight markers (→ NA if none survive — a
   0-weight marker contributes nothing, and this also closes the all-zero
   CCT abort at the cell level) and clamps p at `.Machine$double.xmin`, so
   extreme rows report tiny-nonzero instead of exactly 0.

## New fast-tier assertions (6, all passing)

`test-unit-acatv.R` (new): zero-weight markers dropped (2 assertions:
equals CCT of survivors; weight-0 never influences the value); NA when
nothing survives (all-zero weights, all-NA p, empty input); underflow
survival (`exp(-900)` → finite 0<p<1e-300; denormal 1e-315 → finite >0;
NA-p still dropped). `test-unit-groupfile.R` (appended): K=1 negative
weights still accepted (legacy preserved); K=2 negative weight fails fast
naming set CADD / region GENE1.

## RESOLUTION (appended after the snapshot below)

The relaunched full tier came back **FAIL 0 | WARN 0 | SKIP 0 | PASS 496**,
goldens byte-identical → **C6 committed as 02bb7832 and pushed to `slim`**
(steps 2–4 below done: plan file + memory updated; CI watch running in the
background). Step 5 done: the C7 implementation agent (BURDEN writer Weight
column) is in flight. The section below is the historical snapshot.

## Gate status

- Fast tier after fixes: **FAIL 0 | WARN 0 | SKIP 19 | PASS 380** (new tests
  executed and passed; installed-namespace suite, install done first).
- Canonical gate handle **brns6hb6c**: build half completed
  ("clean install OK (fb122646)") but the full tier was **killed at launch —
  no verdict**. Relaunched as a fresh background full-tier run on the same
  tree (build not repeated: clean build already succeeded on this exact tree).
- Pre-fix baseline for reference: the original C6 gate (b0amkhqm9) was green
  at FAIL 0 | PASS 486 before this fix round.

## Remaining steps (ordered, with dependencies)

2. **Commit C6** ← full-tier verdict green. Files: the four R/test files
   above. Message cites the gate result + the review round (4 defects fixed).
3. **Push to `slim` + check CI** ← step 2 (numeric-mode standard tier,
   ~3m20s historically).
4. **Update plan file + memory** (`rosy-munching-pudding.md`,
   `saige-strip-project.md`: C6 landed incl. review-fix round) ← step 2.
5. **Launch C7** (BURDEN writer/header Weight column; lifts C6's K>1+BURDEN
   stop; K=1 emits old bytes) ← step 2. Ultracode is off — single edit agent,
   review workflows only where risk warrants.
