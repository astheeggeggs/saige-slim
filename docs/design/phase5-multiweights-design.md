<!-- Extracted from the design session's final message; the raw JSONL
     transcript is preserved at harness/work/docs-raw/ (gitignored).
     This is the D1-D10 specification referenced by the working plan. -->

# Design: porting multi-weight region tests ("flexRV") onto the `strip` branch

## 0. Ground truth recap (verified in both trees)

**Target invariants that must survive at defaults** (all verified in code/tests):

- `checkGroupFile` returns exactly `list(nRegions, is_weight_included)`; pinned by `expect_identical` in `tests/testthat/test-unit-groupfile.R:19-32` for both the unweighted and 3-line weighted formats.
- `SAIGE.getRegionList_new` returns per-region `WEIGHT` as a character vector; unweighted sentinel `WEIGHT = c(0)` (`R/SAIGE_SPATest_Region.R:1393-1395`), pinned at `test-unit-groupfile.R:56`.
- Four group goldens (P2 fastTest, P3 conditional, P5 BURDEN, P6 nofast) are all **unweighted**; `test-invariants-goldens.R:78` asserts exactly **1 Cauchy row per region** and recomputes ACAT per p-column to 1e-12; `:86-105` pins SKAT-O sandwich bounds and Cauchy-row NA layout.
- P1Mat/P2Mat/VarMat rows for real (non-UR) markers are **unweighted** (`src/Main.cpp:1176-1177`); all SKAT/SKAT-O weighting is R-side (`AnnoWeights`/`weightMat`/`wStatVec`/`wadjVarSMat`, `R/SAIGE_SPATest_Region.R:389-455`), duplicated near-verbatim at `:673-901` for the fastTest re-run.
- But the C++ side is **not** weight-free: `genoSumMat` accumulation (`Main.cpp:1214`, `:1384-1397`), the URV collapse `max(w0*GVec)` + `genoURMat_noweights` (`:1282-1293`), the BURDEN j/m loop (`:1676-1836`), the conditional cross-term `G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * w0` (`:1171`, `:1465`), and ER gating `(!isWeightCustomized || isEqualWeights)` (`:1442`, `:1734`) all consume weights.
- R hardcodes `dbeta(MAFVec, 1, 25)` at `:440`/`:720` regardless of `weights.beta` (which only reaches C++ via `setAssocTest_GlobalVarsInCPP`). This upstream asymmetry must be preserved verbatim on the default path.
- Both bundled weighted fixtures are unusable on the PLINK path: `group_new_chrposa1a2_withWeights.txt` keys on `chr:pos:A:C`, and `groupFile_geneBasedtest_withWeights.txt` is the retired `;weight`-suffix format. **A new rsID weighted fixture is required.**

**Fork semantics to adopt** (from `flexrv-fork`, branch `multi-weights`): weight lines tagged `weight` (set name "WEIGHT") or `weight:NAME`/`weight,NAME` (token 2 split on `[,:]`); `checkGroupFile` returns `numberofWeightlists`/`nameofWeightlists`; per-region result rows become (anno × maxMAF × weightSet) with jmr layout `jmr = j*q_maf*q_weight + m*q_weight + r` (weight fastest); `weights.beta` becomes a list of pairs producing sets named `Beta_a_b`; `is_no_weight_in_groupTest=TRUE` appends an all-1s set `No_Weights`; per-row ACAT-V + ACAT-O; single pooled Cauchy row.

**Fork defects catalogued during exploration — none to be replicated** (beyond those in the brief, I confirmed three more):
1. `genoSumcount_noweight(jmr) += genoURVec(...)` in the UR collapse accumulates *weighted* dosages as the "unweighted" allele count (fork `Main.cpp:~1590`), so BURDEN pseudo-marker MAF/MAC is wrong for custom sets — the target's separate `genoURMat_noweights` (jm-indexed) must be kept.
2. `get_SKAT_pvalue_Burden_SKAT_ACATV` has a precedence bug: `SE_BURDEN = abs(BETA_BURDEN/qnorm(Pvalue_BURDEN)/2)` (fork `Region_Func.R:~433`) vs the correct `qnorm(P/2)`; plus it duplicates two `Met_SKAT_Get_Pvalue` calls where target's `get_SKAT_pvalue` (r.corr=0) already returns all three p-values in one call.
3. Positional WEIGHT-column extraction after the `merge()` (`RegionData[posSNP, 4:(nline_per_gene+1)]`) survives only because merge happens to move `SNP` to column 1; a user weight-set named `CHROM`/`POS`/`genoIndex` would collide with markerInfo columns. Extraction must be name-based with sanitized internal names.
Also confirmed: fork's `w0G2Mat_cond` per-set loop is immediately clobbered by the old single-vector outer product two lines later (fork `Main.cpp:1141-1148`); the fastTest re-run lives in a never-invoked `.f = function()` closure and still contains a stale 20-arg (VCF-era) `mainRegionInCPP` call; `saveRDS(outList, "outList.rds")` after every region.

---

## 1. Core architectural decisions

### D1 — jmr extension lives in C++ (like the fork), not R-only. Argued:
An R-only extension (loop over weight sets around the existing post-processing) would suffice **iff** all weight influence were the R-side `AnnoWeights` outer product. It is not:
- The URV collapse `max(w_r * GVec)` produces a **different pseudo-genotype vector per custom set**, whose single-variant score/SPA test, P1Mat/P2Mat rows, and VarMat entries differ per set.
- `genoSumMat` (binary Phi_cc rescale `g.sum`, and the whole BURDEN path) is weight-baked per column.
So the fork's layout is right in kind: `q_anno_maf_weight = q_anno*q_maf*q_weight` columns for `genoSumMat`/`genoURMat`/group vectors, `q = q0 + q_anno_maf_weight` rows for per-marker outputs, `annoMAFIndicatorMat` of `q × q_anno_maf_weight` (real-marker row i sets 1 at every jmr of its (j,m); UR row `q0+jmr` sets 1 only at its own jmr). Real markers are still read and score-tested **once** — no marker re-read, P1Mat/P2Mat rows for real markers shared across sets. The degeneracy requirement is the design keystone: **at `q_weight==1`, `jmr ≡ jm` and every loop, index, and float operation must be the exact current sequence** (this is what makes byte-identical goldens provable rather than hoped).

### D2 — Per-kind w0 semantics table (the fork got beta/UR interplay wrong vs the target; this table is normative):

| weight-set kind | rare marker w0 | UR collapse | UR contribution to genoSumMat | R-side AnnoWeights |
|---|---|---|---|---|
| custom (file) set r | `w_r(i)` | `max(w_r(i)*GVec)` | `* 1` | `c(w_r[kept], rep(1, nUR))` |
| beta set r `(a_r,b_r)` | `pdf(beta_r, MAF_i)` | `max(GVec)` (unweighted) | `* pdf(beta_r, MAF_collapsed)` | `dbeta(MAFVec, a_r, b_r)` |
| `No_Weights` | 1 | `max(GVec)` | `* 1` | `rep(1, n)` |

Single-set default (K=1 implicit beta set) keeps the literal hardcoded `dbeta(MAFVec,1,25)` in R and `g_cfg.weights_beta` in C++ untouched — byte-identical. In multi-weight mode beta sets use their true params (correct on the *new* path; the default-path asymmetry is preserved bug-for-bug per branch rules — add a LEDGER.md note documenting the asymmetry rather than fixing it).

UR bookkeeping: keep `genoURMat_noweights` and `genoSumcount_noweight` indexed by **jm** (weight-independent), consumed by every r — fixes fork defect (1) and keeps counts identical at K=1.

### D3 — Output schema (consistent across both writers)
Multi-weight mode is active iff **K > 1** (or the new `--is_output_weight_column=TRUE`, needed by the equivalence invariant test). Then:

- SKAT-O path (R `fwrite`): `Region  Group  max_MAF  Weight  Pvalue  Pvalue_Burden  Pvalue_SKAT  Pvalue_ACATV  Pvalue_ACATO  BETA_Burden  SE_Burden  [Pvalue_cond  Pvalue_Burden_cond  Pvalue_SKAT_cond  Pvalue_ACATV_cond  BETA_Burden_cond  SE_Burden_cond]  MAC  [MAC_case MAC_control]  Number_rare  Number_ultra_rare`.
  - `Weight` sits **after max_MAF** in *both* writers (first three columns unchanged from today; fixes the fork's R-vs-C++ position inconsistency, where C++ used `Region Group Weight max_MAF`).
  - `Pvalue` **keeps its SKAT-O meaning** (never silently replaced); ACAT-V and ACAT-O are additive columns. ACAT-V = `CCT(exp(adjPVec[tempPos]), AnnoWeights[tempPos])` via the existing `CCT(pvals, weights)` in `R/CCT_modified.R:19`; ACAT-O = `get_CCT_pvalue(c(ACATV, Burden, SKAT))`. `Pvalue_ACATV_cond` is computable from the already-exported `pval_cVec`.
- BURDEN path (C++ writer): insert `Weight` after `max_MAF` in `openOutfile` header and `writeOutfile_BURDEN`, emitted only when `q_weight>1` (pass K to `openOutfile`; R knows K before the `openOutfile` call at `Region.R:176` via `checkGroupFile(detail=TRUE)`).
- Row order: for each (anno, maf), K consecutive weight rows (weight-innermost — matches the fork writer and jmr memory layout). The `m > maxMAF0` carry-forward reuses a small cache keyed `(j, r)` (as the fork did, correctly).
- **Cauchy: one pooled row per region** (recommended over per-set rows). Rationale: preserves `test-invariants-goldens.R:78` semantics (exactly 1 Cauchy row, each p-column = ACAT over that column across all rows), matches the fork's intent, and gives the flexRV omnibus for free: the Cauchy row's `Pvalue_ACATO` = ACAT over per-row ACAT-Os. The row is built as a **named one-row data.frame** (NA-filled, then assign p-columns by name), then `rbind` — eliminating the fork's hand-ordered `cctVec` misalignment defect class permanently. Per-set summaries remain derivable by filtering on `Weight`.

### D4 — checkGroupFile extension without breaking pinned tests
Add optional arg: `checkGroupFile(groupFile, detail = FALSE)`. Default return stays the exact 2-element list (both `expect_identical` pins pass verbatim). `detail=TRUE` adds `nWeightSets` and `weightSetNames`. `SAIGE.Region` calls with `detail=TRUE`. Parsing: token 2 split on `[,:]`; bare `weight` → name `"WEIGHT"`; the first region defines the global stride (2 + K lines) — same global-stride model as today — but pass 2 additionally verifies each region's weight-line **names match region 1's names in order** (the fork only checks the `weight` prefix → silent set mislabeling on permuted lines; we hard-error). Fork group files parse identically.

### D5 — WEIGHT vector → matrix through R
`SAIGE.getRegionList_new` gains `weightSetNames`; for weighted files it builds `RegionData` with internal sanitized columns `.W1..WK` (name-based extraction post-merge — fixes fork defect (3)), and `region$WEIGHT` becomes an n×K numeric matrix with `colnames = weightSetNames`. Unweighted sentinel `WEIGHT = c(0)` unchanged (pinned test intact). The annotation-filter subsetting at `:1381-1383` becomes `WEIGHT[-rmind, , drop=FALSE]`; the `setorderv(CHROM, POS)` sort at `:1292` already carries weight columns inside the data.table, so alignment is structural. `SAIGE.Region:292-296` builds the composed weight structure: `customMat (n×Kc)` from the file, `betaPairs (Kb×2)` from `weights.beta`, `includeNoWeights` flag; `weightSetNames = c(fileNames, paste0("Beta_",a,"_",b), "No_Weights")` in that order (fork-compatible).

### D6 — mainRegionInCPP signature (all three call sites via one wrapper; see D8)
Current arg 12 `arma::vec & t_weight` (of 19) becomes:
```
arma::mat & t_weight,              // n x Kc custom sets (0-matrix sentinel when none)
arma::mat & t_beta_params,         // Kb x 2 beta pairs (0-matrix sentinel)
bool t_includeNoWeights,
std::vector<std::string> & t_weightNames   // length q_weight, for the BURDEN writer
```
(21 args total; regenerate `RcppExports` per branch rule, `Main.hpp` was deleted 2026-09-14 (audit F), so there is nothing to align.) `q_weight = Kc + Kb + includeNoWeights`; degeneracy contract: today's calls map to (n×1 custom, empty beta, false) for weighted files, (0-matrix, empty beta... ) — **no**: the K=1 default (unweighted file) maps to `t_weight = 0-matrix sentinel, t_beta_params = 0-matrix, includeNoWeights=false` and C++ then treats `q_weight=1` with the *implicit* current `g_cfg.weights_beta` beta set — i.e., `isWeightCustomized=false` exactly as today. `w0_vec(r)` filled per marker per the D2 table; at K=1 `w0_vec(0)` reproduces today's scalar `w0` computation (`Main.cpp:1106-1110`) exactly.

### D7 — Conditional analysis: single conditioning-weight vector; staged delivery
Decision: conditioning-marker weights stay a **single vector** across all sets (user `weights_for_condition`, else `pdf(beta(1,25))` — exactly today's `w0G2Vec_cond`), and `SAIGE_test.cpp`/`m_G2_Weight_cond` (arma::vec) is untouched. Rationale: the fork's per-set `m_G2_Weight_cond` matrix was never wired (`.col(ci)` OOB, result clobbered); a per-set G2 weighting is a modeling refinement with no working reference. What *must* vary per set is the G1-side factor `w0`: `G1tilde_P_G2tilde_Weighted_Mat.row(i) = Vec % w0G2^T * w0` bakes `w0` in C++ (`:1171`,`:1465`). Design: when `q_weight>1`, C++ exports the **raw** `G1tilde_P_G2tilde_Mat` (no `*w0`; still `% w0G2^T`) plus skips the C++ `TstatAdjCond`/`VarMatAdjCond` computation (`:1860-1869`); R computes, per set r, `G1w_r = raw * AnnoWeightsMat[,r]` (row-scale), then `AdjCondMat_r`, `TstatAdjCond_r`, `VarMatAdjCond_r` from the already-exported `m_VarInvMat_cond`-derived pieces. At `q_weight==1` the current baked path runs unchanged (byte-identical P3 golden). Staging: **commit C9 ships a hard, clear error** (`isCondition && K>1` → "conditional analysis supports a single weight set; got K=…") satisfying constraint (4); **commit C11** (last, separable) implements the R-side per-set adjustment and adds the conditional multi-weight golden. BURDEN conditional per-jmr works as-is (per-set `genoSumVec` through `getMarkerPval`'s internal conditional machinery, weight-independent `m_*_cond`).

### D8 — De-duplicate before extending (prep commit)
The ~250-line post-processing block (`Region.R:383-628`) and its near-verbatim twin (`:667-901`) are extracted into one function, e.g. `region_postprocess(outList, ctx)` returning `list(pval.Region, cctpval)`; the three `mainRegionInCPP` call sites collapse into a local `call_mainRegion(useSparseGRM)` closure. This is a pure refactor gated byte-exact against all four group goldens, and it is precisely what prevents the fork's failure mode (their two copies diverged: one updated, one left with a dead closure and a stale 20-arg call → fastTest broken). All subsequent multi-weight R changes then land exactly once. The two blocks have three deliberate diffs to fold as parameters (the block-2 re-derivation of `noNAIndices`/`annoMAFIndicatorMat`, and comment noise); confirmed no semantic differences otherwise.

### D9 — ER gating per set
Replace the global `(!isWeightCustomized || isEqualWeights)` at `:1442`/`:1734` with a per-set predicate `erAllowed(r)`: true for `No_Weights`, true for a custom set whose weights are all 1, false for other custom sets, true for beta sets **only in the exact case it is true today** — note today `isWeightCustomized=false` (beta path) makes the gate pass, so `erAllowed(beta set) = true` preserves K=1 behavior; document that ER on a beta-weighted collapsed pseudo-marker is the current (questionable) upstream behavior, preserved.

### D10 — CLI surface
- `--weights.beta="1,25"` unchanged meaning; **new syntax** `--weights.beta="1,25;1,1"` — semicolon-separated pairs, comma within a pair (no shell-quoting hazards, backward-compatible parse: `strsplit(x,";")` then per-pair `strsplit(...,",")`, validate 2 finite numerics each). In `SPAGMMATtest`, `weights.beta` becomes the Kb×2 matrix; row 1 (or the single pair) continues to feed `setAssocTest_GlobalVarsInCPP` (marker path + cond-weight default) unchanged.
- Beta sets join the group test only when (file has no weight lines) — then the single default pair is the implicit current set — or when the user passes >1 pair. A weighted file + the default single pair behaves exactly as today (custom set only). This is the rule that makes the fork's always-append-Beta_1_25 non-compatibility go away.
- `--is_no_weight_in_groupTest=TRUE`: unweighted file → today's behavior verbatim (single all-1s set, no new columns). Weighted file → **appends** `No_Weights` (replaces today's hard error at `Region.R:207-209`; safe: previously-erroring invocations gain behavior, nothing silent changes).
- New `--is_output_weight_column=TRUE` (default FALSE): forces the multi-weight schema at K=1 — exists for users wanting stable schemas across runs and for the equivalence invariant test.

---

## 2. Staged commit sequence (each gated: `bash harness/build.sh` from clean objects + `bash harness/run_tests.sh standard` in byte mode; full tier where flagged)

**C1 — prep: de-duplicate region post-processing.**
Files: `R/SAIGE_SPATest_Region.R` (extract `region_postprocess()`; `call_mainRegion()` wrapper absorbing the 3 call sites at `:318`, `:349`, `:646`). No C++/schema change. Gate: standard byte (P2/P3/P5/P6 byte-identical is the whole point). Risk: subtle env capture in the closure — keep it a plain function taking explicit args, not a closure over mutables.

**C2 — checkGroupFile: `weight:NAME` parsing, `detail=TRUE`, validation.**
Files: `R/SAIGE_SPATest_Region_Func.R:73-188`; `tests/testthat/test-unit-groupfile.R` (new cases only — the two pinned `expect_identical`s untouched). Validation: duplicate set names, reserved names (`Cauchy`), per-region name/order match, weight-count-vs-marker-count per line. Transcribe the fork's two-pass skeleton (~70%), rewrite validation. Gate: standard.

**C3 — getRegionList: WEIGHT matrix + matrix-aware subsetting; SAIGE.Region consumes K=1 identically.**
Files: `R/SAIGE_SPATest_Region.R:1191-1425` (+ `:292-296`). `WEIGHT=c(0)` sentinel preserved; n×K matrix with sanitized internal columns; `as.numeric` path for K=1 verified equal to old vector path. Unit tests: filter-subset alignment, POS-sort alignment, name-based extraction with adversarial set names (`CHROM`). Gate: standard.

**C4 — fixtures (no goldens yet).**
New files: `extdata/input/group_new_snpid_withWeights.txt` (single unnamed weight line, rsID-keyed — reuse the numeric weights from `group_new_chrposa1a2_withWeights.txt` mapped rs1..rs100) and `extdata/input/group_new_snpid_multiWeights.txt` (two named lines, e.g. `weight:CADD` = same values, `weight:REVEL` = a fixed deterministic transform, committed literally — never generated at test time). Extend unit tests to parse both. Gate: standard (fixtures inert to goldens).

**C5 — C++ jmr core.** *(the risk commit; full tier)*
Files: `src/Main.cpp` (`mainRegionInCPP` signature per D6; `w0_vec` fill at `:1106`; rare-marker jmr accumulation around `:1198-1229`; UR collapse `:1269-1300` per D2 with jm-indexed `genoURMat_noweights`; UR-collapse chunk `:1344-1520` (P1/P2 resize `m1new = max(m1, q_anno_maf_weight)`, UR pseudo-marker names gain `:<WeightName>` suffix only when K>1); BURDEN loop `:1676-1836` gains inner r; ER predicate per D9; conditional raw-export branch stub for K>1), `src/Main.hpp` (aligned decls), regenerated `src/RcppExports.cpp`/`R/RcppExports.R`; `R/SAIGE_SPATest_Region.R` wrapper passes the new args (K=1 sentinels). K>1 allowed through C++ but R still errors before use (guard added here, lifted in C6) — so behavior at defaults provably unchanged. Gate: **build.sh (clean objects — signature change), standard byte, plus one full run**; additionally a temporary assertion test: run P2's exact config through `SPAGMMATtest` with the C4 single-weight fixture and diff against the same run made at C3's commit? Not possible read-only across commits — instead assert internal degeneracy directly: unit-test that `checkGroupFile`+`getRegionList`+`mainRegionInCPP` on the weighted single-set fixture produces `annoMAFIndicatorMat` with `q_anno_maf` columns and the same `MAC_GroupVec` as the K=1 vector path (callable from testthat via the installed package).

**C6 — R multi-weight post-processing + ACAT-V/ACAT-O + pooled Cauchy.**
Files: `R/SAIGE_SPATest_Region.R` (inside the single `region_postprocess()`: `AnnoWeightsMat` per D2; per-set `wStatVec`/`wadjVarSMat` (list); r-loop with `(j,r)` carry-forward cache; resultDF gains `Weight`/`Pvalue_ACATV`/`Pvalue_ACATO` **only when K>1 or forced**; named-data.frame Cauchy row per D3), `R/SAIGE_SPATest_Region_Func.R` (small `get_ACATV_pvalue(pvals, weights)` helper delegating to `CCT`; reuse existing `get_SKAT_pvalue` — do not port the fork's duplicate/buggy helper), `R/SAIGE_Test_main.R` (set composition per D10 rules; pass-through). K=1 path exercises literally the current code lines (branch on K, not a rewrite of the K=1 arm). Gate: standard byte + new unit tests over the fixtures (run `SPAGMMATtest` in-process against the multiWeights fixture; assert row counts = anno×maf×K + 1 Cauchy, `Weight` values, ACAT-O = CCT of its three columns).

**C7 — BURDEN C++ writer/header Weight column.**
Files: `src/Main.cpp` (`openOutfile:2267` gains `nWeightSets` param; `writeOutfile_BURDEN:2564` gains `q_weight`+names, emits `Weight` after `max_MAF` iff K>1; pooled Cauchy row with the extra NA slot), `R/SAIGE_SPATest_Region.R:176` (pass K), RcppExports regen. Gate: build.sh + standard byte (P5 golden byte-identical since K=1 emits the old header/rows verbatim).

**C8 — CLI + docs.**
Files: `extdata/step2_SPAtests.R` (`--weights.beta` pair-list parse, `--is_output_weight_column`, help text incl. `weight:NAME` syntax), `R/SAIGE_Test_main.R` (roxygen for `groupFile`/`weights.beta`/`is_no_weight_in_groupTest`), `CLAUDE.md` (one paragraph: multi-weight mode, schema trigger, fixture names), `man/` regen if the repo does. Gate: standard.

**C9 — conditional guard + negative tests.**
Files: `R/SAIGE_SPATest_Region.R` (clear `stop()` when `isCondition && K>1`), test asserting the exact message via CLI. Gate: standard.

**C10 — goldens + invariant tests.** *(full tier)*
Files: `harness/make_goldens.sh` (P8: SKAT-O fastTest binary on `group_new_snpid_multiWeights.txt`, same anno/maxMAF flags as P2; P9: BURDEN r.corr=1 on the same file; optionally P8b adding `--weights.beta="1,25;1,1" --is_no_weight_in_groupTest=TRUE` to exercise composed sets), committed goldens under `harness/goldens/step2/`, `tests/testthat/test-step2-goldens.R` (two new cases), `tests/testthat/test-invariants-goldens.R` additions:
- Cauchy pooled: still exactly 1 row/region; per-column ACAT recompute extended to `Pvalue_ACATV`/`Pvalue_ACATO` (1e-12).
- ACAT-O recompute: per row, `Pvalue_ACATO == acat(c(ACATV, Burden, SKAT))`.
- ACAT-V recompute: from `.singleAssoc.txt` p-values + fixture weights (+1 for UR rows), matched by marker membership per (anno,maf) via the markerList output or MAF cutoffs.
- Burden linearity (two-ways): fixture contains a third region (or a dedicated small fixture) with sets `w` and `2w` → assert `Pvalue_Burden` identical and `BETA_Burden` halves, `SE_Burden` halves (score-test scale invariance).
- Single-set equivalence: CLI run of P2's config with `--is_output_weight_column=TRUE` on the *unweighted* file → drop the `Weight`/ACATV/ACATO columns → identical (numeric 1e-12, or byte after column excision) to the P2 golden.
- `No_Weights` set equivalence: in the composed run, rows with `Weight=="No_Weights"` match a plain `is_no_weight_in_groupTest=TRUE` run's rows.
Gate: full (golden regeneration by definition), and old goldens byte-identical.

**C11 — conditional multi-weight (stretch, separable).**
Files: `src/Main.cpp` (K>1 branch exports raw `G1tilde_P_G2tilde_Mat`, skips baked `TstatAdjCond`/`VarMatAdjCond`), `R/SAIGE_SPATest_Region.R` (per-set adjustment in `region_postprocess`), lift the C9 error, new conditional multi-weight golden P10 + invariant (conditional run leaves unconditional columns untouched — mirror of the existing `:107-113` test). Gate: full.

Total: 11 commits (10 if C11 deferred), each independently green.

---

## 3. Fork commits: transcribe vs rewrite

| Fork material | Verdict |
|---|---|
| `checkGroupFile` two-pass parser + `[,:]` split (`fc41e71a`, port commit) | **Transcribe ~70%**; add name/order validation, wrap in `detail=` |
| `getRegionList` multi-column WEIGHT build | Transcribe the geneData `cbind` loop; **rewrite** extraction name-based |
| C++ jmr layout, `w0_vec` fill, per-set loops (`34ae5a7a`, `4a862dc5`, `3e44ae56`, `4a0a15f8`) | **Transcribe ~50-60%** (layout, loop skeletons, BURDEN `WeightName_Vec`); rewrite UR/beta interplay per D2, drop `genoURMat_noweights` removal, drop diagnostics |
| R jmr post-loop + `AnnoWeightsMat` (`31ac3b39`, `01b28be4`) | Semantics transcribable ~40%; **rewrite** into the deduped helper; drop silent `next` guards, `saveRDS`, prints; keep the `(j,r)` cache idea |
| `writeOutfile_BURDEN` + header (`ee84e03f`) | Transcribe skeleton; fix column position + conditional emission |
| `get_SKAT_pvalue_Burden_SKAT_ACATV` | **Rewrite thin** (reuse target `get_SKAT_pvalue` + `CCT(p,weights)`; fork has SE precedence bug + duplicated SKAT calls) |
| Conditional multi-weight (`m_G2_Weight_cond` matrix) | **Discard**; replace with D7 |
| PCG relative-tolerance commit `7c039806` | **Out of scope** (collides with deliberately-deferred LEDGER #6) |
| `.dockerignore`, zero-genotype catches (`968565a2`, `91428eb4`) | Skip; the zero-genotype guard is superseded by C13-style validation (all-zero weight columns rejected at parse; all-zero genotype regions already handled by target's `NumRare+NumUltraRare==0` gate at `Region.R:400`) |

Net: roughly 35-40% of the fork's diff is transcribable; the conditional path, the helper functions, and everything touching output assembly must be rewritten.

---

## 4. Risks

1. **Index skew C++↔R** (jmr formula divergence) — single normative formula in one comment on each side; unit test cross-checks `ncol(annoMAFIndicatorMat) == anno*maf*K` and MAC vector lengths; C5's degeneracy test.
2. **Byte drift at K=1 from incidental reordering** in C5/C6 — mitigate by branching on K around existing lines rather than rewriting the K=1 arm; the standard byte gate catches any slip on 4 goldens × ~250 invariant assertions.
3. **BURDEN header/writer column-count coupling** (cond × trait × K combinations) — the Cauchy row NA-slot math is hand-counted in C++; add a golden for BURDEN multi-weight (P9) and, in the invariant suite, a header-vs-row field-count check for every group golden.
4. **Memory**: P1/P2Mat UR-resize to `max(m1, anno*maf*K)` rows and `genoURMat` at n×(anno*maf*K) — for the bundled data trivial; document scaling (K multiplies UR columns, not marker reads) in CLAUDE.md.
5. **Restart/append (.index)**: schema decided per-run from K — a restart with a different weight config would append mismatched columns; add a guard (existing header sniffed on `isappend`, error on mismatch) in C6.
6. **`format(genoIndex)` leading-space quirk** and other data.table idiosyncrasies interacting with the new matrix columns — covered by C3 unit tests including the pinned `" 0"` behavior.
7. **fastTest cctpval semantics at K>1**: pooled Cauchy now includes ACAT-O-influencing columns? No — keep the re-run trigger `cctpval` = ACAT over `Pvalue` (SKAT-O) rows exactly as today, K rows simply widen the pool; document.

---

### Critical Files for Implementation
- /Users/dpalmer/Repositories/SAIGE/R/SAIGE_SPATest_Region.R
- /Users/dpalmer/Repositories/SAIGE/src/Main.cpp
- /Users/dpalmer/Repositories/SAIGE/R/SAIGE_SPATest_Region_Func.R
- /Users/dpalmer/Repositories/SAIGE/R/SAIGE_Test_main.R
- /Users/dpalmer/Repositories/SAIGE/tests/testthat/test-unit-groupfile.R (plus harness/make_goldens.sh for the new golden runs)
