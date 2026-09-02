<!-- Extracted from the architecture-panel agent's result field; the raw
     JSON is preserved at harness/work/docs-raw/ (gitignored). -->

# Step-2 threading verdict

> Judged architecture panel for an in-process multithreaded step-2 marker loop with byte-identical output

## Status (updated 2026-09-02) -- READ THIS BEFORE USING THE PLAN BELOW

The document below is the original judgment, verified against `64b0f677`. Two of
its commits have since landed and one of its stated gates was wrong.

- **S0 DONE** (`b87fcaba`). Goldens P27-P33 close the F8 coverage gap: the
  marker-path sparse `scoreTest` + PCG, the `scoreTestFast` arm, the fastTest
  re-test block, the Firth arm and the ER **resampling** arm each have a golden
  now, with the per-arm counts measured on an instrumented build and recorded in
  `harness/make_goldens.sh`. Two corrections to S0 as specified: (a) the sparse
  arm needs a step-1 model with BOTH `sparse` and `null_noXadj` variance-ratio
  rows, which only `--useSparseGRMforVarRatio=TRUE` produces, hence the new
  step-1 golden `binary_sparseGRM_vr`; (b) **`--max_MAC_for_ER=25` gates
  nothing** -- it is byte-identical to P1 on the bundled fixture, as are 20, 50
  and 100. The ER golden that gates S2 is **P33, at `--max_MAC_for_ER=200`**,
  the smallest value on this fixture that admits a marker with >= 21 carriers
  and so reaches the resampling arm at all.
- **S2 DONE** (this commit). `SL_setseed(1)` in `ComputeExact::Run` is now
  conditional on some `m_IsExact[i] != 1`, i.e. on a resampling arm actually
  running. Measured by tracing `base::set.seed` with both builds installed side
  by side: at the default `max_MAC_for_ER=4` the call count on golden P1's
  configuration goes **2 -> 0**, and at 200 (golden P33) **3 -> 1**, with all
  four runs reproducing their goldens byte for byte. F4 -- "the single highest-
  risk fact in the whole exercise", the R API entered from every ER marker -- is
  therefore closed for the default configuration and for the resampling one:
  the only remaining R-API entry is on markers that genuinely draw, and S6's
  deferral of those is still required.
- **F8's premise is spent**: the two per-marker allocation items C3 excluded for
  want of coverage landed as `8959b5bd` (paired A/B: 1.9%, p = 0.023, on the new
  `s2_marker_all_sparse` bench workload).
- **S1 DONE** (this commit). `SPA.cpp`/`SPA_binary.cpp` carry their root and
  saddle results in POD structs; the exported functions are thin wrappers that
  still return `Rcpp::List`, and `RcppExports.{cpp,R}` regenerate byte-identical.
  Measured with `Rprofmem(threshold = 0)` over a whole-genome scan (128,868
  markers, 5,572 of them SPA-corrected): R-heap page allocations **2,958 ->
  1,736**, a 41% cut, with every large-allocation bucket identical between the
  arms -- the change removes small-pool traffic and touches nothing else. It is
  also, unexpectedly, the largest step-2 speedup measured so far: **11.8% on the
  whole-genome marker scan** (16 paired reps in two blocks with the arm order
  swapped between them, t = 7.53 on 15 df, 95% CI 8.4-15.1%), i.e. ~74 us per
  SPA-corrected marker. Consistent with the profile in `harness/bench/NOTES.md`,
  which put `SPA_fast` at 15% of the marker loop.
- **S3 DONE** (this commit). `TestContext` carries the three per-marker values
  that were mutable members (`m_varRatioVal`, `m_flagSparseGRM_cur`,
  `m_isnoadjCov_cur` -- exactly F2's inventory); `getMarkerPval`, the three
  score tests, `getadjGFast`, `getPCG1ofSigmaAndGtilde` and
  `fast_logistf_fit_simple` are `const`, and the marker loop builds a context per
  marker instead of writing the model. The region path still drives the members
  through the setters and snapshots them with `currentContext()`, so it is
  unchanged by construction. **F3 (LEDGER #32) is fixed** at the landmine:
  `selectVarianceRatio` throws on an uncovered MAC instead of leaving the
  previous marker's ratio in place. The compiler was the audit -- const-ifying
  surfaced every remaining write, and the only ones were the SPA kernels and
  `fast_logistf_fit_simple` taking `arma::vec &` for read-only inputs, now
  `const arma::vec &`. One dead parameter fell out: `getMarkerPval`'s
  `t_isSparseGRM` was passed by every call site and never read (the body used the
  member), so a caller could disagree with the code and never know.
- **Line references in the plan below predate S0-S3 and have drifted.** Anchor
  on names, not numbers. Current as of `90f21c77`: the marker loop is
  `Main.cpp:244-576` inside `mainMarkerInCPP` (`:155`), the writer call that
  follows it is `:578`, and the per-marker slot vectors S4 folds into
  `MarkerResult` (`isSPAConvergeVec` among them) are declared from `:197`.
- **Not yet done**: S4-S6 (kernel extraction, reader clone, threading), S7. The
  marker path is now free of member writes; what remains before threading is
  S4's `MarkerResult`/`MarkerScratch` (the `vector<bool>` race and the
  loop-hoisted buffers) and S5's reader `clone()`, which still throws.


All verification complete. Delivering the judgment.

---

# Judgment: threading the SAIGE step-2 marker loop

All claims below verified against the working tree at `64b0f677` (branch `strip`) unless marked otherwise.

## 1. Factual disagreements resolved (verified against live code)

**F1 — Output writes: once per chunk. All three designs correct.** `writeOutfile_single` has one call site, `src/Main.cpp:558`, after the `for(i=0..q)` loop (:240-555); results accumulate in q-length slot vectors indexed by `i`; the writer streams slots in index order, skipping `pvalVec[k]=="NA"` and computing `numtest` itself. In-loop side effects: stdout progress line (:241-243), `std::cerr` on `stod` failure (:429-432), `break` on `!isReadMarker` (:281-284, unreachable — `PlinkReader::getOneMarker` returns `true` unconditionally, `src/GenoReaderRegistry.cpp:45`). Slot-indexed parallel compute + unchanged serial writer is byte-safe by construction.

**F2 — SAIGEClass members mutated per marker: exactly three.** `m_varRatioVal` (written by `assignVarianceRatio`/`assignSingleVarianceRatio`, `src/SAIGE_test.cpp:942-1042`), `m_flagSparseGRM_cur` (:1189), `m_isnoadjCov_cur` (:1194, read back only in `Main.cpp:452,454,463`, never inside SAIGE_test.cpp). `m_varRatio` (:935 etc.), `m_XVG` (:413), `m_mu21`/`m_res1` are misnamed function locals. ER path copies `res_er=m_res, pi1_er=m_mu, resout_er=m_resout` per call before `SKATExactBin_Work` mutates `res_out` — members untouched. Conditional factors are written once pre-loop (`assign_conditionMarkers_factors`, called from `R/SAIGE_Test_main.R:~455` before `SAIGE.Marker`; R forces `is_fastTest=FALSE` under `--condition`). All three designs agree and are correct.

**F3 — Cross-marker state read: one exists; only Design 1 found it.** `assignVarianceRatio` returns `false` leaving `m_varRatioVal` stale when the MAC-category vectors don't cover; `hasVarRatio` is computed and never checked (`Main.cpp:397,452`); `m_varRatioVal` is never ctor-initialized. Unreachable with default category vectors (they cover (0,∞)), but it is a genuine previous-marker read and thread-order-dependent under threading. Design 1's resolution (explicit error, LEDGER entry) is required before threading.

**F4 — THE DECISIVE DISAGREEMENT — ER exact test calls the R API on every ER marker. Design 3 is right; Designs 1 and 2 are wrong.** `ComputeExact::Run` calls `SL_setseed(1)` unconditionally at entry (`src/Binary_ComputeExact.cpp`, ~:352, before the IsExact loop); `SL_setseed` in the R build (`src/Binary_resampling.cpp:50-66`, the LEDGER #2 fix) calls `Rcpp::Environment("package:base")` + `Rcpp::Function("set.seed")` + `GetRNGstate()` — full R heap/GC API. Call chain: `getMarkerPval` ER branch → `SKATExactBin_Work` → `SKAT_Exact` (`Binary_global.cpp:69`) → `exact.Run()`. This runs for **every** binary marker with MAC ≤ `max_MAC_for_ER` (default 4) whose `StdStat > SPA_Cutoff` — common in any rare-variant scan.
- Design 2's claim "ER exact test is thread-clean" is **false**; its plan contains no mitigation, and its threaded gate (P1, minMAC=1) would crash/UB before a redesign.
- Design 1 correctly located where random *draws* happen (resampling requires `n_total=2^k > 2e6` ⇒ k≥21 carriers, then `C(k,i) > ExactMax=1e4`, satisfied at k≥21) but its H4 defers only k≥21 markers — leaving `set.seed` via Rcpp on worker threads for every k<21 ER marker. As written, broken; its own TSAN/gctorture gates would catch it, but the design's central safety claim fails.
- Design 3's fix (seed only when some `m_IsExact[i]==0`) is correct at defaults but leaves the resampling branch (`--max_MAC_for_ER ≥ 21`, legal config) calling R API on workers. **Neither fix alone suffices; the union does** (see synthesis). Note: post-fix, each resampling marker re-seeds with `set.seed(1)` at `Run()` entry, so its draws are order-independent; deferral ordering matters only for `.Random.seed` end-state parity, not output bytes.

**F5 — RcppParallel semantics: Design 1 right, Design 3's rationale spurious.** Installed RcppParallel 5.1.11.2; `parallelFor(begin,end,worker,grainSize=1,numThreads=-1)` exists, and `resolveValue` (Common.h:16-42) returns the explicit `numThreads` whenever `>0` and non-default — the explicit argument **overrides** `RCPP_PARALLEL_NUM_THREADS`, so `harness/env.sh`'s pin does not neuter threaded test runs. Design 3's stated reason for raw `tbb::task_arena` ("worker count independent of the env var") is therefore moot; `parallelFor` with explicit `numThreads` achieves the same with the wrapper already used in step 1.

**F6 — Reader.** mmap path of `PlinkClass::getOneMarker` (`src/PLINK.cpp:182-318`) writes no members (only the fread fallback touches `m_fin`/`m_OneMarkerG4`); `closegenofile` (:175-179) munmaps eagerly and `fclose`s non-idempotently — clone lifetime must own the map (shared_ptr) and close must be hardened, as Designs 2/3 argue. `clone()` currently throws (`GenoReaderRegistry.cpp:~52`).

**F7 — Silent no-op CLI path confirmed.** `extdata/step2_SPAtests.R`: `nThreads>1` branch acts only inside `if(opt$idstoIncludeFile != "")`; empty include file ⇒ nothing runs. Defaults confirmed: `markers_per_chunk=10000`, `pCutoffforFirth=0.01`, `is_fastTest=FALSE`, `is_noadjCov=TRUE`, `max_MAC_for_ER=4`.

**F8 — Golden coverage gap confirmed.** Marker-path goldens are P1 (binary, minMAC=1/ER, moreDetails, default noadjCov ⇒ `scoreTestFast_noadjCov` only), P4 (quant conditional), P7 (chunked). Nothing gates the fast-test re-test (`Main.cpp:436-475`), marker-path sparse `scoreTest`+PCG, `scoreTestFast` (noadjCov=FALSE), or Firth. C0-first is mandatory.

**F9 — SpMat lazy cache (Design 1 H6).** Verified in the pixi build env's armadillo: `mutable MapMat cache; mutable state_type sync_state; mutable std::mutex cache_mutex` (SpMat_bones.hpp:696-703); the ctor family sets `sync_state=0` (SpMat_meat.hpp:4577-4592); `ARMA_USE_STD_MUTEX` defined (config.hpp:139-140). Concurrent const CSC reads are safe; Design 1's version number ("Armadillo 15.4.2") is wrong but the underlying facts hold.

**F10 — Design 2's dead-code deletion is a live-path break.** It proposes deleting the 1-arg `assignSingleVarianceRatio` as "only called by uncalled `set_varianceRatio`" — but `mainRegionInCPP` calls it at `Main.cpp:1029`. Compile-time catch, but an error in the design.

**F11 — Plan alignment (Design 3's claim).** Verified: `~/.claude/plans/rosy-munching-pudding.md:69` reserves "TestContext + const-ified getMarkerPval/scoreTest*, clone() with shared metadata" for Phase 4. LEDGER has 16 entries; #17/#18 are free as Design 1 assumes. `SL_out()` is indeed never called.

## 2. Scores (1-10)

| Criterion | D1 minimal | D2 clean | D3 throughput |
|---|---|---|---|
| (a) byte-identity + determinism safety | **7** — sound skeleton (slot writes, serial writer, explicit-numThreads verified, break/progress/Firth handled, drills); but H4 ships R-API calls on workers for every k<21 ER marker | **6** — same sound skeleton and the best gate-controls discipline, but declares the ER path thread-clean; no mitigation exists in the plan | **8** — only design that found the unconditional `SL_setseed`; fix correct at defaults, incomplete for `max_MAC_for_ER≥21`; optional per-thread formatting adds gated risk |
| (b) mutable-state isolation | **9** — complete inventory incl. the latent VR-carry (the only real cross-marker read) with a resolution; region path untouched via wrappers | **7** — clean API (`firstPassParams`/`retestParams`, private members) and a good two-pass carry-over audit, but misses the VR-carry consequence and deletes a live-path function (F10) | **7** — correct three-member core and const overloads; thinner analysis (no VR-carry, wrappers read members) |
| (c) speedup at biobank N | **7** — parallelizes the dominant cost; no bandwidth/allocator thinking | **7** — same; explicit Amdahl envelope stated | **8** — contiguous-block streaming, madvise, allocator plan, ER copy trim, DRAM analysis; only design reasoning at N=500k scale |
| (d) incrementality/reviewability | **9** — 7 gated commits, byte-argument per commit, LEDGER entries, files listed | **9** — equally strong; cleanly quarantines the CLI-behaviour change as a separate decision (C7) | **8** — good sequence; raw TBB `task_arena` is novel review surface vs. the already-used RcppParallel wrapper, on a moot rationale (F5) |
| (e) effort | **8** — smallest diff: shared_ptr-wrapped whole PlinkClass, no member migration, no writer change | **6** — PlinkShared split + private-member migration + API redesign | **6** — BedShared split + ER trim + seed hoist + optional formatting commit |
| **Overall** | **8.0** | **7.0** | **7.4** |

Design 1 wins on structure, inventory depth, and reviewability despite the H4 error, because the error is repairable by a one-commit graft from Design 3 without touching the skeleton; Design 2's equivalent gap has no in-plan repair and its extra API surface buys little; Design 3's unique contributions are exactly the graftable parts.

## 3. Recommended design: D1 skeleton + grafts

**Skeleton (D1):** `MarkerResult`/`MarkerScratch`, `testOneMarker` as verbatim extraction of `Main.cpp:244-546`, `TestContext {useSparseGRM, noadjCov, varRatioVal}` + const overloads with legacy mutating wrappers (region path zero-diff), unchanged serial `writeOutfile_single`, `RcppParallel::parallelFor` with explicit `numThreads` (F5), serial path literally serial (no TBB init in mclapply children), default OFF, `--nThreads` contract untouched.

**Grafts:**
- **G1 (from D3, mandatory):** conditional `SL_setseed` in `ComputeExact::Run` — seed only when some `m_IsExact[i]==0`. **Plus D1's own serial deferral** for markers that would resample (`binary && MAC≤cutoff && nNonZero≥21` — exact criterion verified, F4), run ascending on the main thread. Union closes the ER hazard at defaults *and* at `--max_MAC_for_ER≥21`.
- **G2 (from D2):** the gate's own positive/negative controls in C6 — a deliberately-raced scratch build (shared workspace) must fail the byte comparator and TSAN; 1-thread-through-TBB must pass. (Matches the repo's drill-verified-suite culture.)
- **G3 (from D2):** name the param builders `firstPassParams(MAC, singleVR)` / `retestParams(MAC, singleVR)` as const methods rather than open-coded locals — same semantics, reviewable transcription of `Main.cpp:389-398/440-455`. Do **not** take D2's member-privatization/deletion (F10 risk class); defer to post-threading cleanup.
- **G4 (from D2/D3):** shared_ptr ownership of the mmap so `closeGenoFile` cannot dangle clones (structural, not convention); idempotent `closegenofile`. Keep D1's minimal "share the whole PlinkClass" (no PlinkShared split) — the mmap path is member-write-free (F6), so per-clone state is only the fallback `FILE*`/G4 buffer.
- **G5 (from D3):** contiguous grain (~64) with a simple partitioner for sequential .bed streaming; `madvise` optional and output-invariant. Skip D3's optional parallel formatting (irrelevant at biobank N, byte-risk at review).
- **G6 (from D1):** LEDGER #17 (VR-category gap → explicit error before threading; stale read is thread-order-dependent) and #18 (`--nThreads>1` no-op).

**Staged commit sequence** (each gated: `compileAttributes && rm -f src/*.o && R CMD INSTALL . && bash harness/run_tests.sh standard`, byte mode; header-change commits always clean `.o` per repo CLAUDE.md):

| # | Commit | Content | Extra gate |
|---|---|---|---|
| S0 | goldens first (no src change) | union of D1-C0/D2-C0: marker goldens for fastTest re-test + sparse (needs `example_binary_sparseGRM` model with sparse VR rows), Firth (`--is_Firth_beta=TRUE --pCutoffforFirth=0.05`), `--is_noadjCov=FALSE`, quant unconditional, `--max_MAC_for_ER=25` | double-run each for run-to-run stability |
| S1 | de-Rcpp SPA kernels | `Rcpp::List` → POD structs in `SPA.cpp`/`SPA_binary.cpp` (values round-trip exactly); exported List functions become thin wrappers; unreachable `Rcpp::stop`s → `std::logic_error` | — |
| S2 | ER thread-safety (G1) | conditional `SL_setseed`; LEDGER #2 annotation | ER25 golden byte-identical |
| S3 | TestContext + const core (G3, G6) | const overloads of `getMarkerPval`/`scoreTest*`/`getadjGFast`/`getPCG1ofSigmaAndGtilde`; legacy wrappers; `firstPassParams`/`retestParams`; VR-gap → error (LEDGER #17) + characterization test | region goldens P2/P3/P5/P6 unchanged |
| S4 | kernel extraction | `testOneMarker` + `MarkerResult` (fixes the `vector<bool>` race by construction) + `MarkerScratch`; writer takes `const vector<MarkerResult>&`; still serial; progress line kept | — |
| S5 | reader clone (G4) | `shared_ptr<PlinkClass>` (or equivalent) shared by clones; map dies with last owner; idempotent close; const `getOneMarker` | unit test: clone vs primary byte-equal incl. LEDGER #16 mid-stream row-0 pattern; close-primary-then-read-clone |
| S6 | threading, default OFF | `nThreads_marker` knob (R arg, no CLI change); explicit-numThreads `parallelFor`, grain ~64 contiguous (G5); resampling-ER deferral post-pass; break-truncation; progress re-emit; Firth int sums | all marker goldens at 1/2/4/8 threads × 5 repeats byte-identical; TSAN build; gctorture run; G2 positive/negative control drills; mclapply path unchanged; bench + maxRSS at 1 vs 4 |
| S7 | separate decision | route the silent `--nThreads>1` no-op path to in-process threads (LEDGER #18); revisit default-on | CLI behavioural-change review |

**Member-relocation table (final):**

| Member | Current writes (per marker) | New home |
|---|---|---|
| `m_varRatioVal` (hpp:34) | `SAIGE_test.cpp:942-1042` | `TestContext::varRatioVal`, produced by pure `selectVarianceRatio(MAC, sparse, noXadj, bool& has) const` (body of the 3-arg selector, same comparison order); wrapper `assignVarianceRatio` stores to member for region path; gap ⇒ error (S3) |
| `m_flagSparseGRM_cur` (hpp:53) | `:1189` via setter, from `Main.cpp:390/392/443/447` | `TestContext::useSparseGRM`; reads at `SAIGE_test.cpp:488,503,508,635,844,920` (+ `scoreTest` :165 region-only paths) → `ctx`; region setters keep writing the member, wrappers forward |
| `m_isnoadjCov_cur` (hpp:56) | `:1194`, from `Main.cpp:449` | `TestContext::noadjCov` (member retained solely for wrapper parity; read only in Main.cpp) |
| loop-hoisted buffers `t_GVec/gtildeVec/t_P2Vec/G1tilde_P_G2tilde_Vec/index vecs` (`Main.cpp:197-217`) | reused across markers, never read-before-write | `MarkerScratch`, one per thread |
| `std::vector<bool> isSPAConvergeVec` (`Main.cpp:199`) | bit-packed race | `MarkerResult::isSPAConverge` (per-slot bool) |
| reader `m_fin`, `m_OneMarkerG4` | fread fallback only | per-clone; everything else shared-immutable behind shared ownership of the map |
| conditional factors `m_*_cond`, model matrices, `m_spSigmaMat`/`m_diagSigma` | pre-loop only | stay on model, accessed const (assert `sync_state==0` at construction) |

## 4. Claims I could not verify

- `obj_cc$res.out` is length-0 in real models (D3 commit 5's rationale) — depends on step-1 artifacts; verify before the ER copy-trim commit or drop the trim.
- OpenBLAS concurrent-single-threaded-caller safety and version claims (D1 H7, D3 risk 2) — empirical gate only.
- Rmath `warning()` unreachability for in-domain `pchisq/qnorm/pnorm` (all three assert) — call sites verified guarded, but nmath sources not audited; D2's `options(warn=2)` + TSAN gate is the right check.
- Exact exception-type propagation through TBB (D1 H9 / D2).
- RcppParallel-bundled TBB version "2023.0 / interface 12180" (D3) — not checked.
- Bench figures (9.5 s / ~40 µs/marker / profile ordering) — taken from the brief as given.
- D1's "Armadillo 15.4.2" — version wrong (env has 15.2.x-line RcppArmadillo); the SpMat facts it cites are nonetheless verified (F9).

**Bottom line:** adopt Design 1's skeleton and commit discipline, graft Design 3's ER-seed fix (the one hazard only it found) combined with Design 1's resampling deferral, Design 2's gate controls and param-builder naming, and shared-ownership reader teardown. The single highest-risk fact in the whole exercise — R-API entry from the ER exact test on the default path — was misdiagnosed by two of three designs; commit S2+S6's deferral closes it in both default and non-default configurations.
