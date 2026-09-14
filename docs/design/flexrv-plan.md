# FlexRV on the `strip` branch

<!-- Plan written 2026-09-06/07. Status 2026-09-09 (evening): F0, F1, F2, F3
     and F4 LANDED -- kernels (src/FlexRV.cpp, test-unit-flexrv.R), the engine
     (SET_TRANSFORM in mainRegionInCPP, R/FlexRV.R, the five --flexRV_* flags,
     test-step2-flexrv{,-guards}.R), the fixture
     (extdata/input/group_new_snpid_flexrv.txt) and goldens P58-P61 with their
     invariants. See docs/state/handoff-2026-09-09-flexrv-f2-f4.md and the
     "Implementation notes" section of docs/design/flexrv-design.md, which
     record what this plan got wrong (the fixture's in-sample N is 1,000, not
     10,000) and what was found on the way (LEDGER #117: a zero weight NAs a
     binary SKAT-O cell). 2026-09-09 (late): the end-to-end suite landed
     (test-step2-flexrv-{reconcile,cond,paths,multimodel}.R, see
     docs/state/handoff-2026-09-09-flexrv-e2e.md) and opened LEDGER #118 (a
     set's burden is tested with the cell's unweighted MAC). NOT done:
     --flexRV_output_sets (every set's row is
     always written), F5 (the N = 50k measurement, which is where the RAM/CPU
     question is actually answered -- the N = 1,000 fixture is startup-bound
     and cannot see it), F6's roxygen for the R API is done but the CLAUDE.md
     paragraph is the one in this commit series. The multi-weight engine this
     builds on is Phase 5 (C1-C11), landed. -->

## How to read the references in this file

**Every code reference here is anchored on a SYMBOL, with a line number as a
dated hint only** (`~:1956` means "near line 1956 at commit `97db5a88`,
2026-09-08"). The first version of this plan cited bare line ranges and they had
drifted by up to ~145 lines within a day -- `SAIGEWT_k` was written as `:1579`
and was at `:1745` -- which is exactly the failure CLAUDE.md's "re-resolve every
line number before using it" rule exists to prevent. Grep for the symbol; treat
the number as a starting point, not an address.

## Context

**Is FlexRV implemented? No, not as a feature.** What the repo calls "flexRV"
(`docs/design/phase5-multiweights-design.md`, `docs/state/c6-gate.md`) is the
Phase 5 multi-weight region engine, fully landed (C1-C11, commits
`0fc41971`..`94a95dd7`): K weight sets per region from `weight:NAME` group-file
lines plus `--weights.beta` pairs and `No_Weights`, rows per
(anno x max_MAF x weight set), per-row ACAT-V/ACAT-O, and ONE pooled Cauchy row
per region whose each p-column is ACAT over every row (both writers: the `CCT_cpp` call in the
BURDEN writer, `src/Main.cpp` ~:3551-3587, for BURDEN, and the `cctVec`
assembly in `R/SAIGE_SPATest_Region.R` ~:563 for SKAT-O). That engine is the substrate; FlexRV itself is absent.

**What FlexRV is** (Schwartzentruber, Fiziev, McRae, Ulirsch, Farh; bioRxiv
2025.11.04.686320, Illumina; STAAR implementation, code promised at
github.com/Illumina/flexrv which does not exist yet). Per gene, 192 fixed
weight vectors = 16 pathogenicity transforms x 12 MAF transforms, a burden test
per vector (SKAT optional), Cauchy-combined into `p_FlexRV`. No training.
Per-variant score `s` in [0,1]: LoF = 1, missense = PrimateAI-3D percentile
within the gene's possible missense variants, combined by max. Variants with
MAF < 0.1%. Transforms as quoted from the Methods (via a lossy fetch; **F0 re-pinned
them from the PDF in `docs/design/flexrv-design.md` -- read that, not this
table, before writing a kernel**):

| family | count | weight |
|---|---|---|
| score: unweighted | 1 | 1 |
| score: LoF-only mask | 1 | 1 if variant is LoF else 0 (annotation-based, not score-based) |
| score: cubic root, t in {0.2..0.9} | 8 | `(sign(s-t)*abs(s-t)^(1/3) + t^(1/3)) / ((1-t)^(1/3) + t^(1/3))` |
| score: beta, b in {1.5,2,3,5,10,20} | 6 | `dbeta(1-s,1,b)/dbeta(0,1,b)` = `s^(b-1)` |
| MAF: unweighted | 1 | 1 |
| MAF: beta, b in {2,5} | 2 | `dbeta(MAF/0.001,1,b)/dbeta(0,1,b)` = `(1-MAF/0.001)^(b-1)` for MAF<0.001, else 0 |
| MAF: hard MAC, t in {1,3,5,10,25,50,100,200,500} | 9 | 1 if MAC <= t else 0 |

`w = w_score * w_MAF`. Paper's large-gene collapse (>= 400 variants: AC <= 5
collapsed into <= 50 pseudo-variants by AC and score decile) is NOT
replicated; SAIGE's own ultra-rare collapse (`--MACCutoff_to_CollapseUltraRare`,
default 10, per set `max(w_r*G)`) plays that role and is documented as the
difference.

**Why not just write 192 weight lines?** It works today (custom sets) and is
kept as the reconciliation control, but MAC thresholds need the in-sample MAC
that only C++ knows after QC (the `MACVec.at(i)`/`MAFVec.at(i)` fill,
`src/Main.cpp` ~:2463), and group files grow
192x. Decision (user): build the sets inside SAIGE from one `score:NAME` line,
and include the N=50k measurement stage.

## Design

**Group file.** One optional line per gene, `GENE score:NAME v1 v2 ...`
(same `[,:]` split as weight lines -- `parse_weight_type`,
`R/SAIGE_SPATest_Region_Func.R` ~:86, and `parse_score_type` beside it since F1),
after the var/anno lines and before/among weight lines; values numeric in
[0,1], NA rejected. Carried like weight columns through
`SAIGE.getRegionList_new` (sanitized `SAIGESC_1`, cf. `SAIGEWT_k` at
the `paste0("SAIGEWT_", ...)` line in
`SAIGE.getRegionList_new`, `R/SAIGE_SPATest_Region.R` ~:1745) and subset with
the annotation filter. DONE in F1 as `SAIGESC_1`.
Region 1 defines the layout; every region must carry the same score name
(the same global-stride model as weights -- the "must appear in the same order
in every region" check, `R/SAIGE_SPATest_Region_Func.R` ~:294).

**CLI.** `--flexRV_score=NAME` enables FlexRV sets from that score line
(error if the file has no such line). `--flexRV_scoreTransforms` and
`--flexRV_MAFTransforms` are spec strings with the paper's defaults
(`unw;lof;cbrt:0.2,...,0.9;sbeta:1.5,2,3,5,10,20` and
`unw;mafbeta:2,5;mac:1,3,5,10,25,50,100,200,500`); `--flexRV_lofAnno=lof`
names the annotation label(s) that the LoF mask keys on. Parsed in
`R/SAIGE_Test_main.R` next to `parse_weights_beta` (~:11, unchanged), into a
`Kt x 4` transform table (scoreKind, scoreParam, mafKind, mafParam).

**Set composition** (the "weight-set composition (design D5/D10)" block,
`R/SAIGE_SPATest_Region.R` ~:862 onwards): order becomes
custom, **transform**, beta, No_Weights; names
`<NAME>_<scoreT>_<mafT>` e.g. `PAI3D_cbrt0.5_MAC10`, `PAI3D_unw_unw`.
`weightSetDescriptor` (`R/WeightSets.R`) gains kind `transform`; its R-side
AnnoWeights come from a C++-exported `transformWeightMat` (q0 x Kt, real
markers only) with `rep(1, numofUR)` appended, i.e. **custom-set semantics**
(collapse `max(w_r*G)`, UR pseudo-marker `w0 = 1` unscaled, LEDGER #17 shape).
Weights are computed ONCE, in C++, and exported; R never re-derives them.
Exception to the two-ways rule is deliberate and replaced by a kernel test
(below).

**C++** (`src/Main.cpp` `mainRegionInCPP`): new args `arma::vec t_score`,
`arma::mat t_transformTable`, `std::vector<std::string> t_lofAnno` -> new
`SET_TRANSFORM` kind in the set tables (the `SET_CUSTOM`/`SET_BETA`/`SET_NOWEIGHT`
constants and `setKindVec` fill, ~:1956-1969). At the `pr.w0(r) = t_weight(...)`
fill (~:2331) and every custom-branch site (the weighted `genoURMat ... std::max`
collapse ~:2577-2599, the UR emission from `genoURMat.col(jmr)` ~:2891-2931,
conditional raw export) read from `transformW(i, c)` instead of
`t_weight(i, c)` through one accessor, so the two kinds share code paths.
`transformW(i, c) = flexrvScoreWeight(kind, param, s_i, isLof_i) *
flexrvMafWeight(kind, param, MAF_i, MAC_i)` evaluated right after MAF/MAC are
final (the `MACVec`/`MAFVec` fill, ~:2463); `isLof_i` comes from the marker's
anno label. Kernels live in
a new `src/FlexRV.cpp` with an `[[Rcpp::export]]` wrapper for tests (add to
`OBJECTS` in `src/Makevars`, `compileAttributes`, align `Main.hpp`).
`erAllowedVec[r]` is true only for the `unw_unw` set (weights identically 1),
false otherwise, matching the custom rule in the `erAllowedVec` fill (~:1963-1969). At `Kt == 0` nothing in
the existing arms changes: every branch is on kind, the K=1 legacy arm is
untouched, so all committed goldens stay byte-identical (predicted movers:
none).

**Beta-weight independence (owner's concern, verified in code).** The
default `beta(1,25)` MAF weighting must NOT touch FlexRV sets: they replace
it. Today's custom path already guarantees this on the unconditional path
(`pr.w0(r) = t_weight(i, setColVec[r])` at `src/Main.cpp` ~:2331, R AnnoWeights
= file weights + ones for UR rows at `R/SAIGE_SPATest_Region.R` ~:218, and
the implicit beta set is composed only when no other set exists,
the "composed ONLY when more than one pair is supplied" rule in the composition
block, ~:883). Transform sets inherit that by construction (same accessor, same
branches). Two guards make it an invariant rather than a hope:
- Test: a FlexRV run with `--weights.beta="1,1"` and one with the default
  `"1,25"` produce identical FlexRV rows (numeric 1e-12); the single pair
  must be inert on transform sets. `unw_unw` == `No_Weights` rows is the
  second witness (any beta leak would show as `dbeta(MAF,1,25)` scaling).
- With `--weights.beta` given as >1 pairs, `Beta_a_b` sets are appended
  beside the transform sets as separate rows, never multiplied in; the
  pooled Cauchy then pools them too, so document that a pure `p_FlexRV`
  needs the default single pair.

The one place beta(1,25) DOES reach a custom-set run is the conditioning
markers: `w0G2Vec_cond` defaults to `pdf(beta(1,25), MAF)` when
`--weights_for_condition` is absent (the `t_weight_cond.is_zero()` fallback,
`src/Main.cpp` ~:2013 and ~:3998), one vector
shared by every set (design D7). Conditioning markers live outside the group
file and have no score, so a transform cannot be applied to them. Rule:
`--flexRV_score` + `--condition` is a hard error unless
`--weights_for_condition` is supplied (message names both flags); with it,
the D7 raw-export path (C11) serves transform sets like custom sets, gated by
`test-step2-multiweight-cond.R`'s reconciliation extended to one transform
set. Documented in F0 and CLAUDE.md as a known asymmetry, not silently
defaulted.

**Zero-weight cells** (already characterised): BURDEN omits the row and the
Cauchy pool skips it (the BURDEN writer's row-emission guard and its
`nonMissingPvalVec` filter, `src/Main.cpp` ~:3551); SKAT-O reports NA via the SKAT
failure path. Fine for MAC thresholds below a gene's smallest MAC; document.

**Memory.** `genoURMat`, `genoSumMat` (N x anno*maf*K), `P1Mat`/`P2Mat`
(max(m1, anno*maf*K) x N): at K=192 and the default 12 cells that is
2304 x N x 8 B x ~4 arrays (9 GB per array at N=500k). Print an estimate once
at region start when Kt > 0, and document the recommended FlexRV run: one
anno group (`lof:missense`), one max MAF (0.001), `--r.corr=1`. Whether UR
storage must go sparse is decided by the F5 measurement, not now.

`p_FlexRV` = the pooled Cauchy row's `Pvalue_Burden` when exactly one
(anno, max_MAF) cell is run. With several cells the row is a super-omnibus;
per-cell Cauchy-over-sets rows are a follow-up decision, not in scope.

## Commits (each gated: `bash harness/build.sh` when C++/headers change,
`bash harness/run_tests.sh standard` in byte mode, full tier where marked)

- **F0 design doc.** `docs/design/flexrv-design.md`: the transform table
  re-pinned VERBATIM from the PDF (bioRxiv blocks fetches; the user supplies
  the PDF), closed forms with their derivations, the collapse difference, the
  ER rule, the memory formula. No code.
- **F1 score line. DONE 2026-09-08.** Two decisions taken during
  implementation, both recorded because they narrow the design above:
  (1) **the score line's position is FIXED** -- immediately after `anno`,
  before any weight line -- where this plan said "before/among weight lines".
  Relaxing a position later is backward compatible for existing files and
  tightening later is not, and it keeps every offset in the parser and in
  `getRegionList_new` a constant instead of a per-region layout vector. A
  second score line, or one after a weight line, is a named error.
  (2) **`scoreName` is appended to the `detail = TRUE` return ONLY when the
  file has a score line.** Three existing tests pin that return as exactly
  four elements, not two as this plan assumed, and they had to stay verbatim;
  callers read `region_list$scoreName`, which is NULL when absent, and
  `length(NULL) == 0` is the arithmetic they want. The same rule applies to
  `RegionList[[r]]$SCORE` and to the `nline_per_gene` mismatch message, which
  keeps its historical wording when no score is in play (`test-unit-regionlist.R`
  pins it). Line references below had drifted by ~145 lines and were
  re-resolved.
- **F1 (original text).** `checkGroupFile` (`R/SAIGE_SPATest_Region_Func.R` ~:79)
  parses `score:NAME`, `detail=TRUE` adds `scoreName`; the two pinned
  `expect_identical` returns stay verbatim. `SAIGE.getRegionList_new` carries
  `SCORE`. Unit tests in `test-unit-groupfile.R` (parse, [0,1] validation, name
  mismatch across regions, score with and without weight lines). Inert at
  defaults.
- **F2 kernels.** `src/FlexRV.cpp` + `flexrvWeightsInCPP()` export;
  `tests/testthat/test-unit-flexrv.R` pins each transform to the R closed
  form on a grid (fast tier, like `test-unit-spa.R`): s=0, s=1, s=t,
  MAF=0.001 boundary, MAC=t boundary, and the 16x12 product ordering.
- **F3 engine + CLI.** `SET_TRANSFORM` through `mainRegionInCPP` (signature
  change: regenerate RcppExports, align `Main.hpp`, delete `src/*.o`),
  `transformWeightMat` export, descriptor kind, composition/naming, the three
  CLI flags in `extdata/step2_SPAtests.R` and `SPAGMMATtest`, the memory
  print, the `--condition` without `--weights_for_condition` hard error,
  restart-schema guard already keyed on K. Gate: standard byte, every
  existing golden byte-identical.
  **Decisions taken 2026-09-09**, all of them narrowing the CLI rather than
  adding machinery (rationale in `docs/design/flexrv-design.md`):
  - **One cell per FlexRV run.** `--flexRV_maxMAF` (default 0.001) replaces
    `--maxMAF_in_groupTest` for the run, and exactly ONE annotation group is
    allowed. Both dimensions must be pinned or the region's pooled Cauchy row is
    a super-omnibus over anno x maxMAF x set rather than `p_FlexRV`. Internally
    this is `maxMAFVec = c(M)`, `q_maf = 1` -- no new machinery, and it retires
    the deferred "per-cell Cauchy row" question. Passing
    `--maxMAF_in_groupTest` as well, with a different value, is an error naming
    both flags: silently overriding a knob is the LEDGER #30 shape.
  - **M != 0.001 warns loudly and renames the sets.** The `mafbeta*` transforms
    rescale to `(1 - MAF/M)^(b-1)`, because the paper defines them relative to
    the inclusion boundary ("the most common variants included in this test will
    have the smallest scores"). The MAC ladder does NOT rescale: `MAC <= 1,3,5`
    are absolute and interpretable, and rescaling them would invent something
    the paper does not have. The warning says what actually changes: at
    M = 0.01 every sub-0.1% variant lands in (0.9, 1] for b=2, so the beta sets
    stop discriminating among exactly the variants FlexRV is for; and the MAC
    ladder then spans only the bottom sliver of the tested range (at N = 500k,
    M = 0.01 admits MAC up to 10,000). M is encoded in the set name
    (`PAI3D_sbeta5_mafbeta2@0.005`) so the OUTPUT is self-describing -- a stdout
    warning does not survive into a golden, a column label does.
  - **A MAC rung at or above `2*N*M` is dropped before the sets are composed
    (LEDGER #142, decided 2026-09-13).** Ingest admits a marker only if its
    folded MAF is at or below M, so every admitted MAC is at most `2*N*M`, and
    a rung at or above that bound gives every marker weight 1: it is `unw`
    under another name, evaluated again and counted again in the pooled row.
    At the paper's N (~424k, M = 0.001, bound ~847) no default rung is
    degenerate and the table is the full 16 x 12; at N = 50,000 the rungs
    100/200/500 go (144 sets), at N = 1,000 with M = 0.01 the rungs 25..500
    (112 sets). Measured on the sim50k study before deciding
    (`harness/sim/flexrv_study/pool_rules.R`): the null ECDF and the power are
    unchanged under every candidate rule and the per-gene `p_FlexRV` moves by
    up to 5x in either direction, so the choice was one of DEFINITION; this
    rule equals the shipped arithmetic bit for bit wherever no rung is
    degenerate and makes `p_FlexRV` independent of N through duplication. The
    ladder's VALUES are still absolute and still not rescaled; the drop is
    computed from N and the flags before a marker is read
    (`flexrvTransformTable(..., N = )`), and the run prints a warning naming
    N, the bound and the rungs dropped. The rejected alternatives -- dropping
    value-identical rows per region, or a weighted CCT -- are on the ledger
    row and in `docs/state/handoff-2026-09-13-night-135-139-142.md`.
  - **`--flexRV_lofAnno` takes a list**, separators `,` `:` `;`
    interchangeable (one flat dimension, unlike `--annotation_in_groupTest`
    where `,` separates groups and `:`/`;` join within one). Matching is exact
    and case-sensitive like every other annotation -- SAIGE has NO annotation
    vocabulary (`RegionData$ANNO[posSNP] %in% ...`,
    `R/SAIGE_SPATest_Region.R` ~:1944), so `pLoF`, `plof`, `ptv` are all equally
    valid and none is special. A label matching nothing in the group file's ANNO
    column is a hard error naming the labels present (the #113(a) precedent);
    otherwise the mask silently zeroes and 12 of the 192 sets vanish from the
    pool with no NA and no message. `isLof_i` is per-marker, evaluated once,
    reused across cells. Graded pLoF severity (LOFTEE HC/LC) belongs in the
    score line, not in the mask, which is binary by construction.
  - **Warn when a cell's score column is constant** (all 16 score transforms
    coincide; 192 sets collapse to 12) and when it takes only the values 0 and
    1 (16 collapse to <= 3). A property of the input, so no code path is made
    data-dependent.
  - **Write every set's row by default**, `--flexRV_output_sets=FALSE` to emit
    only the pooled row. The CCT is already computed on the fly from the
    in-memory p-vector, so this is verbosity, not computation order. Measured
    cost: 70 B/row in a burden-only golden, so 192 sets x 18,790 genes is
    ~325 MB per trait and ~10 s of writer time (2-3 us/row). Every analysis
    worth doing twice -- which transform won, the effective threshold, the
    score-vs-MAF optimum, why a `p_FlexRV` surprised you -- needs those rows.
    Do NOT make row emission conditional on a p-value threshold: it puts a
    statistics-dependent branch in both writers (#51).
- **F4 fixture, goldens, invariants** (full tier). Fixture
  `extdata/input/group_new_snpid_flexrv.txt` (2 genes, literal deterministic
  scores, LoF rows scored 1). Goldens (next free P numbers in
  `harness/make_goldens.sh`, plus PGEN mirrors, proved byte-reproducible twice):
  BURDEN full default 192 sets on one cell; SKAT-O on a reduced spec (2x2 sets)
  so the R writer is covered without 192 Davies calls. Invariants extend
  `tests/testthat/test-invariants-goldens.R` (the Weight-column block, ~:130) to
  the new files. Reconciliations, all
  numeric at ~1e-12, never byte across paths:
  - `unw_unw` rows == `No_Weights` rows of a plain run (null control for the
    machinery).
  - Internal expansion == external route: a checked-in weight file carrying
    the same score-only transforms (MAF `unw`, so no in-sample dependence) run
    as custom sets gives the same rows.
  - `lof_unw` rows == an unweighted `--annotation_in_groupTest=lof` run at the
    same max MAF (quantitative model, so ER gating cannot differ).
  - MAC threshold above the fixture's largest MAC == `unw`; below the
    smallest MAC: row omitted, Cauchy pool unaffected.
  - Beta inert: `--weights.beta="1,1"` vs default gives identical FlexRV rows.
  - Conditional: transform set + `--weights_for_condition` reconciles the
    way `test-step2-multiweight-cond.R` does for custom sets; without the
    flag the exact error message is asserted.
  **The fixture arithmetic is a constraint, not a detail (2026-09-09;
  CORRECTED later the same day).** The fam file of `genotype_100markers` has
  10,000 rows, but step 2 subsets to the MODEL's samples and every bundled
  model holds 1,000, so the in-sample 2N is 2,000 -- an earlier version of this
  paragraph said 20,000 and drew the cells one decade too high. Measured with
  `plink2 --freq counts --keep <the 1,000 IDs>`: the 0.001 cell (MAC <= 2)
  holds no real marker at the default collapse cutoff (every MAC <= 2 marker
  collapses), and the 0.01 cell (MAC <= 20) holds 21 polymorphic markers in
  GENE1 (MAC 1..20) and 18 in GENE2 (MAC 1..13), of which rs4/rs17/rs31/rs33
  and rs60/rs70/rs91 survive the cutoff as real markers. So at 0.01:
  `MAC25`..`MAC500` would be all-ones (five rungs that duplicate `unw`, which
  with `unw` itself makes six of twelve transforms coincide) and are dropped
  before composition since LEDGER #142 (2026-09-13), `MAC1/3/5/10`
  select only collapsed markers at cutoff 10 and are fully live at cutoff 0,
  and `mafbeta@0.01` is live on every marker. A FlexRV golden on the paper's
  0.001 would gate nothing (the "a golden that reproduces P1 gates nothing"
  trap), which is why F4 pins `--flexRV_maxMAF=0.01` (P58-P61) and asserts how
  many of the 192 weight columns are distinct on the fixture.
- **F5 N=50k stage** (SAIGE-sim worktree, fixture outside git). Derive a
  score line from `variants.tsv.gz` (`a_lat` is not written; use the anno
  class plus a seeded deterministic jitter, or add `a_lat` to the table and
  rebuild the fixture, 65 min). Run y_null -> KS on `p_FlexRV` (negative
  control, expect uniform), causal-gene ranking (positive control), and wall
  time and peak RSS at K=192 on one cell and on 12 cells; record in
  `harness/bench/NOTES.md`. **First, before any p-value** (2026-09-09): the
  fraction of individuals carrying >= 2 qualifying ultra-rare variants per gene,
  straight from the genotypes -- it bounds a priori how far the `max` collapse
  can diverge from the paper's additive pseudo-variants, and if it is a fraction
  of a percent the `sum_geno` question closes without a single test being run.
  Then `p_FlexRV` at `--MACCutoff_to_CollapseUltraRare` 0 vs 10, null and signal
  arms; check the null by the ECDF at small alpha with a binomial interval, not
  by lambda_GC. Timing must separate the ultra-rare pass, which is serial per
  cell and now runs 192 times instead of once -- the suspected new hot spot, and
  the input to the deferred duplicate-column memoisation decision. Decide from the numbers whether sparse UR storage
  or the BURDEN-loop threading (`docs/state/region-perf-handoff.md`) becomes a
  next.md item.
- **F6 docs.** CLAUDE.md paragraph (feature, flags, recommended run, the
  collapse and ER differences), roxygen, a dated handoff, next.md pointer.

## Verification

- Fast tier after F2 (`bash harness/run_tests.sh fast`): kernel closed forms.
- Standard tier byte mode after F3 and F4; `bash harness/prepush.sh` before
  any push (both compare modes). Predicted golden movement: none until F4 adds
  new goldens; any movement of an existing golden is a bug.
- Full tier at F4 (golden regeneration, double-run reproducibility).
- F5 controls reported together: null KS, positive ranking, timing; a
  misbehaving control halts, not interprets.
- TSAN not needed (region loop stays single-threaded).

## Critical files
- `R/SAIGE_SPATest_Region_Func.R` (parser), `R/SAIGE_SPATest_Region.R`
  (composition ~:862, `region_postprocess` ~:218, `SAIGE.getRegionList_new`
  ~:1700),
  `R/WeightSets.R`, `R/SAIGE_Test_main.R`, `extdata/step2_SPAtests.R`
- `src/Main.cpp` (set tables ~:1956, MAF/MAC fill ~:2463, w0 fill ~:2331,
  collapse ~:2577, UR emission ~:2891, BURDEN writer ~:3551),
  new `src/FlexRV.cpp`, `src/Main.hpp`, `src/Makevars`
- `harness/make_goldens.sh`, `tests/testthat/test-invariants-goldens.R`,
  `tests/testthat/helper-multiweight.R`, `harness/sim/02_build_dataset.py`
