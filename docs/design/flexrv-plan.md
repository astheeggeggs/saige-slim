# FlexRV on the `strip` branch

<!-- Plan written 2026-09-06/07, NOT implemented. Status: awaiting the owner's
     go-ahead; step F0 needs the FlexRV PDF (bioRxiv blocks fetches). The
     multi-weight engine it builds on is Phase 5 (C1-C11), landed. -->

## Context

**Is FlexRV implemented? No, not as a feature.** What the repo calls "flexRV"
(`docs/design/phase5-multiweights-design.md`, `docs/state/c6-gate.md`) is the
Phase 5 multi-weight region engine, fully landed (C1-C11, commits
`0fc41971`..`94a95dd7`): K weight sets per region from `weight:NAME` group-file
lines plus `--weights.beta` pairs and `No_Weights`, rows per
(anno x max_MAF x weight set), per-row ACAT-V/ACAT-O, and ONE pooled Cauchy row
per region whose each p-column is ACAT over every row (both writers:
`src/Main.cpp:2888-2909` for BURDEN, `R/SAIGE_SPATest_Region.R:606-651` for
SKAT-O). That engine is the substrate; FlexRV itself is absent.

**What FlexRV is** (Schwartzentruber, Fiziev, McRae, Ulirsch, Farh; bioRxiv
2025.11.04.686320, Illumina; STAAR implementation, code promised at
github.com/Illumina/flexrv which does not exist yet). Per gene, 192 fixed
weight vectors = 16 pathogenicity transforms x 12 MAF transforms, a burden test
per vector (SKAT optional), Cauchy-combined into `p_FlexRV`. No training.
Per-variant score `s` in [0,1]: LoF = 1, missense = PrimateAI-3D percentile
within the gene's possible missense variants, combined by max. Variants with
MAF < 0.1%. Transforms as quoted from the Methods (via a lossy fetch; **step F0
re-pins them from the PDF**):

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
that only C++ knows after QC (`src/Main.cpp:2012-2016`), and group files grow
192x. Decision (user): build the sets inside SAIGE from one `score:NAME` line,
and include the N=50k measurement stage.

## Design

**Group file.** One optional line per gene, `GENE score:NAME v1 v2 ...`
(same `[,:]` split as weight lines, `R/SAIGE_SPATest_Region_Func.R:83-92`),
after the var/anno lines and before/among weight lines; values numeric in
[0,1], NA rejected. Carried like weight columns through
`SAIGE.getRegionList_new` (sanitized `SAIGESC_1`, cf. `SAIGEWT_k` at
`R/SAIGE_SPATest_Region.R:1579-1584`) and subset with the annotation filter.
Region 1 defines the layout; every region must carry the same score name
(the same global-stride model as weights, `:129-165`).

**CLI.** `--flexRV_score=NAME` enables FlexRV sets from that score line
(error if the file has no such line). `--flexRV_scoreTransforms` and
`--flexRV_MAFTransforms` are spec strings with the paper's defaults
(`unw;lof;cbrt:0.2,...,0.9;sbeta:1.5,2,3,5,10,20` and
`unw;mafbeta:2,5;mac:1,3,5,10,25,50,100,200,500`); `--flexRV_lofAnno=lof`
names the annotation label(s) that the LoF mask keys on. Parsed in
`R/SAIGE_Test_main.R` next to `parse_weights_beta` (`:11-50`), into a
`Kt x 4` transform table (scoreKind, scoreParam, mafKind, mafParam).

**Set composition** (`R/SAIGE_SPATest_Region.R:792-931`): order becomes
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
`SET_TRANSFORM` kind in the set tables (`:1805-1826`). At the `w0_vec(r)` fill
(`:2036-2044`) and every custom-branch site (collapse `:2238-2267`, UR emission
`:2397-2409`, conditional raw export) read from `transformW(i, c)` instead of
`t_weight(i, c)` through one accessor, so the two kinds share code paths.
`transformW(i, c) = flexrvScoreWeight(kind, param, s_i, isLof_i) *
flexrvMafWeight(kind, param, MAF_i, MAC_i)` evaluated right after MAF/MAC are
final (`:2016`); `isLof_i` comes from the marker's anno label. Kernels live in
a new `src/FlexRV.cpp` with an `[[Rcpp::export]]` wrapper for tests (add to
`OBJECTS` in `src/Makevars`, `compileAttributes`, align `Main.hpp`).
`erAllowedVec[r]` is true only for the `unw_unw` set (weights identically 1),
false otherwise, matching the custom rule at `:1810`. At `Kt == 0` nothing in
the existing arms changes: every branch is on kind, the K=1 legacy arm is
untouched, so all committed goldens stay byte-identical (predicted movers:
none).

**Beta-weight independence (owner's concern, verified in code).** The
default `beta(1,25)` MAF weighting must NOT touch FlexRV sets: they replace
it. Today's custom path already guarantees this on the unconditional path
(`w0_vec(r) = t_weight(i, col)` at `src/Main.cpp:2036-2044`, R AnnoWeights
= file weights + ones for UR rows at `R/SAIGE_SPATest_Region.R:216-226`, and
the implicit beta set is composed only when no other set exists,
`:917-927`). Transform sets inherit that by construction (same accessor, same
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
`--weights_for_condition` is absent (`src/Main.cpp:1852-1860`), one vector
shared by every set (design D7). Conditioning markers live outside the group
file and have no score, so a transform cannot be applied to them. Rule:
`--flexRV_score` + `--condition` is a hard error unless
`--weights_for_condition` is supplied (message names both flags); with it,
the D7 raw-export path (C11) serves transform sets like custom sets, gated by
`test-step2-multiweight-cond.R`'s reconciliation extended to one transform
set. Documented in F0 and CLAUDE.md as a known asymmetry, not silently
defaulted.

**Zero-weight cells** (already characterised): BURDEN omits the row and the
Cauchy pool skips it (`:2807`, `:4004`); SKAT-O reports NA via the SKAT
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
- **F1 score line.** `checkGroupFile` (`R/SAIGE_SPATest_Region_Func.R:76-239`)
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
- **F4 fixture, goldens, invariants** (full tier). Fixture
  `extdata/input/group_new_snpid_flexrv.txt` (2 genes, literal deterministic
  scores, LoF rows scored 1). Goldens (next free P numbers in
  `harness/make_goldens.sh`, plus PGEN mirrors, proved byte-reproducible twice):
  BURDEN full default 192 sets on one cell; SKAT-O on a reduced spec (2x2 sets)
  so the R writer is covered without 192 Davies calls. Invariants extend
  `test-invariants-goldens.R:66-89` to the new files. Reconciliations, all
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
- **F5 N=50k stage** (SAIGE-sim worktree, fixture outside git). Derive a
  score line from `variants.tsv.gz` (`a_lat` is not written; use the anno
  class plus a seeded deterministic jitter, or add `a_lat` to the table and
  rebuild the fixture, 65 min). Run y_null -> KS on `p_FlexRV` (negative
  control, expect uniform), causal-gene ranking (positive control), and wall
  time and peak RSS at K=192 on one cell and on 12 cells; record in
  `harness/bench/NOTES.md`. Decide from the numbers whether sparse UR storage
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
  (composition :792-931, `region_postprocess` :216-226, getRegionList),
  `R/WeightSets.R`, `R/SAIGE_Test_main.R`, `extdata/step2_SPAtests.R`
- `src/Main.cpp` (:1760-1827, :2012-2044, :2238-2267, :2383-2439, :2975),
  new `src/FlexRV.cpp`, `src/Main.hpp`, `src/Makevars`
- `harness/make_goldens.sh`, `tests/testthat/test-invariants-goldens.R`,
  `tests/testthat/helper-multiweight.R`, `harness/sim/02_build_dataset.py`
