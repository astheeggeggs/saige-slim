# Plan: mine the upstream issue tracker for defects live on this branch (2026-09-04)

Copied verbatim from the session-local plan file
`~/.claude/plans/can-you-take-a-dreamy-breeze.md` on 2026-09-05, because
`docs/state/next.md` cites it for the raw-evidence pointers and the "do not
re-run the sweep" decision, and plan files outside the repo are unversioned.

---

# Mining the upstream SAIGE issue tracker for bugs live on `strip`

## Context

`LEDGER.md` now has 62 entries, every one found by reading this branch's own code.
Nobody has mined the **upstream issue tracker** — 186 issues, 108 open, 321 comments —
which is a corpus of *symptoms observed on real biobank data*, at N and in
configurations the bundled N=1,000 fixtures cannot reach. Upstream reporters routinely
misattribute and maintainers routinely close threads with a workaround, so the tracker
is full of symptoms whose mechanism was never located.

Two facts make it worth mining now:

1. **`git merge-base HEAD origin/main` is `365deeff` (2026-06-25) and
   `git rev-list --count $MB..origin/main` is 0.** This branch is level with upstream
   HEAD and *nothing has been fixed upstream since the fork*. Every unfixed upstream
   issue is live here, minus what left with BGEN/VCF/SAV, step 3, RareEffect and survival.
2. Several of these bugs the strip already removed, and knowing which is as valuable as
   knowing which survived — it says what the refactor bought.

The deliverable is a ranked collection of candidate defects, each with a mechanism in
code terms, a named output column or file it corrupts, the flags that reach it, why no
golden catches it, and the cheapest experiment that settles it.

## Method

All 186 issues pulled with `gh` (body + comments) to `<scratch>/issues/NNNN.md`, then:

- **14 triage agents**, one per batch of ~14 issues, reading every thread and grepping
  this branch for the error strings and function names it names → 73 in-scope candidates.
- **~70 locate agents**, one per candidate, tracing the mechanism in source and querying
  `harness/tools/ledger.py` → **44 PRESENT**, 11 already fixed here, 7 absent, 4 already
  ledgered, 3 not-a-bug.
- Clustering into 25 distinct defects, then **75 adversarial skeptics** (three lenses per
  candidate: code-reality, consequence, prior-art), each instructed to default to
  *refuted* → **21 survived, 4 killed**.
- A completeness critic over the whole corpus, which found the single largest item.

Findings marked **[verified here]** were re-derived by hand, independently of the agents.

> **Note the tree moved mid-sweep.** Another session landed LEDGER #58–#62 today
> (`1055706b`). The prior-art skeptics ran against that tree, so the surviving list is
> already checked against 62 entries — but **line numbers below are a mix of pre- and
> post-commit and must be re-resolved against HEAD before use.**

## Findings

### Tier 0 — RETRACTED

**N1. `conditional-spa-scale-mismatch` (#153, #37) — the stated mechanism is FALSE.**

The claim (mine and the completeness critic's, independently) was that
`q_c = t_Tstat_c/sqrt(t_varT_c/t_var2) + m1` ([src/SAIGE_test.cpp:1015](src/SAIGE_test.cpp#L1015))
divides a projected numerator by an unprojected denominator, inflating `q_c` by
`sqrt(t_var1/t_varT_c)` and diverging as the conditioning marker approaches complete LD.

**`varT_c` cancels.** `T_c/sqrt(varT_c/var2) = (T_c/sqrt(varT_c))*sqrt(var2) = z_c*sqrt(var2)`,
and the marginal line at [:784](src/SAIGE_test.cpp#L784) reduces to `z*sqrt(var2)` the same
way. The two constructions are algebraically parallel; both map a standardised statistic
onto the raw-score scale. Verified numerically: holding `z_c` fixed and sweeping
`var1/varT_c` over 1…57 leaves `q_c` constant at 56.653711. There is no inflation and
nothing diverges.

*Both of us made the same error by reading the numerator and denominator as independent
quantities without reducing the expression. Two agents agreeing is not a check.*

**What survives, measured, and demoted to Tier 3.** The conditional statistic IS
saddlepointed against the **marginal** marker's CGF: [:1036](src/SAIGE_test.cpp#L1036)
and [:1039](src/SAIGE_test.cpp#L1039) pass `t_gtilde`, `gNA`, `gNB`, `muNA`, `muNB` and
`m1` — all un-residualised — because the residual vector `g̃ − P_cond g̃` is never formed.
So the case-control-imbalance correction applied to the conditional score is the
*marginal* marker's. Measured through the exported `SPA_pval` kernel at N=20,000, 2%
prevalence, comparing the marginal-CGF call against a residualised-CGF call at matched
`z_c` (negative control at LD r²=0 gives exactly 0.000 decades, so the probe does not
fabricate signal):

| LD r² | variance kept | decades @ z_c=5 |
|---|---|---|
| 0.000 | 1.000 | **0.000** |
| 0.378 | 0.622 | −0.007 |
| 0.854 | 0.146 | −0.445 |
| 0.982 | 0.018 | −1.643 |

Anti-conservative, sub-decade until extreme LD, and *growing* with `z_c` — the wrong
shape to explain #153/#37's collapse to 1e-20…1e-103. Worth a LEDGER row as an
approximation error; **not** worth Tier 0.

**#153 and #37 therefore remain UNEXPLAINED.** Their discriminating detail is still good
evidence of something — "p becomes 0 after conditioning *when using SPA*, while the
p-value without SPA still looks normal" — and the next place to look is where `t_var2`
and `m1` are assigned on the conditional path (are they this marker's, or stale from the
marginal block?), which is the LEDGER #40 family and was not checked.

**The bundled fixture cannot reach this regime by any flag combination** — established,
and worth keeping:
- In golden P12 the conditional SPA fires on exactly 4 of 70 markers, all with
  `var/var_c` between 1.000 and 1.007.
- The two markers with real leverage are the conditioning markers themselves
  (rs79 at 57×, rs13 at 49×); conditioning a marker on itself drives `Tstat_c → 0`, so
  their `stat_c` is 0.003–0.005 against a gate of `SPAcutoff²`.
- `--SPAcutoff` is floored at 0.5 by `checkArgNumeric`, so the lowest reachable gate is
  `stat_c > 0.25`, and the best `var/var_c` among markers clearing it is **1.0085**.
- Re-running P12's command at `--SPAcutoff=2` reproduces the golden byte-for-byte
  (negative control passes).

So P12 — added *specifically* to cover this block with the LEDGER #40 fix — exercises it
only where conditioning is negligible. Any real test of the conditional SPA needs a
constructed high-LD pair, which is Phase 2 item 2.

### Tier 1 — silent wrong numbers (survived 3/3 skeptics unless noted)

**A. `step0-grm-drops-one-pair` (#160)** — **[verified here]**, 0/3 refuted.
`long int totalCombination = Ntotal2*(Ntotal2-1)/2 - 1` handed to a **half-open**
`RcppParallel::parallelFor(0, totalCombination, ...)`. I enumerated the k→(i,j) triangle
fold for N = 2,4,5,6,8,10,11,604,1000: it is an exact bijection onto all C(N,2) pairs
over `k = 0..C(N,2)-1`, so the `-1` drops exactly one **real** pair. `geno.indiceVec` is
the only pair list `refineKin` ever re-scores, so a pair never screened can never reach
the `.mtx`.

**The smoking gun is already committed** — no new run needed:

```
harness/goldens/logs/step0_sparseGRM.log:99   Ntotal: 1000
harness/goldens/logs/step0_sparseGRM.log:103  totalCombination: 499499   # C(1000,2) = 499500
```
and at N=10,000 the sim log prints `49994999` against C(10000,2) = 49995000.

*Correction to my first draft*: "two adjacent `.fam` rows" holds only for **even** N
(0-based `(N/2-1, N/2)`); for odd N = 2m+1 the dropped pair is `(m-1, N-1)`. Even N is
the common case and is the bad one — adjacent rows are exactly where family-structured
PLINK files put relatives.
*Gate gap*: on this fixture the dropped pair (500, 501) has kinship −0.0087, below the
0.05 cutoff, so **fixing it leaves the golden byte-identical**. The LEDGER #44 shape.

**B. `se-not-rolled-back-on-spa-underflow` (#94, #145, #46)** — 0/3 refuted.
The marginal SPA arm back-calculates `t_seBeta = |Beta|/t_qval` from the saddlepoint p,
then discovers that p is unusable and rolls **only the p-value** back to the score p,
never the SE it already destroyed ([src/SAIGE_test.cpp:814-838](src/SAIGE_test.cpp#L814)).
Reported upstream as `SE == 0` (#94) and "z = BETA/SE does not reproduce the p" (#145).
Distinct from LEDGER #42, which fixed the log-scale halving in the same block.
*Gate gap*: `test-invariants-goldens.R` checks the SE↔p identity only on `Is.SPA=="true"`
rows, and **no committed marker golden has `SE == 0`** (0 rows across all 16 — checked).

**C. `empty-ultra-rare-collapse-wrong-region-p` (#12, #30, #60, #67)** — 0/3 refuted.
Four independent reports converge on the binary region path's case-control rescaling when
a cell's ultra-rare collapse is empty: wrong or all-`NA` `Pvalue` / `Pvalue_Burden` /
`Pvalue_SKAT` / `BETA_Burden` / `SE_Burden`. **The four do not agree on the mechanism** —
treat "a defect here" and "this defect" as separate questions.
*Gate gap*: CLAUDE.md records that the ultra-rare collapse **always** appends a chunk on
the bundled fixtures, so an empty collapse is precisely what no golden has.

**D. `ai-reml-guards-are-scale-dependent` (#23, #47, #55, #180, #82)** — 0/3 refuted.
`--tol` is documented and used as a *relative* tolerance
([R/SAIGE_fitGLMM_fast.R:470](R/SAIGE_fitGLMM_fast.R#L470)) but also acts as an absolute
floor on tau; the phenotype-variance guard tests `Var(Y)` against 0.1; the divergence test
is `max(tau) > tol^(-2)`. All absolute — so **multiplying a quantitative phenotype by 100
changes which exit fires**: "first variance component is 0" (#23, #47), the invNormalize
abort (#55), or `solve.default(XVX): system is computationally singular` (#180).
The critic added three sharp details from #82: the guard tests `Y = eta - offset + …`,
which under the step-1 **CLI default** `isCovariateOffset=TRUE` is the *covariate
residual*, not the phenotype the message names; it is a `stop()` whose text begins
"WARNING"; and it sits inside `if (sum(tauInit[fixtau==0]) == 0)`, so `--tauInit`
silently bypasses it.

**E. `conditional-variance-cancellation-inf-se` (#24, #40)** — 0/3 refuted.
Conditional variance cancels to ~1e-17 or negative when the tested marker is collinear
with the conditioning set → `SE_c = inf`, `BETA_c = nan`, `p.value_c = 1`. LEDGER #44
added exactly this guard — an absolute `var1 <= DBL_MIN` test plus a relative cancellation
test — **to the marginal score tests only**.

**F. `duplicate-pheno-ids-deflate-grm-freqs` (#135)** — 0/3 refuted.
A repeated `--sampleIDColinphenoFile` value multiplies rows through the merge at
[R/SAIGE_fitGLMM_fast.R:1012](R/SAIGE_fitGLMM_fast.R#L1012); the inflated sample set is
what GRM marker selection and standardisation see. Step 1 *notes* the duplicates in its
log and neither dedupes nor stops. Exit 0. LEDGER #52 added duplicate-ID refusal to the
PGEN psam path only.

**G. `covariate-names-not-quoted-in-formula` (#112)** — **[verified here]**, 0/3 refuted.
The formula is pasted together without backtick quoting
([R/SAIGE_fitGLMM_fast.R:966-977](R/SAIGE_fitGLMM_fast.R#L966)); `fread` is called with
`check.names=FALSE`. `Gram+` is the loud half. The silent half is worse:
`attr(terms(as.formula("y~x1+x2+x1^2")),"term.labels")` returns `"x1" "x2"` — **a
covariate named `x1^2` is dropped with no error**, and `x1*x2` becomes an interaction.

**H. `male-x-rewrite-is-chromosome-blind` (#159)** — **[verified here]**, 1/3 refuted.
`processMale_XnonPAR` takes only the genotype vector, the position and the PAR matrix —
**no chromosome argument**. `applyMaleXnonPAR` is called for every marker gated solely on
the flag. A whole-genome run with `--is_rewrite_XnonPAR_forMales=TRUE` silently doubles
male genotypes on every autosome. The docstring warns the user instead of the code
checking. (One skeptic argued the docstring makes it documented-behaviour; the flag also
STOPs for >1 model in multimodel, so it is a live option.)

**I. `quantitative-missing-varianceratiofile-silent` (#178)** — 1/3 refuted.
`Get_Variance_Ratio` handles an absent `--varianceRatioFile` with a bare `cat()` and a
ratio of 1 ([R/readInGLMM.R:348-362](R/readInGLMM.R#L348)) — no `warning()`, no `stop()`.
For a quantitative trait r ≈ 1/Var(y|X,GRM), so omitting it scales `var`, `BETA`, `SE`
and `p.value` in **every row**. Upstream reports λ ≈ 0.4.

**J. `skato-integrate-compact-support-p1` (#209)** — 1/3 refuted.
SKAT-O's outer integral over the rho grid can return 0 on a compact support, reporting
`Pvalue` exactly 1 on a row whose per-rho p-values are small.

**K. `cate-varratio-bin-edges-disagree` (new, from the critic)** — **[verified here]**
Step-1 API default is `cateVarRatioMinMACVecExclude = c(10, 20.5)`
([R/SAIGE_fitGLMM_fast.R:710](R/SAIGE_fitGLMM_fast.R#L710)); step-2 API default is
`c(10.5, 20.5)` ([R/SAIGE_Test_main.R:131](R/SAIGE_Test_main.R#L131)). **Both CLI drivers
default to `"10,20.5"`**, so the CLI is self-consistent and only the R API disagrees with
itself across steps: a marker with 10 < MAC ≤ 10.5 is *estimated* in one category and
*applied* in another. The step-1 roxygen documents 10.5 while the code sets 10. Adjacent
to the just-landed LEDGER #61; check whether that already covers it.

### Tier 2 — crashes and silent data loss

**L. `burden-cct-stod-out-of-range` (#164)** — 0/3 refuted.
When a p-value underflows to exactly 0 the formatters emit `"%.1fE%d"` with exponent
≤ −324, a literal outside double range. The BURDEN Cauchy aggregation calls a **bare**
`std::stod` on it at [src/Main.cpp:2793](src/Main.cpp#L2793) (and :2796 for conditional).
The single-marker path does the identical conversion at [:561](src/Main.cpp#L561) and
**is** wrapped in try/catch at :559-568; the R side has a hardened `convert_str_to_log`
([R/Util.R:319](R/Util.R#L319)). This one site kept the naive parse. Quantitative traits
have **no SPA rescue arm**, so an underflowed score p goes straight to the throw.
Reproduces upstream's `Error in mainRegionInCPP(...): stod` exactly.
*Gate gap*: six goldens run `--r.corr=1` and execute the line, but no golden p-value has
a three-digit exponent, and **there is no quantitative group/BURDEN golden at all.**

**M. `conditioning-markers-have-no-qc-gate` (#88)** — 0/3 refuted.
`assign_conditionMarkers_factors` is the one `getMarkerPval` call site with no QC gate; a
dropped conditioning marker gives `arma::mean()` an empty object →
`mean(): object has no elements`, whole run aborts. **30 of the 100 bundled markers are
dropped by the ordinary marker QC**, so this is reachable on the shipped fixture today.

**N. `bim-id-type-inference-join-abort` (#137)** — **[verified here]**, 0/3 refuted.
[R/Geno.R:92](R/Geno.R#L92) reads the bim with `fread(..., select=c(1,2,4,5,6))` and no
`colClasses`, so an all-numeric ID column infers **integer** while group-file and
`--condition` IDs are always character. I reproduced the exact upstream message:
`Incompatible join types: x.ID (integer) and i.SNP (character)`. Every bundled fixture ID
is `rs<N>`, so no golden can reach it.

**O. `varratio-cursor-na-abort` (#32, #151)** — 0/3 refuted.
`extractVarianceRatio`'s cursor is incremented in two places; past the end the read yields
NA and `if (genoInd == 0)` throws. The `.rda` has already been written but the
`.varianceRatio.txt` has not, so the model is unusable. Reached by any non-autosomal
marker; every bundled step-1 fixture is autosome-only.

**P. `chrom-mismatch-zero-byte-output` (#117, #26)** — 0/3 refuted.
`--chrom` filters `markerInfo` unconditionally with no emptiness check
([R/Geno.R:102-104](R/Geno.R#L102), :125-126), and `Check_OutputFile_Create` has already
deleted the prior output → **zero-byte file at exit 0**, indistinguishable from "no hits".
The critic tied in closed issue #26 ("use 09 not 9"): the comparison is
`markerInfo[,1] == chrom` against an fread-inferred type, so both `"09" ≠ 9` and
`"chr9" ≠ 9` silently yield zero rows.

**Q. `sexcol-dropped-before-sex-filter` (#186)** — **[verified here]**, 0/3 refuted.
[R/SAIGE_fitGLMM_fast.R:908](R/SAIGE_fitGLMM_fast.R#L908) subsets `data` to
`c(phenoCol, covarColList, sampleIDColinphenoFile)`; `sexCol` is not in that list, so the
check at :938 always fires. **`--MaleOnly`/`--FemaleOnly` cannot work at all** unless
sexCol happens to be a covariate, and the error blames the phenotype file.

**R. `sparse-grm-nnz-integer-overflow` (#188)** — 0/3 refuted.
`getsubGRM` reads the `.mtx` via `Matrix:::readMM` with no size check; past 2^31 nonzeros
R reports `cannot have negative length vector` *after* reading 39 GB. One skeptic
corrected the location: the doubling is not inside `getsubGRM` (which preserves half-
storage throughout) — re-locate before writing the entry.

**S. `qcovar-cardinality-unbounded` (#73)** — 1/3 refuted.
`as.factor` on every `--qCovarColList` column with no cardinality or type check. Upstream:
n=43,790 → 216 GB RSS → `protect(): protection stack overflow`, no output. The refuting
skeptic showed the proposed cheap test does **not** reproduce at fixture scale, so the
consequence needs re-deriving before this is written up.

### Tier 3 — diagnostics, performance, defaults

**T. `sparse-sigma-refactorised-per-solve` (#109, #119, #124, #141, #183)** — 0/3 refuted.
The direct-solve branch rebuilds Sigma and runs a full SuperLU factorisation **per
right-hand side** although Sigma is unchanged across them. Five upstream reports (spsolve
memory failure, "solution not found", "gets stuck") are the biobank-N consequence.
Caching the factorisation is numerically identical, so the byte gate is not in tension.

**U. `is-spa-true-for-er-markers` (#196, #20)** — **[verified here]**, 1/3 refuted.
`t_isSPAConverge` starts false, the only SPA call sites are :802/:806, but the **ER arm**
sets it true at [src/SAIGE_test.cpp:921](src/SAIGE_test.cpp#L921). In the committed
goldens `rs44` (AC=1) and `rs56` (AC=2) both read `Is.SPA=true` — MAC-1 markers do not go
through a saddlepoint approximation. `git blame` puts :921 at upstream `626e7d71f`, added
as an undocumented line inside a commit about spsolve/PCG, so the column's meaning
silently flipped between releases. **Diagnostic only — no statistic is affected**;
the critic is right that it should not rank beside the variance defects.

**V. `region-index-off-by-one` (#56)** — **[verified here]**, 1/3 refuted.
`plink_group_burden.txt.index` reads `nEachChunk = 0` / *"Have completed the analysis of
chunk 3"* for a **two-region** run and has no *"all chunks"* line; the marker index has
all three correct plus byte offsets. A skeptic checked all 20 group indexes: same pattern,
0 completion lines. LEDGER #25 covers only the `nEachChunk` mismatch, and
`compare_goldens.sh` excludes `*.index` from the byte gate.

**W. `is_noadjCov` CLI/API default divergence (new)** — **[verified here]**
CLI defaults TRUE ([extdata/step2_SPAtests.R:557](extdata/step2_SPAtests.R#L557)), API
defaults FALSE ([R/SAIGE_Test_main.R:167](R/SAIGE_Test_main.R#L167)). These select
*different score tests*, so the two entry points give different p-values from identical
arguments. (The critic framed this as "the CLI defaults into LEDGER #28" — that is stale,
#28 is fixed; the divergence itself is still real.)

**X. Two more from the critic, unverified.** `R/readInGLMM.R:390-398` — the
`null_noXadj` cardinality check sits in the *rows-absent* arm, so it can never fire when
the rows are present (adjacent to LEDGER #39). And `R/SAIGE_fitGLMM_fast.R:2547` reads
`varRatio_sparseGRM_vec` in an arm where it is only ever assigned in the *other* arm —
unbound at k=1, stale afterwards; reachable only via the R API.

### Refuted — do not pursue

- **`plink-fam-map-sibling-of-ledger-52` (#35, #43, #103)** — killed 3/3. A `.fam` is
  headerless with a fixed six-column layout, so `line_elements[1]` **is** the IID by
  definition; there is no name-lookup ambiguity for #52 to be a sibling of. The row-count
  check the claim said was missing is **LEDGER #58**, landed today. A skeptic ran the
  proposed test and got a loud named abort, not the claimed silent shift.
- **`group-output-written-without-header` (#70, #81, #85)** — killed 3/3. This is
  **LEDGER #59 verbatim, fixed today** in `6e64b1e7`, and gated by
  `tests/testthat/test-step2-file-integrity.R:88`. A skeptic re-ran that file: 20 pass.
- **`varratio-mac-bin-staging` (#11)** — killed 2/3 **empirically**: the staging predicate
  reads only MAC bounds, and a probe at `--minMAFforGRM` 0.01/0/0.2 moved the GRM pool
  (0/40/0 markers) while the VR pool was identical in all three.
- **`cli-vs-api-covariate-offset-default` (#64)** — killed 2/3. Each entry point documents
  its own default; step 2 discards `modglmm$coefficients` outright
  ([R/readInGLMM.R:55](R/readInGLMM.R#L55), :202), so nothing reaches an output. Subsumed
  by W as a *class* of divergences rather than a defect in this flag.

### Already fixed — record, don't act

- **#197 / #206** (variance inflated for AF > 0.5 under the default fast score test) =
  LEDGER #28, fixed in Phase 6.
- **#202 / #212** (SE overestimation from log-p mishandling) = LEDGER #42, fixed (S11).
- **#191** (`ratioVec_null_noXadj` not found on legacy files) = LEDGER #13, fixed — but
  see X for a second, still-live defect in that function.
- **#201** (malformed output rows under `--nThreads=5`) — structurally eliminated: S7
  retired the forked shard-merge writer upstream still uses. A skeptic checked the one
  residual concern (`R/SkatoParallel.R:56` forking with unflushed C++ `ofstream`s) and
  showed fork children terminate without flushing inherited buffers. **Eliminated, not
  moved.**
- **#182** (`extractNglmm` dies: `Covariate_Transform` missing `traitType`) —
  **accidentally fixed**. Upstream forces the promise at `if(traitType == "survival")`;
  Phase 3a deleted survival, so `traitType` is now an unused formal. But
  [R/SAIGE_extractNeff.R:459](R/SAIGE_extractNeff.R#L459) still calls it without the
  argument — it survives on R's lazy evaluation alone and **re-arms the moment anyone
  uses that formal.** Worth a one-line fix or a LEDGER row.

## Plan

**Phase 0 — DONE for the free ones.** Nothing here needs a new run:
- **A** is settled from a committed artifact: `harness/goldens/logs/step0_sparseGRM.log:103`
  reads `499499` against C(1000,2) = 499500.
- **N1 is retracted** (see Tier 0). The surviving CGF-mismatch is measured and sub-decade.
- **G, N, U, V** are settled by the probes recorded against them.

The first thing to do now is **Phase 1**, not more probing. The one open question left
over from N1 is narrow: check whether `t_var2` and `m1` on the conditional path are this
marker's or stale from the marginal block (the LEDGER #40 family) — that is the remaining
candidate explanation for #153/#37.

**Phase 1 — write LEDGER entries, do not fix.** This branch's rule is bug-for-bug
preservation with fixes in a dedicated phase. Each survivor gets a row in the established
WHERE / WHAT / PROVENANCE / STATUS format with the upstream issue number in PROVENANCE.
**Re-resolve every line number against HEAD first** — the sweep spans commits `72051f5a`
to `1055706b` and several citations are 15–40 lines stale.

**Phase 2 — characterization tests before any fix**, asserting what the code does *today*.
Order by **gate gap, not severity** — the paths with no golden at all are where the next
#53 is hiding:
1. a **quantitative group/BURDEN golden** — currently *zero* coverage, reaches L;
2. a conditioning marker in high LD with a tested marker — reaches E, and is the
   only way to exercise the conditional SPA at non-negligible conditioning at all;
3. a group cell with an empty ultra-rare collapse — reaches C;
4. a non-autosomal step-1 fixture, made by rewriting one `.bim` column (the `.bed` is
   indexed by row, so no regeneration) — reaches O;
5. an all-numeric-ID bim — reaches N;
6. a phenotype file with a duplicated ID — reaches F.

**Phase 3 — fix, one commit per ledger entry, goldens regenerated.** Expect A, C, F
and H to move goldens. Expect B, L, M, N, O, P, Q, R to leave every golden byte-identical
— which per CLAUDE.md is a statement about coverage, not a pass.

## Verification

- `bash harness/run_tests.sh standard` green before and after every step.
- `bash harness/make_goldens.sh harness/work/candidate` +
  `bash harness/compare_goldens.sh harness/goldens harness/work/candidate` for anything
  touching step 0/1/2 output.
- New tests must pass in **both** compare modes
  (`SAIGE_COMPARE_MODE=numeric bash harness/run_tests.sh standard`) before pushing.
- Any header (`.hpp`) change: `rm src/*.o` before `R CMD INSTALL`.
- **Each new characterization test must be shown to fail at unpatched HEAD.** A test that
  passes both before and after gates nothing — which is exactly how the two refuted
  candidates above got as far as they did.

## Files this touches

`LEDGER.md`, `docs/state/next.md`, `tests/testthat/` (new characterization files),
`harness/make_goldens.sh` (new goldens), and — Phase 3 only —
`src/SAIGE_test.cpp` (N1, B, E, M, U), `src/SAIGE_fitGLMM_fast.cpp` (A, T),
`src/Main.cpp` (L, H), `R/SAIGE_fitGLMM_fast.R` (D, F, G, O, Q, R, X),
`R/SAIGE_SPATest_Region.R` (C, V), `R/Geno.R` (N, P), `R/readInGLMM.R` (I, X),
`R/SAIGE_Test_main.R` (K, W).

## Where the raw evidence lives (for a session that starts cold)

- **Issue corpus**: 186 markdown files (body + all comments), one per issue. Regenerate
  in ~30 s — the session scratchpad does not survive:
  ```sh
  gh issue list --repo saigegit/SAIGE --state all --limit 500 \
     --json number,title,state,createdAt,body,labels,author,url > issues_full.json
  gh api "repos/saigegit/SAIGE/issues/comments?per_page=100" --paginate > issue_comments.json
  # then join on issue_url and write one file per issue
  ```
- **Agent transcripts and per-agent return values** persist under
  `~/.claude/projects/-Users-dpalmer-Repositories-SAIGE/f5bd846e-4ac8-4ad4-99f6-925f4edff49e/subagents/workflows/`:
  - `wf_652dab26-b83` — triage + locate (the 44 PRESENT write-ups, with full mechanism,
    gate_gap and cheapest_test per issue). `journal.jsonl` has one `{"type":"result"}`
    line per agent; the locate records are the ones with an `issue` key.
  - `wf_7e8abe08-d79` — the 75 adversarial skeptics and the completeness critic.
- **Do not re-run the sweep.** It cost ~160 agents and ~13M subagent tokens. Everything
  worth keeping is in this file; the journals are only for chasing a specific citation.

## One standing lesson from this work

The retracted Tier 0 finding was produced independently by two agents that agreed with
each other, and it survived a three-lens adversarial pass, because every one of them read
`T_c/sqrt(varT_c/var2)` as a ratio of independent quantities instead of reducing it.
**Agreement between agents is not verification; reducing the expression is.** Before any
finding in this file is written into `LEDGER.md`, do the arithmetic by hand.
