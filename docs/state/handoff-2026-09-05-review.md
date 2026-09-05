# Handoff: the code-quality review of the Phase 3 sweep, executed (2026-09-05, after the upstream-compare handoff)

One session. `/thermo-nuclear-code-quality-review` over the sweep range
(`f8cbe746^..c0688cc9`: 24 source files, +636/-183) produced the findings
below; the owner then fixed the order of execution and added LEDGER #97 at
the front and the N = 50,000 strip-vs-upstream comparison as an extra gate on
the two refactors. Every commit followed the branch rules: one commit per
ledger row, the row's test red at unpatched HEAD (scratch libraries built
from `c0688cc9`, and from `2c43edf3` for the fold) and green after, fast tier
after each R-only step, `harness/build.sh` + the standard tier after the C++
step, goldens byte-identical throughout (no golden moved this session).

## The review's findings, as written before execution

1. **`#65` rescaled files it did not write** (structural, blocker). A
   `1:22` glob filtered by `file.exists` picked up any `<prefix>_chr<j>.rda`
   left by an earlier run and multiplied it by this run's s. -> LEDGER #98.
   Its second half, save-in-phenotype-units instead of save-then-rescale, is
   an owner decision (below).
2. **Dead mirrored code in the binary branch**: the quantitative branch's
   reload-scale block pasted into `if (traitType == "binary")`, unreachable.
   Symptom of the two 270-line halves with 89 differing lines. -> #98 deleted
   the block; the fold is the last commit of the session.
3. **Five guards inlined into a 1,126-line function, and a stale copy**:
   `getNeff` carried an older version of the phenotype block that got none
   of the sweep's guards. -> `preparePhenotypeFrame` (`R/PhenoPrep.R`) and
   LEDGER #100.
4. **`parsePrintedPval` duplicated the marker path's inline parse.** -> the
   review-items commit.
5. **`implicit_beta` was the same thing as `beta` after #19.** -> retired;
   `weightSetDescriptor` (`R/WeightSets.R`) is the one constructor and the
   `weight_sets = NULL` contract is gone.
6. **`scaleQuantFit` reuse by field injection.** -> explicit arguments, and
   the injection turned out to double-scale y into every LOCO entry's
   `obj.noK`: LEDGER #99.
7. Minor: `#24` check-then-convert order, a redundant reset, a residue
   comment, `ledger.py`'s duplicated listers, `#38`'s parser as a named
   function, the stale `--nThreads` help (LEDGER #94) and the README sentence
   that documented the staleness instead of fixing it.

Not flagged, on purpose: the `#79` static factoriser (house style), `#84`'s
`m_G2tilde_cond` (set once, read const).

## What landed, in order

| commit | what | gate |
|---|---|---|
| `4fb1a975` | LEDGER #97: sparse variance-ratio rows are scale-free, not rescaled | test red at HEAD (2 failures: rows scaled 1e-4, unscaled row 1.47 = 1/s^2), 10 pass after; fast 0/1162 |
| `72ad18d7` | LEDGER #98: only the files this run wrote are rescaled; dead binary block gone | red at HEAD (planted `_chr5.rda` came back rescaled), 31 pass after; low-memory arm on a two-chromosome fixture derived at run time |
| `2d8444e3` | LEDGER #94: `--nThreads` help; README note deleted; `--help` pinned by a test | text only |
| `c70c0421` | LEDGER #99: LOCO entries' `obj.noK` rebuilt from y scaled once | red at HEAD (14 failures incl. `S_a` closed form), 43 pass after; standard 0/2661 after `build.sh` |
| `5ed4be9e` | review items 4 and 7 (C++ + R), behaviour unchanged | standard 0/2661 |
| `33a470c5` | `implicit_beta` retired; `weightSetDescriptor`; NULL contract gone | standard 0/2661, every region golden byte-identical |
| `766895f3` | `preparePhenotypeFrame` extracted (fitNULLGLMM 1,126 -> 940 lines) + unit test | fast 0/1181, standard 0/2680; **50k gate below** |
| `2c43edf3` | LEDGER #100: `getNeff` uses it (a duplicated ID is refused; Nglmm unchanged at 142.2884) | red at HEAD (accepted onto 1001 rows, died in `spsolve()`), 6 pass after; standard 0/2686 |
| `994c16c5` | the binary/quantitative halves folded into one body (fitNULLGLMM ~590 lines; the file 2,682 -> 2,265) | **cmp**: all 8 model files (both traits, both layouts) byte-identical to the pre-fold build; standard 0/2719; **50k gate below** |

`ledger.py open` is the count of what is open; rows #97-#100 are FIXED.

## The N = 50,000 gate (extra, on the two refactors)

Run from `harness/sim/upstream_compare/` against a scratch install of the
commit under test (`R_LIBS_USER`, recorded in every log's `lib=` line), so
the working tree could move on while it ran: frozen-model `s2single` and
`s2group` for `y_binary` and `y_quant`, `compare_all.sh 1 1`, then own-model
step 1 for both traits.

- **After `766895f3`** (the extraction): the frozen tables carry only
  `AF>0.5` (#28), 8 strip-only markers (#29) and the 965 / 305 strip-only
  single-in-group rows (#50); every gene cell identical; no `UNPREDICTED` in a
  frozen table (the one in the summary is the own-model fast-test straddle the
  comparison handoff already records). `cmp_summary.txt` came back
  byte-identical to the committed one. Own models: the binary `.rda`
  bit-identical to `sim50k/step1_y_binary.rda`; the quantitative `sparse` rows
  0.999996 (the #97 fix, against 1.11 before) and its `null` rows within 8e-7
  of the pre-#65 fixture's.
- **After `994c16c5`** (the fold): the same result on every count. Frozen
  tables: `AF>0.5` / `none` on the single-variant tables, every gene cell
  identical, 8 and 965 / 305 strip-only rows; no `UNPREDICTED` in a frozen
  table. Own models: the binary `.rda` bit-identical; the quantitative rows
  the same digits as after `766895f3`. One artefact to know about: this
  second run's `compare_all.sh` found the OWN-model quantitative ratio file
  the first run had refitted (sparse rows 1.00 instead of 1.11) and re-tagged
  the earlier session's own-model step-2 outputs -- produced with the 1.11
  models -- against it, so the own-model `gene y_quant` section of
  `cmp_summary.txt` and three own-model plots changed (`VR-mix` 1330 -> 0,
  `small` 23 -> 1703, `UNPREDICTED` 0 -> 119). That is a mixture of two model
  generations, not a regression; the tracked files were restored with `git
  checkout` and the own-model step 2 was not re-run (out of scope). Re-running
  `MODELS=own` step 2 for `y_quant` would make the own tables consistent
  again.

## What is NOT done: the one decision left

**`#65` save-in-phenotype-units** (review finding 1, second half). Today the
run saves the standardised quantitative model, reloads it for the ratio
estimation and rescales exactly the files it wrote at the end. The
alternative saves in phenotype units at the three save sites and standardises
on the reload: it deletes `rescaleSavedQuantModels`, the `savedModelFiles`
bookkeeping and the CLAUDE.md sentence about standardised files. Its cost:
the ratio estimation then reads a `x s x 1/s` round-tripped model (ULP
level), so P4, P16, P32, P36, P37 (+ `singleAssoc`) and the quantitative
step-1 goldens (+ pgen) would need a second rounding-level re-baseline, and
the #97 test's 1e-8 tolerances would be re-checked. Not implemented; the
owner's word decides it. Also in `next.md`: the binary LOCO warm-start
divergence and the two by-trait asymmetries the fold made visible.

## Traps met, for the next session

- **A scratch library for a drill can reuse compiled objects**: `git archive
  <commit>` into a fresh tree, copy `src/*.o` from an earlier scratch tree,
  `touch` them, delete the `.o` of any `.cpp` that changed between the two
  commits, then `R CMD INSTALL -l`. The pre-fold library built in ~1 min this
  way instead of ~6.
- **`cmp` of model files is the right drill for a step-1 refactor**, not a
  tolerance: same inputs, same code path, the `.rda` bytes must match. The
  standing test that compares the two LOCO layouts is a tolerance test and
  found a pre-existing 1e-5 disagreement (next.md) that the cmp drill was
  immune to.
- **A sparse-GRM fit never does LOCO** (`LOCO <- FALSE` at
  `R/SAIGE_fitGLMM_fast.R` under `useSparseGRMtoFitNULL`), and a dense LOCO fit
  on the full 22-chromosome fixture costs ~15 s per chromosome. The
  two-chromosome subset in `helper-pheno-scale.R` (first 1,000 records of chr 1
  and 2, the `.bed` being SNP-major) fits in 3 s and is what every LOCO test
  now uses.
- **`y` in a saved model is a one-column matrix**; `X * (y - mu)` is
  non-conformable, `as.vector` first (as `ScoreTest_NULL_Model` does).
- **`run_tests.sh` from a shell left in `extdata/`** fails on
  `source harness/env.sh`; `cd` inside a compound command persists.
- **compare_all.sh rewrites tracked files** (`cmp_summary.txt`, `plots/`).
  This session they came back byte-identical; check `git status` after a gate.

## Gate numbers

`harness/prepush.sh` at `994c16c5` (working tree = the docs of this handoff):
ledger lint OK; standard tier byte mode `[ FAIL 0 | WARN 2 | SKIP 4 | PASS 2719 ]`;
numeric mode `[ FAIL 0 | WARN 2 | SKIP 4 | PASS 2894 ]`. Not pushed: the owner
did not ask for a push this session.
