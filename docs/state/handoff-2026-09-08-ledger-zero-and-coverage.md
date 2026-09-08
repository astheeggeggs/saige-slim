# Handoff: the ledger reaches 0 open, and a coverage audit closes the quantitative sparse-GRM hole (2026-09-08, second session of the day)

Branch end commit: `1350394a`. The first session of the day is
`docs/state/handoff-2026-09-08-burden-variance-fixes.md` (#103, #106, #104/#105
and the N = 50,000 solver benchmark); read it first if you need that story, but
this file is the current pointer and carries everything after it.

`python3 harness/tools/ledger.py lint` reports **0 open of 106, 4 decided**.

## NEXT UP, and it is a named job, not a suggestion

**Run the permuted-phenotype (null) REGION control at N = 50,000 and check
p-value uniformity after the day's fixes.** This is the one calibration check
nobody has run, and it is the honest gap left by everything below.

Why it is the next thing:

* Three fixes committed today deliberately MOVED p-values -- #104/#105 by up to
  12.6 % on every binary SKAT-O cell, #106 across 43 goldens, and the new
  quantitative sparse-GRM goldens by up to 23 %. Each is justified individually
  and drilled both ways, but "the variance estimate is now closer to the exact
  model variance" is an argument, not a calibration measurement.
* The bundled negative control DOES pass at HEAD (see "What was verified"), but
  it gates a KS test on **32 markers** and only on SINGLE-VARIANT output, while
  every change today was on the REGION path. `harness/check_null.R` says so in
  its own header. It rules out gross breakage and almost nothing else.
* The sim50k fixture gives the KS test real power: the existing null region run
  has **3,339 region rows and 8,795 single-variant rows**.

**Everything needed is already on disk. Checked 2026-09-08, do not re-derive:**

* Fixture `/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k`, **3.1 GB,
  present**, gitignored. `SAIGE-sim` is a git WORKTREE sharing
  `/Users/dpalmer/Repositories/SAIGE/.git`, so `strip` is checked out there too.
  Rebuild is ~65 min if it is ever lost (`CONFIG=bench50k bash
  harness/sim/run_all.sh`), but it does not need rebuilding.
* `pheno.txt` / `pheno_pcs.txt` carry a **`y_null`** column -- a phenotype with
  no simulated effect, which is the control you want. (Permuting is not needed;
  the fixture was built with a null arm.)
* `step1_y_null.rda` + `step1_y_null.varianceRatio.txt` already exist AND the
  ratio file **already has `sparse` rows** (2 null, 2 null_noXadj, 2 sparse), so
  the region run can be given `--sparseGRMFile`. This is the exact enabler that
  was MISSING on the bundled fixture and that this session had to add there as
  `quant_sparseGRM_vr`; at 50k it is already in place, so there is no step-1
  work to do.
* `step2_group_y_null.txt` **already exists from 2026-09-01**, produced WITH
  `--sparseGRMFile=.../pedigree_truth.sparseGRM.mtx`. It therefore predates
  #103/#104/#105/#106 and is a genuine BEFORE baseline: diffing a fresh run
  against it quantifies the day's fixes at biobank scale on a null, which
  complements #106's N = 50,000 accuracy work (that was on the alternative, and
  measured variances against an exact Cholesky solve, not p-value uniformity).
* Driver: `harness/sim/04_run_saige.sh:82` builds `step2_group_$col.txt`. The
  40-gene benchmark harness `harness/work/bench106/run_arm.sh` is the other
  reusable piece, and `harness/work/lib-old` is still the pre-#106 scratch
  install if a three-way comparison is wanted.

Suggested shape, not prescribed: re-run the y_null region arm at HEAD, then (a)
KS / lambda_GC on the region `Pvalue`, `Pvalue_Burden` and `Pvalue_SKAT`
columns and on the singleAssoc p-values, stratified the way
`harness/check_null.R` stratifies (MAC <= 4 is ATOMIC on a binary trait and a
continuous KS against it is invalid -- that lesson is already paid for), and (b)
diff against the 2026-09-01 file to see what the fixes did. State the expected
sign first: the fixes REMOVED an inflation of the burden variance, so p-values
should get SMALLER and lambda_GC should move UP toward 1 -- if the null was
conservative before, it should be less so now, and it must not overshoot.

**Two fixture defects are known and do NOT block this job** -- they block the
POWER controls, not the null one: the NCP model is wrong (binary detects 16/20
causal genes, quantitative 5/20 at identical per-gene NCP) and no directional
genes exist (`|sum(b)|/sum(|b|)` spans only 0.02-0.50). `y_null` was clean on
both gene and single-variant tests on 2026-09-01 and the ladder statistic was
correctly null (0.085, p = 0.71). See the `saige-bench50k-fixture` memory.

## What was done in this session

Seven commits, `9e9f848c..1350394a`, each gated by `bash harness/prepush.sh`
green in BOTH compare modes before it landed.

**The ledger sweep, four rows, one commit each.** All four were "dead surface
left by the strip", and each was bigger than its row said.

* **#92 `9e9f848c`** -- `extdata/cmd.sh`. The dead `--bgenFile`/`--vcfFile`
  flags were only the first layer. Rewriting for PLINK/PGEN found three more
  ways the blocks were dead: a `--varianceRatioFile` two of them named
  (`...varianceRatio.col1.txt`) **does not exist in the repo at all**; no
  committed `extdata/output/` artifact can drive a single-variant test any more
  (they carry a `null` row only, and `prepareNullModel` requires
  `null_noXadj`); and the `--cateVarRatio*` overrides two blocks passed abort on
  a count check. So cmd.sh is now SELF-CONTAINED -- three step-1 blocks fit
  their own models under `./output/example_*_cmd.*` and seven step-2 blocks test
  against them. Measured, not asserted: the whole file runs end to end under
  `set -e` with every output redirected to scratch, and block 2 is block 1 with
  only the reader swapped, byte-identical output included.
* **#95 `c902edd2`** -- six help strings and four roxygen twins. One was not
  merely stale but WRONG about current code: `--is_rewrite_XnonPAR_forMales`
  still told users "the program does not check the chromosome value", which
  LEDGER #69 made false; the roxygen had been corrected then and the CLI help
  had not.
* **#96 `22770bfd`** -- the 51 bgen/vcf/sav fixtures (18 MB) MOVED, not deleted,
  to `harness/goldens_retired/inputs/`, which is what the row itself proposed:
  the retired bgen/vcf GOLDENS already live under `harness/goldens_retired/`, so
  inputs and expected outputs stay paired. `extdata/input/readme` went with them
  (every line was a recipe for a moved file). The `.log` files stayed, because
  those runs also produced the live bed/bim/fam.
* **#93 `af24ece8`** -- the man pages. The row assumed roxygen2 needed adding to
  the env; it is already there at 8.1.0. **The trap is NAMESPACE**: it is
  hand-written around `exportPattern("^[[:alpha:]]+")` while only 6 functions
  carry `@export`, so the default roclets would have silently narrowed the
  package's exports to those six. Use `roxygen2::roxygenise(".", roclets = "rd")`
  and `cmp` NAMESPACE afterwards, every time. Regeneration fixed only one
  direction, so 28 `@param` were added to the roxygen source and 8 wrong names
  deleted; all six pages now match `names(formals())` exactly and `checkRd` is
  clean. Two pre-existing defects surfaced: a stray `#'` inside a function body
  generating a garbage page, and an **unescaped `%`** in the `nThreads` `@param`
  -- `%` is Rd's comment character, so it truncated the line, swallowed its
  closing brace and broke the whole `\arguments` block.

**#91 `b09d5eb4` -- the last open row.** Sigma's diagonal was floored at an
absolute `1e-4` in phenotype units at three live sites, and **the fix differs by
consumer**, which only tracing every caller showed. `gen_sp_Sigma`'s clamped
diagonal IS the matrix -- its only caller is `gen_spsolve_v4`, the direct sparse
solve the sparse-GRM fits take at default flags -- so its floor is DELETED with
no replacement: a singular Sigma now reaches the
`"spsolve(): solution not found"` that function already threw, and the clamp was
exactly what stopped that being reached. `getDiagOfSigma` and
`getDiagOfSigma_LOCO` are only ever `1/getDiagOfSigma(...)`, the Jacobi
PRECONDITIONER, which cannot change a CONVERGED PCG's answer -- but SAIGE's PCG
stops on an absolute residual (#6), so it did; those now clamp only
non-positive/non-finite entries, to 1. They keep a fallback rather than throwing
because the AI-REML step sets a `tau` below `tol` to exactly zero
(`tauVec.elem(find(tauVec < tol)).zeros()`), so a zero diagonal is REACHABLE and
throwing would turn a limping fit into a hard error. No golden moved, predicted
beforehand from the goldens: the minimum diagonal over all five step-1
configurations is 0.216 (quantitative) and 4.0 (binary), far above 1e-4.

**The coverage audit, and `4e47b760`.** Answering "what do the goldens not
gate": parse the 38 step-2 golden commands out of `harness/make_goldens.sh`,
cross them against trait type and test type, sweep all 133 CLI options for ones
no golden touches, check those against the test suite. Full result is in
`docs/state/next.md`, section "Golden coverage audit (2026-09-08)".

The headline: **23 of 33 binary step-2 goldens passed `--sparseGRMFile` and 0 of
5 quantitative ones did**, so nothing put a quantitative model through
`setSparseSigma_new` or `getPCG1ofSigmaAndGtilde`. That is the regime where
LEDGER #106 mattered MOST by its own WHAT cell ("on a QUANTITATIVE trait ... the
missing triangle is first-order"), so #106 moved 43 binary goldens and could not
move one quantitative byte -- not because that path was right, but because
nothing looked. Measured against the pre-#106 build on the new BURDEN golden
(which isolates #106, since #104 provably does not touch `r.corr=1`):
`BETA_Burden` up to **64.5 %** out, `Pvalue_Burden` up to 11.9 %, against ~0.5 %
on the binary goldens. Closed with `quant_sparseGRM_vr` (step 1, + PGEN mirror)
and goldens P38-P41.

## What was verified at the end

`bash harness/run_tests.sh full`: **FAIL 0, WARN 2, SKIP 0, PASS 2998**. `SKIP 0`
is the point -- golden regeneration, the step-1 determinism double-run and the
permuted-phenotype control all actually executed. The control passes through
both readers (70 markers, min p 2.62e-02 against a 1.4e-05 flag, lambda_GC
0.505 / 1.078 / 0.783 over MAC <= 4 / 5-19 / >= 20, KS D = 0.202 p = 0.126 on
the 32 near-continuous markers), and all four committed control goldens came
back BYTE-IDENTICAL -- predicted beforehand, and a real confinement check, since
the control's step 2 carries no `--sparseGRMFile` and its step-1 Sigma diagonal
is >= 4.

## Traps met, each already paid for

* **A guard whose predicate rests on an ACCIDENT of another component's
  tolerance is a coincidence with a test around it.** #106 RE-OPENED LEDGER #44:
  `scoreTest`'s sparse arm reported `var == 0` exactly only because the OLD
  absolute PCG tolerance was satisfied at `x = 0`, so the solver never stepped
  and #44's `var1 <= DBL_MIN` guard fired on that. The converged solve returns
  the true tiny variance, the guard stops firing, and P27 printed
  `BETA=-4.4e13 p=0.72`. Caught by the GATE, not by the movement prediction,
  which said "P27, all 70 `var`" and was right while missing what a column
  moving off EXACTLY ZERO implied. Fixed in the same commit with
  `projectionCollapsed` (`||gtilde||^2 <= (n*eps)^2 ||g||^2`, 29.5 orders of
  separation on the fixture). When predicting movement, a column moving off an
  exact special value gets the question "what fired on that value?".
* **A new golden gates nothing until a test re-runs it.** The standard tier
  re-runs goldens from an EXPLICIT list, not by globbing. P38-P41 would have sat
  in `harness/goldens/` gating nothing; the pass count going 2949 -> 2972 (byte)
  is the evidence they now run. And a new golden must be checked to DIFFER from
  the twin it was derived from, or it gates nothing either.
* **A tolerance tighter than the printed precision measures the writer's `%g`.**
  `test-step2-sigma-triangle.R` compared a golden `var` (6 significant digits)
  against an exact solve at `1e-6`; the quantisation is ~3e-6.
* **Cloning a test to a new regime: every hard-coded NUMBER is a claim about
  that regime.** The quantitative twin of a binary invariant inherited
  `abs(sparse_ratio - 1) < 1e-3`; the quantitative ratios are ~3.2e-3 out,
  because the GRM's share of `Sigma = tau0/W + tau1*GRM` is 2.2x the diagonal
  term at W = 1 against <= 0.08 when `tau0/W >= 4` -- about 26x, the same
  ordering that made #106 first-order for quantitative traits. Re-measure, do
  not copy. And run the FAST TIER before the 9-minute gate on any new test; that
  rule is in memory and was still violated once today.
* Do NOT `rm src/*.o.tmp` while a build is running -- it deletes the file the
  compiler is about to rename and fails the build. Self-inflicted, once.

## Still open, in the order I would take them

1. The N = 50,000 null region control above.
2. **A PLINK/PGEN missingness fixture.** Verified by decoding the .bed:
   `genotype_100markers` has ZERO missing calls, and `impute_method`,
   `maxMissing`, `is_imputed_data`, `minInfo`, `dosage_zerod_cutoff` and
   `dosage_zerod_MAC_cutoff` have no golden AND no test. The only missingness
   fixtures the repo ever had were bgen/vcf, which #96 retired this morning.
   Generate one from the existing bed with plink2.
3. `LOCO=TRUE` appears in ZERO golden runs (three tests cover pieces).
4. The `sum(Phi)` divide in `get_newPhi_scaleFactor_traitType` has no zero
   guard -- a decision, not a defect on any bundled input.
5. Gate timing: `prepush.sh` is ~9 min against the ~5 CLAUDE.md implies, never
   measured idle. The open question in `next.md` is whether numeric mode needs
   to re-run every step-1 fit, since the modes differ only in how step-2 TABLES
   are compared. That is close to a halving of every future gate.
6. FlexRV: planned in `docs/design/flexrv-plan.md`, nothing built; its F5 stage
   is blocked on `a_lat` not being written by `02_build_dataset.py`.

Scratch material from both of the day's sessions is under
`harness/work/resume-2026-09-08/` (untracked): the apply scripts, `setcell.py`,
`coldiff.R`, the composed ledger cells and every commit message.
