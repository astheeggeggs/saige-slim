# Handoff: LEDGER #103, #106, #104/#105 fixed (690894f1, 4ebc7c38, e338bdbd), #44 re-fixed, and the solver benchmarked at N = 50,000 (2026-09-08)

One session, following the owner's instruction to apply the three fixes the
2026-09-07 diagnosis left as decisions, in the order #103 (beta denominator),
#106 (sparse Sigma and its PCG), #104/#105 (the burden rescale), and to run
the solver at biobank scale with relatedness present, comparing time and
accuracy with the current one. One commit per ledger row, characterization
inverted in the same commit, golden movement predicted before each run and
matched, goldens re-baselined only after two `make_goldens.sh` runs were
`cmp`-identical, the pre-push gate in both compare modes before each commit.

Read `python3 harness/tools/ledger.py 103` / `106` / `104` for the rows; the
measurement that decided the order is
`docs/state/diag-2026-09-07-region-vs-burden-exact.md` (both burden
constructions against the exact model variance).

**One thing the plan did not anticipate, and it is the session's main lesson:
the #106 fix RE-OPENED LEDGER #44.** `scoreTest`'s sparse arm reported
`var == 0` EXACTLY on a covariate-collinear marker only because the OLD
absolute PCG tolerance was already satisfied at `x = 0` (`||b||^2 ~ 1e-29 <
0.02`), so the solver never took a step and #44's `var1 <= DBL_MIN` guard was
firing on an accident of the stopping rule. The converged solve returns the
true tiny variance, the guard stops firing, and P27's two rows printed
`BETA=-4.4e13 SE=1.2e14 p=0.72`. This was caught by the pre-push gate, NOT by
the golden movement prediction, which had said "P27, all 70 `var`" and was
right about that while missing what it implied for `BETA`/`SE`/`p`. See the
"Traps met" section.

## #103 -- BETA_Burden divides by 1'Phi 1

One line in `get_SKAT_pvalue`. Predicted and measured: 15 group goldens moved,
only `BETA_Burden`/`SE_Burden` (and the `_cond` twins in the three conditional
ones), multi-marker cells only; every p-value, BURDEN, marker, step-0 and
step-1 output byte-identical; 60 step-2 files identical between the two
regeneration runs. `test-step2-burden-trace.R` inverted;
`test-unit-region-cond.R`'s #33 block re-pointed at the diagonal identity it
rests on. The golden `.index` files were left as they were: they have not
been re-baselined since Phase 3 commit 0 and are excluded from the byte gate.

## #106 -- one-triangle Sigma, absolute PCG tolerance

`setSparseSigma_new` expands the subset to a general triplet matrix (both
triangles) and `getPCG1ofSigmaAndGtilde` stops on `||r||^2 <= 1e-8 ||b||^2`
(file-scope constants at the top of `src/SAIGE_test.cpp`; the constants must
sit ABOVE their first use at :227, the first build failed on that).
Characterization `test-step2-sigma-triangle.R`: green at HEAD (10 pins), red
against `harness/work/lib-old` once inverted (8 failures), green after the
golden re-baseline (its P27 pins read the golden file). Movement predicted and
matched: every binary group golden with `--sparseGRMFile` and its singleAssoc,
the PGEN/multi-weight/conditional/BURDEN mirrors, P27 (all 70 `var`), P29 (2
re-tested markers); the quantitative group goldens, other marker goldens,
step 0 and step 1 byte-identical. On the SKAT-O path SE_Burden and the
p-values move but BETA_Burden does not: #104's rescale cancels a scale on Phi
along the burden direction, as the diagnosis predicted.

## #104/#105 -- the burden rescale to the naive variance

Fix A, one change in `get_newPhi_scaleFactor_traitType`: the burden statistic
is standardised by `1'Phi 1` before the `|z| < 2` cutoff and the saddlepoint
(`var_b <- sum(Phi)`; `q_adj <- m1 + (q.sum - m1) * sqrt(var1 / var_b)`), the
device the C++ marker SPA already uses. Below the cutoff `VarQ_2` then equals
`1'Phi 1`, so `r = min(VarQ/VarQ_2, 1)` is a no-op and the whole of `Phi` stops
being inflated by `var1 / 1'Phi 1`. #105 closed with it, one mechanism.

Characterization `test-step2-phi-rescale.R` green at HEAD (pinning the defect:
adjusted burden variance `== var1`, and a scalar on `Phi` cancelling entirely),
then inverted; red against `harness/work/lib-old` on exactly the three inverted
assertions (8.93 vs 8.04, 8.9 vs 4.5), green after the rebuild. The negative
control -- `Phi` the naive covariance -- holds both ways.

Movement predicted before the run and matched: 14 step-2 goldens, all of them
the binary SKAT-O group runs and no others; `Pvalue`, `Pvalue_Burden`,
`Pvalue_SKAT` on every row, `Pvalue_ACATO` where the multi-weight schema
carries it, `BETA_Burden`/`SE_Burden` on every multi-marker cell, and the
`_cond` twins in the conditional files. `Pvalue_ACATV` does NOT move (it is the
ACAT over the single-variant p-values, which never pass through `Phi`) and no
`.singleAssoc.txt` moves. Sign and magnitude control, stated first: the burden
variance falls from `var1` to `VarQ`, so implied variance ratio 0.94148-0.99999
over P6's cells (the row's `r` = 0.945-0.953, shifted slightly by #106),
`SE_Burden` up to +3.06 % = `1/sqrt(0.94148)`, `BETA_Burden` up to +6.2 %, and
every p-value on every cell falls (max 12.6 % relative on `Pvalue_Burden`).
Byte-identical: both quantitative group goldens, every `--r.corr=1` BURDEN
golden, every marker golden, step 0, and the 10 step-1 comparisons PASS.

**Watch item, not acted on:** `var_b <- sum(Phi)` divides without a guard. `Phi`
is PSD so `1'Phi 1 >= 0`, and it is 0 only if the burden direction lies in its
null space -- no bundled cell comes near it, but a region of entirely
monomorphic markers would give `NaN` rather than a message. Left alone because
widening #104's scope unasked is not this branch's rule; recorded here so the
decision is visible.

## The N = 50,000 solver benchmark

Fixture: `/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k` (not in
git; 42,500 singletons + 1,875 families of 4, pedigree-truth sparse GRM with
59,375 entries, binary trait `y_binary` fitted with the sparse GRM). Harness
under `harness/work/bench106/` (not in git): `run_arm.sh <tag> [R_LIBS]` runs
SAIGE-GENE+ on the first 40 genes (`group40.txt`) at one thread under
`/usr/bin/time -l`, loading SAIGE from a scratch library if one is given;
`acc.R <arm> [<arm2>]` computes the exact `vr x g~' Sigma_s^-1 g~` by sparse
Cholesky for 300 sampled sparse-arm markers (MAC 11-20) and compares each
arm's printed `var` (Tstat reproduced as the control); `iters.R` runs the R
replica of the C++ PCG (old rule and new rule) on the same markers.
The old-solver arm ran from `harness/work/lib-old` (the tree at #103, before
#106, installed with `R CMD INSTALL -l`).

Old solver (before #106): 40 genes, 58.8 s wall / 45.8 s user (machine load
3-5, provisional); `var/exact - 1` median +1.4e-4, max +2.9e-3, 37 of 300
markers beyond 1e-3; one PCG iteration on every marker. The R replica
reproduces those errors to the digit, and the converged symmetric replica
needs 2-3 iterations for errors <= 1.2e-8.

New solver (after #106, main env, goldens regenerating concurrently):
40 genes, 42.6 s wall / 37.7 s user, max RSS 726 MB (old arm 722 MB);
`var/exact - 1` median -1.4e-8, mean +6.0e-8, min -4.5e-6, max +4.8e-6 -- the
printed precision of `var` -- 0 of 300 beyond 1e-3. Over all 1,631 markers in
the 40 genes, `new/old var - 1` ranges -2.9e-3 to +4.5e-5, median 0 (the
common-arm markers, MAC > 20.5, do not touch the solver and are identical).
The region output `group.txt` differs between the arms (same 638 rows), as it
must: the sparse-arm variances feed every Phi.

Timing, like-for-like (arms alternated back to back, two repeats each, load
3-4, outputs byte-identical across repeats): old 40.5 / 40.5 s wall
(35.7 / 36.1 user), new 39.9 / 40.8 s wall (35.6 / 36.3 user). No measurable
cost: 2-3 PCG iterations replace 1, and the solve is a small share of step 2.

## Traps met

- **A guard whose predicate rests on an ACCIDENT of another component's
  tolerance is not a fix -- it is a coincidence with a test around it.** #44's
  absolute `var1 <= DBL_MIN` guard passed every gate for five days because the
  PCG's absolute tolerance happened to return exactly 0. Changing the tolerance
  in #106 broke it. When predicting golden movement, a column that moves off an
  EXACT special value (0, NA, inf) deserves the question "what fired on that
  value?", not just "by how much did it move?".
- The consequence was worse than the original bug: a `4.4e13` entry raises the
  numeric comparator's column-scaled absolute floor by 14 orders and would have
  silently masked the other 68 markers of P27's `BETA` column -- a detection-
  power loss that byte mode would have caught and numeric mode would not.
- Two `expect_equal(golden, exact)` comparisons in the new sigma-triangle test
  were toleranced at 1e-6, below the 6-significant-digit precision `var` is
  WRITTEN with (`1.62837` for `1.6283745`, a 2.8e-6 gap). A tolerance tighter
  than the writer's own quantisation measures the `%g`, not the estimator.

- The `.bed` bundled with `genotype_100markers` holds 10,000 samples; the
  model's 1,000 are matched by IID (the new sigma-triangle test does this).
- plink2's `--export A` counts bim column 6; SAIGE's Allele2 under
  `--AlleleOrder=alt-first` is bim column 5, so every exported column is
  flipped before use.
- `expect_identical(x, 1L)` on a counter built by `it <- it + 1` compares a
  double with an integer; use `expect_equal`.

## RESUME HERE (kept current as the session advances)

**All three rows are done and committed.** #103 at `690894f1`, #106 at
`4ebc7c38`, #104/#105 at `e338bdbd`. LEDGER #44's
STATUS carries the third shape; #104's and #105's STATUS cells are filled and
`ledger.py lint` reports 5 open (from 7).

The pre-push gate was run in BOTH compare modes before each commit. #106's
first gate run was RED -- five failures, and both were real: the
sigma-triangle test's 1e-6 tolerance sat below the golden's printed precision,
and the marker-paths characterization had caught #44 re-opening. Both are
fixed in the #106 commit; see "Traps met".

Goldens: 43 files re-baselined for #106 (plus P27's two degenerate rows once
`projectionCollapsed` landed), 14 for #104. Each set proved byte-reproducible
by two `make_goldens.sh` runs into separate trees before being copied.

**Verified at `4e47b760` (end of day): the FULL tier and the negative control,
which had last run on 2026-08-30 -- before every fix in this handoff.**
`bash harness/run_tests.sh full` gives FAIL 0, WARN 2, **SKIP 0**, PASS 2998;
SKIP 0 is the point, since it means the golden regeneration, the step-1
determinism double-run and the permuted-phenotype control all actually ran
rather than being skipped as they are in the standard tier. The control passes
through both readers: 70 markers, min p = 2.62e-02 against a 1.4e-05 flag,
lambda_GC 0.505 / 1.078 / 0.783 over the MAC<=4 / 5-19 / >=20 strata, KS on the
32 near-continuous markers D = 0.202, p = 0.126.

All four committed control goldens came back BYTE-IDENTICAL, which was predicted
before the run and is a real confinement check rather than a formality: the
control's step 2 carries no `--sparseGRMFile`, so #103/#104/#105/#106 cannot
reach it, and its step 1 is a binary fit whose Sigma diagonal is >= 4, far above
#91's old 1e-4 floor. Any of those leaking outside its intended path would have
moved these bytes.

Read the control for what it is: 32 markers at the KS gate has power against
GROSS miscalibration and very little against subtle, and `harness/check_null.R`
says so in its own header. It also only exercises the SINGLE-VARIANT output,
while everything this handoff changed is on the REGION path. The strongest
calibration evidence for the day remains #106's N = 50,000 comparison against an
exact sparse-Cholesky solve (median relative error 1.4e-4 -> 1.4e-8), which is
ground truth for variances rather than uniformity for p-values. A permuted
phenotype through a REGION run on the `sim50k` fixture would be the honest next
step there, and is the one calibration check nobody has run.

What is NOT done, for whoever comes next:

- The `sum(Phi)` divide in `get_newPhi_scaleFactor_traitType` has no guard for
  a zero burden variance (see the #104 section's watch item). A decision, not
  a defect on any bundled input.
- `docs/state/next.md`'s open items are otherwise unchanged by this session;
  re-read `python3 harness/tools/ledger.py open` for the 5 that remain.
- The scratch material is still under `harness/work/resume-2026-09-08/`
  (untracked): the apply scripts, `setcell.py`, `coldiff.R`, the composed
  ledger cells and all three commit messages. `harness/work/lib-old` is the
  pre-#106 scratch install the drills used; `harness/work/bench106/` is the
  N = 50,000 solver benchmark harness.
