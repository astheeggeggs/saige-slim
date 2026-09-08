# The BINARY SKAT-O path at N = 50,000: where #104/#105 live (2026-09-08)

Follow-on to `docs/state/null50k-2026-09-08.md`, which covered the quantitative
arm and could not reach #104/#105 at all. Predictions pre-registered in
`harness/work/null50k-20260908/PREDICTIONS-binary.md` before any arm was read;
raw tables in the same directory.

## The binary null had to be constructed, and the construction is checkable

The fixture has none: `y_binary` carries 21 causal genes, and **permuting is ruled
out by LEDGER #104's own provenance** -- a permuted phenotype drives tau2 to 0,
which switches off the sparse-GRM machinery under test and makes the control
vacuous. So

    y_null_bin = 1[ y_null > quantile(y_null, 1 - K) ],  K = 0.0100 from y_binary

-- the same pedigree-transmitted liability background, the same 500 cases in
50,000, zero chr20 heritability. The fixture directory is untouched; the
phenotype file and the step-1 fit live in the scratch dir.

**The check that construction had to pass.** Step 1 fits `theta = (1, 0.6810)`
against `y_binary`'s `(1, 0.8464)`. tau2 is NOT zero, so the sparse GRM is
genuinely active; and it is smaller than `y_binary`'s, which is the right
ordering because `y_null_bin` lacks the chr20 gene effects `y_binary`'s null
model does not absorb. The variance-ratio file carries `sparse` rows.

## Four arms, verified at source level, not assumed

| arm | commit | #103 `sum(Phi)` | #104 `q_adj`/`var_b` | #106 `generalMatrix` |
|---|---|---|---|---|
| `bin-pre103` | `cca6b497` | no | no | no |
| `bin-lib103` | `harness/work/lib-old` | YES | no | no |
| `bin-pre104` | `4ebc7c38` | YES | no | YES |
| `bin-HEAD` | HEAD | YES | YES | YES |

Read out of each installed library with `lazyLoad` + `deparse`, because "which
build did that arm actually load" is the failure mode a version string cannot
catch.

## Every falsifiable prediction held

* **#103 alone: ZERO p-values move.** All three region columns on all 3,128
  cells, and all 9,049 `.singleAssoc.txt` rows, byte-identical.
  `BETA_Burden`/`SE_Burden` move on exactly the 1,644 multi-marker cells and on
  none of the 1,484 single-marker ones. `kappa_w` median 0.99928, max **4.791**,
  and `1/4.791 = 0.2087` is exactly the minimum `BETA_Burden` ratio observed.
  Same shape and nearly the same numbers as the quantitative arm, as expected for
  a trait-agnostic fix.
* **#106 with #104 still present: `SE_Burden` and the p-values move, `BETA_Burden`
  does not.** 1,453 cells move `SE_Burden` against **366** for `BETA_Burden`, and
  those 366 move by at most 2e-4. That is #106's STATUS claim -- "#104's rescale
  cancels a scale on `Phi` along the burden direction" -- tested where it is
  testable, and confirmed. The cancellation is not exact because #106 changes
  `Phi` element-wise rather than by a scalar.
* **#104/#105: `.singleAssoc.txt` is BYTE-IDENTICAL.** 0 of 9,049 markers moved
  on `var`, `p.value` or `p.value.NA`. The marker path never passes through
  `Phi`, and this is the measurement of that.

## #106 is ~20x smaller on binary than on quantitative, exactly as its row says

Burden-variance factor `c = (Z_old/Z_new)^2`:

| trait | median `c` | min `c` |
|---|---|---|
| quantitative (`y_null`) | 0.9988 | **0.9387** |
| binary (`y_null_bin`) | 1.000000 | **0.99685** |

#106's WHAT cell says the missing triangle is "small here only because `1/W ~ 11`
dominates `tau2 x kinship ~ 0.16`; on a quantitative trait (diagonal `tau1 ~ 1`)
... the missing triangle is first-order". Two arms at the SAME N, the SAME GRM
and the SAME genotypes differ by a factor of ~20 in exactly that direction. The
row's own reasoning is now measured rather than argued.

## #104/#105 is the dominant binary fix, and it is smaller at 50k than the goldens showed

Isolated (`4ebc7c38` -> HEAD), on 3,128 cells:

* `Pvalue` moved on 2,615 cells and got **SMALLER on 99.9 %**; `Pvalue_Burden`
  2,717 cells, **100.0 % smaller**; `Pvalue_SKAT` 2,699 cells, 99.9 % smaller.
* `BETA_Burden` **rises**, median +1.29 %, max +2.82 %; `SE_Burden` rises,
  median +0.64 %, max +1.40 %. That is the correct direction and not a typo: with
  `B = 1'Phi 1`, `BETA = S/B` and `SE = 1/sqrt(B)`, a FALLING `B` raises both and
  raises `Z = S/sqrt(B)`, so p falls. The three numbers reconcile exactly --
  `1/c = 1.0124` against a measured `BETA` ratio of 1.0129/1.0123, and
  `1/sqrt(c) = 1.0062` against a measured `SE` ratio of 1.0064/1.0060.
* `c` = `r`, the clamp's own quantity: median **0.9878**, min 0.868.

**Predicted before the run: `r` should sit closer to 1 than the bundled
fixture's 0.945-0.953, because both of its drivers are weaker here** (variance
ratios 0.988 null / 1.0000000 sparse, and only 7,500 of 50,000 samples carry any
GRM off-diagonal). Measured median `r` = 0.9878. **Held.** The consequence is
that the p-value movement is a median of **0.4 %**, not the 12.6 % the bundled
goldens showed -- though the tail still reaches **12.9 %**, so the golden headline
is the maximum, not the typical case.

## Calibration: conservative, and the fixes move it slightly toward nominal

Sign held for the third time. Every lambda rises from `cca6b497` to HEAD and none
rises far: Cauchy 0.481 -> 0.521, `hml/0.01 Pvalue` 0.813 -> 0.820,
`Pvalue_SKAT` 0.886 -> 0.931, all-cells `Pvalue` 0.675 -> 0.687.

**A binary rare-variant scan at K = 0.01 is heavily DISCRETE, and that dominates
the shape.** With 500 cases, a marker at MAC 10 has 0.1 expected case carriers, so
most rare variants have ZERO and their p piles up near 1. Real markers at
MAC 5-19 read `P(p<=0.5) = 0.128` against 0.5 -- lambda 0.317 -- while
`P(p<=0.05) = 0.058` and `P(p<=0.01) = 0.0052`. The body is not a calibration
failure; it is the discreteness, and it is why a continuous KS (D = 0.193) is
uninformative here. Read the ECDF at alpha.

**One thing to flag, and it is NOT the fixes.** The 210 independent Cauchy
omnibus rows give **6 at p <= 0.01 against 2.1 expected** (one-sided binomial
p = 0.020), and the count is **6 at `cca6b497` too** -- identical, so no part of
it is attributable to #103/#104/#105/#106. For contrast, on the same genotypes
and GRM the quantitative null gives **0 of 210** and `y_binary`'s 189 non-causal
genes give **1 of 189**. Three ~200-gene realisations reading 6, 0 and 1 against
an expectation of ~2 is suggestive at best: 0.020 one-sided before any correction
for having looked at three arms, and the smallest p (7.28e-04 over 210 genes) is
unremarkable on its own. **A single 210-gene realisation cannot settle binary tail
calibration at K = 0.01; that needs replicate phenotypes**, which is a separate
job and is recorded in `next.md` rather than answered here.

## Positive control: the signal survives, and slightly strengthens

`y_binary` run at `cca6b497` and at HEAD, against `truth_genes.tsv`:

| | `sig-pre103` | `sig-HEAD` |
|---|---|---|
| causal genes in top 20 | 16 / 20 | 16 / 20 |
| lambda over 21 causal genes | 26.01 | 26.14 |
| min p causal | 7.16e-21 | 3.98e-21 |
| LADDER Spearman(NCP, -log10 p) | 0.682 (p = 6.6e-04) | 0.682 (p = 6.6e-04) |
| non-causal `P(p<=0.05)` / `P(p<=0.01)` | 0.0317 / 0.0053 | 0.0317 / 0.0053 |

The 16/20 and the 0.682 reproduce the fixture's own 2026-09-01 `RESULTS.txt`
exactly, which cross-validates this script against `06_analyse.R`. **A fix that
improved the null by destroying the signal would have shown here, and none did**
-- the causal minimum p gets slightly SMALLER, the same direction the null moved.
