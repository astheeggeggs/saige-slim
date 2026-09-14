# Handoff 2026-09-14: the FlexRV power fixture p2, and the answer it gives

The 2026-09-13 sweep left one open scientific question
(`handoff-2026-09-13-flexrv-sweep-rerun.md`, section 5): FlexRV beat the
unweighted burden by the same margin with the score permuted as with it
informative. Two readings were possible -- the gain is the multi-set Cauchy
combination, or the old fixture was too weak to separate the arms -- and 6e
set the design constraint for a fixture that could tell them apart: causal
signal where a burden test can resolve it. That fixture was built today,
its outcome predicted before the sweep (`prediction-flexrv-p2-2026-09-14.md`),
and the sweep run. **Answer: the second reading. On a fixture with resolvable
signal the gain does not survive permuting the score.** Owner's decisions
taken today, both recorded: the collapse cutoff stays at 10 for both traits
(upstream's default, `git show main:extdata/step2_SPAtests.R:100`), and the
fixture is this one.

## 1. What was built, and what was not touched

`harness/sim/flexrv_study/make_power_pheno.py` (docstring is the spec).
Phenotypes only: the 50k genotypes, pedigree-truth GRM, PCs, group file and
the three score group files (`harness/work/flexrv_study/scores/`) are the
2026-09-10 fixture's, byte-untouched. New files beside the fixture
(`/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k`, gitignored):
`pheno_pcs_p2.txt` (every old column plus `y_quant_p2`, `y_binary_p2`,
`y_null_bin_p2`), `truth_genes_p2.tsv` (analyse.R's columns plus per-gene
predicted NCPs), `causal_p2.tsv`, `sim_meta_p2.json`, and the three step-1
models `step1_y_{quant,binary,null_bin}_p2.*` (`step1.sh` now takes suffixed
trait names). The design in one line each:

- 60 causal genes of 210, drawn among the 179 whose visible cell
  (`high:medium`, MAF <= 0.001 -- the cell every run tests) has >= 30 minor
  alleles.
- P(causal) = 0.9 * sigmoid(10 * (within-gene percentile of a_lat - 0.5)),
  `high` at the top. a_lat is what the informative score ranks, so the score
  is informative by construction (mean cor(score, causal) 0.67); the
  realistic score (a_lat + N(0, 1.5^2)) partly; the permuted score not at all.
- One beta per gene, same sign. Per-gene NCP of the causal-only burden
  log-uniform on [10, 80], delivered exactly on the realised burden.
- Background = the fixture's `y_null` scaled by sqrt(1 - h2): same pedigree
  polygenic component the null models were built on.
- Binary at prevalence 0.10 (5,000 cases) on its own liability at
  h2 / 0.342, so NCP_binary = NCP_quant on paper. (At K = 0.01 the equal-NCP
  liability h2 is 0.55 and does not fit the budget.)
- Realised chr20 h2 0.0397; 1,182 causal variants; causal minor alleles per
  gene median 91, q10 20, min 3.

Generation controls, all passed before the sweep: ladder delivered with zero
error; sum of per-gene h2 = total var (ratio 0.999); cor(y_quant_p2, y_null)
= 0.980 = sqrt(1 - h2); bim/tsv and fam/pheno orders asserted; 0 missing
calls; the three fits converged with tau2 > 0 (binary 0.126, binary null
0.127) and carry `sparse` VR rows.

## 2. The sweep

72 runs, collapse 10 only, at `841e8d7b` clean: 4 traits (`y_null`,
`y_null_bin_p2`, `y_quant_p2`, `y_binary_p2`) x (3 score arms x 4 FlexRV
methods + unw/beta/skato/grid + the synonymous pair). Results
`harness/work/flexrv_study/results-p2/`, analysis
`harness/sim/flexrv_study/results-p2-2026-09-14.txt` (`analyse.R <results>
<fixture> _p2`). Every log ends `exit=0`; every FlexRV run composes 144 sets
(the #142 warning fires, MAC100/200/500 dropped). A first launch lost 55 runs
to `run.sh`'s load guard (the box carried an external load of 6-16 from
Chrome and Spotlight; the guard refuses before writing a log, so a refused
run leaves no trace); the relaunch with `MINLOAD_REFUSE=999` completed and the
driver's idempotence kept the 17 finished runs. Wall times are meaningless
today and are not reported.

## 3. Prediction against outcome (quantitative trait, 60 causal genes, Bonferroni 0.05/210)

| arm / method | predicted E[det] | observed | verdict |
|---|---|---|---|
| unw | 33.7 (band 30-37) | **36** | on prediction; per-gene Spearman(pred NCP, chi2) 0.83, mean ratio 1.09 |
| beta(1,25) burden | -- | 37 | |
| SKAT-O / grid (6-cell Cauchy) | -- | 37 / 40 | |
| flex192 informative | 39.0-43.5 | **47** | ABOVE the stated upper bound, see 3a |
| flex192 realistic | 37.7-40.3 | **45** | above, same reason |
| flex192 uninformative | 32.8-36.7 | **39** | inside |
| flexscore inf / uninf | -- | 47 / 39 | the score axis carries all of the gain |
| flexmaf (any arm) | -- | 39 | the MAF axis alone is +3 over unw in every arm |
| flex1 (single unweighted set) | -- | 36 | = unw, as it must |
| oracle (min over sets, Bonferroni) | <= 49.5 | 42 (inf) / 37 (uninf) | Cauchy BEATS the Bonferroni oracle on the informative arm |
| synonymous (low class only) | 0 | 0 / 0 | |
| false detections, 150 null genes | ~0 | 0 on every method | |

Paired per gene (sign test over the 60 causal genes, `y_quant_p2`):

| comparison | smaller p | sign p | median dlog10 p |
|---|---|---|---|
| informative vs unw | 47 : 13 | 1.2e-5 | +1.19 |
| **uninformative vs unw** | **29 : 31** | **0.90** | **-0.02** |
| informative vs uninformative | 55 : 5 | 1.0e-11 | +1.13 |
| realistic vs uninformative | 47 : 13 | 1.2e-5 | +0.64 |
| flexscore uninf vs unw | 28 : 32 | 0.70 | -0.08 |
| flexmaf vs unw | 23 : 37 | 0.09 | -0.21 |

**So the permuted score gives exactly the unweighted burden, gene for gene.**
The 2026-09-13 `uninf<=unw=FALSE` was the old fixture's resolution, not a
property of the method. CONTROL 3 still prints `uninf<=unw=FALSE` today
(39 vs 36) because it counts detections, and the +3 is the MAF ladder, which
`flexmaf` shows in every arm and which the design doc's control anticipated
("if it does, the gain is the MAF axis, which the MAF-only arm isolates").
Monotone holds: 47 >= 45 >= 39.

### 3a. Why the informative arm beat the prediction's upper bound

The bound was the max over my own family s^e (e in 0.2..20) of the
informative percentile score. FlexRV's `sbeta` transforms (a beta density in
the score) are sharper than any power of the rank on this link, and the
Cauchy combination of 144 correlated, individually significant sets returned
MORE detections (47) than the best single set Bonferroni-corrected (42). The
first is a wrong choice of family in the prediction; the second is the known
behaviour of CCT when many inputs are jointly small. Neither is a defect.
The bound was written as an upper bound and it was wrong by 3.5 genes; the
sign and the ordering of every arm were right.

## 4. Three things run down before this was written

**4a. Null calibration: 7 of 144 lines outside the binomial interval, all
the same line.** `y_null` (quantitative) at alpha 0.05 on the pooled Cauchy
rows: flex192 and flexmaf 0.0905-0.0952 in every arm, grid 0.0952 (which
was outside on 2026-09-13 too). At alpha 0.01, on `y_null_bin_p2`, and on the
150 non-causal genes of both signal traits, every line is inside; unw, beta,
skato, flex1 are at 0.048 on `y_null`. Side by side with 2026-09-13 (192
sets): flexmaf 0.081 -> 0.0905, flex192 0.081/0.0905 -> 0.0952/0.0905 -- the
#142 rung drop removed 48 duplicate-`unw` sets whose weight in the CCT had
been pulling the pooled row toward the calibrated `unw` row, as the row said
it would. This is the body-of-distribution excess of a Cauchy combination
over positively dependent inputs, already the open item in
`handoff-2026-09-12-126-null-calibration.md` ("Cauchy row 0.06-0.08 at
0.05"); it is 0.09 now on a quantitative null with the 144-set ladder, and
it is not FlexRV-specific (grid has it). Nothing in the tail (alpha 0.01,
Bonferroni) moves. NOT resolved; see 6.

**4b. The binary arm detects 22, not the 34 designed.** NCP_binary =
NCP_quant by construction, yet `unw` finds 22/60 against 36/60, Spearman
(predicted, observed) 0.33 against 0.83. It is monotone in the effect size:

| liability |beta| (x1.71 on the binary scale) | n | chi2 quant / binary | ratio |
|---|---|---|---|
| < 0.3 | 5 | 12.4 / 11.4 | 1.20 |
| 0.3-0.6 | 26 | 19.0 / 20.5 | 1.22 |
| 0.6-1.0 | 16 | 25.1 / 22.0 | 0.63 |
| > 1.0 | 13 | 25.8 / 13.9 | 0.22 |

Below 0.6 sd per allele the binary arm is AT the quantitative one; above 1 sd
it has a fifth of it. That is the liability-threshold linearisation
(`phi(t)^2 / K(1-K)`, the same `BIN_NCP_FRAC` the generator has used since
2026-09-01) overstating a large effect: a carrier whose liability is shifted
by 1.7 sd is a case with probability ~0.5 at K = 0.10, not 0.1 + a linear
increment, so the score test's non-centrality saturates. It is the same
defect the bench50k memory records as "equal NCP, unequal power" (there in
the other direction, at K = 0.01 with the ladder on the whole gene). It is
the FIXTURE's NCP column, not the engine: the per-gene comparisons in
section 3 hold on the binary trait too (informative vs uninformative 54 : 6,
p 1e-10; uninformative vs unw 20 : 40, i.e. the permuted score is slightly
WORSE than unw there, median -0.34). [ASSUMPTION: the saturation mechanism
is inferred from the monotone table, not from an exact recomputation; the
exact threshold-model NCP is a one-line integral and is the fix, see 6.]

**4c. Two ways to one number: the pooled row.** `p_FlexRV` recomputed from
the per-set `Pvalue_Burden` rows of every FlexRV file (10,920 pooled rows,
16 files): with the Bonferroni fallback I first wrote, max |dlog10 p| = 1.5
-- MY rule was wrong, not the file. The rule in `R/CCT_modified.R:56-69` is
ACAT's: `1/(p*pi)` for p < 1e-16 in place of `tan`, and `1/(stat*pi)` for
stat > 1e15. With it, every row agrees to 1e-3 decades except 16 rows whose
printed p is 4e-16 to 3e-15, where `1 - pcauchy(stat)` sits on the double
floor (the recomputation returns exact multiples of 2.2e-16). Reconciled.
CLAUDE.md's "`min(1, n*min(p))` Bonferroni fallback of #86" describes the
sidecar's rule label for a DIFFERENT case (all-1 or degenerate inputs), not
the small-p branch; do not read it as ACAT's small-p rule.

## 5. What this settles

- The FlexRV engine on this branch delivers the paper's mechanism: with a
  score that carries information about which variants are causal, the
  weighted-set Cauchy combination finds 47 of 60 genes where the unweighted
  burden finds 36, monotone in score quality, with zero false detections
  among 150 null genes at Bonferroni, and the permuted score buys nothing
  gene for gene. The 2026-09-13 open question is closed in the method's
  favour.
- The MAF ladder alone is worth +3 detections here but LOSES to unw gene by
  gene (23 : 37) -- it helps the genes whose causal alleles are rarest and
  hurts the rest. Not a defect, a property of the transform family.
- The collapse question does not arise on this fixture (at MAC <= 10 of
  100,000 alleles the pseudo-marker is a no-op, 6d of the sweep handoff) and
  the cutoff stays at 10 by decision.

## 6. Open, in order of cost

1. **Binary NCP column (4b).** Replace the linearised `BIN_NCP_FRAC` scaling
   in `make_power_pheno.py` (and `02_build_dataset.py`, same constant) with
   the exact threshold-model non-centrality per gene, or scale the binary
   liability so that the EXACT NCP matches the quantitative one. A
   prediction-first check then repeats 3 on `y_binary_p2`. Until then the
   binary counts in `results-p2-2026-09-14.txt` are a fixture artefact and
   the paired comparisons are the only binary statements to quote.
2. **Cauchy body excess at alpha 0.05 (4a).** 0.09 on the quantitative null
   at 144 sets, 0.095 on the 6-cell grid; the tail is clean. Whether a
   paper-style permutation on real null variants shows the same is the
   question the paper's own calibration answered by construction; on this
   branch the number to watch is the 210-gene ECDF at 0.05, and a claim
   about FlexRV at nominal alpha 0.05 should not be made from the pooled row
   without it. A dedicated null series (several `y_null` draws) would turn
   the 0.09 into an interval.
3. The 2026-09-12 decisions memo's recommendation ("stratified + class-only
   + flex24 at collapse 0") is superseded for its collapse-0 half by the
   decision and for its fixture by this one; the memo is dated and stays.
4. `harness/sim/flexrv_study/study.sh` still defaults `COLL="0 10"`; the
   collapse-0 column is no longer run by decision. Left as a default so the
   old sweep remains reproducible; pass `COLL=10`.

## 7. Reproduce

```sh
# fixture (2 s), models (~1 min), sweep (~25 min at load 1)
python3 harness/sim/flexrv_study/make_power_pheno.py --dir $D            # miniconda python
bash harness/sim/flexrv_study/step1.sh y_quant_p2 y_binary_p2 y_null_bin_p2
OUT=harness/work/flexrv_study/results-p2 TRAITS="y_null y_null_bin_p2 y_quant_p2 y_binary_p2" \
  COLL=10 MINLOAD_REFUSE=999 bash harness/sim/flexrv_study/study.sh 3
Rscript harness/sim/flexrv_study/analyse.R harness/work/flexrv_study/results-p2 $D _p2
```
`--dry-run` prints the prediction table without writing. The paired per-gene
statistics of section 3 were computed by a scratch script on the same files;
the numbers in the tables are what it printed.
