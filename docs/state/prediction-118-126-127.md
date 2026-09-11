# Golden-movement prediction for LEDGER #118, #126, #127

Written 2026-09-11, **before any gate was run**, against the working tree that
implements the three rows (scratch build, `R_LIBS=<scratch>/newlib`). Method
copied from the #124 re-baseline's `prediction-regen.txt`: name the file, the
column and the mechanism first, then let the gate falsify it.

Handoff: `docs/state/handoff-2026-09-11-percell-mac.md`. **Its "Prediction to
beat" is wrong and is corrected below** -- movement is NOT confined to
P58/P59/P61, because option (b) changes the pseudo-marker's units on every
custom set, zero weights or not. Section 6 says why, and what that does to the
gate-5 matrix the same handoff predicted.

---

## 1. The five mechanisms that can move a byte

| # | Site (`src/Main.cpp` unless noted) | What changes | Bites when |
|---|---|---|---|
| M1 | `:2604`-`:2618`, the five non-collapsed accumulations, now gated on `pr.w0(r) != 0` | `MAC`, `MAC_case`, `MAC_control`, `Number_rare` | a set weights a NON-collapsed marker exactly 0 |
| M2 | `:2690` `NumUltraRare_GroupVec`, gated on `wUR != 0` | `Number_ultra_rare` | a set weights an ULTRA-RARE marker exactly 0 |
| M3 | `genoSumcount_noweight` widened to `jmr` and weight-gated (`:2111`, `:2637`, `:3014`, read at `:3517`) | the burden path's variance-ratio band, ER eligibility and min-group-MAC gate, hence `Pvalue_Burden`/`BETA`/`SE` on `--r.corr=1` | a set zeroes any marker AND the per-set MAC lands in a different variance-ratio band |
| M4 | `:3147`/`:3212`, the pseudo-marker's contribution taken from the UNWEIGHTED twin `genoURMat_noweights` (option b) | `MAC`, `MAC_case`, `MAC_control` | **any custom or transform set whose cell has an ultra-rare collapse** -- no zero weight needed |
| M5 | `R/SAIGE_SPATest_Region.R`, the min-group-MAC gate applied on the SKAT-O path (#126) + default 5 -> 1 | rows present/absent | a cell's per-set MAC falls below the gate |

Plus two non-golden effects: the FlexRV memory line's constant (4.0 -> 5.0
columns per cell, stdout only), and the #127 sidecar (section 5).

**M4 is the wide one and it is deliberate.** `MAC` used to be
`(sum of the cell's non-collapsed marker MACs) + (the set's WEIGHTED collapsed
pseudo-marker AC)` -- mixed units, which is what the fractional MACs in the
goldens are. Under (b) it is an allele count throughout. LEDGER #118 records
this as a knowing divergence from upstream on the user-weights path.

---

## 2. Files predicted UNCHANGED (byte-identical at gate 1)

* every `step0/` and `step1/` golden -- no code on those paths was touched;
* every `plink_marker_*` / `pgen_marker_*` golden -- the marker path is untouched;
* every `*.singleAssoc.txt` and `*.singleAssoc.txt_temp`. **Mechanism, not
  hope:** a zero weight already leaves `genoURMat` alone (`max(x, 0*g) == x`),
  so no pseudo-marker VECTOR changes; the UR pseudo-marker's own test context
  still reads the WEIGHTED `MAC` at `:3085`, which M4 does not touch; and
  `AC_Allele2` is `altCounts`, still weighted.
* every group golden whose only weight sets are beta / No_Weights sets:
  `plink_group_burden{,_quant,_quant_sparse}`, `plink_group_composed`,
  `plink_group_cond_binary`, `plink_group_cond_burden`, `plink_group_cond_quant`,
  `plink_group_fastTest{,_escalated}`, `plink_group_loco_chr1`,
  `plink_group_loco_chr2`, `plink_group_nofast{,_quant,_quant_sparse}`, and the
  pgen mirrors `pgen_group_burden`, `pgen_group_composed`,
  `pgen_group_cond_binary`, `pgen_group_cond_burden`, `pgen_group_fastTest`,
  `pgen_group_nofast`, `pgen_group_nofast_quant_sparse`.
  Why: a beta set's ultra-rare collapse is already unweighted, so the twin holds
  the same vector and M4 is the identity; `dbeta(MAF, a, b)` is never exactly 0
  on this fixture, so M1/M2/M3 never fire; and the smallest cell MAC in these
  files is 5, so M5 at gate 1 refuses nothing.

## 3. Files predicted to MOVE at gate 1 -- custom weight sets, no zero weights

Mechanism: **M4 only.** These group files carry no weight of exactly 0
(`group_new_snpid_withWeights.txt` and `group_new_snpid_multiWeights.txt` have
none), so M1/M2/M3 are inert and **no p-value, BETA or SE moves**. Only
`MAC`, `MAC_case`, `MAC_control` move; `Number_rare` and `Number_ultra_rare`
do not.

| file | rows whose MAC moves |
|---|---|
| `plink_group_multiweight.txt` | 24 of 24 cells |
| `plink_group_multiweight_burden.txt` | 24 |
| `plink_group_multiweight_quant.txt` | 24 (no MAC_case/MAC_control: quantitative) |
| `plink_group_cond_multiweight.txt` | 24 |
| `plink_group_cond_customweight.txt` | 12 (K = 1 custom set) |
| `pgen_group_multiweight.txt` | 24 |
| `pgen_group_multiweight_burden.txt` | 24 |
| `pgen_group_cond_multiweight.txt` | 24 |

**The predicted VALUES, exactly.** `group_new_snpid_multiWeights.txt` and
`group_new_snpid_withWeights.txt` have `var`/`anno` lines byte-identical to
`group_new_snpid.txt` (checked), and all eight plain goldens agree on these
columns to the digit. So every weight set of a cell must land on the plain
golden's numbers -- the CADD row and the REVEL row of a cell become EQUAL in
MAC, which is the sharpest single check here:

| Region | Group | max_MAF | MAC | MAC_case | MAC_control |
|---|---|---|---|---|---|
| GENE1 | lof | 0.001 | 5 | 2 | 3 |
| GENE1 | lof | 0.01 | 98 | 7 | 91 |
| GENE1 | missense;lof | 0.001 | 10 | 3 | 7 |
| GENE1 | missense;lof | 0.01 | 119 | 10 | 109 |
| GENE1 | missense;lof;synonymous | 0.001 | 10 | 3 | 7 |
| GENE1 | missense;lof;synonymous | 0.01 | 119 | 10 | 109 |
| GENE2 | lof | 0.001 | 5 | 0 | 5 |
| GENE2 | lof | 0.01 | 22 | 2 | 20 |
| GENE2 | missense;lof | 0.001 | 15 | 2 | 13 |
| GENE2 | missense;lof | 0.01 | 77 | 10 | 67 |
| GENE2 | missense;lof;synonymous | 0.001 | 15 | 2 | 13 |
| GENE2 | missense;lof;synonymous | 0.01 | 77 | 10 | 67 |

(Worked example: `GENE1 lof 0.001` today reads 38.4 for CADD and 3.34 for
REVEL, `Number_rare = 0`, `Number_ultra_rare = 4` -- pure pseudo-marker, so the
whole value is the weighted AC. Both become 5, the plain golden's value, because
the unweighted collapse of those 4 markers does not depend on the weight set.)

## 4. Files predicted to MOVE at gate 1 -- FlexRV (zero weights present)

| file | mechanisms | predicted |
|---|---|---|
| `plink_group_flexrv_burden.txt` (P58) | M1, M2, M3, M4 | `MAC`, `MAC_case`, `MAC_control`, `Number_rare`, `Number_ultra_rare` move on most of the 384 cell rows; `Pvalue_Burden`/`BETA`/`SE` move on the subset whose per-set MAC crosses a variance-ratio band (LEDGER #118 measured 79 of 384 band changes, SE ~0.29% on this model's two `null` ratios 0.9475/0.9421), and the Cauchy row follows |
| `pgen_group_flexrv_burden.txt` (P61) | same | must move **identically** to P58: same cells, same numbers |
| `plink_group_flexrv_burden_quant.txt` (P59) | M1, M3 only (collapse off, `--MACCutoff_to_CollapseUltraRare=0`, so there is no pseudo-marker and M4 cannot fire; quantitative, so no MAC_case/MAC_control) | `MAC` and `Number_rare` move on the zero-weight sets; `Pvalue_Burden`/`BETA`/`SE` move where the band changes |
| `plink_group_flexrv_skato.txt` (P60) | **M4 only** | 8 cell rows: `MAC`, `MAC_case`, `MAC_control` move; **no p-value moves**, because the reduced spec (`unw;sbeta:2` x `unw;mafbeta:2` at max MAF 0.02) contains no weight of exactly 0 -- which is why it has zero NA rows (LEDGER #117) |
| `plink_group_flexrv_skato.txt.singleAssoc.txt` | none | unchanged |

**Two invariants inside P58/P61 that must hold:**
* a set whose weight vector is identically 1 on the cell (`unw` score transform
  x any MAF transform that is 1 for every in-cell MAC -- `unw`, and the `MAC<=k`
  rungs at or above the cell's largest MAC of 20) must keep its `MAC` **exactly**:
  weighted and unweighted collapse coincide there. Those rows are the negative
  control inside the moving file.
* `Number_rare + Number_ultra_rare` can only DECREASE, never increase, for every
  set; and it must stay at its old value for every all-ones set.

## 5. New files the candidate will contain (#127 sidecar)

`<SAIGEOutputFile>.pooledTests.txt` is written only when some region has
something to disclose -- a refused cell, an NA that dropped out of the Cauchy
combination, a rule other than the plain Cauchy one, or nothing tested at all.
At gate 1 nothing is refused, so the trigger is the rule. Predicted, from the
committed goldens' own p-values (a cell p of exactly 1 turns CCT into
`min(1, n*min(p))`, LEDGER #86):

| run | regions disclosed | rule |
|---|---|---|
| `plink_group_nofast_quant.txt` | GENE1, GENE2 | min-p Bonferroni |
| `plink_group_nofast_quant_sparse.txt` | GENE1, GENE2 | min-p Bonferroni |
| `pgen_group_nofast_quant_sparse.txt` | GENE1, GENE2 | min-p Bonferroni |
| `plink_group_cond_quant.txt` | GENE1, GENE2 | min-p Bonferroni |
| `plink_group_multiweight_quant.txt` | GENE1, GENE2 | min-p Bonferroni |

Every other run prints its one-line summary and writes no file.
`compare_goldens.sh` will report these five as **EXTRA in candidate** until they
are committed as goldens -- that is the expected first-run failure, not a defect,
and the five are exactly the runs whose pooled row was already produced by a rule
the output never named.

**No pooled p-value changes value.** The `NA`-instead-of-0.5 change (#127, the
CCT of no inputs) needs a region in which nothing was tested; no committed
golden has one (checked: every region of every group golden has at least one
non-NA cell row).

## 6. Gate 5 (the second run of #126's matrix)

The handoff predicted "the SKAT-O path loses exactly the cells below MAC 5,
which is 2 rows in each of 5 goldens". **That prediction is void once #118 lands
in the same commit sequence**, and the reason is section 3: the candidate rows
were the `REVEL` sets at `lof/0.001` with MAC 3.34 and 3.55, and those are
weighted-unit values. Under option (b) they become 5 and 5, which the gate at 5
admits. (They are also 2 rows in **7** goldens, not 5.)

So the prediction for gate 5 is:

* every non-FlexRV golden **byte-identical to its own gate-1 output**. The
  smallest per-set MAC anywhere else is exactly 5.
* `plink_group_flexrv_skato.txt` byte-identical to gate 1 (its 4 sets sit at
  MAC ~80-180).
* `plink_group_flexrv_burden.txt` / `pgen_group_flexrv_burden.txt` /
  `plink_group_flexrv_burden_quant.txt` lose rows: **at most 32 of 384**
  (LEDGER #118's count of sets that pass on the cell MAC and fail on the per-set
  one), the same rows in P58 and P61, and their Cauchy rows move because the
  combination shrinks. A `.pooledTests.txt` sidecar appears for these runs with
  a non-zero `Cells_refused_minGroupMAC`.
* the number is "at most 32" rather than 32 because that count was taken on a
  per-set MAC whose pseudo contribution was not yet pinned; option (b)'s
  unweighted pseudo-marker can only make a set's MAC larger than the weighted
  one where weights are below 1, so refusals can only be fewer.

## 7. What would falsify the implementation rather than the prediction

* any `plink_marker_*`, `pgen_marker_*`, `step0/` or `step1/` file moving;
* any `*.singleAssoc.txt` moving;
* a p-value moving in section 3's files (would mean M3 fired without a zero
  weight, i.e. the per-set count is wrong);
* a section 3 MAC that is not the plain golden's value (would mean the twin is
  not the collapse, or the non-collapsed part was gated when it should not be);
* `Number_rare`/`Number_ultra_rare` moving in section 3's files;
* a row appearing or disappearing anywhere at gate 1 (M5 must be inert there --
  `indexNonZeroVec_arma.n_elem > 0` implies an unweighted MAC of at least 1, so
  gate 1 is inert **by construction** once the MAC is an allele count);
* the two paths disagreeing on WHICH cells they refuse at gate 5, compared on
  cell identity alone (never on p-values -- the two burden estimators differ by
  2.76x at MAC 2).

## 8. The reconciliation that now has to hold

Under (b), `MAC_GroupVec(jmr)` and the burden path's
`genoSumcount_noweight(jmr)` are the same quantity computed two ways --
`sum of pr.MAC` against `arma::accu` over the non-zero genotypes, both
unweighted, both per set. They must agree to within the ultra-rare collapse's
`imputeGenoAndFlip` residual (#109 dosage zeroing acts on the WEIGHTED
pseudo-marker only, so on the unweighted twin there is none). Any gap is a
located error, not rounding.

---

# OUTCOME (appended after the gate; nothing above was edited)

## Gate 1 -- the prediction held, file for file and value for value

`bash harness/make_goldens.sh <cand-gate1>` then
`bash harness/compare_goldens.sh harness/goldens <cand-gate1>`, scratch build.

* **12 files DIFFER, exactly the 12 named in sections 3 and 4**, and nothing
  else: every step-0, step-1, marker and plain-group golden is byte-identical,
  and so is **every `*.singleAssoc.txt`**.
* **5 files EXTRA, exactly the 5 named in section 5** -- the sidecars, on the
  five quantitative SKAT-O runs whose pooled row takes CCT's `min(1, n*min(p))`
  fallback. They now read `Rule = min-p Bonferroni`, which is LEDGER #86 made
  visible in the output for the first time.
* **Row counts unchanged in all 12** -- no row appeared or disappeared, so the
  gate at 1 is inert, as predicted from `n_elem > 0 => MAC >= 1`.
* **Columns that moved, per file:**

| file | moved |
|---|---|
| the 7 binary custom-weight goldens (section 3) | `MAC`, `MAC_case`, `MAC_control` |
| `plink_group_multiweight_quant.txt` | `MAC` (quantitative: no case/control columns) |
| `plink_group_flexrv_skato.txt` (P60) | `MAC`, `MAC_case`, `MAC_control` -- **no p-value**, as predicted from "no in-cell weight is exactly 0" |
| `plink_group_flexrv_burden.txt` / `pgen_group_flexrv_burden.txt` (P58/P61) | `Pvalue_Burden`, `BETA_Burden`, `SE_Burden`, `MAC`, `MAC_case`, `MAC_control`, `Number_rare`, `Number_ultra_rare` |
| `plink_group_flexrv_burden_quant.txt` (P59) | `Pvalue_Burden`, `BETA_Burden`, `SE_Burden`, `MAC`, `Number_rare` -- and NOT `Number_ultra_rare`, because the collapse is off |

* **Section 3's value table: 0 mismatches over all 8 files.** Every weight set of
  every cell landed on the plain golden's `MAC`/`MAC_case`/`MAC_control`.
* **P58/P61 are byte-identical to each other**, as they must be.
* **P58's negative control fired**: the 12 rows whose `MAC` did not move are
  exactly the all-ones sets -- `REVEL_unw_unw` and `REVEL_unw_MAC{25,50,100,200,500}`
  in both genes -- where the weighted and unweighted collapses coincide.
* **`Number_rare + Number_ultra_rare` never increased** (0 rows), and decreased
  on **174 of 384** P58 rows. That is the same 174 LEDGER #118 measured
  independently, from the weights, before any of this was written.

## Gate 5 -- the prediction held, including the correction to the row's own

Step-2 re-run with `--minGroupMAC_in_BurdenTest=5` appended to every step-2
command (optparse keeps the LAST value), compared against the gate-1 candidate.

* **Exactly 3 result files differ from gate 1**: `plink_group_flexrv_burden.txt`
  (P58), `pgen_group_flexrv_burden.txt` (P61) and
  `plink_group_flexrv_burden_quant.txt` (P59). **Every other group golden --
  including `plink_group_flexrv_skato.txt` and all seven custom-weight files --
  is byte-identical to its gate-1 output**, which is section 6's correction of
  this row's own gate-5 prediction: the `REVEL lof 0.001` rows it named are 5
  after #118, not 3.34, and 5 passes a gate of 5.
* **32 of 384 cells lost in each of the three**, which is LEDGER #118's
  independently measured count exactly, not merely the bound section 6 claimed.
  The lost cells' MACs run 1 to 4 and the smallest surviving MAC is exactly 5,
  so the gate cuts where it says it does.
* **Every surviving row is byte-equal to its gate-1 twin** -- the gate refuses
  cells, it does not perturb the ones it admits.
* P58 and P61 lose the same 32 and stay byte-identical to each other.
* Three `.pooledTests.txt` sidecars appear for those runs, each reading
  `Cells_refused_minGroupMAC = 16` per gene against `Rows = 176`, which is the
  disclosure #127 exists for: without it the pooled p-value would have shrunk
  its denominator by 8% with nothing in the output to say so.
* Why "at most 32" became exactly 32: option (b)'s unweighted pseudo-marker is
  larger than the weighted one wherever weights are below 1, so it could only
  have moved sets UP across the boundary -- and on this fixture none of the 32
  crossed it, because the refused sets are MAC-ladder rungs whose whole burden
  is a two- to four-allele pseudo-marker either way.

**One artefact to know about**: the gate-5 arm regenerates step 2 only (its
step-1 models are symlinked from the gate-1 candidate), so `make_goldens.sh`'s
own PGEN-vs-PLINK check reports 2 failures for the step-0 sparse-GRM files it
was never asked to produce. The full gate-1 run reports zero.
