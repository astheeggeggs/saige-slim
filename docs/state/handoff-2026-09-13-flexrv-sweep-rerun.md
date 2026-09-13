# The FlexRV sweep re-run at fixed HEAD (2026-09-13)

The sweep the #129/#130 handoff left owed has been re-run. This is what was
run, what moved, why it moved, and the one scientific question the results
raise that this fixture cannot answer.

**Read `docs/state/next.md`'s 2026-09-13 correction block first** -- it records
that the previous handoff's reason for doing #134 before the sweep was false,
which is why the sweep could be run now.

## 1. What was run

Only the **72 collapse-10 runs**, not all 144. The 72 collapse-0 runs are
provably unaffected by #129/#130: every non-comment hunk of both fixes is inside
`if(i2 > 0)` at `src/Main.cpp:3308`, and the committed c0 logs carry **16,992
lines all reading "0 markers are ultra-rare"**. That is a measurement, not an
argument: with `--MACCutoff_to_CollapseUltraRare=0` nothing is ever collapsed, so
the block never executes.

The pre-fix c10 outputs were MOVED, not deleted, to
`harness/work/flexrv_study/results-prefix-38eee5f5/` (232 files). They are the
BEFORE state for section 3 and should be kept until someone decides otherwise.

**Provenance of the new runs: all 72 at `head=976e73ce dirty=0`**, one
`so_md5=07423e2f...` matching the installed build, exit 0 on every run. The
previous sweep ran at `38eee5f5` with `dirty=1..6` -- an uncommitted tree -- so
the new results are the first with clean provenance.

**Two operational notes.**
- No rebuild was needed: the installed `.so` (Sep 12 17:21) postdates every
  tracked `src`/`R` source, the last of which is `a49c2cec` at 17:19.
- **The load guard fired and cost one pass.** `run.sh` refuses when the 1-minute
  load exceeds `MINLOAD_REFUSE` (default 2.5). This box idles at 3-5 from desktop
  apps alone, so the first attempt raised it to 6 and still lost **56 of 72 runs**
  to `REFUSE: load 7.58 > 6` -- the sweep's own two workers on top of the
  baseline. `study.sh` is idempotent, so the second pass at
  `MINLOAD_REFUSE=100` picked them up with 0 failures. **Every timing recorded
  in these logs is void by construction** and must not be compared with
  `harness/bench/NOTES.md`. If a future run needs timings, it needs a quiet box,
  not a raised guard.

## 2. A method error worth recording, because it inverted a conclusion

The first before/after comparison zipped the two files POSITIONALLY. The fix
EMITS ROWS THAT DID NOT EXIST BEFORE (section 3), so the row counts differ and
every row after the first insertion was compared against its neighbour. That
reported **535,236 of 679,168 rows differing (73%) and a max |dlog10 p| of
51.7**, and led to a confident, wrong hypothesis: that the same weight set gave
different answers depending on how many sibling sets were in the run, i.e.
cross-set contamination through the shared worker scratch.

The hypothesis was killed by its own control. `SIM_unw_unw` is a set of
identically-1 weights that exists in a 1-set run (`flex1`), a 16-set run
(`flexscore`) and a 192-set run (`flex192`). Comparing it ACROSS those runs, keyed
on (Region, Group, max_MAF, Weight): **2,520 cells, 0 differences, in the pre-fix
tree AND at HEAD.** There was never any cross-set dependence. The correct
comparison is keyed, and it is section 3.

The general lesson is the project's own: a positional comparison of two files
that may differ in length is not a comparison. Key it, and check
`only_old`/`only_new` explicitly -- here `only_old` is 0 everywhere and
`only_new` is 4-45 per file, which is itself the signal that the length changed.

## 3. What actually moved

Keyed on (Region, Group, max_MAF, Weight), over all 72 c10 runs:

| | |
|---|---|
| files with any change | **28 of 72** |
| rows compared by key | 728,632 |
| rows differing | **87,267 (12.0%)** |
| rows present only in the NEW output | 4-45 per changed file; **0 removed** |
| max abs difference in log10 p | **5.77** (`GENE0181 high;medium 0.001 SIM_sbeta10_MAC25`) |

**The movement is confined to sets with FRACTIONAL weights, and is monotone in
how hard the transform shrinks them.** By score transform, share of shared rows
that moved:

```
sbeta20  36.0%   cbrt0.9 19.2%   cbrt0.6 11.3%   cbrt0.3  5.4%
sbeta10  30.6%   cbrt0.8 15.8%   cbrt0.5  9.4%   cbrt0.2  4.0%
sbeta5   22.2%   sbeta3  14.0%   cbrt0.4  7.4%   sbeta1.5 0.9%
                 cbrt0.7 13.5%   sbeta2   6.1%
unw       0.0%  (73,008 rows)    lof      0.0%  (31,752 rows)
```

**`unw` and `lof` are the internal negative control and they are exactly zero
over 104,760 rows.** `unw` has `w == 1`; `lof` has `w` in {0, 1}. Neither ever
produces a fractional `w*g`, and neither moved a single row. Any explanation of
this diff has to predict that, and the one below does.

**The mechanism is LEDGER #129, and the ordering is predictable from one
number.** #129 was the dosage filter reaching the collapsed pseudo-marker;
`--dosage_zerod_cutoff` defaults to **0.2** (`extdata/step2_SPAtests.R:482`), and
the collapse is `max(w * g)`. The score-beta transform is `w = |s|^(b-1)`, so at
a mid score `s = 0.5`:

```
sbeta1.5 w=0.707  above 0.2      sbeta5  w=0.0625   BELOW 0.2
sbeta2   w=0.500  above 0.2      sbeta10 w=0.00195  BELOW 0.2
sbeta3   w=0.250  above 0.2      sbeta20 w=1.9e-06  BELOW 0.2
```

The observed movement jumps between sbeta3 (14.0%) and sbeta5 (22.2%) -- exactly
where the weighted genotype crosses the cutoff and the old code began zeroing it.
The MAF axis agrees: `MAC1` moves most (29.9%), because it keeps only MAC <= 1
markers and multiplies an already-small score weight.

So #129 dominates and #130 (the fold, which needs the pseudo-marker's weighted
alt frequency above 0.5) contributes the smaller part. **This was not predicted
before the run and should have been** -- the direction was predictable from the
two fixes plus the cutoff.

## 4. What did NOT move: the null calibration

`harness/sim/flexrv_study/results-2026-09-13.txt` against
`results-2026-09-12.txt`: **the CONTROL 1 and CONTROL 1b null-calibration lines
are identical except one**, and the exceptional one moves within its CI
(`y_null c10 flex192 uninformative` at alpha 0.01, 0.0190 -> 0.0238).

**8 lines of ~290 fall outside a 95% binomial CI, before and after, the same 8.**
At a 95% interval ~14 of 290 would be expected outside by chance, so 8 is if
anything conservative and is not evidence of miscalibration.

That is the expected result and is the point of the re-run: the fixes change
individual pseudo-marker statistics without biasing the null distribution. **The
null calibration is now established at a clean HEAD**, which it was not before.

## 5. What improved, and the question that did not

CONTROL 3 (`results-2026-09-13.txt`):

```
y_quant   c0  flex192: informative  2, realistic  2, uninformative  2, unw 1  monotone=TRUE  uninf<=unw=FALSE
y_quant   c10 flex192: informative  2, realistic  2, uninformative  2, unw 1  monotone=TRUE  uninf<=unw=FALSE
y_binary  c0  flex192: informative 11, realistic 11, uninformative 10, unw 8  monotone=TRUE  uninf<=unw=FALSE
y_binary  c10 flex192: informative 11, realistic 11, uninformative 10, unw 8  monotone=TRUE  uninf<=unw=FALSE
```

Two real gains. **All four lines are now `monotone=TRUE`**; two were FALSE before,
both because the uninformative arm over-detected (`y_quant` c10 3 -> 2,
`y_binary` c10 uninformative 11 -> 10). And **c0 and c10 now agree exactly on
both traits**, where before they did not -- the collapse cutoff should not change
the answer much, and after the fixes it does not.

**`uninf<=unw=FALSE` on every line, unchanged.** This is the finding that
outlives the re-run and it is stated carefully because the fixture cannot settle
it: the uninformative arm has its scores PERMUTED WITHIN GENE (cor with the true
class ~0), and it still detects 10 genes against the unweighted burden's 8. So at
this fixture's resolution, FlexRV's advantage over an unweighted burden does not
depend on the score carrying information. That is consistent with the advantage
coming from the 192-set Cauchy combination itself -- which a permuted score gets
in full -- rather than from the weighting. It is ALSO consistent with the fixture
being too underpowered to separate the two. **Distinguishing those two is the
power question, and it needs a different fixture.**

## 6. The power arm, RESOLVED 2026-09-13: the fixture is behaving correctly

Section 6 originally said the power arm showed a "10-40x shortfall" against the
fixture's h2 budget that was "NOT explained", and named three candidate causes.
**That was wrong, it has been run down, and the answer is that there is no
anomaly.** The corrected account, with each of the three candidates resolved:

### 6a. The variance decomposition (a sum over the truth table, no model)

Summing `2*MAF*(1-MAF)*beta^2` over `variants.tsv.gz` gives 0.02248, against
`sim_meta.json`'s `h2_realised` 0.02200 -- that reconciliation is the control
that the formula and the beta scale are right. Split by what the sweep can see
(`--annotation_in_groupTest=high:medium`, `--flexRV_maxMAF=0.001`):

| bucket | variants | share of causal variance |
|---|---|---|
| **visible: mask AND MAF <= 0.001** | 498 | **30.2%** |
| outside mask, MAF > 0.001 | 105 | 37.5% |
| outside mask, MAF <= 0.001 | 972 | 26.4% |
| in mask, MAF > 0.001 | 15 | **6.0%** |

**CANDIDATE 1 IS REFUTED. The MAF window costs 6 points; the annotation mask
costs 64.** `--flexRV_maxMAF=0.001` is not what starves the test -- the
simulation put most causal variance in `low` and `noncoding`, which a gene-based
test on a `high:medium` mask is not supposed to see. So this is not a defect, it
is the wrong denominator: the right one for this sweep is 30.2% of 0.0225 =
0.0068, not `h2 * frac_h2_coding` = 0.0132.

### 6b. Candidates 2 and 3 are both CONFIRMED, and both are real

Per causal gene, over visible variants only, unweighted burden
(`NCP = N * Cov^2 / Var`, `Cov = sum b*2p(1-p)`, `Var = sum 2p(1-p)`):

- **22 of 62 genes have mixed-sign visible betas**, and directional cancellation
  takes expected detections from 9.1 (no-cancellation bound) to **4.9**.
- **The ladder is heavily skewed**: median gene NCP **0.62**, with only a handful
  appreciable (34.3, 34.2, 15.5, 12.5, ...).

### 6c. What was actually wrong, and it was MY ESTIMATOR, not the fixture

The 4.9 figure still did not match the observed 1-2 on `y_quant`, so the
estimator was validated against the data before being trusted further -- and it
FAILED: Spearman(predicted NCP, observed -log10 p) = 0.243, and the two
top-predicted genes, **both at NCP ~34, gave observed -log10 p of 0.13 and
9.01**. Same prediction, one null and one a 1e-9 hit, is a wrong model rather
than a mis-scaled one.

The cause, found by checking the one thing never checked -- whether the variants
being summed over are variants a test can actually resolve: **GENE0190's entire
NCP of 34.3 comes from ONE variant at MAF 2e-05, i.e. MAC 2, two carriers, with
beta 4.14.** The chi-square NCP approximation is meaningless at two carriers,
which is exactly what SPA and the exact test exist for. GENE0210's NCP comes from
MAC 9/29/74 variants and duly detects at p = 1e-9.

Quantified, splitting the naive 4.9 by whether a gene's signal is resolvable:

| gene set | genes | E[detections] |
|---|---|---|
| all 62 | 62 | 4.9 |
| total visible MAC >= 20 | 35 | **1.8** |
| total visible MAC >= 50 | 19 | 1.4 |
| total visible MAC >= 100 | 9 | 1.1 |

**3.1 of the naive 4.9 comes from 27 genes whose entire visible signal sits below
MAC 20.** Against resolvable genes only the prediction is **1.8**, and the
observed `y_quant` detections are **1 (`unw`) and 2 (`flex192`)`. That
reconciles. `y_binary`'s 8-11 follows from `h2_binary_mult` 13.9 lifting every
beta above threshold.

**So there is no unexplained shortfall.** The fixture is doing what its
parameters say it should; the earlier "10-40x" and the intermediate "~3x" were
both artefacts of crediting signal that no test can resolve.

### 6d. The collapse is a no-op here, which also refutes a fourth guess

`c0` (nothing collapsed) against `c10`, same arm: **193-194 of 210 genes have
IDENTICAL p-values**, median |dlog10 p| 0.0000, max 0.264. At MAC <= 10 out of
2N = 100,000 essentially no sample carries two of a gene's ultra-rare variants,
so `max()` and `sum()` coincide. The collapse cutoff cannot explain anything on
this fixture.

### 6e. The design constraint this hands to a power fixture

**Control the MAC of the causal variants, not just their h2.** A beta of 4.14 on
a MAC-2 variant generates heritability on paper and nothing detectable in
practice, and 27 of this fixture's 62 causal genes are built that way. A fixture
intended to separate weighting schemes needs its causal signal in variants the
test can resolve -- concentrated in fewer genes at higher MAC -- or the oracle
ceiling stays at a handful of genes however much h2 is nominally present.

Three hypotheses were refuted in the course of this (the MAF window, the
collapse, and the estimator itself). The missing control, in hindsight, was the
cheap one: **reconcile the variants being summed against the variants a test can
resolve, before predicting anything from them.**

## 7. What to do next on this

- **Do not re-run the sweep for power.** It is re-run, its null half is valid and
  current, and its power half is at the oracle ceiling.
- **Item 1 is DONE (section 6) and the answer was no**: the MAF window costs 6
  points of causal variance and the annotation mask costs 64, so a different
  `--flexRV_maxMAF` fixes nothing. What a new fixture needs is in 6e -- causal
  signal at resolvable MAC, not more h2.
- **The power question, stated so it can be answered**: does FlexRV's advantage
  over the unweighted burden survive permuting the score? Today's answer is
  "yes", which would be a statement about the method rather than about this
  branch's code. A fixture sized so the oracle detects tens of genes rather than
  two would settle it.
- LEDGER **#142** (the `mac:` ladder degeneracy) bears directly on that question,
  and **this sweep confirms it at N = 50,000**, which the row could previously
  only show on the bundled N = 1,000 fixture. At `--flexRV_maxMAF=0.001` the
  ingest bound is MAC <= 2*N*maxMAF = 100, so every rung at or above 100 is
  identically 1. Measured on `y_binary_informative_flex192_c10_r1.txt`, share of
  cells where the rung's row is identical to its own score transform's `unw` row
  on every value column:

  ```
  MAC100 100.0%   MAC50 74.3%   MAC5  8.9%    mafbeta2 0.0%
  MAC200 100.0%   MAC25 48.5%   MAC3  5.0%    mafbeta5 0.0%
  MAC500 100.0%   MAC10 20.5%   MAC1  1.6%
  ```

  So **three of the twelve MAF transforms are fully degenerate here (48 of the
  192 sets per gene) and three more are degenerate in half the cells or more**,
  and all of them are combined into the pooled Cauchy row. The duplicate
  structure is part of what the permuted-score arm is picking up, which is the
  concrete reason #142 and the power question have to be settled together
  rather than separately.
