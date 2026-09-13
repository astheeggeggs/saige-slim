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

## 6. The power arm is not interpretable, and the reason is not (only) low power

Stating the expected magnitude first, from the fixture's own `sim_meta.json`:
coding h2 = `h2_realised` 0.0220 x `frac_h2_coding` 0.6 = 0.0132 over
`n_causal_genes` 63 -> mean per-gene h2 2.1e-4 -> **mean NCP = N x h2 = 10.5**. At
Bonferroni alpha = 0.05/210 = 2.38e-4 the 1-df critical value is 13.50, giving
power 0.332 at the mean NCP. **Expected detections: 21 if every gene sat at the
mean, 48 under the `ncp_ladder = 1,60` spread.**

Observed: **`y_quant` 1-3 of 63**, `y_binary` 8-11 of 63.

That is a 10-40x shortfall and it is NOT explained. Per this project's own rule a
result that lands that far off is a located error until shown otherwise, so the
power arm must not be read as "underpowered by design" until one of these is
tested:

1. **The tested MAF window is not the h2 window.** The sweep tests
   `--flexRV_maxMAF=0.001` (MAC <= 100 at N = 50,000) while the gene's coding h2
   spans all frequencies. If most coding h2 sits above MAF 0.001 the test never
   sees it. **Check this first** -- it is a one-pass calculation over the
   fixture's `variants.tsv.gz` and needs no run.
2. **`frac_directional = 0.5`**: half the causal genes have mixed-sign effects,
   which a burden cancels by construction, and `link_mix` is step/lof/graded so
   the genotype-to-burden map is not linear either.
3. The NCP ladder may be concentrated at the low end rather than uniform on
   [1, 60].

The `oracle` column bounds where the ceiling is and is the sharpest number in the
file: on `y_quant` c0 the oracle -- best single set per gene, Bonferroni over
sets -- detects **2**. flex192's 2 is AT that ceiling. No weighting scheme can do
better on this fixture, which is why CONTROL 2 (rho > 0.5 expected, observed
-0.04 to 0.25 on `y_quant`) and CONTROL 3 are being evaluated at a resolution of
one or two genes.

## 7. What to do next on this

- **Do not re-run the sweep for power.** It is re-run, its null half is valid and
  current, and its power half is at the oracle ceiling.
- **Resolve item 1 above before designing anything**, because if the h2 is mostly
  outside the tested MAF window then the fixture is not underpowered, it is
  mis-targeted, and the fix is a different `--flexRV_maxMAF` rather than a bigger
  effect.
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
