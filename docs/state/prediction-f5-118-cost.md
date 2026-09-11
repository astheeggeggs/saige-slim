# Prediction: what LEDGER #118 cost the FlexRV memory profile (written before measuring)

Arms: `harness/work/flexrv50k/run.sh k192`, sim50k `y_null`, 210 genes, ONE cell,
192 transform sets, N = 50,000, 1 thread. Control = `35654e19` (the commit
before #127/#126/#118); candidate = `0c1edfb2`.

Baseline from `harness/bench/NOTES.md` (2026-09-10, same recipe): **35.1 s user,
805 MB peak RSS**.

## The mechanism, and the number it implies
#118 widened `genoURMat_noweights` from `q_anno_maf` (per cell) to
`q_anno_maf_weight` (per cell x weight set), a dense N x cells double matrix
allocated and zero-filled per region:
  50,000 x 192 x 8 B = **76.8 MB**, against 50,000 x 1 x 8 B = 0.4 MB before.

## Predictions
1. **Peak RSS: +70 to +80 MB**, i.e. ~805 -> ~880 MB (+9.5%). Sign positive,
   one array's worth, no more: nothing else #118 touched is N-sized.
2. **User time: +0 to +2 s** (35.1 -> 35-37 s). The extra cost is one
   `fill::zeros` of 76.8 MB per region x 210 regions = ~16 GB of memset; at
   ~10 GB/s that is ~1.6 s, ~5% of a 35 s run. It is NOT free, but it is small
   against the 28 s the 191 extra sets already cost.
3. **The k1 and plain arms must NOT move** (805 -> unchanged at 506 MB): at
   q_weight = 1, `q_anno_maf_weight == q_anno_maf`, so the array is the same
   size it always was. This is the negative control and it is what separates
   "#118 widened one array" from "something else regressed".

## What would falsify the implementation rather than the prediction
* RSS up by much more than one array (~154 MB+) -> something else became per-set.
* k1 or plain moving at all -> the widening is not confined to q_weight > 1.
* User time up by more than ~10% -> the cost is not just the zero-fill.

## Bonus: the first look at #118 on realistic data
The same run is a 210-gene, N = 50,000 instance of #118 with real zero weights
(the MAC ladder and the `high` lof mask). Expect `MAC`, `Number_rare` and
`Number_ultra_rare` to move on most of the 192 x 210 rows, and `Pvalue_Burden`
to move only slightly: this model's two `null` variance ratios are 1.0892 and
1.0951, 0.54% apart, so a band change is worth ~0.3% in SE. A LARGE p-value
move here would mean the band is not the only route and would need explaining.

---

# OUTCOME (appended after measuring; nothing above was edited)

Paired, alternated, 2 reps, control `35654e19` vs HEAD `0c1edfb2`, same fixture,
libraries the only difference. `prepare` is page-cache noise (5.2-11.2 s) and
totals with it; `urPass`/`emit` are stable within arm. Rep 2 of each arm ran at
comparable load and is the estimate; rep 1 of head is a high-load outlier.

| bucket | ctl | head | delta |
|---|---|---|---|
| peak RSS | 797 / 801 MB | 844 / 877 MB | **+47 to +76 MB** |
| urPass | 13.61 / 13.49 | 20.21 / 15.93 | **+2.4 s** (rep 2) |
| emit | 19.17 / 18.88 | 32.74 / 22.78 | +3.9 s (rep 2) |
| .. of which emitTest | 13.15 / 12.87 | 24.31 / 16.91 | **+4.0 s** (rep 2) |
| emit - emitTest | 6.02 / 6.01 | 8.43 / 5.87 | ~0 |
| k1 arm RSS | 507 / 507 MB | 505 / 505 MB | **0** |
| plain arm RSS | 504 / 504 MB | 503 / 505 MB | **0** |

**Prediction 1 (RSS +70-80 MB): held at the top of the range, 2 reps spanning
+47 to +76.** One array's worth, nothing more.

**Prediction 3 (k1 and plain unchanged): held exactly.** The negative control
fired, so the growth is confined to q_weight > 1, which is what separates
"#118 widened one array" from "something else regressed".

**Prediction 2 (user time +0 to +2 s): FALSIFIED.** User CPU went 41.4 s -> 52.0 s
(+26%). The buckets split it in two, and only one half is overhead:

* **emitTest +4.0 s is the cost of being CORRECT, not a defect.**
  `ctx.useSparseGRM = (MAC > m_cateVarRatioMinMACVecExclude.back())` -- back() is
  20.5 at the default -- so a set whose true MAC is small now takes the
  sparse-GRM score test that the cell's MAC used to let it skip. MEASURED on the
  same two output files: rows with MAC <= 20.5 went **2,722 -> 4,287 (+57%)**.
  That is the fix working: a set with 3 alleles should not be tested with the
  common-variant shortcut because the CELL around it has 93.
* **urPass +2.4 s IS overhead I introduced and it is removable.** The
  `genoSumcount_noweight` accumulation moved INSIDE the r loop, so what ran once
  per (anno x max_MAF) cell now runs 192 times: an `arma::find` over an N-vector
  plus a 400 KB `.col()` copy, 210 regions x 192 sets = 40,320 times instead of
  210. Order of magnitude checks out (~100 us each -> ~4 s).
  FIX: `arma::accu(genoURMat_noweights.col(jmr))` over the subview -- the same
  sum, since adding exact zeros is exact in IEEE, with no `find`, no uvec
  allocation and no copy. Expect to recover most of the 2.4 s (~5% of the run).

## A CORRECTION this forces to what was said earlier today

I told the owner the variance-ratio consequence of #118 was "nearly invisible on
every fixture we have", from the two `null` ratios being 0.54% apart at N=50k.
That is right about the BAND and wrong about the total, because there are TWO
routes off the per-set MAC and I had only separated one:
  (a) which variance-ratio band is selected -- worth ~0.3%, as stated;
  (b) whether `useSparseGRM` is on at all -- a different variance ESTIMATOR.
Measured over 40,196 cell rows at N = 50,000: `MAC` moved on **96.4%** of rows,
`Pvalue_Burden` on **8.8%**, median relative move 0.8% but **max 20%**. A 20%
move cannot come from (a); it comes from (b). So #118 at realistic N is a
statistical change to ~9% of burden p-values, not the reporting-and-gating fix
the bundled fixture made it look like.
