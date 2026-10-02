# Handoff 2026-10-02: LEDGER #176 explained -- the variance-ratio rescaling moves `q` off the lattice

Read with `python3 harness/tools/ledger.py 176`. **The row's own NEXT STEP was executed and it
answered the question.** The row says "the 3.12 SE factor has no account at all"; it has one now,
confirmed in the source and reproduced arithmetically on four cells. The row's `OBSERVATION ONLY`
framing of the 2026-10-01 counter is right about the counter and wrong about the defect: the
cells it counts are **mis-tested**, not merely un-routed.

Nothing in this document changes code. The repair is costed below and not written.

## The mechanism, from source

`src/SAIGE_test.cpp:507`  `double var1 = var2 * t_ctx.varRatioVal;`  -> `var1/var2 == varRatio`
`src/SAIGE_test.cpp:1039` `m1 = dot(m_mu, t_gtilde);`                -> the null mean
`src/SAIGE_test.cpp:1050` `q = t_Tstat/sqrt(t_var1/t_var2) + m1;`    -> `q = Tstat/sqrt(varRatio) + m1`

So the **observed score is divided by `sqrt(varRatio)`** and then compared against a CGF built
from the **unrescaled** `gtilde` and `mu`. That is the standard location-scale trick and it is
correct for a smooth distribution: shifting the quantile is equivalent to scaling the variance.
It is wrong for a distribution with atoms, because it moves `q` **off the lattice the statistic
can actually attain**.

Verified against the four dumps, `Tstat/sqrt(vr) + m1` against the header's `q`:

| dump | predicted | q (header) | rel diff |
|---|---|---|---|
| 183_136 | 0.999807 | 0.999807 | 3.1e-07 |
| 183_137 | 0.985358 | 0.985358 | 4.1e-07 |
| 183_138 | 0.943254 | 0.943254 | 3.0e-07 |
| 183_139 | 0.999805 | 0.999805 | 2.4e-07 |

## Why it is catastrophic only on concentrated cells

The rescaling inflates `q` by about 1% of the statistic (0.0091 to 0.0097 on these four). The
discrete part absorbs that using its non-atom carriers -- unless they hold nothing.

On the four cells nearly all of the weight sits in the largest term, so everything else among
the carriers sums to ~0.003, less than the ~0.0096 needed. The only
remaining route to `q` is an excursion of the Gaussian remainder, whose sd is ~0.0012:

```
183_136: q - obs = 0.00965 = 8.27 sd of the remainder  ->  normal tail 6.6e-17
183_137: q - obs = 0.00951 = 8.27 sd                   ->  6.5e-17
183_138: q - obs = 0.00910 = 8.27 sd                   ->  6.6e-17
183_139: q - obs = 0.00965 = 8.27 sd                   ->  6.6e-17
```

**Exactly 8.27 sd on all four.** Both the inflation and the remainder's sd scale with the
vector, so the ratio is invariant -- a signature, not a coincidence. 6.6e-17 is the reported
1.1e-17 before the saddlepoint's correction.

### The two contrast cells, which discriminate this from the alternatives

- **183_0**, `nz` `<20` -- BELOW the 13-carrier cap, so #172's routing fired and sent it to the
  exact test, **which ignores the variance ratio**. Its header has `q` equal to the observed
  statistic to the last digit, and its p is 8.7e-3. The fix that exists already works; the cap
  is what withholds it.
- **183_130**, less concentrated than the broken four. Same ~1% inflation, but 0.22 of slack
  in its non-atom weight, so it absorbs it: p 2.1e-3.

## The true p, by four independent routes

| route | value |
|---|---|
| exact convolution at the **unrescaled** `q` | 5.99e-03 |
| the closed form of the same model (`q` sits at the mean of the atom-fired branch) | ~6.1e-03 |
| the `sbeta5` cells of the same gene, which are less concentrated and were not broken | 9.4e-03 to 1.08e-02 |
| binomial tail implied by `MAC_case` (`<20` of a cell MAC of 84, 9.4 expected) | >= 2.1e-03 |

Against a reported 1.1e-17. **Fifteen orders**, from a 1% misalignment, because the
distribution had no middle.

## The evidence chain, and what a reviewer can re-run

| step | artefact | reproducible outside the TRE |
|---|---|---|
| `q = Tstat/sqrt(varRatio) + m1` | `src/SAIGE_test.cpp:507,1039,1050` | YES, source |
| the four-row arithmetic check | table above | YES, the numbers are here |
| exact-convolution oracle | `harness/tools/spa_exact_ref.py` | YES |
| oracle's controls (closed-form binomial to 5e-34, null at 0.5, MC agreement at 1e-2 and 1e-3) | `harness/tools/spa_exact_ref_selftest.py` | YES |
| Fix 1 leak sweep | same file, phase 2 | YES |
| the dumps themselves | `amr/flexrv_abovecap/*.spaDump.tgz` | NO -- individual level, stays in the TRE |
| `provenance()` output (concentration measures, z) | printed by `spa_exact_ref.py` inside the TRE, not written here | NO, but its inputs are the dumps |

## Three reversals, and the reason for them

I concluded, in order: the saddlepoint is broken; the saddlepoint is right; the variance ratio
band is wrong. All three were wrong, and the pattern is the useful part.

**Every check was downstream of `q`.** The saddlepoint, the exact convolution and the
BETA/SE-to-p reconciliation all take `q` and `gtilde` as given and ask whether the tail is
computed correctly. They agreed because they are the same question. In particular the
convolution's agreement with the kernel was **guaranteed, not informative** -- I handed it the
same inflated `q`, so it was a correct answer to the wrong question.

The provenance check (`spa_exact_ref.py::provenance`) was the first thing that looked upstream
of `q`, comparing `gtilde`, the raw `g`, and `sum(gtilde*(y-mu))` against the header, and it
found the inconsistency on the first run.

**A two-way reconciliation only tests what the two ways do not share.** Both ways shared the
inputs, so they could only ever validate the arithmetic. The same shape as the repo's standing
trap that a check which can return clean from a wrong path is not a check.

The second lesson: `n_eff` is a concentration measure, not a count of attainable values. My
first mechanism ("the distribution is two-point, so Lugannani-Rice has no justification") was
wrong because carriers whose weights span 0.0027 to 1 give a rich lattice. The thing that mattered was
not the number of values but whether the non-atom mass could cover the inflation.

## The repair, costed

**Fix 1 -- stop refusing these cells. Enumerate the heaviest 13, Gaussian for the rest.** The
routing predicate (`n_eff <= --max_MAC_for_ER` and concentrated) already identifies them
correctly; the 13-carrier cap exists only because the exact test enumerates `2^n`. Make the cap
a SPLIT POINT rather than a refusal: the exact test over the top 13 by weight, the remaining
material carriers folded into the Gaussian remainder -- the same carrier/remainder split the
fast kernel already uses, one level deeper, and sound here for the reason the bug exists (the
low-weight carriers are the smooth part). The exact test ignores the variance ratio, which is
what these cells need.

Measured against the exact oracle on a SYNTHETIC near-one-atom shape (one heavy weight, many
negligible ones; round values, `spa_exact_ref_selftest.py`, rerun 2026-10-02):

```
  q = OBSERVED statistic              oracle 5.1155e-03   fix1 5.1268e-03   ratio 1.002
  q = observed / sqrt(varRatio) [BUG] oracle 6.5210e-26   fix1 6.6464e-26   ratio 1.019
```

**0.2% on that shape.** Over an adversarial sweep (one atom plus k equal secondary
carriers, k 2 to 60, secondary mass 0.05 to 1.0, restricted to shapes the routing would admit)
the worst error was 519x deep in the tail and 2.3x in the range p >= 1e-3 (the review handoff
disputes that this is the decision-relevant range). The sweep used equal secondary weights and a single `mu`;
**that family is a plausible adversary, not a proven worst case.**

**Fix 2 -- convolution instead of enumeration.** A weighted sum of independent Bernoullis is a
convolution: DP over a value grid, O(n x grid), not `2^n`. That deletes `kExactByWeightCarriers`
rather than reinterpreting it. `spa_exact_ref.py` is a working, validated reference
implementation. Only worth it if Fix 1 is measured to leak on real cells.

**Fix 3 -- the guard, either way.** A cell the routing WANTS and cannot take must not silently
return a saddlepoint p. That is #163/#164's principle and this is a new instance of it. The
counter built on 2026-10-01 already identifies exactly the right cells; it needs to gate rather
than observe.

**The leak rate should be measured, not estimated.** Build Fix 1, then run the oracle beside it
on every above-cap cell in the rollout and compare. Any cell disagreeing by more than a few
percent gets Fix 2.

## Open questions for the owner

1. **Should the exact test ignore the variance ratio at all?** It does, by precedent (#172, and
   before that for tiny unweighted sets). Fix 1 depends on that. It means two cells in one gene
   are corrected differently for relatedness, which is worth deciding rather than inheriting.
2. **Is `q = Tstat/sqrt(varRatio)` wrong, or only wrong here?** It is fine on smooth
   distributions. The repair might belong in the routing rather than the formula -- or the
   formula may be wrong for every discrete test and this is merely where it became visible.
3. **Is 13 still the right number** once it is a split point rather than an enumeration budget?
   It would want a justification from how much weight the top-k holds, not from `2^13`.

## Blast radius: not yet measured, and measurable without reruns

How many above-cap cells exist across the rollout, in how many traits? The
`Concentrated cells above the exact test's enumeration cap: N ...` line is in every v7 job log,
so the count is recoverable from logs already in the bucket. The tidied notebook's `counters()`
cell scrapes it. That decides whether this is fix-and-move-on or rerun-and-splice. **Until it is
measured, do not read any above-cap cell's p-value as real** -- which is the row's existing
advice and still stands.

## What else came out of the day, and is NOT in the repo

- **`.exactByWeightAboveCap.txt` was never delivered.** `R/ExactByWeight.R:71` writes it inside
  the job; the AoU step-2 wrapper's sidecar loop had no entry for it and `submit()` declared no
  output (it declared `OUT_FILE_EXACTBYWEIGHT` twice instead). Every run before 2026-10-02
  discarded it. Fixed in the workbench notebook, NOT in this repo -- the wrapper lives in a
  `%%writefile` cell.
- **`.spaPinned.txt` (#174) has the same gap** and was fixed the same way.
- **`SAIGE_SPA_DUMP_BELOW` was unreachable from a dsub job** -- `src/SAIGE_test.cpp:131` reads
  it, the wrapper exported only `SAIGE_SPA_DUMP_DIR`. Without that, this investigation could not
  have run at all.
- The notebook rebuild that carries those repairs is `harness/tools/tidy_04b_v2.py`.

## Proposed STATUS amendment for #176, not applied

Applying it is a judgement call: the convention is one commit per row with the characterization
test inverted, and nothing is fixed yet. Suggested wording to append:

> **EXPLAINED 2026-10-02, and the defect is worse than the counter implied.** The row's own NEXT
> STEP (dump with `SAIGE_SPA_DUMP_BELOW`, replay) was run on amr EFRMB `ENSG00000213030`.
> `src/SAIGE_test.cpp:1050` sets `q = Tstat/sqrt(varRatio) + m1` while the CGF is built on the
> unrescaled `gtilde`/`mu`; the ~1% inflation is absorbed by a cell's non-atom carriers unless
> they hold nothing, and on the most concentrated cells it is 8.27 sd of the Gaussian remainder -- identically on
> all four cells -- giving 1e-17 where 6.0e-3 is right. Four independent routes agree on
> 6.0e-3. So the above-cap cells are MIS-TESTED, not merely un-routed, and `OBSERVATION ONLY`
> describes the counter, not the defect. The 3.12 SE factor is accounted for: it is the same
> inflation seen through `SE_Burden`. Contrast cells discriminate it -- a cell below the cap was
> routed to the exact test and its `q` equals the observed statistic exactly, and a less
> concentrated cell with the same inflation reads 2.1e-3. Repair costed in
> `docs/state/handoff-2026-10-02-176-variance-ratio-atom.md`: make the 13-carrier cap a SPLIT
> POINT (top 13 exact, the rest Gaussian), measured at 0.2% against an exact convolution on the
> real shape. Blast radius unmeasured; the per-job log line makes it countable without reruns.

## A note on counts

Counts of people, alleles or carriers below 20 are written `<20` throughout, per the AoU rule.
**One disclosure defect found and not fixed**: the region results file masks `MAC_case` at `<20`
while printing `MAC` and `MAC_control` in full, so the masked value is recoverable by
subtraction in every row. Wherever that table is emitted for export, the mask needs to cover the
combination, not the column.
