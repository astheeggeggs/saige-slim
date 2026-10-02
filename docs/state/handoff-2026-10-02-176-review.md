# Handoff 2026-10-02 (later): independent review of LEDGER #176 -- mechanism confirmed, the proposed repair replaced

Read after `handoff-2026-10-02-176-variance-ratio-atom.md` (the account under review) and
`python3 harness/tools/ledger.py 176`. **Nothing here changes code, a golden or the ledger row.**
The row's STATUS still describes the earlier state; amending it is the owner's call.

Counts in this document are counts of CELLS, GENES or JOBS. Carrier counts are binned so that
nothing below 20 is written.

## Verdict in five lines

1. The mechanism in the earlier handoff is right. Read independently from source, before the
   handoff: `q = Tstat/sqrt(varRatio) + m1` (`src/SAIGE_test.cpp:1050`) is compared against a
   CGF built on the unscaled `gtilde` and `mu`.
2. Its proposed repair ("Fix 1": enumerate the heaviest 13, fold the rest into the Gaussian
   remainder) should not be built. It is anti-conservative by orders on shapes most real
   above-cap cells resemble.
3. The damage is a continuum in ONE quantity: the stretch divided by the sd of everything except
   the heaviest carrier. It is not a property of "concentrated" cells as the routing defines them.
4. On real data the gross damage is rare: seven unfiltered chromosome runs, about 732,000 cells,
   five wrong cells, one gene.
5. Proposal: gate on that quantity, exact convolution above the gate, variance ratio applied on
   the p-value scale. Four decisions are the owner's and are listed at the end.

## Part 1: the source, read before the handoff

- `q` is the observed sum `gtilde'y` with its deviation from the null mean `m1 = mu'gtilde`
  divided by `sqrt(varRatio)` (`var1 = var2 * varRatio`, `:507`).
- The CGF is the independent-Bernoulli one on unscaled `gtilde`/`mu`. The comparison is exactly
  `P(sqrt(r) (S - m1) >= T)`: the lattice is shrunk by `sqrt(r)` and the observation is not.
- The shift is `delta = |T| (1/sqrt(r) - 1)`. For `r > 1` the shift moves `q` into the gap BELOW
  the atom, so a one-atom cell gets essentially no relatedness correction at all (by structure;
  not measured).

## Where the earlier handoff is wrong or overstated

| claim | finding |
|---|---|
| "four independent routes agree on 6.0e-3" | one route: route 2 is the closed form of route 1, route 3 is a different statistic 1.6-1.8x away, route 4 is a bound |
| contrast cell 183_0: "`q` equals the observed statistic because the exact test ignores the variance ratio" | its `varRatio` is 0.999999, a different band; `q` is unstretched for that reason. Its reported saddlepoint p is 1.41e-2; the handoff's 8.7e-3 is the Monte Carlo value |
| contrast cell 183_130: "p 2.1e-3" | that is the Monte Carlo value; the run reported 9.4e-3 |
| "the `Concentrated cells above ... cap` line is in every v7 job log" | false: no `fxv7-` log prints it, the image predates the counter. Only jobs on `1.5.2-dev-7b7b447a` do |
| Fix 1 "measured at 0.2%" | measured on a kernel that does not exist in the repo (one-sided convolution on `gtilde` plus a remainder). `SKATExactBin_Work` is two-sided on Q, mid-p, on raw `G`, with no remainder. Fix 1 needs new kernel code either way |
| `fasttest="auto"` is what v4 used (`aou_spa_rerun.md`, `submit_abovecap`) | literally false, no effect: amr and eur records read `fasttest FALSE`, and on a group run the wrapper's `auto` passes no flag, so the commands are the same. Recipe wording corrected in the follow-up commit |

The mechanism survives all of these. The two contrast cells were offered as discriminating
evidence and each quoted the reference in place of the reported p.

## Evidence 1: Fix 1 against a worse adversary

`harness/tools/spa_176_fix1_adversary.py` (numpy; scratch venv). The earlier sweep used equal
secondary weights, one `mu`, and `q = 1.0`, which is not an attainable observed value. This one
uses low-`mu` secondaries (0.001-0.003) and `q` at attainable lattice points, `n_eff <= 4`, z > 2:

| true p | worst fix1/oracle | share off by >2x |
|---|---|---|
| 1e-2 to 1 | 0.52 | 0% |
| 1e-3 to 1e-2 | 0.14 | 33% |
| 1e-4 to 1e-3 | 0.028 | 69% |
| 1e-6 to 1e-4 | 4.5e-4 | 99% |
| 1e-8 to 1e-6 | 4.5e-4 | 100% |

Always anti-conservative (largest ratio 1.10). Every constant in the script is a round synthetic value. The earlier handoff's "decision-relevant is
p >= 1e-3" is the wrong band for a gene-level threshold near 2.5e-6.

The oracle itself (`spa_exact_ref.py`) held: 0 bracket violations in 1,200 comparisons against
`2^n` enumeration (signed weights, sub-step weights, `mu` 1e-4 to 0.5, both tails, remainder sd
0 / 1e-3 / 0.05). Two limits: with no remainder and `q` on an atom the bracket is the whole atom,
and its selftest reports bracket midpoints (12% wide on the inflated-`q` case, so the printed
1.019 is inside its own noise).

## Evidence 2: real data, seven unfiltered chromosome runs

Sampling jobs `fxrate-<anc>-<trait>-c<chr>` on image `1.5.2-dev-7b7b447a`, production FlexRV
flags copied from the `fxv4-` records, no region filter, outputs in
`<anc>/flexrv_caprate/`. They change no production output.

| run | cells | above cap | share | burden cells with `n_eff` <= 1.1, all genes | cells > 3 orders below their gene's ordinary median |
|---|---|---|---|---|---|
| amr EFRMB chr19 | 202,016 | 4,189 | 2.1% | 144 | 5, one gene (`ENSG00000213030`) |
| amr EFRMB chr20 | 76,593 | 2,081 | 2.7% | 67 | 0 |
| amr Asth chr20 | 85,353 | 1,640 | 1.9% | 66 | 0 |
| amr CAD chr20 | 85,359 | 1,596 | 1.9% | 66 | 0 |
| eur EFRMB chr20 | 94,278 | 409 | 0.43% | 6 | 0 |
| eur Asth chr20 | 94,387 | 243 | 0.26% | 3 | 0 |
| eur CAD chr20 | 94,388 | 251 | 0.27% | 3 | 0 |

- About half of each above-cap count is the ultra-rare pass, which the `BURDEN cells`
  denominator does not include, so the shares are overstated by about 2x for burden cells.
- chr19: 422 genes carry an above-cap cell; `n_eff` is 2-4 for 73% and <= 1.1 for 6.5%; 92% come
  from `sbeta5`/`sbeta10`/`sbeta20`. Material carriers: 34% below 20, 46% at 20-49, 20% at
  50-499, none at 500 or more, in any run.
- Controls: ordinary cells sit near nominal on chr19 (224 of 197,900 below 1e-3); the known gene
  is recovered; routed (exact-test) cells show zero extremes everywhere.
- Above-cap cells are RARER in eur than in amr. The prediction that they would grow with N was
  wrong (mine, and the row's withdrawn N-scaling story before it).
- Apart from the five, none of about 10,500 concentrated cells is below 1e-3 (about 10 nominal):
  conservative overall.
- **Blind spot**: this screen sees errors above three orders. A cell 5x or 50x out passes unseen.

Routed share of burden cells from the existing `fxv7-` logs, for scale: sas 36.4% (6 of 6 jobs),
eas 18.0% (15 of 15), afr about 1.5%, amr about 1.6%, eur about 0.6% (samples).

## Evidence 3: the dumpcheck of the known gene at a higher cutoff

`dumpcheck-amr-22`, 140 converged-below dumps, 78 distinct cells, replay equal to the reported p
in every one. Saddlepoint against the Monte Carlo reference AT THE SAME `q`:

| concentration of the weight (five groups, least to most) | cells | reported / MC |
|---|---|---|
| 1, least | 14 | 0.75-1.12 |
| 2 | 5 | 0.84-1.45 |
| 3 | 43 | 0.52-0.67 |
| 4 | 8 | 4.4-6.1 |
| 5, most | 4 | 1e-17 against nothing in 2e6 draws |

So the saddlepoint alone is good to about a factor of 2 on this gene. It is a hard configuration
for a saddlepoint, chosen because it is the broken one, so these ratios are not a rate. The
grouping is by how much of the weight the largest term holds; the values, and anything about
the carriers behind them, are deliberately not written here.

**Checked, and it did NOT fire**: dumps 0 and 60 ran on the DENSE kernel with a NEGATIVE
`NAsigma` (-9.7e-4, -2.8e-4), which looked like LEDGER #173's negative-remainder fallback.
The `fxcap-amr-efrmb` log (one gene, 144 cells, the run the dumps came from) finished, prints
`SPA fallbacks: 0 ...` and has NO `SPA remainder:` line, which is printed only when the count
is above zero. So CLAUDE.md's "never observed" stands. Why those two cells were dense is not
confirmed; the likely reason is that they sit in the low-MAC band that uses the sparse GRM
(their `varRatio` is 0.999999), where the dense kernel always runs and the dumped `NAsigma`
is computed but never read.

## Evidence 4: synthetic cells, with the real table as the control

`harness/tools/spa_176_stretch_synth.py`. Every input is a round SYNTHETIC value: one atom
(weight 1, fitted probability 0.01, set to fire), 81 secondaries of total weight T with fitted
probabilities drawn log-uniformly, a variance ratio of 0.98, three weight shapes.

Controls against Evidence 3: the well-spread regime reproduces in magnitude (saddlepoint over
exact 0.56-0.60 against the real 0.52-0.67); the moderately concentrated one in direction only
(1.6-2.2 against 4.4-6.1); the collapse in kind (8.5 sd, many orders). A first parametrisation
failed two of the three controls and was discarded. Read the SHAPE of the table, not its digits.

| T | `n_eff` | stretch / sd of the rest | cost of the stretch alone (truth / exact at the stretched `q`) | reported / truth |
|---|---|---|---|---|
| 0.003 | 1.01 | 8.5 | about 6e6 | about 1.6e-7 |
| 0.03 | 1.06 | 4.9-6.6 | 3e2 to 2e4 | 6e-5 to 3e-3 |
| 0.1 | 1.21 | 1.8-3.1 | 5-24 | 0.05-0.30 |
| 0.22 | 1.49 | 0.8-1.5 | 2.0-3.9 | 0.41-1.1 |
| 0.5 | 2.2 | 0.37-0.66 | 1.3-1.8 | 1.1-1.9 |
| 1.0 | 3.9 | 0.18-0.33 | 1.2-1.3 | 0.63-1.0 |
| 2.3 | 9-10 | 0.08-0.15 | 1.04-1.10 | 0.52-0.58 |

Two errors are in the last column: the stretch (fourth column) and the saddlepoint's own error at
a lumpy point, which is 0.56-3.0x here with or without the stretch. Below a ratio of about 0.7 the
stretch costs under 2x and the saddlepoint's own error is the larger term.

The variance ratio applied AFTER the exact tail, on the normal-score scale
(`p' = 2 Phibar(Phibar^-1(p/2) / sqrt(r))`), gives 0.93-0.94 of the unstretched answer in every
row and every shape. Synthetic only; uncalibrated under real relatedness.

## Evidence 5: the R-side saddlepoint has the same stretch, with a bounded effect

`get_newPhi_scaleFactor_traitType` (`R/SAIGE_SPATest_Region_Func.R`) sets
`q_adj = m1 + (q.sum - m1) * sqrt(var1 / sum(Phi))` on the RAW burden vector, deliberately
modelled on the C++ line (#104/#105). `harness/tools/spa_176_rside_rescale.R`, installed package:

| vector | VR | burden p it uses | Phi out / in |
|---|---|---|---|
| one atom, no stretch (control) | 1.00 | 3.2e-2 | 21.65 |
| one atom, stretched | 0.98 | 9.1e-24 | 15.22 |
| moderate secondary weight, stretched | 0.98 | 8.2e-5 | 15.15 |
| well spread, stretched | 0.98 | 3.7e-3 | 13.64 |

The burden p collapses the same way, but it only feeds `r = min(VarQ / VarQ_2, 1)`, and a
too-small p drives `r` to its cap. The consequence is a case-control correction that is NOT
applied (Phi 1.42x too small on the one-atom vector), not a p that is orders out. This path is
the SKAT-O arm and the conditioning markers; FlexRV's `--r.corr=1` cells do not use it. One
crude single-marker fixture; treat it as "same device, bounded", not as a measured size.

## Proposal (not implemented)

1. **Gate** where a binary test reaches the saddlepoint: `delta / sd_rest`, `sd_rest` being the
   sd of everything but the heaviest carrier. Below about 0.5 the saddlepoint runs as today.
2. **Kernel** above the gate: the convolution at the UNSTRETCHED score, replacing the
   above-cap give-up branch. The enumeration at or below 13 material carriers stays, so no
   routed cell moves.
3. **Variance ratio on the p-value scale** for cells that take the convolution.
4. **Counter and sidecar**, as every other route has.
5. **AoU**: a per-cell screen of the existing output in a dsub job (extreme p, far below the
   gene's ordinary cells), then rerun only the flagged genes. No full rerun.

Cost: the numpy reference takes 52 ms at 82 carriers and a 2e5-point grid, 323 ms at 500; the
worst case of every above-cap cell on amr chr19 is about 3.5 minutes, the realistic case seconds.
The grid cannot be coarsened for speed (bracket 7% wide at 2e4 points, 82 carriers).

Before building: the above-cap fixture (P77/P78) exists, so predict which goldens move first;
the change is under `src/`, so the full gate, one commit against row #176. The R-side function
needs the same decision.

## Decisions for the owner

1. **The threshold.** In the synthetic table the stretch alone costs about 1.6x at a ratio of 0.5
   and about 2x at 0.7. A fixture must exercise the chosen value before it is adopted.
2. **Sidedness and mid-p.** The existing exact test is a mid-p on a two-sided statistic; the
   convolution as written is a tail. They should agree by design.
3. **Whether the EXISTING exact test also gets the p-scale variance ratio.** It moves about a
   third of sas burden cells; its own row.
4. **Whether a factor-of-2 saddlepoint error on well-spread cells matters.** The gate does not
   address it, and it is a property of the method, not of #176.

## What is still unmeasured

- Errors between 3x and three orders on real data. Only the exact kernel beside the saddlepoint
  finds them.
- The exact reference at the unstretched `q` on the 78 real cells (needs the vectors; an R port
  of `spa_exact_ref.py` would run in the image).
- Any ancestry beyond amr and eur, any chromosome beyond 19 and 20, any trait beyond three.
- Why so few concentrated cells are hit. The account (the observed score must sit at the
  dominant term and the cell must pass the |Z| cutoff) is reasoning, not a measurement.

## Notebook notes (the workbench, not this repo)

- Fetch logs into the notebook VM sequentially and by sample; a thread pool is the shape of the
  earlier egress alerts.
- `replay_lines` breaks if a later cell reassigns the global `KEEP`; a helper named `_cat`
  shadows the one `fetch_log` uses. Both happened in this session.
- `next.md` was modified by another session and is NOT updated to point here.
