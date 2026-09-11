# Handoff 2026-09-12 (small hours): the #126 null-calibration check, and #128

Read this first, then `docs/state/next.md`, then
`python3 harness/tools/ledger.py open`.
Previous: `handoff-2026-09-11-evening-117-115-114.md` (still the place for
what #118 cost and bought, and the disk clear-out).

## The job and its answer

LEDGER #126 moved `--minGroupMAC_in_BurdenTest` from 5 to 1, so burden now
tests cells at unweighted MAC 1-4 that were refused before, and those cells
take the sparse-GRM score test (`src/Main.cpp:3221`, the ternary
`(MAC > 20.5) ? false : flagSparseGRM` -- the handoff's abbreviated quote of
that line read as if the opposite). Nobody had checked they were calibrated.

**They are.** Prediction written before any run, outcome appended after:
`docs/state/prediction-126-null-calibration.md`. Quantitative null: ECDF(0.05)
0.059 / 0.057 / 0.054 / 0.050 at MAC 1 / 2 / 3 / 4 over 1,186 primary cells,
pooled 0.056 [0.043, 0.070], the MAC 5-20 reference 0.059, the same on the
SKAT-O path and on the 192-set FlexRV arm. Binary null (K = 0.01):
conservative exactly as discreteness predicts (pooled 0.023 [0.015, 0.033],
strata 0.020-0.025 against 1 - 0.99^m; never liberal). Gate 5 removes exactly
the sub-5 rows and nothing else moves in 22 run pairs. Positive control on the
63 causal genes: detections identical at both gates on burden and FlexRV
(14, 11), 16 -> 15 on the SKAT-O grid, the one loss being CCT dilution of a
gene 9% inside the threshold by three copies of one MAC-2 cell, within the
pre-registered bound. All on the #126 row.

**The fixture has 63 causal genes, not 21.** The job description carried the
2026-09-01 number; the 2026-09-10 rebuild (`sim_meta.json`) has 63.

## Three things that make the check non-vacuous, kept for next time

* **The gene-level file yields 11-16 sub-5 cells per stratum** (140 in all),
  a 95% interval of [0, 0.2] at alpha 0.05. The powered arm cuts each gene
  into 25-marker windows (1,614 regions, 2,249 sub-5 rows, 392 / 353 / 240 /
  201 primary cells at MAC 1-4); a null phenotype cannot see the grouping.
  `harness/sim/calib126/make_windowed.py`.
* **Predict the cell MAC with the collapse.** The engine's MAC is the sum over
  non-ultra-rare kept markers plus the per-sample MAX over the ultra-rare ones
  (#118's unweighted twin); the plain allele sum is off by the co-carrier
  count on 387 of 1,875 gene cells (by up to 80). Under the right rule: 0
  mismatches over 11,054 cells. `predict_collapsed.py` (needs numpy; the
  miniconda python has it, the system one does not).
* **One cell per (region, stratum) is the primary test**; the max-MAF ladder
  and the annotation masks nest, and FlexRV's 192 sets share supports.

## What it found that it was not looking for

**LEDGER #128 (new, OPEN).** A SKAT-O run dies at the first region with two or
more cells when the gate refuses every one of them: the pooled row of a region
with no surviving cell is built by the degenerate `rbind(NULL, cctVec)` branch
as a name-less character matrix, and `regionFinish`'s cross-region `rbind`
fails in `match.names`. One cell escapes (no pooled row is built), the BURDEN
path escapes (C++ pooled row, `NA`). Unreachable at the default gate 1; any
raised gate on a small gene reaches it. Reproduced 8 of 8 on the bundled
fixture; `tests/testthat/test-step2-region-allrefused.R` pins the crash and
both escapes. Bug-for-bug: not fixed; two cheap fixes on the row, neither
moving a golden.

**#116 re-observed**: two windowed regions are analysed and produce no cell;
the burden path writes an `NA` Cauchy row, the SKAT-O path writes nothing.
The two paths already disagree on an empty region; the decision on #116
should make them agree, and #128's fix belongs in the same place.

**The pooled Cauchy row at alpha 0.05 reads 0.06-0.08 on the quantitative
null in every arm at both gates**, on target at alpha 0.01. Not a gate
effect. [ASSUMPTION: the Cauchy combination of nested, positively dependent
cells is mildly liberal at a large alpha and accurate in the tail.] Untested
here; it is a `next.md` item, not a finding about #126.

## Traps this session added

* **`echo "exit=$?"` after a `$(date)` substitution reports date's status.**
  `harness/sim/flexrv_study/run.sh` had it (every `exit=` it ever logged was
  0) and is fixed; the calibration driver captures `rc` first. Combined with
  `run.sh`'s row count including the header, a dead run could look green
  twice over.
* **The Bash tool's shell is zsh, which does not word-split an unquoted
  variable.** A flag string in `$COMMON` reaches `Rscript` as ONE argument
  and the CLI reports the model file as missing (an empty path). Put a
  multi-flag driver in a `bash` script file.
* **The pinned library was at HEAD but nothing said so**: the scratch build
  from `git archive HEAD` had the same `SAIGE.so` md5. Building it was still
  right -- another session's tier was running on the pinned env, and
  `pixi run` can re-link it -- and the conda compilers can be activated by
  hand (`CONDA_PREFIX`, `PATH`, source `etc/conda/activate.d/*.sh`; drop
  `set -u`, the gfortran script reads an unset variable).
* **The #127 sidecar counts EVALUATED refusals**, and a refused cell's
  identical successor across the max-MAF ladder is sometimes re-evaluated and
  sometimes carried: 1 region in 53 disagreed with the distinct-cell count on
  the gene file, 0 of 1,029 on the windowed one. Bookkeeping, not a result.

## Where things are

`harness/sim/calib126/` (driver, analysis, group-file generator, the
collapse-rule predictor, the raw results). Runs and the scratch library were
in a session scratchpad and are gone; `README.md` there regenerates everything
in ~25 minutes of CPU at one thread. LEDGER: #116 and #128 open.
