# Handoff, 2026-09-09: the 50k benchmark, and five ledger rows closed

Session summary. Eight commits, all on `strip`, NONE PUSHED. The gate
(`harness/prepush.sh`, both compare modes) is green at `0ddd0f85`.

```
0ddd0f85 LEDGER #113(a): a group file matching no marker stops instead of succeeding
a18c346f LEDGER #115: two more step-2 help strings contradict their own defaults
c1c4bf0b LEDGER #109: the dosage zeroing is documented, not gated
2c110894 LEDGER #108: --is_imputed_data refuses rather than trading a column for a constant
f067d49d LEDGER #107: --minInfo refuses a threshold it cannot apply
dcfab2d6 Benchmark strip vs upstream at N=50k, with a Firth-off arm
e779fda4 LEDGER #112: LOCO with a sparse GRM is refused, not silently overridden
04216280 LEDGER #110: the missing --chrom message names no test type
6f66e43b LEDGER #111: step 2's --LOCO help states its actual default
```

**Open ledger rows: 8 -> 3.** #114 (the owner's call, unchanged), and two
opened today: #115 and #116.

## 1. The benchmark

Full record: `harness/sim/upstream_compare/RESULTS-2026-09-09.md`. `RESULTS.md`
(2026-09-05) gained a supersession banner rather than being edited, because its
figures remain correct for the code they were taken on.

Headline, paired 3 reps, frozen models, 1 thread, upstream/strip on user CPU:
single-variant **1.78x binary / 1.99x quant** (was 0.94x / 1.28x), step 1
1.32x/1.34x and gene tests 1.38x/1.37x, both unchanged. The cause is `b11e6b07`
(workspace kernels, the `countOneMarker` QC pre-screen, closed-form Firth) and
the fact that step 1 and the region path are unchanged is the specificity check
on that attribution.

Firth is a p-value no-op, proved on printed bytes in both builds. All
strip-vs-upstream p-value differences reconcile; the gene-test cells the 09-05
tagger could not explain are #103-#106, pinned by re-running strip at
`cca6b497`, which reproduces the 09-05 figure of 0.018 to the printed digit.

New harness: `run_arm.sh` gained `VARIANT`; `compare_firth.R`, `bench_firth.sh`,
`run_pvalue_phase.sh` are new; `compare_all.sh` and `analyse.R` are
variant-aware.

## 2. The LOCO group (#110, #111, #112)

Mechanism chosen by the owner: **tri-state**. `--LOCO` defaults to NULL in the
step-1 CLI and in `fitNULLGLMM()`; absent resolves to `!useSparseGRMtoFitNULL`,
which is exactly what the silent override produced, so no existing command or
golden changes; an EXPLICIT `LOCO=TRUE` with a sparse GRM is refused.

Two traps this hit, both recorded in the rows:

- **The resolution must sit above `modelOut`.** `if (LOCO & isLowMemLOCO)` reads
  LOCO earlier than the sparse block, and NULL there is `logical(0)`, which `if`
  rejects.
- **Step 2's rewritten message must lead with the GENERAL condition.** The guard
  fires for a DENSE fit run with `--LOCO=FALSE` too (golden `binary_denseGRM`),
  so a draft that led with the sparse cause was #110's error committed again.
  The sibling test caught it.

**Two claims in #112 were WRONG and are corrected in the row.** (1) The override
was never silent: an unconditional `cat("Leave-one-chromosome-out is not
applied")` sits at line 203 of `harness/goldens/logs/step1_binary_sparseGRM_a.log`.
The row's evidence was `grep -i loco`, which cannot match that string. (2) "All
51 golden runs pass `--LOCO` explicitly" conflated the step-2 golden runs, which
do, with the step-1 ones, which mostly do not -- so the default combination had
been running on every sparse step-1 golden all along, which is why the refusal
could not break golden generation.

## 3. The dosages/INFO group (#107, #108, #109)

The owner's reasoning changed the repair on every row, and the reasoning is the
part worth carrying:

- **#107 `--minInfo`, #108 `--is_imputed_data`: REFUSED, not removed.** Upstream's
  behaviour is defensible where BGEN/VCF supply a real R2; this is a
  consequence of stripping those readers, not an upstream defect. So the
  options, their declarations and all four `t_isImputation` writer branches
  STAY, and a restored reader re-enables both features. The failure `--minInfo`
  actually had is worse than "inert": `--minInfo=0.8` reads `1 < 0.8` = FALSE
  and filters NOTHING, silently.
- **#109: the zeroing is KEPT and the DOCUMENTATION corrected.** An ultra-rare
  variant should not be recorded as fractionally present where it was merely
  missing, and the region path's ultra-rare collapse would accumulate those
  fractions. The counter-argument (mean imputation is unbiased under MCAR;
  zeroing biases allele counts down) is in the row, because this is a judgement
  and not a theorem.

**Gating #109 would have destroyed coverage, which is why it was not done.**
`is_imputed_data` defaults FALSE and no golden sets it, so gating would have
made P43 byte-identical to P46 and P49/P51 byte-identical to P50 -- three
goldens silently reduced to duplicates of two others, taking the only coverage
the dosage clean has. **Golden movement across the whole group: P48 RETIRED
(with `--is_imputed_data=TRUE` refused there is no output to compare), nothing
re-baselined.**

## 4. #113(a), and a third route found while fixing it

`SAIGE.getRegionList_new` now stops when the group-file/variant-file ID join
leaves zero markers. The condition is "no marker matched AT ALL"; individually
empty regions are normal and the positive control (20 real markers among 2
fakes) proves they still run.

**Unlike #76 the guard cannot prevent the empty file existing**: `openOutfile()`
runs at `R/SAIGE_SPATest_Region.R:952`, the group file is ingested at `:1364`,
and deleting the file would be unsafe because under
`--is_overwrite_output=FALSE` it can hold earlier chunks' real results. The
guard removes the SUCCESS, not the file.

**#116 opened**: a region that IS analysed but yields no qualifying cell also
leaves a zero-byte file at exit 0. Reproduced on the bundled fixture; the
MECHANISM IS NOT DIAGNOSED and the row says so -- do not quote a cause from it
without checking. The #113(a) control deliberately uses 20 markers, not 2, to
avoid that path, and says so in a comment so the avoidance is not mistaken for
coverage.

## 5. #115 opened: a sweep, and its method

Two more step-2 help strings contradict their declared default, the #111 shape:
`--impute_method` (declared `best_guess`, help says `[default=minor]`) and
`--is_output_markerList_in_groupTest` (declared `FALSE`, help says
`[default=TRUE]`). Recorded, not fixed in passing -- fixing in passing is how
#111's sibling stayed invisible.

**The sweep's method matters more than the two findings.** A naive string
comparison over all 130 `make_option` blocks reports ~30 mismatches and all but
these two are formatting (`1L` vs `1`, `1e-05` vs `1e-5`, quoting), plus
`--maxMAF_in_groupTest` and `--tauInit`, whose declared values contain commas a
lazy regex truncates and which are correct. Step 1's `--LOCO` also shows as a
mismatch and IS correct since #112. Compare values, not strings.

## Traps from this session

- **`pgrep -f <name>` matches your own watcher.** It bit three times: two
  waiters deadlocked on each other and idled 33 minutes, and later a gate was
  reported finished while it was still running, because `pgrep -f | head -1`
  returned the watcher's own shell. Wait on a PID captured from `ps`.
- **Timing on this box needs a load gate.** A run under load 33 took 2050 s
  wall against 198 s, at 20% CPU utilisation, and its USER CPU was itself
  inflated 2.2x. `bench_firth.sh` refuses above load 2.5;
  `harness/work/upstream_compare/await_quiet_then_bench.sh` waits for it.
  Spotlight indexing the outputs is a recurring tax -- `.metadata_never_index`
  markers are now in both `harness/work/` trees (gitignored).
- **The R suite loads the INSTALLED package**, so an R-only change needs
  `pixi run R CMD INSTALL .` before the tests see it. No header changed in this
  session, so the object clean was not needed; it IS needed after any `.hpp`.
- **`git add` aborts the whole staging on one bad pathspec.** A commit went in
  with only a `git rm`'d file in it and had to be amended. Check
  `git show --stat` after committing with explicit paths.

## NEXT UP

Three open rows, and the shape of each:

* **#114** -- the owner's call, unchanged from the last handoff: deleting the two
  unreachable C++ globals forecloses whatever the preconditioner was for.
* **#115** -- two one-line help fixes, whenever wanted. Worth doing as a sweep of
  all three CLIs rather than another single row.
* **#116** -- needs the mechanism diagnosed FIRST, then a decision: write a header
  even when no cell qualifies (fixes the downstream ambiguity for every route at
  once, but is schema-shaped and moves every golden), or refuse only when NO
  region in the run produced a row (closest to #76 and #113(a), moves nothing).

Still blocked on the owner: **FlexRV F0/F2 need the paper PDF** (bioRxiv blocks
fetches; the 16x12 transform table must be pinned verbatim).

Still open on performance: **batching SKAT-O cells into fewer, larger tasks**,
which closes the measured 1.30x -> 1.43x gap (`docs/state/region-perf-handoff.md`).

**Not mine, still in the working tree**: `harness/sim/README.md` is modified by
another session sharing this checkout and was deliberately left untouched.
