# Handoff, 2026-09-09: identical region cells evaluated once

One change on `strip` after `ab047839`, R only (no C++ touched). Previous
handoff: `docs/state/handoff-2026-09-09-benchmark-and-ledger-sweep.md`.

## The question it answers

A region's cells are (annotation x max_MAF x weight set), and the annotation
list is a list of UNIONS (`lof`, `missense;lof`, `missense;lof;synonymous`).
On real data, and more so at low N, many cells hold exactly the same markers:
a union equals the smaller union in every gene without a marker of the extra
class under the cutoff, a max_MAF band that adds no marker repeats the band
below. Each such cell was evaluated again in full -- ~445 `SKAT_davies` calls
a cell, which is the 45% of a SKAT-O run measured in
`harness/bench/NOTES.md`. The C++ carry-forward (`q_maf_for_anno`,
`src/Main.cpp:3324`) covers only cutoffs ABOVE an annotation's largest MAF;
nothing covered an empty intermediate band or two annotations coinciding.

In the bundled fixture (golden P6, `plink_group_nofast.txt`) 4 of the 12
tested cells are such twins: `missense;lof;synonymous` == `missense;lof` in
both genes at both non-empty cutoffs.

## What landed

`R/SkatoParallel.R`, section "Identical cells evaluated once":

- `skatoJobKey(job)` -- exact key (`%a` hex of Score, Phi, r.corr, plus
  shapes and `regionTestType`); `skatoDedupe(jobs)` -- first-occurrence
  index per job, a key hit confirmed with `identical()` before it counts
  (the key is a pre-filter; the contract is R's own equality, attributes
  included). Keys go through `match`, never environment names (a 50-marker
  cell's key is ~60 KB, past R's symbol limit).
- `skatoEvaluate(cl, jobs, dedupe)` sends only the distinct jobs to the fork
  pool and fans the answers back out. The batch line now ends
  `(N distinct)`.
- `skatoMemoKernel()` -- a `skat_fun` memo for the drives that do not use the
  pool (`--nThreads=1`, and the BURDEN drive, where it is inert): one memo
  per `region_postprocess` drive, which is where the twins live.
- Per-run counters printed once before "Analysis done!":
  `SKAT-O cells: N asked for, M evaluated, K served from an identical cell.`
  Always printed, so a run where the memo never fired says 0.

Why it cannot move a byte: a cell is served another cell's answer only when
its inputs are `identical()`, and `get_SKAT_pvalue` is a pure function of
those inputs (the premise the fork pool already rests on). Whether the C++
produces bit-identical (Score, Phi) for identical marker sets is NOT part of
the argument -- if it did not, the memo would not fire. The duplicate rows are
still emitted, so the Cauchy row combines them exactly as before (they are
double-weighted there; upstream behaviour, deliberately preserved).

## The coverage finding

P6 predicts 4 served cells from its 4 twin rows; the run reports **3**. The
fourth pair (GENE2, `missense;lof` vs `missense;lof;synonymous` at max_MAF
0.01, a 4-marker cell) has identical Score and a Phi that differs in ONE
symmetric off-diagonal pair, entries (3,4)/(4,3), by 2 ULP
(`-0x1.8c593c4324454p+3` vs `...452p+3`). The golden row says the cell is
3 rare markers + 15 ultra-rare ones, so column 4 IS the collapsed
pseudo-marker and the differing entry is its cross-term with a real marker.
[ASSUMPTION: a P1Mat/P2Mat product whose summation order differs with the
cell's chunk layout; not traced.] The
printed digits agree, so the golden row pair is identical. Nothing to fix:
the memo's refusal is the negative control working. It does mean the served
fraction on real data is a measurement, not a count of equal rows.

## Gate

- `test-unit-skato-batch.R` +3 tests (dedupe pairs identical cells; a
  `(1 + 2 eps)` perturbation and an attribute-only difference each stay
  distinct; `skatoEvaluate` serves twins their answer and counts them; the
  memo kernel calls the kernel once per distinct cell, fresh per drive).
- `test-step2-skato-parallel.R` +2 tests: P6 at `--nThreads=1` and `=4`,
  golden bytes AND the printed counts `12 asked for, 9 evaluated, 3 served`,
  the batch line `(9 distinct)` -- a memo that never fired cannot pass.
- All region goldens exercise the memo (every one has twin rows), so the
  standard tier's byte gate is the positive control that it fired without
  changing anything.
- `harness/prepush.sh` at this change on `ab047839`: FAIL 0 in both modes
  (byte 3170 / numeric 3370 passes).

## The 50k measurement: design and predictions (written BEFORE the run)

Driver: `harness/work/region_profile/drivers/region_s6_ab.sh` (gitignored
dir, like the earlier drivers). Paired A/B on the sim50k fixture, `y_binary`
SKAT-O with `harness/sim/04_run_saige.sh`'s flags and frozen models: a
scratch install of `ab047839` (the parent; same `.o` files, so `SAIGE.so`
md5 must match the installed one) against the installed build, arms
alternated per rep, 3 reps, `--nThreads=1` first and `4` second. Every run
refuses to start at load > 2.5 and records `lib=`/`so_md5=` (the null
control: which build actually loaded). Outputs cmp'd a-vs-b per config.

The fixture is dense (every gene carries all three classes, ~1,000 markers a
gene), so its own twin fraction is expected to be ~0 -- that arm is the
negative control for the memo's overhead, not a measurement of the saving.
Two more configs, through two env hooks added to `04_run_saige.sh`
(`GROUP_ANNO`, `GROUP_FILE`; defaults unchanged):

| config | what | twin fraction (predicted) | new/old wall, 1 thread | 4 threads |
|---|---|---|---|---|
| native | driver's 5 groups x 3 MAF | ~0 | 1.00 +- 0.02 | 1.00 +- 0.02 |
| planted | the 5 groups plus each unioned with a label no marker carries (`:zz`) -- 5 of 10 groups twins by construction | ~0.5 (less the ULP misses) | 0.67-0.75 (cells double, so Davies is ~2/3 of the old run and half of it goes) | ~0.85 |
| thinned | `high` relabelled `medium` in a random half of the genes (seed 20260909; `high:medium` == `medium` there) -- the real-data mechanism, no refit | ~0.11 (105 genes x 3 twin cells of ~2,835) | 0.94-0.96 | ~0.97 |

Basis: R-side SKAT-O is 45-54% of a 1-thread `y_binary` run
(`harness/bench/NOTES.md`), the pool absorbs about half of that at 4
threads, and a served cell costs one key + one `identical()`. A result far
from these is a bug to locate (units, arm swap, the wrong library loaded),
not a discovery. The twin fraction is read from the `SKAT-O cells:` line of
each log, so it is measured, not assumed.

## The 50k measurement: results (2026-09-09 19:21-21:00)

Log `harness/work/region_profile/ab_s6.log`, TSV
`harness/work/region_profile/results_s6_ab.tsv`, outputs under
`harness/work/upstream_compare/strip/out/s2group_y_binary_r[ab]*_frozen_{native,planted,thinned}.txt`.
Both libraries' `SAIGE.so` md5 `28d522d1...` (equal); the 18 old-arm logs
all record the `oldlib` path and the 18 new-arm logs none of them (the
which-library null control). **Every a/b output pair is byte-identical:
54 of 54 files** (region rows, singleAssoc, markerList). Spotlight held
the load at 4-12 throughout (the 2.5 rule was waived deliberately: paired,
alternated arms; ratios only).

| config | twin fraction, measured | new/old wall, 1 thread (3 reps) | predicted | 4 threads (3 reps) |
|---|---|---|---|---|
| native | 14 / 2527 = 0.6% | 0.999, 1.000, 1.005 -> **1.00** | 1.00 +- 0.02 | 0.99, 1.09, 1.07 -> 1.05 |
| planted | 2391 / 5054 = 47% | 0.709, 0.729, 0.684 -> **0.71** | 0.67-0.75 | 0.91, 0.84, 0.85 -> 0.87 |
| thinned | 236 / 2411 = 9.8% | 0.952, 0.948, 0.954 -> **0.95** | 0.94-0.96 | 0.87, 0.98, 0.92 -> 0.92 |

Every 1-thread ratio is inside its pre-registered band; the served counts
are identical across reps (deterministic). Of the 2527 planted twins,
136 (5.4%) were not bit-identical -- the same ULP class as P6's fourth
pair -- so the served fraction on real data is (twin rows) x ~0.95.

**Rule of thumb from this:** saving = twin fraction x ~0.5 at one thread
(the Davies share), about half that in wall at 4 threads where the pool
absorbs part of it. The memo's own cost: 1-thread user time equal to
within 0.4 s on 2527 cells; a direct in-process timer put `skatoDedupe`
at 0.78 s per 4-thread run over 24 MB of keys.

**The 4-thread numbers are noisy at +-10% on this box and only the
1-thread ratios should be quoted.** native's 1.05 at 4 threads is not
overhead: rep 1 was 0.99, and the thinned OLD arm alone swung 178-198 s
of user CPU between its own reps. Four workers plus four TBB threads on
8 cores under Spotlight's ~2 cores is the reason.

## Not done

- The C++ per-cell work (ultra-rare collapse, pseudo-marker test,
  `genoSumMat` residualisation, BURDEN emitter: ~25-30% of a BURDEN run's
  C++) is still per cell. Deduping it needs set-equality detection before
  the collapse and column sharing, and its byte-identity would rest on
  iteration-order reasoning rather than value equality. Worth it only if a
  BURDEN profile on real data shows a high twin fraction; measure first.
- A quiet-box 4-thread measurement, if anyone needs the 4-thread ratio
  to better than +-10%.
