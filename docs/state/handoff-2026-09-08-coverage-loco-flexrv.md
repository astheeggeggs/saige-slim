# Handoff: three coverage holes closed, eight decisions opened, FlexRV started (2026-09-08, third session of the day)

Branch end commit: `fd54eb04`. Eleven commits, `3bb270c5..fd54eb04`, each gated
by `bash harness/prepush.sh` green in BOTH compare modes before it landed. The
two earlier sessions of the day are `handoff-2026-09-08-burden-variance-fixes.md`
(#103/#104/#105/#106) and `handoff-2026-09-08-ledger-zero-and-coverage.md`; this
file supersedes both as the pointer.

`python3 harness/tools/ledger.py lint` reports **8 open of 114, 4 decided**.
Gate at the end: **FAIL 0, PASS 3110 byte / 3310 numeric** (from 2972 at the
start of the day).

## NEXT UP, and it is a decision, not a job

**Eight open ledger rows, and every one of them is yours to decide.** This branch
records rather than fixes, and coverage work has now produced more decisions than
it has questions. They group into three commits:

* **#110, #111, #112 -- the LOCO surface.** One subject, blocked only on
  choosing #112's mechanism. #112 is the one with blast radius: `--LOCO` defaults
  to TRUE in BOTH CLIs, so a hard error on `useSparseGRMtoFitNULL & LOCO` rejects
  the DEFAULT invocation of the recommended workflow until every script adds
  `--LOCO=FALSE`. The alternatives are a GRM-dependent default for `--LOCO`, or a
  loud warning. **The owner has already decided the DIRECTION** (incompatible
  combinations should be refused, not silently reconciled); what is open is which
  mechanism, i.e. who breaks. #110 and #111 are one-line message fixes that ride
  along.
* **#107, #108, #109 -- missingness and dosages.** `--minInfo` is a switch at 1;
  `--is_imputed_data` REPLACES `MissingRate` with a constant rather than adding a
  column; `--dosage_zerod_cutoff` silently turns `--impute_method=mean` into
  `minor` below the MAC cutoff.
* **#113(a) and #114.** #113's half (b) is already done. #114's deletion
  forecloses whatever the preconditioner was for, so it is a decision.

**FlexRV F0 and F2 are blocked on you supplying the PDF** (bioRxiv blocks
fetches; the 16x12 transform table must be pinned verbatim). F1 is landed.

## What was done

**Two calibration studies at N = 50,000** (`docs/state/null50k-2026-09-08.md`,
`null50k-binary-2026-09-08.md`). The null is calibrated at alpha = 0.05 and
conservative at alpha = 0.01 on both traits; the day's burden-variance fixes are
correctly signed and calibration-neutral at this scale. Predictions were
pre-registered before any result was read and held, including the one that said
lambda_GC could NOT resolve the effect and it had to be a paired per-row ratio.

**Three coverage holes closed**, each of which had "no golden AND no test":
missing calls (P42-P48), dosages (P49-P51), LOCO end to end (P52-P55) plus the
quantitative region shapes (P56-P57).

**Eight ledger rows opened**, #107-#114, all from those three jobs.

**FlexRV F1 landed** and the plan re-resolved.

## What each job cost, and the traps worth carrying

* **A confounded baseline looks like a clean one.** The handoff before this
  called the 2026-09-01 `step2_group_y_null.txt` "a genuine BEFORE baseline"; ~80
  commits separate it from HEAD, and the diff showed it at once -- 255 extra
  singleAssoc rows in six whole gene blocks, which are exactly the six genes that
  escalate under `--is_fastTest`, i.e. LEDGER #50 fixed 2026-09-03. Attribution
  needed scratch installs at the actual parent commits. Three exist and are
  reusable under `harness/work/null50k-20260908/`.
* **`lambda_GC` is invalid over a point mass.** The SKAT-O columns carry
  `p == 1` exactly (26.7 % of the Cauchy rows; LEDGER #86). Excluding the atom
  reports lambda 1.47 on an acceptable null; including it reports 0.73. Use the
  ECDF at alpha with a binomial interval.
* **Size a QC fixture against the ANALYSED subset.** `missingRate` is
  `numMissing / N` over the model's samples (fam rows 1..1000), not the .fam's
  10,000. A first fixture sized against 10,000 made SAIGE report 0.032/0.096/0.161
  where 0.05/0.10/0.15 were meant, putting the boundary marker on the wrong side
  of `--maxMissing`.
* **Check what a change REMOVED, not only what it added.** `--is_imputed_data`
  looked like it ADDED a column; it REPLACES one. Caught by a FILE SIZE -- the
  golden came out smaller than its twin, and an added column cannot shrink a
  file. The check made was `setdiff(new, old)` and never `setdiff(old, new)`.
* **A comparator can silently compare nothing.** `read_saige_table` returns
  character columns by design; a helper gating on `is.numeric` compared no column
  and reported "no differences" for every pair. Three tests passed vacuously.
* **So can a GOLDEN comparator.** `compare_step1.R` descended into a list only
  when it was named `obj.noK`, so `LOCOResult` was skipped SILENTLY -- a LOCO
  step-1 golden would have gated theta and the base fields while comparing not one
  value of the thing it exists for. Fixed and drilled.
* **Measure before writing a golden, not after.** The quantitative `fastTest`
  region run is BYTE-IDENTICAL to P38 (the region path never reads the flag,
  #51). Written first and checked after, the repo would have gained a file that
  looks like coverage and is none. It became a test instead.
* **A plan's line numbers rot fast.** Every one of the 18 references in
  `flexrv-plan.md` was stale within a day, by up to ~145 lines. They are now
  SYMBOL-anchored with the number as a dated hint.

## Still open, in the order I would take it

1. The eight decisions above.
2. FlexRV F2, once the PDF exists. F3 is the commit to review: it changes
   `mainRegionInCPP`'s signature, which pulls in `RcppExports`, the stale
   `Main.hpp`, `src/Makevars` and the delete-`src/*.o` rule. F5's memory number
   (~2 GB at N = 500k, extrapolated) is a genuine go/no-go.
3. Smaller coverage: PGEN mirrors of P56/P57, quantitative BURDEN twins,
   dosages AND missingness together, the 14 untested step-1 numerical knobs.
4. The binary tail at K = 0.01 (6 of 210 Cauchy rows at p <= 0.01 against 2.1
   expected, one-sided p = 0.020, and IDENTICAL before the fixes -- so not caused
   by them). Settling it needs replicate phenotypes, and it is a question about
   SAIGE-GENE+ rather than about this branch.
5. `covr` baseline is stale (R 55.2 %, step-2 C++ 72.2 % line / 42.9 % branch, at
   C8). Today added a lot of test surface, so it is wrong in an unknown direction.
6. Halving the pre-push gate -- and note the modes are NOT redundant: numeric
   carries 200 expectations byte does not, so it means classifying 95 test files.

Scratch material, all untracked: `harness/work/null50k-20260908/` (the three
scratch installs and every arm), `harness/work/loco-20260908/`,
`harness/work/missing-20260908/`, `harness/work/quantregion/`.
