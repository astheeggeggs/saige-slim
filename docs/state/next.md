# What is worth doing next (2026-09-03, updated 2026-09-04)

**Start with `docs/state/handoff-2026-09-04.md`** -- it says what landed, which
gates passed with their numbers, where CI lives, and what is a fresh start.

Short list, kept deliberately short. `LEDGER.md` says what is broken;
`docs/state/*.md` says what was done and why; this file says what is NOT done
and is worth keeping in mind. It is not a plan and has no phases.

**Provenance is marked on every item, because it varies.** Two came from a
read-only assessment written by another session on 2026-09-03 and are recorded
here because they are expensive to rediscover -- neither is backed by a
measurement taken for it. Two came out of the S12a / LEDGER #50 work in the same
session. Nothing here has been profiled at realistic N, which is itself the
point: see "the gate" below.

## Resolved: the four known-wrong goldens WERE re-baselined (2026-09-04)

[Grounded; LEDGER #53/#56] Both step-0 sparse GRMs and both `binary_denseGRM`
step-1 fits were regenerated after the `geno2` fix and committed. Each was
proved byte-reproducible first, as CLAUDE.md requires: two independent
`make_goldens.sh` runs into separate candidate trees, `cmp` on all six files.
`compare_goldens.sh` is back to "all comparisons PASS" and the standard tier is
green again. Nothing else moved -- every step-2 golden and every sparse-GRM
step-1 golden was byte-identical throughout.

## The gate on all of it: profile at realistic N first -- NOW DONE

Every step-2 profile so far is N = 1,000, where `harness/bench/NOTES.md:352-359`
shows the serial writer -- not the work -- sets the ~2x ceiling at 4 threads.
That ranking is not expected to hold at biobank N, so choosing among the items
below before measuring is guessing. The sim fixtures are the instrument
(`harness/sim/`; the N=50k `bench50k` build lives outside the repo).

**MEASURED 2026-09-03, so this section is now evidence rather than a worry.**
See `harness/bench/NOTES.md` (the group-test section) and LEDGER #51. At
N = 50,000 over 210 genes: `--is_fastTest=TRUE` costs **+16.3%** and returns
BYTE-IDENTICAL output, because the region path never reads the flag. Escalation
rates run 2.9% (`y_null`, against a 5% nominal -- a clean negative control) to
14.8% on signal arms.

**An earlier draft of this file said an underpowered fixture would understate
group cost "by up to ~2x". That was wrong as stated** and is corrected here: the
2x applies per ESCALATED gene, but at a ~10% escalation rate the aggregate bias
is ~10%. Still worth reporting the rate alongside any per-gene cost -- it is
signal-dependent across ~5x -- but it is not the large effect claimed.

**The measurement moved the target, and then the profile moved it again.**
Swapping the PCG score test for `scoreTestFast` across a whole pass saves only
**8.3%**, so ~90% of gene-test time is NOT the per-marker score test. That was
written as an open question -- ultra-rare collapse, P1Mat/P2Mat assembly, or
R-side SKAT-O post-processing? **Answered 2026-09-04 by `Rprof` at N = 50,000**
(`harness/bench/NOTES.md`, "Where the OTHER 90% of a gene test goes"):

| bucket | share |
|---|---|
| C++ `mainRegionInCPP` | **43.9%** |
| R `region_postprocess` | **53.8%** |
| .. of which `SKAT_davies` under `integrate()` | **44.9%** |
| .. R glue | 5.0% |

So the largest single item in a SKAT-O group run is the SKAT package's Davies
routine under R adaptive quadrature, and **threading the C++ region loop caps at
1.49x on 4 threads** for SKAT-O. For BURDEN (`--r.corr=1`) the same fixture is
93.7% C++ and the same threading would give 3.36x. Replicates agree to 0.3 pp and
the BURDEN arm is the negative control -- it predicts `SKAT_davies` = 0 and
measures 0.0%.

**DONE (2026-09-04): the cells are evaluated in parallel.** `R/SkatoParallel.R`
drives `region_postprocess` three times -- collect the cells, evaluate them on a
persistent fork pool, serve the answers -- so its 175-line triple-nested loop did
not have to be restructured. `--nThreads` on a group run drives it (it used to be
discarded there). **1.30x at N = 50,000, byte-identical**, 8 workers no better
than 4. The remaining ceiling is 1.43x, so most of the available win is taken.

Still open on this path, in rough order of value:
* the C++ region loop is the other 43.9%; threading it would compose with this
  (serial fraction ~11% if both), but `qfc` thread-safety is unchecked and the
  SKAT-O side would have to be reimplemented to run in-process -- see the
  rejected option in NOTES.md;
* per-task overhead still costs the gap between 1.30x and 1.43x; batching cells
  into fewer, larger tasks would narrow it and has not been tried.

## From the 2026-09-04 silent-failure audit

A 26-agent fan-out hunted the bug class this branch keeps finding by accident
(silent failure that every golden reproduces). 15 raw findings -> 11 unique ->
10 adversarially verified by two skeptics each, both instructed to default to
REFUTED. Five survived; five were killed, including a plausible `raw_variant_ct`
claim and a `dosage_zerod` one.

**Two were reproduced by hand and are now FIXED (LEDGER #58, #59)** -- the
unbounded mmap `.bed` read (a `.bim` claiming one more marker than the `.bed`
holds gave byte-identical output at exit 0, while the pread path stopped), and a
region results file written with data rows and no header. Both have positive and
negative controls in the ledger entries, and both are gated by
`tests/testthat/test-step2-file-integrity.R`.

**Three more survived verification but are NOT independently reproduced.**
Recorded here rather than in the ledger, because the standard for a ledger entry
on this branch is evidence, and these have agent citations rather than a run:
* `src/SAIGE_fitGLMM_fast.cpp:101` (`BedByteSource::fill`) -- `seekg`/`read`
  status discarded in the GRM ingest, and failbit is sticky, so one short read
  would silently zero EVERY subsequent marker. Same family as #58; likely the
  same one-line size check fixes both.
* `R/SAIGE_Test_main.R:602` -- `maxMAF_in_groupTest` is extended by the
  MAC-derived cutoffs AFTER `g_cfg.maxMAFLimit` was fixed from the un-extended
  vector at `:337`, so a bin whose implied MAF exceeds the largest
  `--maxMAF_in_groupTest` would be emitted with all its markers already
  discarded. Conditional on `max(maxMAC)/(2N) > max(maxMAF_in_groupTest)`, i.e.
  small N or large MAC.
* `R/Util.R:38` (`updateChrStartEndIndexVec`) -- LOCO windows are `[min, max]`
  of each chromosome's marker indices, and the guard below tests only
  monotonicity of consecutive starts, so a handful of misplaced markers in an
  otherwise sorted `.bim` silently depletes a LOCO GRM. Over-exclusion, not
  contamination: the leave-one-out property survives.

Each needs the same treatment #58 and #59 got -- reproduce it or drop it.

## Outstanding DECISION for the user, not a task to pick up

[LEDGER #57, the part left open] `get_SKAT_pvalue`'s `try(..., silent = TRUE)`
(`R/SAIGE_SPATest_Region_Func.R:363`) turns any SKAT failure into
`Pvalue <- c(NA, NA, NA)` with an `error.code <- 2` that is assigned and never
read, printed or returned. The #57 fix broke the I/O path into it -- a failed
chunk-file read now stops the run instead of arriving here -- but a genuine
NUMERICAL failure (a singular or near-singular Phi, say) still turns a whole gene
into NA rows at exit 0 with no signal.

Three defensible answers and they are materially different operationally: warn
and keep the NA, abort the scan, or carry a count and print an end-of-run summary
line. Whichever is chosen changes what a long exome run does when one gene
misbehaves, so it was NOT defaulted. Ask before implementing.

## Items

**Format step-2 output rows on the workers.**
[Now MORE valuable than when this was written: a multi-model scan formats P rows
per marker on the serial writer, so the Amdahl serial fraction under
`--nThreads` rises with the number of models -- see
`docs/state/multimodel-step2.md`.]
[From the 2026-09-03 assessment; NOT independently measured.] `testOneMarker`
would build the row into a `std::string` on `Marker::MarkerResult`, leaving
`writeOutfile_single` (`src/Main.cpp:3072`, called at `:811` after the
`parallelFor`) as one `<<` per row on the main thread. The predicted effect is
on the serial fraction measured above (s = 0.305 at N = 1,000), so it should
shrink as the writer stops formatting. Byte-exactness rests on the
`std::ostringstream` carrying the same stream state as `g_out.single`
(`src/Main.cpp:89`, where no precision or flags are set) -- check that rather
than assume it, and gate on the marker goldens P1/P4/P7/P11/P12/P27-P33 plus
their PGEN mirrors.
`[ASSUMPTION: the prediction of a specific speedup is the assessment's
arithmetic from NOTES.md's Amdahl fit, not a measurement of this change.]`

**If phenotype batching is ever attempted, one trap is already known.**
[From the same assessment; recorded because it is cheap to write down and
expensive to hit.] Computing all P traits' score statistics as a single
`G' x [res_1..res_P]` matrix product CHANGES THE REDUCTION ORDER and therefore
breaks the byte gate; the per-trait dot products have to keep their serial
order. This is the same tension as the step-1 GetTrace GEMM decision, which was
resolved the same way. Nothing else about batching is designed, and the natural
unit for it -- a per-trait `TestContext` -- only became possible on the region
path with S12.

**`--is_fastTest` should warn on group runs.** [Measured; LEDGER #51] It is
inert for `--groupFile` and costs +16.3%. Making it work is bounded at 8% AND
would break byte-identity with every group golden and with upstream (a cheap
screen uses a different variance ratio and a different test, so genes near the
cutoff could fail to escalate). The defensible change is a warning, not an
implementation.

**`pval_cutoff_for_fastTest` is not reachable from the CLI.**
[Grounded, S12a] It is a `SPAGMMATtest()` parameter with no `--` option here or
upstream, so the fastTest mechanism cannot be tuned by anyone using the scripts,
and the escalated region path it controls cannot be reached by a golden --
which is why LEDGER #50 survived and why its regression test has to call
`SPAGMMATtest` directly (`tests/testthat/test-step2-region-escalation.R`).
Exposing it is additive and would let that path become a golden. Note the
CLI-option rules in `CLAUDE.md` (lint-clean, `option_list` metadata verified).

**Resolved: the psam sample map keyed on FID, not IID.** [LEDGER #52, FIXED
2026-09-04] `readPsamFile` now finds the `IID` column by name from the header,
strips `\r`, and refuses a duplicate IID instead of overwriting it. All 13 PGEN
goldens stayed byte-identical, exactly as predicted -- which is why the gate had
to be a fixture: `tests/testthat/test-step2-pgen-psam.R` derives psams at run
time (one FID shared by all 10,000 samples; `#IID SEX` with no FID at all; a
duplicated IID) and drives golden P13 through each. 6 failures at unpatched
HEAD, 0 after.

**Dosage pgens are hard-called for the GRM, and that is a real limitation, not
just a warning.** [Grounded, 2026-09-03] `genoClass`'s store is 2 bits per
sample, so `PgenByteSource` rounds any dosage to the nearest of {0,1,2} and
warns once with the count. For imputed data that discards information the GRM
and the variance ratios could have used. Carrying dosages properly means a
second storage path through `genoClass` and touches the GRM kernels, LOCO and
every step-1 golden -- deliberately out of scope, and worth measuring (does the
GRM actually move?) before being attempted.

## Deliberately not listed

The consolidation work (README is still the upstream stub, `pixi.lock` is
untracked, bench logs are in git) is real but needs a decision about who this
branch is for, which is not a decision to make from inside the code. The
remaining speedup candidates named in the older plan -- step-1 variance-ratio
estimation blocking, `setgeno` ingest -- are left where they are until the
profile above exists.
