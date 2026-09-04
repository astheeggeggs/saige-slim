# What is worth doing next (2026-09-03)

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

## The gate on all of it: profile at realistic N first

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

**The measurement moved the target.** Swapping the PCG score test for
`scoreTestFast` across a whole pass saves only **8.3%**, so **~90% of gene-test
time is NOT the per-marker score test.** A profile that goes looking there will
find nothing. The unaccounted 90% is some mix of the ultra-rare collapse, the
P1Mat/P2Mat gene covariance assembly, and the R-side SKAT-O post-processing
(`R/SAIGE_SPATest_Region_Func.R`) -- and WHICH of those it is decides whether
threading the C++ region loop is worth anything, because if the time is in R
then it is not. That is the next measurement.

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

**The psam sample map keys on FID, not IID.** [Grounded, LEDGER #52; found
2026-09-03] `PgenShared::readPsamFile` (`src/PGEN.cpp:158`) does
`sampleInPgen[line_elements[0]] = idx;  // IID`, but column 0 of a
`#FID IID SEX` psam is FID. Step 2 would fail outright on any psam with
FID != IID, and silently undercount samples on one with repeated FIDs. Both
committed PGEN fixtures have FID == IID, which is the only reason all 13 PGEN
step-2 goldens pass. The steps-0/1 path added on 2026-09-03 does not go through
this map (it uses the identity-mapping `PgenShared` overload and counts psam
rows directly), so the bug is contained but still live for step 2. The fix is a
header-driven `IID` column lookup plus the `\r` stripping the pvar parser has
and this one lacks. **The gate is the missing fixture, not the fix**: every
existing golden will stay byte-identical, so a psam with FID != IID has to be
added or the fix proves nothing.

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
