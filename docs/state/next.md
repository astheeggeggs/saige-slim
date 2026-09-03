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

## The gate on all of it: profile at realistic N first

Every step-2 profile so far is N = 1,000, where `harness/bench/NOTES.md:352-359`
shows the serial writer -- not the work -- sets the ~2x ceiling at 4 threads.
That ranking is not expected to hold at biobank N, so choosing among the items
below before measuring is guessing. The sim fixtures are the instrument
(`harness/sim/`; the N=50k `bench50k` build lives outside the repo).

**A group-test profile MUST report the ESCALATION RATE, not just per-gene time.**
[Grounded, S12a] Under `--is_fastTest` a region whose Cauchy p falls below
`pval_cutoff_for_fastTest` runs `mainRegionInCPP` a SECOND time; the escalated
pass is numerically inert on the bundled fixture (P2 and P6 are byte-identical)
but costs a full extra pass. So per-gene cost is bimodal in a way that depends on
how much signal the fixture carries, and a fixture whose causal genes are
underpowered will understate it -- by up to ~2x on exactly the genes that decide
whether region threading is worth its risk. Quote per-gene cost split by
escalated/not, and the rate alongside.

## Items

**Format step-2 output rows on the workers.**
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

**`pval_cutoff_for_fastTest` is not reachable from the CLI.**
[Grounded, S12a] It is a `SPAGMMATtest()` parameter with no `--` option here or
upstream, so the fastTest mechanism cannot be tuned by anyone using the scripts,
and the escalated region path it controls cannot be reached by a golden --
which is why LEDGER #50 survived and why its regression test has to call
`SPAGMMATtest` directly (`tests/testthat/test-step2-region-escalation.R`).
Exposing it is additive and would let that path become a golden. Note the
CLI-option rules in `CLAUDE.md` (lint-clean, `option_list` metadata verified).

## Deliberately not listed

The consolidation work (README is still the upstream stub, `pixi.lock` is
untracked, bench logs are in git) is real but needs a decision about who this
branch is for, which is not a decision to make from inside the code. The
remaining speedup candidates named in the older plan -- step-1 variance-ratio
estimation blocking, `setgeno` ingest -- are left where they are until the
profile above exists.
