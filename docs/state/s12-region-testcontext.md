# S12: the region path gets a TestContext (handoff, 2026-09-03)

Branch `strip`, on top of `2d733b42`. The "Where to start next" of
`docs/state/s8-pgen-clone.md` named group-test threading and said the region
path needs S3's treatment first. This is that step -- plus a live defect that
fell out of doing it.

Two commits: `6b857ff5` (LEDGER #50) and `2914336c` (the refactor). In that
order deliberately, so the refactor's gate could be plain byte-identity rather
than "byte-identical except for a bug I am also fixing".

## The finding: LEDGER #50, the twin of #27

`m_flagSparseGRM_cur` was doing two unrelated jobs. R raises it between region
calls to mean *"this is the accurate re-run, DO write your output"*; the marker
and ultra-rare loops inside `mainRegionInCPP` also assign it, from each marker's
MAC, to mean *"use the sparse-GRM variance ratio for this marker"*. By the time
a write decision read it, the caller's intent was gone.

**#27 fixed that at the BURDEN write decision** (`Main.cpp:2240`, an
`isEscalatedPass` snapshot taken at function entry) **and left the non-BURDEN
one a hundred lines below still reading the member.** That is #50.

The single-variant-in-group rows have two routes to the output file and exactly
one fires per region:

| region | route |
|---|---|
| not escalated | C++ parks the rows in `.singleAssoc.txt_temp`; R copies them across (`copy_singleInGroup`) |
| escalated | C++ must write them directly, from the misread decision |

So an escalated region parked its rows too -- and the temp file is opened
`std::ofstream::out`, truncating, so the next region overwrote them. Measured on
the bundled binary fixture with golden P2's configuration:

| `pval_cutoff_for_fastTest` | regions escalated | `.singleAssoc.txt` |
|---|---|---|
| 1 (forced) | 2 | **header only -- all 19 rows lost**, exit 0, no warning |
| 0.05 (default, positive control) | 0 | 19 rows, byte-identical to golden P2 |

**No golden could have caught it, and the reason is worth keeping.**
`pval_cutoff_for_fastTest` is not a CLI option -- not here and not upstream --
and at its default neither bundled region escalates. So P2 (`--is_fastTest=TRUE`)
exercised only the copy route, P6 (`FALSE`) only the unconditional write, and
the escalated writer was gated by nothing. P2 and P6 are in fact **byte-identical
to each other in both files**, which says the flag they differ in reaches no
number on this fixture; the pair looks like coverage and is not.

The gate is `tests/testthat/test-step2-region-escalation.R`: a direct
`SPAGMMATtest` call in a child process (the C++ singletons are
one-dataset-per-process -- same reason as `test-step2-er-seeding.R`), which
**asserts the escalation actually happened -- 2 regions -- before it looks at
the rows**. Without that assertion a run that quietly failed to escalate would
pass as a green test, which is the LEDGER #30 shape. Drill-verified: restoring
the member read fails it two ways (0 rows, dims differ). The unescalated run is
the file's own negative control and still matches golden P2 in both compare
modes.

## The refactor

Mechanically what S3 did for the marker path:

- The three per-marker selection sites in `mainRegionInCPP` -- the rare-marker
  loop, the ultra-rare loop, the BURDEN collapse -- each build a local
  `SAIGE::TestContext` and hand it down. `sqrt(g_model->m_varRatioVal)` at the
  `P1Mat`/`P2Mat` rows becomes `sqrt(ctx.varRatioVal)`: the same marker's ratio,
  read from a local rather than off the model. The arithmetic is unchanged by
  construction because it goes through the pure selectors S3 had already
  extracted.
- `assign_conditionMarkers_factors` likewise.
- `SAIGEClass` **loses** `m_varRatioVal`, `m_flagSparseGRM_cur`,
  `m_isnoadjCov_cur`, `currentContext()`, both `assignVarianceRatio` overloads,
  both `assignSingleVarianceRatio` overloads,
  `assignSingleVarianceRatio_withinput`, `set_flagSparseGRM_cur` and
  `set_isnoadjCov_cur`. None of the three members was initialised by
  `setSAIGEobjInCPP`, so they also carried an uninitialised read for any caller
  that got the order wrong.

**The escalation signal deliberately did NOT become a TestContext field.**
Putting it back next to the per-marker sparse-GRM choice is exactly the
collision that #27 and #50 were. It is now `g_regionEscalatedPass` in
`Main.cpp`: written only by `set_flagSparseGRM_cur_SAIGE` / `_org` (R-facing
names and semantics unchanged), read once at the top of `mainRegionInCPP`, and
both write decisions read that. The comment on `struct TestContext` says so, so
the next person to want it there is warned.

Three dead things went with it, each verified unreachable rather than assumed:
`set_varianceRatio` (defined, declared, called by nothing); the
`assignSingleVarianceRatio` at the top of `mainRegionInCPP`, whose written
member no read in that function can reach because all three sites reassign
first; and the same call left behind in `mainMarkerInCPP` by S3.

**Gates**: 54 step-2 golden files and step0 byte-identical, five step-1 goldens
at `max_rel_err 0.000e+00`, standard tier **FAIL 0 | PASS 1839** in both compare
modes (+10 on 1829, exactly the new test file). For a pure refactor
byte-identity is the whole claim, which is why #50 was landed first.

## Where this leaves group-test threading

**The region path no longer writes to the model. That is one of the three
things the marker path needed, and the smallest.** The other two:

- **S4's shape** -- per-marker state owned rather than loop-hoisted. The region
  path still holds the two remaining `std::vector<bool>`s in `Main.cpp` and the
  parallel q-length vectors the marker path shed.
- **A parallel unit that pays.** This is the real obstacle and it is not in
  C++ at all: **the region loop lives in R.**
  `R/SAIGE_SPATest_Region.R:1108-1130` calls `run_mainRegionInCPP()` once per
  region, and the escalation decision round-trips through R --
  `region_postprocess` -> `cctpval` -> `set_flagSparseGRM_cur_SAIGE(TRUE)` ->
  re-run. So threading **across** regions (independent, seconds each, the unit
  that would actually pay) needs `region_postprocess`, the CCT and the SKAT-O
  assembly moved into C++; threading **within** a region (the loop this commit
  just cleaned) is a `parallelFor` over tens to a few hundred markers against a
  per-region fixed cost.

**So the honest next step is a measurement, not an implementation**: where does
a group run's wall time actually go -- the C++ marker loops, `region_postprocess`
in R, or per-region fixed cost? Nobody has profiled it. If it is dominated by
the R postprocess, C++ threading buys nothing and the work above is not worth
its risk.

Still open otherwise: #6 (deferred by choice), #14 (unverified), #18-#26,
#35-#39, #45, #46, #49. The PGEN per-handle arena trim still wants a large-N
fixture.
