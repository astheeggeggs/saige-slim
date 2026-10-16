# LEDGER #181: the region path's R-side saddlepoint returned NaN and killed chromosomes (2026-10-16)

**Found by** the v6 summary job's chromosome scan (section 19.5 of the runbook, moved into
`aou_flexrv_summary.sh` that morning after the notebook scan killed a session): afr HTN,
`sk6-afr-htn-b2` (chr 5 and 20) and `am6-afr-htn-b0/b1/b2` (chr 6, 12, 5), every other afr
job clean. The jobs delivered as SUCCESS. The failed chromosomes' logs (stderr, after
`*** chrN FAILED; full log follows ***`):

    Error in if (VarQ_2 == 0) { : missing value where TRUE/FALSE needed
    Calls: do.call ... get_newPhi_scaleFactor_traitType -> SPA_ER_kernel_related_Phiadj_fast_new

**Mechanism.** The region path's burden rescale calls the C++ saddlepoint kernel `SPA_pval`
from R and hands its value to `VarQ_2 <- Q_b / qchisq(p, log.p = TRUE)` and `if (VarQ_2 == 0)`.
That `if` is upstream's; a NaN reaching it is an R error outside LEDGER #57's per-cell catch,
so the chromosome's R process dies, the wrapper delivers the rest and exits 0. The NaN is new
since the v4 image (f2fb24da, 2026-09-23): the LEDGER #168 early returns (09-26) hand the
kernel's value back with no finiteness check, and the LEDGER #174 pinned guard (09-29) can
pass the kernel a vector with carriers removed. Which of the two, and on what vector, is
NOT yet known -- the fix dumps the vector so the reruns will say.

**The kernel half (added the same afternoon, after the pinned guard was ruled out: both HTN
models report `Pinned samples: 0`).** The NaN is born in the kernel: the #168 early returns in
`SPA_pval`, `SPA_ws` and `SPA_fast_ws` returned `getSaddle.pval` on `isSaddle` alone. All
three now require `std::isfinite` there and check the both-roots sum before returning it; a
non-finite value is a non-convergence (normal p, flag FALSE), which is what these cells got
before #168 and what v4 shipped. On the C++ test paths that flag is what triggers the #164
exact-test fallback, so there the behaviour is better than v4's. The R net below stays as
defence in depth and as the instrument. Warnings ratchet at baseline; `test-unit-spa.R` (129),
`test-step2-spa-dump.R` (53) and the new plant test (20) pass.

**What the value is used for.** Only the rescale `r = min(VarQ / VarQ_2, 1)` of Phi (the
case-control imbalance correction); the cell's `Pvalue`, `Pvalue_Burden`, `Pvalue_SKAT` come
from the SKAT kernel on the rescaled Phi. Falling back to the normal-approximation p for the
rescale is what `SPA_pval` itself already does on non-convergence (`pval = pval_noadj;
Isconverge = false`), so the fallback is upstream's own behaviour, now applied to the NaN
case too and DISCLOSED rather than silent: `Region saddlepoint NA: N ...` at the end of every
region run, `<SAIGEOutputFile>.spaRegionNA.txt` naming each cell (kind, region, group, MAF,
weight, carriers, pinned dropped, q, the normal log p used, the kernel's value, the dump
path). `is.na(VarQ_2)` is guarded as the second net. `R/SpaRegionNA.R`.

**Gate.** `test-step2-spa-region-na.R`: as shipped, golden P2 unchanged and the line at 0;
planted (`SAIGE_SPA_PLANT_NA_R=1`), the run completes, every rescale that reached the
saddlepoint is counted, named and dumped with the model's full vector, every cell's p is
finite. Standard tier, both comparators, every golden byte-identical (no bundled cell makes
the NaN). Lint clean.

**Open (on the row).** The cause, from the dumps of the four reruns (run with `SPA_DUMP=1` so
`SAIGE_SPA_DUMP_DIR` is set; the R dumps `spaRegionNA_*.tsv` land in the same tarball as the
C++ ones). And whether, for small cells, a rescale derived from the exact test's p should
replace the normal one -- a design change to make with the dumps in hand.

**The lesson, now three times.** #170 (a `stop()` on a missing label), the IMG14 log10 defect
(an `if` on NaN in `CCT_log10`), #181 (an `if` on NaN in the rescale): a value computed in R
outside the per-cell catch that can be NA costs a whole chromosome, silently. The CLAUDE.md
trap now says so in one line.

## Addendum, same afternoon: the instrument was blind, and why

The five HTN reruns on IMG16 all delivered with `failed none` and `Region saddlepoint NA: 0`.
The zero was NOT "no NaN": the kernel guard caught the non-finite value first and returned the
normal p with `isSPAConverge = false`, and R never sees that flag -- `SPA_pval` takes it by C++
reference but R passes it by value, so the `isSPAConverge <- TRUE` R sets beforehand is never
updated. The rescale's non-convergence has therefore been invisible to R for SAIGE's whole life,
and the net only fired on a value that reached R non-finite, which after the kernel fix nothing
does. The dumps never fired either, so the cause stayed unread.

Repair: `SPA_pval_status` (a new export over the same kernel) returns the value with the flag
and whether the non-finite guard fired; `get_newPhi_scaleFactor_traitType` calls it and records
every non-converged rescale with its reason (`nonconverged`, `nonfinite`, or
`nonfinite_reached_R` for the plant), dumping when `SAIGE_SPA_DUMP_DIR` is set. The line now
reads `Region saddlepoint NA: N burden rescale saddlepoint(s) did not converge ... K of them
returned a non-finite value`. Image IMG17; one HTN bin rerun with `SPA_DUMP=1` reads the cause.
