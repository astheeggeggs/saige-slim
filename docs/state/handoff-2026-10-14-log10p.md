# --is_output_log10P: log10 p-value columns, exact in log space (2026-10-14)

**Why.** On AoU the strongest genes print a p-value of 0: PNPLA3/ALT sits at
1.48e-323 on eur chr20 (the denormal floor; one more order and it is 0), and
on the BURDEN path the C++ writer already prints an underflowed p as the
log-derived string `"1.0E-400"` -- which the pooled-row Cauchy then parses
back as 0 (`parsePrintedPval`), so the gene's Cauchy row is 0 by CCT's
`is.zero` rule. A 0 is a statement about the representation, not the
evidence; the owner wants those genes estimated, not floored.

**What.** `--is_output_log10P` (default FALSE in the package; the AoU wrapper
passes TRUE via `LOG10P`) adds, right after `SE_Burden`:
- SKAT-O path: `log10_Pvalue`, `log10_Pvalue_Burden`, `log10_Pvalue_SKAT`
  (plus `log10_Pvalue_ACATV` / `log10_Pvalue_ACATO` at K > 1);
- BURDEN path: `log10_Pvalue_Burden`;
- the pooled Cauchy row, combined in log space (`CCT_log10` in R,
  `CCT_log10_cpp` in C++, twins to 2e-15).
Every pre-existing column and byte is unchanged with the flag off (the goldens)
and with it on (`test-step2-log10p.R` asserts column-for-column identity).

**Where the log is exact.** The SKAT-O saddlepoint tier (both layers:
`skatoSaddlepointSurvival(log = TRUE)` is `pnorm(..., log.p = TRUE)` on the
same Lugannani-Rice value, finite at q = 5000 where the linear p is 0); the
C++ score test's underflow spelling, read from its digits
(`log10FromPrinted`); the Cauchy combination (log-sum-exp over the 1/(pi p)
terms). Elsewhere it is log10(p) of a p above the floor by construction: the
Davies path only runs above the tail threshold, so the only way to a 0 with no
log is `--SKATO_tail_method=davies`, which gives NA there.

**A side finding.** `CCT()` computes the Cauchy tail as `1 - pcauchy(stat)`,
which cancels at large statistics: measured 2e-3 in log10 p against the
cancellation-free `pcauchy(lower.tail = FALSE)` on random vectors. The log
versions use `atan(1/t)/pi`, exact for t > 0. Not changed in `CCT()` (it
would move every pooled row); recorded here and asserted as the control in
the test (the reference is NEEDED).

**Gate.** Goldens P79/P80 (P2 and P58 with the flag), each proved
byte-reproducible twice; `test-step2-log10p.R` (53 checks: additivity and the
printed report line, column position and agreement with log10 of the printed
p on both paths including the pooled row, the two CCT implementations against
the reference with underflowed entries, the saddlepoint log twin, the digit
parser); the standard tier in both compare modes; `man/` regenerated
(`SPAGMMATtest.Rd`, new `CCT_log10.Rd`). `complementIndices` was deleted on
the way: #179's third site had been its last caller and the warnings ratchet
flagged it.

**AoU.** The wrapper's new `LOG10P` (default FALSE) passes the flag; the
production rerun sets it. A file with the columns cannot be SPLICED into a v4
file (the splice control refuses a header that differs), which is fine for a
complete rerun and is the reason not to turn it on for a partial one.

## Addendum 2026-10-15: the first AoU run with the flag lost a chromosome

IMG14 (`1.5.2-dev-524a299f`, the commit above) ran the SpliceAI v1-vs-v2 (0.40)
comparison with `LOG10P=TRUE`. `sai2v040-eur-alt-b5` delivered as SUCCESS with
`chrStatus` `failed 8`: 450 chr8 regions absent, the last region analysed GPT
(the ALT gene, p ~1e-146 on the other image). The log:

    Error in if (lin != 0 && log10S < 300) { : missing value where TRUE/FALSE needed
    Calls: do.call ... <Anonymous> -> get_CCT_log10pvalue -> CCT_log10

Two defects, one chain. (1) In `skatoSaddlepointTail` a per-rho saddlepoint p
that UNDERFLOWS to 0 is replaced by Davies' value (also 0) by the linear rule,
and the first log twin followed that branch, so `logpval` was `log(0) = -Inf`
and the pooled rules propagated it to `log.p.value` -- exactly where the column
exists to be finite. (2) `CCT_log10` had no rule for a -Inf input: its
log-sum-exp computed `Inf - Inf = NaN` and the `if` on it errored. The error is
raised in `region_postprocess`, outside the per-cell SKAT catch (LEDGER #57), so
it killed the chromosome's R process; the wrapper delivers the surviving
chromosomes and exits 0 (the deliberate fail-open, LEDGER #170's shape).

Fixes (one commit, gate green both comparators, every golden byte-identical
because no bundled cell reaches the tail): the twin keeps the saddlepoint's log
whenever it is finite and the saddlepoint did not FAIL (underflow is not
failure); `CCT_log10` and `CCT_log10_cpp` return -Inf on any -Inf input, which
is `CCT()`'s `is.zero` rule in log space, and check it BEFORE the `is.one` rule
as `CCT()` does; `log(area)` is clamped at 0 (a quadrature over a ~0 integrand
can return a tiny negative, which the linear sum absorbed and the log could
not) and a NaN combination falls to -Inf. Synthetic ladder measured before and
after: linear p identical at every scale; log10 p -17.2 / -67.6 / -268.1 at
scales 2 / 4 / 8 on both, then -Inf (before) against -602 / -1069 / -1670 /
-6677 / -26704 (after) at 12 / 16 / 20 / 40 / 80, scaling as the square of the
score as a burden statistic must. `test-step2-log10p.R` pins the ladder and
the -Inf rule in both twins.

Why the gate could not see it: no bundled cell is below 1e-5, so the tail's
log path ran end to end for the first time on AoU. The fixture bounds a
PHASE, not the regime (the repo's own trap, again). LEDGER #180 records a
related pre-existing divergence found while mirroring the rule: `CCT_cpp`
does not return on an exact zero as `CCT()` does.

Consequence for the AoU runs: every `LOG10P=TRUE` job on IMG14 must have its
`chrStatus` read; a failed chromosome is rerun on IMG15 (the fix commit).
Production had not started on IMG14.
