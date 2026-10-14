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
