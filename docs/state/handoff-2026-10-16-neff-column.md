# --is_output_nEff: Kish's effective carrier count per cell (2026-10-16)

**Why.** Reading the v6 FlexRV results (runbook 19.8), three eur Height FlexRV-only hits at
1e-6..1e-8 were carried by sparse rungs (`MAC1`/`MAC3`/`MAC5`) with steep score weights and,
by Kish's count over the carriers' weights, rested on FEWER THAN FIVE effective carriers each,
while their plain counts (MAC >= 20, 20+ variants) showed nothing. On an inverse-normalised
quantitative trait one carrier contributes at most z ~ 4.9 at eur's N, so two or three tall
carriers of the top-scored singletons reach 1e-6 with no gene-level effect; binary traits are
protected by the exact-test routing of one-atom sets (LEDGER #172), quantitative ones are not.
The owner's decision: make the number a column, so the meta-analysis stage can filter on it.

**What.** `--is_output_nEff` (default FALSE in the package; `NEFF=TRUE` in the AoU wrapper)
adds `nEff_Burden` after `SE_Burden` (after the log10 columns when both are on) for every
cell and weight set on both region paths, NA on the pooled row. The value is
(sum x)^2 / sum x^2 over the per-sample weighted genotype sums the burden test saw, the
ultra-rare pseudo-marker included, taken from the carrier store once the ultra-rare pass has
completed every cell's vector (`NEff_GroupVec`, `src/Main.cpp`), exported to R for the SKAT-O
path and written by the BURDEN writer on the other; O(carriers x cells). `kishNEffInCPP` is
the arithmetic as an export, pinned to closed forms. The count is at most the number of
carriers and at most MAC; on an unweighted set (every carrier holds 1 or 2) it is at least
MAC/2, and exactly MAC when every carrier is a single het.

**Gate.** Goldens P81 (P2 + flag) and P82 (P58 + flag), each byte-reproducible twice; with the
flag off every golden is byte-identical; `test-step2-neff.R` (additivity and the report line
on both paths, the column position, the pooled NA, the bounds, the unweighted-set bounds and
the single-het identity, the closed forms); standard tier in both comparators.

**Reading it.** A FlexRV or custom-weight hit whose `nEff_Burden` is below ~10 rests on a
handful of people; the random-score null (runbook 19.9) sizes how often that reaches 1e-6
under no signal and is the number that sets the filter.
