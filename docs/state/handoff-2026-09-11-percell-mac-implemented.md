# Handoff 2026-09-11 (evening): #118, #126 and #127 implemented and gated

Read this first, then `docs/state/next.md`, then
`python3 harness/tools/ledger.py open`.
Previous: `handoff-2026-09-11-percell-mac.md`, which is the DECISIONS handoff
for this work and is still the place to read why each choice was made. This one
records what was built, what moved, and the two places where that handoff's own
predictions were wrong.

**Three ledger rows closed.** #118 (the per-set MAC), #126 (the min-group-MAC
gate on both region paths, default 1) and #127 (the pooled row's provenance).
Nothing was left half-done; what remains open is #114, #115, #116 and #117.

## What was built

**#118, `src/Main.cpp`.** `genoURMat_noweights` and `genoSumcount_noweight` are
both `q_anno_maf_weight` (jmr) now rather than `q_anno_maf` (jm). The five
non-collapsed accumulations (`MAC_GroupVec`, `MACCase_`, `MACControl_`,
`NumRare_` at :2604-:2618) are gated on `pr.w0(r) != 0`; `NumUltraRare_` at
:2690 is gated on the set's ultra-rare weight, at a separate site. The
pseudo-marker's contribution to `MAC_GroupVec` and to the case/control halves
comes from the unweighted twin (option b): `MAC_nw` is computed from the same
pre-`imputeGenoAndFlip` vector the weighted `MAC` is, with the same
`(1 - missingRate)` factor.

**Deliberately NOT changed**: the ultra-rare pseudo-marker's OWN single-variant
test context at :3085 still reads the WEIGHTED `MAC`, because that row tests the
weighted pseudo-marker. That is the reason every `*.singleAssoc.txt` came out
byte-identical, and it is the line to reconsider first if anyone ever wants the
UR row's variance ratio to follow the allele count instead.

**#126.** The default moved 5 -> 1 in `extdata/step2_SPAtests.R`,
`R/SAIGE_Test_main.R` and `man/SPAGMMATtest.Rd`. The SKAT-O half is enforced in
`region_postprocess`'s cell loop on `outList$MAC_GroupVec[jmr]`; a refused cell
emits no row and clears the carry-forward flag, exactly as `c.tested == false`
does on the burden path. The two gates read the same quantity computed two ways
(`MAC_GroupVec` against `genoSumcount_noweight`), which is what makes the
cross-path invariant sound.

**#127.** `R/PooledRow.R` plus a small record store in `src/Main.cpp`
(`getPooledRecordsInCPP`). Sidecar `<SAIGEOutputFile>.pooledTests.txt`, columns
`Region Model Path Rows Pvalues_combined Cells_refused_minGroupMAC Rule
Pvalue_pooled`, written only when a region has something to disclose. The
arithmetic is unchanged except one thing: `CCT_cpp` of no inputs returns 0.5,
and the BURDEN emitter writes `NA` instead, with the fastTest write decision
explicitly preserved so that only the printed number changes.

## The gate, and the two predictions that were wrong

`docs/state/prediction-118-126-127.md` was written BEFORE the gate and carries
the outcome appended after it. Read it rather than re-deriving any of this.

**Gate 1 (the new default).** 12 files move, 5 sidecars appear, nothing else --
every step-0, step-1, marker and plain-group golden byte-identical, every
`*.singleAssoc.txt` byte-identical, no row added or removed anywhere.

**The decisions handoff's "Prediction to beat" was wrong**, and the reason is
worth carrying: it said movement was confined to P58/P59/P61 "because neither
bundled custom-weight fixture contains an exact 0". That is the right test for
the COUNTER gating, but option (b) also changes the pseudo-marker's UNITS, and
that bites on every custom or transform set whose cell has an ultra-rare
collapse, zero weights or not. Eight more goldens moved for that reason, in
`MAC`/`MAC_case`/`MAC_control` only.

**The second wrong prediction was #126's own gate-5 matrix**, and it was wrong
BECAUSE #118 landed with it: it named the `REVEL lof 0.001` rows at MAC 3.34 and
3.55 as the cells the SKAT-O path would lose at gate 5 (and said 5 goldens where
it is 7). Those are weighted-unit values; under (b) they are 5, which a gate of
5 admits. At gate 5 the only files that move are the three FlexRV BURDEN
goldens, losing exactly 32 of 384 cells each, MAC 1 to 4, smallest survivor
exactly 5.

**Two checks worth knowing exist**, because they are cheap and they caught
nothing only because the code was right:
* every weight set of a cell must now print the PLAIN unweighted run's
  `MAC`/`MAC_case`/`MAC_control` -- the group files' `var`/`anno` lines are
  byte-identical, so this is one quantity computed two ways. 0 mismatches over
  8 files.
* inside P58, the sets whose weight vector is identically 1 on the cell
  (`REVEL_unw_unw` and `REVEL_unw_MAC{25,50,100,200,500}`) must keep their MAC
  exactly. They do -- 12 rows of 384, the negative control inside a moving file.

## Traps this work leaves behind

* **`--minGroupMAC_in_BurdenTest` at its default is inert by construction.** A
  cell that survives `indexNonZeroVec_arma.n_elem > 0` carries at least one
  allele, so its unweighted MAC is at least 1. Any test that wants to exercise
  the gate must NAME a value; `test-step2-region-threads.R`'s stale-record
  fixture was rewritten to pass `--minGroupMAC_in_BurdenTest=5` rather than
  inherit it, and it is the pattern to copy.
* **Never compare `Pvalue_Burden` across `--r.corr=1` and `--r.corr=0`.** Still
  true, and now easier to do by accident because the two paths share a gate.
  The one legitimate cross-path assertion is the SET of refused cells, compared
  on cell identity.
* **`region_postprocess` runs TWICE on the parallel SKAT-O path** (collect, then
  serve). Anything that RECORDS rather than computes must be suppressed on the
  collect drive -- `skatoCollect` passes `record_pooled = FALSE`. The same
  hazard is why `skatFailureRecord` is safe: the collect stub returns no `error`
  field. A third recorder added without this will double-count silently.
* **The FlexRV memory line grew by one column per cell** (4.0 -> 5.0 x N x
  cells), because the unweighted twin is now per set. It is still the formula,
  not a measurement.
* **Five goldens now carry a `.pooledTests.txt`** -- the quantitative SKAT-O
  runs, whose saturated cells send CCT down LEDGER #86's Bonferroni fallback. A
  NEW sidecar appearing in a candidate means a region started combining
  something differently, which is exactly what it is for.

## Still open

`#114` (dead preconditioner globals), `#115` (two help strings state the
opposite of their declared default), `#116` (a region analysed but producing no
row leaves a 0-byte file; mechanism NOT diagnosed) and `#117` (a weight of
exactly 0 inside a binary SKAT-O cell returns NA). #117 now has one more reason
to be taken: with #118 done, a zero weight is honest everywhere EXCEPT in which
markers reach the SKAT-O kernel -- `annoMAFIndicatorVec` is deliberately left
ungated, because that is #117's decision and not this one's.
