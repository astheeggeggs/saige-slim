# Handoff 2026-09-11: the per-set MAC (#118), the SKAT-O gate (#126), the pooled row (#127)

Read this first, then `docs/state/next.md`, then `python3 harness/tools/ledger.py open`.
Previous: `handoff-2026-09-10-flexrv-sweep-124.md` (the FlexRV sweep and the
fitter rows #121-#125, still correct).

**No source file changed this session.** Everything below is decisions,
measurements and ledger. The implementation is fully specified and NOT started.

## First, the tree was committed and pushed

The working tree the previous sessions left uncommitted is in: LEDGER #119 and
#120 fixed (`c8589476`, `d4c073e2`), the FlexRV end-to-end suite and #118's
characterization (`3206f359`), saige-ae's skato/docdrift/CLAUDE.md unit
(`b750458a`), the sim study and its S2-S4 sweep (`0b53d7d7`), plus the
scaffolding move and a `.gitignore` for `flexrv.pdf`. Pushed as
`a7d8af30..606bb3f4`; the pre-push gate ran green in both compare modes at
`606bb3f4` and CI run 34547342979 followed.

## The four decisions, all the owner's

1. **#118: the per-set MAC, EVERYWHERE** -- variance-ratio band, ER
   eligibility, the min-group-MAC gate and the printed columns. Not the
   printed column alone, because the gate consequence is statistical.
2. **#118 units: FULLY UNWEIGHTED (option b).** `MAC` becomes a genuine
   allele count, so a gate threshold of 1 or 5 literally means alleles.
3. **#126: ONE flag, default 1, applied to BOTH paths.** The owner's reason is
   the load-bearing one: the SKAT-O path does not merely contain a burden
   component, it REPORTS `Pvalue_Burden`, `Pvalue_SKAT` and `Pvalue`, so a
   gate on burden testing that skips the path printing `Pvalue_Burden` is
   incoherent. Consequence: burden becomes more permissive on real data
   (refusing below MAC 1 where it refuses below 5 today); no bundled fixture
   can show that, since every committed burden golden is byte-identical at
   gate 0, 1 and 5.
4. **#127: a sidecar, disclosure only.** The pooled row's arithmetic is the
   paper's own and does not change.

## The measurements that drove them

- **The burden STATISTIC is already correct under a zero weight.** Weighting a
  marker 0 gives output identical to deleting it from the group file, to every
  printed digit, on both the non-collapsed arm (BETA -0.0654086, Pb
  2.490047E-01) and the ultra-rare arm (Pb 2.556854E-01). Positive control:
  the marker at its real weight differs. So #118 is a REPORTING AND GATING
  defect, statistical only where the wrong MAC crosses a band or the gate.
- **#117 under SKAT-O is worse than "one NA row".** The failed set drops
  silently out of the Cauchy combination and the gene reports the surviving
  set's p verbatim (0.6787 against 0.5319 for the control). GENE2, with no
  zero weight, combines normally -- that is the control. The fix target is
  therefore exact and needs NO new golden: zeroing must equal deleting, the
  identity burden already satisfies.
- **Impact of #118 over P58's 384 sets:** 45.3% of printed MACs count markers
  the set zeroed (median 2.83x); 20.6% would change variance-ratio band, but
  this model's two `null` ratios are 0.9475 and 0.9421 so the SE moves 0.29%;
  8.3% pass the gate today and would be refused. Worst case
  `REVEL_lof_MAC1` in GENE2: printed 37, true 1.
- **The two paths do NOT share a burden estimator.** Same cells: 1.010x at MAC
  336.4, 1.002x at 74.73, **2.76x at MAC 2**. No tolerance spans that.
- **Cell-level MAC is identical across paths**, and the gene-level (Cauchy)
  row carries no MAC at all in either. So the gate will refuse the same cells
  in both paths -- which is what makes #126's cross-path invariant sound.

## Three corrections to claims committed EARLIER THE SAME DAY

Each was committed wrongly first and is corrected in the history rather than
silently overwritten. Read them before trusting an older sentence.

1. **The SKAT-O served-cell discriminator is NOT the runner image** (`606bb3f4`
   corrects `b750458a`). Image 20260907.300.1 served 4 in run 34403657864 and
   3 in run 34538797394. The cause is UNIDENTIFIED and no doc may name one.
2. **The dosage-zeroing refutation was testing the wrong object.** It is
   #109's zeroing acting on the weighted PSEUDO-marker, not on the
   non-collapsed markers; the pseudo-marker's own MAC is what must be at or
   below `dosage_zerod_MAC_cutoff`. Perfect necessary condition: of 318 sets
   with pseudo MAC > 10, ZERO show the residual; of 66 at or below, 39 do.
3. **#118 is NOT "upstream's semantics"**, which CLAUDE.md and the row both
   said. Upstream 365deeff has NO per-set dimension: `jmr` appears zero times
   in its `src/Main.cpp` and all three group counters are `q_anno_maf`, per
   CELL, and mutually consistent. Phase 5 of THIS branch widened them to
   `q_anno_maf_weight` without making the accumulations per-set. **That
   removes the bug-for-bug objection**: the rule protects upstream behaviour
   and there is none here to protect. Guarded by docdrift rule
   `flexrv-percell-mac-upstream`.

## Primary source: SAIGE-GENE+ (PMC9534766, open access)

The owner supplied it; four rows were re-grounded on it (`c27b16b4`).
- `"Ultra-rare variants with MAC <= 10 are collapsed to a single marker"` and
  the per-sample **maximum** are the published method, so `std::max` at
  :2664/:2686 is not a quirk.
- User weights are applied **before** the collapse -- our custom-set order.
- **No minimum-MAC filter appears anywhere in the paper.** Timeline: received
  2021-07-11, `minGroupMAC_in_BurdenTest` first appears upstream 2022-03-12
  (v0.99.3), accepted 2022-07-29. Added DURING review, Methods never updated.
  Undescribed, not outdated. NOT checked: the online docs and release notes.
- The Bonferroni fallback of #86 is **the paper's own prescription**, which
  strengthens the earlier decision to leave it and shrinks #127 to disclosure.

## Implementation shape (line numbers resolved 2026-09-11, nothing written)

- Widen `genoURMat_noweights` (`src/Main.cpp:2069`, accumulated :2687) from
  `q_anno_maf` to `q_anno_maf_weight`, gating its accumulation on
  `w0(r) != 0`. This is the unweighted twin upstream already maintains.
- Take the pseudo-marker contribution to `MAC_GroupVec` (:3147, :3212) and to
  `MACCase_`/`MACControl_GroupVec` (:3155-:3157 and its twin site) from that
  unweighted matrix instead of the weighted `genoURVec`.
- Gate the five non-collapsed accumulations on `w0(r) != 0`: `MAC_GroupVec`
  (:2604), `MACCase_` (:2606), `MACControl_` (:2607), `NumRare_` (:2617), and
  `NumUltraRare_` (:2690, a SEPARATE site).
- Make `genoSumcount_noweight` per-set and weight-gated so the gate, the
  variance-ratio band and ER read it (:3470 is where the gate reads it today).
- `pr.w0(r)` is already in scope two lines below :2604.

## First moves

1. Scratch build (`R_LIBS_USER=<scratch>`), not the pinned env.
2. **Predict per-file golden movement BEFORE running the gate**, in a file, as
   the #124 re-baseline did (`prediction-regen.txt` is the method worth
   copying). Prediction to beat: movement is confined to P58/P59/P61; every
   SKAT-O golden and every non-FlexRV golden stays byte-identical, because
   neither bundled custom-weight fixture contains an exact 0 and P60 has zero
   NA rows (which by #117 proves no in-cell zero weight).
3. Run at **gate 1**, then at **gate 5** (#126's two-run matrix).

## Traps specific to this work

- **Gate 1 does NOT mean nothing moves.** It neutralises only the gate;
  #118's per-set MAC still moves 174 MAC values and 79 bands.
- **Never compare `Pvalue_Burden` across `--r.corr=1` and `--r.corr=0`.** See
  the 2.76x above. The ONE legitimate cross-path assertion is that the SET of
  refused cells agrees, compared on cell identity only.
- **The gate on the burden path wraps more than the test** (:3480 also guards
  the `TestContext` and ER setup). Transplanting that shape to SKAT-O could
  skip something else while still passing at gate 1. That is the most likely
  way the byte-identity prediction fails, and it would fail loudly.
- **Option (b) diverges from upstream** on the user-weights path, where
  upstream's pseudo contribution stays weighted. Unlike the rest of #118 this
  one IS upstream behaviour, so the divergence is deliberate and recorded.
- **A free reconciliation falls out of (b):** `MAC_GroupVec` and
  `genoSumcount_noweight` become the same quantity computed two ways
  (`pr.MAC` against `arma::accu` of the non-zero genotypes). They must agree.
- **docdrift patterns are matched against FLATTENED text** -- `flatten()`
  blanks backticks, commas and apostrophes to spaces, so a pattern containing
  any of them silently never fires. Selftest catches it; mine did.

## Still open, not decided

`#114` (dead preconditioner globals), `#115` (two help strings state the
opposite of their declared default -- a two-line fix, no decision needed),
`#116` (a region analysed but producing no row leaves a 0-byte file; mechanism
NOT diagnosed), and `#117` itself (drop weight-0 markers before the kernel, or
leave binary FlexRV a `--r.corr=1` run) -- the last now has an exact,
golden-free gate if it is taken.
