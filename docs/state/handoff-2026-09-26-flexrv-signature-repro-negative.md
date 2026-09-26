# Handoff 2026-09-26 (evening): the FlexRV ultra-rare-only signature does NOT reproduce in the repo

Follows `handoff-2026-09-26-oldnew-comparison-and-flexrv-signature.md`, whose
OPEN item asked for a repo repro of the afr/BenIntNeo signature (three genes
at p_FlexRV 1e-120..1e-256, every extreme set row `Number_rare = 0`,
`BETA_Burden` 260588 / `SE_Burden` 7619 on MAC 86). Three drills on the
bundled fixture, all negative; the scan of every sim50k FlexRV output on
disk, negative. Not a LEDGER row: nothing in the tree has been shown wrong.

## What was tested (bundled fixture, binary sparse-GRM model + `$SG`, `--r.corr=1`)

Drill 1, pure rescale. One custom weight set with every weight = c on GENE1
(`GENE1 weight:S c c ...`), c in {1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10}, at
`--maxMAF_in_groupTest=0.005` (MAC <= 10: the cell is ONLY the collapsed
pseudo-marker, `Number_rare = 0`, MAC 54) and 0.01 (4 rare + 17 UR, MAC
119); binary and quantitative. `Pvalue_Burden` BYTE-IDENTICAL across all six
scales in all four arms; `BETA_Burden` and `SE_Burden` scale exactly by 1/c.
|Z| stayed below `--SPAcutoff` (2), so this arm exercised the score test only.

Drill 2, saddlepoint forced. Same sets with `--SPAcutoff=0.5` (0.1 is refused:
the minimum is 0.5), plus the MAC <= 1 cell (`--maxMAF_in_groupTest=0.0005`,
4 singletons, |Z| 2.4, SPA fires at the default cutoff), plus an
`--is_fastTest=TRUE` arm. p identical from c = 1e-2 down to 1e-10 in every
arm. Two steps between c = 1 and c = 1e-2, both explained and both small:
the MAC <= 1 cell goes 2.940e-2 -> 1.843e-2 because ER fires only on a set
whose weights are all exactly 1 (`sets[r].erAllowed`, `src/Main.cpp:2636`);
the MAC 54 cell goes 5.9425e-1 -> 5.9445e-1 (3e-4 relative), which is the
parked #143 finding, the SPA root's ABSOLUTE tolerance `eps^0.25`
(`src/SAIGE_test.cpp`, `tol1`), and it vanishes once the root is large.
Neither is the AoU shape, which is ~250 orders of magnitude.

Drill 3, collapse cutoff raised. Golden P58's command (the full ladder, 112
sets) with `--MACCutoff_to_CollapseUltraRare` 10 / 20 / 50 / 100. At 20 and
above every one of the 224 rows is `Number_rare = 0`, and the largest |Z| is
2.444 at every cutoff (GENE1 `REVEL_lof_MAC1`, unchanged from 10); the pooled
rows move from 0.0956982 / 0.93903 to 0.0965245 / 0.938945. Raising the
cutoff on this fixture inflates nothing.

Also negative: the P58 golden already mixes weights of 1 (LoF, score 1) and
~1e-38 (`sbeta20` on REVEL 0.01) in one pseudo-marker and its `sbeta20_MAC*`
rows converge to the `lof_MAC*` rows as they should. And every sim50k FlexRV
output on disk (164 files under `harness/work/flexrv_study/`, cutoff 10, both
traits, null series and planted signal): the largest |Z| on a
`Number_rare = 0` row is 4.13 on the null replicates and 9.3 on phenotypes
with planted causal genes, i.e. signal.

So on this build the burden test is scale-invariant to 1e-10 on every path a
FlexRV set can reach, and no cell shape the fixture can make -- ultra-rare
only, mixed, singletons, everything collapsed -- produces the signature.

## What the recorded numbers imply, and what to pull next

The 2026-09-26 handoff records `AM_lof_unw` at beta 792 with a SANE p. The
score builder (`harness/tools/flexrv_score_from_alphamissense.py`) gives every
LoF variant score 1.0, so `lof_unw` is a UNIT-weight burden. For a unit
burden, beta = S/var1 and Z^2 = S^2/var1, so var1 = (Z/beta)^2 ~ 1.6e-6 at
Z ~ 1 -- against MAC x mu(1-mu) x VR ~ O(1) for a hundred carriers at any
plausible prevalence. That is a 1e6 shortfall in var1 (or a 1e-3 scale on the
vector) on a set whose p looks fine, and it is a sharper anomaly than the
extreme rows, because it cannot come from any weight. It also cannot be
checked from here. The cheapest next measurements, all aggregate output:

1. One extreme gene, ALL its set rows: `Weight`, `MAC`, `Number_rare`,
   `Number_ultra_rare`, `Pvalue_Burden`, `BETA_Burden`, `SE_Burden`. Compute
   var1 = (Z/beta)^2 per row and compare it with MAC x prevalence x (1 - prev).
   If var1 is ~1e-6 on the unit-weight rows too, the vector or the variance
   is scaled for the whole gene, not for the weighted sets.
2. The model's `varianceRatio.txt` (every row) and the run's
   `--cateVarRatioMinMACVecExclude` / `MaxMACVecInclude`; a MAC band whose
   ratio is tiny inflates every cell in that band.
3. The wrapper's exact step-2 flags. The owner reports the burden arm did
   NOT use the default "max MAC cutoff" (which flag is not yet on record:
   `--MACCutoff_to_CollapseUltraRare`, `--max_MAC_for_ER` or
   `--maxMAC_in_groupTest`). If it is the collapse cutoff, a value >= 50 is
   what makes the MAC25/MAC50 rungs `Number_rare = 0`; the fixture says the
   raise itself is inert, so the value alone does not explain the p.
   SAIGE-GENE+ collapses at MAC <= 10 and every calibration this branch has
   (`docs/state/handoff-2026-09-13-flexrv-sweep-rerun.md`, sim50k) is at 10;
   nothing in SAIGE-GENE+ licenses another value.
4. The single-variant rows of the same gene (`.singleAssoc.txt`, if written):
   a marker whose per-variant p is itself absurd would move this off FlexRV.

Drill scripts: `run.sh`, `run2.sh`, `run3.sh` in this session's scratchpad;
each is P58's command with one flag changed, and the weight-line group files
are two `head -2` lines of `group_new_snpid_flexrv.txt` plus one
`GENE1 weight:S` line of fifty equal values.

## Later the same evening: the owner located the mechanism in the AoU output

The owner read the FlexRV arm's `.singleAssoc.txt`: the extreme genes'
pseudo-marker rows (e.g. `ENSG00000182759:pLoF;damaging_missense;other_missense:0.001:AM_lof_unw`)
carry `Is.SPA = false` with `p.value == p.value.NA == 3.943736e-174`, while
the SKAT-O arm's row for the SAME pseudo-marker has `Is.SPA = true`, the
same `p.value.NA`, and `p.value = 0.00063`. So the score test agrees between
the arms (same Tstat, same var1) and the burden arm REPORTED THE UNADJUSTED
NORMAL-APPROXIMATION P, which on a rare binary trait is the inflation the
saddlepoint exists to remove. That also explains the betas: 792 on the
unit-weight set is Tstat/var1 with the score-test SE back-derived from a
p of 1e-174. Every "impossible |Z|" is p.NA restated.

In code, `p == p.NA` with `Is.SPA == false` is `getMarkerPval`'s fallback
(`src/SAIGE_test.cpp` after the SPA block, `t_pval = t_pval_noSPA` when
`t_isSPAConverge` is false). It is reached when (a) |Z| never exceeded
`--SPAcutoff` (impossible at 1e-174 unless the cutoff was raised), (b) a
root finder returned non-convergence -- `q >= gpos || q <= gneg` (LEDGER
#22's unattainable-q rule), a NaN Newton step, or 1000 iterations -- (c)
`Get_Saddle_Prob_Binom_core` returned `isSaddle = false` (non-finite k1/k2,
negative `temp1` or k2, `w == 0`), (d) the linear-scale SPA p underflowed to
exactly 0 (LEDGER #64's rollback), or (e) `R::qnorm` overflowed. The ER
branch is NOT a candidate: it sets `Is.SPA = true` and reports the exact p,
and non-unit weight sets never take it (`erAllowed`).

Since the SKAT-O arm converged on the same vector, the SPA INPUTS differed
between the arms, which means a flag (or build) differed. Fixture drill 4/5,
`--r.corr=1` vs SKAT-O x `--max_MAC_for_ER` 4/100 x `--is_fastTest` x with/
without `--sparseGRMFile`: the pseudo-marker rows are byte-identical between
the BURDEN and SKAT-O arms in every pairing, the SPA fires on all 16 rows
above the cutoff and converges in every arm (the no-GRM arms move p by a few
percent through var2 and the fast kernel, and still converge). So nothing the
fixture can vary makes the burden arm fall back. What is needed from the
workbench, in order:

1. The two step-2 command lines verbatim (FlexRV arm, SKAT-O arm) and the
   image tag each ran on. The diff between them is the determinant.
2. Systematic or sporadic: over the FlexRV arm's `.singleAssoc.txt`, the
   count of pseudo-marker rows with `p.value.NA < 1e-4` split by `Is.SPA`.
   All false = a flag that disables or breaks the SPA on that arm (a raised
   `--SPAcutoff`, or a kernel/GRM choice); a fraction false = data-dependent
   non-convergence at N = 50,000 (rule b or c above), which the fixture cannot
   reach and sim50k with the same flags can.
3. For one failing row: `Tstat`, `var`, and the model's case count, so `q`
   against `gpos` (the sum of positive adjusted genotypes, roughly the
   carrier count) can be checked by hand for rule (b).

Consequence for the results regardless of the cause: any FlexRV-arm gene
whose driving set has `Is.SPA = false` in the single-variant sidecar is not a
result. The SKAT-O arm's pseudo-marker p (0.00063 here) is the right number
for that set, and the SKAT-O arm is unaffected by construction.

## The wrapper's flags (owner pasted `saige_step2_brava_main_original.sh`, same evening)

Both group arms get: `--minMAF=0 --minMAC=0.5 --LOCO=FALSE
--is_output_moreDetails=TRUE --is_single_in_groupTest=TRUE
--is_output_markerList_in_groupTest=TRUE`, `--is_Firth_beta=FALSE`
(the notebook's `submit()` defaults `firth="FALSE"`, which overrides the
wrapper's `auto`; the owner corrected an earlier draft of this paragraph
that said Firth was on), the sparse GRM and
`--relatednessCutoff` (RELCUTOFF is required when a GRM is given), and NO
`--is_fastTest` (FASTTEST auto = variant runs only). The wrapper has no way
to pass `--max_MAC_for_ER`, `--MACCutoff_to_CollapseUltraRare` or
`--SPAcutoff`, so those are at the build defaults on both arms. The arms
differ ONLY in: `--r.corr=1` + `--flexRV_score/maxMAF/lofAnno=pLoF` + one
`--maxMAF_in_groupTest=0.001` (FlexRV, pan-ancestry group file `am_grp`)
against no `--r.corr`, three max MAFs and the per-ancestry group file
(SKAT-O), `MIN_GROUP_MAC` if set, and the IMAGE TAG, which `submit()`
records per run in `RUNS[name]["image"]`.

Reviewed for `--r.corr=1`-only mechanisms and found none that can reach the
sidecar row: the BURDEN twin memo (`src/Main.cpp:4325`) serves the REGION
burden record, not the ultra-rare pass's single-variant row, and the UR pass
has no memo; `is_gtilde` is assigned inside `getMarkerPval` before it is
read; `SPAWorkspace` views its buffers at the current `n` (`ws.ensure`), so
a previous vector's length cannot leak. With identical `p.value.NA` the
score test saw identical Tstat and var1 in both arms, hence identical q,
gtilde and mu, so on one build the two SPA calls are the same call. The
first thing to establish is therefore whether the arms ran on the same
image; the second is the failing rows' position in set order within a gene
(the first set, `unw_unw`, is never served or preceded; a pattern by set
order points at per-worker scratch state, every set failing points at the
inputs).

## sim50k planting, both arms, the wrapper's flags: NEGATIVE (same evening)

Both AoU arms ran on ONE image (`saige-slim:1.5.2-dev-f2fb24da`, from
`RUNS[...]["image"]`), Firth off on both (`submit()` defaults
`firth="FALSE"`), `--nThreads` at the build default 1 (the cell tests are
serial). So the divergence is inside one build on identical inputs, and the
last local lever was N. Planted on the sim50k fixture
(`/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k`, model
`harness/work/null50k-20260908/step1_y_null_bin`, 500 cases of 50,000,
`sparse` VR 1.0000): five genes of 86 SINGLETONS each with distinct
carriers, 0 / 7 / 14 / 26 / 40 of them cases (`scratchpad/plant50k/plant.py`
picks them from `allcounts.acount` + the bed), annotation `pLoF`, score line
all 1. FlexRV arm (`--r.corr=1 --flexRV_score=AM --flexRV_maxMAF=0.001
--flexRV_lofAnno=pLoF --maxMAF_in_groupTest=0.001`, 144 sets) and SKAT-O arm
(`--maxMAF_in_groupTest=0.0001,0.001,0.01`), both with the wrapper's
`--minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_output_moreDetails=TRUE
--is_single_in_groupTest=TRUE --is_output_markerList_in_groupTest=TRUE
--is_Firth_beta=FALSE` and the sparse GRM at cutoff 0.

Result: the pseudo-marker sidecar rows are IDENTICAL between the arms for
every gene (AC, Tstat, var, p, p.NA, Is.SPA), the SPA fires on every gene
with signal and converges every time -- PLANT07 p.NA 9.0e-9 -> p 3.8e-5,
PLANT14 8.7e-30 -> 5.2e-12, PLANT26 4.3e-92 -> 3.3e-28, PLANT40 9.9e-199 ->
8.5e-50 -- and of 720 FlexRV-arm UR rows the 144 with `Is.SPA = false` are
all PLANT00 (|Z| 0.96, below the cutoff). `var` for 86 carriers at 1%
prevalence is 0.87-1.62, i.e. MAC x mu(1-mu) as it should be.

So a fallback at p.NA 1e-174 that the SKAT-O arm does not share is not
reproducible with this build on any input this repo can make. What is left
is in the two AoU rows themselves, and the next step is to read them whole.

## RESOLVED (late 2026-09-26): LEDGER #163 and #164, reproduced on the kernel with no AoU data

The owner pasted the two sidecar rows whole (all columns; the cell's
counts are below the AoU disclosure floor and are not reproduced here).
Two columns settle it. `Tstat` and `var` together give the carriers' fitted
probabilities, which are of order 1e-3 and below on a 7% trait (age in the
model; young carriers), and `var` = sum mu(1-mu) is CORRECT, not deflated.
The normal approximation then gives |Z| ~28, p 3.9e-174, for a cell whose
exact tail is of order 1e-3. The SKAT-O arm's 6.3e-4 is right and came from
the EXACT RESAMPLING test (the cell is at or below `--max_MAC_for_ER`; a
beta-weighted set is ER-eligible; the ER branch sets `Is.SPA = TRUE`, which
is what misled the arm comparison). The FlexRV arm's `AM_lof_unw` set is
ER-INELIGIBLE (only `unw_unw` is), so it went to the saddlepoint, whose
Newton root finder starts at t = 0 with a step of (q - K1(0))/K2(0) of
several hundred, past the exp underflow, and returns non-convergence for a
root that exists (#163); SAIGE then reports the normal p as `p.value` with
`Is.SPA = FALSE` (#164). `BETA` 792 is `Tstat/var` and `SE` 28 is
back-derived from the 1e-174, so the region file carries no signal of the
fallback.

Kernel repro, `tests/testthat/test-unit-spa.R` (#163/#164 blocks, run at
standard: 62 pass): a SYNTHETIC two-carrier shape at mu = 1e-3 converges
(root 6.907); at 6e-4 the finder returns `-Inf`, `Isconverge = FALSE`, and `SPA_pval` hands
back `pval_noadj` = 6.5e-183 against an exact 1.2e-3. Threshold for that
shape between 7.0e-4 and 6.8e-4 (first step ~730). Why no earlier drill saw
it: the bundled and sim50k models have min mu 0.005, so the first step is
at most ~100 there, and it converges. The wrapper's flags, the image, Firth,
threads and the raised cutoff were all red herrings, and "Number_rare = 0"
was a coincidence of the shape (MAC <= 4 cells are collapsed cells).

Consequences and decisions (owner's): #163 has a candidate fix on the row
(bracket then safeguarded Newton, both kernels; no golden moves, predicted
from the models' min mu). #164 is design: whether the fallback should be
visible on the region path, replaced by ER when MAC <= max_MAC_for_ER
regardless of weights, or NA. For the AoU FlexRV results as they stand:
a hit whose driving set has `Is.SPA = FALSE` in the sidecar is not a
result; for MAC <= 4 driving sets the SKAT-O arm's pseudo-marker row for the
same mask is the correct p; and the arm-wide count of `Is.SPA = FALSE` rows
with `p.value.NA < 1e-4` is the size of the problem (still to be measured).

## #163 FIXED (same day), gated; #164 awaits the owner's choice

`rescueRootK1` (`src/SPA_binary.cpp`): after the plain Newton fails on an
attainable q, bracket the sign change of K1 - q by doubling from 0, then
safeguarded Newton with bisection fallback. Both root finders. Zero blast
radius by construction (only runs where the old code failed) and confirmed:
no golden reaches the fallback today, and `harness/prepush.sh` (standard
tier, `SAIGE_COMPARE_MODE=both`) is green with every golden byte-identical.
`test-unit-spa.R` #163 block inverted (68 pass). Two findings on the way:
the synthetic vector must be CENTRED to be realistic -- an adjusted genotype's
negative entries sum to about minus the carrier count, which is what makes
the mirrored q attainable; with no negative entries (the `noadjCov` raw-G
path) the mirrored tail is unattainable by construction and still falls into
#164. And the composed SPA p for the two-carrier shape is 1.6e-4 for an
exact tail of 1.2e-3: right order, ~7x anticonservative, the saddlepoint's
own limit at a discrete extreme, not the root's. NOT committed: the owner
commits (one commit for the row: `src/SPA_binary.cpp`, the test file, the
LEDGER row, this handoff). #164's three options are on its row; the owner
has not chosen. `docs/state/next.md` has another session's uncommitted
edits and was not touched; its FlexRV OPEN item is resolved by this file.

## #164 FIXED, #165 raised and DECIDED (late 2026-09-26); owner chose both halves and rejected NA

Owner's decisions in order: "let's go with both" (visibility AND a conservative
fallback), "why would we want NA?" (so the fallback is the EXACT test, not
NA), and "we shouldn't worry too much about super high MAC for ER; this
should be low by default; simply reject MACs that are very high for this
option" (so `--max_MAC_for_ER` is capped at 200, the largest value ever
verified, golden P33, and the fallback applies the same ceiling).

What landed (four files under `src/`, three under `R/`, one new test file,
CLI help, roxygen and `man/`, CLAUDE.md, LEDGER #164 FIXED and #165 DECIDED):
- `getMarkerPval`: the ER block is a lambda; on an attempted saddlepoint that
  does not converge it runs in place of the normal p, any weight set, MAC
  <= 200; above the ceiling or if the exact test throws, p is NA with the
  reason; `Rcpp::LongjumpException` and `NeedsMainThread` are rethrown (the
  existing worker deferral re-runs the call on the main thread).
- `MarkerResult::tested`: a tested marker whose p is NA keeps its row (the
  writer's `pval != "NA"` rule meant "QC-dropped"; it now reads the flag too).
- Records: `KernelWorkspace::spaFallback` + `spaFallbackNote` set per call,
  `SAIGE::spaFallbackCount()` atomic, `recordSpaFallback` at the marker,
  single-in-group, ultra-rare and BURDEN sites; `R/SpaFallbacks.R` prints
  `SPA fallbacks: N ...` at the end of every step-2 run and writes
  `<SAIGEOutputFile>.spaFallbacks.txt` when N > 0.
- Plant `SAIGE_SPA_PLANT_FAIL=1`; `tests/testthat/test-step2-spa-fallback.R`
  (48 pass): negative control, planted marker run, two-ways against native
  ER at 200, the region path, the CROSS-ARM control (planted FlexRV
  `unw_MAC1` == SKAT-O arm's native exact p on the same vector: the AoU
  discrepancy closed end to end), and the refusal of `--max_MAC_for_ER=201`.
- #165: the native exact test throws a non-standard exception on this
  fixture at rs23 (folded MAC 607, |Z| above the cutoff) -- `--max_MAC_for_ER
  =2000 --idstoIncludeFile=<rs23>` reproduces it, 500 runs; identity not
  established (no `throw` in `src/` makes it). Mitigated by the cap.

Two lessons paid for today: NEVER edit `extdata/`, `R/`, `src/` or `tests/`
while `harness/prepush.sh` runs -- it tests the working tree, and a help-text
edit mid-run made three golden cases read a half-written CLI script and fail
to parse (the gate was rerun clean); and `harness/build.sh --incremental`
refuses after a header edit, so a `.hpp` change means the clean build.

No golden moves (no bundled saddlepoint fails; the plant is off by default).
Not committed: the owner commits, one commit per row is the convention
(#163: `src/SPA_binary.cpp` + `test-unit-spa.R`; #164/#165: the rest).

## What this means for the AoU production results

The defect is upstream's too, so the 2024 comparison could not see it. Every
binary run on the old image is affected wherever a saddlepoint failed: the
FlexRV arm (transform sets never take the exact test natively, so a tiny-mu
ultra-rare cell reports the normal p, and `p_FlexRV` inherits the minimum),
the SKAT-O arm's single-variant sidecar rows for real markers above MAC 4,
and the VARIANT-level runs for any binary marker above MAC 4 with tiny-mu
carriers. The region-level SKAT-O `Pvalue` does not go through the
saddlepoint and is unaffected. Cheapest path: rebuild and push the image
(`make docker`, then tag), rerun the FlexRV arm for every binary trait
(the `SPA fallbacks:` line in each log is now the size of the problem per
trait), and for the variant-level arm count `Is.SPA = FALSE` rows with
`p.value.NA < 1e-4` per trait from the existing sidecars to decide which
traits to rerun -- a trait with zero such rows needs nothing.
