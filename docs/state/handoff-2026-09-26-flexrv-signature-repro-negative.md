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

## Night: the first AoU rerun, LEDGER #166 and #167, and what is still open

The SAS FlexRV rerun on `1.5.2-dev-adfe6266` (filtered to the scan's regions,
`REGIONS` in the wrapper) worked as designed for five traits (Asth: 107
regions, 21 chromosomes in 5-20 s each, 9,554 rerun pseudo-marker rows, none
with the old signature, none NA) and was SLOW on HTN's chr16: 12m11s for 27 fallbacks including a
209-carrier cell that the fallback sent to the exact test's resampling
regime, against 5-15 s for every other chromosome -- it was read as a hang
and a `ddel` was attempted, but it finished normally (SAIGE_SECONDS=912).
On the kernel and on sim50k the exact test at 13-209 carriers takes under
5 s at N = 50,000, with weights of 1 or 2e-6, and with carriers chosen at
the extremes of mu (max/mean 4.7); so the ~x10-100 on real data is a
data property not pinned down (`aou_spa_rerun.md` section 4 has a
per-chromosome `timeout` if it is ever worth naming the cell, and timing
`erExactTestInCPP` on that cell's (Z, y, mu) INSIDE the workspace is the
probe). The fallback now runs the exact test only in its enumeration regime
(13 carriers) and reports NA with the reason above that, which bounds it.

Owner's decisions, in order: both halves of #164; NOT NA when an exact test
is possible; cap `--max_MAC_for_ER` at 200 rather than chase the MAC-607
crash (#165); and "why can't we do similar things on the burden path?" --
which became #166: every weight set is ER-eligible at MAC <=
`--max_MAC_for_ER`, as beta sets always were. On the way, #167: the exact
test's ABSOLUTE tie epsilon makes every configuration a tie under weights
of 1e-3 and below (mid-p exactly 0.5, silently); the vector is now scaled to
a largest entry of 1 first, which leaves every 0/1/2 vector, and so every
existing exact-test row, byte-identical.

Robustness, as asked: `tests/testthat/test-unit-er.R` is an R REFERENCE
implementation of the exact test (enumeration, odds within a case count,
Q = S^2, the kernel's mid-p and tie rule) fed the kernel's own
weight-independent P(j) through the test export `erExactTestInCPP`; it
matches the kernel to 1e-10 on 80 unit-vector cases (the control) and 150
weighted ones, on weights spanning six orders in one vector, and is exactly
scale-invariant from 1e-8 to 1e3. The cross-arm test asserts, natively and
without any plant, that the FlexRV `unw_MAC1` row prints the SKAT-O arm's
exact p. Goldens: exactly the predicted movers moved (P58 and P61, 16 set
rows + GENE1's pooled row each), reproducible across two trees, every
control identical; re-baselined.

Also learned: the exact test's p is a MID-p (`pval - pval_same/2`), which
the reference had to adopt before it matched; and the AoU sidecar's
`AC_Allele2` for a weighted pseudo-marker is the WEIGHTED sum, so the
carrier count must be read from the results file's per-set `MAC`.

Next: commit (rows #163-#167), `make docker` + push, rerun the SAS six on the
new image (most of the previously failing rows now take the exact test
natively, so `SPA fallbacks:` should be near zero and the sidecar's Note
column says why any remaining saddlepoint failed), splice, then the other
ancestries from their scans. Wrapper: consider delivering
`.spaFallbacks.txt` (a `--output OUT_FILE_SPAFALLBACKS` next to the other
sidecars) so the reasons are readable without a diagnostic job.

## STATE FOR THE NEXT SESSION (written 2026-09-27 ~17:15, before a context clear)

**Where the code is.** Committed on `strip`: #163-#168 (last: `8ca5d0e4`
clause (b) + the dump hook; `08641f7f` the dump-check script). Working tree:
`src/SPA_binary.hpp`, `src/SPA.cpp`, `src/SAIGE_test.cpp` -- the sidecar Note
now also carries `[root1=<t>/<niter> root2=<t>/<niter> qsaddle=valid|INVALID|n/a
bound=<chernoff>]`, and the unused `kExactTestMACCeiling` constant is gone
(the ratchet caught it). Built, `test-unit-spa.R` 88 and `test-unit-er.R`
471 pass; the full gate on it was launched at 17:11 (`scratchpad/gate_diag.log`
of the session that wrote this) -- it PASSED (6573, both modes, 10m08s) and
the change is committed with this section. `docs/state/next.md`
has another session's uncommitted edits; leave them.

**Images on gcr.io/weighty-elf-452116-c7/saige-slim**: `1.5.2-dev-adfe6266`
(#163-#165), `-7767bfee` (+#166/#167), `-a6213623` (+#168 a), `-8ca5d0e4`
(+#168 b, dump hook) -- the current production one. The Note diagnostics
above are NOT in any image yet.

**AoU state (SAS, the pilot ancestry).** Scans done for all five ancestries
(SKAT-O arms clean everywhere; FlexRV arms: every binary trait fails, 30-150
regions each). SAS: Asth/EFRMB/HTN/PAD/T2Diab rerun on `a6213623` with ZERO
fallbacks, spliced into `sas/flexrv_v4_fixed/` (splice controls green).
BenIntNeo rerun on `8ca5d0e4`: STILL 24 fallbacks, all `root(mirrored q) not
found` BY THE FINDER, 20 above 13 carriers -> NA -> 102 rows omitted (twins
of NA primaries go unrecorded); NOT spliced. A dump run
(`fxv5-sas-benintneo-dump`, `SPA_DUMP=1`, output under `sas/diag/`) was
launched ~17:00 to capture those 24 vectors; the analysis job is
`harness/tools/aou_spa_dumpcheck.sh` (recipe section 5) on the SAIGE image;
its log's four lines per dump say whether root(q) converges, whether its
saddle is `isSaddle=TRUE`, whether the mirror is attainable, where K1-qinv
changes sign vs the bracket cap 2^20/max|g|, and the Chernoff minimum. On
every synthetic vector at realistic N the mirrored root is within the
bracket (~-6e4 at N = 8,000, scaling with N), so the real cells have a
property the constructions lack; the two candidates the numbers will
separate are "q-side saddle invalid at its root" and "bound not small".
Whatever it is: fix, gate, image, rerun BenIntNeo, splice; then the other
four ancestries follow the recipe (scan summaries already delivered under
`spa_scan/<anc>/flexrv_v4/<stamp>/`; afr/BenIntNeo's three 1e-174 genes are
the positive control and should land near 1e-3).

**Notebook cells that must exist** (a fresh kernel loses them): env
(`GOOGLE_CLOUD_PROJECT`, `PET_SA_EMAIL`, `WORKSPACE_BUCKET`, `BUCKET`, `GRM`,
`STEP1`, `REG`, `RUNS`), the helpers cell (`submit` with the
`OUT_FILE_SPAFALLBACKS` and `OUT_FILE_SPADUMP` outputs), `load_runs()`, the
rerun cell (`IMG5`, `STAMP=20260926-170535` for sas, `regions_for`,
`submit_rerun`), `submit_splice`, and on disk the patched wrapper
(`grep -ci "spafallbacks\|REGIONS\|spadump"` >= 9), `aou_spa_splice.sh`,
`aou_spa_dumpcheck.sh`. Job names are dsub job names; `job_id` in RUNS is
always None (the parse never worked) -- status by `dstat --age`, ids by
`dstat --full --format json` for `ddel`.

**Addendum 17:40.** The dump run delivered 20 files for 24 fallbacks: each
chromosome's SAIGE process numbers its dumps from 0 into the wrapper's single
`spadump/`, so later processes overwrote same-numbered earlier files. Fixed in
the working tree (`src/SAIGE_test.cpp`: the name carries `getpid()`), BUILT
but NOT gated and NOT committed -- run `bash harness/prepush.sh` and commit it
with the next source change. The 20 vectors are enough: the analysis job
(`spadumpcheck-sas-benintneo`, `harness/tools/aou_spa_dumpcheck.sh` on
`1.5.2-dev-8ca5d0e4`) is the next step and its log is what to read first.

## LEDGER #169 (2026-09-27, evening): the dumps said "q-side saddle invalid", and why

The `spadumpcheck-sas-benintneo` log answered the question the addendum left:
in every one of the 20 dumps `root(q)` converged (2.9e3..2.5e5, 5-9
iterations, `K1 - q` ~1e-15), the mirror was UNATTAINABLE (`root(qinv)` Inf,
0 iterations -- the #168 (a) rule, not the finder), the q-side saddle was
`isSaddle=FALSE` with `Korg(root1)=Inf` and a finite `K2(root1)` (2.6e-11..
7.9e-8), and the Chernoff bound on the mirror was 0. So: **q-side saddle
invalid at its root**, and the bound is not the problem. The sidecar's
`root(mirrored q) not found` was the final else's label (reason 2 chosen on
`outuni1.Isconverge` alone) after the (a) branch had run and been refused.

**Mechanism.** `Korg_*_ws` sums `log(1 - mu + mu e^{g t})` as written and
`exp` overflows past `g t` ~709; K1 and K2 are written in `e^{-g t}` and stay
finite, so Newton is happy where K is Inf. The root is that large when the
cell's positive mass is dominated by carriers with tiny fitted probabilities
and S is within a hair of gpos: in the single-dominant-carrier model of a
covariate-adjusted vector, the rest contributes its expectation plus a
leverage term, and `root * max|g| = (Delta/sd) * max|g| / sd(rest)`, roughly
1/mu of the dominant carrier -- 730 to 6,000 on the dumps. My synthetic
vectors had dominant carriers at mu 1e-2 and up, so `root * max|g|` was ~10.

**Fix (LEDGER #169, one commit).** When the plain sum is non-finite, K is
recomputed term by term as `x + log((1 - mu) e^{-x} + mu)` for `x > 0`, in
both kernels; the finite path is untouched (bit-identical). The final else
labels a refused q-side saddle reason 3. Test: 20-sample vector, dominant
carrier at mu 1e-3, exact enumeration over 2^20 configurations as the
reference; the fixed p is 2.16e-4 against exact 1.57e-4, and the same shape
one notch inside exp's range is 1.77x -- r* is 1.3-1.8x conservative on this
boundary shape on BOTH sides of the overflow, so the fix changes
evaluability, not the approximation (recorded on the row, not changed). Red
at unpatched HEAD (9 assertions), controls green. Gate: see the commit.

**What the BenIntNeo cells will now report.** A p of order mu of the
dominant carrier times the rest's tail: 1e-3 to 1e-4, not NA and not
1e-100. Commit `4ff9576b`, image `1.5.2-dev-4ff9576b`; rerun BenIntNeo on it (regions
list unchanged), expect `SPA fallbacks: 0`, splice, then the other four
ancestries per `harness/tools/aou_spa_rerun.md`. The dump-name collision
(pid in the file name) ships in the same commit.

**Image pushed** (same evening): `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-4ff9576b`,
digest `sha256:9a44f52a...`, built by `make docker` from `4ff9576b`. Rollout
order and the per-stage convergence check are `harness/tools/aou_spa_rerun.md`
section 6 (owner's order: sas/BenIntNeo, then eas + afr, then amr, then eur;
`aou_spa_stagecheck.sh` over each stage's sidecars before the next launch).

**sas/BenIntNeo rerun on `4ff9576b` (same evening, `sas/flexrv_v6/`): SUCCESS,
sidecar header-only (0 fallbacks against 24 on `8ca5d0e4`), dump tarball
empty, and the results file equals the v4 file ROW FOR ROW in every rerun
region** -- the `8ca5d0e4` shortfall was confined to a single region and is
fully recovered. So the "102 rows missing for 20 NA cells" written above and
in LEDGER #168 was wrong (corrected on the row): the sidecar counts TESTS,
and the fastTest escalation pass re-tests a region and records its fallback
again, which the dump log shows as the same vector dumped several times
(identical q, N, root). The wrapper delivers the sidecars gzipped
(`.spaFallbacks.txt.gz`, ~92 bytes when header-only); `aou_spa_stagecheck.sh`
reads both forms since `f228e86b`. Next: splice (old v4, new v6), then
eas + afr.

## LEDGER #170 (2026-09-27, night): the region the #169 rerun could not produce

The splice of sas/BenIntNeo (old v4, new `flexrv_v6` on `4ff9576b`) REFUSED
85 rows. Read against v4 they were one region's 91 rows minus its recorded
failures and its Cauchy row: the region was absent from v6 altogether, and
from the `8ca5d0e4` run too. `chrStatus` said `failed 5` on both while both
jobs reported SUCCESS (the wrapper tolerates a dead chromosome). The job
log's chromosome-5 SAIGE output ended in `--flexRV_lofAnno label(s) 'pLoF'
match no annotation in the group file` -- the rerun wrapper writes one group
file per chromosome SUBSET to the rerun regions, chromosome 5 had one region
in the list, and it has no pLoF variant. In v4 the whole-chromosome file had
pLoF elsewhere, the file-level check passed, and that region simply had no
`lof` sets (90 of 96), which is what a region without the label always gets:
those sets keep no marker, have MAC 0, and the min-group-MAC gate refuses
them. The error was stricter than the code beneath it.

**Fixed (LEDGER #170)**: the `stop` is a `warning` naming the absent labels,
the labels present, and what the `lof` transform will do (key on the present
labels, or be refused in every region). Test inverted in
`test-step2-flexrv-guards.R` (red at unpatched HEAD, 8 assertions): the
present-plus-absent case is byte-identical to the default run; the
all-absent case equals the run that never asked for `lof` (keys identical, p
within 1e-10). Doc surfaces updated (roxygen, CLI help, CLAUDE.md,
`cli_status.Rmd`) with a docdrift rule (`170-lofanno-must-occur`, 0 hits at
the drill commit). `aou_spa_stagecheck.sh` reads `chrStatus` and flags a
lost chromosome. The five spliced sas traits all read `failed none`
(checked), so nothing in `flexrv_v4_fixed` is affected.

**Two facts for the rollout.** (1) On the other ancestries this WILL recur
wherever a chromosome's rerun list holds only regions without a pLoF
variant, so the reruns must go on the #170 image, and the stage checker's
chromosome line is now part of the per-stage check. (2) A splice refusal
whose unexplained rows are one whole region is this signature; read
`chrStatus` before anything else.

Gate 6611 pass both modes, zero goldens moved; commit `3e92d89d`; image
`1.5.2-dev-3e92d89d`. Next: rerun sas/BenIntNeo into `flexrv_v6` on it (only chromosome 5 should differ from today's run), splice,
then eas + afr.

**sas COMPLETE (2026-09-27, ~19:30 UTC).** BenIntNeo rerun on `1.5.2-dev-3e92d89d`
into `flexrv_v6`: `chrStatus` `failed none`, sidecar header-only, dump empty,
every rerun region present, chromosome 5's region back with its 91 rows.
Splice green: results 2,840 changed rows all explained (1,420 cells: the
1,414 counted before plus that region's recorded failures and Cauchy row),
totals unchanged; singleAssoc likewise; pooledTests lost the 8 rerun regions'
rows because the v4 artefact (p underflowed to 0) had fired the CCT fallback
rule and a plain Cauchy discloses nothing (the ROWS wording is per file since
`60e3297c`). Six of six sas traits are in `flexrv_v4_fixed/`. Next: eas + afr
on `IMG7` (section 6 of the recipe), stage checker over sidecars AND
chrStatus, then amr, then eur.

## eas + afr stage (2026-09-27 night, into 2026-09-28)

Reruns on `1.5.2-dev-3e92d89d`: 128 jobs (every binary trait of eas, and of
afr in its chromosome bins; the other submissions had empty region lists),
all delivered, `aou_spa_stagecheck.sh`: 0 fallbacks in every sidecar, no
lost chromosome, `STAGE OK`. So #169 lands on two more ancestries and the
#170 guard stayed quiet.

The SPLICES then found three structures of larger cohorts that the sas stage
(N = 2,645 in-sample, every marker ultra-rare) could not show, each first
refused by the splice's control as designed and then given an exact rule
(`aou_spa_splice.sh`, commits `e4606698`, `d69e85cd`, `5febb6b4`,
`0d0452cd`, `0df0707a`; recipe section 3):
1. **Real-marker singleAssoc rows** (MAC above the collapse cutoff, which
   exists once `2 N maxMAF > 10`): keyed by MarkerID, explained by a recorded
   real-marker failure or MAC <= the exact-test cutoff; a marker spanning
   overlapping regions (protocadherin, ZNF clusters) has one row per region,
   replaced as a multiset.
2. **Burden cells with `Number_rare > 0`**: a test the scan could not
   observe (its saddlepoint can fail while every component's own test
   converged; v4 recorded nothing). Accepted at region level, counted, with
   the direction reported (p rose in nearly all, as an artefact that was too
   small predicts).
3. **Pseudo-marker rows changed by #166** (every set now ER-eligible by its
   own unweighted count): bounded by the 0/1-weight twin's count.
Classifier over the 128 splice logs: 98 clean first time, 27 of the 30
resubmitted clean, 3 refused on the `lof x mafbeta` twin case, fixed in
`0df0707a`; those 3 are the next resubmission. Then afr/BenIntNeo's positive
control (scan-tail genes near 1e-174 should read ~1e-3 in `flexrv_v4_fixed`),
then amr, then eur. No AoU count below 20 is written anywhere in the repo
(checked 2026-09-27 over every change of the session, commit messages and
memory notes included; the one such number written was scrubbed).

**eas + afr COMPLETE (2026-09-28).** The last three splices passed on
`0df0707a`. afr/BenIntNeo positive control over the scan's tail regions:
the three artefact genes went 2.3e-256 -> 1.4e-4, 7.8e-120 -> 9.8e-8 and
3.9e-174 -> 6.3e-4; the six genuine tail regions (8e-7 .. 8.6e-6) are
unchanged to the printed digits -- only the failed tests moved. The gene at
9.8e-8 is a minimum over ~100 weight sets, not a stop; worth reading (set,
MAC, SKAT-O arm) when afr's results are read. Next: amr on `IMG7`
(`STAMPS["amr"]` from the scan listing), stage checker, splices; then eur.

**amr COMPLETE (2026-09-28).** 118 reruns on `IMG7` (the other submissions
had empty region lists), stage checker `STAGE OK` (0 fallbacks, no lost
chromosome), 118 splices clean on the first pass with the four
larger-cohort rules in place -- nothing new. amr runs WITH its sparse GRM
like the others (recipe corrected, `52fbaec3`). Next: eur on `IMG7`, on
demand (its bins run longer), then the plotting pass from per-ancestry
gene-level summary tables written in the bucket (the Saturday egress
warning coincided with plotting from the raw per-trait files).

## OPEN: a submission system for thousands of traits (owner, 2026-09-28)

The rollout's mechanics do not scale: eur was 320 serial `aou_dsub` calls
from the notebook (seconds each), then 320 serial `fetch_log` reads to
classify verdicts, then a second wave of 320 splice jobs. At thousands of
traits that is hours of submission before any compute. Levers, in order:
1. `dsub --tasks` with a task table (trait, bin, regions, outputs): one API
   call per stage. Already identified for step 1 in
   `docs/state/` (the step-1 scaling handoff).
2. Verdicts as files, not logs: the wrapper writes a one-line status per task
   (fallbacks, failed chromosomes, splice verdict) under one prefix; one
   `gcloud storage cat` over the prefix replaces per-job log fetching, and
   displays kilobytes (the egress-safe pattern).
3. The splice runs at the END of the rerun task, in the same container, so
   the second wave disappears and the controls run where the data is.
4. Several traits per genotype pass (`--GMMATmodel_varianceRatio_multiTraits_File`
   on the region path) -- only for traits fitted on ONE sample set; measure
   how many AoU binary traits share one before designing around it.
Not started; the current notebook loop is fine for the last stage of this
rollout and not for production.

## ROLLOUT COMPLETE (2026-09-28)

**eur:** 198 reruns on `IMG7` (of 320 binary-trait bins; the rest had empty
region lists), stage checker `STAGE OK` (0 fallbacks, no lost chromosome),
198 splices SUCCESS. **All five ancestries** are now rerun on
`1.5.2-dev-3e92d89d` (LEDGER #163-#170) and spliced into
`<anc>/flexrv_v4_fixed/`: sas 6 traits, eas + afr 128 jobs, amr 118, eur 198
-- every FlexRV binary-trait job the scans flagged, with zero saddlepoint
fallbacks anywhere on the final image and every splice passing the
"only failed tests moved" control (with the four larger-cohort rules of
`aou_spa_splice.sh`, self-tested in `aou_spa_splice_selftest.sh`).
The afr/BenIntNeo positive control held (artefact genes to 1e-4..1e-7,
genuine tails unchanged). The SKAT-O arms never needed rerunning (scans
clean). The quantitative traits have no saddlepoint and were never
affected.

**Next:** the plotting pass -- what FlexRV added over SKAT-O/burden -- from
per-ancestry gene-level SUMMARY tables written by a job into the bucket
(kilobytes; figures only displayed), never from the raw per-trait files:
the Saturday egress warning coincided with plotting from those. The summary
job needs the SKAT-O arm's results path pattern and the columns to plot,
which the owner has and this handoff does not. After that, the OPEN item
above (a submission system for thousands of traits) before any production
launch.

**eur checks (2026-09-28).** Over the 198 splice logs: no WARNING, no
`ROWS .txt.gz` line, every unobservable burden cell rose and none fell,
real-marker rows in every job (eur's MAC range). Positive/negative control
over the scan's eight most extreme v4 burden tails among the rerun binary
traits: BRCA1 (5e-78) and BRCA2 (3e-61) for breast cancer, PKD1 (6e-43)
for CRF and APC (1e-31) for BenIntNeo are unchanged to the printed digits;
the other four (2e-88, 9e-70, 2e-38, 7e-32) now read 4e-4 to 1e-3 -- the
artefact, gone; the biology, untouched.

## Plotting layer (2026-09-28)

`harness/tools/aou_flexrv_summary.sh` (one dsub job per ancestry, cloud-cli
image, `summary2-<anc>`) writes two gene-level tables under
`SAIGE_SLIM_TESTING/summary_v4/`: `<anc>.flexrv_vs_skato.v4.tsv.gz` (one row
per trait x gene: FlexRV pooled and best set, SKAT-O pooled and best cell,
best burden and SKAT, fixed-vs-v4 provenance) and
`<anc>.skato_cells.v4.tsv.gz` (every SKAT-O cell wide, burden and SKAT-O p
per mask@maxMAF). A separate plotting notebook (`05_flexrv_power.ipynb`,
workbench) draws from those alone: the per-(ancestry, trait) QQ grid -- six
mask panels with the three max MAF series and the Beta order-statistic 95%
band, a seventh with the six masks at max MAF 0.001 plus FlexRV pooled --
one PNG per pair at 110 dpi under `summary_v4/qq/<anc>/`, bulk thinned to
~1,500 ranks below expected 2 with every tail point kept. The four
power/calibration figures (stacked counts below a threshold, per-trait
counts, -log10 scatter, QQ by trait type) read the first table. sas verified
end to end; the other four follow their jobs. Notebook tidy-up:
`tidy_04b.py` (scratchpad of this session; a builder that copies the kept
cells verbatim by marker) produced `04b_production_tidy.ipynb`.

## universal-saige (BRaVa-genetics) moved to this build (2026-09-28)

Branch `saige-slim-aou-parameters` in `~/Repositories/universal-saige`, one
commit `e5c2223`, pushed; draft PR https://github.com/BRaVa-genetics/universal-saige/pull/7 awaiting the owner's review and merge: the analyst pipeline now pulls
`astheeggeggs/saige-slim:<tag>` from Docker Hub (the namespace is a
placeholder until the owner's Docker Hub push; `SAIGE_IMAGE`/`SAIGE_VERSION`
in its `download_resources.sh`), refuses VCF (convert once with plink2),
takes PLINK 2 as the recommended input, carries the AoU step-2 choices and
`--tol 0.02`, ADDS `--useSparseGRMforVarRatio=TRUE` to step 1 (without it
step 2 given the GRM refuses the model on this build), runs FlexRV via
`--flexRVscore NAME` with `04_flexrv_groupfile.sh` around a copy of
`harness/tools/flexrv_score_from_alphamissense.py`, and passes the fit gates
through. Checked end to end on this repo's fixture with the local image:
gate refusal control, both step-1 arms byte-identical, variant / SKAT-O /
FlexRV runs, PGEN vs PLINK FlexRV byte-identical. Its
`docs/saige-slim-migration.md` records the decisions. Found on the way:
LEDGER #171 (a refused step-1 fit leaves 0-byte outputs). The three
`harness/tools/{flexrv_score_from_alphamissense.py,build_flexrv_groupfiles.sh,check_groupfile_overlap.py}`
here are still UNTRACKED (another session's); universal-saige holds a copy of
the first.

**Docker Hub (2026-09-28):** `astheeggeggs/saige-slim:1.5.2-dev-3e92d89d` pushed, public, digest identical to the gcr copy; universal-saige pins it. Registries and push commands: `docker/readme.md`.

**LEDGER #172 (2026-09-28, OPEN): three FlexRV cells with a CONVERGED but
impossible saddlepoint p.** Found on the FlexRV-vs-SKAT-O QQ grids (amr/eas):
p 1e-53..1e-93 with `Is.SPA true`, one of them 91 orders below its NORMAL p at
|Z| = 2.3. The rows are the v4 image's: the scan flags only failures, so a
converged wrong value was never rerun, and every kernel change since
`f2fb24da` is on a failure path -- the current image is predicted to reproduce
it. Five synthetic shapes through both kernels give sane answers, so the real
vector is needed: `SAIGE_SPA_DUMP_BELOW=<p>` now dumps converged saddlepoints
below `<p>`, every dump carries the kernel/roots/remainder line and an `nz`
carrier flag, `aou_spa_dumpcheck.sh` replays the run's own kernel and prints a
Monte Carlo reference (`aou_spa_rerun.md` section 8; test
`test-step2-spa-dump.R`). Decision recorded on the row: do NOT rerun every
SPA cell until the bug is found; binary traits only afterwards; the
egress-safe lower bound is the count of rows with `log10(p.value.NA) -
log10(p.value) > 3`.
Commit `f14d47e5` (gate green, pushed); image **IMG8**
`gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f14d47e5` = Docker Hub
`astheeggeggs/saige-slim:1.5.2-dev-f14d47e5`. Next action is the owner's: `submit_below` for the
three (trait, gene) pairs, then `aou_spa_dumpcheck.sh` on each `.spaDump.tgz`.

**#172 MECHANISM (2026-09-29, from the replay of the amr dumps; owner's
decision on the fix pending).** Not the kernel. `q = Tstat/sqrt(var1/var2) +
m1` rescales the score by the band's variance ratio; at 0.9893 / 0.9907 a raw
score of 0.998 (the weight-1 LoF atom, MAC <20, observed in cases) becomes 1.003, ABOVE
`gmax` = 1, and the saddlepoint charges the remainder (sd 3.6e-4) 13 sigma for
the excess: `t*q-K` 115 where ~6 is right, both kernels agreeing. Planted on a
synthetic vector: VR 1.0 -> 1.4e-2, 1.01 -> 2.7e-3, 0.9893 -> 0, 0.98 -> 4e-46
(reference 3.6e-3); `test-unit-spa.R` #172 block characterizes it. Options
laid out for the owner on the row: route atom-dominated sets to the exact test
by the weighted effective carrier count (recommended), refuse at the crossing,
or rescale the CGF's support. The eas MatHem cell is NOT yet explained: on IMG8
its `AM_sbeta10_MAC5` prints 0.86 with its MAC above the exact-test cutoff, so
something changed a saddlepoint that v4 accepted -- its singleAssoc row and
fallback sidecar are the next read.

**Correction (2026-09-29):** the eas pair is MatHem x `ENSG00000128607`; an earlier note carried `ENSG00000136802`, a gene id copied from a splice-count line, and the first IMG8 eas run tested that wrong gene (it reproduced v4 trivially, empty dump). Every mention above is corrected; the eas rerun on the right gene is the next step.
