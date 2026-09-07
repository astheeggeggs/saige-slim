# What is worth doing next (2026-09-03, updated 2026-09-07)

**Current handoff: `docs/state/handoff-2026-09-07-overnight.md`** (then this
file, then `python3 harness/tools/ledger.py open`). Handoffs are dated and
never renamed; when several share the newest date `ls | tail -1` sorts them by
name, not by time, so this line is the pointer and is updated with every new
handoff. `docs/state/handoff-2026-09-04.md` is where CI lives.

Short list, kept deliberately short. `LEDGER.md` says what is broken;
`docs/state/*.md` says what was done and why; this file says what is NOT done
and is worth keeping in mind. It is not a plan and has no phases.

**Provenance is marked on every item, because it varies.** Two came from a
read-only assessment written by another session on 2026-09-03 and are recorded
here because they are expensive to rediscover -- neither is backed by a
measurement taken for it. Two came out of the S12a / LEDGER #50 work in the same
session. Nothing here has been profiled at realistic N, which is itself the
point: see "the gate" below.

## Resolved: the four known-wrong goldens WERE re-baselined (2026-09-04)

[Grounded; LEDGER #53/#56] Both step-0 sparse GRMs and both `binary_denseGRM`
step-1 fits were regenerated after the `geno2` fix and committed. Each was
proved byte-reproducible first, as CLAUDE.md requires: two independent
`make_goldens.sh` runs into separate candidate trees, `cmp` on all six files.
`compare_goldens.sh` is back to "all comparisons PASS" and the standard tier is
green again. Nothing else moved -- every step-2 golden and every sparse-GRM
step-1 golden was byte-identical throughout.

## The gate on all of it: profile at realistic N first -- NOW DONE

Every step-2 profile so far is N = 1,000, where `harness/bench/NOTES.md:352-359`
shows the serial writer -- not the work -- sets the ~2x ceiling at 4 threads.
That ranking is not expected to hold at biobank N, so choosing among the items
below before measuring is guessing. The sim fixtures are the instrument
(`harness/sim/`; the N=50k `bench50k` build lives outside the repo).

**MEASURED 2026-09-03, so this section is now evidence rather than a worry.**
See `harness/bench/NOTES.md` (the group-test section) and LEDGER #51. At
N = 50,000 over 210 genes: `--is_fastTest=TRUE` costs **+16.3%** and returns
BYTE-IDENTICAL output, because the region path never reads the flag. Escalation
rates run 2.9% (`y_null`, against a 5% nominal -- a clean negative control) to
14.8% on signal arms.

**An earlier draft of this file said an underpowered fixture would understate
group cost "by up to ~2x". That was wrong as stated** and is corrected here: the
2x applies per ESCALATED gene, but at a ~10% escalation rate the aggregate bias
is ~10%. Still worth reporting the rate alongside any per-gene cost -- it is
signal-dependent across ~5x -- but it is not the large effect claimed.

**The measurement moved the target, and then the profile moved it again.**
Swapping the PCG score test for `scoreTestFast` across a whole pass saves only
**8.3%**, so ~90% of gene-test time is NOT the per-marker score test. That was
written as an open question -- ultra-rare collapse, P1Mat/P2Mat assembly, or
R-side SKAT-O post-processing? **Answered 2026-09-04 by `Rprof` at N = 50,000**
(`harness/bench/NOTES.md`, "Where the OTHER 90% of a gene test goes"):

| bucket | share |
|---|---|
| C++ `mainRegionInCPP` | **43.9%** |
| R `region_postprocess` | **53.8%** |
| .. of which `SKAT_davies` under `integrate()` | **44.9%** |
| .. R glue | 5.0% |

So the largest single item in a SKAT-O group run is the SKAT package's Davies
routine under R adaptive quadrature, and **threading the C++ region loop caps at
1.49x on 4 threads** for SKAT-O. For BURDEN (`--r.corr=1`) the same fixture is
93.7% C++ and the same threading would give 3.36x. Replicates agree to 0.3 pp and
the BURDEN arm is the negative control -- it predicts `SKAT_davies` = 0 and
measures 0.0%.

**DONE (2026-09-04): the cells are evaluated in parallel.** `R/SkatoParallel.R`
drives `region_postprocess` three times -- collect the cells, evaluate them on a
persistent fork pool, serve the answers -- so its 175-line triple-nested loop did
not have to be restructured. `--nThreads` on a group run drives it (it used to be
discarded there). **1.30x at N = 50,000, byte-identical**, 8 workers no better
than 4. The remaining ceiling is 1.43x, so most of the available win is taken.

Still open on this path, in rough order of value. **Both now have a standalone
brief: `docs/state/region-perf-handoff.md`**, which carries the profile, the
Amdahl table, the two traps that each cost a measurement cycle, what was
measured and REJECTED so it is not retried, and the one design question that
must be answered before item 2 below is coded.
* per-task overhead still costs the gap between 1.30x and 1.43x; batching cells
  into fewer, larger tasks would narrow it and has not been tried. Small, low
  risk, and the gate already exists.
* [DECIDED 2026-09-05, LEDGER #51] `--is_fastTest` on group runs is pure cost
  (+16.3% at N = 50,000, byte-identical output) and stays so by decision: the
  region path never reads the flag, the R side still escalates, and a warning
  says so. The potential speed-up is bounded at the **8.3%** measured for
  `scoreTestFast` throughout, and taking it would move every fastTest region
  golden (P2, P3, P8, P8b, P10, P14, P15, P21, P23, P24, P34, P35) and retire
  P35's reason to exist. Only worth revisiting once the two items above are
  taken and the region loop is the bottleneck again.
* the C++ region loop is the other 43.9%; threading it would compose with the
  above (serial fraction ~11% if both), but `qfc` thread-safety is unchecked and
  the SKAT-O side would have to be reimplemented to run in-process. The brief
  proposes doing **BURDEN only** first -- 93.7% C++, 3.36x at 4 threads, and
  none of the pipelining tension.

## From the 2026-09-04 upstream-issue-tracker mine (LEDGER #63-#85)

The 186-issue upstream tracker was mined for defects live on this branch; the
plan and the raw-evidence pointers are in
`docs/state/plan-tracker-mine.md` (**do not re-run the sweep** --
it cost ~160 agents). Phase 1 is DONE: **23 entries, LEDGER #63-#85**, one row
each, no fixes. Every line number was re-resolved against HEAD and every
arithmetic claim re-derived before the row was written, which is what the plan's
Tier 0 retraction was about.

**#63 is the worked example and is the only one with a test**:
`tests/testthat/test-characterization-ledger63.R` (standard tier), detection
power drilled BOTH ways -- 20/20 at HEAD, 8 failures with the `- 1` removed and
the package rebuilt into a scratch library.

**Reproduced by hand this session, not merely traced**: #63 (the dropped GRM
pair, with the committed `.mtx` reproduced exactly as the positive control),
#68 (`x1^2` silently dropped from a formula), #72 (the `std::stod` boundary is
DBL_MIN, i.e. every subnormal -- not the reported 1e-324), #73 (conditioning on
a monomorphic marker: exit 0, and the `_c` columns identical to the marginal
ones on 70/70 rows, against 67/70 changed for a real conditioning marker), #74
(the data.table join-type message reproduced verbatim), #76 (a region run at
`--chrom=2` writes a ZERO-BYTE results file at exit 0 and says "Analysis
done!"), #80/#81/#82/#84 (verified against the committed goldens and `git
blame`). The rest are traced-at-HEAD with the consequence explicitly marked as
not reproduced -- read each row's PROVENANCE, it says which.

**Four candidates were REFUTED during this pass and deliberately have no row.**
Recorded here so nobody re-opens them:
- **empty ultra-rare collapse** (upstream #12/#30/#60/#67) -- the C++ handles it
  coherently, checked at four points: the pseudo-marker slot keeps its NA/NaN
  initialisation, `annoMAFIndicatorMat.row(i)` stays all-zero so the cell is
  excluded from the region test, `i1`/`i1InChunk` are NOT incremented so VarMat
  stays square and aligned with `mPassCVVec`, and `writeOutfile_singleInGroup`
  skips any row whose p is `"NA"`. If a defect lives here it is in the
  `noNAIndices`-vs-`VarMat` alignment in `region_postprocess`, not in the
  collapse.
- **sparse-GRM nnz overflow** (upstream #188) -- the proposed half-to-full
  storage doubling does not happen: in the pinned env (Matrix 1.7.6) a
  `dsCMatrix` survives `*`, `diag<-` and `[i, i]` with half storage intact
  (checked). What remains is only that `Matrix:::readMM` is called unguarded at
  `R/SAIGE_fitGLMM_fast.R:2109` and `:2160`, with no mechanism located.
- **SKAT-O compact-support p = 1** (upstream #209) -- the outer integral is
  inside `SKAT:::Met_SKAT_Get_Pvalue`, a dependency, not in this tree.
- **`null_noXadj` cardinality check in the wrong arm** -- this is already
  **LEDGER #39**, which names the unreachable `stop()` explicitly. #39's line
  numbers are stale by ~17 lines (`:405-414` -> `:388-399`).

**Phase 2 is COMPLETE; Phase 3 (fixes) has STARTED** -- #74 and #66 are done
(`22be2a69`, `c9f7e9ae`), each with its characterization test INVERTED in the
same commit. The plan orders Phase 2 by gate
gap rather than severity. Item 1 is DONE: goldens **P36/P37**
(`plink_group_burden_quant.txt`, `plink_group_nofast_quant.txt`) close the
QUANTITATIVE region-path gap -- every previous group golden is binary -- and
`tests/testthat/test-step2-group-quant.R` drives them. Both were proved
byte-reproducible over three independent runs before being accepted.

**That golden immediately paid for itself: LEDGER #86.** P37 reports SKAT-O
`Pvalue` **exactly 1** on 5 of its 12 cells while their own `Pvalue_Burden` and
`Pvalue_SKAT` are 0.82-0.99, which the suite's existing invariant
(`test-invariants-goldens.R:106-110`, SKAT-O lies between its component tests)
forbids -- it was asserted on one BINARY golden only. Measured on the kernel
directly: `SKAT:::Met_SKAT_Get_Pvalue` returns exactly 1 for 55 of 400 random
inputs, all with `min(p.val.each) > 0.777`, and never below 0.763. And it is not
inert -- both CCT implementations fall back to `min(1, n*min(p))` when any input
is 1, so one saturated cell silently turns the region's `Cauchy` row into a
Bonferroni (confirmed to the last digit on GENE2). **This corrects the
refutation recorded above**: the sweep's SKAT-O finding was dismissed as living
in the `SKAT` dependency and therefore unlocatable. The integral does live
there, but the consequence is here and is now pinned.

**Item 2 is DONE**: `tests/testthat/test-step2-cond-highld.R` builds a
conditioning fixture AT RUN TIME (seeded, no golden -- the genotypes come from
R's RNG, so a golden would pin the RNG rather than the code) over the same
10,000 fam rows as `genotype_100markers`, so the committed
`binary_sparseGRM_a` model is reused unchanged. It carries an LD ladder against
one conditioning marker (r^2 = 1.000 / 0.90 / 0.49 / 0.000, MEASURED on the
1,000 analysed samples inside the test) plus a marker built to correlate 0.35
with the null model's own residual, so a large conditional statistic survives
conditioning.

Both target regimes are reached, and both controls behave: at r^2 ~ 0
conditioning is a no-op (var/var_c = 1.001), at r^2 ~ 0.5 the variance really
drops (1.89, against 1/(1-r^2) = 2.0).
- **#84**: the associated marker gives `stat_c = 91` against a gate of 4 with
  `var/var_c = 1.51`, and `p.value_c` (2.8e-18) differs from `p.value.NA_c`
  (1.3e-21) by three decades -- the conditional SPA block executes at real
  conditioning for the first time on this branch. DRILL: on golden P12, **zero**
  of 70 rows satisfy both `stat_c > 4` and `var/var_c > 1.3` (its three
  high-ratio rows are the conditioning markers themselves, `stat_c` ~ 0.003).
- **#66 -- FIXED in Phase 3 (2026-09-05, `c9f7e9ae`); what follows is the
  measurement that located it.** Under `--is_noadjCov=FALSE` a perfectly
  collinear marker gives
  `var_c = 2.5e-14` -- cancelled, and yet ~294 orders ABOVE the `DBL_MIN` guard,
  which is exactly the geometry LEDGER #44 had to fix with a relative test.
  `BETA_c` and `SE_c` are then ratios of two cancelled quantities (SE_c = 6.3e6)
  and `p.value_c` is forced to 1. DRILL: adding a relative cancellation guard to
  `src/SAIGE_test.cpp` and rebuilding into a scratch library fails 6 assertions;
  source reverted and reinstalled.

**SETTLED 2026-09-05 (recorded in LEDGER #82).** At the CLI default
(`--is_noadjCov=TRUE`) the same perfectly collinear marker reports a non-zero
`var_c` because the tested marker's variance is the UNADJUSTED, allele-frequency-
centred kernel while the conditioning cross-term is built from the covariate-
ADJUSTED marker. The residual is exactly `vr_noXadj * [centring-mismatch term +
W-weighted covariate term]`: predicted from the model to every printed digit
(0.362754 on the run-time fixture) and, from the two runs' printed columns alone,
`var(def) - (vr_noXadj/vr_null) * var(adj)` to printing precision. The earlier
`var * R^2_OLS` prediction was low because it dropped the centring term (30% of
the residual) and used an unweighted R^2. Pinned in `test-step2-cond-highld.R`.

**Item 4 is DONE**: `tests/testthat/test-step1-varratio-cursor.R`. The
variance-ratio marker loop advances `indexInMarkerList` once unconditionally and
a SECOND time on each skip branch, so a skipped marker consumes the marker after
it, which is never read at all -- and because the exhaustion guard tests
EQUALITY, a double step from the last position jumps over it and the next read is
an NA.

Reproduced by hand, and the arithmetic reconciles two ways. With a run-time
fixture of M markers at MAC 40 (all in the single VR pool) over the step-1
fixture's own 1,000 fam rows:
- NEGATIVE CONTROL, all autosomal: exit 0, exactly 30 markers read (the
  `--numRandomMarkerforVarianceRatio` default), a valid ratio file.
- **All non-autosomal, EVEN pool (M = 40)**: exactly **20 of 40** markers read --
  the cursor visits 1, 3, 5, ..., which is the arithmetic signature; a
  single-increment loop reads all 40. Lands exactly on the guard, exits the loop
  cleanly, then dies downstream in `calCV(NULL)` with no marker accepted.
- **All non-autosomal, ODD pool (M = 41)**: 22 reads, and the failure MODE
  changes -- `Error in if (genoInd == 0) { : missing value where TRUE/FALSE
  needed`, which is #75's predicted abort and upstream #32/#151 verbatim. **One
  marker of difference flips the outcome**, which no single-increment loop can do.
- Both leave a 139 KB `.rda` beside a **0-byte** `.varianceRatio.txt` -- the
  unusable model pair #75 describes.
- The SILENT case is the partial one: 20 of 40 non-autosomal gives **exit 0**
  with the ratio estimated from 12 markers instead of 30, and 14 pool markers
  never examined. Reconciled: 12 accepted + 14 skipped = 26 reads, cursor
  = 1 + 12 + 2*14 = 41, exactly the guard value for a pool of 40.

DRILL: removing the second increment and rebuilding into a scratch library fails
5 assertions (read counts become 40 and 41, and the NA abort disappears); the
negative control correctly stays green. Source reverted and reinstalled.

**Item 5 is DONE**, and **#74 is now FIXED** (Phase 3, `22be2a69`; the test
below was inverted with it, and a pgen arm added because the pvar half of the
fix was otherwise ungated). What follows is the characterization as first
written. The bim WAS
read with no `colClasses`, so data.table INFERS the ID column's type; every
bundled fixture uses `rs<N>` (character), and an all-numeric ID column infers
INTEGER. The fixture is derived at run time and reuses the committed `.bed` and
`.fam` UNCHANGED -- a `.bed` is indexed by row, so rewriting bim column 2 changes
no genotype, which is what makes the control exact.

- NEGATIVE CONTROL: a plain marker scan is UNAFFECTED -- same row count, and
  `p.value` identical string-for-string between the two bims. The ID column is
  not joined on that path.
- `--condition=13` -> exit 1, `Incompatible join types: x.ID (integer) and
  i.SNP (character)`, upstream #137's message verbatim.
- a group file with numeric IDs -> exit 1, the same `bmerge` failure.
- POSITIVE CONTROL: both commands with the committed rs-spelled bim exit 0.

DRILL: adding `colClasses = list(character = 2)` to the `fread` and rebuilding
into a scratch library makes both numeric-ID runs SUCCEED, failing 5 assertions;
the plain-scan control stays green. That one argument is also the Phase 3 fix.
Source reverted and reinstalled. **Phase 3 then took exactly that argument**,
plus `list(character = 3)` on the pvar, and strengthened the positive control:
the numeric-ID and rs-spelled runs must now agree on EVERY column but
`MarkerID`, which they do.

**Item 6 is DONE**: `tests/testthat/test-step1-duplicate-pheno-id.R`. Nothing in
the step-1 phenotype path checks for duplicate sample IDs -- the only
`duplicated()` call in `R/SAIGE_fitGLMM_fast.R` is on `SampleIDIncludeFile` -- so
a repeated ID multiplies that sample's rows through the merge at `:1012-1014`.

A phenotype file of 1,001 rows with 1,000 distinct IIDs runs to **exit 0** and
reports "1001 samples will be used for analysis"; the unmodified file reports
1000 (negative control).

The consequence was DERIVED FROM THE CODE BEFORE IT WAS MEASURED, which is what
makes it a test rather than a fishing expedition. In `setGenoObjFromSource`,
`Nnomissing = subSampleInGeno.size()` counts the DUPLICATED list, while the
allele-count accumulation is gated on `indicatorGenoSamplesWithPheno_in[ind]`
with `ind` running over the FILE's samples, so each distinct sample is counted at
most once. `altFreq = alleleCount / float(Nnomissing * 2)`
(`src/SAIGE_fitGLMM_fast.cpp:677`) therefore puts the duplicated sample in the
DENOMINATOR but not the NUMERATOR. Predicted gap: exactly that sample's genotype.
Measured, on a fixture where the duplicated sample carries the ALT allele in half
the markers: store holds 41 alleles and the frequency implies 40 where it is a
carrier, 39 vs 39 where it is not -- `store_AC - reported_AC == g_dup` on all 40
markers to float32 precision. Both arms in one run, so the probe cannot be
fabricating the signal.

DRILL, two patches in one scratch build because the two tests pin different code:
deduplicating the merge in R breaks the sample-count assertion, and changing the
`altFreq` divisor to the indicator count breaks the three frequency-gap
assertions. 4 failures, exactly the intended ones. Both reverted and reinstalled.

**Phase 2 is COMPLETE.** See the handoff for the gate numbers.

## From the 2026-09-04 silent-failure audit

A 26-agent fan-out hunted the bug class this branch keeps finding by accident
(silent failure that every golden reproduces). 15 raw findings -> 11 unique ->
10 adversarially verified by two skeptics each, both instructed to default to
REFUTED. Five survived; five were killed, including a plausible `raw_variant_ct`
claim and a `dosage_zerod` one.

**Two were reproduced by hand the same day and are FIXED (LEDGER #58, #59)** -- the
unbounded mmap `.bed` read (a `.bim` claiming one more marker than the `.bed`
holds gave byte-identical output at exit 0, while the pread path stopped), and a
region results file written with data rows and no header. Both have positive and
negative controls in the ledger entries, and both are gated by
`tests/testthat/test-step2-file-integrity.R`.

**The other three were reproduced on 2026-09-04 and are now FIXED
(LEDGER #60, #61, #62).** None was dropped; all three were real, and two were
worse than the finding had claimed. Each carries a positive and a negative
control in its ledger entry:
* **#60**, `src/SAIGE_fitGLMM_fast.cpp` (`BedByteSource::fill`) -- the step-1
  twin of #58. `.bed` truncated to 60 of 100 markers, `.bim` left at 100: exit
  0, no warning, GRM from 18 QC-passing markers instead of 33. Truncating
  MID-record is worse than the finding predicted -- half a real marker plus
  fabricated homozygotes PASSED QC and entered the GRM, a marker the intact file
  rejects. Fixed with #58's size check plus a `gcount` check in `fill()`.
* **#61**, `R/SAIGE_Test_main.R:602` -- and the effect is not "all its markers
  discarded" as the finding read it, it is that the bin comes out a
  byte-for-byte DUPLICATE of the bin below, which looks like an answer.
  `--maxMAC_in_groupTest=100` on the 1,000-sample fixture (implied MAF 0.05
  against a largest cutoff of 0.01) gives `Number_rare` 3/4/4/1/3/3 where the
  same bin asked as a MAF gives 4/7/7/2/7/7, and moves a p-value 3x. Fixed by
  re-issuing `setRegion_GlobalVarsInCPP` with the extended vector.
* **#62**, `R/Util.R` (`updateChrStartEndIndexVec`) -- confirmed exactly as
  described: one chr1 marker inside chr2 passes both arms of the guard and
  chr1's window then spans 16 markers for its 11, dropping five chr2 markers
  from both LOCO GRMs. Fixed by asserting contiguity, which is the invariant the
  windows actually need.

Detection power drill-verified in both directions for all three: the new tests
are red at unpatched HEAD (1 / 5 / 4 failures) and green after.

## From the strip-vs-upstream comparison at N = 50,000 (2026-09-05)

[Measured; `harness/sim/upstream_compare/RESULTS.md`] Upstream `main` (365deeff)
built beside `strip` in its own pixi env pinned to the same R 4.4.3 / OpenBLAS
0.3.30 builds, both arms on `sim50k`, paired and alternated.

- **LEDGER #97 -- FIXED 2026-09-05** (`R/SAIGE_fitGLMM_fast.R`, the `sparse` rows
  of `varRatioTable` no longer multiplied by `varRatioScale`). The #65 rescale had
  multiplied the scale-free `sparse` variance-ratio rows by 1/s^2, so every
  quantitative model fitted with `--useSparseGRMforVarRatio=TRUE` between
  `66a65bb9` and the fix carried `sparse` rows 11-12% high on the N = 50,000
  fixture (1/s^2 ~ 1.5 on the bundled one); refit those models. Gate:
  `tests/testthat/test-step1-vr-sparse-scale.R` (standard tier) -- a quantitative
  step-1 with `--useSparseGRMforVarRatio=TRUE` whose `sparse` rows are invariant
  to y * 100 while its `null` rows scale by 1e-4, and whose `sparse` rows are ~1.
  No committed golden moved (none has quantitative `sparse` rows).
- **RESOLVED 2026-09-05 (evening): the binary single-variant scan's 0.94x was not
  like-for-like, and its kernels were profiled and de-allocated.**
  (`docs/state/handoff-2026-09-05-binary-profile.md`; `harness/bench/NOTES.md`.)
  Firth is 56% of the binary scan's CPU on both arms; upstream's #28 made it skip
  ~1,000 Firth fits and ~450 SPA corrections on AF > 0.5. A per-worker
  `KernelWorkspace` (byte-identical at 1 and 4 threads, standard tier green) took
  the 1-thread user CPU to 205 s (upstream 236) and the 4-thread sys time from
  55 s to 11 s. `--firth_solver=closed` (-27% more, output byte-identical, P31 unmoved) is
  now the DEFAULT with `qr` selectable; the warm start (118 non-converged fits)
  was benchmarked and deleted. The PLINK reader now
  decides QC gate 1 from 2-bit code tallies before expanding a marker
  (quantitative 50k scan 62 -> 45 s, byte-identical).
- Strip's step-2 peak RSS at N = 50,000 reads 2.3-2.8 GB against upstream's
  0.35-0.38 GB. Measured to be resident pages of the mmap'd 2.6 GB `.bed`, not
  allocations: `SAIGE_PLINK_NO_MMAP=1` gives 312 MB with byte-identical output.
  Worth a sentence in the user docs, not a fix.

## DECISIONS for the owner (2026-09-05): all twelve answered

The eleven rows the Phase 3 sweep left open, plus #18's second half, were
each decided and implemented the same day
(`docs/state/handoff-2026-09-05-decisions.md`). `python3 harness/tools/ledger.py
open` lists what is still open (documentation rot, #92/#93/#95/#96, and the
Sigma floors #91); `ledger.py decided` the four left-as-is rows.

## Time the gate at the next run (2026-09-05)

The standard tier took ~4 min 20 s per compare mode on the last pre-push
(`harness/work/prepush-{byte,numeric}.log` mtimes 21:28:34 and 21:32:54, on a
machine still at load 10-14 from the N = 50,000 gate), so `prepush.sh` is ~9
min against the ~5 CLAUDE.md implies. Not yet measured idle. What this session
added to the standard tier, all step-1 fits (wall time is bounded by the
slowest FILE): `test-step1-nglmm-prep.R` (two extractNglmm runs, ~20 s each),
`test-step1-lowmem-loco.R` (four 3 s fits), `test-step1-stale-model-file.R`
(three fits), `test-step1-vr-sparse-scale.R` (two ~5 s fits), and the
morning's `test-step1-pheno-scale.R` (three fits + three step-2 runs).
At the next gate: check `uptime` first, run the tier once with
`reporter = "progress"` and collect the per-file `[Ns]` markers (the parallel
result object carries no timings), then split or demote the long poles
(`nglmm-prep` is the likely one; its second arm could share the first's run).
Also worth asking whether numeric mode needs to re-run every step-1 fit, since
the two modes differ only in how step-2 tables are compared.

## From the 2026-09-05 code-quality review of the Phase 3 sweep

[`docs/state/handoff-2026-09-05-review.md`; measured on the two-chromosome
dense-GRM fixture that `tests/testthat/helper-pheno-scale.R` derives at run time]
The review's refactors landed (LEDGER #97-#100 fixed on the way; the fold of
`fitNULLGLMM`'s binary/quantitative halves proved byte-identical by `cmp` of
every model file against the pre-fold build). Three things it made visible and
deliberately did NOT change:

- **Binary LOCO refits are start-dependent, and the two layouts start
  differently.** `glmmkin.ai_PCG_Rcpp_Binary`'s LOCO loop passes `alpha`/`eta`
  that it overwrites with each chromosome's result, so chromosome j warm-starts
  from chromosome j-1's refit; `fitNULLGLMM`'s `--isLowMemLOCO` loop starts every
  chromosome from the full model's `alpha0`/`eta0`. At `tol.coef = 0.02` they stop
  at different points: chr1 agrees to 1e-10, chr2 to ~1e-5 (coefficients 8.5e-6,
  fitted values 2.4e-5, cov 7.7e-6) and the LOCO entry's `obj.noK$S_a`, the score
  at the refit, not at all. The quantitative refit is one linear step and agrees
  to 1e-10 on both. Pinned as found in `test-step1-lowmem-loco.R`; which start
  upstream intends is not known, and either choice moves every binary LOCO model.
  [ASSUMPTION: the warm start is the whole cause -- derived from the code, and the
  chr1-exact/chr2-loose pattern matches it; not isolated by forcing the start.]
- **Two by-trait asymmetries in the low-memory refit, preserved.** The binary
  refit passes the glm's offset (none unless `isCovariateOffset`) and
  `tol.coef = tol`; the quantitative refit passes the model's covariate offset
  (`X[, -1] %*% coef[-1]`, alongside the same X and alpha) and the default
  `tol.coef = 0.1`. Both mirror what each trait's fitter does in its own LOCO
  loop, so the layouts agree per trait; whether the quantitative offset is
  double-counted is a derivation nobody has done. Marked "by trait" in the code.
- **The N = 50,000 own-model quantitative tables are a mixture and want one
  re-run.** The review's gates refitted the strip own models (sparse rows now
  1.00, the #97 fix) but did not re-run own-model step 2, so
  `harness/sim/upstream_compare/RESULTS.md`'s own-model `y_quant` section
  describes outputs from the pre-#97 models. From
  `harness/sim/upstream_compare/`: `MODELS=own bash run_arm.sh strip s2single
  y_quant 1`, the same for `s2group`, then `bash compare_all.sh 1 1`; ~10 min.
  Expected: the own `gene y_quant` tags lose their `VR`/`VR-mix` mass.
- **`#65`'s save-then-rescale is a design the owner may still change** (review
  finding 1, second half): the run saves the standardised quantitative model,
  reloads it for the ratio estimation, and rescales the files it wrote at the
  end (`savedModelFiles`, #98). The alternative -- save in phenotype units at the
  three save sites, standardise on the reload -- deletes `rescaleSavedQuantModels`
  and the CLAUDE.md sentence about standardised files, at the cost that the ratio
  estimation reads a `x s x 1/s` round-tripped model (ULP level), so P4, P16, P32,
  P36, P37 (+`singleAssoc`) and the quantitative step-1 goldens would need a
  second rounding-level re-baseline. Not done; owner's call.

## Resolved DECISION: SKAT kernel failures are counted and summarised

[LEDGER #57, the open half; decided 2026-09-05] `get_SKAT_pvalue`'s
`try(..., silent = TRUE)` still turns a kernel error into NA rows, by choice: the
user picked count-and-summarise over warn-and-NA and abort, so that what drives
the NAs can be determined. The run now prints `SKAT failures: N ...` with a
by-message tally before `Analysis done!` and writes
`<SAIGEOutputFile>.skatFailures.txt` naming each failed cell when N > 0. Result
files are byte-identical. Tests: `test-unit-skat-failures.R`,
`test-step2-skat-failures.R` (planted kernel failure in a child process).

## Items

**Resolved: formatting step-2 output rows on the workers is CLOSED by
measurement (2026-09-06).** [Measured; `harness/bench/NOTES.md`, the writer
section; `docs/state/handoff-2026-09-06-writer-bench.md`] The item predicted
that P rows per marker on the one writer thread cap `--nThreads` as P grows.
At N = 50,000, HEAD `b11e6b07`, the writer is 2-3 us per row and takes
w = 0.2-0.9% of the main thread's `mainMarkerInCPP` time at 1 thread and
0.5-2.1% at 8, on 1, 3 and 5 models alike (worst case 3 quantitative models,
w = 0.021, ceiling 1.02x). The thread-scaling loss the Amdahl fit called a
serial fraction is CPU inflation (1.3-2.7x more user CPU for the same work at
8 threads: allocation churn and hyperthreading on 4 physical cores); the main
thread is 91-95% busy as a worker at T=8. Do not re-open without a fixture
where marker work is under ~20 us (N well below 1,000), which is not the
target regime. Driver: `harness/bench_writer.sh`; reader:
`harness/bench/writer_profile.py`; table: `harness/bench/writer_analyse.py`.

**If phenotype batching is ever attempted, one trap is already known.**
[From the same assessment; recorded because it is cheap to write down and
expensive to hit.] Computing all P traits' score statistics as a single
`G' x [res_1..res_P]` matrix product CHANGES THE REDUCTION ORDER and therefore
breaks the byte gate; the per-trait dot products have to keep their serial
order. This is the same tension as the step-1 GetTrace GEMM decision, which was
resolved the same way. Nothing else about batching is designed, and the natural
unit for it -- a per-trait `TestContext` -- only became possible on the region
path with S12.

**Resolved: `--is_fastTest` warns on group runs, and
`--pval_cutoff_for_fastTest` is a CLI option.** [LEDGER #51 / S12a; done
2026-09-04] The warning sits beside the existing Firth/r.corr one, where
`is_fastTest` is final -- after the conditional path has forced it FALSE. #51
itself stays OPEN: the asymmetry is upstream's and unchanged, only the silence
is fixed. Exposing the cutoff made golden **P35** possible (P2's command with
the cutoff at 1, so both regions take the ESCALATION pass) -- **its bytes are
identical to P2's and that is #51 restated**, not a defect: the second pass
recomputes what the first already had. It pins the other WRITER, and under
LEDGER #50 it lost all 19 single-variant rows while P2 stayed green.

**Resolved: the psam sample map keyed on FID, not IID.** [LEDGER #52, FIXED
2026-09-04] `readPsamFile` now finds the `IID` column by name from the header,
strips `\r`, and refuses a duplicate IID instead of overwriting it. All 13 PGEN
goldens stayed byte-identical, exactly as predicted -- which is why the gate had
to be a fixture: `tests/testthat/test-step2-pgen-psam.R` derives psams at run
time (one FID shared by all 10,000 samples; `#IID SEX` with no FID at all; a
duplicated IID) and drives golden P13 through each. 6 failures at unpatched
HEAD, 0 after.

**PCG stops on an absolute residual, and that is a recorded hazard, not a
fix.** [DECIDED 2026-09-05, LEDGER #6] All eight PCG loops in
`src/SAIGE_fitGLMM_fast.cpp` stop on `sum(r^2) > tolPCG` with no right-hand-side
normalisation, so `--tolPCG=1e-5` is a different relative tolerance at every N
and for every phenotype scale. On the bundled fixtures it only constrains the
harness's step-1 tolerances (theta/VR 1e-2). At biobank N with a DENSE GRM the
stop is relatively looser than at N = 1,000, which is where it could bite; the
sparse-GRM fits never enter a PCG loop at the default flags (direct `spsolve`),
so a normalised rule would move only `binary_denseGRM` and its pgen mirror
[ASSUMPTION: traced, not instrumented]. If revisited: perturb `tolPCG` first to
see which goldens actually move, then normalise by `dot(b, b)` and re-baseline.

**Dosage pgens are hard-called for the GRM, and that is a real limitation, not
just a warning.** [Grounded, 2026-09-03] `genoClass`'s store is 2 bits per
sample, so `PgenByteSource` rounds any dosage to the nearest of {0,1,2} and
warns once with the count. For imputed data that discards information the GRM
and the variance ratios could have used. Carrying dosages properly means a
second storage path through `genoClass` and touches the GRM kernels, LOCO and
every step-1 golden -- deliberately out of scope, and worth measuring (does the
GRM actually move?) before being attempted.

## Deliberately not listed

The consolidation work (README is still the upstream stub, `pixi.lock` is
untracked, bench logs are in git) is real but needs a decision about who this
branch is for, which is not a decision to make from inside the code. The
remaining speedup candidates named in the older plan -- step-1 variance-ratio
estimation blocking, `setgeno` ingest -- are left where they are until the
profile above exists.
