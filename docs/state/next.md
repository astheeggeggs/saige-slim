# What is worth doing next (2026-09-03, updated 2026-09-04)

**Start with `docs/state/handoff-2026-09-04.md`** -- it says what landed, which
gates passed with their numbers, where CI lives, and what is a fresh start.

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
* the C++ region loop is the other 43.9%; threading it would compose with the
  above (serial fraction ~11% if both), but `qfc` thread-safety is unchecked and
  the SKAT-O side would have to be reimplemented to run in-process. The brief
  proposes doing **BURDEN only** first -- 93.7% C++, 3.36x at 4 threads, and
  none of the pipelining tension.

## From the 2026-09-04 upstream-issue-tracker mine (LEDGER #63-#85)

The 186-issue upstream tracker was mined for defects live on this branch; the
plan and the raw-evidence pointers are in
`~/.claude/plans/can-you-take-a-dreamy-breeze.md` (**do not re-run the sweep** --
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

**Phase 2 is STARTED; Phase 3 (fixes) is not.** The plan orders Phase 2 by gate
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

Remaining Phase 2 items, in the plan's order:
2. a **constructed high-LD conditioning pair** -- the only way to exercise the
   conditional SPA at non-negligible conditioning at all; reaches #66 and #84.
3. (dropped -- the empty ultra-rare collapse was refuted, see above)
4. a **non-autosomal step-1 fixture**, made by rewriting one `.bim` column (the
   `.bed` is indexed by row, so no regeneration) -- reaches #75.
5. an **all-numeric-ID bim** -- reaches #74.
6. a **phenotype file with a duplicated ID** -- reaches #67.

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

## Outstanding DECISION for the user, not a task to pick up

[LEDGER #57, the part left open] `get_SKAT_pvalue`'s `try(..., silent = TRUE)`
(`R/SAIGE_SPATest_Region_Func.R:363`) turns any SKAT failure into
`Pvalue <- c(NA, NA, NA)` with an `error.code <- 2` that is assigned and never
read, printed or returned. The #57 fix broke the I/O path into it -- a failed
chunk-file read now stops the run instead of arriving here -- but a genuine
NUMERICAL failure (a singular or near-singular Phi, say) still turns a whole gene
into NA rows at exit 0 with no signal.

Three defensible answers and they are materially different operationally: warn
and keep the NA, abort the scan, or carry a count and print an end-of-run summary
line. Whichever is chosen changes what a long exome run does when one gene
misbehaves, so it was NOT defaulted. Ask before implementing.

## Items

**Format step-2 output rows on the workers.**
[Now MORE valuable than when this was written: a multi-model scan formats P rows
per marker on the serial writer, so the Amdahl serial fraction under
`--nThreads` rises with the number of models -- see
`docs/state/multimodel-step2.md`.]
[From the 2026-09-03 assessment; NOT independently measured.] `testOneMarker`
would build the row into a `std::string` on `Marker::MarkerResult`, leaving
`writeOutfile_single` (`src/Main.cpp:3072`, called at `:811` after the
`parallelFor`) as one `<<` per row on the main thread. The predicted effect is
on the serial fraction measured above (s = 0.305 at N = 1,000), so it should
shrink as the writer stops formatting. Byte-exactness rests on the
`std::ostringstream` carrying the same stream state as `g_out.single`
(`src/Main.cpp:89`, where no precision or flags are set) -- check that rather
than assume it, and gate on the marker goldens P1/P4/P7/P11/P12/P27-P33 plus
their PGEN mirrors.
`[ASSUMPTION: the prediction of a specific speedup is the assessment's
arithmetic from NOTES.md's Amdahl fit, not a measurement of this change.]`

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
