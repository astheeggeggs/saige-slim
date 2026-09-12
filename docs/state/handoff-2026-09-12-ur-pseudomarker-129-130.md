# Handoff 2026-09-12 (night): the ultra-rare pseudo-marker's two defects

`LEDGER #129` and `#130` are FIXED, one commit each, plus a goldens commit
before them and a docs commit after. Read this, then `docs/state/next.md`, then
`python3 harness/tools/ledger.py open` (two rows left: #116 and #128).

Both were at one site -- the ultra-rare pseudo-marker pass in `src/Main.cpp` --
and both were DECISIONS, put to the owner before anything was implemented. The
decisions are the load-bearing part of this handoff.

## 1. The two decisions, and the reasoning that settled them

**#130: never fold the pseudo-marker.** Four options were put up (minor-allele
coding with a dense `sumFill`; the same statistic implemented as a
constant-shifted sparse column; ALT-coding the burden while folding the test;
never folding at all). The recommendation was the first, on the grounds that
every other marker folds and the SKAT-O path's `BETA_Burden` already reads the
folded coding. **The owner reversed it with one question: the annotation is
based on the ALT allele, so what does that do to "rarity" and to the mapping?**
That is decisive, but be precise about WHAT it settles -- the adversarial review
in section 8 caught this row, the commit message and three comments overclaiming
it. What it settles is DIRECTION: `max_m(w_m g_m)` reads "carries at least one
of them", its fold reads "carries NONE of them", and only the first is a burden
of anything, so inverting a collapse is wrong whatever the allele convention.
What it does NOT settle is that the burden counts the ALT allele, because every
constituent marker has already gone through `imputeGenoAndFlip` in the
per-marker prepare (`src/Main.cpp:2746`) and both the collapse and the burden
column are built from its post-fold values. The collapse is over MINOR-coded
dosages. The ALT-vs-minor question for the constituents is the same question one
level down; it is NOT fixed, and it is now the larger half of LEDGER #131.
Never folding also turned out to be the smallest change (the un-folded branch
already existed inline), the only one that keeps the burden column sparse and
non-negative, and the one that makes both region paths agree.

**#129: skip the dosage filter for the pseudo-marker** (fix (b) of the two on
the row, over re-pointing it at the unweighted twin). Its entries are `w*g`, not
dosages; the threshold 0.2 and the gate MAC 10 are dosage and allele-count
quantities. At the defaults the two cutoffs are co-extensive, so every
constituent marker has already been through this same filter as a dosage before
the collapse. The two candidates differ only for a marker whose MAC falls
between `--dosage_zerod_MAC_cutoff` and a RAISED
`--MACCutoff_to_CollapseUltraRare`, where a sub-cutoff entry now survives into
the pseudo-marker. That entry is NOT only mean-imputation noise: PGEN decodes
real dosages (`PgrGet1D` with `dosage_main`, `src/PGEN.cpp:485`, and the
`pgen_marker_dosage*` goldens exist to prove it), so above the defaults this is
a knowing choice rather than a free one. At the defaults it is free, because the
two cutoffs are co-extensive and every constituent was filtered as a dosage
before the collapse.

## 2. What the fixes bought, measured

**#130 gave an exact three-way reconciliation** where there had been a
contradiction. A cell holding ONE marker has, by definition, that marker's
burden, so the region row's `Pvalue_Burden` and the run's own
`.singleAssoc.txt` p are one quantity computed twice -- three times, with the
SKAT-O path's `sum(Score)/sum(Phi)` in R. On the drill (goldens P62/P63) all
six cells now agree on both arms:

| cell | before (BURDEN arm) | after | single-marker p |
|---|---|---|---|
| GENE1 lof (AF 0.374, never folded) | −651.793 / 2.516356E-01 | unchanged | 2.516356E-01 |
| GENE1 missense;lof | +2.34563 / 5.758082E-01 | **−0.860935 / 8.334812E-01** | 8.334812E-01 |
| GENE2 lof | +383.633 / 2.634589E-01 | **−445.642 / 1.663821E-01** | 1.663821E-01 |
| GENE2 missense;lof | +1.86926 / 7.021056E-01 | **−1.76402 / 7.175196E-01** | 7.175196E-01 |

`GENE1 lof` is the in-fixture control: it is byte-identical before and after on
both arms, which is what says the agreement is a property of the statistic and
not of the arithmetic. The two arms' `BETA_Burden` agree to printed precision
(−0.860935 against −0.860934561967558), and both `.singleAssoc.txt` files came
out byte-identical, which is stronger than the "to printed precision" that was
predicted.

**#129's fix lands on numbers the row had recorded before it.** Golden P64
(weights 0.01, collapse 10) now reads `Pvalue_Burden 2.522646E-01` and
`BETA_Burden 232.633`. The row's PROVENANCE had already measured `2.522646E-01`
for the SAME data at weights 0.5 and at collapse 0, and 232.633 is exactly 50x
the weights-0.5 coefficient 4.65266. A burden p cannot depend on the weight
scale; the fixed corner is where that invariance had been broken, and it now
holds over the whole 2x2 of (0.01, 0.5) x (collapse 10, collapse 0).

**The FlexRV fixture's discriminating power stopped depending on the collapse
cutoff.** `test-invariants-goldens.R` records the number of DISTINCT
(Pvalue_Burden, BETA_Burden, SE_Burden) triples per gene among the 192 sets. It
moved 100 -> 111 (GENE1) and 104 -> 111 (GENE2) -- and 111/111 is what the
collapse-OFF golden has always had. The filter had been cleaning several sets'
pseudo-markers down to one shared support: `REVEL_cbrt0.4..0.9_MAC1` all printed
the same `1.453839E-02`. That is the same statement as the row's "collapse 0 and
10 do not even agree on which sets exist", arriving from the other side.

## 3. Two things the rows had wrong, now corrected on them

- **#130 said the SKAT-O path was not affected. It was, by another route.** On a
  BINARY trait `region_postprocess` rescales Phi through `g.sum`, which IS the
  burden column the emitter builds (`get_newPhi_scaleFactor_traitType`,
  `R/SAIGE_SPATest_Region.R`). So that arm moved too, and improved the same way:
  `GENE2 lof` went 0.166847187 -> 0.1663821, its own single-marker p. A
  quantitative trait never reads `g.sum`, so there only the BURDEN arm moves.
- **#130 said the fix was "one line". Carrier space had made that false.**
  `scatterSum` hardcodes fill 0.0 and the UR pass cannot insert carrier slots
  (insertion is serial-phase only; the cell prepares run on the pool), so
  emitting over the post-fold list -- the row's own suggestion -- would have
  needed a `sumFill` twin of `urFill` threaded through six sites. The decision
  made that moot, but a future row at this site should not trust "one line".

## 4. Blast radius: predicted before each gate, confirmed exactly

- **#130**: P62/P63 only. Predicted, and the gate agreed (`FAIL 0 | PASS 4926`,
  7m45s, both comparators on one generation).
- **#129**: **18 files moved, 1 deleted, nothing else** -- predicted file by
  file before the build. The instrument was `--dosage_zerod_cutoff=0` over every
  weight-bearing group golden, an upper bound turned exact by two controls: P59
  (the same markers with the collapse OFF) and the two runs whose smallest
  weight is 0.5 (`withWeights`, `composed`) are byte-identical including
  siblings, so marker-level zeroing is inert on this fixture and every remaining
  difference is the pseudo-marker's. Moved: `plink`/`pgen_group_flexrv_burden`,
  the four multiWeights runs and their `.singleAssoc` siblings (that fixture's
  REVEL set carries weights 0.01-0.07), and P64. Deleted:
  `plink_group_ur_smallweights.txt.pooledTests.txt` -- the region has nothing
  left to disclose, which is the #127 sidecar doing its job in both directions.

Every new golden was proved byte-reproducible by two generations into separate
trees before being accepted or re-baselined, and the four PGEN twins were
re-checked byte-identical to their PLINK runs.

## 5. Dead surfaces, and where the line was drawn

Both fixes killed code. The line taken: **delete what has no callers AND no API
entanglement; note what has entanglement.**
- `CarrierStore::cleanDense` -- DELETED. Four lines, this was its only caller,
  no signature elsewhere depended on it.
- `CarrierStore::urFill` and the dense `clear` path -- KEPT, noted at the three
  places a reader would tidy them away, the way LEDGER #114's arms are. A folded
  pseudo-marker was their only non-zero writer, so they are dead, but removing
  them touches the store's constructor, `scatterUr`, `gatherUr` and the
  `ScratchGuard` predicate. That is a removal candidate for the next
  dead-surface audit, not a passing change. CLAUDE.md's "a flipped pseudo-marker
  keeps the dense routines" is corrected.
- The dense case/control arm in the UR pass -- DELETED as unreachable, and it
  was NOT "value-identical" as this handoff and the commit first said. Half of
  it was: `mean(2 - G_case)/2` then `1 - that` IS `sum(G_case)/(2 n_case)`. The
  other half computed `N_case - het - hom` from the FUNCTION-SCOPE
  `N_case`/`N_ctrl` of `mainRegionInCPP`, declared at `:2463` and **never
  assigned** -- the identically-named locals in the per-marker prepare are a
  different, shadowing declaration -- so under `--is_output_moreDetails` a
  folded pseudo-marker's hom counts came from an uninitialised read. The
  deletion removes it, and the two declarations went with it (the build
  confirms nothing else read them). Noted, not fixed: the surviving arm's
  hom/het predicates are dosage-scaled, so a custom-weighted collapse's `w*g`
  above 2 counts as neither hom nor het.

## 6. Opened by the same question: LEDGER #131 (OPEN, a decision)

The owner's question had a second layer that #130 does not touch. **The MAF
filter and the beta weight still fold, while the annotation is ALT-based.** A
site whose ALT allele is COMMON (AF 0.778) is admitted to a
`--maxMAF_in_groupTest=0.5` cell as MAF 0.222, handed `dbeta(0.222, 1, 25)` =
0.0606 -- the weight of a 22 %-frequency variant -- and printed with `MAC` 444
rather than 1556. So after #130 the collapse and the burden describe one
vector -- still minor-coded per constituent, see section 8.1 -- while the
weight and the admission use the folded frequency. That is upstream's behaviour
and the published SKAT / SAIGE-GENE+ convention, so it is a modelling question,
not an arithmetic defect; the row carries three candidate answers and
deliberately has no characterization test yet, because what it would pin is the
choice. Reachable on bundled data only via P62/P63's raised collapse cutoff or
`--AlleleOrder=ref-first`.

## 7. What is now STALE because of this work

**The FlexRV sweep's collapse-10 columns.** Every collapse-10 number in
`results-2026-09-12.txt` and in `flexrv-study-decisions-2026-09-12.md` was
measured with #129 in, and the standing "read the collapse-0 column" caveat is
now MOOT rather than satisfied -- only a re-run says what the collapse-10
numbers are. Cheapest check that the fix reached the sweep and not just the
bundled fixture: the per-trait counts of sets present at collapse 0 and absent
at 10 (20 / 45 / 44 / 21, `harness/sim/flexrv_study/collapse_presence.py`)
should all be ZERO.

## 7b. P65, and why rolling back to minor-allele coding was rejected

The owner asked, correctly, whether the corrected understanding in section 8.1
argues for rolling back to option A (fold the pseudo-marker like every other
marker). Three measurements say no, and none of them depends on the annotation
argument that drove the original decision:

1. **In a mixed cell, folding sets a quarter to a half of the cell against the
   rest.** Rare markers enter the burden as `+w g`; a folded pseudo-marker
   enters as `w(2 - g) = const - w g`. So a cell holding both would measure
   `rare - ultra-rare + const`, a signed contrast rather than a burden. That
   configuration is REACHABLE and is now golden P65: at
   `--MACCutoff_to_CollapseUltraRare=700`, FIVE cells of the bundled fixture
   hold both halves, and the pseudo-marker carries 28.3 / 28.3 / 28.7 / 28.7 /
   52.5 % of the cell's MAC. Found by sweeping the cutoff between the default 10
   and P62's 2000 (50/100/200 reach no cell, 400 reaches two, 1000 leaves
   nothing rare). P62/P63 cannot discriminate here -- they collapse everything,
   so a fold there only flips a sign.
2. **For a weighted set the fold is not even computable.** The predicate was
   `altFreq > 0.5` on `Sum max(w g) / 2N` -- a weighted sum compared against a
   frequency -- and the operation is `2 - G`. On the committed multiweight
   golden, `GENE1 lof 0.001 CADD` has weighted AC 38.4 over an unweighted MAC
   of 5 with CADD weights 0.5-10, so entries run ~7.7 to ~20 and `2 - G` gives
   -5.7 to -18: not an allele complement. Scaling those weights by 30 would
   make the same genotypes cross 0.5 and be folded. That is #129's category
   error at the same site, which is the sense in which the two rows are one
   mistake.
3. **Folding a collapse inverts its meaning** -- `max(w G)` reads "carries at
   least one", its fold reads "carries none". Independent of 1 and 2.

The one variant left on the table, if per-marker consistency is wanted where
the fold IS well-defined: fold beta / No_Weights sets only (their collapse is
unweighted `max(g)` in {0,1,2}, so both the predicate and `2 - G` mean
something). Its cost is that two weight sets over the same markers would then
disagree in sign, and it brings back the `sumFill` machinery for that half.

## 8. The adversarial review, which found six real things

The diff was put through a six-lens adversarial review before pushing (the
same device that found three defects in the 2026-09-12 lever diff, one of them
live). Every finding below was verified against the code by hand afterwards --
two were re-checked because an agent's framing was wrong even though its
finding was right. What it caught, all now corrected above:

1. **#130's stated rationale was overclaimed** (two lenses found this
   independently, and it is the most important one). "The burden counts the ALT
   allele" is false: the constituents are folded before the collapse. Corrected
   on the row, in the code comment, in the flip test's header and in section 1.
   The consequence is that #131 is WIDER than it was written -- the genotype
   fold belongs to it, and that is its larger half.
2. **The deleted case/control arm read uninitialised memory**, so
   "value-identical" was half wrong. Verified: `N_case`/`N_ctrl` at
   `src/Main.cpp:2463` are declared and never assigned, and the build accepted
   their removal, which is the proof nothing else read them.
3. **PGEN is not a hard-call reader.** `PgrGet1D` decodes dosages and there are
   committed `pgen_marker_dosage*` goldens. The #129 justification survives at
   the defaults but the sentence was false.
4. **Two claims in the P62-P64 golden block were falsified by my own
   measurement**: the UR-row denominator is 359, not 269, and "no group golden
   has a custom weight small enough" was contradicted by the 17 pre-existing
   goldens #129 moved (`group_new_snpid_multiWeights.txt` carries REVEL weights
   of 0.01-0.07). Both corrected in `make_goldens.sh`.
5. **A computed p-value was pinned as an exact STRING** in the inverted #129
   test (`"2.522646E-01"`). That is the byte-vs-numeric trap CLAUDE.md says has
   fired twice; it passed locally and would have been a platform bet on CI. Now
   numeric with a tolerance.
6. **The inverted #129 test was VACUOUS about the collapse** -- every assertion
   in it was also satisfied by a build where the ultra-rare collapse never
   happened, since the two singletons would then be tested as rare markers with
   the same MAC and the same p. It now pins `Number_ultra_rare` == 2 at
   collapse 10 and 0 at collapse 0. It also spent 7 CLI runs on 4
   configurations while its header claimed four; a per-process cache makes the
   header true.

It also produced **LEDGER #132 and #133**, both OPEN, both pre-existing, and
both located by controls rather than argued:
- **#132**: ACAT-V's combination weights are `c(w, rep(1, numofUR))`, so the
  pseudo-marker's fixed 1 does not scale with the user's weights and CCT
  normalises by `sum(w)` -- doubling every weight in a group file, which is
  exactly a no-op for `Pvalue`, `Pvalue_Burden` and `Pvalue_SKAT` (byte-
  identical, measured) and exactly halves BETA/SE, shifts `Pvalue_ACATV` by up
  to 10.5 %. The control that localises it: the same doubling with the collapse
  OFF leaves every column identical. Found by inverting the #129
  characterization in `test-step2-multiweight-behaviour.R`, whose comment had
  attributed the old inexactness to the wrong mechanism.
- **#133**: a cell's `MAC` is folded while `MAC_case`/`MAC_control` are not, so
  `MAC != MAC_case + MAC_control` on any cell above AF 0.5 -- 444 against
  151 + 1405 in P62/P63. `test-invariants-goldens.R` asserts that identity, but
  only on a golden where no cell exceeds 0.5, so it had never been tested where
  it fails. Now pinned in the flip test.

Three findings were REFUTED on verification and are recorded here so they are
not re-opened: that the fix breaks the multi-model or threaded region paths
(the UR pass's per-model loop and per-cell pool are untouched), that
`urFill` being dead is itself a hazard (it is inert, and noted), and that the
new goldens leave an output ungated (`compare_goldens.sh` globs, and would flag
an EXTRA file).

## Open, and needing the owner

- **#131** -- the rarity half of the ALT-vs-minor question, above.
- **#116 + #128** -- unchanged; #128's fix is best decided with #116, whose
  mechanism is still undiagnosed, so that package starts with a diagnosis.
- **The FlexRV sweep re-run** at the fixed HEAD, which is also what any power
  claim now depends on.
- **The binary FlexRV arm at 50k** -- still the only untested scale arm.
- **The pooled row's size AT a discovery threshold**, if anyone wants it pinned
  rather than argued (`prediction-cct-dependence-2026-09-12.md`).
- **`urFill`'s removal**, at the next dead-surface audit.
- **#132 and #133**, both opened by this work (section 8).

## 9. The second half of the session: what was AUDITED, not changed

After the two fixes, the same session audited surfaces rather than changing
them. No source was touched; everything below is a ledger row plus a document,
and all of it is pushed. A later session should read the documents rather than
these summaries.

- **`docs/cli-options-and-group-file.md`** -- all 139 CLI options (step 0: 16,
  step 1: 58, step 2: 65) traced from `make_option()` to the code that consumes
  them, with a live / refused / inert / partly-live status and this branch's own
  additions marked. The sections the help strings cannot give you: 2 refused
  outright, 3 genuinely inert, the trait-gated set (on a quantitative trait
  `--is_output_moreDetails` emits NO columns), the path-gated set, and the ones
  that are live but cannot move a result byte. Appendix: the group-file format
  without weights, with K weight sets, and with FlexRV's score line, read off
  `checkGroupFile`, plus what the annotation mask excludes and why LoF scores of
  1 are a convention rather than an enforcement.
- **`docs/state/open-allele-and-weight-issues-2026-09-12.md`** -- #131 to #135
  on one page with "cost of ignoring it" and "moves goldens if fixed" per row,
  a suggested order, and the ONE root cause they share: SAIGE's single-variant
  idioms (the fold, MAF, `2 - G`, the dosage filter) reused on quantities that
  are not one variant's dosage vector -- which is also what #129 and #130 were.
- **`docs/flag-surprises.md`** -- twelve flag combinations that silently do not
  do what a user expects, from a ten-class adversarial hunt (40 candidates, 14
  verified). Three were hand-checked and have rows: **#139** (chrX doubling
  plus any missing call gives two different wrong answers by reader -- PLINK
  drops the marker, PGEN writes p = 1 beside `nan`), **#140** (a region restart
  can LOSE regions, because #47's byte-count truncation is marker-path only and
  the region writers are never flushed) and **#141** (a duplicated `.fam` IID
  is accepted by all three steps and step 2 then reads the wrong genotype row;
  PGEN's psam reader refuses it). **The other nine are LEADS, workflow-verified
  only** -- confirm before acting.
- Rows **#136** (six false help strings, each re-verified by hand) and **#137**
  (the manifest's promised refusal of `--varianceRatioFile` is missing) and
  **#138** (a group-file `weight` line contaminates `p_FlexRV`, the one
  contamination route of four that is not refused).

**Where to start next, and the reason is sequencing rather than severity:**
**#134**. It is the only one of the five that is wrong in a DEFAULT
configuration, and it moves the same FlexRV numbers that #129's fix already
made stale -- so fixing it before the sweep is re-run means running the sweep
once instead of twice. #139/#140/#141 are the ones that produce wrong results
rather than merely ignoring the user, and all three are cheap.

