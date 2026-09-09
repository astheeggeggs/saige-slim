# Handoff, 2026-09-09: FlexRV F0, and the decisions that came with the paper

Docs only, four commits, no code and no gate movement. Previous handoff:
`docs/state/handoff-2026-09-09-skato-identical-cells.md` (a different session
was working in this checkout at the same time -- see "The checkout was shared"
below, it matters).

`722c3610` F0, `ab047839` the F3/F4/F5 decisions, `93f41251` next.md,
this commit the handoff plus one design-doc correction.

## What changed

**The owner supplied the FlexRV PDF, which was the blocker on F0 and F2.**
`docs/design/flexrv-design.md` now pins the 16x12 transform table VERBATIM from
the Methods, derives each closed form, and records what it is safe to assume.
It SUPERSEDES the table in `flexrv-plan.md`, which came from a lossy fetch.
The PDF is not committed (5.7 MB); its sha256 is in the design doc.

`pdftotext -layout` extracts the Methods cleanly and Extended Data Fig. 1 is on
PDF page 13, readable as a rendered page. Neither needed a web fetch.

## The one resolution that changes the program

The Methods say the MAF beta transform "followed the definition above, but we
mapped input scores s to w as f(MAF / 0.001)". Substituted literally into their
own score-beta formula that gives `(MAF/0.001)^(b-1)`, which INCREASES with
MAF. It is contradicted by their next sentence and by Extended Data Fig. 1b,
whose axis is allele count 1-1000 and whose two beta panels start at weight 1
and fall to 0.

Pinned as `(1 - MAF/0.001)^(b-1)`, zero at and above the boundary. At MAC = 100
the pinned form gives 0.900 (b=2) and 0.656 (b=5); the literal reading gives
0.100 and 1e-4, and the b=5 panel reads about 0.65 there, so the literal
reading is excluded rather than merely disfavoured. Their figure axis coincides
with `MAF/0.001` only because 2N is about 1e6 at their N; we key on MAF, which
is N-independent.

Four more ambiguities are resolved in the doc (A2 the LoF mask is
annotation-based though the figure draws it as `s == 1`; A3 their thresholds are
MAF values chosen to match MACs at their fixed N, ours read in-sample MAC; A4
their typeset CCT drops both its parentheses and its 1/k, so we use the repo's
`CCT_cpp` unchanged; A5 which AC their pseudo-variants' thresholds read).

Every closed form was computed two ways before being written down: score beta
`dbeta(1-s,1,b)/dbeta(0,1,b)` vs `s^(b-1)` agree to 1.1e-16, MAF beta to
1.7e-16, cubic root hits `w(0)=0`, `w(1)=1`, `w(t)` = 0.3865 / 0.5 / 0.6753 at
t = 0.2 / 0.5 / 0.9. Those are the values F2 pins. The `b=2` score panel being a
straight line is a free check on the derivation, since `s^1` is the identity.

## The finding that shapes F3

**SAIGE's ultra-rare collapse collides with four of the twelve MAF
transforms.** `MAC<=1,3,5,10` sit at or below `--MACCutoff_to_CollapseUltraRare`
(default 10), so every variant they select has already been collapsed: every
real marker gets weight 0 and the entire test is carried by the pseudo-marker.
The masking ORDER is already right -- the collapse is weighted and per set
(`genoURMat` is `N x q_anno_maf_weight`, `src/Main.cpp` ~:2598), so each cutoff
gets its own pseudo-marker -- but the operator is `max`, an "any qualifying
carrier" indicator, not the additive weighted count the paper's corresponding
column reports. `w0 = 1` does NOT bias those four sets: with no real markers
left the pseudo-marker IS the test and a burden p-value is scale-invariant.

We do NOT port the paper's AC x decile pseudo-variant binning. It is a STAAR
speed workaround whose own validation is that it approximates the unbinned
additive burden, and the unbinned test is exactly available with
`--MACCutoff_to_CollapseUltraRare=0` (one condition,
`if(!(MAC > g_cfg.region_minMAC_cutoff))`, `src/Main.cpp` ~:2345; `checkArgs.R`
rejects only negatives). Recommendation by trait: cutoff 0 for quantitative
(the paper's own 62 traits), keep the collapse for binary and drop the four
transforms from the spec string. Upstream's inactive `sum_geno` stub is the
cheap middle option if F5 says the operator matters.

## Decisions recorded in the plan (owner's, taken in session)

- **One MAF cell and one annotation group per FlexRV run**
  (`--flexRV_maxMAF`, default 0.001). Both dimensions must be pinned or the
  region's pooled Cauchy row is a super-omnibus over anno x maxMAF x set and is
  not `p_FlexRV`. Internally `maxMAFVec = c(M)`, `q_maf = 1` -- no new
  machinery, and it retires the deferred per-cell-Cauchy-row question.
- **M != 0.001 warns and renames the sets.** `mafbeta*` rescales to
  `(1 - MAF/M)^(b-1)` because the paper defines it relative to the inclusion
  boundary; the MAC ladder does NOT rescale, because those thresholds are
  absolute and interpretable. M goes in the set name so the OUTPUT is
  self-describing.
- **`--flexRV_lofAnno` is an exactly-matched list** (`,` `:` `;`
  interchangeable). SAIGE has no annotation vocabulary, so `pLoF`, `plof`,
  `ptv` are all equally valid; a label matching nothing is a hard error naming
  what IS present, or 12 of 192 sets vanish from the pool silently.
- **Warn on a constant or {0,1} score column**; write every set's row by
  default (`--flexRV_output_sets=FALSE` to suppress), measured at 70 B/row.

## The checkout was shared, and it changed something under us

Five commits from another session landed between `ab047839` and `93f41251`.
Two matter here:

- **`1b5bf4c1` implements the compute-once-and-copy memo this session had just
  argued for, on the SKAT-O path** (`skatoDedupe`/`skatoMemoKernel`,
  `R/SkatoParallel.R`), keyed on `(Score, Phi, r.corr, regionTestType)`. Two
  FlexRV sets with identical weights should hit it, which would remove the
  192-Davies-call problem for free -- EXPECTED, not verified. `skatoDedupeSummary()`
  prints the counts; F5 reads them rather than assuming. The C++ burden and
  ultra-rare side is still undeduped and is where the `N`-length work is.
- **`056b3d93` adds `harness/tools/docdrift.py`**, a blocklist of falsified
  statements that runs before EVERY tier and scans `CLAUDE.md`,
  `docs/state/next.md`, `R/`, `man/`, `extdata/` (not `docs/design/`, not dated
  handoffs). Checked after the next.md edit: `selftest` 6 rules OK, scan 0 hits
  over 34 files.

`git log --oneline` interleaves the two sessions. Anyone timing anything in
this checkout should re-read the shared-checkout trap first.

## NEXT UP

**F2, and it is fully specified by the design doc.** `src/FlexRV.cpp` with the
16 score and 12 MAF kernels plus a `flexrvWeightsInCPP()` export, and
`tests/testthat/test-unit-flexrv.R` pinning each transform to the closed form on
a grid: `s = 0`, `s = 1`, `s = t`, the `MAF = 0.001` boundary from both sides,
`MAC = t` and `MAC = t+1`, and the product ordering `k = c*12 + d`
(score-major, MAF-minor, normative and identical in R and C++). Fast tier.
Remember `OBJECTS` in `src/Makevars`, `compileAttributes`, and that a header
change means deleting `src/*.o`.

**Three decisions are the owner's and none blocks F2**: whether the
quantitative/binary collapse split is advice or an enforced default; whether the
constant-score warning earns its surface; whether per-set rows are written by
default (recommended yes).

**F4 carries a constraint that costs a golden-generation cycle if rediscovered
late**: `genotype_100markers` is N = 10,000, so at the goldens' default max MAFs
most of the MAF axis is degenerate -- the 0.0001 cell has no real markers at
all, the 0.001 cell makes six of nine MAC thresholds duplicates of `unw`. Choose
`--flexRV_maxMAF` and the collapse cutoff deliberately, and assert how many of
the 192 weight columns are actually distinct.

**F5's first measurement comes before any p-value**: the fraction of individuals
carrying >= 2 qualifying ultra-rare variants per gene, straight from the
genotypes. It bounds a priori how far the `max` collapse can diverge from an
additive one, and may close the `sum_geno` question without a test being run.
