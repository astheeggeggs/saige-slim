# Handoff: the twelve decisions (2026-09-05, after the Phase 3 sweep handoff)

The Phase 3 sweep left eleven open ledger rows plus the second half of #18,
each a decision rather than a fix. This session laid the options out for the
owner (what the code does, what each option moves, a recommendation with its
one reason and its gate), the owner answered each with a word, and the answers
were implemented under the branch rules: one commit per row, the row's test
red at unpatched HEAD (a scratch install of HEAD, `R_LIBS=<scratch>/headlib`)
and green after, goldens that move re-baselined only after two independent
regenerations agreed byte for byte, the standard tier green in byte mode
before the push, and `harness/prepush.sh` before the push itself.

## The decisions and what landed

| row | decision | what changed | goldens |
|---|---|---|---|
| #86 | leave, document | `--r.corr` help + roxygen say the Cauchy row is a Bonferroni when any cell's SKAT-O p is exactly 1 (SAIGE's own edit of STAAR's CCT, checked against both sources) | none |
| #19 | honour | the implicit single beta set uses the first `--weights.beta` pair on the R side too; one statistic two ways vs the explicit two-pair path | none (all goldens at 1,25 or two pairs) |
| #82 | API | `SPAGMMATtest` defaults `is_noadjCov = TRUE` like the CLI; the inverted help rewritten; a roxygen added | none |
| #65 | standardise | `fitNULLGLMM` fits, and estimates variance ratios, on the phenotype divided by the sd of its covariate-adjusted values, then scales the saved model and the ratio file back (`R/PhenoScale.R`) | quant step-1 (+pgen), P4, P16, P32, P36, P37 (+singleAssoc), rounding level |
| #84 | residualise | the conditional saddlepoint uses the residualised genotype's CGF, full SPA (`docs/design/conditional-spa-residualised.md`); wall time unchanged (1.9 s P12 before and after), p.value_c moves only on rows that clear the gate | P12, P17, the single-variant files of P3, P10, P15, P24, P34, and (unpredicted, tier-caught) the BURDEN conditional P25/P26 on one cell |
| #80 | leave, document | `--SPAcutoff`/`--max_MAC_for_ER` help + roxygen: Is.SPA means the SE was back-derived from the p (saddlepoint OR exact test) | none |
| #79 | batch | one SuperLU factorisation of sparse Sigma reused across right-hand sides while (W, tau) are unchanged; 167 solves / 11 factorisations on the binary fixture | none, verified identical() |
| #26 | flag | the Firth cap is reported as not converged; BETA/SE unchanged | none |
| #18b | API | `SPAGMMATtest` defaults `is_single_in_groupTest = FALSE` like the CLI | none |
| #51 | leave, document | the bounded 8.3% speed-up recorded in `next.md` | none |
| #6 | leave, document | the absolute PCG stop recorded as a hazard in `next.md`, with the drill to run first | none |
| #20 | drop | `fast_logistf_fit` loses `col_fit`/`maxhs`/`lconv`; `maxit = 0` returns `init` | none |

`ledger.py` gained a `DECIDED` status (and a `decided` subcommand) for the
four leave-and-document rows, so `ledger.py open` lists only real work.

## Findings worth keeping

- **#65 had to wrap the variance ratios, not just the fit.** With only the
  AI-REML loop standardised, y * 100 was end-to-end consistent (step-2
  p-values to 1.5e-7) but y / 100 gave variance ratios 33% off and step-2
  p-values 57% off: two absolute `1e-4` floors on Sigma's diagonal in
  `src/SAIGE_fitGLMM_fast.cpp` (`gen_sp_Sigma`, `getDiagOfSigma`) bite when
  tau ~ 2e-5. That is the LEDGER #6 shape on the variance-ratio path, and it
  is why the standardisation sits in `fitNULLGLMM`, not in the fitter.
- **The variance ratio scales as 1/s^2 by convention**, not 1: its numerator
  variance is in 1/y^2 units, its denominator in genotype units only. A
  scale-invariance test that assumed a scale-free ratio was wrong and said so
  at once (the x100 arm agreed to 5e-8 once the convention was applied).
- **`spsolve()` and `spsolve_factoriser` are the same four SuperLU calls** at
  the default options (`refine = REF_NONE`), read in `sp_auxlib_meat.hpp`
  before coding #79; the byte-identity of all four sparse step-1 goldens
  confirmed it.
- **SKAT-O's p is bounded below by its smallest component, not above by its
  largest**; the #86 defect is the point mass at exactly 1, not the ordering
  the ledger row had claimed. Corrected in the row.
- **The BURDEN region path saddlepoints its conditional burden score through
  the single-marker block** (`mainRegionInCPP` calls `getMarkerPval` on the
  burden sum with `isCondition`), so #84 moved P25/P26 too. Predicted "every
  region file unchanged"; the tier said otherwise. SKAT-O's conditional burden
  p is the R-side `SPA_pval` and did not move.
- **A stray `// [[Rcpp::export]]` above a comment block silently drops the
  export** (compileAttributes attaches it to whatever follows and warns about
  an "invalid parameter"). It cost one clean build.

## Traps met, for the next session

- The drill library: `git archive HEAD | tar -x -C <scratch>/head-src`,
  symlink `.pixi` and `plink-ng` into it, and install FROM THE REPO DIRECTORY
  with `pixi run R CMD INSTALL -l <scratch>/headlib <scratch>/head-src`. Run
  from elsewhere, `pixi run` does not put the conda compiler on PATH and the
  build fails with `x86_64-apple-darwin13.4.0-clang++: command not found`.
- Two rows that both change `RcppExports.*` cannot be committed separately
  from one regeneration: stash the other row's sources, `compileAttributes`,
  commit, pop, `compileAttributes` again (#20 then #79).
- `R CMD INSTALL .` without cleaning is fine for R-only rows; any header
  change goes through `harness/build.sh`.

## Gate numbers

Standard tier, byte mode: 2536 pass after the R-side rows; 2592 pass after
#65 with its re-baselined goldens. The final `harness/prepush.sh` numbers are
in the push that follows this handoff.
