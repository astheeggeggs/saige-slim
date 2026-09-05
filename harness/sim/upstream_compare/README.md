# strip vs upstream comparison harness

Drives the same `sim50k` commands (`harness/sim/04_run_saige.sh`) through two builds
and compares speed and p-values. Results and their reading: `RESULTS.md`.

## The upstream arm

Worktree `/Users/dpalmer/Repositories/SAIGE-main` (detached at 365deeff = `origin/main`),
its own pixi env, package installed in `SAIGE-main/Rlib`. Upstream needs `savvy`,
`r-metaskat`, `r-rsqlite`, `r-qlcmatrix`, `r-dplyr`, none in strip's `pixi.toml`, so it
cannot share strip's env. Four traps, each of which broke the build once:

1. Unpinned, `pixi install` resolves R 4.5.2 / OpenBLAS 0.3.34 (strip has 4.4.3 /
   0.3.30) -- a timing confounder -- and the newer conda R's `FLIBS` carries
   `-lheapt_w` whose directory does not exist (`ld: library not found for -lheapt_w`).
   Pinning the VERSION is not enough (r-base 4.4.3 build `h07b34e9_7` has the same
   flag); pin the BUILD STRINGS strip's env has (`pixi list` there):
   `r-base = { version = "4.4.3", build = "h46d4a84_2" }`,
   `openblas = { version = "0.3.30", build = "openmp_h30af337_0" }`,
   `gfortran_osx-64 = "13.4.0"`, `libgfortran5 = "14.2.0"`.
2. Pin `r-base`, not the `r` metapackage (`r = "4.4.3.*"` does not solve).
3. Upstream's `platforms` lists osx-arm64, for which that r-base has no candidate, so
   the whole solve fails: set `platforms = ["osx-64"]` in the worktree copy.
4. `LinkingTo` includes `lintools` (not on conda-forge): `install.packages("lintools",
   lib = "Rlib")` first, then `R_LIBS_USER=$PWD/Rlib pixi run R CMD INSTALL -l Rlib .`
   from the worktree root (run elsewhere, `pixi run` does not put the conda clang on PATH).

plink-ng: copy `SAIGE/plink-ng` (v2.0.0-a.6.16) into the worktree, compile
`libplink2_includes.a` as `harness/build.sh` does, keep a copy at the worktree ROOT
(upstream's Makevars links `-L.. -lplink2_includes`) and in `.pixi/envs/default/lib`.

Upstream CLI differences that matter here: PGEN is `--pgenPrefix` and refuses without
`--AlleleOrder=ref-first`; `--rangestoIncludeFile` did not select the window this
harness used ("No markers on chrom 20 are found"), so the smoke test uses
`--idstoIncludeFile`; `--nThreads` is accepted and ignored (LEDGER #30).

## Running

```sh
# one run, timed, appended to harness/work/upstream_compare/results.tsv
bash run_arm.sh <strip|main> <step1|s2single|s2group> <trait> <rep> [nthreads]
#   MODELS=frozen (default): step 2 reads the fixture's committed strip models
#   MODELS=own: step 2 reads the models THIS build's step1 run wrote

# paired timing block: arms adjacent, order alternated per rep; refuses if throttled or load > 2.5
bash bench.sh 3 "step1 s2single s2group" "y_binary y_quant"
Rscript analyse.R                       # per-cell means, paired t, chunk-regression work_s

# p-values: every (models x trait x kind) pair both arms have; TSVs + plots
bash compare_all.sh [rep] [nthreads]
```

`run_arm.sh` points `04_run_saige.sh` at a build through `RS`, `R_LIBS_USER`,
`SAIGE_EXTDATA`, and at prefixes outside the fixture dir through `MODELPREFIX`,
`S1OUTPREFIX`, `S2OUT`, `GROUP_NT` (the overrides added to that script for this).
Each log starts with the library `find.package` resolved and the `SAIGE.so` md5.

Compatibility smoke (both directions, 500-marker window via `--idstoIncludeFile`):
each build reads the other's `.rda` + `.varianceRatio.txt` and returns the same 213 rows.

## Reading the comparator's tags

`compare_pvalues.R` tags every departing row with its cause and reconciles, rather than
bands: `VR` = a chi-square rescale by the two arms' reported `var` predicts strip's p to
1e-3 in log10 (single) or the `SE_Burden²` quotient equals a variance-ratio row quotient
(gene); `VR-mix` = between the row quotients (mixed MAC categories); `AF>0.5` = LEDGER
#28; `Cauchy`/`Cauchy=1` = the per-gene aggregate row (#86 when exactly 1); `tiny` <
1e-5; `small` < 0.05 unreconciled; `UNPREDICTED` = stop and explain.
