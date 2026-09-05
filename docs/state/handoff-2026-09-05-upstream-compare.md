# Handoff: strip vs upstream at N = 50,000 (2026-09-05, after the twelve-decisions handoff)

One session, no R/C++ change. The owner asked for a whole-pipeline speed comparison
against upstream `main` (365deeff) and a -log10 p before/after plot. Everything is in
`harness/sim/upstream_compare/` (`RESULTS.md` first, then `README.md` for how the
arms were built and run), with a section appended to `harness/bench/NOTES.md` and two
items in `docs/state/next.md`.

## What landed

- `harness/sim/upstream_compare/`: `run_arm.sh` (one timed run of one build, records
  the library loaded and the `SAIGE.so` md5), `bench.sh` (paired, alternating, refuses
  a throttled or loaded machine), `analyse.R`, `compare_pvalues.R` (marker-set diff
  first, then per-row tags that RECONCILE rather than band), `plot_pvalues.R`,
  `compare_all.sh`; `results.tsv`, `cmp_summary.txt`, `plots/`, `cmp/` (gene-level
  per-row TSVs; the 70k-row single-variant ones stay under `harness/work`), `logs/`.
- `harness/sim/04_run_saige.sh` / `00_env.sh`: env overrides (`RS`, `SAIGE_EXTDATA`,
  `MODELPREFIX`, `S1OUTPREFIX`, `S2OUT`, `GROUP_NT`) so one script drives either build.
  Defaults unchanged.
- LEDGER #97 (OPEN): the #65 rescale applied to the scale-free `sparse` variance-ratio
  rows. Found by the comparison; fix and test named in the row and in `next.md`.
- Upstream worktree `/Users/dpalmer/Repositories/SAIGE-main` with its own pixi env
  (pinned to strip's conda builds) and `Rlib`. Not in git; ~2 GB; rebuild recipe in
  the README (four traps, each of which cost a build).

## The numbers, one line each

Step 1 strip 1.34-1.37x; single-variant scan 1.28x on quant and **0.94x (n.s.) on
binary**; gene tests 1.33-1.37x; strip at 4 threads 2.35x / 2.72x wall on the scan,
byte-identical. Same-model step 2 is bit-identical to upstream on every gene cell and
every AF <= 0.5 marker; the departures are upstream's #28 (2-7x variance inflation on
AF > 0.5, hand-verified), #29, #50. Own-model departures are #53 and #97.

## What is NOT done

- **#97 is not fixed.** It is a one-line change plus a test; it moves no committed
  golden. Until then every quantitative model fitted at HEAD with
  `--useSparseGRMforVarRatio=TRUE` is conservative on the region path.
- **The binary single-variant scan is not faster than upstream and nobody knows why.**
  Profile at N = 50,000 (`harness/profile.sh` against `sim50k`, `y_binary`).
- The quantitative `null` variance-ratio rows differ 0.2-0.6% between arms; attributed
  to #53 plus the absolute PCG stop (#6) acting in standardised units after #65 --
  consistent, not isolated.
- Timing was 3 reps per arm. Enough for the CIs quoted, not for anything under ~5%.

## Traps met

- Upstream's `--rangestoIncludeFile` did not select a window on this fixture ("No
  markers on chrom 20 are found"); `--idstoIncludeFile` did.
- `compare_pvalues.R` first left stale `.only_*.tsv` files from a partial run, which
  drew a phantom rug on one plot. It now unlinks them first.
- `run_arm.sh`'s step-1 output prefix is per build, not per rep, so the timing reps
  overwrote the own-model files (bit-identical for binary; 1e-7 VR noise for quant,
  LEDGER #56). The own-model step-2 outputs were produced before that.
- The `sim50k` step-1 models on disk are the 09-01 strip fit: identical to HEAD for
  binary, pre-#65 (and hence pre-#97) for quantitative.
