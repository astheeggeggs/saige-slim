# Handoff 2026-09-10 (night): the FlexRV sweep ran, LEDGER #124 is measured and held, the comparator floor

Read this first, then `docs/state/next.md`, then `python3 harness/tools/ledger.py open`.
Previous: `handoff-2026-09-10-fitter-guard.md` (the fitter rows #121-#123, still
correct), then `handoff-2026-09-09-flexrv-e2e.md`.

## Done tonight

1. **The FlexRV study sweep, S2-S4 (144 runs) and its analysis.** Results
   `harness/sim/flexrv_study/results-2026-09-10.txt` (the `analyse.R` output),
   the write-up in `harness/bench/NOTES.md` ("FlexRV simulation study S2-S4"),
   the design doc's status line points at both. Read the controls before the
   tables: null calibration is fine (8 of 288 lines outside their 95%
   interval, chance rate 5%, both directions, never at both alphas of one
   cell); **the positive control FAILED on `y_quant`** (Spearman(NCP, -log10 p)
   -0.16..0.25 for every method, 0-3 of 63 causal genes detected), so nothing
   about power is read on the quantitative trait at this N; on `y_binary`
   rho 0.26-0.38 and 7-13 of 63 detected. **Two control expectations in the
   design are mis-specified**: the "uninformative" arm randomises only the
   SIM score while the `lof` set (annotation) and the 12 MAF transforms stay
   informative, so flex192 on an uninformative score beats the unweighted
   burden (12 vs 8) by construction; and the `low` class carries causal
   variants in this simulation (the synonymous-style mask detects 7/63 causal
   genes on `y_binary`), so "expected null on the signal traits" does not
   hold as built. Both are [DECISION]s for the owner before the power question
   can be answered. Three driver bugs fixed in `run.sh`/`study.sh` (details in
   NOTES); the printed "N rows" of `run.sh` counts the output's lines including
   the header, so a header-only file reports as a success -- checked by hand
   tonight, worth a guard.
2. **LEDGER #124 measured end to end and HELD for the owner's decision** (the
   row is authoritative). Characterization `test-step1-ai-projection.R`
   committed in its TODAY form (green at HEAD, standard tier, two child fits);
   the inverted form and the three-token patch are in
   `harness/work/flexrv_study/ledger124/` (gitignored: re-derive from the row
   if that directory is gone). The fix's gate: goldens move at 1e-4 (theta)
   and <= 2.4e-4 of scale (every field), iteration counts unchanged, two fixed
   trees byte-identical -- but the byte-mode standard tier pins the step-1
   goldens byte-exactly (multipheno, spsolve-cache, the guard test's trace at
   1e-6 against the golden LOG), so landing it regenerates every step-1 golden
   and log and every step-2 golden (step 2 consumes the models). Options in
   the row. The installed package is back at HEAD (unpatched), md5 `db0fcf54`
   was the pre-session install; verify with `md5 -q` on `SAIGE.so` before
   trusting any number.
3. **`harness/compare_step1.R` floors its elementwise ratio at the field's
   scale** for the model's vector/matrix fields (theta and the variance ratios
   keep the floorless per-component ratio). Reason and the five-plant drill in
   the file's comment and the commit; the drill caught my first version
   (flooring theta too let a 3% tau2 move through because tau1 = 1 sets the
   scale). No golden or test result changes at HEAD.

## The environment incident (cost ~40 min, three background jobs)

Another session (saige-6d) had a git worktree whose `.pixi` SYMLINKED to this
checkout's env; every `pixi run` from there re-solved against the worktree's
manifest and re-linked the shared env (19:41:51, 19:54:00), during which
`Rscript` and `libtbb` were absent and every SAIGE process loading died
(seven golden fits, three sweep windows -- all re-run). The worktree is gone.
Separately I stopped an `R CMD INSTALL` mid-compile while the sweep was still
running (the handoff's own warning; I had misread the sweep as finished): R
had already moved the package into `00LOCK-SAIGE/SAIGE`; moving it back
restored it byte-identically. Memory note `pixi-run-relink-hazard`. Rules
that follow: never symlink a worktree's `.pixi`; call
`.pixi/envs/default/lib/R/bin/Rscript` directly in unattended runs (the
sweep driver does; a `pixi run`-free copy of the golden script is trivial);
confirm a background sweep is finished from its task notification, not its
log's last line.

## Handed to me, not gated

saige-6d: `src/Main.cpp:2988` `genoURVec_noweights = genoURMat_noweights.col(jm)`
is built once per weight set inside the r loop and read only under
`q_weight == 1`, so at K = 192 it is 191 dead N-length copies per cell per
region (~16 GB of memcpy over the 210-gene 50k run). One-line fix (a view),
byte-identity expected trivially, NOT gated by anyone.

## Uncommitted state of the tree

Mine, committed tonight: the #124 row + test, the comparator floor, this
handoff and the next.md pointer. Left uncommitted on purpose because they
sit inside other sessions' uncommitted files or the FlexRV e2e session's
untracked work (the owner's to commit, as the previous handoff says):
`harness/sim/flexrv_study/{run,study}.sh` fixes and
`results-2026-09-10.txt`, the NOTES.md section, the design doc's status
line, and rows 118-120. saige-ae holds `test-step2-skato-parallel.R`,
`docdrift.py` (rule `skato-served-count`) and a CLAUDE.md hunk; do not
`git add` those paths.

## Next, in order

1. Owner: #124 (row), the two study [DECISION]s, and whether `y_quant`'s
   planted NCP should be raised for a usable positive control.
2. The estimator test (design in `handoff-2026-09-09-flexrv-e2e.md`) -- the
   in-process kernel-call pattern in `test-step1-ai-projection.R` is the
   reusable piece (intercept `getAIScore*` in the namespace, dense inputs).
3. #114-#118 as before.
