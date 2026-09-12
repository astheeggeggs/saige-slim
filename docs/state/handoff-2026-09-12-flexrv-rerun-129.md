# Handoff 2026-09-12 (morning): the FlexRV sweep re-run at HEAD, and #129

Read this first, then `docs/state/next.md`, then
`python3 harness/tools/ledger.py open`.
Previous: `handoff-2026-09-12-126-null-calibration.md` (the #126 null
calibration, #128), then `handoff-2026-09-11-evening-117-115-114.md`.

## Why the re-run

The owner asked whether the FlexRV simulation study was done. It was not: the
144-run sweep of 2026-09-10 had run at `5dc679c3`, before #118 / #126 / #127 /
#117 landed, its raw outputs were deleted on 2026-09-11, and the study was
blocked on three owner decisions. The owner said to re-run. Same grid, same
frozen step-1 models, HEAD `38eee5f5`, prediction written first
(`docs/state/prediction-flexrv-rerun-2026-09-12.md`, outcome appended after).

## What it said

- **Pipeline control held to the digit**: all 40 `skato`/`grid` summary lines
  (beta-weight SKAT-O, untouched by the four rows) reproduce 2026-09-10
  exactly; `unw`/`beta` moved only through the two genes the gate now admits.
- **Every FlexRV detection count fell by 1-2 and none rose** (flex192 13 ->
  11 on `y_binary`), inside the pre-registered bounds. #118's per-set MAC
  routes small-support sets to the sparse-GRM estimator and their own
  variance-ratio band, which is more conservative here. **At HEAD FlexRV-192
  equals the standard grid (11 = 11) and beats the single burden by 3**; the
  2026-09-10 margin over the grid was a pre-#118 number.
- **The three decisions are unchanged** (predicted: they are structural).
  The memo for the owner is `docs/state/flexrv-study-decisions-2026-09-12.md`
  (drafted by three independent perspectives, merged, each option
  adversarially refereed; every claim not traceable to a file is marked).
- **LEDGER #129 (new, OPEN).** `src/Main.cpp:3127` passes the WEIGHTED
  ultra-rare pseudo-marker through `imputeGenoAndFlip` with
  `--dosage_zerod_cutoff` (0.2) and `--dosage_zerod_MAC_cutoff` (10);
  `src/UTIL.cpp:52-55` cleans every entry <= 0.2 when the weighted MAC is
  <= 10. Consequences, all measured: transform sets whose support is all
  ultra-rare with weights <= 0.2 vanish at collapse 10 (20 / 45 / 44 / 21 of
  ~40k rows per trait across the four score arms, never the reverse); and a
  partially-zeroed set is tested on the surviving carriers while its MAC
  columns count everyone -- `GENE0181 sbeta10`: 1.7e-3 at collapse 0, 2.9e-9
  at collapse 10, 1.710212E-03 again at collapse 10 with the zeroing off.
  Bundled reproduction and `test-step2-region-ur-zerod.R` (weights 0.01
  vanish, 0.5 survive, the burden p equal in every tested case). Not fixed;
  two local fixes on the row. **Read the collapse-0 column of the study
  until it is** -- it is the only #129-free one.

## Traps this session added

* **A prediction that "support does not depend on the collapse" was wrong**,
  and the first collapse-0 file showed it within a minute of the sweep
  starting (20 rows more than collapse 10). Compare the first output of a
  sweep against the prediction before the sweep finishes.
* **`--dosage_zerod_cutoff=0` switches the zeroing off entirely**
  (`if(t_dosage_zerod_cutoff > 0)`), which is the one-flag way to test whether
  a collapse-dependent p-value is #129.
* **`grep -P` is ugrep's in the tool shell and does not exist in a bash
  script on this Mac**; use awk in scripts.
* **The 2026-09-10 "reference" detection counts are like-for-like only with
  the same cell grid**: `unw` is one cell, `grid` six; last night's
  nine-cell calibration arms were not comparable to them (corrected in
  `prediction-126-null-calibration.md`).

## Where things are

`harness/sim/flexrv_study/`: `results-2026-09-12.txt`, `compare_results.py`
(line-by-line diff of two analyses, drilled), `collapse_presence.py` (#129's
set-presence count), `zerod_check.sh` (the GENE0181 proof). Raw run outputs
under `harness/work/flexrv_study/results/` (gitignored, 176 files, ~600 MB;
`study.sh` is idempotent). LEDGER: #116, #128, #129 open.
