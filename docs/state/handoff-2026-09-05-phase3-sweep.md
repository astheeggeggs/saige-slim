# Handoff: the Phase 3 sweep (2026-09-05, after the workflow handoff)

One session, working the open ledger rows in the order judged most useful:
silent-wrong-answer rows that reproduce cheaply first, then the parser and
guard rows, then the latent guard-audit rows, leaving the design decisions.
Every fix is one commit per row with its test red at unpatched HEAD (a clean
rebuild into the same library, stashing only that row's files) and green
after, and the push went through `harness/prepush.sh`.

## What landed

`python3 harness/tools/ledger.py open` is the authoritative count. Fixed
today, in commit order: #76, #77, #73 (twice), #70, #68, #71, #87, #45, #69,
#49, #43, #85, #88, #89, #83, #25+#81, #14, #78, #72, #23, #22, #24, #21, #35,
#36, #37, #38, #39, #46, #18 (as documentation), #64, #67, #90. Plus #57's
open half and the settled #82 mechanism from earlier in the day. One new row,
#90, was found and then fixed.

## What the pre-push gate caught, twice

- **#73, first draft**: refusing any conditioning marker the scan's own QC
  would drop. Golden P4 conditions on rs13 (MAC 2) under `--minMAC=20` --
  conditioning on a marker rarer than the scanned ones is legitimate. The
  guard is now: monomorphic stops, scan thresholds warn.
- **#73, again**: `test-characterization-ledger16.R` conditioned on rs80,
  itself monomorphic; the construction only needed a mid-stream cursor, so it
  is on rs79 now.

Both drafts had passed every test file run by hand. The gate is the thing.

## Traps met this session, for the next one

- **Never edit `tests/` or `extdata/` while the gate runs.** They are read
  from the tree by the running tests; an appended test went red against the
  old install and cost a gate cycle. `R/` and `src/` are safe (the gate tests
  the INSTALLED package).
- **A HEAD drill that stashes a file stashes every uncommitted row in it.**
  Group rows per file, or drill before starting the next row in that file.
- **`pkill -f <test name>` matches the drill script too** if the script's
  command line contains it. Kill the R child by another pattern.
- **Do not drill #78's refusal arm at HEAD**: without the guard it asks for
  a 1,000-dummy fit. `SAIGE_DRILL_HEAD=1` skips it; the warning arm is the
  drill.
- **The zsh Bash tool does not word-split `$var`.** A multi-token argument
  string must go through `bash -c` or a script.
- **Adjacent ledger rows share one `-U0` diff hunk**; stage `LEDGER.md` per
  row by line-mixing, not by patch (see `commit_rows.sh` in the session
  scratchpad, reproduced in spirit here: HEAD's file with one row's line
  taken from the final file).

## What remains open, and why: decisions, not fixes

- **#6** PCG tolerance: deliberately deferred since Phase 0.
- **#19** `weights.beta` reaches C++ but not the R-side ACAT-V: a design
  asymmetry inherited from upstream; fixing it moves every multi-weight golden.
- **#20** `fast_logistf_fit`'s unused `col_fit`/`maxhs`/`lconv` and the
  `maxit = 0` loop: honouring the parameters is implementing logistf's
  step-halving; dropping them changes an exported signature.
- **#26** Firth drifts to +-750 on null markers: cap, report NA, or leave.
- **#51** `--is_fastTest` is inert for group tests: making the region path
  read it changes every fastTest region golden.
- **#65** three AI-REML guards absolute where documented relative: changing
  them can move every step-1 golden.
- **#79** a cost, not an answer; belongs to the perf track.
- **#80** `Is.SPA` true on ER markers: renaming the column's meaning moves a
  column in every marker golden and the `SE == 1/sqrt(var)` invariant.
- **#82** CLI and API default `is_noadjCov` oppositely: which one changes?
- **#84** the conditional SPA saddlepoints against the marginal CGF: a
  statistical change needing its own derivation and a fixture at real
  conditioning (the high-LD one exists).
- **#86** SKAT-O saturating to exactly 1 lives in the SKAT package; the
  Cauchy fallback here is a policy choice.

## Gate numbers

Recorded by the push that follows this handoff; see the CI run on
`astheeggeggs/saige-slim` for the same commit.
