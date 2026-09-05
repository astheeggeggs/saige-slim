# Handoff: workflow hardening (2026-09-05, after the Phase 3 start handoff)

Three commits, no R, C++ or golden change. Prompted by an assessment of the
workflow: the evidence discipline is strong, the prose that carries it had
started to drift (a hand-written open count wrong by five; two files naming
different "newest" handoffs; a ledger row invisible to its own tool; a
compare-mode trap documented but unenforced).

## What landed

| commit | what |
|---|---|
| `fe6a665a` | `ledger.py lint`, run before every tier and in CI; 15 ledger rows repaired |
| (next) | `harness/prepush.sh` + `.githooks/pre-push`: lint + standard tier in BOTH compare modes before any push to `slim` |
| (next) | `CLAUDE.md` 896 -> 376 lines; one entry-point rule; `history-claude-md.md` holds the old file whole |

## The finding

**The true open-row count was 44 of 89, not 39.** The previous handoff found
one row (#66) whose prose pipes broke the markdown table and hand-corrected it,
inferring the true count was 40. The lint finds FIFTEEN such rows, and five of
them are OPEN rows `ledger.py open` had never listed: **#26, #45, #64, #65,
#77**. The old parser kept the first five cells of a broken row, so STATUS was a
prose fragment and the row was neither open nor fixed. Two rows (#27, #41) also
carried an "upstream report" link in a sixth cell that was never displayed.

Repair was content-neutral after unescaping (checked row by row with difflib):
prose pipes are now `\|`, the GFM table escape. Rule recorded in CLAUDE.md.

## Gate numbers

`bash harness/prepush.sh` at `fe6a665a`: byte 2263 pass / 0 fail, numeric 2438
/ 0 -- each +3 on the previous handoff's last row, which is what the #63 commit
added. 8 min wall at load average 5; expect ~4 min on a quiet machine.

Detection drills, each red then green: a planted bare pipe, a renumbered row, a
STATUS of "maybe fixed" (the inherited test accepted "FIXED" anywhere in the
first 16 characters), a planted "N open of M" in `next.md`, a source file newer
than the installed `SAIGE.so` (the gate refuses), the hook on `origin` (not
gated) and on `slim` with `SAIGE_SKIP_PREPUSH=1` (skips, loudly).

## What is NOT done

- The hook is armed per clone by `git config core.hooksPath .githooks`; that
  is not versioned and every fresh clone must run it (Build section).
- `CLAUDE.md` is 376 lines against a 200-line target. It is rules and traps
  now, not narrative; further cuts are judgement calls left to the owner.
- `docs/state/next.md` still carries per-item narrative that overlaps the
  handoffs; it was not restructured, only re-pointed.
- Phase 3 fixes continue from `python3 harness/tools/ledger.py open`.
