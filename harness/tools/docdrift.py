#!/usr/bin/env python3
r"""Fail when a doc repeats a statement that a later commit falsified.

  docdrift.py            check the current-state docs; exit 1 and name each hit
  docdrift.py selftest   the controls: every rule fires on its own canary,
                         and no rule fires on benign text

This is a REGRESSION guard, not a drift detector. It cannot tell that a doc
has gone stale; it can only tell that a sentence known to be false is back.
That is the failure mode actually observed on this branch: on 2026-09-09 five
copies of two false claims were found, all dating from 2026-09-06, and they
were false in the same way because they had been COPIED from each other --
CLAUDE.md into the roxygen, the roxygen into man/, man/ into the CLI help.
One of them had already been corrected once in `next.md` (the step-3 bullet at
line 110 was struck through) while its twin twenty lines away was not.

So each rule below is one sentence that was true when written, the commit that
made it false, and the evidence you can re-run. A rule earns its place only
after the claim has been written down wrongly at least once.

SCANNED is the set of surfaces that describe the CURRENT state. Dated
handoffs, `LEDGER.md` and `docs/state/history-claude-md.md` are deliberately
NOT scanned: a handoff that says what was true on its date is correct, and
rewriting history to please a linter would destroy the branch's record of how
it got here. A line that must quote a falsified claim on purpose can carry the
marker `docdrift-allow` and is skipped.

`harness/run_tests.sh` runs this before every tier, beside `ledger.py lint`.
"""
import re, sys, pathlib

REPO = pathlib.Path(__file__).resolve().parents[2]

# The current-state surfaces. Anything dated or historical stays out.
SCANNED = ["CLAUDE.md", "docs/state/next.md", "R/*.R", "man/*.Rd", "extdata/*.R"]

ALLOW = "docdrift-allow"

# Each rule: a pattern matched against the file with newlines flattened to
# spaces (these claims wrap across lines and one of the five was missed by a
# line-wise grep for exactly that reason), the commit that falsified it, and
# how to re-check. `canary` is the sentence as it was actually written, and is
# the rule's positive control in `selftest`.
RULES = [
    dict(
        id="region-loop-serial",
        pattern=r"marker loop stays on one thread",
        why="the C++ region marker loop has been threaded since 2026-09-06",
        evidence="3eeb554e; R/SAIGE_SPATest_Region.R prints 'Region marker loop will use "
                 "N in-process threads'; tests/testthat/test-step2-region-threads.R",
        canary="For a GROUP test the marker loop stays on one thread, but `nThreads` is "
               "not wasted there: it drives the SKAT-O post-processing",
    ),
    dict(
        id="group-threads-postprocess-only",
        pattern=r"threads drive\s+the SKAT-O post-processing instead",
        why="--nThreads drives three things on a group run, not one: the C++ region "
            "marker loop, the per-cell tests of passes B and D, and the SKAT-O pool",
        evidence="3eeb554e + e7b346ba; the new test asserts both printed lines "
                 "(tests/testthat/test-step2-pgen.R, the group-at-4-threads case)",
        canary="For group tests the marker loop stays on one thread and the threads drive "
               "the SKAT-O post-processing instead (about half of a group run",
    ),
    dict(
        id="cells-serial",
        pattern=r"(ultra-rare pass|BURDEN emitter)[^.]{0,120}still serial per cell",
        why="step 3b put the per-cell TESTS of both passes on the pool; what stays "
            "serial in cell order is the prepare and the write-back either side",
        evidence="e7b346ba; runCellTests in src/Main.cpp, called from the ultra-rare "
                 "pass and the BURDEN emitter; docs/state/handoff-2026-09-06-region-cells.md",
        canary="The ultra-rare pass and the BURDEN emitter are still serial per cell, "
               "which is what caps BURDEN at ~1.6x on 4 threads.",
    ),
    dict(
        id="burden-1.6x",
        pattern=r"caps BURDEN at ~?1\.6x",
        why="that was step 3's number and step 3b superseded it: ~1.84x (binary) / "
            "~1.71x (quantitative) at 4 threads by the paired ratios",
        evidence="harness/bench/NOTES.md, the step-3b section",
        canary="which is what caps BURDEN at ~1.6x on 4 threads",
    ),
    dict(
        id="multimodel-marker-only",
        pattern=r"multiTraits_File[^.]{0,240}marker path\s+only",
        why="a model manifest runs on the region path too, when every model shares "
            "one sample set",
        evidence="1e02a154 + 2ed00165; the guard is R/SAIGE_Test_main.R (isGroupTest, "
                 "the sample-set-group check); tests/testthat/test-step2-region-multimodel.R",
        canary="**Step 2 tests several null models in ONE genotype pass** "
               "(`--GMMATmodel_varianceRatio_multiTraits_File`, marker path only, each "
               "model's file byte-identical to its own single-model run)",
    ),
    dict(
        id="manifest-single-variant-only",
        pattern=r"Single-variant tests only",
        why="a group test takes a manifest as well; only more than ONE sample-set "
            "group among the models is an error",
        evidence="2ed00165; R/SAIGE_Test_main.R stops with 'needs every model fitted on "
                 "the SAME sample set'; tests/testthat/test-step2-region-multimodel.R",
        canary="Single-variant tests only: a groupFile with more than one model is an error.",
    ),
]

BENIGN = (
    "The region marker loop threads, the SKAT-O pool threads, and a model "
    "manifest runs on both paths. BURDEN reads ~1.84x at 4 threads. Nothing "
    "here is stale, so no rule may fire on it."
)


def targets():
    seen, out = set(), []
    for pat in SCANNED:
        for p in sorted(REPO.glob(pat)):
            if p.is_file() and p not in seen:
                seen.add(p)
                out.append(p)
    return out


# Blanked before matching, WITHOUT changing any offset, so the reported line
# number is still exact. Newlines because these claims wrap across lines
# (`marker path\n only` is one sentence, and a line-wise grep for it finds
# nothing -- which is how one of the five survived a search that had already
# found the other four). Quotes, backticks and commas because the same
# sentence in a CLI help string is a run of concatenated R literals:
# `"...threads drive ",\n      "the SKAT-O post-processing..."`, where the
# gap between two words is `",` newline spaces `"`. The first version of this
# tool missed exactly that hit while catching the other nine.
BLANK = set("\n\"'`,")


def flatten(text):
    """Blank the punctuation above to spaces, keeping a per-character line number."""
    line_of, buf, line = [], [], 1
    for ch in text:
        buf.append(" " if ch in BLANK else ch)
        line_of.append(line)
        if ch == "\n":
            line += 1
    return "".join(buf), line_of


def check():
    hits = 0
    for path in targets():
        text = path.read_text(encoding="utf-8", errors="replace")
        flat, line_of = flatten(text)
        lines = text.splitlines()
        for rule in RULES:
            for m in re.finditer(rule["pattern"], flat):
                n = line_of[m.start()]
                # A line (or the one before it, for a wrapped claim) may opt out.
                window = " ".join(lines[max(0, n - 2):n + 1])
                if ALLOW in window:
                    continue
                hits += 1
                rel = path.relative_to(REPO)
                print(f"DOCDRIFT: {rel}:{n}: [{rule['id']}] {m.group(0)[:80]!r}")
                print(f"          false since: {rule['why']}")
                print(f"          evidence:    {rule['evidence']}")
    if hits:
        print(f"docdrift: {hits} falsified statement(s); fix the doc, or mark the line "
              f"{ALLOW} if it quotes one on purpose")
        return 1
    print(f"docdrift OK: {len(RULES)} rules over {len(targets())} files, 0 hits")
    return 0


def selftest():
    """Positive control: every rule fires on its own canary. Negative control:
    no rule fires on benign text -- which is what catches a pattern written so
    loosely that it matches everything, the way a rule would fail silently."""
    problems = []
    for rule in RULES:
        flat, _ = flatten(rule["canary"])
        if not re.search(rule["pattern"], flat):
            problems.append(f"[{rule['id']}] does not match its own canary")
        if re.search(rule["pattern"], BENIGN):
            problems.append(f"[{rule['id']}] fires on benign text (pattern too loose)")
    ids = [r["id"] for r in RULES]
    if len(set(ids)) != len(ids):
        problems.append("duplicate rule ids")
    for p in problems:
        print("DOCDRIFT SELFTEST:", p)
    if problems:
        return 1
    print(f"docdrift selftest OK: {len(RULES)} rules, each fires on its canary "
          f"and none on benign text")
    return 0


def main(argv):
    if argv and argv[0] == "selftest":
        return selftest()
    if argv:
        print(__doc__)
        return 2
    return check()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
