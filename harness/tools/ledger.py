#!/usr/bin/env python3
r"""Query LEDGER.md without reading 180 KB of it.

The ledger is two markdown tables (#1-#39, then #40 on, split by an audit
section) whose rows are single lines of up to ~2,600 characters. Reading it
whole is expensive and scrolling it is useless. This prints one line per entry,
or one entry in full, wrapped.

  ledger.py                 every entry, one line each
  ledger.py open            only the OPEN ones
  ledger.py decided         the DECIDED ones: owner chose to leave, documented
  ledger.py 56              entry #56 in full, wrapped
  ledger.py grep ingest     entries whose text matches (case-insensitive)
  ledger.py lint            check every row parses; exit 1 and name the row if not

A row is one line of exactly five cells separated by UNESCAPED pipes. Any `|`
in the prose -- `|Beta|`, `a || b`, `x |= y`, even inside backticks -- must be
written `\|` (the GFM table escape, which also renders correctly). Row #66
once wrote a bare `|Beta_c|`: the row split into extra cells, the parser kept
the first five, so STATUS read as the wrong cell and `open` never listed it.
`harness/run_tests.sh` runs `lint` before every tier.
"""
import re, sys, textwrap, pathlib

REPO = pathlib.Path(__file__).resolve().parents[2]
LEDGER = REPO / "LEDGER.md"
COLS = ["#", "WHERE", "WHAT", "PROVENANCE", "STATUS"]
ROW = re.compile(r"^\| \d+ \|")
# Files that must not carry a hand-written open count: the count drifts (it
# said "36 open of 86" while the tool said 39 of 89). Quote `ledger.py open`.
NO_COUNT_FILES = ["CLAUDE.md", "docs/state/next.md"]
COUNT_PAT = re.compile(r"\b\d+ open of \d+\b")


SEP = re.compile(r"(?<!\\)\|")  # a cell separator is an UNESCAPED pipe


def split(line):
    s = line.strip()
    return [x.strip() for x in SEP.split(s[1:-1])]


def unescape(cell):
    return cell.replace("\\|", "|")


def rows():
    for line in LEDGER.read_text().splitlines():
        if ROW.match(line):
            c = split(line)
            if len(c) == 5:
                yield c


def status_of(cell):
    s = cell.lstrip("*").upper()
    if s.startswith("OPEN"):
        return "OPEN"
    if s.startswith("REMOVED"):
        return "REMOVED"
    if s.startswith("DECIDED"):
        return "DECIDED"
    return "FIXED" if s.startswith("FIXED") else "?"


def one_line(c):
    where = re.sub(r"`([^`]*)`", r"\1", unescape(c[1]))
    return f"#{c[0]:<3} {status_of(c[4]):<5} {where[:88]}"


def show(c):
    for name, cell in zip(COLS, c):
        print(f"\n=== {name} ===")
        print(textwrap.fill(unescape(cell), 96))


def lint():
    """Every row: exactly 5 cells, a parseable STATUS; IDs 1..N contiguous in
    file order; no hand-written open count in the files that load each session."""
    problems = []
    ids = []
    for n, line in enumerate(LEDGER.read_text().splitlines(), 1):
        if not ROW.match(line):
            continue
        c = split(line)
        rid = c[0]
        if len(c) != 5:
            problems.append(f"LEDGER.md:{n} row #{rid}: {len(c)} cells, expected 5 "
                            f"(a bare '|' in the prose? write it '\\|')")
            continue
        if status_of(c[4]) == "?":
            problems.append(f"LEDGER.md:{n} row #{rid}: STATUS does not start with OPEN, FIXED, REMOVED or DECIDED: "
                            f"{c[4][:60]!r}")
        ids.append(int(rid))
    expected = list(range(1, len(ids) + 1))
    if ids != expected:
        missing = sorted(set(expected) - set(ids))
        dups = sorted({i for i in ids if ids.count(i) > 1})
        order = [i for i, e in zip(ids, expected) if i != e][:3]
        problems.append(f"LEDGER.md: {len(ids)} parseable rows but ids are not 1..{len(ids)} "
                        f"in order (missing {missing}, duplicated {dups}, first out of place {order})")
    for rel in NO_COUNT_FILES:
        p = REPO / rel
        if not p.exists():
            continue
        for n, line in enumerate(p.read_text().splitlines(), 1):
            m = COUNT_PAT.search(line)
            if m:
                problems.append(f"{rel}:{n}: hand-written count {m.group(0)!r}; "
                                f"quote `ledger.py open` instead")
    for p in problems:
        print("LEDGER LINT:", p)
    if problems:
        return 1
    n_open = sum(1 for c in rows() if status_of(c[4]) == "OPEN")
    n_dec = sum(1 for c in rows() if status_of(c[4]) == "DECIDED")
    print(f"ledger lint OK: {len(ids)} rows, {n_open} open, {n_dec} decided")
    return 0


def main(argv):
    if argv and argv[0] == "lint":
        return lint()
    all_rows = list(rows())
    if not argv:
        for c in all_rows:
            print(one_line(c))
    elif argv[0] == "open":
        hits = [c for c in all_rows if status_of(c[4]) == "OPEN"]
        for c in hits:
            print(one_line(c))
        print(f"\n{len(hits)} open of {len(all_rows)}")
    elif argv[0] == "decided":
        hits = [c for c in all_rows if status_of(c[4]) == "DECIDED"]
        for c in hits:
            print(one_line(c))
        print(f"\n{len(hits)} decided of {len(all_rows)}")
    elif argv[0] == "grep":
        pat = " ".join(argv[1:]).lower()
        for c in all_rows:
            if pat in " ".join(c).lower():
                print(one_line(c))
    elif argv[0].isdigit():
        for c in all_rows:
            if c[0] == argv[0]:
                show(c)
                return 0
        print(f"no entry #{argv[0]}", file=sys.stderr)
        return 1
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
