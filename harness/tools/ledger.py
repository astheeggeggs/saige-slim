#!/usr/bin/env python3
"""Query LEDGER.md without reading 104 KB of it.

The ledger is two markdown tables (#1-#39, then #40 on, split by an audit
section) whose rows are single lines of up to ~2,600 characters. Reading it
whole is expensive and scrolling it is useless. This prints one line per entry,
or one entry in full, wrapped.

  ledger.py                 every entry, one line each
  ledger.py open            only the OPEN ones
  ledger.py 56              entry #56 in full, wrapped
  ledger.py grep ingest     entries whose text matches (case-insensitive)
"""
import re, sys, textwrap, pathlib

LEDGER = pathlib.Path(__file__).resolve().parents[2] / "LEDGER.md"
COLS = ["#", "WHERE", "WHAT", "PROVENANCE", "STATUS"]


def rows():
    for line in LEDGER.read_text().splitlines():
        if re.match(r"^\| \d+ \|", line):
            c = [x.strip() for x in line.strip().strip("|").split("|")]
            if len(c) >= 5:
                yield c[:5]


def status_of(cell):
    s = cell.lstrip("*").upper()
    return "OPEN" if s.startswith("OPEN") else ("FIXED" if "FIXED" in s[:16] else "?")


def one_line(c):
    where = re.sub(r"`([^`]*)`", r"\1", c[1])
    return f"#{c[0]:<3} {status_of(c[4]):<5} {where[:88]}"


def show(c):
    for name, cell in zip(COLS, c):
        print(f"\n=== {name} ===")
        print(textwrap.fill(cell, 96))


def main(argv):
    all_rows = list(rows())
    if not argv:
        for c in all_rows:
            print(one_line(c))
    elif argv[0] == "open":
        hits = [c for c in all_rows if status_of(c[4]) == "OPEN"]
        for c in hits:
            print(one_line(c))
        print(f"\n{len(hits)} open of {len(all_rows)}")
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
