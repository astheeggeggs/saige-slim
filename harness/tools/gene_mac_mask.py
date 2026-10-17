#!/usr/bin/env python3
"""No gene is named beside a MAC rung below 20 (owner, 2026-10-16).

All of Us forbids reporting a count below 20 from real data, and a FlexRV set label such as
`sbeta20_MAC1` names an allele-count bin below 20; written beside a gene it says that gene has
such a cell. The rule (CLAUDE.md, Traps): genes in one list, set labels in another, the pairing
stays in the in-workspace table. This tool enforces it on the Markdown surfaces:

  rewrite [PATH...]   in place: in every *.md under the paths (default: the repo), in every
                      PARAGRAPH that names an Ensembl gene, a backticked set label's rung
                      `MAC1` / `MAC3` / `MAC5` / `MAC10` becomes `MACsparse`. Nothing outside
                      Markdown is touched (goldens, tests and R code carry the rung names and
                      must stay byte-identical), nothing in a paragraph without a gene, and
                      the rungs at or above 20 (`MAC25`, `MAC50`, ...) are left alone.
  scan [REV...]       read-only: report every (file, passage) in the working tree, or in each
                      given commit, where a gene and a sub-20 rung share a paragraph. Exit 1
                      if any. `scan $(git rev-list A..B)` checks a range of commits.
  selftest            the rules on fabricated text, both ways.

Run `scan` before pushing a Markdown change that names a gene.
"""
import os, re, subprocess, sys

GENE = re.compile(r"\bENSG\d{6,}\b")
RUNG = re.compile(r"(?<=MAC)(?:1|3|5|10)(?=\b|`)")          # the rung suffix inside a set label
# a rung is `MAC1` etc. whether it follows `_` (a set label) or a space; `\b` would miss `_MAC1`
_R = r"(?<![A-Za-z0-9])MAC(?:1|3|5|10)(?![0-9])"
PAIR = re.compile(r"\bENSG\d{6,}\b[^\n]{0,400}" + _R + r"|" + _R + r"[^\n]{0,400}\bENSG\d{6,}\b")


def paragraphs(text):
    """(start, end) of every blank-line-separated paragraph."""
    out, start = [], 0
    for m in re.finditer(r"\n\s*\n", text):
        out.append((start, m.start()))
        start = m.end()
    out.append((start, len(text)))
    return out


def flat(s):
    return re.sub(r"(?<!\n)\n(?!\n)", " ", s)


def rewrite_text(text):
    pieces, last, n = [], 0, 0
    for a, b in paragraphs(text):
        para = text[a:b]
        if GENE.search(para) and RUNG.search(para):
            new = RUNG.sub("sparse", para)
            n += 1
        else:
            new = para
        pieces.append(text[last:a]); pieces.append(new); last = b
    pieces.append(text[last:])
    return "".join(pieces), n


def offending(text):
    return [m.group(0) for a, b in paragraphs(text) for m in PAIR.finditer(flat(text[a:b]))]


def md_files(paths):
    for p in paths:
        if os.path.isfile(p):
            if p.endswith(".md"): yield p
            continue
        for root, dirs, files in os.walk(p):
            dirs[:] = [d for d in dirs if d not in (".git", ".pixi", "node_modules", "site", "work")]
            for f in files:
                if f.endswith(".md"): yield os.path.join(root, f)


def cmd_rewrite(paths):
    changed = 0
    for f in md_files(paths or ["."]):
        with open(f, encoding="utf-8", errors="surrogateescape") as fh: t = fh.read()
        new, n = rewrite_text(t)
        if n:
            with open(f, "w", encoding="utf-8", errors="surrogateescape") as fh: fh.write(new)
            changed += 1
            print(f"rewrote {f}: {n} paragraph(s)", file=sys.stderr)
    return 0


def cmd_scan(revs):
    bad = 0
    if not revs:
        for f in md_files(["."]):
            with open(f, encoding="utf-8", errors="surrogateescape") as fh: t = fh.read()
            for o in offending(t):
                bad += 1; print(f"{f}: {o[:160]}")
    for r in revs:
        files = subprocess.run(["git", "ls-tree", "-r", "--name-only", r], capture_output=True, text=True).stdout.split("\n")
        for f in files:
            if not f.endswith(".md"): continue
            t = subprocess.run(["git", "show", f"{r}:{f}"], capture_output=True, text=True, errors="surrogateescape").stdout
            for o in offending(t):
                bad += 1; print(f"{r[:8]}:{f}: {o[:160]}")
        msg = subprocess.run(["git", "log", "-1", "--format=%B", r], capture_output=True, text=True).stdout
        for o in offending(msg):
            bad += 1; print(f"{r[:8]} COMMIT MESSAGE: {o[:160]}")
    print(f"{bad} offending passage(s)" + (f" over {len(revs)} commit(s)" if revs else ""))
    return 1 if bad else 0


def selftest():
    t = ("intro\n\nENSG90000000003\n`AM_sbeta20_MAC1` p 2.2e-53, eas MatHem ENSG90000000004 `AM_sbeta10_MAC5`\n\n"
         "| MatHem | ENSG90000000005 | 3.7e-14 | `AM_sbeta5_MAC1` |\n| X | ENSG90000000007 | 8.4e-9 | `AM_sbeta1.5_MAC5` |\n\n"
         "ENSG90000000002 `AM_sbeta1.5_MAC25` stays, ENSG90000000008 `AM_cbrt0.9_MAC10` goes\n\n"
         "no gene here: `sbeta20_MAC1` and `MAC3` are the sparse rungs, MAC10 too\n\n"
         "the genes, Height: ENSG90000000006; three of the four best sets `sbeta20_MAC1`\n")
    assert len(offending(t)) == 4, offending(t)
    new, n = rewrite_text(t)
    assert n == 4, n
    assert offending(new) == [], offending(new)
    assert "`AM_sbeta20_MACsparse`" in new and "`AM_cbrt0.9_MACsparse`" in new and "`AM_sbeta1.5_MAC25` stays" in new
    assert "no gene here: `sbeta20_MAC1` and `MAC3` are the sparse rungs, MAC10 too" in new   # untouched: no gene
    assert "MAC250" not in new
    # idempotent
    assert rewrite_text(new) == (new, 0)
    # a count that is not a rung is untouched
    t2 = "ENSG90000000001 has MAC 15 on the set and `MAC25`\n"
    assert rewrite_text(t2) == (t2, 0) and offending(t2) == []
    print("selftest OK: 9 checks")
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or a[0] not in ("rewrite", "scan", "selftest"):
        sys.exit(__doc__)
    sys.exit({"rewrite": lambda: cmd_rewrite(a[1:]), "scan": lambda: cmd_scan(a[1:]), "selftest": selftest}[a[0]]())
