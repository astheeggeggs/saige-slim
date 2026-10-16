#!/bin/bash
# Self-test for flexrv_rand_scores.py on a FABRICATED group file: var and anno lines are
# copied byte for byte, every score line is replaced by score:RAND with one value in (0, 1)
# per variant, the same variant gets the same value in a second file and a second run (seeded
# by the variant ID), a different seed gives different values, and a score line whose length
# disagrees with its var line is refused.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
python3 - "${T}" <<'PY'
import sys, gzip
T = sys.argv[1]
lines = ["G1 var chr1:10:A:C chr1:20:G:T chr1:30:C:A", "G1 anno pLoF damaging_missense synonymous", "G1 score:AM 1.0 0.912 0.031",
         "G2 var chr1:40:T:C chr1:50:A:G", "G2 anno other_missense synonymous", "G2 score:AM 0.2 0.0"]
gzip.open(f"{T}/a.txt.gz", "wt").write("\n".join(lines) + "\n")
# the same first variant in another file (another chromosome's file would share no variants; a
# re-run of the same file must reproduce every value)
gzip.open(f"{T}/b.txt.gz", "wt").write("G9 var chr1:10:A:C chr9:1:A:T\nG9 anno pLoF synonymous\nG9 score:AM 0.5 0.5\n")
gzip.open(f"{T}/bad.txt.gz", "wt").write("G1 var chr1:10:A:C chr1:20:G:T\nG1 anno pLoF pLoF\nG1 score:AM 1.0\n")
PY
cd "${T}"
python3 "${HERE}/flexrv_rand_scores.py" a.txt.gz out1.txt.gz > log1 && python3 "${HERE}/flexrv_rand_scores.py" a.txt.gz out2.txt.gz > log2
python3 "${HERE}/flexrv_rand_scores.py" b.txt.gz outb.txt.gz > logb
python3 "${HERE}/flexrv_rand_scores.py" a.txt.gz out3.txt.gz --seed 7 > log3
fail=0
python3 - <<'PY' || fail=1
import gzip
a = gzip.open("a.txt.gz", "rt").read().splitlines(); o = gzip.open("out1.txt.gz", "rt").read().splitlines(); o2 = gzip.open("out2.txt.gz", "rt").read().splitlines()
ob = gzip.open("outb.txt.gz", "rt").read().splitlines(); o3 = gzip.open("out3.txt.gz", "rt").read().splitlines()
assert [l for l in a if " var " in l or " anno " in l] == [l for l in o if " var " in l or " anno " in l], "var/anno lines changed"
sc = {l.split()[0]: l.split()[2:] for l in o if " score:" in l}
assert all(l.split()[1] == "score:RAND" for l in o if " score:" in l), "score name"
assert len(sc["G1"]) == 3 and len(sc["G2"]) == 2, "value counts"
assert all(0 < float(v) < 1 for g in sc for v in sc[g]), "values in (0,1)"
assert len({float(v) for v in sc["G1"]}) == 3, "values vary"
assert o == o2, "a rerun reproduces every value"
scb = {l.split()[0]: l.split()[2:] for l in ob if " score:" in l}
assert scb["G9"][0] == sc["G1"][0], "the same variant gets the same value in another file"
sc3 = {l.split()[0]: l.split()[2:] for l in o3 if " score:" in l}
assert sc3["G1"] != sc["G1"], "a different seed gives different values"
print("ok   var/anno unchanged, score:RAND in (0,1), reproducible, variant-keyed, seed-sensitive")
PY
grep -q "2 regions, 2 score lines replaced by score:RAND" log1 && echo "ok   summary line" || { echo "FAIL summary: $(cat log1)"; fail=1; }
if python3 "${HERE}/flexrv_rand_scores.py" bad.txt.gz outbad.txt.gz > logbad 2>&1; then echo "FAIL a mismatched score line was accepted"; fail=1; else echo "ok   a score line of the wrong length is refused"; fi
(( fail )) && { echo "SELFTEST FAILED"; exit 1; }
echo "selftest OK: 8 checks"
