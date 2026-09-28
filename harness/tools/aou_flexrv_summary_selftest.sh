#!/bin/bash
# Self-test for aou_flexrv_summary.sh on FABRICATED files (no AoU data): the
# fixed file must win over v4 for a spliced trait, an NA set p is skipped, a gene
# in one arm only carries NA on the other, sidecars are ignored, and the row
# count and both-arms count are what the fixtures say.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
mkdir -p "${T}/base/xx/flexrv_v4" "${T}/base/xx/flexrv_v4_fixed" "${T}/base/xx/skato_v4" "${T}/out"
cd "${T}"
python3 - <<'PY'
import gzip
H = "Region\tGroup\tmax_MAF\tWeight\tPvalue_Burden\tBETA_Burden\tSE_Burden\tMAC\tMAC_case\tMAC_control\tNumber_rare\tNumber_ultra_rare"
def fx(reg, w, p): return f"{reg}\tgrp\t0.001\t{w}\t{p}\t1\t1\t9\t1\t8\t0\t9"
def cau(reg, p): return f"{reg}\tCauchy\tNA\tNA\t{p}\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
v4a = [H, fx("G1","W1","0.5"), fx("G1","W2","1e-100"), cau("G1","1e-99"), fx("G2","W1","0.2"), cau("G2","0.3")]
fxa = [H, fx("G1","W1","0.5"), fx("G1","W2","0.001"), cau("G1","0.002"), fx("G2","W1","0.2"), cau("G2","0.3")]
v4b = [H, fx("G1","W1","0.04"), fx("G1","W2","NA"), cau("G1","0.05"), fx("G3","W1","0.7"), cau("G3","0.7")]
SH = "Region\tGroup\tmax_MAF\tPvalue\tPvalue_Burden\tPvalue_SKAT\tBETA_Burden\tSE_Burden\tMAC"
def sk(reg, grp, maf, p, pb, ps): return f"{reg}\t{grp}\t{maf}\t{p}\t{pb}\t{ps}\t1\t1\t9"
def skc(reg, p): return f"{reg}\tCauchy\tNA\t{p}\tNA\tNA\tNA\tNA\tNA"
ska = [SH, sk("G1","pLoF","0.001","0.01","0.02","0.03"), sk("G1","pLoF","0.01","0.005","0.5","0.004"), skc("G1","0.008"), sk("G4","pLoF","0.001","0.9","0.9","0.9"), skc("G4","0.9")]
skb = [SH, sk("G1","pLoF","0.001","0.3","0.3","0.3"), skc("G1","0.3")]
for path, lines in [("base/xx/flexrv_v4/fxv4-xx-traita.TraitA.xx.chr1-22.txt.gz", v4a), ("base/xx/flexrv_v4_fixed/fxv4-xx-traita.TraitA.xx.chr1-22.txt.gz", fxa),
                    ("base/xx/flexrv_v4/fxv4-xx-traitb.TraitB.xx.chr1-22.txt.gz", v4b),
                    ("base/xx/skato_v4/skv4-xx-traita.TraitA.xx.chr1-22.txt.gz", ska), ("base/xx/skato_v4/skv4-xx-traitb.TraitB.xx.chr1-22.txt.gz", skb)]:
    gzip.open(path, "wt").write("\n".join(lines) + "\n")
gzip.open("base/xx/flexrv_v4/fxv4-xx-traita.TraitA.xx.chr1-22.singleAssoc.txt.gz", "wt").write("x\n")
PY
ANC=xx BASE="${T}/base" OUT_DIR="${T}/out" bash "${HERE}/aou_flexrv_summary.sh" > run.log 2>&1 || { cat run.log; echo "SELFTEST FAILED: script exited non-zero"; exit 1; }
gzip -dc out/xx.flexrv_vs_skato.v4.tsv.gz > out.tsv
fail=0
check(){ if awk -F'\t' "$1" out.tsv | grep -q .; then echo "ok   $2"; else echo "FAIL $2"; fail=1; fi; }
check 'NR>1 && $2=="TraitA" && $3=="G1" && $4=="fixed" && $5=="0.002" && $6=="0.001" && $7=="W2" && $8=="2"' "fixed file wins over v4, best set found"
check 'NR>1 && $2=="TraitB" && $3=="G1" && $4=="v4" && $6=="0.04" && $7=="W1" && $8=="2"' "NA set p skipped, v4 used when no fixed file"
check 'NR>1 && $2=="TraitA" && $3=="G1" && $9=="0.008" && $10=="0.005" && $11=="pLoF@0.01" && $12=="0.02" && $13=="0.004"' "SKAT-O cauchy, best cell, best burden, best SKAT"
check 'NR>1 && $2=="TraitA" && $3=="G4" && $4=="NA" && $5=="NA" && $9=="0.9"' "gene in the SKAT-O arm only carries NA on the FlexRV side"
check 'NR>1 && $2=="TraitA" && $3=="G2" && $9=="NA"' "gene in the FlexRV arm only carries NA on the SKAT-O side"
[[ $(($(wc -l < out.tsv) - 1)) -eq 5 ]] && echo "ok   5 rows" || { echo "FAIL row count"; fail=1; }
grep -q "rows: 5 (trait, gene) pairs; 2 with both arms; traits: 2" run.log && echo "ok   summary line" || { echo "FAIL summary line"; fail=1; }
# the wide cell table: TraitA G1 has two cells (pLoF@0.001, pLoF@0.01), G4 one; TraitB G1 one
gzip -dc out/xx.skato_cells.v4.tsv.gz > cells.tsv
head -1 cells.tsv | grep -qx $'anc\ttrait\tgene\tburden:pLoF@0.001\tburden:pLoF@0.01\tskato:pLoF@0.001\tskato:pLoF@0.01' && echo "ok   cells header (cell list from the files)" || { echo "FAIL cells header: $(head -1 cells.tsv)"; fail=1; }
awk -F'\t' 'NR>1 && $2=="TraitA" && $3=="G1" && $4=="0.02" && $5=="0.5" && $6=="0.01" && $7=="0.005"' cells.tsv | grep -q . && echo "ok   cells row: both cells of G1, burden and skato" || { echo "FAIL cells row G1"; fail=1; }
awk -F'\t' 'NR>1 && $2=="TraitA" && $3=="G4" && $4=="0.9" && $5=="NA" && $7=="NA"' cells.tsv | grep -q . && echo "ok   cells row: a gene with one cell carries NA in the other" || { echo "FAIL cells row G4"; fail=1; }
[[ $(($(wc -l < cells.tsv) - 1)) -eq 3 ]] && echo "ok   3 cell rows (genes in the SKAT-O arm)" || { echo "FAIL cells row count"; fail=1; }
grep -q "cells table: 3 (trait, gene) rows x 2 cells" run.log && echo "ok   cells summary line" || { echo "FAIL cells summary line"; fail=1; }
(( fail )) && { echo "SELFTEST FAILED"; exit 1; }
echo "selftest OK: 12 checks"
