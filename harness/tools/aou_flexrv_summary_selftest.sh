#!/bin/bash
# Self-test for aou_flexrv_summary.sh on FABRICATED files (no AoU data): the
# verified v7 file must win over fixed, and fixed over v4, the exactByWeight
# sidecar must not be read as a results file, an NA set p is skipped, a gene
# in one arm only carries NA on the other, sidecars are ignored, and the row
# count and both-arms count are what the fixtures say. A second, v6 layout (2026-10-15):
# log10 columns, the am_moderate arm into the cells table only, the floor table.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
mkdir -p "${T}/base/xx/flexrv_v4" "${T}/base/xx/flexrv_v4_fixed" "${T}/base/xx/flexrv_v7_checked" "${T}/base/xx/skato_v4" "${T}/out"
cd "${T}"
python3 - <<'PY'
import gzip
H = "Region\tGroup\tmax_MAF\tWeight\tPvalue_Burden\tBETA_Burden\tSE_Burden\tMAC\tMAC_case\tMAC_control\tNumber_rare\tNumber_ultra_rare"
def fx(reg, w, p): return f"{reg}\tgrp\t0.001\t{w}\t{p}\t1\t1\t9\t1\t8\t0\t9"
def cau(reg, p): return f"{reg}\tCauchy\tNA\tNA\t{p}\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
v4a = [H, fx("G1","W1","0.5"), fx("G1","W2","1e-100"), cau("G1","1e-99"), fx("G2","W1","0.2"), cau("G2","0.3")]
fxa = [H, fx("G1","W1","0.5"), fx("G1","W2","0.001"), cau("G1","0.002"), fx("G2","W1","0.2"), cau("G2","0.3")]
v4b = [H, fx("G1","W1","0.04"), fx("G1","W2","NA"), cau("G1","0.05"), fx("G3","W1","0.7"), cau("G3","0.7")]
v7a = [H, fx("G1","W1","0.5"), fx("G1","W2","0.01"), cau("G1","0.02"), fx("G2","W1","0.2"), cau("G2","0.3")]   # the complete rerun (LEDGER #172)
v4c = [H, fx("G1","W1","0.6"), fx("G1","W2","1e-50"), cau("G1","1e-49")]
fxc = [H, fx("G1","W1","0.6"), fx("G1","W2","0.003"), cau("G1","0.006")]
SH = "Region\tGroup\tmax_MAF\tPvalue\tPvalue_Burden\tPvalue_SKAT\tBETA_Burden\tSE_Burden\tMAC"
def sk(reg, grp, maf, p, pb, ps): return f"{reg}\t{grp}\t{maf}\t{p}\t{pb}\t{ps}\t1\t1\t9"
def skc(reg, p): return f"{reg}\tCauchy\tNA\t{p}\tNA\tNA\tNA\tNA\tNA"
ska = [SH, sk("G1","pLoF","0.001","0.01","0.02","0.03"), sk("G1","pLoF","0.01","0.005","0.5","0.004"), skc("G1","0.008"), sk("G4","pLoF","0.001","0.9","0.9","0.9"), skc("G4","0.9")]
skb = [SH, sk("G1","pLoF","0.001","0.3","0.3","0.3"), skc("G1","0.3")]
for path, lines in [("base/xx/flexrv_v4/fxv4-xx-traita.TraitA.xx.chr1-22.txt.gz", v4a), ("base/xx/flexrv_v4_fixed/fxv4-xx-traita.TraitA.xx.chr1-22.txt.gz", fxa),
                    ("base/xx/flexrv_v7_checked/fxv4-xx-traita.TraitA.xx.chr1-22.txt.gz", v7a),
                    ("base/xx/flexrv_v4/fxv4-xx-traitb.TraitB.xx.chr1-22.txt.gz", v4b),
                    ("base/xx/flexrv_v4/fxv4-xx-traitc.TraitC.xx.chr1-22.txt.gz", v4c), ("base/xx/flexrv_v4_fixed/fxv4-xx-traitc.TraitC.xx.chr1-22.txt.gz", fxc),
                    ("base/xx/skato_v4/skv4-xx-traita.TraitA.xx.chr1-22.txt.gz", ska), ("base/xx/skato_v4/skv4-xx-traitb.TraitB.xx.chr1-22.txt.gz", skb)]:
    gzip.open(path, "wt").write("\n".join(lines) + "\n")
gzip.open("base/xx/flexrv_v4/fxv4-xx-traita.TraitA.xx.chr1-22.singleAssoc.txt.gz", "wt").write("x\n")
gzip.open("base/xx/flexrv_v7_checked/fxv4-xx-traita.TraitA.xx.chr1-22.exactByWeight.txt.gz", "wt").write("Kind\tRegion\tGroup\tmax_MAF\tWeight\tMAC\tn_eff\tn_material\nburden\tG1\tgrp\t0.001\tW2\t9\t1.0\t1\n")
PY
ANC=xx BASE="${T}/base" OUT_DIR="${T}/out" bash "${HERE}/aou_flexrv_summary.sh" > run.log 2>&1 || { cat run.log; echo "SELFTEST FAILED: script exited non-zero"; exit 1; }
gzip -dc out/xx.flexrv_vs_skato.v4.tsv.gz > out.tsv
fail=0
check(){ if awk -F'\t' "$1" out.tsv | grep -q .; then echo "ok   $2"; else echo "FAIL $2"; fail=1; fi; }
check 'NR>1 && $2=="TraitA" && $3=="G1" && $4=="v7" && $5=="0.02" && $6=="0.01" && $7=="W2" && $8=="2"' "verified v7 file wins over fixed and v4, best set found"
check 'NR>1 && $2=="TraitC" && $3=="G1" && $4=="fixed" && $5=="0.006" && $6=="0.003" && $7=="W2" && $8=="2"' "fixed file wins over v4 when there is no v7"
check 'NR>1 && $2=="TraitB" && $3=="G1" && $4=="v4" && $6=="0.04" && $7=="W1" && $8=="2"' "NA set p skipped, v4 used when no fixed file"
check 'NR>1 && $2=="TraitA" && $3=="G1" && $9=="0.008" && $10=="0.005" && $11=="pLoF@0.01" && $12=="0.02" && $13=="0.004"' "SKAT-O cauchy, best cell, best burden, best SKAT"
check 'NR>1 && $2=="TraitA" && $3=="G4" && $4=="NA" && $5=="NA" && $9=="0.9"' "gene in the SKAT-O arm only carries NA on the FlexRV side"
check 'NR>1 && $2=="TraitA" && $3=="G2" && $9=="NA"' "gene in the FlexRV arm only carries NA on the SKAT-O side"
[[ $(($(wc -l < out.tsv) - 1)) -eq 6 ]] && echo "ok   6 rows (the exactByWeight sidecar is not a results file)" || { echo "FAIL row count"; fail=1; }
grep -q "rows: 6 (trait, gene) pairs; 2 with both arms; traits: 3" run.log && echo "ok   summary line" || { echo "FAIL summary line"; fail=1; }
grep -q "1 v7" run.log && echo "ok   one v7 file listed" || { echo "FAIL v7 listing"; fail=1; }
# the wide cell table: TraitA G1 has two cells (pLoF@0.001, pLoF@0.01), G4 one; TraitB G1 one
gzip -dc out/xx.skato_cells.v4.tsv.gz > cells.tsv
head -1 cells.tsv | grep -qx $'anc\ttrait\tgene\tburden:pLoF@0.001\tburden:pLoF@0.01\tskato:pLoF@0.001\tskato:pLoF@0.01' && echo "ok   cells header (cell list from the files)" || { echo "FAIL cells header: $(head -1 cells.tsv)"; fail=1; }
awk -F'\t' 'NR>1 && $2=="TraitA" && $3=="G1" && $4=="0.02" && $5=="0.5" && $6=="0.01" && $7=="0.005"' cells.tsv | grep -q . && echo "ok   cells row: both cells of G1, burden and skato" || { echo "FAIL cells row G1"; fail=1; }
awk -F'\t' 'NR>1 && $2=="TraitA" && $3=="G4" && $4=="0.9" && $5=="NA" && $7=="NA"' cells.tsv | grep -q . && echo "ok   cells row: a gene with one cell carries NA in the other" || { echo "FAIL cells row G4"; fail=1; }
[[ $(($(wc -l < cells.tsv) - 1)) -eq 3 ]] && echo "ok   3 cell rows (genes in the SKAT-O arm)" || { echo "FAIL cells row count"; fail=1; }
grep -q "cells table: 3 (trait, gene) rows x 2 cells" run.log && echo "ok   cells summary line" || { echo "FAIL cells summary line"; fail=1; }
# ---- v6 layout (2026-10-15): log10 columns, the am_moderate arm, the floor table ------
mkdir -p base/xx/skato_v6 base/xx/am_v6 base/xx/flexrv_v6 out6
python3 - <<'PY'
import gzip
SH = "Region\tGroup\tmax_MAF\tPvalue\tPvalue_Burden\tPvalue_SKAT\tBETA_Burden\tSE_Burden\tlog10_Pvalue\tlog10_Pvalue_Burden\tlog10_Pvalue_SKAT\tMAC"
def sk(reg, grp, maf, p, l10): return f"{reg}\t{grp}\t{maf}\t{p}\t{p}\t{p}\t1\t1\t{l10}\t{l10}\t{l10}\t9"
def skc(reg, p, l10): return f"{reg}\tCauchy\tNA\t{p}\t{p}\t{p}\tNA\tNA\t{l10}\t{l10}\t{l10}\tNA"
# G1: an ordinary cell and one PAST THE FLOOR (printed 0, log10 -602.1); G5 ordinary
sk6 = [SH, sk("G1","pLoF","0.001","0.01","-2"), sk("G1","pLoF;damaging_missense;spliceai2_gt040","0.001","0","-602.1"), skc("G1","0","-601.2"),
       sk("G5","synonymous","0.001","0.5","-0.30103"), skc("G5","0.5","-0.30103")]
am6 = [SH, sk("G1","am_moderate","0.001","0.002","-2.69897"), skc("G1","0.002","-2.69897"), sk("G5","am_moderate","0.01","0.4","-0.39794"), skc("G5","0.4","-0.39794")]
FH = "Region\tGroup\tmax_MAF\tWeight\tPvalue_Burden\tBETA_Burden\tSE_Burden\tlog10_Pvalue_Burden\tMAC"
def fx(reg, w, p, l10): return f"{reg}\tgrp\t0.001\t{w}\t{p}\t1\t1\t{l10}\t9"
def cau(reg, p, l10): return f"{reg}\tCauchy\tNA\tNA\t{p}\tNA\tNA\t{l10}\tNA"
fx6 = [FH, fx("G1","W1","0.5","-0.30103"), fx("G1","W2","0","-350.5"), cau("G1","0","-349.8"), fx("G5","W1","0.2","-0.69897"), cau("G5","0.3","-0.52288")]
for path, lines in [("base/xx/skato_v6/sk6-xx-traita.TraitA.xx.chr1-22.txt.gz", sk6), ("base/xx/am_v6/am6-xx-traita.TraitA.xx.chr1-22.txt.gz", am6),
                    ("base/xx/flexrv_v6/fx6-xx-traita.TraitA.xx.chr1-22.txt.gz", fx6)]:
    gzip.open(path, "wt").write("\n".join(lines) + "\n")
gzip.open("base/xx/am_v6/am6-xx-traita.TraitA.xx.chr1-22.singleAssoc.txt.gz", "wt").write("# no singleAssoc output: SINGLE_IN_GROUP=FALSE\n")
# a STRAY from another run tag in the v6 directory: excluded by its basename, named in the log
gzip.open("base/xx/flexrv_v6/fxv5-xx-traita.TraitA.xx.chr1-22.txt.gz", "wt").write("\n".join([FH, fx("G1","W1","1e-30","-30"), cau("G1","1e-29","-29")]) + "\n")
# a LEDGER #181 sidecar on the skato arm: two named cells (never a results file)
gzip.open("base/xx/skato_v6/sk6-xx-traita.TraitA.xx.chr1-22.spaRegionNA.txt.gz", "wt").write(
    "Kind\tRegion\tGroup\tmax_MAF\tWeight\tn_carriers\tn_pinned_dropped\tq_adj\tlog_p_normal_reported\tkernel_value\tdump\n"
    "region_rescale\tG1\tpLoF\t0.001\tNA\t9\t0\t3.5\t-4.2\tNA\t\nregion_rescale\tG5\tpLoF\t0.001\tNA\t3\t0\t2.1\t-1.1\tNaN\t\n")
# chrStatus sidecars: two clean, one that lost chr8 (the job still delivered)
gzip.open("base/xx/skato_v6/sk6-xx-traita.TraitA.xx.chr1-22.chrStatus.txt.gz", "wt").write("requested\t1 22\nfailed\tnone\nrefused_by_gate\tnone\nskipped_no_region\tnone\n")
gzip.open("base/xx/am_v6/am6-xx-traita.TraitA.xx.chr1-22.chrStatus.txt.gz", "wt").write("requested\t1 22\nfailed\tnone\nrefused_by_gate\tnone\nskipped_no_region\tnone\n")
gzip.open("base/xx/flexrv_v6/fx6-xx-traita.TraitA.xx.chr1-22.chrStatus.txt.gz", "wt").write("requested\t1 22\nfailed\t8\nrefused_by_gate\tnone\nskipped_no_region\tnone\n")
PY
ANC=xx BASE="${T}/base" OUT_DIR="${T}/out6" RUNTAG=v6 bash "${HERE}/aou_flexrv_summary.sh" > run6.log 2>&1 || { cat run6.log; echo "SELFTEST FAILED: v6 run exited non-zero"; exit 1; }
gzip -dc out6/xx.skato_cells.v6.tsv.gz > cells6.tsv
head -1 cells6.tsv | tr '\t' '\n' | grep -qx 'burden:am_moderate@0.001' && echo "ok   v6: am_moderate cells in the wide table" || { echo "FAIL v6 am cell column: $(head -1 cells6.tsv)"; fail=1; }
awk -F'\t' 'NR>1 && $2=="TraitA" && $3=="G1"' cells6.tsv | grep -q $'\t0.002\t' && echo "ok   v6: am_moderate value on the gene row" || { echo "FAIL v6 am value"; fail=1; }
gzip -dc out6/xx.flexrv_vs_skato.v6.tsv.gz > out6.tsv
awk -F'\t' 'NR>1 && $3=="G1" && $10=="0" && $11=="pLoF;damaging_missense;spliceai2_gt040@0.001"' out6.tsv | grep -q . && echo "ok   v6: SKAT-O best cell is the floored one (p printed 0)" || { echo "FAIL v6 best cell"; fail=1; }
awk -F'\t' 'NR>1 && $3=="G1" && $9!="NA"' out6.tsv | grep -q . && echo "ok   v6: am_moderate did NOT enter the SKAT-O per-gene best (G1 has its own skato row)" || { echo "FAIL v6 am leaked"; fail=1; }
gzip -dc out6/xx.floor.v6.tsv.gz > floor6.tsv
[[ $(($(wc -l < floor6.tsv) - 1)) -eq 4 ]] && echo "ok   v6: floor table has the 4 planted rows (2 skato incl. Cauchy, 2 flexrv incl. Cauchy)" || { echo "FAIL floor rows: $(cat floor6.tsv)"; fail=1; }
head -2 floor6.tsv | tail -1 | grep -q $'skato\tTraitA\tG1\tpLoF;damaging_missense;spliceai2_gt040@0.001\t-602.1' && echo "ok   v6: floor table sorted by log10 p, the deepest first" || { echo "FAIL floor order: $(head -2 floor6.tsv)"; fail=1; }
grep -q "floor table: 4 row(s)" run6.log && echo "ok   v6: floor summary line" || { echo "FAIL floor summary"; fail=1; }
grep -q "am_moderate: 1 files" run6.log && echo "ok   v6: am listing line (the singleAssoc placeholder is not a results file)" || { echo "FAIL am listing"; fail=1; }
[[ $(gzip -dc out/xx.floor.v4.tsv.gz | wc -l | tr -d " ") == 1 ]] && echo "ok   v4: floor table is header-only (no log10 columns)" || { echo "FAIL v4 floor"; fail=1; }
grep -q "^chrStatus: 3 jobs scanned, 1 with a failed chromosome: flexrv/fx6-xx-traita(chr 8)" run6.log && echo "ok   v6: chromosome scan names the job that lost chr8" || { echo "FAIL chrStatus line: $(grep chrStatus run6.log)"; fail=1; }
gzip -dc out6/xx.chrstatus.v6.tsv.gz | awk -F'\t' 'NR>1 && $3=="fx6-xx-traita" && $5=="8"' | grep -q . && echo "ok   v6: chrstatus table row" || { echo "FAIL chrstatus table"; fail=1; }
grep -q "^chrStatus: 0 jobs scanned, 0 with a failed chromosome$" run.log && echo "ok   v4: no sidecars, scan reports 0/0" || { echo "FAIL v4 chrStatus line: $(grep chrStatus run.log)"; fail=1; }
grep -q "^EXCLUDED 1 file(s)" run6.log && grep -q "excluded (not fx6-\*): .*fxv5-xx-traita" run6.log && echo "ok   v6: the stray fxv5 file is excluded and named" || { echo "FAIL stray exclusion: $(grep -i exclud run6.log)"; fail=1; }
awk -F'\t' 'NR>1 && $2=="TraitA" && $3=="G1" && $5=="0"' out6.tsv | grep -q . && echo "ok   v6: the real fx6 file's Cauchy p (0) was used, not the stray's (1e-29)" || { echo "FAIL stray leaked into the join"; fail=1; }
grep -q "^flexrv: 1 v4 files" run6.log && echo "ok   v6: one FlexRV results file counted" || { echo "FAIL v6 flexrv count: $(grep '^flexrv' run6.log)"; fail=1; }
grep -q "^spaRegionNA: 2 rescale saddlepoint(s) returned nothing (LEDGER #181); skato/TraitA 2" run6.log && echo "ok   v6: the #181 sidecar is scanned and counted by arm/trait" || { echo "FAIL spaRegionNA line: $(grep spaRegionNA run6.log)"; fail=1; }
gzip -dc out6/xx.spaRegionNA.v6.tsv.gz | awk -F'\t' 'NR>1 && $3=="sk6-xx-traita" && $6=="G5" && $14=="NaN"' | grep -q . && echo "ok   v6: spaRegionNA table row" || { echo "FAIL spaRegionNA table"; fail=1; }
grep -q "^spaRegionNA: 0 rescale" run.log && echo "ok   v4: no #181 sidecar, 0 reported" || { echo "FAIL v4 spaRegionNA"; fail=1; }
grep -q "^flexrv: 1 v4 files, 0 fixed, 0 v7; skato: 1 files" run6.log && echo "ok   v6: the #181 sidecar is not counted as a results file" || { echo "FAIL sidecar counted"; fail=1; }
# a control fleet named <tag>-... in <anc>/<tag>/ (the random-score arms): FX_DIR + FX_PREFIX
mkdir -p base/xx/fxr_syn outr
cp "base/xx/flexrv_v6/fx6-xx-traita.TraitA.xx.chr1-22.txt.gz" "base/xx/fxr_syn/fxr_syn-xx-traita.TraitA.xx.chr1-22.txt.gz"
ANC=xx BASE="${T}/base" OUT_DIR="${T}/outr" RUNTAG=fxr_syn FX_DIR=fxr_syn FX_PREFIX=fxr_syn- bash "${HERE}/aou_flexrv_summary.sh" > runr.log 2>&1 || { cat runr.log; echo "SELFTEST FAILED: fxr_syn run exited non-zero"; exit 1; }
grep -q "^flexrv: 1 v4 files, 0 fixed, 0 v7; skato: 0 files" runr.log && [[ -s outr/xx.flexrv_vs_skato.fxr_syn.tsv.gz ]] && echo "ok   control tag: FX_DIR/FX_PREFIX read the <tag>/ directory, SKAT-O absent" || { echo "FAIL control tag: $(grep '^flexrv' runr.log)"; fail=1; }
ANC=xx BASE="${T}/base" OUT_DIR="${T}/outr0" RUNTAG=fxr_syn bash "${HERE}/aou_flexrv_summary.sh" > runr0.log 2>&1 || true
grep -q "^flexrv: 0 v4 files" runr0.log && echo "ok   control tag without the overrides reads nothing (the 2026-10-16 mistake)" || { echo "FAIL control tag w/o overrides: $(grep '^flexrv' runr0.log)"; fail=1; }
# both arms under one control tag (the simulated null): SK_DIR/SK_PREFIX beside FX_DIR/FX_PREFIX
mkdir -p base/xx/sim_sk base/xx/sim_fx outs
cp "base/xx/flexrv_v6/fx6-xx-traita.TraitA.xx.chr1-22.txt.gz" "base/xx/sim_fx/sim_fx-xx-traita.TraitA.xx.chr1-22.txt.gz"
for f in base/xx/skato_v6/sk6-xx-traita.TraitA.xx.chr1-22*.txt.gz; do cp "$f" "base/xx/sim_sk/$(basename "$f" | sed 's/^sk6-/sim_sk-/')"; done
ANC=xx BASE="${T}/base" OUT_DIR="${T}/outs" RUNTAG=simnull FX_DIR=sim_fx FX_PREFIX=sim_fx- SK_DIR=sim_sk SK_PREFIX=sim_sk- bash "${HERE}/aou_flexrv_summary.sh" > runs.log 2>&1 || { cat runs.log; echo "SELFTEST FAILED: simnull run exited non-zero"; exit 1; }
if grep -q "^flexrv: 1 v4 files, 0 fixed, 0 v7; skato: 1 files" runs.log && gzip -dc outs/xx.flexrv_vs_skato.simnull.tsv.gz | awk -F'\t' 'NR>1 && $3=="G1" && $5!="NA" && $9!="NA"' | grep -q .; then echo "ok   simnull tag: both arms read through the overrides, joined on the gene"; else echo "FAIL simnull tag: $(grep '^flexrv' runs.log)"; fail=1; fi
# all three arms under one control tag (the October fleet): AM_DIR/AM_PREFIX beside the other two
mkdir -p base/xx/oct_am outo
for f in base/xx/am_v6/am6-xx-traita.TraitA.xx.chr1-22*.txt.gz; do cp "$f" "base/xx/oct_am/$(basename "$f" | sed 's/^am6-/oct_am-/')"; done
ANC=xx BASE="${T}/base" OUT_DIR="${T}/outo" RUNTAG=oct FX_DIR=sim_fx FX_PREFIX=sim_fx- SK_DIR=sim_sk SK_PREFIX=sim_sk- AM_DIR=oct_am AM_PREFIX=oct_am- bash "${HERE}/aou_flexrv_summary.sh" > runo.log 2>&1 || { cat runo.log; echo "SELFTEST FAILED: oct run exited non-zero"; exit 1; }
if grep -q "^am_moderate: 1 files" runo.log && gzip -dc outo/xx.skato_cells.oct.tsv.gz | head -1 | grep -q "am_moderate@"; then echo "ok   oct tag: AM_DIR/AM_PREFIX read the am_moderate arm into the cells table"; else echo "FAIL oct tag am arm: $(grep '^am_moderate' runo.log)"; fail=1; fi
(( fail )) && { echo "SELFTEST FAILED"; exit 1; }
echo "selftest OK: 35 checks"
