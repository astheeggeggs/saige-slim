#!/bin/bash
# Self-test for aou_v6_vs_v4.py on FABRICATED summary tables (no AoU data). Planted:
# pLoF identical between old and new except THREE cells moved past 1e-5 (the cross-image
# control must report exactly those three moves and a Spearman near 1); synonymous
# uniform in both (fraction below 0.05 near 0.05); the new splice mask equals the new
# missense mask except FIVE genes made significant (the splice contrast must report
# +5 hits, only-second 5); am_moderate with two hits; FlexRV pooled with four; one floor
# row. Needs a python with pandas, matplotlib and scipy (the pyplot image, or a venv);
# pass it as $1, default python3.
set -euo pipefail
PY="${1:-python3}"; HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
"${PY}" - "${T}" <<'PYEOF'
import sys, gzip, numpy as np
T = sys.argv[1]; rng = np.random.default_rng(7); G = 400; traits = ["QuantX", "Asth"]
def u(): return rng.uniform(0.001, 1)
def write_cells(path, masks, pfun):
    mafs = ["0.0001", "0.001", "0.01"]
    cols = [f"burden:{m}@{f}" for m in masks for f in mafs] + [f"skato:{m}@{f}" for m in masks for f in mafs]
    with gzip.open(path, "wt") as fh:
        fh.write("anc\ttrait\tgene\t" + "\t".join(cols) + "\n")
        for tr in traits:
            for g in range(G):
                vals = [f"{pfun(tr, g, c):.6g}" for c in cols]
                fh.write(f"xx\t{tr}\tG{g}\t" + "\t".join(vals) + "\n")
def write_sum(path, tag, pflex):
    with gzip.open(path, "wt") as fh:
        fh.write("anc\ttrait\tgene\tsource\tp_flexrv\tp_flexrv_best\tset_best\tn_sets\tp_skato\tp_skato_best\tcell_best\tp_burden_best\tp_skat_best\n")
        for tr in traits:
            for g in range(G):
                fh.write(f"xx\t{tr}\tG{g}\t{tag}\t{pflex(tr, g):.6g}\tNA\tNA\t112\tNA\tNA\tNA\tNA\tNA\n")
# one fixed draw per (trait, gene, mask-kind) so old and new agree where they should
base = {(tr, g, k): u() for tr in traits for g in range(G) for k in ("plof", "syn", "dmg", "om", "dm", "all")}
old_masks = ["pLoF", "damaging_missense", "other_missense", "synonymous", "pLoF;damaging_missense", "pLoF;damaging_missense;other_missense;synonymous"]
new_masks = ["pLoF", "damaging_missense", "pLoF;damaging_missense", "pLoF;damaging_missense;spliceai2_gt040", "synonymous", "other_missense",
             "pLoF;damaging_missense;spliceai2_gt040;other_missense;synonymous", "am_moderate"]
kind = {"pLoF": "plof", "synonymous": "syn", "pLoF;damaging_missense": "dmg", "other_missense": "om", "damaging_missense": "dm",
        "pLoF;damaging_missense;other_missense;synonymous": "all", "pLoF;damaging_missense;spliceai2_gt040": "dmg",
        "pLoF;damaging_missense;spliceai2_gt040;other_missense;synonymous": "all", "am_moderate": "dm"}
MOVED = {("QuantX", 0), ("QuantX", 1), ("Asth", 2)}          # pLoF cells that moved past 1e-5 between images
SPL = {("QuantX", 10), ("QuantX", 11), ("QuantX", 12), ("Asth", 13), ("Asth", 14)}   # genes the splice label makes significant
AM = {("QuantX", 20), ("Asth", 21)}
def p_old(tr, g, c):
    m = c.split(":", 1)[1].rsplit("@", 1)[0]; p = base[(tr, g, kind[m])]
    if m == "pLoF" and (tr, g) in MOVED: p = 1e-6
    return p
def p_new(tr, g, c):
    m, maf = c.split(":", 1)[1].rsplit("@", 1); p = base[(tr, g, kind[m])]
    if m == "pLoF" and (tr, g) in MOVED: p = 1e-9             # three decades moved
    if m == "pLoF;damaging_missense;spliceai2_gt040" and maf == "0.001" and (tr, g) in SPL: p = 1e-8
    if m == "am_moderate" and maf == "0.001" and (tr, g) in AM: p = 1e-7
    return p
write_cells(f"{T}/xx.skato_cells.v4.tsv.gz", old_masks, p_old)
write_cells(f"{T}/xx.skato_cells.v6.tsv.gz", new_masks, p_new)
write_sum(f"{T}/xx.flexrv_vs_skato.v4.tsv.gz", "v4", lambda tr, g: u())
write_sum(f"{T}/xx.flexrv_vs_skato.v6.tsv.gz", "v6", lambda tr, g: 1e-7 if (tr, g) in {("QuantX", 10), ("QuantX", 30), ("Asth", 31), ("Asth", 32)} else u())
with gzip.open(f"{T}/xx.floor.v6.tsv.gz", "wt") as fh:
    fh.write("anc\tarm\ttrait\tgene\tcell\tlog10_p\nxx\tskato\tQuantX\tG0\tpLoF@0.001\t-602.1\n")
PYEOF
ANC=xx IN_SUM_NEW="${T}/xx.flexrv_vs_skato.v6.tsv.gz" IN_CELLS_NEW="${T}/xx.skato_cells.v6.tsv.gz" IN_FLOOR_NEW="${T}/xx.floor.v6.tsv.gz" \
  IN_SUM_OLD="${T}/xx.flexrv_vs_skato.v4.tsv.gz" IN_CELLS_OLD="${T}/xx.skato_cells.v4.tsv.gz" OUT_DIR="${T}/out" \
  "${PY}" "${HERE}/aou_v6_vs_v4.py" > "${T}/run.log" 2>&1 || { cat "${T}/run.log"; echo "SELFTEST FAILED: script exited non-zero"; exit 1; }
fail=0
chk(){ if grep -qE "$1" "${T}/out/xx.v6_vs_v4.txt"; then echo "ok   $2"; else echo "FAIL $2"; grep -nE "$3" "${T}/out/xx.v6_vs_v4.txt" || true; fail=1; fi; }
chk "^1\. synonymous below 0\.05: old 0\.0[3-6][0-9]* \(n=800\), new 0\.0[3-6]" "synonymous near 0.05 in both runs" "synonymous"
chk "^2\. pLoF old vs new .*n=800 shared, Spearman (0\.9[89]|1\.0)[0-9]*, moved >1 decade: 3 " "pLoF control: 800 shared, Spearman ~1, exactly the 3 planted moves" "pLoF"
chk "^2\. pLoF .*hits old 3 new 3 both 3 only-old 0 only-new 0; distinct genes old 3 new 3" "pLoF hit counts: the 3 moved cells are hits on both sides (1e-6 -> 1e-9), nothing else" "pLoF"
chk "^3\. splice label: .*hits 0 vs 5 \(both 0, only-first 0, only-second 5\); distinct genes 0 vs 5" "splice contrast: +5 hits, all new" "splice"
chk "^4b\. reannotation, missense only: .*Spearman 1\.000, moved 0, hits 0 vs 0" "missense-only contrast identical by construction" "4b"
chk "^5\. AlphaMissense: .*hits 0 vs 2 " "am_moderate: the 2 planted hits" "AlphaMissense"
chk "^6\. FlexRV pooled .*hits 5 vs 4 \(both 1, only-first 4, only-second 3\)" "FlexRV vs splice mask: 5 vs 4 hits, 1 shared" "FlexRV"
chk "^== floor: 1 cell\(s\)" "floor table read" "floor"
[[ -s ${T}/out/xx.controls.png && -s ${T}/out/xx.contrasts.png ]] && echo "ok   two PNGs" || { echo "FAIL PNGs"; fail=1; }
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $h["trait"]=="QuantX" && $h["hits_new_spl"]==3 && $h["hits_new_flexrv"]==2' "${T}/out/xx.v6_vs_v4.tsv" | grep -q . \
  && echo "ok   per-trait table: QuantX 3 splice hits, 2 FlexRV hits" || { echo "FAIL per-trait table"; cat "${T}/out/xx.v6_vs_v4.tsv"; fail=1; }
grep -q "^done$" "${T}/run.log" && echo "ok   done" || { echo "FAIL no done"; fail=1; }
(( fail )) && { echo "SELFTEST FAILED"; exit 1; }
echo "selftest OK: 11 checks"
