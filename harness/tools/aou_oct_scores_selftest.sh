#!/bin/bash
# Self-test for aou_oct_scores.py on FABRICATED summary tables (no AoU data). Planted:
# AMPCT and SCALLION share one uniform draw per (trait, gene) (Spearman near 1) except, in
# QuantX only, AMPCT hits at four genes and SCALLION at three with two shared (so hits 4 vs
# 3, both 2, only-AMPCT 2, only-SCALLION 1) and one gene AMPCT moves to 1e-3 from a fixed
# 0.5 -- exactly four cells moved by more than a decade; RAND uniform on both traits with ONE
# gene of Asth planted at 1e-10, so Asth is leaky (hits_rand 1) and the clean pooled null is
# QuantX alone (400 genes, no hit); synonymous SKAT-O burden uniform. A second run without the
# optional tables must say so.
# Needs a python with pandas, matplotlib and scipy (the pyplot image, or a venv); pass it
# as $1, default python3.
set -euo pipefail
PY="${1:-python3}"; HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
"${PY}" - "${T}" <<'PYEOF'
import sys, gzip, numpy as np
T = sys.argv[1]; rng = np.random.default_rng(11); G = 400; traits = ["QuantX", "Asth"]
base = {(tr, g): rng.uniform(1e-4, 1) for tr in traits for g in range(G)}
base[("QuantX", 50)] = 0.5
AM = {0, 1, 2, 3}; SC = {0, 1, 5}
def p_am(tr, g): return 1e-8 if (tr == "QuantX" and g in AM) else (1e-3 if (tr == "QuantX" and g == 50) else base[(tr, g)])
def p_sc(tr, g): return 1e-8 if (tr == "QuantX" and g in SC) else base[(tr, g)]
def write_sum(path, tag, pf, trs):
    with gzip.open(path, "wt") as fh:
        fh.write("anc\ttrait\tgene\tsource\tp_flexrv\tp_flexrv_best\tset_best\tn_sets\tp_skato\tp_skato_best\tcell_best\tp_burden_best\tp_skat_best\n")
        for tr in trs:
            for g in range(G):
                fh.write(f"xx\t{tr}\tG{g}\t{tag}\t{pf(tr, g):.6g}\tNA\tNA\t112\tNA\tNA\tNA\tNA\tNA\n")
write_sum(f"{T}/xx.flexrv_vs_skato.oct_fxa.tsv.gz", "oct_fxa", p_am, traits)
write_sum(f"{T}/xx.flexrv_vs_skato.oct_fxs.tsv.gz", "oct_fxs", p_sc, traits)
write_sum(f"{T}/xx.flexrv_vs_skato.oct_rand.tsv.gz", "oct_rand", lambda tr, g: 1e-10 if (tr == "Asth" and g == 7) else rng.uniform(1e-4, 1), traits)
masks = ["pLoF", "synonymous"]; mafs = ["0.0001", "0.001", "0.01"]
cols = [f"burden:{m}@{f}" for m in masks for f in mafs] + [f"skato:{m}@{f}" for m in masks for f in mafs]
with gzip.open(f"{T}/xx.skato_cells.oct_fxa.tsv.gz", "wt") as fh:
    fh.write("anc\ttrait\tgene\t" + "\t".join(cols) + "\n")
    for tr in traits:
        for g in range(G):
            fh.write(f"xx\t{tr}\tG{g}\t" + "\t".join(f"{rng.uniform():.6g}" for _ in cols) + "\n")
PYEOF
ANC=xx IN_AMPCT="${T}/xx.flexrv_vs_skato.oct_fxa.tsv.gz" IN_SCALLION="${T}/xx.flexrv_vs_skato.oct_fxs.tsv.gz" \
  IN_RAND="${T}/xx.flexrv_vs_skato.oct_rand.tsv.gz" IN_CELLS="${T}/xx.skato_cells.oct_fxa.tsv.gz" OUT_DIR="${T}/out" DPI=60 \
  "${PY}" "${HERE}/aou_oct_scores.py" > "${T}/run.log" 2>&1 || { cat "${T}/run.log"; echo "SELFTEST FAILED: script exited non-zero"; exit 1; }
fail=0
R="${T}/out/xx.scores.txt"
chk(){ if grep -qE "$1" "${R}"; then echo "ok   $2"; else echo "FAIL $2"; grep -nE "$3" "${R}" || true; fail=1; fi; }
chk "^xx: 2 traits with both scores \(AMPCT 2, SCALLION 2, RAND 2\)" "header: 2 traits, RAND on 2" "^xx"
chk "^1\. RAND traits with a gene below 1e-08 .*: Asth \(min p 1\.0e-10, 1 gene\)$" "leaky: Asth named with its min p and one gene" "^1\."
chk "^   RAND on synonymous, pooled over the 1 clean traits .*n=400; .*hits 0; below 0\.05 0\.0[2-8]" "clean RAND null: QuantX's 400 genes, no hit" "clean"
chk "^   RAND on synonymous, pooled over ALL its traits .*n=800; .*hits 1;" "all-traits RAND pool: 800 genes, the planted hit" "ALL"
chk "^2\. SKAT-O synonymous burden at 0\.001.*n=800; .*hits 0;" "synonymous SKAT-O null: 800 cells, no hit" "^2\."
chk "^3\. AMPCT vs SCALLION: n=800 shared, Spearman (0\.9[5-9]|1\.0)[0-9]*, moved >1 decade: 4, hits AMPCT 4 SCALLION 3 both 2 only-AMPCT 2 only-SCALLION 1; distinct genes AMPCT 4 SCALLION 3" \
    "contrast: Spearman ~1, the 4 planted moves, hits 4 vs 3 with 2 shared" "^3\."
chk "^4\. AMPCT pooled over traits: n=800; .*hits 4;" "AMPCT pooled tail: 4 hits" "^4\."
chk "^5\. SCALLION pooled over traits: n=800; .*hits 3;" "SCALLION pooled tail: 3 hits" "^5\."
[[ -s ${T}/out/xx.scores.png ]] && echo "ok   PNG" || { echo "FAIL PNG"; fail=1; }
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $h["trait"]=="QuantX" && $h["hits_ampct"]==4 && $h["hits_scallion"]==3 && $h["hits_both"]==2 && $h["hits_rand"]==0' "${T}/out/xx.scores.tsv" | grep -q . \
  && echo "ok   per-trait table: QuantX 4 / 3 / 2 shared / RAND 0" || { echo "FAIL per-trait table"; cat "${T}/out/xx.scores.tsv"; fail=1; }
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $h["trait"]=="Asth" && $h["hits_ampct"]==0 && $h["hits_rand"]==1 && $h["rand_leaky"]=="True" && $h["rand_min_p"]+0==1e-10' "${T}/out/xx.scores.tsv" | grep -q . \
  && echo "ok   per-trait table: Asth no score hit, RAND 1 hit, leaky, min p 1e-10" || { echo "FAIL Asth row"; cat "${T}/out/xx.scores.tsv"; fail=1; }
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR>1 && $h["trait"]=="QuantX" && $h["rand_leaky"]=="False"' "${T}/out/xx.scores.tsv" | grep -q . \
  && echo "ok   per-trait table: QuantX not leaky" || { echo "FAIL QuantX leaky flag"; cat "${T}/out/xx.scores.tsv"; fail=1; }
grep -q "^done$" "${T}/run.log" && echo "ok   done" || { echo "FAIL no done"; fail=1; }
# ---- without the optional tables: the controls say they are not shown, the contrast is unchanged
ANC=xx IN_AMPCT="${T}/xx.flexrv_vs_skato.oct_fxa.tsv.gz" IN_SCALLION="${T}/xx.flexrv_vs_skato.oct_fxs.tsv.gz" OUT_DIR="${T}/out2" DPI=60 \
  "${PY}" "${HERE}/aou_oct_scores.py" > "${T}/run2.log" 2>&1 || { cat "${T}/run2.log"; echo "SELFTEST FAILED: second run exited non-zero"; exit 1; }
grep -q "^1\. RAND table not given" "${T}/out2/xx.scores.txt" && grep -q "^2\. cells table not given" "${T}/out2/xx.scores.txt" \
  && echo "ok   optional tables absent: both controls say so" || { echo "FAIL absent-table lines"; fail=1; }
grep -q "^3\. AMPCT vs SCALLION: n=800 shared" "${T}/out2/xx.scores.txt" && echo "ok   contrast still computed" || { echo "FAIL contrast without tables"; fail=1; }
(( fail )) && { echo "SELFTEST FAILED"; exit 1; }
echo "selftest OK: 15 checks"
