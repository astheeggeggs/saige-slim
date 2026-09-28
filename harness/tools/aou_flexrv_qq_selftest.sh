#!/bin/bash
# Self-test for aou_flexrv_qq.py on FABRICATED tables (no AoU data): two traits,
# 18 cells, planted hits; expects one PNG per trait, binary trait first, and the
# 'done' line. Needs a python with pandas, matplotlib and scipy (the pyplot
# image, or a venv); pass it as $1, default python3.
set -euo pipefail
PY="${1:-python3}"; HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
"${PY}" - "${T}" <<'PYEOF'
import sys, gzip, numpy as np
T = sys.argv[1]; rng = np.random.default_rng(1); G = 400
masks = ["pLoF","damaging_missense","other_missense","synonymous","pLoF;damaging_missense","pLoF;damaging_missense;other_missense;synonymous"]
mafs = ["0.0001","0.001","0.01"]; cols = [f"burden:{m}@{f}" for m in masks for f in mafs] + [f"skato:{m}@{f}" for m in masks for f in mafs]
with gzip.open(f"{T}/xx.skato_cells.v4.tsv.gz", "wt") as fh:
    fh.write("anc\ttrait\tgene\t" + "\t".join(cols) + "\n")
    for tr in ["QuantX", "Asth"]:
        for g in range(G):
            vals = [f"{(10**-rng.uniform(4,8) if (g < 3 and c.startswith('burden:pLoF@')) else rng.uniform()):.6g}" for c in cols]
            fh.write(f"xx\t{tr}\tG{g}\t" + "\t".join(vals) + "\n")
with gzip.open(f"{T}/xx.flexrv_vs_skato.v4.tsv.gz", "wt") as fh:
    fh.write("anc\ttrait\tgene\tsource\tp_flexrv\tp_flexrv_best\tset_best\tn_sets\tp_skato\tp_skato_best\tcell_best\tp_burden_best\tp_skat_best\n")
    for tr in ["QuantX", "Asth"]:
        for g in range(G):
            fh.write(f"xx\t{tr}\tG{g}\tv4\t{(1e-7 if g < 2 else rng.uniform()):.6g}\tNA\tNA\t96\tNA\tNA\tNA\tNA\tNA\n")
PYEOF
ANC=xx IN_SUMMARY="${T}/xx.flexrv_vs_skato.v4.tsv.gz" IN_CELLS="${T}/xx.skato_cells.v4.tsv.gz" OUT_DIR="${T}/out" DPI=60 "${PY}" "${HERE}/aou_flexrv_qq.py" > "${T}/run.log" 2>&1 || { cat "${T}/run.log"; echo "SELFTEST FAILED: script exited non-zero"; exit 1; }
fail=0
[[ -s ${T}/out/xx_Asth.png && -s ${T}/out/xx_QuantX.png ]] && echo "ok   one PNG per trait" || { echo "FAIL PNGs"; fail=1; }
grep -q "xx: 2 traits with both tables; cells columns: 18" "${T}/run.log" && echo "ok   18 burden cells recognised" || { echo "FAIL header line"; cat "${T}/run.log"; fail=1; }
grep -q "^done$" "${T}/run.log" && echo "ok   done" || { echo "FAIL no done"; fail=1; }
(( fail )) && { echo "SELFTEST FAILED"; exit 1; }
echo "selftest OK: 3 checks"
