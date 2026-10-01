#!/bin/bash
# Self-test for aou_flexrv_calibration.py on FABRICATED summaries (no AoU data).
# Controls before conclusions: a NULL arm that must stay quiet, a PLANTED arm
# that must fire, and the arithmetic checked against a closed form.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
cd "${T}"
python3 - <<'PY'
import gzip, random
random.seed(20261001)
H = ["anc","trait","gene","source","p_flexrv","p_flexrv_best","set_best","n_sets",
     "p_skato","p_skato_best","cell_best","p_burden_best","p_skat_best"]
M = 18000
def row(anc, trait, i, pf, ps):
    return [anc, trait, f"ENSG{i:08d}", "v7", f"{pf:.6g}", f"{pf:.6g}",
            "AM_unw_unw", "112", f"{ps:.6g}", f"{ps:.6g}", "pLoF@0.001",
            f"{ps:.6g}", f"{ps:.6g}"]
rows = []
# 20 NULL traits: both arms uniform. These are the negative control.
for t in range(20):
    for i in range(M):
        rows.append(row("tst", f"null{t:02d}", i, random.random(), random.random()))
# PLANTED A: one absurd FlexRV gene, SKAT-O left null -> the #172/#176 shape
for i in range(M):
    pf = 1e-40 if i == 7 else random.random()
    rows.append(row("tst", "plantFlexRVonly", i, pf, random.random()))
# PLANTED B: both arms extreme on the same gene -> real biology, must be
# flagged on calibration but marked CONCORDANT by the two arms agreeing
for i in range(M):
    e = (i == 11)
    rows.append(row("tst", "plantBoth", i, 1e-30 if e else random.random(),
                                        1e-28 if e else random.random()))
with gzip.open("tst.flexrv_vs_skato.v4.tsv.gz", "wt") as fh:
    fh.write("\t".join(H) + "\n")
    for r in rows: fh.write("\t".join(r) + "\n")
print("fixture written:", len(rows), "gene rows")
PY
python3 "${HERE}/aou_flexrv_calibration.py" tst.flexrv_vs_skato.v4.tsv.gz > calib.tsv
fail=0
chk(){ if eval "$2"; then echo "ok   $1"; else echo "FAIL $1"; fail=1; fi; }

# --- arithmetic against the closed form --------------------------------------
python3 - <<'PY' > arith.txt
import csv, math
r = {(x["trait"], x["arm"]): x for x in csv.DictReader(open("calib.tsv"), delimiter="\t")}
x = r[("null00", "flexrv")]
m = int(x["n_genes"]); mp = float(x["min_p"])
med_expect = -math.expm1(math.log(0.5)/m)
pcal_expect = -math.expm1(m*math.log1p(-mp))
print("n_genes", m)
# the table is written to 6 significant figures, so every comparison here is
# RELATIVE at 1e-5 -- a tolerance tighter than the digits on disk tests the
# formatter, not the arithmetic (the first version of this test used 1e-12
# absolute and failed for exactly that reason)
def rel(a, b): return abs(a - b) <= 1e-5 * abs(b)
print("median_ok", rel(float(x["null_median"]), med_expect))
print("pcal_ok",   rel(float(x["p_cal"]),       pcal_expect))
print("ratio_ok",  rel(float(x["ratio"]),       mp / med_expect))
PY
chk "n_genes is every usable gene"      "grep -q '^n_genes 18000$' arith.txt"
chk "null_median matches the closed form" "grep -q '^median_ok True$' arith.txt"
chk "p_cal matches the closed form"       "grep -q '^pcal_ok True$' arith.txt"
chk "ratio is min_p / null_median"        "grep -q '^ratio_ok True$' arith.txt"

# --- NEGATIVE control: 20 null traits x 2 arms, p_cal must be ~uniform --------
nlow=$(awk -F'\t' 'NR>1 && $2 ~ /^null/ && $11+0 < 0.01' calib.tsv | wc -l | tr -d ' ')
chk "null traits: at most 1 of 40 arms below p_cal 0.01 (expect 0.4)" "[ ${nlow} -le 1 ]"
nhuge=$(awk -F'\t' 'NR>1 && $2 ~ /^null/ && $11+0 < 1e-6' calib.tsv | wc -l | tr -d ' ')
chk "null traits: NONE below p_cal 1e-6"  "[ ${nhuge} -eq 0 ]"

# --- POSITIVE control: the planted traits must fire --------------------------
chk "planted FlexRV-only fires on the flexrv arm" \
    "awk -F'\t' '\$2==\"plantFlexRVonly\" && \$3==\"flexrv\" && \$11+0 < 1e-20' calib.tsv | grep -q ."
chk "planted FlexRV-only stays QUIET on the skato arm" \
    "awk -F'\t' '\$2==\"plantFlexRVonly\" && \$3==\"skato\" && \$11+0 > 0.01' calib.tsv | grep -q ."
chk "planted FlexRV-only names the right gene" \
    "awk -F'\t' '\$2==\"plantFlexRVonly\" && \$3==\"flexrv\" && \$7==\"ENSG00000007\"' calib.tsv | grep -q ."
chk "planted BOTH fires on both arms (the biology shape)" \
    "[ \$(awk -F'\t' '\$2==\"plantBoth\" && \$11+0 < 1e-20' calib.tsv | wc -l) -eq 2 ]"

# --- the discriminator the instrument exists for -----------------------------
# one arm extreme = the defect shape; both = probably real. The table must make
# that readable without recomputation.
chk "the two shapes are distinguishable from the table alone" \
    "[ \$(awk -F'\t' 'NR>1 && \$11+0 < 1e-20 {print \$2}' calib.tsv | sort -u | wc -l) -eq 2 ]"

if (( fail )); then echo "SELFTEST FAILED"; exit 1; fi
echo "calibration selftest OK: 4 arithmetic checks, 2 negative controls (40 null arms), 4 positive controls, 1 discriminator"
