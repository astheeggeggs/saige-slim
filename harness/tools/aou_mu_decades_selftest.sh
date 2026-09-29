#!/bin/bash
# Self-test of aou_mu_decades.sh on the bundled PINNED fixture (simulated data,
# so its counts may be printed). Controls: the fitted model (a separated group
# of 60 straddling 1e-8: 13 below, 47 above), the planted model (all 60 at
# 1e-12), a model with no pinned sample, a quantitative model (skipped), and
# the masking itself (the same table at MASK=20 shows no count from 1 to 19).
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
RS="${RSCRIPT:-${REPO}/.pixi/envs/default/bin/Rscript}"
W="$(mktemp -d)"; trap 'rm -rf "${W}"' EXIT
mkdir -p "${W}/models/a" "${W}/models/b"
G="${REPO}/harness/goldens/step1"
cp "${G}/binary_sparseGRM_pinned.rda" "${W}/models/a/"
cp "${G}/binary_sparseGRM_pinned_planted.rda" "${W}/models/a/"
cp "${G}/binary_sparseGRM_a.rda" "${W}/models/b/"
cp "${G}/quant_sparseGRM.rda" "${W}/models/b/"
cp "${G}/binary_sparseGRM_a.rda" "${W}/models/b/binary_sparseGRM_a_chr1.rda"   # a per-chromosome file: skipped
fails=0; checks=0
check(){ checks=$((checks+1)); if [[ "$2" == "$3" ]]; then echo "  ok   $1"; else echo "  FAIL $1: got '$2', want '$3'"; fails=$((fails+1)); fi; }
cell(){ # cell <table> <model> <column>
  awk -F'\t' -v M="$2" -v C="$3" 'NR==1{for(i=1;i<=NF;i++) if($i==C) c=i; next} $1==M{print $c}' "$1"; }
cd "${W}"
MODELS="${W}/models" OUT_DIR="${W}/o1" TAG=t MASK=1 RSCRIPT="${RS}" bash "${REPO}/harness/tools/aou_mu_decades.sh" > "${W}/log1" 2>&1 || { cat "${W}/log1"; exit 1; }
gunzip -c "${W}/o1/mu_decades.t.tsv.gz" > "${W}/t1.tsv"
check "three binary models read (quantitative and per-chromosome files skipped)" "$(awk 'END{print NR-1}' "${W}/t1.tsv")" "3"
check "fitted: at or below 1e-8" "$(cell "${W}/t1.tsv" binary_sparseGRM_pinned.rda n_le_1e-08)" "13"
check "fitted: at or below 1e-6" "$(cell "${W}/t1.tsv" binary_sparseGRM_pinned.rda n_le_1e-06)" "60"
check "fitted: decade 1e-8..1e-7" "$(cell "${W}/t1.tsv" binary_sparseGRM_pinned.rda 1e-8_1e-7)" "45"
check "fitted: a gap of empty decades above the group" "$(cell "${W}/t1.tsv" binary_sparseGRM_pinned.rda empty_decades_in_widest_gap)" "3"
check "planted: all 60 in 1e-15..1e-12" "$(cell "${W}/t1.tsv" binary_sparseGRM_pinned_planted.rda 1e-15_1e-12)" "60"
check "planted: mass at or below 1e-8" "$(cell "${W}/t1.tsv" binary_sparseGRM_pinned_planted.rda mass_le_1e-08)" "6.00e-11"
check "no pinned sample: at or below 1e-5" "$(cell "${W}/t1.tsv" binary_sparseGRM_a.rda n_le_1e-05)" "0"
check "no pinned sample: no gap" "$(cell "${W}/t1.tsv" binary_sparseGRM_a.rda empty_decades_in_widest_gap)" "0"
check "log: the straddle is counted" "$(grep -c 'not all at or below 1e-8: 1' "${W}/log1")" "1"
MODELS="${W}/models" OUT_DIR="${W}/o2" TAG=t MASK=20 RSCRIPT="${RS}" bash "${REPO}/harness/tools/aou_mu_decades.sh" > "${W}/log2" 2>&1
gunzip -c "${W}/o2/mu_decades.t.tsv.gz" > "${W}/t2.tsv"
check "masked: 13 is written <20" "$(cell "${W}/t2.tsv" binary_sparseGRM_pinned.rda n_le_1e-08)" "<20"
check "masked: its mass is withheld" "$(cell "${W}/t2.tsv" binary_sparseGRM_pinned.rda mass_le_1e-08)" "masked"
check "masked: 60 is printed" "$(cell "${W}/t2.tsv" binary_sparseGRM_pinned.rda n_le_1e-06)" "60"
# no bare count from 1 to 19 anywhere in the masked table's count columns (N and cases aside)
bad=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i ~ /^(n_le_|le1e|1e-|gt1e|contradicting)/) c[i]=1; next} {for(i in c) if($i ~ /^[0-9]+$/ && $i+0 > 0 && $i+0 < 20) n++} END{print n+0}' "${W}/t2.tsv")
check "masked: no count from 1 to 19 survives" "${bad}" "0"
# the log prints the same table: the same check on its tab-separated lines (the
# widest-gap column is a number of decades, not of people, and is not a count)
badlog=$(grep "$(printf '\t')" "${W}/log2" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i ~ /^(n_le_|le1e|1e-|gt1e|contradicting)/) c[i]=1; next} {for(i in c) if($i ~ /^[0-9]+$/ && $i+0 > 0 && $i+0 < 20) n++} END{print n+0}')
check "masked: the log prints none either" "${badlog}" "0"
# positive control of the check itself: the UNMASKED table trips it
tripped=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i ~ /^(n_le_|le1e|1e-|gt1e|contradicting)/) c[i]=1; next} {for(i in c) if($i ~ /^[0-9]+$/ && $i+0 > 0 && $i+0 < 20) n++} END{print (n+0 > 0) ? "yes" : "no"}' "${W}/t1.tsv")
check "the masking check detects an unmasked table" "${tripped}" "yes"
echo "${checks} checks, ${fails} failed"
exit $(( fails > 0 ))
