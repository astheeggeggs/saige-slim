#!/bin/bash
# Per-stage saddlepoint convergence check for the AoU FlexRV reruns (LEDGER
# #169 rollout, harness/tools/aou_spa_rerun.md section 6). Reads every
# <pre>.spaFallbacks.txt in a directory -- the COMPLETE record of saddlepoint
# failures, since a region-path cell with p NA is not written to the results
# file at all (LEDGER #164) -- and prints one line per trait:
#   <file>: fallbacks F (exact p E, NA X) | reasons: <spa: reason> xN ...
# then every NA row's cell and Note, then a verdict. A fallback that ran the
# exact test is a valid p and is counted, not flagged; an NA is flagged.
# Exit 0 on STAGE OK (no NA cell), 1 otherwise. Prints, never writes.
#   usage: bash aou_spa_stagecheck.sh <dir-with-sidecars>
set -euo pipefail
dir="${1:?usage: aou_spa_stagecheck.sh <dir-with-sidecars>}"
shopt -s nullglob
# the AoU wrapper delivers the sidecars gzipped; both forms are read
files=("${dir}"/*.spaFallbacks.txt "${dir}"/*.spaFallbacks.txt.gz)
if (( ${#files[@]} == 0 )); then echo "no *.spaFallbacks.txt[.gz] under ${dir}" >&2; exit 2; fi
work=$(mktemp -d); trap 'rm -rf "${work}"' EXIT
total=0; total_na=0; na_notes=()
for src in "${files[@]}"; do
  base=$(basename "${src}"); base=${base%.gz}; base=${base%.spaFallbacks.txt}
  if [[ ${src} == *.gz ]]; then f="${work}/${base}.spaFallbacks.txt"; gzip -dc "${src}" > "${f}"; else f="${src}"; fi
  # a placeholder (empty table: header only, or not this sidecar's header) is zero fallbacks
  if [[ $(wc -l < "${f}") -lt 2 ]] || ! head -1 "${f}" | grep -q "p_exact_reported"; then
    echo "${base}: fallbacks 0"; continue
  fi
  # columns by header name, not position
  read -r -a hdr < <(head -1 "${f}" | tr '\t' ' ')
  ip=0; inote=0; ireg=0; igrp=0; imaf=0; iw=0
  for i in "${!hdr[@]}"; do
    case "${hdr[$i]}" in
      p_exact_reported) ip=$((i+1));; Note) inote=$((i+1));; Region) ireg=$((i+1));;
      Group) igrp=$((i+1));; max_MAF) imaf=$((i+1));; Weight) iw=$((i+1));;
    esac
  done
  if (( ip == 0 || inote == 0 )); then echo "ERROR ${base}: header lacks p_exact_reported/Note" >&2; exit 2; fi
  summary=$(awk -F'\t' -v ip="${ip}" -v inote="${inote}" 'NR > 1 {
      n++; if ($ip == "NA") na++
      r = $inote; sub(/^spa: /, "", r); sub(/ \[.*$/, "", r); sub(/;.*$/, "", r); reasons[r]++
    } END {
      printf "fallbacks %d (exact p %d, NA %d)", n, n - na, na
      if (n > 0) { printf " | reasons:"; for (r in reasons) printf " <%s> x%d", r, reasons[r] }
      printf "\n" }' "${f}")
  echo "${base}: ${summary}"
  n=$(awk -F'\t' 'NR > 1 {n++} END {print n+0}' "${f}")
  na=$(awk -F'\t' -v ip="${ip}" 'NR > 1 && $ip == "NA" {na++} END {print na+0}' "${f}")
  total=$((total + n)); total_na=$((total_na + na))
  if (( na > 0 )); then
    while IFS= read -r line; do na_notes+=("${base}: ${line}"); done < <(
      awk -F'\t' -v ip="${ip}" -v inote="${inote}" -v ireg="${ireg}" -v igrp="${igrp}" -v imaf="${imaf}" -v iw="${iw}" \
        'NR > 1 && $ip == "NA" { printf "%s %s %s %s | %s\n", $ireg, $igrp, $imaf, (iw ? $iw : "-"), $inote }' "${f}")
  fi
done
if (( ${#na_notes[@]} > 0 )); then
  echo "--- NA cells (no usable p; the Note says why the saddlepoint failed and why the exact test could not run):"
  printf '%s\n' "${na_notes[@]}"
fi
echo "stage: ${#files[@]} sidecar(s), fallbacks ${total}, NA ${total_na}"
if (( total_na == 0 )); then echo "STAGE OK"; exit 0; fi
echo "STAGE: ${total_na} NA cell(s) -- read their Notes; a 'saddle(q) invalid' with qsaddle=INVALID on 1.5.2-dev-4ff9576b or later is a NEW mechanism, stop the rollout"
exit 1
