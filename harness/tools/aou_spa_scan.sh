#!/bin/bash
# Scan delivered SAIGE step-2 outputs for saddlepoint fallbacks (LEDGER #163/#164)
# and list what a filtered rerun must cover. Runs INSIDE a dsub job on an image
# that has gcloud (google/cloud-sdk:slim); streams every file with `gcloud
# storage cat`, so nothing is localized and nothing leaves the region.
#
# Inputs (env): SCAN_DIR   gs:// directory holding *.singleAssoc.txt.gz and the
#                          matching results files (<pre>.txt.gz)
#               OUT_DIR    local dir dsub delocalizes (--output-recursive)
#               TAIL_P     burden-cell tail threshold, default 1e-5
#               PATTERN    optional glob on basenames, default '*'
#
# Signature of a failed saddlepoint in the sidecar: Is.SPA == FALSE and
# p.value.NA < 0.0455 (|Z| > 2, the default --SPAcutoff: attempted, not
# converged; a row that took the exact test reads Is.SPA = TRUE). Reported p on
# such a row is the normal-approximation p.
#
# Outputs in OUT_DIR:
#   summary.tsv          one line per sidecar: rows, UR rows, failures, regions
#   failures.tsv         every failing row (file, MarkerID, Region, CHR, p, p.NA)
#   tail_regions.tsv     regions with any burden cell below TAIL_P (results file)
#   rerun/<pre>.regions.txt   union: failure regions + tail regions (group arms)
#   rerun/<pre>.markers.txt   real-marker failures (variant arm, --idstoIncludeFile)
set -euo pipefail
: "${SCAN_DIR:?}"; : "${OUT_DIR:?}"
TAIL_P="${TAIL_P:-1e-5}"; PATTERN="${PATTERN:-*}"
mkdir -p "${OUT_DIR}/rerun"
SUM="${OUT_DIR}/summary.tsv"; FAIL="${OUT_DIR}/failures.tsv"; TAIL="${OUT_DIR}/tail_regions.tsv"
printf 'file\trows\tur_rows\tfail_rows\tfail_ur_rows\tfail_marker_rows\tfail_regions\ttail_regions\trerun_regions\n' > "${SUM}"
printf 'file\tMarkerID\tRegion\tCHR\tp.value\tp.value.NA\tIs.SPA\n' > "${FAIL}"
printf 'file\tRegion\tmin_Pvalue_Burden\n' > "${TAIL}"

cat_any() {  # stream a gs:// file, gunzipping if needed
  if [[ $1 == *.gz ]]; then gcloud storage cat "$1" | gzip -dc; else gcloud storage cat "$1"; fi
}

mapfile -t SIDES < <(gcloud storage ls "${SCAN_DIR}/${PATTERN}.singleAssoc.txt.gz" 2>/dev/null || true)
echo "sidecars: ${#SIDES[@]}"
for side in "${SIDES[@]}"; do
  base=$(basename "${side}")
  pre="${base%.singleAssoc.txt.gz}"
  res="${SCAN_DIR}/${pre}.txt.gz"
  # --- the sidecar: failing rows ---
  cat_any "${side}" | awk -F'\t' -v f="${pre}" -v OFS='\t' '
    NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
    { rows++; ur = ($h["CHR"]=="UR"); if(ur) urrows++
      if (tolower($h["Is.SPA"])=="false" && ($h["p.value.NA"]+0) < 0.0455) {
        id=$h["MarkerID"]; reg = ur ? substr(id, 1, index(id, ":")-1) : "-"
        fails++; if(ur) failur++; else failmk++
        if(ur) regs[reg]=1; else mks[id]=1
        print f, id, reg, $h["CHR"], $h["p.value"], $h["p.value.NA"], $h["Is.SPA"] >> "'"${FAIL}"'"
      } }
    END { nr=0; for(r in regs){ nr++; print r > "'"${OUT_DIR}/rerun/${pre}.failregions.tmp"'" }
          for(m in mks) print m > "'"${OUT_DIR}/rerun/${pre}.markers.txt"'"
          printf "%s\t%d\t%d\t%d\t%d\t%d\t%d", f, rows, urrows, fails, failur, failmk, nr > "'"${OUT_DIR}/rerun/${pre}.sum.tmp"'" }'
  # --- the results file: tail regions (group arms only; a variant file has no Region column) ---
  ntail=0
  if gcloud storage ls "${res}" >/dev/null 2>&1; then
    cat_any "${res}" | awk -F'\t' -v f="${pre}" -v tp="${TAIL_P}" -v OFS='\t' '
      NR==1 { for(i=1;i<=NF;i++) h[$i]=i; if(!("Region" in h) || !("Pvalue_Burden" in h)) exit; next }
      $h["Group"]!="Cauchy" && $h["Pvalue_Burden"]!="NA" && ($h["Pvalue_Burden"]+0) < tp {
        r=$h["Region"]; p=$h["Pvalue_Burden"]+0; if(!(r in mn) || p < mn[r]) mn[r]=p }
      END { for(r in mn){ print f, r, mn[r] >> "'"${TAIL}"'"; print r > "'"${OUT_DIR}/rerun/${pre}.tailregions.tmp"'" } }'
    [[ -s "${OUT_DIR}/rerun/${pre}.tailregions.tmp" ]] && ntail=$(wc -l < "${OUT_DIR}/rerun/${pre}.tailregions.tmp")
  fi
  cat "${OUT_DIR}/rerun/${pre}.failregions.tmp" "${OUT_DIR}/rerun/${pre}.tailregions.tmp" 2>/dev/null \
    | sort -u > "${OUT_DIR}/rerun/${pre}.regions.txt" || true
  nrerun=$(wc -l < "${OUT_DIR}/rerun/${pre}.regions.txt")
  printf '%s\t%d\t%d\n' "$(cat "${OUT_DIR}/rerun/${pre}.sum.tmp")" "${ntail}" "${nrerun}" >> "${SUM}"
  rm -f "${OUT_DIR}/rerun/${pre}".*.tmp
  [[ -s "${OUT_DIR}/rerun/${pre}.regions.txt" ]] || rm -f "${OUT_DIR}/rerun/${pre}.regions.txt"
  echo "$(tail -1 "${SUM}")"
done
echo "done: $(($(wc -l < "${SUM}") - 1)) files, $(($(wc -l < "${FAIL}") - 1)) failing rows"
