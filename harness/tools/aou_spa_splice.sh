#!/bin/bash
# Splice a filtered FlexRV rerun (LEDGER #164) over the original delivered files.
# Runs INSIDE a dsub job (gcr.io/google.com/cloudsdktool/google-cloud-cli:slim);
# streams with gcloud storage cat. Env:
#   OLD_PREFIX  gs://.../<pre>          (v4: .txt.gz, .singleAssoc.txt.gz, .pooledTests.txt.gz)
#   NEW_PREFIX  gs://.../<pre>          (v5 rerun of the listed regions only)
#   REGIONS     gs://.../rerun/<pre>.regions.txt
#   FAILURES    gs://.../failures.tsv   (the scan's; rows of this file's failures are used)
#   OUT_DIR     local dir dsub delocalizes (--output-recursive)
set -euo pipefail
: "${OLD_PREFIX:?}"; : "${NEW_PREFIX:?}"; : "${REGIONS:?}"; : "${FAILURES:?}"; : "${OUT_DIR:?}"
mkdir -p "${OUT_DIR}"; pre=$(basename "${OLD_PREFIX}")
cat_any(){ if [[ $1 == *.gz ]]; then gcloud storage cat "$1" | gzip -dc; else gcloud storage cat "$1"; fi; }
gcloud storage cat "${REGIONS}" | sed '/^$/d' | sort -u > regions.txt
cat_any "${FAILURES}" | awk -F'\t' -v f="${pre}" 'NR>1 && $1==f && $3!="-" {print $3}' | sort -u > failregions.txt
echo "regions to splice: $(wc -l < regions.txt); with failed rows: $(wc -l < failregions.txt)"

# region of a row: results/pooled = column Region; singleAssoc = MarkerID prefix (UR rows)
splice_one(){  # suffix  region-column-rule (awk expr on $0 giving the region)
  local sfx="$1" rule="$2" old="${OLD_PREFIX}${sfx}" new="${NEW_PREFIX}${sfx}"
  local out="${OUT_DIR}/${pre}${sfx%.gz}"
  cat_any "${old}" > old.tsv; cat_any "${new}" > new.tsv
  head -1 old.tsv > hdr.old; head -1 new.tsv > hdr.new
  cmp -s hdr.old hdr.new || { echo "ERROR ${sfx}: header differs between old and new" >&2; exit 1; }
  # new rows, keyed by region (only the listed regions are expected in new)
  awk -F'\t' -v rule="${rule}" 'NR==FNR{keep[$1]=1; next} FNR==1{next} {reg=('"${rule}"'); if(!(reg in keep)){bad++} } END{if(bad){print "ERROR: " bad " new rows outside the region list" > "/dev/stderr"; exit 1}}' regions.txt new.tsv
  # CONTROL 1: regions WITHOUT failed rows must be byte-identical old vs new
  awk -F'\t' 'NR==FNR{f[$1]=1; next} !($1 in f)' failregions.txt regions.txt > ctrl.txt
  for r in $(cat ctrl.txt); do
    awk -F'\t' -v R="$r" 'FNR>1 && ('"${rule}"')==R' old.tsv > o.r; awk -F'\t' -v R="$r" 'FNR>1 && ('"${rule}"')==R' new.tsv > n.r
    cmp -s o.r n.r || { echo "ERROR ${sfx}: control region ${r} (no failed row) CHANGED -- splice refused" >&2; exit 1; }
  done
  # CONTROL 2: every region with failed rows must change somewhere
  for r in $(cat failregions.txt); do
    awk -F'\t' -v R="$r" 'FNR>1 && ('"${rule}"')==R' old.tsv > o.r; awk -F'\t' -v R="$r" 'FNR>1 && ('"${rule}"')==R' new.tsv > n.r
    if cmp -s o.r n.r && [[ ${sfx} == ".txt.gz" ]]; then echo "WARNING ${sfx}: failed region ${r} came back identical" >&2; fi
  done
  # the splice: old rows not in the list, then new rows, in old-file order of first appearance
  { head -1 old.tsv
    awk -F'\t' -v rule="${rule}" 'NR==FNR{keep[$1]=1; next} FNR==1{next} {reg=('"${rule}"'); if(!(reg in keep)) print}' regions.txt old.tsv
    awk -F'\t' 'FNR>1' new.tsv
  } > "${out}"
  gzip -f "${out}"
  echo "${sfx}: old $(($(wc -l < old.tsv)-1)) rows, new $(($(wc -l < new.tsv)-1)) rows, spliced $(( $(gzip -dc "${out}.gz" | wc -l) - 1 )) rows"
}
splice_one ".txt.gz"             '$1'
splice_one ".singleAssoc.txt.gz" 'substr($3, 1, index($3, ":") ? index($3, ":")-1 : length($3))'
if gcloud storage ls "${OLD_PREFIX}.pooledTests.txt.gz" >/dev/null 2>&1; then
  splice_one ".pooledTests.txt.gz" '$1'
fi
{ echo "spliced $(date -Is)"; echo "old: ${OLD_PREFIX}"; echo "new: ${NEW_PREFIX}"
  echo "regions replaced: $(wc -l < regions.txt)"; echo "regions with failed saddlepoints: $(wc -l < failregions.txt)"
  echo "--- regions:"; cat regions.txt; } > "${OUT_DIR}/${pre}.splice.txt"
echo done
