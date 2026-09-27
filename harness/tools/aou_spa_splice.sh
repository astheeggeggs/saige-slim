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
cat_any(){  # gs:// or a local path, gzipped or not
  if [[ $1 == gs://* ]]; then src=(gcloud storage cat "$1"); else src=(cat "$1"); fi
  if [[ $1 == *.gz ]]; then "${src[@]}" | gzip -dc; else "${src[@]}"; fi; }
exists(){ if [[ $1 == gs://* ]]; then gcloud storage ls "$1" >/dev/null 2>&1; else [[ -s $1 ]]; fi; }
cat_any "${REGIONS}" | sed '/^$/d' | sort -u > regions.txt
cat_any "${FAILURES}" | awk -F'\t' -v f="${pre}" 'NR>1 && $1==f && $3!="-" {print $3}' | sort -u > failregions.txt
echo "regions to splice: $(wc -l < regions.txt); with failed rows: $(wc -l < failregions.txt)"

# region of a row: results/pooled = column Region; singleAssoc = MarkerID prefix (UR rows)
splice_one(){  # suffix  region-column-rule (awk expr on $0 giving the region)
  local sfx="$1" rule="$2"
  local old="${OLD_PREFIX}${sfx}" new="${NEW_PREFIX}${sfx}"   # a separate `local`: bash expands a
  local out="${OUT_DIR}/${pre}${sfx%.gz}"                     # local's words before it assigns them
  cat_any "${old}" > old.tsv; cat_any "${new}" > new.tsv
  # A sidecar the run had nothing to write is delivered as a one-line '#'
  # placeholder (the wrapper's doing: dsub fails delocalization on a declared
  # output that does not exist). Treat a placeholder as an EMPTY table with the
  # other side's header; both placeholders means nothing to splice.
  local oph=0 nph=0
  [[ $(head -c 1 old.tsv) == "#" ]] && oph=1
  [[ $(head -c 1 new.tsv) == "#" ]] && nph=1
  if (( oph && nph )); then
    cp old.tsv "${out}"; gzip -f "${out}"; echo "${sfx}: placeholder on both sides, carried over"; return 0
  elif (( oph )); then head -1 new.tsv > old.tsv; echo "${sfx}: old is a placeholder (no region had anything to disclose); treating as empty"
  elif (( nph )); then head -1 old.tsv > new.tsv; echo "${sfx}: new is a placeholder; treating as empty"
  fi
  head -1 old.tsv > hdr.old; head -1 new.tsv > hdr.new
  cmp -s hdr.old hdr.new || { echo "ERROR ${sfx}: header differs between old and new" >&2; exit 1; }
  # The singleAssoc holds TWO kinds of row (found on the eas/afr stage,
  # 2026-09-27). A pseudo-marker row's MarkerID is Region:Group:maxMAF:Weight
  # and names its region. A REAL marker -- rare but not ultra-rare, MAC above
  # the collapse cutoff, which exists whenever 2*N*maxMAF > 10, i.e. N above
  # ~5,000 at max MAF 0.001: eas and afr, never sas -- is tested singly under
  # its own chr:pos:ref:alt ID, which names no region. Real-marker rows are
  # keyed by MarkerID: the rerun's rows replace the old rows with the same ID
  # (the rerun processed only the listed regions, so those ARE the rerun
  # regions' markers), and a changed one must be a recorded real-marker
  # failure or a marker at MAC <= ER_MAC. A row is a real-marker row when the
  # MarkerID's second ':'-field is a number (a position).
  local single=0; [[ ${sfx} == ".singleAssoc.txt.gz" ]] && single=1
  local isreal='(split($3, a_, ":") >= 2 && a_[2] ~ /^[0-9]+$/)'
  : > newmarkers.txt
  if (( single )); then awk -F'\t' 'FNR>1 && '"${isreal}"' {print $3}' new.tsv | sort -u > newmarkers.txt; fi
  # new pseudo-marker/region rows, keyed by region (only the listed regions are expected in new)
  awk -F'\t' -v rule="${rule}" -v single="${single}" 'NR==FNR{keep[$1]=1; next} FNR==1{next} single && '"${isreal}"' {next} {reg=('"${rule}"'); if(!(reg in keep)){bad++} } END{if(bad){print "ERROR: " bad " new rows outside the region list" > "/dev/stderr"; exit 1}}' regions.txt new.tsv
  # CONTROL 1 (row-level, since LEDGER #166): a rerun changes exactly two
  # kinds of row -- a cell whose saddlepoint failed under v4 (in failures.tsv)
  # and a cell at MAC <= --max_MAC_for_ER (ER_MAC, default 4), which now takes
  # the exact test whether or not it failed. Any OTHER changed row is a
  # refusal: it would mean the rerun differs from v4 for a reason we have not
  # named. Cells are keyed Region|Weight (results/pooled: the row's own
  # columns; singleAssoc pseudo-marker rows: the MarkerID's first and last
  # fields); a real-marker row is keyed by its MarkerID. MAC is read from the
  # OLD file: per set from the results, per marker from AC_Allele2 folded by
  # AF_Allele2 in the singleAssoc.
  ERMAC="${ER_MAC:-4}"
  if [[ ${sfx} == ".txt.gz" ]]; then
    awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["Group"]!="Cauchy"{print $h["Region"]"|"$h["Weight"]"\t"$h["MAC"]}' old.tsv > cellmac.tsv
  fi
  : > markermac.tsv
  if (( single )); then
    awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++)h[$i]=i; next} '"${isreal}"' {ac=$h["AC_Allele2"]+0; af=$h["AF_Allele2"]+0; mac=ac; if(af>0.5 && af<1) mac=ac*(1-af)/af; print $3"\t"mac}' old.tsv > markermac.tsv
  fi
  cat_any "${FAILURES}" | awk -F'\t' -v f="${pre}" 'NR>1 && $1==f {n=split($2,a,":"); if(!(n>=2 && a[2] ~ /^[0-9]+$/)) print a[1]"|"a[n]}' | sort -u > failcells.txt
  cat_any "${FAILURES}" | awk -F'\t' -v f="${pre}" 'NR>1 && $1==f {n=split($2,a,":"); if(n>=2 && a[2] ~ /^[0-9]+$/) print $2}' | sort -u > failmarkers.txt
  # changed rows = old rows the rerun replaces (listed regions' rows; for the
  # singleAssoc also the rerun markers' rows) that do not appear verbatim in new, and vice versa
  awk -F'\t' -v rule="${rule}" -v single="${single}" 'FILENAME=="regions.txt"{keep[$1]=1; next} FILENAME=="newmarkers.txt"{nm[$1]=1; next} FNR==1{next} { if(single && '"${isreal}"'){ if($3 in nm) print; next } reg=('"${rule}"'); if(reg in keep) print }' regions.txt newmarkers.txt old.tsv | sort > o.all
  awk -F'\t' 'FNR>1' new.tsv | sort > n.all
  comm -3 o.all n.all | sed 's/^\t//' > changed.rows
  nbad=$(awk -F'\t' -v sfx="${sfx}" -v single="${single}" 'FILENAME=="failcells.txt"{fc[$1]=1; next}
        FILENAME=="failmarkers.txt"{fm[$1]=1; next}
        FILENAME=="cellmac.tsv"{mac[$1]=$2; next}
        FILENAME=="markermac.tsv"{mmac[$1]=$2; next}
        { ok=0
          if(single && '"${isreal}"'){ key=$3; if(key in fm) ok=1; if((key in mmac) && mmac[key]+0 <= '"${ERMAC}"') ok=1 }
          else {
            if(single){ id=$3; n=split(id,a,":"); key=a[1]"|"a[n]; reg=a[1] }
            else if(sfx==".pooledTests.txt.gz" || $2=="Cauchy"){ key=$1; reg=$1 }
            else { key=$1"|"$4; reg=$1 }
            if(key in fc) ok=1
            if((key in mac) && mac[key]+0 <= '"${ERMAC}"') ok=1
            if(sfx==".pooledTests.txt.gz" || $2=="Cauchy"){ for(k in fc){ split(k,b,"|"); if(b[1]==reg) ok=1 } for(k in mac){ split(k,b,"|"); if(b[1]==reg && mac[k]+0 <= '"${ERMAC}"') ok=1 } }
          }
          if(!ok){ bad++; if(bad<=5) print "  unexplained change: " substr($0,1,120) > "/dev/stderr" } }
        END{print bad+0}' failcells.txt failmarkers.txt cellmac.tsv markermac.tsv changed.rows)
  if [[ ${nbad} -gt 0 ]]; then echo "ERROR ${sfx}: ${nbad} changed row(s) are neither a recorded failure nor a MAC <= ${ERMAC} cell/marker -- splice refused" >&2; exit 1; fi
  echo "${sfx}: $(wc -l < changed.rows) changed row(s), all explained (failed saddlepoint or MAC <= ${ERMAC})"
  # CONTROL 2: every recorded failure must have changed. The recorded failure
  # is the PSEUDO-MARKER's own single test (the scan reads the singleAssoc),
  # so it is checked exactly in the singleAssoc; in the results file the
  # burden cell of a set with rare, non-collapsed markers (Number_rare > 0)
  # is a different test that need not have failed, so there it is a note.
  if [[ ${sfx} == ".txt.gz" && -s failcells.txt ]]; then
    nmiss=$(awk -F'\t' 'NR==FNR{ch[$1"|"$4]=1; next} !($1 in ch){m++; if(m<=5) print "  failed pseudo-marker whose burden cell is unchanged: " $1 > "/dev/stderr"} END{print m+0}' changed.rows failcells.txt)
    [[ ${nmiss} -eq 0 ]] || echo "NOTE ${sfx}: ${nmiss} recorded pseudo-marker failure(s) sit in a burden cell that did not change (a set with Number_rare > 0 is a different test; the pseudo-marker row itself is checked in the singleAssoc)" >&2
  fi
  if (( single )) && [[ -s failcells.txt || -s failmarkers.txt ]]; then
    nmiss=$(awk -F'\t' 'NR==FNR{ if('"${isreal}"'){ ch[$3]=1 } else { n=split($3,a,":"); ch[a[1]"|"a[n]]=1 } next} !($1 in ch){m++; if(m<=5) print "  recorded failure unchanged in the singleAssoc: " $1 > "/dev/stderr"} END{print m+0}' changed.rows failcells.txt failmarkers.txt)
    [[ ${nmiss} -eq 0 ]] || echo "WARNING ${sfx}: ${nmiss} recorded failure(s) came back identical in the singleAssoc" >&2
  fi
  # ROW ACCOUNTING: a region-path row whose p is NA is not written (LEDGER
  # #127's convention; since #164 a saddlepoint failure above 13 carriers
  # yields such a row), so a rerun region can come back with FEWER rows than
  # it had. Print every region whose count changed, so the drop is visible.
  for r in $(cat regions.txt); do
    no=$(awk -F'\t' -v R="$r" 'FNR>1 && ('"${rule}"')==R' old.tsv | wc -l); nn=$(awk -F'\t' -v R="$r" 'FNR>1 && ('"${rule}"')==R' new.tsv | wc -l)
    if [[ ${no} -ne ${nn} ]]; then
      if [[ ${sfx} == ".pooledTests.txt.gz" ]]; then why="a region with nothing to disclose -- plain Cauchy, no refused cell -- writes no pooled row (LEDGER #127)"
      else why="NA rows are not written on the region path; a WHOLE region missing is a lost chromosome, read chrStatus (LEDGER #170)"; fi
      echo "ROWS ${sfx}: ${r} old ${no} new ${nn} (${why})"
    fi
  done
  if (( single )); then
    nold=$(awk -F'\t' 'NR==FNR{nm[$1]=1; next} FNR>1 && '"${isreal}"' && ($3 in nm)' newmarkers.txt old.tsv | wc -l)
    echo "${sfx}: real-marker rows: $(wc -l < newmarkers.txt) rerun marker(s), ${nold} old row(s) replaced by $(awk -F'\t' 'FNR>1 && '"${isreal}"'' new.tsv | wc -l) new"
  fi
  # the splice: old rows not replaced, then new rows, in old-file order of first appearance
  { head -1 old.tsv
    awk -F'\t' -v rule="${rule}" -v single="${single}" 'FILENAME=="regions.txt"{keep[$1]=1; next} FILENAME=="newmarkers.txt"{nm[$1]=1; next} FNR==1{next} { if(single && '"${isreal}"'){ if(!($3 in nm)) print; next } reg=('"${rule}"'); if(!(reg in keep)) print }' regions.txt newmarkers.txt old.tsv
    awk -F'\t' 'FNR>1' new.tsv
  } > "${out}"
  gzip -f "${out}"
  echo "${sfx}: old $(($(wc -l < old.tsv)-1)) rows, new $(($(wc -l < new.tsv)-1)) rows, spliced $(( $(gzip -dc "${out}.gz" | wc -l) - 1 )) rows"
}
splice_one ".txt.gz"             '$1'
if exists "${OLD_PREFIX}.singleAssoc.txt.gz"; then
  splice_one ".singleAssoc.txt.gz" 'substr($3, 1, index($3, ":") ? index($3, ":")-1 : length($3))'
fi
if exists "${OLD_PREFIX}.pooledTests.txt.gz"; then
  splice_one ".pooledTests.txt.gz" '$1'
fi
{ echo "spliced $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "old: ${OLD_PREFIX}"; echo "new: ${NEW_PREFIX}"
  echo "regions replaced: $(wc -l < regions.txt)"; echo "regions with failed saddlepoints: $(wc -l < failregions.txt)"
  echo "--- regions:"; cat regions.txt; } > "${OUT_DIR}/${pre}.splice.txt"
echo done
