#!/bin/bash
# Gene-level summary of one ancestry's FlexRV and SKAT-O arms, for the power
# comparison plots. Runs INSIDE a dsub job (cloud-cli image); streams every
# results file with gcloud storage cat and writes ONE small TSV: one row per
# (trait, gene), aggregate values only -- the plotting notebook reads this and
# never the per-trait results files (the 2026-09-26 egress warning coincided
# with plotting from those).
#   ANC      ancestry (afr|amr|eas|eur|sas)
#   BASE     gs://.../SAIGE_SLIM_TESTING   (or a local dir for the self-test)
#   RUNTAG   v4 (default)
#   OUT_DIR  local dir dsub delocalizes (--output-recursive)
# FlexRV: the final file of a job is flexrv_<tag>_fixed/<pre>.txt.gz when a
# splice exists (binary traits with rerun regions), else flexrv_<tag>/<pre>.txt.gz
# (jobs with nothing to rerun, every quantitative trait); `source` says which.
# Columns: anc trait gene source p_flexrv p_flexrv_best set_best n_sets
#          p_skato p_skato_best cell_best p_burden_best p_skat_best
set -euo pipefail
: "${ANC:?}"; : "${BASE:?}"; : "${OUT_DIR:?}"; RUNTAG="${RUNTAG:-v4}"
mkdir -p "${OUT_DIR}"
if [[ ${BASE} == gs://* ]]; then
  lsd(){ gcloud storage ls "$1" 2>/dev/null || true; }
  catz(){ gcloud storage cat "$1" | gzip -dc; }
else
  lsd(){ local f; for f in $1; do [[ -e ${f} ]] && echo "${f}"; done; true; }   # $1 unquoted: the glob expands here
  catz(){ gzip -dc "$1"; }
fi
FX="${BASE}/${ANC}/flexrv_${RUNTAG}"; FXF="${BASE}/${ANC}/flexrv_${RUNTAG}_fixed"; SK="${BASE}/${ANC}/skato_${RUNTAG}"
lsd "${FX}/*.txt.gz"  | grep -E '\.txt\.gz$' | grep -vE '\.(singleAssoc|markerList|pooledTests|skatoMethod|skatFailures|spaFallbacks|chrStatus)\.txt\.gz$' > fx.list || true
lsd "${FXF}/*.txt.gz" | grep -E '\.txt\.gz$' | grep -vE '\.(singleAssoc|markerList|pooledTests|skatoMethod|skatFailures|spaFallbacks|chrStatus)\.txt\.gz$' > fxf.list || true
lsd "${SK}/*.txt.gz"  | grep -E '\.txt\.gz$' | grep -vE '\.(singleAssoc|markerList|pooledTests|skatoMethod|skatFailures|spaFallbacks|chrStatus)\.txt\.gz$' > sk.list || true
echo "flexrv: $(wc -l < fx.list) v4 files, $(wc -l < fxf.list) fixed; skato: $(wc -l < sk.list) files"
trait_of(){ basename "$1" | awk -F. '{print $2}'; }   # <job>.<TRAIT>.<anc>.chr<...>.txt.gz
# --- FlexRV: per gene, the Cauchy row's p, the best set and how many sets ----
: > fx.tsv
while read -r f; do
  [[ -n ${f} ]] || continue
  pre=$(basename "${f}"); src=v4
  if grep -qxF "${FXF}/${pre}" fxf.list; then f="${FXF}/${pre}"; src=fixed; fi
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" -v src="${src}" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; if(!("Weight" in h)){print "ERROR: no Weight column in " FILENAME > "/dev/stderr"; exit 1} next}
    { g=$h["Region"]; p=$h["Pvalue_Burden"]
      if($h["Group"]=="Cauchy"){ pc[g]=p; next }
      n[g]++; if(p!="NA" && (!(g in pb) || p+0 < pb[g]+0)){ pb[g]=p; wb[g]=$h["Weight"] } }
    END{ for(g in n) printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\n", tr, g, src, ((g in pc)?pc[g]:"NA"), ((g in pb)?pb[g]:"NA"), ((g in wb)?wb[g]:"NA"), n[g] }' >> fx.tsv
done < fx.list
# --- SKAT-O: per gene, the Cauchy row, the best cell, best burden, best SKAT --
: > sk.tsv
while read -r f; do
  [[ -n ${f} ]] || continue
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    { g=$h["Region"]; p=$h["Pvalue"]
      if($h["Group"]=="Cauchy"){ pc[g]=p; next }
      cell=$h["Group"] "@" $h["max_MAF"]
      if(p!="NA" && (!(g in pb) || p+0 < pb[g]+0)){ pb[g]=p; cb[g]=cell }
      q=$h["Pvalue_Burden"]; if(q!="NA" && (!(g in bb) || q+0 < bb[g]+0)) bb[g]=q
      q=$h["Pvalue_SKAT"];   if(q!="NA" && (!(g in sb) || q+0 < sb[g]+0)) sb[g]=q
      seen[g]=1 }
    END{ for(g in seen) printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", tr, g, ((g in pc)?pc[g]:"NA"), ((g in pb)?pb[g]:"NA"), ((g in cb)?cb[g]:"NA"), ((g in bb)?bb[g]:"NA"), ((g in sb)?sb[g]:"NA") }' >> sk.tsv
done < sk.list
# --- join on (trait, gene); a gene in one arm only keeps NA on the other -----
out="${OUT_DIR}/${ANC}.flexrv_vs_skato.${RUNTAG}.tsv"
awk -F'\t' -v anc="${ANC}" 'BEGIN{OFS="\t"; print "anc","trait","gene","source","p_flexrv","p_flexrv_best","set_best","n_sets","p_skato","p_skato_best","cell_best","p_burden_best","p_skat_best"}
  FILENAME=="fx.tsv"{k=$1"\t"$2; fx[k]=$3"\t"$4"\t"$5"\t"$6"\t"$7; keys[k]=1; next}
  FILENAME=="sk.tsv"{k=$1"\t"$2; sk[k]=$3"\t"$4"\t"$5"\t"$6"\t"$7; keys[k]=1; next}
  END{ for(k in keys){ print anc, k, ((k in fx)?fx[k]:"NA\tNA\tNA\tNA\tNA"), ((k in sk)?sk[k]:"NA\tNA\tNA\tNA\tNA") } }' fx.tsv sk.tsv | (read -r hdr; echo "${hdr}"; sort -t $'\t' -k2,2 -k3,3) > "${out}"
gzip -f "${out}"
n=$(( $(gzip -dc "${out}.gz" | wc -l) - 1 ))
both=$(gzip -dc "${out}.gz" | awk -F'\t' 'NR>1 && $5!="NA" && $9!="NA"{n++} END{print n+0}')
echo "rows: ${n} (trait, gene) pairs; ${both} with both arms; traits: $(gzip -dc "${out}.gz" | awk -F'\t' 'NR>1{t[$2]=1} END{n=0; for(k in t)n++; print n}')"
echo done
