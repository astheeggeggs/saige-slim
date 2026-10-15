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
# FlexRV: the final file of a job is, in this order of precedence,
#   flexrv_v7_checked/<pre>.txt.gz   the complete rerun after LEDGER #172 (binary
#                                    traits), verified by the splice control
#                                    (`source` = v7; FLEXRV_V7_DIR names the dir)
#   flexrv_<tag>_fixed/<pre>.txt.gz  the splice of the filtered reruns (fixed)
#   flexrv_<tag>/<pre>.txt.gz        the original (v4): jobs with nothing to
#                                    rerun, every quantitative trait
# `source` says which.
# Columns: anc trait gene source p_flexrv p_flexrv_best set_best n_sets
#          p_skato p_skato_best cell_best p_burden_best p_skat_best
# A second table, <anc>.skato_cells.<tag>.tsv.gz, holds every SKAT-O cell of
# every gene, wide: anc trait gene then burden:<mask>@<maxMAF> and
# skato:<mask>@<maxMAF> for each cell the ancestry's files contain (the QQ
# grids per mask and max MAF are drawn from it). Gene level, no counts.
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
FX7="${BASE}/${ANC}/${FLEXRV_V7_DIR:-flexrv_v7_checked}"
# 2026-10-15 (v6): a third arm, the am_moderate SKAT-O run (one mask, no single-variant
# file), read into the CELLS table only (its cells are `am_moderate@<maf>`), never into the
# per-gene best of the main SKAT-O arm; absent for v4. And a FLOOR table: every row whose
# log10 column (IMG14+) is below -300, i.e. past the double floor where the printed p is a
# denormal or 0 -- the genes the --is_output_log10P flag exists for. Empty for v4.
AM="${BASE}/${ANC}/am_${RUNTAG}"
# every sidecar a step-2 run or the splice control can leave beside the results
# file; a new one that is not listed here is counted as a results file and the
# `K v7` provenance line reads STALE (2026-10-04: spaPinned, exactByWeightAboveCap
# and stretchGate were missing once flexrv_v8_checked carried them)
SIDECARS='\.(singleAssoc|markerList|pooledTests|skatoMethod|skatFailures|spaFallbacks|chrStatus|exactByWeight|exactByWeightAboveCap|spaPinned|stretchGate|spliceStatId|splice)\.txt(\.gz)?$'
lsd "${FX}/*.txt.gz"  | grep -E '\.txt\.gz$' | grep -vE "${SIDECARS}" > fx.list || true
lsd "${FXF}/*.txt.gz" | grep -E '\.txt\.gz$' | grep -vE "${SIDECARS}" > fxf.list || true
lsd "${SK}/*.txt.gz"  | grep -E '\.txt\.gz$' | grep -vE "${SIDECARS}" > sk.list || true
lsd "${FX7}/*.txt.gz" | grep -E '\.txt\.gz$' | grep -vE "${SIDECARS}" > fx7.list || true
lsd "${AM}/*.txt.gz"  | grep -E '\.txt\.gz$' | grep -vE "${SIDECARS}" > am.list || true
echo "flexrv: $(wc -l < fx.list) v4 files, $(wc -l < fxf.list) fixed, $(wc -l < fx7.list) v7; skato: $(wc -l < sk.list) files"
echo "am_moderate: $(wc -l < am.list | tr -d " ") files"
: > floor.long
trait_of(){ basename "$1" | awk -F. '{print $2}'; }   # <job>.<TRAIT>.<anc>.chr<...>.txt.gz
# --- FlexRV: per gene, the Cauchy row's p, the best set and how many sets ----
: > fx.tsv
while read -r f; do
  [[ -n ${f} ]] || continue
  pre=$(basename "${f}"); src=v4
  if grep -qxF "${FX7}/${pre}" fx7.list; then f="${FX7}/${pre}"; src=v7
  elif grep -qxF "${FXF}/${pre}" fxf.list; then f="${FXF}/${pre}"; src=fixed; fi
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" -v src="${src}" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; if(!("Weight" in h)){print "ERROR: no Weight column in " FILENAME > "/dev/stderr"; exit 1} next}
    { g=$h["Region"]; p=$h["Pvalue_Burden"]
      if(("log10_Pvalue_Burden" in h) && $h["log10_Pvalue_Burden"]!="NA" && $h["log10_Pvalue_Burden"]+0 < -300)
        printf "flexrv\t%s\t%s\t%s\t%s\n", tr, g, (($h["Group"]=="Cauchy")?"Cauchy":$h["Group"] "@" $h["max_MAF"] "@" $h["Weight"]), $h["log10_Pvalue_Burden"] >> "floor.long"
      if($h["Group"]=="Cauchy"){ pc[g]=p; next }
      n[g]++; if(p!="NA" && (!(g in pb) || p+0 < pb[g]+0)){ pb[g]=p; wb[g]=$h["Weight"] } }
    END{ for(g in n) printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\n", tr, g, src, ((g in pc)?pc[g]:"NA"), ((g in pb)?pb[g]:"NA"), ((g in wb)?wb[g]:"NA"), n[g] }' >> fx.tsv
done < fx.list
# --- SKAT-O: per gene, the Cauchy row, the best cell, best burden, best SKAT --
: > sk.tsv; : > cells.long
while read -r f; do
  [[ -n ${f} ]] || continue
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    { g=$h["Region"]; p=$h["Pvalue"]
      if(("log10_Pvalue" in h) && $h["log10_Pvalue"]!="NA" && $h["log10_Pvalue"]+0 < -300)
        printf "skato\t%s\t%s\t%s\t%s\n", tr, g, (($h["Group"]=="Cauchy")?"Cauchy":$h["Group"] "@" $h["max_MAF"]), $h["log10_Pvalue"] >> "floor.long"
      if($h["Group"]=="Cauchy"){ pc[g]=p; next }
      cell=$h["Group"] "@" $h["max_MAF"]
      printf "%s\t%s\t%s\t%s\t%s\n", tr, g, cell, $h["Pvalue_Burden"], p >> "cells.long"
      if(p!="NA" && (!(g in pb) || p+0 < pb[g]+0)){ pb[g]=p; cb[g]=cell }
      q=$h["Pvalue_Burden"]; if(q!="NA" && (!(g in bb) || q+0 < bb[g]+0)) bb[g]=q
      q=$h["Pvalue_SKAT"];   if(q!="NA" && (!(g in sb) || q+0 < sb[g]+0)) sb[g]=q
      seen[g]=1 }
    END{ for(g in seen) printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", tr, g, ((g in pc)?pc[g]:"NA"), ((g in pb)?pb[g]:"NA"), ((g in cb)?cb[g]:"NA"), ((g in bb)?bb[g]:"NA"), ((g in sb)?sb[g]:"NA") }' >> sk.tsv
done < sk.list
# --- am_moderate: its cells into the wide table only (same schema as SKAT-O) --
while read -r f; do
  [[ -n ${f} ]] || continue
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    { g=$h["Region"]; p=$h["Pvalue"]
      if(("log10_Pvalue" in h) && $h["log10_Pvalue"]!="NA" && $h["log10_Pvalue"]+0 < -300)
        printf "am\t%s\t%s\t%s\t%s\n", tr, g, (($h["Group"]=="Cauchy")?"Cauchy":$h["Group"] "@" $h["max_MAF"]), $h["log10_Pvalue"] >> "floor.long"
      if($h["Group"]=="Cauchy") next
      printf "%s\t%s\t%s\t%s\t%s\n", tr, g, $h["Group"] "@" $h["max_MAF"], $h["Pvalue_Burden"], p >> "cells.long" }' 
done < am.list
# --- the floor table: rows past the double floor, by their log10 columns ------
floorout="${OUT_DIR}/${ANC}.floor.${RUNTAG}.tsv"
{ printf 'anc\tarm\ttrait\tgene\tcell\tlog10_p\n'; awk -F'\t' -v anc="${ANC}" 'BEGIN{OFS="\t"}{print anc,$0}' floor.long | sort -t $'\t' -k6,6g; } > "${floorout}"
gzip -f "${floorout}"
echo "floor table: $(wc -l < floor.long | tr -d " ") row(s) past log10 p = -300"
# --- join on (trait, gene); a gene in one arm only keeps NA on the other -----
out="${OUT_DIR}/${ANC}.flexrv_vs_skato.${RUNTAG}.tsv"
awk -F'\t' -v anc="${ANC}" 'BEGIN{OFS="\t"; print "anc","trait","gene","source","p_flexrv","p_flexrv_best","set_best","n_sets","p_skato","p_skato_best","cell_best","p_burden_best","p_skat_best"}
  FILENAME=="fx.tsv"{k=$1"\t"$2; fx[k]=$3"\t"$4"\t"$5"\t"$6"\t"$7; keys[k]=1; next}
  FILENAME=="sk.tsv"{k=$1"\t"$2; sk[k]=$3"\t"$4"\t"$5"\t"$6"\t"$7; keys[k]=1; next}
  END{ for(k in keys){ print anc, k, ((k in fx)?fx[k]:"NA\tNA\tNA\tNA\tNA"), ((k in sk)?sk[k]:"NA\tNA\tNA\tNA\tNA") } }' fx.tsv sk.tsv | (read -r hdr; echo "${hdr}"; sort -t $'\t' -k2,2 -k3,3) > "${out}"
gzip -f "${out}"
# --- every SKAT-O cell, wide: the cell list is whatever the files contain ---
cut -f3 cells.long | sort -u > cells.list
cellsout="${OUT_DIR}/${ANC}.skato_cells.${RUNTAG}.tsv"
sort -t $'\t' -k1,1 -k2,2 -S 1G cells.long | awk -F'\t' -v anc="${ANC}" 'BEGIN{OFS="\t"}
  FILENAME=="cells.list"{nc++; cl[nc]=$1; next}
  FNR==1{ printf "anc\ttrait\tgene"; for(i=1;i<=nc;i++) printf "\tburden:%s", cl[i]; for(i=1;i<=nc;i++) printf "\tskato:%s", cl[i]; printf "\n" }
  { k=$1"\t"$2
    if(k!=cur){ if(cur!=""){ flush() } cur=k; delete b; delete s }
    b[$3]=$4; s[$3]=$5 }
  function flush(   i){ printf "%s\t%s", anc, cur; for(i=1;i<=nc;i++) printf "\t%s", ((cl[i] in b)?b[cl[i]]:"NA"); for(i=1;i<=nc;i++) printf "\t%s", ((cl[i] in s)?s[cl[i]]:"NA"); printf "\n" }
  END{ if(cur!="") flush() }' cells.list - > "${cellsout}"
gzip -f "${cellsout}"
echo "cells table: $(( $(gzip -dc "${cellsout}.gz" | wc -l) - 1 )) (trait, gene) rows x $(( $(wc -l < cells.list) )) cells"
n=$(( $(gzip -dc "${out}.gz" | wc -l) - 1 ))
both=$(gzip -dc "${out}.gz" | awk -F'\t' 'NR>1 && $5!="NA" && $9!="NA"{n++} END{print n+0}')
echo "rows: ${n} (trait, gene) pairs; ${both} with both arms; traits: $(gzip -dc "${out}.gz" | awk -F'\t' 'NR>1{t[$2]=1} END{n=0; for(k in t)n++; print n}')"
echo done
