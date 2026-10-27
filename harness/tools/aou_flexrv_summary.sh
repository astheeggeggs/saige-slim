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
#   FX_DIR, FX_PREFIX   optional: the FlexRV directory under BASE/ANC and the job prefix, for
#            a control fleet that follows neither naming rule (e.g. FX_DIR=fxr_syn FX_PREFIX=fxr_syn-);
#   SK_DIR, SK_PREFIX   the same for the SKAT-O arm (the simulated null: SK_DIR=sim_sk SK_PREFIX=sim_sk-)
#   AM_DIR, AM_PREFIX   the same for the am_moderate arm (the October fleet: oct_sk/oct_am/oct_fxa with oct_*- prefixes)
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
# the fixed/v7 precedence is the v4 era's splice of filtered reruns (sections 5-7); any other
# tag's FlexRV files are complete runs in flexrv_<tag>/ alone (2026-10-16: a v6 run listed the
# v4-era v7 directory -- harmless, the basenames never match, but wrong by construction)
if [[ ${RUNTAG} != "v4" ]]; then FXF="${FX}/__no_fixed_for_${RUNTAG}__"; FX7="${FX}/__no_v7_for_${RUNTAG}__"; fi
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
SIDECARS='\.(singleAssoc|markerList|pooledTests|skatoMethod|skatFailures|spaFallbacks|chrStatus|exactByWeight|exactByWeightAboveCap|spaPinned|stretchGate|spaRegionNA|spliceStatId|splice)\.txt(\.gz)?$'
# Only the run's own job files are results files (2026-10-16: a stray fxv5-* filtered rerun
# sat in sas/flexrv_v6 and would have been summarised as a second BenIntNeo). v4 jobs were
# named fxv4-/skv4-, later tags fx<N>-/sk<N>-/am<N>-; anything else in the directory is listed
# in the log as excluded and never read.
if [[ ${RUNTAG} == "v4" ]]; then PFX_FX="fxv4-"; PFX_SK="skv4-"; PFX_AM="amv4-"; else n="${RUNTAG#v}"; PFX_FX="fx${n}-"; PFX_SK="sk${n}-"; PFX_AM="am${n}-"; fi
# 2026-10-16: a control fleet whose jobs are named <tag>-... and whose files sit in <anc>/<tag>/
# (the FlexRV random-score arms fxr_syn / fxr_dmg / fxa_dmg) fits neither rule above; FX_DIR
# names its directory (under BASE/ANC) and FX_PREFIX its job prefix. The SKAT-O and am arms of
# such a tag do not exist and report 0 files.
[[ -n ${FX_DIR:-} ]] && FX="${BASE}/${ANC}/${FX_DIR}"
[[ -n ${FX_PREFIX:-} ]] && PFX_FX="${FX_PREFIX}"
[[ -n ${SK_DIR:-} ]] && SK="${BASE}/${ANC}/${SK_DIR}"        # the SKAT-O arm's twins (the simulated null: sim_sk / sim_sk-)
[[ -n ${SK_PREFIX:-} ]] && PFX_SK="${SK_PREFIX}"
[[ -n ${AM_DIR:-} ]] && AM="${BASE}/${ANC}/${AM_DIR}"        # and the am_moderate arm's (the October fleet: oct_am / oct_am-)
[[ -n ${AM_PREFIX:-} ]] && PFX_AM="${AM_PREFIX}"
keep(){ grep -E '\.txt\.gz$' | grep -vE "${SIDECARS}" | awk -v pfx="$1" '{ b=$0; sub(/.*\//, "", b); if (index(b, pfx) == 1) print; else print "excluded (not " pfx "*): " $0 > "/dev/stderr" }'; }
lsd "${FX}/*.txt.gz"  | keep "${PFX_FX}" > fx.list  2>> excluded.log || true
lsd "${FXF}/*.txt.gz" | keep "${PFX_FX}" > fxf.list 2>> excluded.log || true
lsd "${SK}/*.txt.gz"  | keep "${PFX_SK}" > sk.list  2>> excluded.log || true
lsd "${FX7}/*.txt.gz" | keep "${PFX_FX}" > fx7.list 2>> excluded.log || true
lsd "${AM}/*.txt.gz"  | keep "${PFX_AM}" > am.list  2>> excluded.log || true
if [[ -s excluded.log ]]; then echo "EXCLUDED $(wc -l < excluded.log | tr -d " ") file(s) whose basename is not this run's:"; cat excluded.log; fi
echo "flexrv: $(wc -l < fx.list | tr -d " ") v4 files, $(wc -l < fxf.list | tr -d " ") fixed, $(wc -l < fx7.list | tr -d " ") v7; skato: $(wc -l < sk.list | tr -d " ") files"
echo "am_moderate: $(wc -l < am.list | tr -d " ") files"
: > floor.long
# --- the chromosome scan (2026-10-16): every job's chrStatus sidecar, all three arms. A job
# that lost a chromosome still delivers and exits 0 (LEDGER #170, the IMG14 log10 defect), so
# this is the first thing read; a failed chromosome means the job is rerun BEFORE its trait
# is summarised. Written as a table and printed as one line; the scan is the VM's, not the
# notebook's (2,418 sidecars at ~2 s each killed a notebook session).
: > chrstatus.long
for arm in flexrv:"${FX}" flexrv_fixed:"${FXF}" flexrv_v7:"${FX7}" skato:"${SK}" am:"${AM}"; do
  name="${arm%%:*}"; dir="${arm#*:}"
  for f in $(lsd "${dir}/*.chrStatus.txt.gz"); do
    job=$(basename "${f}" | awk -F. '{print $1}')
    catz "${f}" | awk -F'\t' -v arm="${name}" -v job="${job}" 'BEGIN{OFS="\t"; req="?"; fl="?"} $1=="requested"{req=$2} $1=="failed"{fl=$2} END{print arm, job, req, fl}' >> chrstatus.long
  done
done
chrout="${OUT_DIR}/${ANC}.chrstatus.${RUNTAG}.tsv"
{ printf 'anc\tarm\tjob\trequested\tfailed\n'; awk -F'\t' -v anc="${ANC}" 'BEGIN{OFS="\t"}{print anc,$0}' chrstatus.long; } > "${chrout}"
gzip -f "${chrout}"
nbad=$(awk -F'\t' '$4!="none" && $4!="?"' chrstatus.long | wc -l | tr -d " ")
echo "chrStatus: $(wc -l < chrstatus.long | tr -d " ") jobs scanned, ${nbad} with a failed chromosome$( (( nbad > 0 )) && awk -F'\t' '$4!="none" && $4!="?"{printf ": %s/%s(chr %s)", $1, $2, $4}' chrstatus.long )"
# --- LEDGER #181 (2026-10-16): the rescale saddlepoints that returned nothing, all arms. The
# sidecar is written only when a run had any; absent for v4 and for IMG15 runs. One line, and
# a table <anc>.spaRegionNA.<tag>.tsv.gz of every named cell with the normal p it used.
: > sparegionna.long
for arm in flexrv:"${FX}" skato:"${SK}" am:"${AM}"; do
  name="${arm%%:*}"; dir="${arm#*:}"
  for f in $(lsd "${dir}/*.spaRegionNA.txt.gz"); do
    job=$(basename "${f}" | awk -F. '{print $1}'); tr=$(basename "${f}" | awk -F. '{print $2}')
    catz "${f}" | awk -F'\t' -v arm="${name}" -v job="${job}" -v tr="${tr}" 'BEGIN{OFS="\t"} NR>1 && $1!~/^#/ {print arm, job, tr, $0}' >> sparegionna.long
  done
done
nna=$(wc -l < sparegionna.long | tr -d " ")
if (( nna > 0 )); then
  naout="${OUT_DIR}/${ANC}.spaRegionNA.${RUNTAG}.tsv"
  { printf 'anc\tarm\tjob\ttrait\tKind\tRegion\tGroup\tmax_MAF\tWeight\tn_carriers\tn_pinned_dropped\tq_adj\tlog_p_normal_reported\tkernel_value\tdump\n'
    awk -F'\t' -v anc="${ANC}" 'BEGIN{OFS="\t"}{print anc,$0}' sparegionna.long; } > "${naout}"
  gzip -f "${naout}"
fi
echo "spaRegionNA: ${nna} rescale saddlepoint(s) returned nothing (LEDGER #181)$( (( nna > 0 )) && awk -F'\t' '{c[$1"/"$3]++} END{for(k in c) printf "; %s %d", k, c[k]}' sparegionna.long )"
trait_of(){ basename "$1" | awk -F. '{print $2}'; }   # <job>.<TRAIT>.<anc>.chr<...>.txt.gz
# --- the three per-file passes, PARALLEL streams (2026-10-16): the job's time is gcloud cat
# one file at a time, so PARALLEL (default 8) files are streamed at once; each worker writes
# its own numbered part and the parts are concatenated in list order, so every table is
# byte-identical to the serial version's (the joins and sorts below never depended on order).
PARALLEL="${PARALLEL:-8}"
mkdir -p parts
do_fx(){   # $1 = index, $2 = file
  local f="$2" pre src=v4
  pre=$(basename "${f}")
  if grep -qxF "${FX7}/${pre}" fx7.list; then f="${FX7}/${pre}"; src=v7
  elif grep -qxF "${FXF}/${pre}" fxf.list; then f="${FXF}/${pre}"; src=fixed; fi
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" -v src="${src}" -v floor="parts/$1.fx.floor" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; if(!("Weight" in h)){print "ERROR: no Weight column in " FILENAME > "/dev/stderr"; exit 1} next}
    { g=$h["Region"]; p=$h["Pvalue_Burden"]
      if(("log10_Pvalue_Burden" in h) && $h["log10_Pvalue_Burden"]!="NA" && $h["log10_Pvalue_Burden"]+0 < -300)
        printf "flexrv\t%s\t%s\t%s\t%s\n", tr, g, (($h["Group"]=="Cauchy")?"Cauchy":$h["Group"] "@" $h["max_MAF"] "@" $h["Weight"]), $h["log10_Pvalue_Burden"] >> floor
      if($h["Group"]=="Cauchy"){ pc[g]=p; next }
      n[g]++; if(p!="NA" && (!(g in pb) || p+0 < pb[g]+0)){ pb[g]=p; wb[g]=$h["Weight"] } }
    END{ for(g in n) printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\n", tr, g, src, ((g in pc)?pc[g]:"NA"), ((g in pb)?pb[g]:"NA"), ((g in wb)?wb[g]:"NA"), n[g] }' > "parts/$1.fx.tsv"
}
do_sk(){   # $1 = index, $2 = file: the per-gene best AND the cells
  local f="$2"
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" -v cells="parts/$1.sk.cells" -v floor="parts/$1.sk.floor" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    { g=$h["Region"]; p=$h["Pvalue"]
      if(("log10_Pvalue" in h) && $h["log10_Pvalue"]!="NA" && $h["log10_Pvalue"]+0 < -300)
        printf "skato\t%s\t%s\t%s\t%s\n", tr, g, (($h["Group"]=="Cauchy")?"Cauchy":$h["Group"] "@" $h["max_MAF"]), $h["log10_Pvalue"] >> floor
      if($h["Group"]=="Cauchy"){ pc[g]=p; next }
      cell=$h["Group"] "@" $h["max_MAF"]
      printf "%s\t%s\t%s\t%s\t%s\n", tr, g, cell, $h["Pvalue_Burden"], p >> cells
      if(p!="NA" && (!(g in pb) || p+0 < pb[g]+0)){ pb[g]=p; cb[g]=cell }
      q=$h["Pvalue_Burden"]; if(q!="NA" && (!(g in bb) || q+0 < bb[g]+0)) bb[g]=q
      q=$h["Pvalue_SKAT"];   if(q!="NA" && (!(g in sb) || q+0 < sb[g]+0)) sb[g]=q
      seen[g]=1 }
    END{ for(g in seen) printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", tr, g, ((g in pc)?pc[g]:"NA"), ((g in pb)?pb[g]:"NA"), ((g in cb)?cb[g]:"NA"), ((g in bb)?bb[g]:"NA"), ((g in sb)?sb[g]:"NA") }' > "parts/$1.sk.tsv"
}
do_am(){   # $1 = index, $2 = file: cells only (same schema as SKAT-O)
  local f="$2"
  catz "${f}" | awk -F'\t' -v tr="$(trait_of "${f}")" -v cells="parts/$1.am.cells" -v floor="parts/$1.am.floor" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    { g=$h["Region"]; p=$h["Pvalue"]
      if(("log10_Pvalue" in h) && $h["log10_Pvalue"]!="NA" && $h["log10_Pvalue"]+0 < -300)
        printf "am\t%s\t%s\t%s\t%s\n", tr, g, (($h["Group"]=="Cauchy")?"Cauchy":$h["Group"] "@" $h["max_MAF"]), $h["log10_Pvalue"] >> floor
      if($h["Group"]=="Cauchy") next
      printf "%s\t%s\t%s\t%s\t%s\n", tr, g, $h["Group"] "@" $h["max_MAF"], $h["Pvalue_Burden"], p >> cells }'
}
export -f do_fx do_sk do_am catz trait_of
export FX7 FXF
progress(){   # $1 = index, $2 = total, $3 = label: one line every 50 files (dsub ships stdout every minute)
  if (( 10#$1 % 50 == 0 || 10#$1 == $2 )); then echo "  $3: $((10#$1))/$2 files streamed ($(date -u +%H:%M:%S))"; fi
}
export -f progress
run_parallel(){   # $1 = list file, $2 = worker; one numbered argument pair per file
  [[ -s $1 ]] || return 0
  local total; total=$(wc -l < "$1" | tr -d " ")
  nl -ba -nrz -w6 "$1" | xargs -P "${PARALLEL}" -n 2 bash -c "$2 \"\$@\" && progress \"\$1\" ${total} $2" _
}
run_parallel fx.list do_fx
run_parallel sk.list do_sk
run_parallel am.list do_am
# concatenate the parts in list order (byte-identical to the serial loops' output)
cat_parts(){ local sfx="$1" f; for f in $(ls parts/ 2>/dev/null | grep -E "\.${sfx}$" | sort); do cat "parts/${f}"; done; }
cat_parts fx.tsv   > fx.tsv
cat_parts sk.tsv   > sk.tsv
{ cat_parts sk.cells; cat_parts am.cells; } > cells.long
{ cat_parts fx.floor; cat_parts sk.floor; cat_parts am.floor; } >> floor.long
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
