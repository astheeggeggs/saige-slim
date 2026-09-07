#!/usr/bin/env bash
# run_quant.sh <name> <r.corr>: golden P37's (r.corr=0) / P36's (r.corr=1)
# command on the quantitative step-1 golden -- the negative control (no
# R-side Phi adjustment on quantitative traits, no sparse GRM).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="${DIAG_OUT:-$REPO/harness/work/diag-burden}"
name=$1; rcorr=$2
run="$OUT/$name"; rm -rf "$run"; mkdir -p "$run"
cd "$REPO/extdata"
source "$REPO/harness/env.sh"
SAIGE_DUMP="$run" SAIGE_REPO="$REPO" pixi run --manifest-path "$REPO/pixi.toml" Rscript "$HERE/drv.R" \
  --bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim \
  --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first \
  --SAIGEOutputFile="$run/out.txt" --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  --GMMATmodelFile="$REPO/harness/goldens/step1/quant_sparseGRM.rda" \
  --varianceRatioFile="$REPO/harness/goldens/step1/quant_sparseGRM.varianceRatio.txt" \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 --r.corr="$rcorr" --is_fastTest=FALSE \
  > "$run/log.txt" 2>&1
echo "$name exit=$?"
