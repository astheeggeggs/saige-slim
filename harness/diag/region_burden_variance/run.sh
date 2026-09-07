#!/usr/bin/env bash
# run.sh <name> <r.corr> <varianceRatioFile> <minMACExclude> <maxMACInclude> [extra flags]
# One instrumented step-2 run of golden P6's (r.corr=0) / P5's (r.corr=1)
# command on the committed binary model, writing to $DIAG_OUT/<name>/.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="${DIAG_OUT:-$REPO/harness/work/diag-burden}"
name=$1; rcorr=$2; vr=$3; minexc=$4; maxinc=$5; shift 5
run="$OUT/$name"; rm -rf "$run"; mkdir -p "$run"
cd "$REPO/extdata"
source "$REPO/harness/env.sh"
SAIGE_DUMP="$run" SAIGE_REPO="$REPO" pixi run --manifest-path "$REPO/pixi.toml" Rscript "$HERE/drv.R" \
  --bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim \
  --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first \
  --SAIGEOutputFile="$run/out.txt" --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  --GMMATmodelFile=./output/example_binary_sparseGRM.rda --varianceRatioFile="$vr" \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr="$rcorr" --is_fastTest=FALSE \
  --cateVarRatioMinMACVecExclude="$minexc" --cateVarRatioMaxMACVecInclude="$maxinc" "$@" \
  > "$run/log.txt" 2>&1
echo "$name exit=$?"
