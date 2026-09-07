#!/usr/bin/env bash
# sim_rep.sh <label> <id> <beta_b> <full|s1>: one replicate -- simulate the
# phenotype, refit step 1 (the binary_sparseGRM_vr configuration), and for
# `full` run the instrumented SKAT-O region test (golden P6's command, no fix
# installed; the fixes are applied offline to the dumps) -- then rep_stats.R.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="${DIAG_OUT:-$REPO/harness/work/diag-burden}"
label=$1; id=$2; beta=$3; mode=$4
prep="${SIM_PREP:-$OUT/sim/prep.rds}"
rep="$OUT/sim/$label/rep_$(printf '%04d' "$id")"
case "$label" in null) seed=$((1000 + id)) ;; planted) seed=$((2000 + id)) ;; *) seed=$((3000 + id)) ;; esac
rm -rf "$rep"; mkdir -p "$rep"
cd "$REPO/extdata"
source "$REPO/harness/env.sh"
R="pixi run --manifest-path $REPO/pixi.toml Rscript"
$R "$HERE/sim_pheno.R" "$prep" "$seed" "$beta" "$rep/pheno.txt" > "$rep/log.txt" 2>&1
$R step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE --useSparseGRMforVarRatio=TRUE \
  --phenoFile="$rep/pheno.txt" --phenoCol=y_binary --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 --sampleIDColinphenoFile=IID --traitType=binary \
  --IsOverwriteVarianceRatioFile=TRUE --isCateVarianceRatio=TRUE \
  --outputPrefix="$rep/m" >> "$rep/log.txt" 2>&1 || { echo "$label $id step1 FAILED"; exit 1; }
if [ "$mode" = "full" ]; then
  SAIGE_DUMP="$rep" SAIGE_REPO="$REPO" $R "$HERE/drv.R" \
    --bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim \
    --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first \
    --SAIGEOutputFile="$rep/out.txt" --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile="$rep/m.rda" --varianceRatioFile="$rep/m.varianceRatio.txt" \
    --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=./input/group_new_snpid.txt \
    --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 --is_fastTest=FALSE >> "$rep/log.txt" 2>&1 \
    || { echo "$label $id step2 FAILED"; exit 1; }
fi
$R "$HERE/rep_stats.R" "$rep" "$prep" >> "$rep/log.txt" 2>&1 || { echo "$label $id stats FAILED"; exit 1; }
rm -f "$rep/m.rda" "$rep/pheno.txt" "$rep/out.txt.index"
echo "$label $id ok"
