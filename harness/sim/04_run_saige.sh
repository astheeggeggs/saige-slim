#!/usr/bin/env bash
# Steps 0/1/2 over the simulated cohort, following the upstream SAIGE recipe.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
export D="${SIMDIR:-${D:-$REPO/harness/work/sim}}"
. harness/sim/00_env.sh
# absolute, via 00_env.sh -- see the note in 03_plink_prep.sh
RS=${RS:-$REPO/.pixi/envs/default/lib/R/bin/Rscript}
D=${D:-harness/work/sim}
# Overrides used by harness/sim/upstream_compare to point one command at another
# build (its extdata/) and at model/output prefixes outside the fixture dir.
E=${SAIGE_EXTDATA:-extdata}
NT="${NTHREADS:-6}"

CUT="${RELCUT:-0.05}"
# chr20 has only ~3,055 markers that are pairwise r2<0.05 within 1 Mb, so the
# recipe's 5000 random markers for the relatedness screen is not attainable.
NRM="${NRANDMARK:-3000}"
# GRMFILE overrides the step-0 output, so steps 1-2 can be pointed at the
# pedigree-truth matrix (03b_pedigree_grm.py) instead of an estimated one.
GRM=${GRMFILE:-$D/sparseGRM_relatednessCutoff_${CUT}_${NRM}_randomMarkersUsed.sparseGRM.mtx}
IDS=${GRMIDS:-$GRM.sampleIDs.txt}

step0() {
  echo "=== step 0: createSparseGRM ==="
  $RS $E/createSparseGRM.R \
    --plinkFile=$D/data_for_grm \
    --nThreads=$NT \
    --outputPrefix=$D/sparseGRM \
    --numRandomMarkerforSparseKin=$NRM \
    --relatednessCutoff=$CUT
}

# step1 <phenoCol> <traitType>
step1() {
  local col=$1 tt=$2 extra=""
  [ "$tt" = "quantitative" ] && extra="--invNormalize=TRUE"
  echo "=== step 1: $col ($tt) ==="
  $RS $E/step1_fitNULLGLMM.R \
    --plinkFile=$D/data_for_vr \
    --sparseGRMFile=$GRM \
    --sparseGRMSampleIDFile=$IDS \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    --phenoFile=$D/pheno_pcs.txt \
    --phenoCol=$col \
    --sampleIDColinphenoFile=IID \
    --traitType=$tt $extra \
    --covarColList=age,sex,PC1,PC2,PC3,PC4 \
    --qCovarColList=sex \
    --skipVarianceRatioEstimation=FALSE \
    --isCateVarianceRatio=TRUE \
    --relatednessCutoff=$CUT \
    --IsOverwriteVarianceRatioFile=TRUE \
    --outputPrefix=${S1OUTPREFIX:-$D/step1_$col}
}

# step2_group <phenoCol>
step2_group() {
  local col=$1
  echo "=== step 2 (SAIGE-GENE+, set + single-in-set): $col ==="
  $RS $E/step2_SPAtests.R \
    --bedFile=$D/chr20_sim.bed \
    --bimFile=$D/chr20_sim.bim \
    --famFile=$D/chr20_sim.fam \
    --chrom=20 \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=${MODELPREFIX:-$D/step1_$col}.rda \
    --varianceRatioFile=${MODELPREFIX:-$D/step1_$col}.varianceRatio.txt \
    --sparseGRMFile=$GRM \
    --sparseGRMSampleIDFile=$IDS \
    --LOCO=FALSE \
    --is_Firth_beta=TRUE \
    --pCutoffforFirth=0.1 \
    --is_output_moreDetails=TRUE \
    --is_fastTest=TRUE \
    --is_output_markerList_in_groupTest=TRUE \
    --is_single_in_groupTest=TRUE \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 \
    --groupFile=$D/groupFile.txt \
    --annotation_in_groupTest="high,medium,low,high:medium,high:medium:low" \
    ${GROUP_NT:+--nThreads=$GROUP_NT} \
    ${S2EXTRA:-} \
    --SAIGEOutputFile=${S2OUT:-$D/step2_group_$col.txt}
}

# step2_single <phenoCol>
step2_single() {
  local col=$1
  echo "=== step 2 (SAIGE-variant, all variants MAC>=20): $col ==="
  $RS $E/step2_SPAtests.R \
    --bedFile=$D/chr20_sim.bed \
    --bimFile=$D/chr20_sim.bim \
    --famFile=$D/chr20_sim.fam \
    --chrom=20 \
    --minMAF=0 \
    --minMAC=20 \
    --GMMATmodelFile=${MODELPREFIX:-$D/step1_$col}.rda \
    --varianceRatioFile=${MODELPREFIX:-$D/step1_$col}.varianceRatio.txt \
    --sparseGRMFile=$GRM \
    --sparseGRMSampleIDFile=$IDS \
    --LOCO=FALSE \
    --is_Firth_beta=TRUE \
    --pCutoffforFirth=0.1 \
    --is_fastTest=TRUE \
    --nThreads=$NT \
    ${S2EXTRA:-} \
    --SAIGEOutputFile=${S2OUT:-$D/step2_single_$col.txt}
}

"$@"
