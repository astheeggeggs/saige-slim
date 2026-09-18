#!/bin/bash
# Profiled sim50k SKAT-O run at HEAD: pgen reader, y_binary, the 04_run_saige.sh
# step2_group flags (fastTest, single-in-group, sparse GRM), 1 thread.
set -uo pipefail
REPO=/Users/dpalmer/Repositories/SAIGE
D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
S=$1; arm=$2
cd $REPO; source harness/env.sh; cd extdata
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
ARGS="--pgenFile=$D/chr20_sim_pgen.pgen --pvarFile=$D/chr20_sim_pgen.pvar --psamFile=$D/chr20_sim_pgen.psam \
 --chrom=20 --minMAF=0 --minMAC=0.5 \
 --GMMATmodelFile=$D/step1_y_binary.rda --varianceRatioFile=$D/step1_y_binary.varianceRatio.txt \
 --sparseGRMFile=$D/pedigree_truth.sparseGRM.mtx --sparseGRMSampleIDFile=$D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt \
 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_output_moreDetails=TRUE --is_fastTest=TRUE \
 --is_output_markerList_in_groupTest=TRUE --is_single_in_groupTest=TRUE \
 --maxMAF_in_groupTest=0.0001,0.001,0.01 --groupFile=$D/groupFile.txt \
 --annotation_in_groupTest=high,medium,low,high:medium,high:medium:low --nThreads=1"
if [ "$arm" = burden ]; then ARGS="$ARGS --r.corr=1"; fi
out=$S/prof/out_$arm.txt; log=$S/prof/log_$arm.txt
{
  echo "arm=$arm start=$(date '+%F %T') load=$(uptime | sed 's/.*load averages: //')"
  echo "so_md5=$(md5 -q $REPO/.pixi/envs/default/lib/R/library/SAIGE/libs/SAIGE.so) head=$(git -C $REPO rev-parse --short HEAD)"
  SAIGE_PROFILE_REGION=1 /usr/bin/time -l $RS step2_SPAtests.R $ARGS --SAIGEOutputFile=$out
  echo "end=$(date '+%F %T')"
} > $log 2>&1
echo "$arm done"
