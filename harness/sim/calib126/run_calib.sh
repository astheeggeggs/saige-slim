#!/bin/bash
# LEDGER #126 null-calibration driver (docs/state/prediction-126-null-calibration.md).
# usage: CALIB_DIR=<dir> [SAIGE_LIB=<lib>] run_calib.sh <trait> <arm> <gate>
#   trait  y_null | y_null_bin | y_binary      gate  the --minGroupMAC_in_BurdenTest value
#   arm    plain | flex | grid  (gene-level: fixture groupFile.txt / flexrv_study informative scores)
#          wplain | wgrid | wgridr (windowed files from make_windowed.py in CALIB_DIR)
# Records head, dirty flag, library and SAIGE.so md5 in the log; exit status captured
# BEFORE the $(date) that would otherwise reset it.
set -o pipefail
S=${CALIB_DIR:?set CALIB_DIR to the working directory (group files, runs/)}
REPO=/Users/dpalmer/Repositories/SAIGE
D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
trait=$1; arm=$2; gate=$3
export R_LIBS="${SAIGE_LIB:-$REPO/.pixi/envs/default/lib/R/library}"
cd $REPO; source harness/env.sh
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
OUT=$S/runs; mkdir -p $OUT
COMMON="--bedFile=$D/chr20_sim.bed --bimFile=$D/chr20_sim.bim --famFile=$D/chr20_sim.fam \
 --chrom=20 --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=FALSE --nThreads=1 \
 --GMMATmodelFile=$D/step1_$trait.rda --varianceRatioFile=$D/step1_$trait.varianceRatio.txt \
 --sparseGRMFile=$D/pedigree_truth.sparseGRM.mtx \
 --sparseGRMSampleIDFile=$D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt \
 --minGroupMAC_in_BurdenTest=$gate"
GF=$D/groupFile.txt; GFLEX=$REPO/harness/work/flexrv_study/scores/groupFile_informative.txt; GW=$S/groupFile_w25.txt
case $arm in
  plain) ARGS="$COMMON --groupFile=$GF --annotation_in_groupTest=high,high:medium,high:medium:low --maxMAF_in_groupTest=0.0001,0.001,0.01 --is_no_weight_in_groupTest=TRUE --r.corr=1" ;;
  flex)  ARGS="$COMMON --groupFile=$GFLEX --annotation_in_groupTest=high:medium --flexRV_score=SIM --flexRV_maxMAF=0.001 --flexRV_lofAnno=high --r.corr=1" ;;
  grid)  ARGS="$COMMON --groupFile=$GF --annotation_in_groupTest=high,high:medium,high:medium:low --maxMAF_in_groupTest=0.0001,0.001,0.01" ;;
  wplain) ARGS="$COMMON --groupFile=$GW --annotation_in_groupTest=high,high:medium --maxMAF_in_groupTest=0.0001,0.001 --is_no_weight_in_groupTest=TRUE --r.corr=1" ;;
  wgrid)  ARGS="$COMMON --groupFile=$GW --annotation_in_groupTest=high,high:medium --maxMAF_in_groupTest=0.0001,0.001" ;;
  wgridr) ARGS="$COMMON --groupFile=$S/groupFile_w25r.txt --annotation_in_groupTest=high,high:medium --maxMAF_in_groupTest=0.0001,0.001" ;;
  *) echo "unknown arm $arm" >&2; exit 2 ;;
esac
tag=${trait}_${arm}_g${gate}
out=$OUT/$tag.txt; log=$OUT/$tag.log
{
  echo "tag=$tag start=$(date '+%F %T') load=$(uptime | sed 's/.*load averages: //')"
  echo "head=$(git rev-parse --short HEAD) dirty=$(git status --porcelain | grep -c .) lib=$R_LIBS so_md5=$(md5 -q $R_LIBS/SAIGE/libs/SAIGE.so)"
  cd extdata
  /usr/bin/time -l $RS step2_SPAtests.R $ARGS --SAIGEOutputFile=$out; rc=$?
  echo "end=$(date "+%F %T") exit=$rc"
} > $log 2>&1
n=$(awk 'NR>1' $out 2>/dev/null | wc -l | tr -d ' ')
echo "$tag: $n data rows, exit=$(grep -o 'exit=[0-9]*' $log | tail -1), regions=$(grep -c 'Analyzing Region' $log), $(grep -E ' real ' $log | awk '{print $3 " s user"}')"
