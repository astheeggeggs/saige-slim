#!/bin/bash
# Null replicate study for the SKAT-O p-value options (2026-09-18): 32 replicate nulls x
# {quantitative, binary K=0.01, binary K=0.002}; per run the SKAT-O cells are dumped and scored offline.
set -uo pipefail
REPO=/Users/dpalmer/Repositories/SAIGE; W=$REPO/harness/work/nullrep-20260918; D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
cd $REPO; . harness/sim/00_env.sh; export MPLBACKEND=Agg
GRM=$D/pedigree_truth.sparseGRM.mtx; PH=$W/pheno_pcs_nullrep.txt
log() { echo "$(date '+%F %T') load=$(uptime | sed 's/.*load averages: //') $*" >> $W/logs/driver.log; }
fit1() { # col tt
  local col=$1 tt=$2 extra=""; [ "$tt" = quantitative ] && extra="--invNormalize=TRUE"
  [ -s $W/models/step1_$col.varianceRatio.txt ] && [ -s $W/models/step1_$col.rda ] && return 0
  local t0=$(date +%s)
  $RS extdata/step1_fitNULLGLMM.R --plinkFile=$D/data_for_vr --sparseGRMFile=$GRM --sparseGRMSampleIDFile=$GRM.sampleIDs.txt \
    --useSparseGRMtoFitNULL=TRUE --useSparseGRMforVarRatio=TRUE --phenoFile=$PH --phenoCol=$col --sampleIDColinphenoFile=IID \
    --traitType=$tt $extra --covarColList=age,sex,PC1,PC2,PC3,PC4 --qCovarColList=sex --skipVarianceRatioEstimation=FALSE \
    --isCateVarianceRatio=TRUE --relatednessCutoff=0.125 --IsOverwriteVarianceRatioFile=TRUE --nThreads=2 \
    --outputPrefix=$W/models/step1_$col > $W/logs/s1_$col.log 2>&1
  echo "$(date '+%F %T') fit $col rc=$? $(( $(date +%s) - t0 ))s" >> $W/logs/driver.log
}
run2() { # col modelprefix
  local col=$1 M=$2; local out=$W/runs/$col.txt
  [ -s $W/runs/$col.cells.rds ] && return 0
  local t0=$(date +%s)
  ( cd extdata; $RS $W/dump_arm.R $out $W/runs/$col.cells.rds $D $M > $W/logs/s2_$col.log 2>&1 )
  echo "$(date '+%F %T') step2 $col rc=$? $(( $(date +%s) - t0 ))s" >> $W/logs/driver.log
}
score() { local col=$1; [ -s $W/runs/$col.variants.rds ] && return 0; mkdir -p $W/runs/$col.d; ln -sf $W/runs/$col.cells.rds $W/runs/$col.d/cells.rds
  local t0=$(date +%s); $RS $W/variants_null.R $W/runs $col.d > $W/logs/sc_$col.log 2>&1; mv -f $W/runs/$col.d/variants.rds $W/runs/$col.variants.rds 2>/dev/null
  echo "$(date '+%F %T') score $col rc=$? $(( $(date +%s) - t0 ))s" >> $W/logs/driver.log; }
export -f fit1 run2 score; export RS W D GRM PH REPO
log "STAGE 1: step-1 fits"
{ for r in $(seq 6 32); do echo "y_null_r$r quantitative"; done; for r in $(seq 1 32); do echo "y_nb01_r$r binary"; echo "y_nb002_r$r binary"; done; } | xargs -P 4 -L 1 bash -c 'fit1 "$@"' _
log "STAGE 2: step-2 dump runs"
{ for r in $(seq 1 32); do if [ $r -le 5 ]; then echo "y_null_r$r $D/step1_y_null_r$r"; else echo "y_null_r$r $W/models/step1_y_null_r$r"; fi; echo "y_nb01_r$r $W/models/step1_y_nb01_r$r"; echo "y_nb002_r$r $W/models/step1_y_nb002_r$r"; done; } | xargs -P 6 -L 1 bash -c 'run2 "$@"' _
log "STAGE 3: offline scoring"
{ for r in $(seq 1 32); do echo y_null_r$r; echo y_nb01_r$r; echo y_nb002_r$r; done; } | xargs -P 6 -L 1 bash -c 'score "$@"' _
log "DONE"
