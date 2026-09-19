#!/bin/bash
# Does splitting the group file across PROCESSES beat threading within one
# process? BURDEN only -- SKAT-O already got its cross-region win (step 5).
set -uo pipefail
REPO=/Users/dpalmer/Repositories/SAIGE
D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
W=$REPO/harness/work/regthread
cd $REPO; source harness/env.sh; cd extdata
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
GENO="--pgenFile=$D/chr20_sim_pgen.pgen --pvarFile=$D/chr20_sim_pgen.pvar --psamFile=$D/chr20_sim_pgen.psam"
COMMON="$GENO --chrom=20 --minMAF=0 --minMAC=0.5 \
 --GMMATmodelFile=$D/step1_y_binary.rda --varianceRatioFile=$D/step1_y_binary.varianceRatio.txt \
 --sparseGRMFile=$D/pedigree_truth.sparseGRM.mtx \
 --sparseGRMSampleIDFile=$D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt \
 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_output_moreDetails=TRUE \
 --is_fastTest=FALSE --is_single_in_groupTest=TRUE --r.corr=1 \
 --maxMAF_in_groupTest=0.0001,0.001,0.01 \
 --annotation_in_groupTest=high,medium,low,high:medium,high:medium:low"
t() { date +%s.%N; }
for rep in 1 2; do
  s=$(t); /usr/bin/time -l $RS step2_SPAtests.R $COMMON --groupFile=$D/groupFile.txt --nThreads=1 \
      --SAIGEOutputFile=$W/serial_r$rep.txt > $W/serial_r$rep.log 2>&1; e=$(t)
  echo "rep$rep serial_1thread   $(echo "$e - $s" | bc)  rss=$(grep 'maximum resident' $W/serial_r$rep.log | awk '{printf "%.0f MB", $1/1048576}')"
  s=$(t); /usr/bin/time -l $RS step2_SPAtests.R $COMMON --groupFile=$D/groupFile.txt --nThreads=4 \
      --SAIGEOutputFile=$W/threads_r$rep.txt > $W/threads_r$rep.log 2>&1; e=$(t)
  echo "rep$rep threads_4        $(echo "$e - $s" | bc)  rss=$(grep 'maximum resident' $W/threads_r$rep.log | awk '{printf "%.0f MB", $1/1048576}')"
  s=$(t)
  for q in 1 2 3 4; do
    /usr/bin/time -l $RS step2_SPAtests.R $COMMON --groupFile=$W/grp_q$q.txt --nThreads=1 \
        --SAIGEOutputFile=$W/proc${q}_r$rep.txt > $W/proc${q}_r$rep.log 2>&1 &
  done
  wait; e=$(t)
  echo "rep$rep procs_4x1thread  $(echo "$e - $s" | bc)  rss_sum=$(cat $W/proc?_r$rep.log | grep 'maximum resident' | awk '{s+=$1} END {printf "%.0f MB", s/1048576}')"
done
echo ALLDONE
