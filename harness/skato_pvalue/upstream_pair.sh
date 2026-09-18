#!/bin/bash
# Paired SKAT-O timing, upstream main (SAIGE-main, its own pixi R and Rlib) vs strip HEAD, pgen reader, y_binary, 1 thread.
set -uo pipefail
S=$1; REPO=/Users/dpalmer/Repositories/SAIGE; MAIN=/Users/dpalmer/Repositories/SAIGE-main; D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
cd $REPO; source harness/env.sh
ARGS="--chrom=20 --minMAF=0 --minMAC=0.5 --GMMATmodelFile=$D/step1_y_binary.rda --varianceRatioFile=$D/step1_y_binary.varianceRatio.txt \
 --sparseGRMFile=$D/pedigree_truth.sparseGRM.mtx --sparseGRMSampleIDFile=$D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt \
 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_output_moreDetails=TRUE --is_fastTest=TRUE \
 --is_output_markerList_in_groupTest=TRUE --is_single_in_groupTest=TRUE --maxMAF_in_groupTest=0.0001,0.001,0.01 \
 --groupFile=$D/groupFile.txt --annotation_in_groupTest=high,medium,low,high:medium,high:medium:low --nThreads=1"
for arm in main strip main strip; do
  if [ $arm = main ]; then RS=$MAIN/.pixi/envs/default/lib/R/bin/Rscript; export R_LIBS_USER=$MAIN/Rlib; SCRIPT=$MAIN/extdata/step2_SPAtests.R; else RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript; unset R_LIBS_USER; SCRIPT=$REPO/extdata/step2_SPAtests.R; fi
  if [ $arm = main ]; then GENO="--pgenPrefix=$D/chr20_sim_pgen --AlleleOrder=ref-first"; else GENO="--pgenFile=$D/chr20_sim_pgen.pgen --pvarFile=$D/chr20_sim_pgen.pvar --psamFile=$D/chr20_sim_pgen.psam"; fi
  out=$S/upcmp/out_${arm}_$RANDOM.txt
  { echo "arm=$arm start=$(date '+%F %T') load=$(uptime | sed 's/.*load averages: //')"
    $RS -e 'cat("lib=", find.package("SAIGE"), " version=", as.character(packageVersion("SAIGE")), "\n")'
    cd $(dirname $SCRIPT); /usr/bin/time -l $RS $SCRIPT $GENO $ARGS --SAIGEOutputFile=$out
    echo "end=$(date '+%F %T')"; } >> $S/upcmp/log_$arm.txt 2>&1
  echo "$arm: $(grep -E 'real|user' $S/upcmp/log_$arm.txt | tail -2 | tr '\n' ' ') rows=$(wc -l < $out 2>/dev/null)"
done
