#!/bin/bash
# LEDGER #157 local reproduction: does a bed SUBSET perturb region cells while
# the PGEN of the SAME markers reproduces the full bed exactly?
set -uo pipefail
REPO=/Users/dpalmer/Repositories/SAIGE
D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
W=$REPO/harness/work/l157
cd $REPO; source harness/env.sh; cd extdata
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
COMMON="--chrom=20 --minMAF=0 --minMAC=0.5 \
 --GMMATmodelFile=$D/step1_y_binary.rda --varianceRatioFile=$D/step1_y_binary.varianceRatio.txt \
 --sparseGRMFile=$D/pedigree_truth.sparseGRM.mtx \
 --sparseGRMSampleIDFile=$D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt \
 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_output_moreDetails=TRUE \
 --is_fastTest=FALSE --is_single_in_groupTest=TRUE \
 --maxMAF_in_groupTest=0.0001,0.001,0.01 --groupFile=$D/groupFile.txt \
 --annotation_in_groupTest=high,medium,low,high:medium,high:medium:low --nThreads=1"
run() { echo "$(date '+%T') START $1"; $RS step2_SPAtests.R $2 $COMMON --SAIGEOutputFile=$W/$1.txt > $W/$1.log 2>&1; echo "$(date '+%T') DONE $1 rc=$?"; }
run fullbed "--bedFile=$D/chr20_sim.bed --bimFile=$D/chr20_sim.bim --famFile=$D/chr20_sim.fam"
run subbed  "--bedFile=$W/sub_bed.bed --bimFile=$W/sub_bed.bim --famFile=$W/sub_bed.fam"
run subpgen "--pgenFile=$W/sub_pgen.pgen --pvarFile=$W/sub_pgen.pvar --psamFile=$W/sub_pgen.psam"
echo ALLDONE
