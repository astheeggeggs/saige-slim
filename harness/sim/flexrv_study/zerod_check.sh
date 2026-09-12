#!/bin/bash
# LEDGER #129 check: GENE0181 (uninformative score arm, y_binary) run four ways --
# collapse 10 and 0, each with the dosage zeroing on (default 0.2) and off (0).
# With the zeroing off the collapse-10 p of every set equals collapse 0's.
# usage: ZEROD_DIR=<scratch> bash zerod_check.sh
S=${ZEROD_DIR:?set ZEROD_DIR to a scratch directory}
REPO=/Users/dpalmer/Repositories/SAIGE; D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
cd $REPO; source harness/env.sh; RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
grep '^GENE0181 ' harness/work/flexrv_study/scores/groupFile_uninformative.txt > $S/group.txt
cd extdata
COMMON="--bedFile=$D/chr20_sim.bed --bimFile=$D/chr20_sim.bim --famFile=$D/chr20_sim.fam --chrom=20 --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=FALSE --nThreads=1 --GMMATmodelFile=$D/step1_y_binary.rda --varianceRatioFile=$D/step1_y_binary.varianceRatio.txt --sparseGRMFile=$D/pedigree_truth.sparseGRM.mtx --sparseGRMSampleIDFile=$D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt --groupFile=$S/group.txt --annotation_in_groupTest=high:medium --flexRV_score=SIM --flexRV_maxMAF=0.001 --flexRV_lofAnno=high --r.corr=1"
for cfg in "c10_zerod0.2:--MACCutoff_to_CollapseUltraRare=10" "c10_zerod0:--MACCutoff_to_CollapseUltraRare=10 --dosage_zerod_cutoff=0" "c0_zerod0.2:--MACCutoff_to_CollapseUltraRare=0" "c0_zerod0:--MACCutoff_to_CollapseUltraRare=0 --dosage_zerod_cutoff=0"; do
  tag=${cfg%%:*}; extra=${cfg#*:}
  $RS step2_SPAtests.R $COMMON $extra --SAIGEOutputFile=$S/out_$tag.txt > $S/log_$tag.txt 2>&1; rc=$?
  echo "== $tag exit=$rc: $(awk -F'\t' '$4=="SIM_sbeta10_MAC10"{print $5, $6, $7, $8, $9, $10}' $S/out_$tag.txt)  | rows $(awk 'NR>1' $S/out_$tag.txt | wc -l | tr -d ' ') | pooled $(awk -F'\t' '$2=="Cauchy"{print $5}' $S/out_$tag.txt)"
done
