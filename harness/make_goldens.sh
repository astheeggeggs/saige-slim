#!/usr/bin/env bash
# Generate golden outputs for the strip/refactor harness.
#
# Run from the repo root, with the SAIGE version under test installed in the pixi env:
#   bash harness/make_goldens.sh <golden-dir>
# e.g. bash harness/make_goldens.sh harness/goldens   (baseline, at main HEAD)
#      bash harness/make_goldens.sh harness/work/candidate   (candidate, on strip branch)
#
# Design:
# - Step-2 runs consume the COMMITTED model/varianceRatio artifacts in extdata/output/
#   as frozen inputs, so step-2 outputs are byte-comparable between versions,
#   isolated from step-1 numeric drift.
# - Step-1 runs write to $GOLD/step1/ and never overwrite extdata/output/.
# - Nothing writes into extdata/output/ at all; the tree stays clean.
# - Every run sources harness/env.sh (step 1 is irreproducible without it; see LEDGER.md #7/#8).

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
GOLD="$(mkdir -p "$1" && cd "$1" && pwd)"
cd "$REPO/extdata"
source "$REPO/harness/env.sh"
mkdir -p "$GOLD/step0" "$GOLD/step1" "$GOLD/step2" "$GOLD/logs"

FAILED_RUNS=()
run() { # run <logname> <cmd...>
  local name="$1" log="$GOLD/logs/$1"; shift
  echo "=== $log"
  "$@" >"$log" 2>&1 || { echo "FAILED: see $log"; FAILED_RUNS+=("$name"); }
}
R2() { pixi run --manifest-path "$REPO/pixi.toml" Rscript "$@"; }

############################ step 0: sparse GRM #############################
# cmd_test_sparseGRM.sh, with nThreads=1 for determinism (72 in the original).
run step0_sparseGRM.log R2 createSparseGRM.R \
  --plinkFile=input/plinkforGRM_1000samples_10kMarkers \
  --nThreads=1 \
  --outputPrefix="$GOLD/step0/test_sparseGRM" \
  --numRandomMarkerforSparseKin=5000 \
  --relatednessCutoff=0.05

######################## step 1: null model fitting #########################
# cmd.sh lines 3-16: binary, sparse GRM to fit null, categorical VR.
run step1_binary_sparseGRM_a.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/step1/binary_sparseGRM_a"

# cmd.sh lines 65-80: binary, sparse GRM, explicit categorical VR MAC bins.
run step1_binary_sparseGRM_b.log R2 step1_fitNULLGLMM.R \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --useSparseGRMtoFitNULL=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2 \
  --qCovarColList=x2 \
  --sampleIDColinphenoFile=IID \
  --cateVarRatioMinMACVecExclude="0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
  --cateVarRatioMaxMACVecInclude="1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
  --traitType=binary \
  --isCateVarianceRatio=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix="$GOLD/step1/binary_sparseGRM_b"

# quantitative null model (no bundled step-1 quantitative command survives at HEAD;
# the committed example_quantitative*.rda artifacts crash step 2 — LEDGER.md #10)
run step1_quant_sparseGRM.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_quantitative \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=quantitative \
  --invNormalize=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/step1/quant_sparseGRM"

# cmd.sh line 177: survival (kept through the strip). nThreads forced 1 (was 4):
# setThreadOptions(4) would override the pinned RCPP_PARALLEL_NUM_THREADS.
# --eventTimeBinSize=1 added: without it survival step 1 crashes at HEAD
# ("object 'minEventTime' not found") — LEDGER.md #9.
run step1_survival.log R2 step1_fitNULLGLMM.R \
  --eventTimeBinSize=1 \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
  --phenoFile=./input/pheno_1000samples_survival.txt \
  --phenoCol=casecontrol \
  --covarColList=X \
  --eventTimeCol=AgeOfEventFinal \
  --sampleIDColinphenoFile=IND_ID \
  --traitType=survival \
  --outputPrefix="$GOLD/step1/survival" \
  --nThreads=1 \
  --LOCO=FALSE \
  --minMAFforGRM=0.01 \
  --skipModelFitting=FALSE \
  --tol=0.01 \
  --isCovariateOffset=FALSE \
  --IsOverwriteVarianceRatioFile=TRUE

##################### step 2: association tests (frozen inputs) #####################
SG="--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

# Frozen model inputs. Two sources, both HEAD-compatible (verified 2026-08-30):
# - committed extdata/output/example_binary_sparseGRM.{rda,varianceRatio.txt}
#   (has sparse+null VR rows, needed by group tests with a sparse GRM in step 2)
# - our own step-1 goldens (null/null_noXadj VR rows -> step 2 WITHOUT sparse args)
# The other committed artifacts (example_binary.rda, example_quantitative*.rda)
# crash step 2 at HEAD (LEDGER.md #10) and are not used.
# BGEN conditional analysis is broken at HEAD (LEDGER.md #12) -> no bgen cond runs.
# --rangestoIncludeFile is broken at HEAD (LEDGER.md #11) -> not used; the ER path
# is exercised by --minMAC=1 (ultra-rare MAC<=4 markers exist in the 100-marker set).
BIN_MODEL="--GMMATmodelFile=$GOLD/step1/binary_sparseGRM_a.rda --varianceRatioFile=$GOLD/step1/binary_sparseGRM_a.varianceRatio.txt"
QUANT_MODEL="--GMMATmodelFile=$GOLD/step1/quant_sparseGRM.rda --varianceRatioFile=$GOLD/step1/quant_sparseGRM.varianceRatio.txt"
GRP_MODEL="--GMMATmodelFile=./output/example_binary_sparseGRM.rda --varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
PLINK_G="--bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first"

# --- PLINK runs (must survive the strip byte-identically). PLINK keys markers on
# bim column 2 (rs IDs): group file is group_new_snpid.txt (verified identical to
# group_new_chrposa1a2.txt for GENE1/GENE2), conditions are rsNN.

# P1: single-variant, binary, minMAC=1 so the ER path is exercised
run step2_plink_marker_binary.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_binary.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE

# P2: group test, binary, fastTest (mirror of cmd.sh:22)
run step2_plink_group_fastTest.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_fastTest.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE

# P3: group test conditional, binary (rs30,rs79)
run step2_plink_group_cond.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_cond_binary.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --condition=rs30,rs79 --is_fastTest=TRUE

# P4: single-variant conditional, quantitative (mirror of cmd.sh:128; rs13,rs79)
run step2_plink_marker_cond.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_cond_quantitative.txt" \
  --chrom=1 --minMAF=0 --minMAC=20 --LOCO=FALSE \
  $QUANT_MODEL \
  --is_output_moreDetails=TRUE \
  --condition=rs13,rs79

# --- bgen/vcf runs (cross-format checks; dropped when the readers are stripped) ---
# B1: group test, bgen, fastTest (cmd.sh:22)
run step2_bgen_group_fastTest.log R2 step2_SPAtests.R \
  --bgenFile=./input/genotype_100markers.bgen \
  --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
  --SAIGEOutputFile="$GOLD/step2/bgen_group_fastTest.txt" \
  --chrom=1 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  --sampleFile=./input/samplelist.txt \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_chrposa1a2.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE

# B2: single-variant, bgen, binary (cross-format vs P1; bgen dosages are 8-bit
# quantized so numeric columns may differ slightly from PLINK/VCF hardcalls)
run step2_bgen_marker_binary.log R2 step2_SPAtests.R \
  --bgenFile=./input/genotype_100markers.bgen \
  --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
  --SAIGEOutputFile="$GOLD/step2/bgen_marker_binary.txt" \
  --chrom=1 --AlleleOrder=ref-first --minMAF=0 --minMAC=1 --LOCO=FALSE \
  --sampleFile=./input/samplelist.txt \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE

# V1/V2 (vcf marker runs) retired to harness/goldens_retired/ when the VCF/SAV
# reader was stripped (Strip 5). Their PLINK counterparts P1/P4 continue; the
# vcf conditional output had been verified identical to P4's.

if [ "${#FAILED_RUNS[@]}" -gt 0 ]; then
  echo "=== COMPLETE WITH FAILURES ($GOLD):"
  printf '  FAILED: %s\n' "${FAILED_RUNS[@]}"
  exit 1
fi
echo "=== all runs complete: $GOLD"
