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

# P12 (added with the LEDGER #40 fix): single-variant conditional, BINARY.
# The conditional SPA block is gated on traitType != "quantitative", so P4
# above cannot reach it and no golden covered it -- which is exactly why #40's
# uninitialised reads went unobserved. minMAC=1 keeps the ultra-rare markers
# that take the ER path (rs44, rs56), which are the ones that hit it.
run step2_plink_marker_cond_binary.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_cond_binary.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE \
  --condition=rs13,rs79

# P5-P7 (added Phase 3 commit 0): coverage for output paths that the
# de-globalization commits relocate — the BURDEN OutFile stream (r.corr=1),
# the non-fastTest sparse-GRM region path, and chunked marker append.

# P5: group test, BURDEN only (r.corr=1) — exercises openOutfile/OutFile
run step2_plink_group_burden.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE

# P6: group test, SKAT-O, is_fastTest=FALSE — exercises the org sparse-GRM path
run step2_plink_group_nofast.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_nofast.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=FALSE

# P7: single-variant, binary, multi-chunk (markers_per_chunk has a floor of
# 1000, so chunking needs a big marker set: the 128k-marker step-1 bed filtered
# to chr1 gives ~6 chunks) — exercises per-chunk append on OutFile_single and
# the index-file chunk records
run step2_plink_marker_chunked.log R2 step2_SPAtests.R \
  --bedFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bed \
  --bimFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim \
  --famFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.fam \
  --AlleleOrder=alt-first \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_chunked.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --markers_per_chunk=1000

# P8-P9/P8b (added Phase 5 C10): the multi-weight region schema. P2's
# configuration throughout, so a diff against P2 isolates the weight dimension.
# The weighted group files carry K = 2 custom sets (weight:CADD, weight,REVEL);
# rows are (anno x max_MAF x weight set) with the weight innermost.

# P8: group test, SKAT-O, K=2 custom weight sets — the Weight column plus the
# additive Pvalue_ACATV/Pvalue_ACATO columns
run step2_plink_group_multiweight.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_multiweight.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE

# P9: same, BURDEN only (r.corr=1) — the C7 Weight column in the C++ BURDEN
# writer, which is a different emitter from P8's R post-processing.
# is_fastTest=FALSE mirrors P5: BURDEN + is_fastTest=TRUE silently drops every
# region but the last (LEDGER.md #27), so that combination is characterized in
# the test suite rather than frozen into a golden.
run step2_plink_group_multiweight_burden.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_multiweight_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE

# P8b: C8 set composition (D10) on an UNWEIGHTED group file — two beta pairs
# plus the appended No_Weights set, i.e. Beta_1_25, Beta_1_1, No_Weights in
# that normative order
run step2_plink_group_composed.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_composed.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --weights.beta="1,25;1,1" --is_no_weight_in_groupTest=TRUE \
  --is_fastTest=TRUE

# P10 (added Phase 5 C11): conditional multi-weight — P3's configuration on
# the K=2 custom-weight group file (design D7: single conditioning-weight
# vector across sets; the per-set G1-side adjustment happens in R from the
# raw exported G1tilde_P_G2tilde matrix). The K=1 conditional path is
# untouched (P3 stays byte-identical).
run step2_plink_group_cond_multiweight.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_cond_multiweight.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --condition=rs30,rs79 --is_fastTest=TRUE

# P11 (added during the AF>0.5 investigation; P10 above is the C11
# conditional multi-weight golden): the SAME bed as P1, read
# with
# --AlleleOrder=ref-first, which makes every marker AF_Allele2 > 0.5 (a
# mirror-image dataset, since the flip is exact -- see the invariance control in
# tests/testthat/test-step2-allele-order.R). It pins the two AF>0.5 bugs at
# once, and is a CHARACTERIZATION golden: it captures buggy output on purpose.
#   LEDGER #28 - variance mis-centred by scoreTestFast_noadjCov (var inflated
#                ~178x median, p driven to 1)
#   LEDGER #29 - 7 of the 70 markers (the MAC=1 singletons) silently dropped by
#                the gate-1 MAC round-trip, hence 63 data rows not 70
# When #28/#29 are fixed this file MUST change and every other golden MUST NOT.
run step2_plink_marker_reffirst.log R2 step2_SPAtests.R \
  --bedFile=./input/genotype_100markers.bed \
  --bimFile=./input/genotype_100markers.bim \
  --famFile=./input/genotype_100markers.fam \
  --AlleleOrder=ref-first \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_reffirst.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE

# P13-P15 (Phase 7): PGEN (plink2 native format) full-format-parity goldens.
# genotype_100markers.{pgen,pvar,psam} was generated once via
# `plink2 --bfile genotype_100markers --make-pgen --out genotype_100markers`
# (deterministic, committed like the bed/bim/fam fixture it mirrors). pvar
# REF/ALT for every marker matches the alt-first bim's REF/ALT exactly (both
# label the same allele ALT), so PGEN and PLINK runs over the SAME phenotype/
# covariate/model inputs should reconcile at float tolerance, not just look
# structurally similar -- see the numeric reconciliation check in
# tests/testthat/test-step2-pgen.R. Same model inputs, chrom/MAF/MAC cutoffs
# as the PLINK runs they mirror (P1/P2/P3) so a diff isolates the reader.
PGEN_G="--pgenFile=./input/genotype_100markers.pgen --pvarFile=./input/genotype_100markers.pvar --psamFile=./input/genotype_100markers.psam"

# P13: single-variant, binary, minMAC=1 (mirrors P1; ER path exercised)
run step2_pgen_marker_binary.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_marker_binary.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE

# P14: group test, binary, fastTest (mirrors P2)
run step2_pgen_group_fastTest.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_fastTest.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE

# P15: group test conditional, binary (mirrors P3; rs30,rs79) -- exercises
# region + conditional logic together through the pgen reader in one golden.
run step2_pgen_group_cond.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_cond_binary.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --condition=rs30,rs79 --is_fastTest=TRUE

# P16-P24 (Phase 7 coverage extension): PGEN mirrors of every DISTINCT code
# path the PLINK goldens exercise beyond P1-P3 (P13-P15's mirrors), except
# P11 (--AlleleOrder=ref-first is a PLINK bim-specific toggle; pvar carries
# fixed REF/ALT columns, so there is no ref-first/alt-first ambiguity to
# mirror). Same model/phenotype/sparse-GRM inputs as each PLINK counterpart,
# so a diff isolates the reader -- exactly the P13-P15 pattern.

# P16: single-variant conditional, QUANTITATIVE (mirrors P4) -- the pgen path
# had ZERO quantitative-trait coverage before this; QUANT_MODEL is fit on the
# same step-1 quant golden P4 uses.
run step2_pgen_marker_cond_quant.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_marker_cond_quantitative.txt" \
  --chrom=1 --minMAF=0 --minMAC=20 --LOCO=FALSE \
  $QUANT_MODEL \
  --is_output_moreDetails=TRUE \
  --condition=rs13,rs79

# P17: single-variant conditional, BINARY, minMAC=1 (mirrors P12) -- the
# LEDGER #40 conditional-SPA fix path (reached via the ER path); confirms the
# pgen reader hands the ER-path markers the same decoded values as PLINK.
run step2_pgen_marker_cond_binary.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_marker_cond_binary.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE \
  --condition=rs13,rs79

# P18: group test, BURDEN only, r.corr=1, is_fastTest=FALSE (mirrors P5) --
# the BURDEN OutFile writer through the pgen reader.
run step2_pgen_group_burden.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE

# P19: group test, SKAT-O, is_fastTest=FALSE (mirrors P6) -- the non-fastTest
# sparse-GRM region path through the pgen reader.
run step2_pgen_group_nofast.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_nofast.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=FALSE

# P20: single-variant, multi-chunk (mirrors P7) -- uses a DIFFERENT, larger
# fixture (the 128k-marker step-1 bed filtered to chr1), so this also
# regenerates a second committed PGEN fixture:
# nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.{pgen,pvar,psam},
# generated once via `plink2 --bfile
# nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr --make-pgen --out
# nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr` the same way as
# the 100-marker fixture. Exercises per-chunk append through the pgen reader.
PGEN_CHUNK_G="--pgenFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.pgen --pvarFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.pvar --psamFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.psam"
run step2_pgen_marker_chunked.log R2 step2_SPAtests.R \
  $PGEN_CHUNK_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_marker_chunked.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --markers_per_chunk=1000

# P21: group test, SKAT-O, K=2 custom weight sets, is_fastTest=TRUE (mirrors
# P8) -- the multi-weight Weight/Pvalue_ACATV/Pvalue_ACATO schema through the
# pgen reader.
run step2_pgen_group_multiweight.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_multiweight.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE

# P22: group test, multiweight BURDEN, r.corr=1, is_fastTest=FALSE (mirrors
# P9) -- the C7 Weight column in the C++ BURDEN writer through the pgen
# reader.
run step2_pgen_group_multiweight_burden.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_multiweight_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE

# P23: group test, "composed" set -- beta weights + appended No_Weights
# (mirrors P8b) -- the C8 set-composition path (D10) through the pgen reader.
run step2_pgen_group_composed.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_composed.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --weights.beta="1,25;1,1" --is_no_weight_in_groupTest=TRUE \
  --is_fastTest=TRUE

# P24: group test, conditional multi-weight (mirrors P10) -- design D7's
# per-set G1-side adjustment through the pgen reader.
run step2_pgen_group_cond_multiweight.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_cond_multiweight.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --condition=rs30,rs79 --is_fastTest=TRUE

# B1/B2 (bgen runs) and V1/V2 (vcf runs) retired to harness/goldens_retired/
# when the BGEN (Strip 6) and VCF/SAV (Strip 5) readers were stripped. Their
# PLINK counterparts P1-P4 continue; before retirement, the bgen group
# fastTest run had byte-stable goldens and the vcf conditional output was
# verified identical to P4's.

if [ "${#FAILED_RUNS[@]}" -gt 0 ]; then
  echo "=== COMPLETE WITH FAILURES ($GOLD):"
  printf '  FAILED: %s\n' "${FAILED_RUNS[@]}"
  exit 1
fi
echo "=== all runs complete: $GOLD"
