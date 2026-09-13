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
# - LEDGER #127: a group run also writes `<out>.pooledTests.txt` when a region's
#   pooled (Cauchy) row was NOT a plain Cauchy combination of every cell it
#   tested. Five runs produce one -- the quantitative SKAT-O goldens, whose
#   saturated cells (Pvalue exactly 1) send CCT down its min(1, n*min(p))
#   fallback, which is LEDGER #86. They are committed like any other output, so
#   compare_goldens.sh gates the disclosure; every other run writes no file, and
#   a new one appearing means a region started combining something differently.
#   (P64 wrote a sixth for one commit, when its only cell vanished to the bug
#   LEDGER #129 records; fixing that took the file away again, which is the
#   disclosure doing its job in both directions.)

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

# Same run through the PGEN reader. extdata/input/plinkforGRM_1000samples_10kMarkers.
# {pgen,pvar,psam} were generated once from the bed/bim/fam with
# `plink2 --bfile <prefix> --make-pgen --out <prefix>` and committed; verified at
# generation time that pvar ALT == bim column 5 and pvar REF == bim column 6 for
# all 10,000 markers, and that psam IID order matches the fam, which is what makes
# the two readers comparable at all.
#
# The output is BYTE-IDENTICAL to the PLINK run above. That is a real gate here
# rather than a claim about two code paths agreeing: the pgen front-end re-encodes
# each variant as a PLINK-1 .bed record and hands it to the SAME
# Get_OneSNP_Geno_atBeginning, so the QC/impute/repack pipeline and every float
# downstream of it is one code path with two front-ends.
# NOTE the seed: extdata/createSparseGRM.R does set.seed(1) before sampling the
# 5,000 markers out of the 9,588 that pass MAF. WITHOUT that seed neither run
# reproduces itself, so a comparison of two unseeded runs measures the RNG, not
# the reader.
run step0_sparseGRM_pgen.log R2 createSparseGRM.R \
  --plink2File=input/plinkforGRM_1000samples_10kMarkers \
  --nThreads=1 \
  --outputPrefix="$GOLD/step0/test_sparseGRM_pgen" \
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

# Binary, FULL (dense) GRM -- the only step-1 golden that is not a sparse-GRM
# fit. Without it, getCrossprodMatAndKin always takes its sparse branch and the
# dense TBB kernel (parallelCrossProd / GetTrace's Au) is exercised by NO
# golden; a change to it would leave every golden byte-identical while being
# completely untested. Config mirrors the s1_full bench workload so bench rows
# and golden runs describe the same code path. ~50 s.
run step1_binary_denseGRM.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --nThreads=1 \
  --LOCO=FALSE \
  --minMAFforGRM=0.01 \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix="$GOLD/step1/binary_denseGRM"

# Binary, sparse GRM to fit the null, PLUS --useSparseGRMforVarRatio=TRUE. The
# ONLY reason this exists: it is the only configuration that emits `sparse`
# variance-ratio rows (R/SAIGE_fitGLMM_fast.R:2467), and the step-2 MARKER path
# needs sparse VR rows AND null_noXadj rows in the same file to take its
# sparse-GRM branch. The committed example_binary_sparseGRM.varianceRatio.txt
# has sparse+null but no null_noXadj, so single-variant tests refuse it
# (SAIGE_Test_main.R:400, the LEDGER #10 guard); our other step-1 goldens have
# null+null_noXadj but no sparse. Hence P27-P29 below.
# Config is binary_sparseGRM_a's with the one flag added, and that is load-bearing:
# the resulting .rda is BYTE-IDENTICAL to binary_sparseGRM_a.rda and the null /
# null_noXadj rows are identical too (the flag only appends sparse rows), so a
# diff of P27-P29 against P1 isolates the step-2 branch and not the model.
# tests/testthat/test-invariants-goldens.R pins that identity.
run step1_binary_sparseGRM_vr.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --useSparseGRMforVarRatio=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/step1/binary_sparseGRM_vr"

# quant_sparseGRM_vr: the QUANTITATIVE counterpart of binary_sparseGRM_vr, and
# the reason the quantitative sparse-GRM step-2 path had no golden at all until
# 2026-09-08. A region test needs `sparse` rows in the variance-ratio file and
# only --useSparseGRMforVarRatio=TRUE writes them; quant_sparseGRM above is a
# --useSparseGRMtoFitNULL fit, so its file has null/null_noXadj only and step 2
# must NOT be given --sparseGRMFile for it (which is exactly why P36/P37 do
# not). Config is quant_sparseGRM's with the one flag added, so -- as for the
# binary pair -- the .rda should come out BYTE-IDENTICAL to quant_sparseGRM.rda
# and a diff of the runs below against P36/P37 isolates the step-2 branch
# rather than the model.
#
# It also puts a quantitative trait through LEDGER #97: #65 rescales the
# null/null_noXadj ratios by 1/s^2 on the way out and leaves the `sparse` ratio
# alone because it is scale-free. Before this config no golden contained a
# quantitative `sparse` row at all.
run step1_quant_sparseGRM_vr.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --useSparseGRMforVarRatio=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_quantitative \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=quantitative \
  --invNormalize=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/step1/quant_sparseGRM_vr"


################## step 1 mirrors through the PGEN reader ##################
# The same five configurations, reading
# extdata/input/nfam_100_..._22chr.{pgen,pvar,psam} (already committed for the
# step-2 PGEN goldens) instead of the bed/bim/fam.
#
# Each is BYTE-IDENTICAL to its PLINK counterpart above -- both the .rda and the
# .varianceRatio.txt -- and tests/testthat/test-step1-pgen.R asserts exactly
# that with cmp. This is not the forbidden "two code paths agree" comparison:
# the pgen front-end re-encodes each variant as a PLINK-1 .bed record and hands
# it to the SAME Get_OneSNP_Geno_atBeginning, so QC, mean imputation, the 2-bit
# repack and every float after them are literally one code path. If this ever
# needs a tolerance, the design has broken -- do not relax the assertion.
#
# binary_denseGRM_pgen is kept for the same reason its PLINK twin is: it is the
# only one of the five on the DENSE GRM path, so without it the TBB kernel is
# reached by no pgen golden at all.

run step1_binary_sparseGRM_a_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
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
  --outputPrefix="$GOLD/step1/binary_sparseGRM_a_pgen"

run step1_binary_sparseGRM_b_pgen.log R2 step1_fitNULLGLMM.R \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --plink2File=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
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
  --outputPrefix="$GOLD/step1/binary_sparseGRM_b_pgen"

run step1_quant_sparseGRM_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
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
  --outputPrefix="$GOLD/step1/quant_sparseGRM_pgen"

run step1_binary_denseGRM_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --nThreads=1 \
  --LOCO=FALSE \
  --minMAFforGRM=0.01 \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix="$GOLD/step1/binary_denseGRM_pgen"

run step1_binary_sparseGRM_vr_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --useSparseGRMforVarRatio=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/step1/binary_sparseGRM_vr_pgen"

run step1_quant_sparseGRM_vr_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --useSparseGRMforVarRatio=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_quantitative \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=quantitative \
  --invNormalize=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/step1/quant_sparseGRM_vr_pgen"

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
############################ step 1: LOCO ###################################
# LOCO had ZERO golden coverage: `--LOCO=TRUE` appeared in 0 of the 51 golden
# runs and `--LOCO=FALSE` in 46 of them. It is three mechanisms (CLAUDE.md) --
# step-1 global-state mutation per chromosome, the save-time file layout, and
# step-2 `ReadModel` overwriting fitted values from `LOCOResult[[chrom]]` -- and
# the third had no end-to-end gate at all. `test-step1-lowmem-loco.R` covers the
# second and `test-unit-readmodel.R` the third at unit level.
#
# ONLY a dense-GRM fit can produce a LOCO model: `useSparseGRMtoFitNULL=TRUE`
# forces `LOCO <- FALSE` at R/SAIGE_fitGLMM_fast.R:836, so four of the five
# step-1 goldens can never carry a LOCOResult.
#
# The fixture is plinkforGRM_1000samples_10kMarkers, not the 22-chromosome one
# the other dense golden uses, and the choice was measured rather than assumed:
#   - 22 chromosomes: 47 s, a 4.7 MB .rda (LOCOResult holds 22 copies of every
#     per-sample vector, against ~268 KB for every other step-1 golden), and
#     chr1-vs-chr2 residuals differing by only 1.5e-3.
#   - this one:        7 s, a 677 KB .rda, and chr1-vs-chr2 residuals differing
#     by 0.19 -- 130x more signal, because it is 9,988 markers on chr1 against
#     12 on chr2, so leaving out chr1 changes the GRM enormously and leaving out
#     chr2 barely at all.
# LOCOResult keeps its length-22 shape either way; here 20 entries are
# `list(isLOCO = FALSE)` stubs. It is also the first step-1 golden on this
# fixture, which until now fed step 0 only.
#
# NOTE the harness had to be fixed to gate this at all: compare_step1.R
# descended into a list only when the field was named `obj.noK`, so LOCOResult
# was skipped SILENTLY and this golden would have compared theta and the base
# fields while checking not one value of the thing it exists for.
run step1_binary_loco_2chr.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --nThreads=1 \
  --LOCO=TRUE \
  --minMAFforGRM=0.01 \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix="$GOLD/step1/binary_loco_2chr"

# A LOCO model that can also drive a REGION test. A group test given
# --sparseGRMFile needs `sparse` variance-ratio rows, and the only model in the
# repo that had them (extdata/output/example_binary_sparseGRM) is LOCO=FALSE,
# while binary_loco_2chr above carries only null/null_noXadj -- so no model could
# do both until this one. The combination that works is
# `--useSparseGRMtoFitNULL=FALSE` (a DENSE fit, so LOCO survives) with
# `--useSparseGRMforVarRatio=TRUE` (which supplies the sparse rows): the two
# flags are independent, and only the FIRST forces LOCO off
# (R/SAIGE_fitGLMM_fast.R:836). 20 s, and the variance-ratio file comes out with
# all three row types over two MAC categories.
run step1_binary_loco_vr_2chr.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/plinkforGRM_1000samples_10kMarkers \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=FALSE \
  --useSparseGRMforVarRatio=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --nThreads=1 \
  --LOCO=TRUE \
  --minMAFforGRM=0.01 \
  --isCateVarianceRatio=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix="$GOLD/step1/binary_loco_vr_2chr"

# PGEN mirror of it; byte-identity asserted in the candidate dir, not committed.
run step1_binary_loco_vr_2chr_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/plinkforGRM_1000samples_10kMarkers \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=FALSE \
  --useSparseGRMforVarRatio=TRUE \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --nThreads=1 \
  --LOCO=TRUE \
  --minMAFforGRM=0.01 \
  --isCateVarianceRatio=TRUE \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix="$GOLD/step1/binary_loco_vr_2chr_pgen"

# PGEN mirror of the LOCO fit; byte-identity is asserted in the candidate dir
# below rather than committed, as for every other step-1 pgen run.
run step1_binary_loco_2chr_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/plinkforGRM_1000samples_10kMarkers \
  --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --nThreads=1 \
  --LOCO=TRUE \
  --minMAFforGRM=0.01 \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outputPrefix="$GOLD/step1/binary_loco_2chr_pgen"

BIN_MODEL="--GMMATmodelFile=$GOLD/step1/binary_sparseGRM_a.rda --varianceRatioFile=$GOLD/step1/binary_sparseGRM_a.varianceRatio.txt"
QUANT_MODEL="--GMMATmodelFile=$GOLD/step1/quant_sparseGRM.rda --varianceRatioFile=$GOLD/step1/quant_sparseGRM.varianceRatio.txt"
VR_MODEL="--GMMATmodelFile=$GOLD/step1/binary_sparseGRM_vr.rda --varianceRatioFile=$GOLD/step1/binary_sparseGRM_vr.varianceRatio.txt"
QUANT_VR_MODEL="--GMMATmodelFile=$GOLD/step1/quant_sparseGRM_vr.rda --varianceRatioFile=$GOLD/step1/quant_sparseGRM_vr.varianceRatio.txt"
GRP_MODEL="--GMMATmodelFile=./output/example_binary_sparseGRM.rda --varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
PLINK_G="--bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first"

# --- PLINK runs (must survive the strip byte-identically). PLINK keys markers on
# bim column 2 (rs IDs): group file is group_new_snpid.txt (verified identical to
# group_new_chrposa1a2.txt for GENE1/GENE2, now in
# harness/goldens_retired/inputs/ -- LEDGER #113), conditions are rsNN.

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
# --is_output_markerList_in_groupTest=TRUE on the six goldens below (LEDGER
# #132, 2026-09-13). It writes a `.markerList.txt` sidecar: one row per
# (region, annotation, max_MAF) with the surviving RARE rsIDs and the folded
# ULTRA-RARE rsIDs (R/SAIGE_SPATest_Region.R:1497-1556, weight-independent
# because the indicator marks a marker folded regardless of its weight).
#
# It is here because that membership is the one input
# tests/testthat/test-invariants-goldens.R could not reconstruct once #132 made
# the pseudo-marker's ACAT-V weight the MEAN of the folded markers' non-zero
# weights: folded membership is per cell (annotation AND max_MAF AND below the
# collapse MAC cutoff) and a folded marker has no single-assoc row, so its MAC
# is in no committed file. With the sidecar the test recomputes ACAT-V from
# first principles again -- members from here, weights from the group file,
# per-set non-zero filter applied by the test -- instead of solving for a free
# parameter and bounding it.
#
# The flag gates only the block at :1497, so it adds a file and changes nothing
# else; that was PREDICTED and then CHECKED against a full regeneration rather
# than assumed. #118 and #134 both wanted this membership pinned and never had it.
run step2_plink_group_multiweight.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_multiweight.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE \
  --is_output_markerList_in_groupTest=TRUE

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
  --condition=rs30,rs79 --is_fastTest=TRUE \
  --is_output_markerList_in_groupTest=TRUE

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
  --is_fastTest=TRUE \
  --is_output_markerList_in_groupTest=TRUE

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
  --condition=rs30,rs79 --is_fastTest=TRUE \
  --is_output_markerList_in_groupTest=TRUE

# P34 (added with the LEDGER #17 fix): K=1 CUSTOM-weighted group file WITH
# --condition -- the one combination that read an uninitialised `w0`. Nothing
# covered it: P3 conditions on an UNWEIGHTED file (no custom set at all), P10
# conditions on the K=2 file (which takes the raw-export path and reads no w0),
# and the C11 reconciliation test deliberately used the BETA route for exactly
# this reason. That is why the UB survived every gate, and it is why this golden
# exists rather than a test alone.
#
# Its weight line is byte-identical to the CADD line of
# group_new_snpid_multiWeights.txt, so P10's CADD rows are the independent
# comparison: the two paths now agree on every conditional column (measured
# 0.000e+00 relative on all 12 cells; before the fix, up to 2.1e-01 on
# BETA_Burden_cond, with the adjustment vanishing entirely on the 6 all-ultra-
# rare cells). No PGEN mirror: the defect is in the region path's weight
# bookkeeping, which no reader can reach differently.
run step2_plink_group_cond_customweight.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_cond_customweight.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_withWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --condition=rs30,rs79 --is_fastTest=TRUE

# P25/P26 (added with the LEDGER #41 fix): BURDEN-only (--r.corr=1) group test
# WITH --condition. This is the combination through which the halved conditional
# p-value reached a GENE-level statistic: mainRegionInCPP parses the printed
# string back with std::stod (src/Main.cpp:2019) and feeds it to CCT_cpp, so the
# region Cauchy row inherited the halving. No golden covered r.corr=1 together
# with --condition before the fix -- which is exactly why #41 survived P3/P10.
# Mirrors P5's is_fastTest=FALSE and P3/P10's rs30,rs79 conditioning set.
run step2_plink_group_cond_burden.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_cond_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE --condition=rs30,rs79

run step2_pgen_group_cond_burden.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_cond_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE --condition=rs30,rs79

# P27-P33 (S0 of the step-2 threading/speedup work): the MARKER-path coverage
# gap. Before these, every single-variant golden (P1/P4/P7/P11/P12) ran the same
# one branch of the getMarkerPval cascade -- scoreTestFast_noadjCov -- so a
# change to any other branch left every golden byte-identical while being
# completely untested. Verified with an instrumented build (counters on each
# branch, run at 61ac4be9 and discarded); the per-run counts are in the commit
# message and are quoted per golden below. All seven share P1's fixture and
# filters, so a diff against P1 isolates exactly one thing.
# No PGEN mirrors: these vary the TEST kernel, not the reader, and the reader is
# already mirrored P13-P26 over the same fixture.

# P27: single-variant with a sparse GRM -- the marker-path scoreTest + PCG.
# 70/70 markers take the sparse branch (SAIGE_test.cpp:522), 70 PCG solves /
# 94 iterations, and 4 SPA calls go to SPA rather than SPA_fast (the
# p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur split at :692).
run step2_plink_marker_sparse.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_sparse.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $VR_MODEL $SG \
  --is_output_moreDetails=TRUE \
  --is_fastTest=FALSE

# P28: single-variant fastTest -- the RE-TEST block at Main.cpp:436-475, which
# no golden reached (the group goldens P2/P3/P8/P10 use fastTest, but the
# re-test is marker-loop code). 73 getMarkerPval calls for 70 markers: the 3
# markers with p < 0.05 and MAC > MACCutoffforER are re-tested, all 3 through
# scoreTestFast. The sparse arm of the re-test's own branch is live here (the
# sparse GRM is supplied) but not selected -- all 3 have MAC above the last
# category bound; P29 selects it.
run step2_plink_marker_fastTest.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_fastTest.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $VR_MODEL $SG \
  --is_output_moreDetails=TRUE \
  --is_fastTest=TRUE

# P29: fastTest re-test, SPARSE arm. --max_MAC_for_ER=0 turns the ER path off,
# which is what lets the two ultra-rare markers (rs44 MAC=1, rs56 MAC=2) reach
# the re-test at all -- it is gated on MAC > MACCutoffforER, so at the default
# cutoff of 4 they are answered by ER and never re-tested. 75 calls (5
# re-tests), of which 2 take the sparse branch and 3 scoreTestFast; ER runs 0
# times, which is itself the only golden coverage of "ER disabled".
run step2_plink_marker_fastTest_sparse.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_fastTest_sparse.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $VR_MODEL $SG \
  --is_output_moreDetails=TRUE \
  --is_fastTest=TRUE --max_MAC_for_ER=0

# P30: --is_noadjCov=FALSE, no sparse GRM -- the third arm of the branch at
# SAIGE_test.cpp:514-528, scoreTestFast, on all 70 markers. P1's model, so the
# diff against P1 is exactly "regress the covariates out of the genotype".
run step2_plink_marker_noadjCov.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_noadjCov.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE \
  --is_noadjCov=FALSE

# P31: the Firth arm -- the last step of the cascade, live in the test suite
# (test-step2-firth.R) since Phase 3 but pinned by no golden. 5 of the 70
# markers qualify at pCutoffforFirth=0.05: 3 through the SPA path
# (SAIGE_test.cpp:777/784) and 2 through the ER path (:835), a site the SPA-path
# reading of the code misses. Do NOT raise the cutoff further -- above it the
# penalised fit stops converging and drifts to maxstep*maxit (LEDGER #26).
run step2_plink_marker_firth.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_firth.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE \
  --is_Firth_beta=TRUE --pCutoffforFirth=0.05

# P32: the quantitative marker path, UNCONDITIONAL. P4 is the only other
# quantitative marker golden and it passes --condition, so the plain
# quantitative single-variant run -- the most ordinary use of step 2 there is --
# was gated by nothing. No ER and no SPA on this path (both are binary-only), so
# it is also the only marker golden where the cascade stops at the score test.
run step2_plink_marker_quantitative.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_quantitative.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $QUANT_MODEL \
  --is_output_moreDetails=TRUE

# P33: the ER RESAMPLING branch. ER itself runs in P1 (2 markers, both
# ultra-rare), but with k <= 2 carriers it always takes SKAT_Exact_Recurse; the
# resampling arm (Binary_ComputeExact.cpp:367, and the SL_setseed(1) that the
# threading verdict doc calls the single highest-risk fact in step 2) needs
# n_total = 2^k > 2e6, i.e. k >= 21 carriers, which no marker can reach while
# ER is capped at MAC <= 4. Raising the cap to 200 pulls in exactly one more
# marker, rs38 (MAC 102, 98 carriers): ER=3, and its p-value moves 1.004E-02
# (SPA) -> 1.367E-02 (ER), same order of magnitude, which is the sanity check
# that the resampling arm is not producing nonsense.
# The verdict doc proposed --max_MAC_for_ER=25 for this; MEASURED, 25 is
# byte-identical to P1 (so are 20, 50 and 100) because the extra markers it
# admits all have |Z| below the SPA cutoff and never take ER. A golden at 25
# would have gated nothing.
run step2_plink_marker_er.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_er.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
  $BIN_MODEL \
  --is_output_moreDetails=TRUE \
  --max_MAC_for_ER=200

# B1/B2 (bgen runs) and V1/V2 (vcf runs) retired to harness/goldens_retired/
# when the BGEN (Strip 6) and VCF/SAV (Strip 5) readers were stripped. Their
# PLINK counterparts P1-P4 continue; before retirement, the bgen group
# fastTest run had byte-stable goldens and the vcf conditional output was
# verified identical to P4's.

# P35 (added when --pval_cutoff_for_fastTest reached the CLI): P2's command with
# the escalation cutoff forced to 1, so BOTH bundled regions take the fastTest
# ESCALATION pass. Until the option existed this path could not be reached by a
# golden at all, which is exactly why LEDGER #50 survived -- at the default 0.05
# neither region escalates, so P2 covers only the copy_singleInGroup route to
# the single-variant-in-group rows and P6 only the unconditional write.
#
# READ THIS BEFORE ASSUMING IT IS REDUNDANT. Its bytes are identical to P2's, in
# both the region file and the .singleAssoc.txt -- that is LEDGER #51: the region
# path never consults m_isFastTest, so the escalated pass recomputes exactly what
# the fast pass already had. The golden is NOT the numbers; it is that the OTHER
# WRITER produces them. Under #50 this run lost all 19 single-variant rows at
# exit 0 while P2 stayed green. So do not "simplify" it away because it
# duplicates P2, and do not assert it against P2's file -- it has its own frozen
# copy, and the equality is a note, not a claim.
run step2_plink_group_fastTest_escalated.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_fastTest_escalated.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE --pval_cutoff_for_fastTest=1

# P36/P37 (Phase 2 of the upstream-tracker mine): the QUANTITATIVE region path,
# which had ZERO golden coverage. Every group golden above is binary, so the
# quantitative arms of region_postprocess and of the C++ region loop -- the
# `traitType == "quantitative"` branches at src/Main.cpp:2486 (N_Vec instead of
# the case/control accounting) and R/SAIGE_SPATest_Region.R's quantitative
# conditional arm -- were gated by nothing at all.
#
# P36 is the one LEDGER #72 needs: --r.corr=1 puts the run through the BURDEN
# Cauchy aggregation at src/Main.cpp:2789-2800, whose `std::stod` is BARE where
# the single-marker twin at :561 is wrapped. It cannot make that parse throw on
# this fixture (no bundled region produces a subnormal p; the two Cauchy rows are
# 0.893901 and 0.472394) -- what it does is put a QUANTITATIVE trait on that
# line, which matters because the quantitative path has no SPA rescue arm, so an
# underflowed score p reaches the writer unfiltered. Six goldens already execute
# the line; all six are binary.
#
# NOTE the missing $SG: quant_sparseGRM is a --useSparseGRMtoFitNULL fit, so its
# variance-ratio file has only null/null_noXadj rows and step 2 must NOT be given
# --sparseGRMFile for it (see the step1->step2 contract in CLAUDE.md). That is
# also why these two use $QUANT_MODEL rather than $GRP_MODEL.
run step2_plink_group_burden_quant.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_burden_quant.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_MODEL \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE

# P37: the same model through the DEFAULT SKAT-O region test, so the
# quantitative region path is covered on its main arm as well as on BURDEN.
# is_fastTest=FALSE mirrors P6.
run step2_plink_group_nofast_quant.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_nofast_quant.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_MODEL \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=FALSE

# P38-P41 (2026-09-08): the QUANTITATIVE sparse-GRM step-2 path, which had NO
# golden of any kind. Measured before they were added: of 38 step-2 goldens, 23
# of the 33 binary ones passed --sparseGRMFile and 0 of the 5 quantitative ones
# did, because no quantitative step-1 golden had `sparse` variance-ratio rows.
# quant_sparseGRM_vr above supplies them.
#
# This is not an abstract gap. LEDGER #106 -- step 2's sparse Sigma carrying one
# triangle, and its PCG stopping on an absolute residual -- says in its own WHAT
# cell that the error is "small here only because 1/W ~ 11 dominates
# tau2 x kinship ~ 0.16; on a QUANTITATIVE trait (diagonal tau1 ~ 1) or a
# strongly related cohort the missing triangle is first-order". The regime where
# that bug was largest is the one no golden looked at: the #106 fix moved 43
# binary goldens and, necessarily, not one quantitative byte.

# P38: group test, quantitative, SKAT-O with a sparse GRM. The quantitative twin
# of P6 and the main one of this group -- it is the run that puts a quantitative
# model through setSparseSigma_new and getPCG1ofSigmaAndGtilde.
run step2_plink_group_nofast_quant_sparse.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_nofast_quant_sparse.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_VR_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=FALSE

# P39: the same model on the BURDEN-only arm (r.corr=1), the quantitative twin
# of P5. Reaches the collapsed-vector burden variance with a sparse GRM, which
# is the construction LEDGER #104/#105 turned on.
run step2_plink_group_burden_quant_sparse.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_burden_quant_sparse.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_VR_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --r.corr=1 --is_fastTest=FALSE

# P40: single-variant with a sparse GRM, the quantitative twin of P27 -- the
# marker-path scoreTest + PCG arm, which is where LEDGER #44's third shape and
# #106's var column both live. Binary had this arm covered; quantitative did not.
run step2_plink_marker_sparse_quant.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_marker_sparse_quant.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_VR_MODEL $SG \
  --is_fastTest=FALSE

# P41: PGEN mirror of P38 -- the only quantitative PGEN region golden, and the
# only place the pgen reader meets a quantitative sparse-GRM region test.
run step2_pgen_group_nofast_quant_sparse.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_nofast_quant_sparse.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_VR_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=FALSE

# P42-P48 (2026-09-08): MISSING GENOTYPE CALLS, which no fixture had.
# Every bundled .bed carries EXACTLY ZERO missing calls -- 7 files, ~150 million
# calls, verified by decoding them -- so `--maxMissing`, `--impute_method`,
# `--dosage_zerod_cutoff`, `--dosage_zerod_MAC_cutoff` and `--is_imputed_data`
# had no golden AND no test. The fixture is generated, not hand-made:
#   python3 harness/tools/make_missing_fixture.py \
#       extdata/input/genotype_100markers extdata/input/genotype_100markers_missing
#   plink2 --bfile <that> --make-pgen --out <that>
# and it is byte-reproducible from that script alone (two runs, cmp).
#
# The ladder is on the marker index, sized against the MODEL's 1,000 samples
# (fam rows 1..1000) because `missingRate = numMissing / N` is computed over the
# analysed SUBSET, not the .fam's 10,000 rows:
#   j%10 :   0     1      2      3     4     5     6      7      8     9
#   rate :   0.0   0.001  0.005  0.01  0.05  0.10  0.15   0.151  0.30  0.60
# `j%10 == 6` sits EXACTLY on the --maxMissing default and must PASS, because
# the gate is a strict `>`. Every marker also carries 500 missing calls OUTSIDE
# the analysed subset which must never reach the output.
#
# All seven runs are P1's command with one thing changed, so the comparisons are
# one-variable. Measured differences, checked BEFORE these were accepted as
# goldens -- a golden that cannot differ from its twin gates nothing:
#   P42 vs P44 : exactly 6 markers (rs23 rs26 rs64 rs66 rs94 rs97), the ones
#                with AF >= 0.25 where best_guess imputes 1 and minor imputes 0
#   P42 vs P45 : 48 markers against 70; the 22 above the cutoff come back
#   P43 vs P46 : exactly the 20 passing markers with MAC <= 10, where
#                --dosage_zerod_cutoff=0.2 zeroes what `mean` just imputed
#   P42 vs P47 : BYTE-IDENTICAL by design (reader parity under missingness)
#   P42 vs P48 : column 8 SWAPPED, MissingRate -> imputationInfo (always 1); the
#                other 22 columns identical. P48 is SMALLER than P42, which is how
#                the swap was noticed at all -- an added column cannot shrink a file
MISS_G="--bedFile=./input/genotype_100markers_missing.bed --bimFile=./input/genotype_100markers_missing.bim --famFile=./input/genotype_100markers_missing.fam --AlleleOrder=alt-first"
MISS_PGEN_G="--pgenFile=./input/genotype_100markers_missing.pgen --pvarFile=./input/genotype_100markers_missing.pvar --psamFile=./input/genotype_100markers_missing.psam --AlleleOrder=alt-first"
MISS_COMMON="--chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE --is_output_moreDetails=TRUE"

# P42: the base -- CLI defaults, so impute_method=best_guess and maxMissing=0.15.
run step2_plink_marker_missing.log R2 step2_SPAtests.R \
  $MISS_G --SAIGEOutputFile="$GOLD/step2/plink_marker_missing.txt" \
  $MISS_COMMON $BIN_MODEL

# P43: mean imputation. imputeG = 2*AF rather than round(2*AF).
run step2_plink_marker_missing_mean.log R2 step2_SPAtests.R \
  $MISS_G --SAIGEOutputFile="$GOLD/step2/plink_marker_missing_mean.txt" \
  $MISS_COMMON $BIN_MODEL --impute_method=mean

# P44: minor imputation, imputeG = 0. Differs from P42 only where round(2*AF)
# is 1, i.e. on the six markers with AF >= 0.25 that pass the missing-rate gate.
run step2_plink_marker_missing_minor.log R2 step2_SPAtests.R \
  $MISS_G --SAIGEOutputFile="$GOLD/step2/plink_marker_missing_minor.txt" \
  $MISS_COMMON $BIN_MODEL --impute_method=minor

# P45: the missing-rate gate switched off, so the 0.151/0.30/0.60 markers come
# back. This is the only golden where a marker above the default cutoff is
# tested at all.
run step2_plink_marker_missing_maxmiss1.log R2 step2_SPAtests.R \
  $MISS_G --SAIGEOutputFile="$GOLD/step2/plink_marker_missing_maxmiss1.txt" \
  $MISS_COMMON $BIN_MODEL --maxMissing=1

# P46: mean imputation with the dosage-zeroing switched OFF. Against P43 this
# isolates `t_GVec.clean(t_dosage_zerod_cutoff)` (src/UTIL.cpp:108-112), which
# on hard calls silently undoes mean imputation for every marker with
# MAC <= --dosage_zerod_MAC_cutoff: 2*AF is ~0.003 there, below the 0.2 default.
run step2_plink_marker_missing_meannozerod.log R2 step2_SPAtests.R \
  $MISS_G --SAIGEOutputFile="$GOLD/step2/plink_marker_missing_meannozerod.txt" \
  $MISS_COMMON $BIN_MODEL --impute_method=mean --dosage_zerod_cutoff=0

# P47: PGEN mirror of P42. Byte-identical to it, which is the point: no earlier
# PGEN mirror had a single missing call to disagree about.
run step2_pgen_marker_missing.log R2 step2_SPAtests.R \
  $MISS_PGEN_G --SAIGEOutputFile="$GOLD/step2/pgen_marker_missing.txt" \
  $MISS_COMMON $BIN_MODEL

# P48 RETIRED 2026-09-09 (LEDGER #108). It ran --is_imputed_data=TRUE, which
# does NOT add a column -- it EXCHANGES `MissingRate` for `imputationInfo`,
# which on this branch can only ever read 1, because both readers assign
# `t_imputeInfo = 1` unconditionally (src/PLINK.cpp:312, src/PGEN.cpp:585). So
# the flag cost a real column and returned a constant, and it is now refused.
# There is no output left to compare; the refusal is pinned by
# tests/testthat/test-step2-missing.R. The option and the writer branch both
# stay, so restoring a BGEN/VCF reader restores this golden with it.

# P49-P51 (2026-09-08): DOSAGES, which only the pgen reader can carry and which
# no golden touched. `--dosage_zerod_cutoff`/`--dosage_zerod_MAC_cutoff` zero
# every dosage below the cutoff when the marker's MAC is at or below the MAC
# cutoff (src/UTIL.cpp:108-112). On hard calls the only sub-threshold values are
# the ones `--impute_method=mean` writes (LEDGER #109, goldens P43/P46); on
# dosages the option reaches the DATA, and the pgen dosage decode path
# (`PgrGet1D` + `Dosage16ToDoubles`, src/PGEN.cpp:485-489) runs for the first
# time in a golden.
#
# extdata/input/dosage_20markers_100samples.* already carried real dosages but
# has 100 samples and is wired only to test-step1-pgen.R; the bundled models
# carry 1,000, so step 2 needed a dosage pgen over THAT sample set. Two steps,
# both deterministic:
#   python3 harness/tools/make_dosage_fixture.py \
#       extdata/input/genotype_100markers <scratch>/dosfix
#   plink2 --import-dosage <scratch>/dosfix.dosage noheader format=1 ref-last \
#          --fam <scratch>/dosfix.fam --map <scratch>/dosfix.map \
#          --hard-call-threshold 0.499 --make-pgen \
#          --out extdata/input/genotype_100markers_dosage
# The .dosage/.map intermediates are NOT committed (2 MB of text that the script
# reproduces byte-identically); the .pgen/.pvar/.psam are, at 228 KB total.
# --hard-call-threshold 0.499 is deliberate: every dosage then also gets a hard
# call, so the fixture introduces no missing genotypes and varies ONE thing.
#
# Every sample keeps its exact hard call except twelve hom-ref samples per
# marker, which take a ladder straddling the 0.2 default:
#   0.02 0.05 0.10 0.15 0.19 0.199 | 0.21 0.25 0.30 0.50 0.80 1.30
# Only hom-ref samples are overwritten, so MAC can only rise, and by at most
# sum = 4.069 -- which is what keeps 59 markers at or below the MAC gate and 41
# above it. All 100 markers are emitted (against 70 for the hard-call fixture)
# because the 30 monomorphic ones become polymorphic.
#
# Measured before these were accepted, and both differ-sets are IDENTICAL to the
# MAC band rather than merely overlapping it:
#   P49 vs P50 : exactly the 59 markers with MAC <= 10
#   P49 vs P51 : exactly the 17 markers with 10 < MAC <= 100
# All 100 AC_Allele2 values are fractional, which is the proof the dosages
# reached the statistics rather than being hard-called on the way in.
DOSAGE_G="--pgenFile=./input/genotype_100markers_dosage.pgen --pvarFile=./input/genotype_100markers_dosage.pvar --psamFile=./input/genotype_100markers_dosage.psam --AlleleOrder=alt-first"
DOSAGE_COMMON="--chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE --is_output_moreDetails=TRUE"

# P49: CLI defaults -- dosage_zerod_cutoff=0.2, dosage_zerod_MAC_cutoff=10.
run step2_pgen_marker_dosage.log R2 step2_SPAtests.R \
  $DOSAGE_G --SAIGEOutputFile="$GOLD/step2/pgen_marker_dosage.txt" \
  $DOSAGE_COMMON $BIN_MODEL

# P50: the zeroing switched off. Against P49 this isolates `clean()` on real
# dosages, on exactly the markers the MAC gate admits.
run step2_pgen_marker_dosage_nozerod.log R2 step2_SPAtests.R \
  $DOSAGE_G --SAIGEOutputFile="$GOLD/step2/pgen_marker_dosage_nozerod.txt" \
  $DOSAGE_COMMON $BIN_MODEL --dosage_zerod_cutoff=0

# P51: the MAC gate widened from 10 to 100, so the clean reaches 17 more
# markers. This is the only golden that moves --dosage_zerod_MAC_cutoff at all.
run step2_pgen_marker_dosage_mac100.log R2 step2_SPAtests.R \
  $DOSAGE_G --SAIGEOutputFile="$GOLD/step2/pgen_marker_dosage_mac100.txt" \
  $DOSAGE_COMMON $BIN_MODEL --dosage_zerod_MAC_cutoff=100

# P52-P53 (2026-09-08): the step-2 half of LOCO, which no golden reached.
# `ReadModel` overwrites the model's fitted values from `LOCOResult[[chrom]]`
# (--chrom is mandatory when LOCO=TRUE), so the SAME model tested on two
# chromosomes uses two different refits. genotype_100markers_2chr carries 90
# markers on chr1 and 10 on chr2 and shares the sample IDs.
#
# Measured before these were accepted, against the same commands with
# --LOCO=FALSE on the same model -- which is the contrast that proves the swap
# reaches the output rather than merely being present in the .rda:
#   chr1: all 63 emitted markers move, median 6.5% and max 52% in p.value
#   chr2: all  7 emitted markers move, median 1.3% and max  6.0%
# The ratio is the point and it was predicted before the run: LOCOResult[[1]]
# leaves out 9,988 markers and LOCOResult[[2]] leaves out 12.
LOCO_G2="--bedFile=./input/genotype_100markers_2chr.bed --bimFile=./input/genotype_100markers_2chr.bim --famFile=./input/genotype_100markers_2chr.fam --AlleleOrder=alt-first"
LOCO_MODEL="--GMMATmodelFile=$GOLD/step1/binary_loco_2chr.rda --varianceRatioFile=$GOLD/step1/binary_loco_2chr.varianceRatio.txt"
LOCO_COMMON="--minMAF=0 --minMAC=1 --LOCO=TRUE --is_output_moreDetails=TRUE"

# P52: chromosome 1 left out of the GRM that the model was refitted on.
run step2_plink_marker_loco_chr1.log R2 step2_SPAtests.R \
  $LOCO_G2 --SAIGEOutputFile="$GOLD/step2/plink_marker_loco_chr1.txt" \
  --chrom=1 $LOCO_COMMON $LOCO_MODEL

# P53: chromosome 2, the same model and command with only --chrom changed, so
# the pair isolates LOCOResult[[chrom]] and nothing else.
run step2_plink_marker_loco_chr2.log R2 step2_SPAtests.R \
  $LOCO_G2 --SAIGEOutputFile="$GOLD/step2/plink_marker_loco_chr2.txt" \
  --chrom=2 $LOCO_COMMON $LOCO_MODEL

# P54-P55 (2026-09-08): the REGION path under LOCO, which no golden reached --
# every group golden runs --LOCO=FALSE, so `ReadModel`'s per-chromosome swap had
# been gated only on the single-variant path (P52/P53).
#
# The group file is `group_new_snpid.txt`, the SAME one the other region goldens
# use, against the two-chromosome genotypes: its GENE1 is entirely on chr1 while
# GENE2 has 40 markers on chr1 and 10 on chr2, so GENE2 SPANS THE BOUNDARY and
# is tested against a different refit in each run. (The repo's
# `group_new_chrposa1a2_2chr.txt` could not be used here: PLINK group files key on
# bim column 2, and that file's chr:pos:A1:A2 ids match nothing in this .bim --
# "100 markers in 'RegionFile' are not in 'GenoFile'", a zero-byte output at
# exit 0. It was built for the retired bgen/vcf fixtures.)
#
# Measured before these were accepted, against the same commands with
# --LOCO=FALSE on the same model:
#   chrom 1 (14 cells): every cell moves; max 39% on the three p-value columns
#                       and 156% on BETA_Burden
#   chrom 2 ( 4 cells): every cell moves; max 0.65% on Pvalue, 1.4% BETA_Burden
# Same asymmetry, and same cause, as the marker path: LOCOResult[[1]] leaves out
# 9,988 GRM markers and LOCOResult[[2]] leaves out 12.
LOCO_R_MODEL="--GMMATmodelFile=$GOLD/step1/binary_loco_vr_2chr.rda --varianceRatioFile=$GOLD/step1/binary_loco_vr_2chr.varianceRatio.txt"
LOCO_R_COMMON="--minMAF=0 --minMAC=0.5 --LOCO=TRUE --is_fastTest=FALSE --groupFile=./input/group_new_snpid.txt --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous --maxMAF_in_groupTest=0.0001,0.001,0.01"

# P54: SKAT-O region test under LOCO, chromosome 1 (GENE1 + GENE2's chr1 half).
run step2_plink_group_loco_chr1.log R2 step2_SPAtests.R \
  $LOCO_G2 --SAIGEOutputFile="$GOLD/step2/plink_group_loco_chr1.txt" \
  --chrom=1 $LOCO_R_COMMON $LOCO_R_MODEL $SG

# P55: the same, chromosome 2 (GENE2's chr2 half), so the pair differs only in
# --chrom and therefore only in which refit ReadModel installs.
run step2_plink_group_loco_chr2.log R2 step2_SPAtests.R \
  $LOCO_G2 --SAIGEOutputFile="$GOLD/step2/plink_group_loco_chr2.txt" \
  --chrom=2 $LOCO_R_COMMON $LOCO_R_MODEL $SG

# P56-P57 (2026-09-08): the QUANTITATIVE region path's two remaining shapes.
# The coverage audit's largest surviving asymmetry: binary had 14 fastTest, 9
# multi-weight and 9 conditional region goldens and quantitative had NONE of the
# three. Both of these use $QUANT_VR_MODEL, the only quantitative model with
# `sparse` variance-ratio rows, so they exercise the same sparse-GRM region path
# LEDGER #106 was first-order on.
#
# There is deliberately NO quantitative fastTest golden. Measured first: the
# quantitative `--is_fastTest=TRUE` run is BYTE-IDENTICAL to P38, exactly as the
# binary pair plink_group_fastTest/plink_group_nofast is, because the region path
# never reads the flag (LEDGER #51). A third file with P38's bytes would gate
# nothing, so that invariant is asserted in tests/testthat/test-step2-group-quant.R
# instead -- where it also states, for the quantitative path, the fact #51
# records for the binary one.
#
# Both were checked to DIFFER from P38 before being accepted, and both bring a
# schema no quantitative golden had: P56 the K>1 multi-weight columns
# (`Weight`, `Pvalue_ACATV`, `Pvalue_ACATO`; 26 rows against P38's 14) and P57
# the five `_cond` columns.

# P56: multi-weight (K=2) SKAT-O on a quantitative model. The quantitative twin
# of plink_group_multiweight.
run step2_plink_group_multiweight_quant.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_multiweight_quant.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_VR_MODEL $SG \
  --groupFile=./input/group_new_snpid_multiWeights.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --is_fastTest=TRUE \
  --is_output_markerList_in_groupTest=TRUE

# P57: conditional region test on a quantitative model, the quantitative twin of
# plink_group_cond_binary. Conditioning markers are rsIDs (bim column 2).
run step2_plink_group_cond_quant.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_cond_quant.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_VR_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.0001,0.001,0.01 \
  --condition=rs30,rs79 --is_fastTest=TRUE

# P58-P61 (2026-09-09, FlexRV F4): the FlexRV transform sets
# (docs/design/flexrv-design.md) on extdata/input/group_new_snpid_flexrv.txt --
# group_new_snpid.txt plus a `score:REVEL` line per gene (LoF markers 1,
# missense markers the REVEL values of the multiWeights fixture).
#
# THE FIXTURE ARITHMETIC PINS THESE FLAGS. The model holds 1,000 samples, so
# max MAF 0.01 admits MAC <= 20: 21 polymorphic markers in GENE1 (MAC 1..20)
# and 18 in GENE2 (MAC 1..13). At the default collapse cutoff the real markers
# are rs4/rs17/rs31/rs33 and rs60/rs70/rs91; MAC25..MAC500 are then duplicates
# of unw, and MAC1/3/5/10 select only collapsed markers. 0.001 (MAC <= 2) would
# leave no real marker at all, so a golden on the paper's boundary would gate
# nothing (the "a golden that reproduces P1 gates nothing" trap).
#
# P58: the default 192 sets, binary BURDEN, one cell (lof:missense at 0.01),
# ultra-rare collapse on -- the recommended binary run. The C++ BURDEN writer
# emits 192 rows + 1 Cauchy row per gene.
run step2_plink_group_flexrv_burden.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_flexrv_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_flexrv.txt \
  --flexRV_score=REVEL --flexRV_maxMAF=0.01 \
  --annotation_in_groupTest=lof:missense \
  --r.corr=1 --is_fastTest=FALSE

# P59: the paper's own configuration -- quantitative trait, NO ultra-rare
# collapse (--MACCutoff_to_CollapseUltraRare=0, the unbinned additive burden),
# so the whole MAC ladder is live on real markers. No $SG: the quantitative
# model's variance-ratio file has only `null` rows.
run step2_plink_group_flexrv_burden_quant.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_flexrv_burden_quant.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $QUANT_MODEL \
  --groupFile=./input/group_new_snpid_flexrv.txt \
  --flexRV_score=REVEL --flexRV_maxMAF=0.01 \
  --annotation_in_groupTest=lof:missense \
  --MACCutoff_to_CollapseUltraRare=0 \
  --r.corr=1 --is_fastTest=FALSE

# P60: SKAT-O on a REDUCED spec (unw;sbeta:2 x unw;mafbeta:2 = 4 sets), so the
# R writer, the transformWeightMat export and the per-set ACAT-V/ACAT-O columns
# are covered without 192 Davies integrations per cell. Max MAF 0.02 rather
# than 0.01, and mafbeta rather than a MAC threshold, so that NO marker in the
# cell gets weight exactly 0: a zero weight inside a binary SKAT-O cell makes
# the kernel return NA (LEDGER #117, pre-existing, reached by custom weights
# too), and a golden that encodes that failure would gate the bug, not the
# feature. At 0.02 the largest in-cell MAF is 0.017 (GENE1 rs37, MAC 34), so
# every mafbeta@0.02 weight is positive.
run step2_plink_group_flexrv_skato.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_flexrv_skato.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_flexrv.txt \
  --flexRV_score=REVEL --flexRV_maxMAF=0.02 \
  --flexRV_scoreTransforms="unw;sbeta:2" --flexRV_MAFTransforms="unw;mafbeta:2" \
  --annotation_in_groupTest=lof:missense \
  --is_fastTest=FALSE \
  --is_output_markerList_in_groupTest=TRUE

# P61: P58 through the pgen reader (the in-sample MAC/MAF the MAC ladder and
# the mafbeta transforms read must not depend on the reader).
run step2_pgen_group_flexrv_burden.log R2 step2_SPAtests.R \
  $PGEN_G \
  --SAIGEOutputFile="$GOLD/step2/pgen_group_flexrv_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid_flexrv.txt \
  --flexRV_score=REVEL --flexRV_maxMAF=0.01 \
  --annotation_in_groupTest=lof:missense \
  --r.corr=1 --is_fastTest=FALSE

# P62-P64 (2026-09-12): the two defects of the ultra-rare pseudo-marker pass,
# LEDGER #129 and #130, both now FIXED. These three were added at unpatched
# HEAD so the fixes would MOVE something, and all three were byte-reproducible
# there (two runs into separate trees, cmp) before being accepted.
#
# TWO CLAIMS THIS BLOCK ORIGINALLY MADE WERE WRONG, and the corrections matter
# more than the claims did:
#  - "0 of the 269 UR rows across every group golden that writes a
#    .singleAssoc.txt has AF > 0.5". The count is 359, not 269 (measured over
#    the pre-P62 goldens). None of them has AF > 0.5, so the coverage claim
#    holds; the denominator was simply wrong, and it had been copied from the
#    ledger row.
#  - "no group golden has a custom weight small enough to be mistaken for an
#    imputed dosage". FALSE, and falsified by the #129 fix itself, which moved
#    17 pre-existing golden files: group_new_snpid_multiWeights.txt carries a
#    REVEL set with weights of 0.01-0.07, every one of them at or below the
#    --dosage_zerod_cutoff of 0.2. P64 is still worth having -- it is the
#    minimal two-marker reproduction, where the whole cell vanishes rather than
#    shifting -- but it is not the only golden that reaches #129.
#
# P62/P63 (LEDGER #130): the flip drill, both region paths. maxMAF 0.5 admits
# every marker of both genes and a collapse cutoff above 2N makes every one of
# them ultra-rare, so each (region x annotation) cell is ONE collapsed
# pseudo-marker and five of the six carry alt frequency > 0.5, which is the
# regime #130 used to fold. A cell of one marker HAS that marker's burden, so
# the region row's Pvalue_Burden and the run's own .singleAssoc.txt p are the
# same quantity computed twice -- and since the fix they AGREE, on all six
# cells and on both arms. `GENE1 lof` (AF 0.374, never folded) is the control:
# it agreed before the fix too and is byte-identical across it.
# --is_single_in_groupTest=TRUE so both halves are in the goldens.
# These two files also record MAC != MAC_case + MAC_control (444 against
# 151 + 1405): the group MAC is folded and its case/control split is not. That
# is LEDGER #133, pre-existing and unrelated to #130 -- these are simply the
# first goldens to exhibit it.
run step2_plink_group_urflip_burden.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_urflip_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.5 --MACCutoff_to_CollapseUltraRare=2000 \
  --is_single_in_groupTest=TRUE --r.corr=1 --is_fastTest=FALSE

# P63: the same drill on the SKAT-O path, which reaches the pseudo-marker
# through the per-marker exports instead of the BURDEN emitter -- and, on a
# BINARY trait, reads the burden column too (`g.sum` in
# get_newPhi_scaleFactor_traitType, R/SAIGE_SPATest_Region.R). So this file is
# the gate on the OTHER route into the same column.
run step2_plink_group_urflip_skato.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_urflip_skato.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.5 --MACCutoff_to_CollapseUltraRare=2000 \
  --is_single_in_groupTest=TRUE --is_fastTest=FALSE

# P64 (LEDGER #129, FIXED): input/group_ur_smallWeights.txt is GENET over rs28
# and rs32 (two singletons, cell MAC 2) with ONE custom weight set at 0.01.
# Both markers are ultra-rare, so the cell is one weighted pseudo-marker whose
# every entry is 0.01. That used to be below --dosage_zerod_cutoff (0.2) at a
# weighted MAC below --dosage_zerod_MAC_cutoff (10), so the filter took it for
# imputation noise, zeroed it, and the cell emitted NO row -- a bare Cauchy NA
# row plus a #127 sidecar reading 0 / 0 / `none`. It is now tested, at
# Pvalue_Burden 2.522646E-01, which is the value the SAME data gives at
# --MACCutoff_to_CollapseUltraRare=0 and at weights 0.5: a burden p cannot
# depend on the weight SCALE, and this golden is where that is pinned. BETA
# 232.633 is the weights-0.01 coefficient, 50x the weights-0.5 one, which is
# the other half of the same statement. test-step2-region-ur-zerod.R holds the
# three configurations against each other.
run step2_plink_group_ur_smallweights.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_ur_smallweights.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_ur_smallWeights.txt \
  --annotation_in_groupTest=lof --maxMAF_in_groupTest=0.01 \
  --MACCutoff_to_CollapseUltraRare=10 --r.corr=1 --is_fastTest=FALSE

# P65 (2026-09-12): the MIXED cell -- the only bundled configuration in which a
# cell holds BOTH surviving rare markers AND a collapsed pseudo-marker above alt
# frequency 0.5. It exists because that is the regime where the coding of the
# pseudo-marker changes the STATISTIC rather than only its sign, and nothing
# gated it: P62/P63 collapse everything (`Number_rare` 0, so the pseudo-marker
# IS the cell) and every other group golden sits below 0.5.
#
# --MACCutoff_to_CollapseUltraRare=700 is the whole fixture arithmetic, found by
# sweeping the cutoff between the default 10 and P62's 2000: at 50/100/200 no
# pseudo-marker reaches 0.5, at 1000 nothing survives as rare, at 400 two cells
# are mixed and at 700 FIVE are (measured -- two pairs of them duplicates, since
# no marker of the fixture is annotated synonymous, so `missense;lof;synonymous`
# selects what `missense;lof` does). Those five carry 28.3 %, 28.3 %, 28.7 %,
# 28.7 % and 52.5 % of their cell's MAC in the pseudo-marker, which is what
# makes the file discriminating: folding it (LEDGER #130's rejected option)
# would enter `w*(2 - g) = const - w*g`, so a quarter to a half of each cell's
# allele content would oppose the rest and the cell's burden would become a
# signed contrast between its ultra-rare and rare halves rather than their sum.
# `GENE1 lof` stays at AF 0.374 with `Number_rare` 0 and is the control row --
# it is byte-identical to its P62 row.
run step2_plink_group_urmixed_burden.log R2 step2_SPAtests.R \
  $PLINK_G \
  --SAIGEOutputFile="$GOLD/step2/plink_group_urmixed_burden.txt" \
  --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
  $GRP_MODEL $SG \
  --groupFile=./input/group_new_snpid.txt \
  --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
  --maxMAF_in_groupTest=0.5 --MACCutoff_to_CollapseUltraRare=700 \
  --is_single_in_groupTest=TRUE --r.corr=1 --is_fastTest=FALSE

############ PGEN-vs-PLINK byte identity, checked in the candidate dir #########
# The step-0 and step-1 PGEN runs above are NOT committed as separate goldens:
# each is byte-identical to its PLINK twin, so committing them would double the
# size of harness/goldens/step1 for no extra signal. The identity is asserted
# HERE instead, between the two runs of THIS candidate directory, which is both
# stronger than compare_goldens.sh's step-1 tolerance and free in git.
#
# (compare_goldens.sh globs $GOLD/step1/*.rda over the COMMITTED goldens, so the
# extra *_pgen files here are simply not looked at by it. This block is what
# stops that being a silent gap. tests/testthat/test-step1-pgen.R makes the same
# assertion from the test suite, against the committed PLINK goldens.)
echo "### byte identity: PGEN reader vs PLINK reader (step0/step1)"
PGEN_IDENTITY_FAILS=0
check_identical() { # check_identical <plink-file> <pgen-file>
  if cmp -s "$1" "$2"; then
    echo "  ok   $(basename "$2")"
  else
    echo "  FAIL $(basename "$2") differs from $(basename "$1")"
    PGEN_IDENTITY_FAILS=$((PGEN_IDENTITY_FAILS+1))
  fi
}
S0SUF="_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
check_identical "$GOLD/step0/test_sparseGRM$S0SUF" "$GOLD/step0/test_sparseGRM_pgen$S0SUF"
check_identical "$GOLD/step0/test_sparseGRM$S0SUF.sampleIDs.txt" \
                "$GOLD/step0/test_sparseGRM_pgen$S0SUF.sampleIDs.txt"
for n in binary_sparseGRM_a binary_sparseGRM_b quant_sparseGRM binary_denseGRM binary_sparseGRM_vr quant_sparseGRM_vr binary_loco_2chr binary_loco_vr_2chr; do
  check_identical "$GOLD/step1/$n.rda"                "$GOLD/step1/${n}_pgen.rda"
  check_identical "$GOLD/step1/$n.varianceRatio.txt"  "$GOLD/step1/${n}_pgen.varianceRatio.txt"
done
if [ "$PGEN_IDENTITY_FAILS" -gt 0 ]; then
  echo "=== PGEN/PLINK byte identity FAILED ($PGEN_IDENTITY_FAILS file(s))"
  FAILED_RUNS+=("pgen_plink_byte_identity")
fi

if [ "${#FAILED_RUNS[@]}" -gt 0 ]; then
  echo "=== COMPLETE WITH FAILURES ($GOLD):"
  printf '  FAILED: %s\n' "${FAILED_RUNS[@]}"
  exit 1
fi
echo "=== all runs complete: $GOLD"
