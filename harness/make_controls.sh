#!/usr/bin/env bash
# Negative control: permuted phenotype through the full step1 -> step2 pipeline.
# The pipeline must NOT fabricate signal: step-2 p-values should be ~Uniform(0,1).
#   bash harness/make_controls.sh <golden-dir>
# Outputs land in <golden-dir>/control/. Also byte-comparable across versions
# (the permutation seed is fixed), so the control doubles as a regression fixture.

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
GOLD="$(mkdir -p "$1" && cd "$1" && pwd)"
cd "$REPO/extdata"
source "$REPO/harness/env.sh"
mkdir -p "$GOLD/control" "$GOLD/logs"

run() { local log="$GOLD/logs/$1"; shift; echo "=== $log"
  "$@" >"$log" 2>&1 || { echo "FAILED: see $log"; tail -20 "$log"; exit 1; }; }
R2() { pixi run --manifest-path "$REPO/pixi.toml" Rscript "$@"; }

# 1. permuted phenotype (fixed seed; regenerated identically every time)
R2 -e '
  set.seed(42)
  d <- read.table("./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
                  header = TRUE, stringsAsFactors = FALSE)
  d$y_binary <- sample(d$y_binary)
  write.table(d, "'"$GOLD"'/control/pheno_permuted.txt",
              quote = FALSE, row.names = FALSE, sep = " ")
'

# 2. step 1 on the permuted phenotype (same flags as golden step1 run "a")
run control_step1.log R2 step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --phenoFile="$GOLD/control/pheno_permuted.txt" \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/control/null_binary"

# 2b. the SAME step-1 fit through the PGEN reader. Two things this buys that
# the goldens do not: the model is fitted on a PERMUTED phenotype, so it is a
# configuration no step-1 golden covers, and the comparison is byte-exact
# rather than toleranced. If the pgen front-end ever stopped re-encoding
# variants identically, this is a second place it would show.
run control_step1_pgen.log R2 step1_fitNULLGLMM.R \
  --plink2File=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
  --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
  --useSparseGRMtoFitNULL=TRUE \
  --phenoFile="$GOLD/control/pheno_permuted.txt" \
  --phenoCol=y_binary \
  --covarColList=x1,x2,a9,a10 \
  --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE \
  --outputPrefix="$GOLD/control/null_binary_pgen"

echo "=== byte identity: step-1 PGEN vs PLINK on the permuted phenotype"
for ext in .rda .varianceRatio.txt; do
  if cmp -s "$GOLD/control/null_binary$ext" "$GOLD/control/null_binary_pgen$ext"; then
    echo "  ok   null_binary_pgen$ext"
  else
    echo "  FAIL null_binary_pgen$ext differs from null_binary$ext"
    exit 1
  fi
done

# 3. step 2, PLINK single-variant, on the permuted-phenotype model
run control_step2.log R2 step2_SPAtests.R \
  --bedFile=./input/genotype_100markers.bed \
  --bimFile=./input/genotype_100markers.bim \
  --famFile=./input/genotype_100markers.fam \
  --AlleleOrder=alt-first \
  --SAIGEOutputFile="$GOLD/control/null_marker.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 \
  --GMMATmodelFile="$GOLD/control/null_binary.rda" \
  --varianceRatioFile="$GOLD/control/null_binary.varianceRatio.txt" \
  --LOCO=FALSE

# 4. uniformity check
R2 "$REPO/harness/check_null.R" "$GOLD/control/null_marker.txt"

# 5. step 2, PGEN single-variant, same permuted-phenotype model -- confirms
# the new reader doesn't fabricate signal either (Phase 7).
run control_step2_pgen.log R2 step2_SPAtests.R \
  --pgenFile=./input/genotype_100markers.pgen \
  --pvarFile=./input/genotype_100markers.pvar \
  --psamFile=./input/genotype_100markers.psam \
  --SAIGEOutputFile="$GOLD/control/null_marker_pgen.txt" \
  --chrom=1 --minMAF=0 --minMAC=1 \
  --GMMATmodelFile="$GOLD/control/null_binary.rda" \
  --varianceRatioFile="$GOLD/control/null_binary.varianceRatio.txt" \
  --LOCO=FALSE

# 6. uniformity check, PGEN path
R2 "$REPO/harness/check_null.R" "$GOLD/control/null_marker_pgen.txt"

echo "=== controls complete: $GOLD/control"
