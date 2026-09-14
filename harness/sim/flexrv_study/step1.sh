#!/bin/bash
# FlexRV simulation study, S1: the four null models on the REBUILT fixture, so
# the model set is consistent with its phenotypes and PCs. y_null_bin is the
# 2026-09-08 recipe (harness/work/null50k-20260908/binnull_step1.sh): y_null
# thresholded at y_binary's own prevalence, a binary null with the same
# pedigree background and zero chr20 heritability; permuting is ruled out
# (LEDGER #104: it drives tau2 to 0). Writes $D/step1_<trait>.* and the
# phenotype file with the y_null_bin column beside the fixture's.
# usage: step1.sh [trait ...]   (default: all four; env D, NT)
set -euo pipefail
REPO=/Users/dpalmer/Repositories/SAIGE
D=${D:-/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k}
NT=${NT:-4}
cd "$REPO"; . harness/sim/00_env.sh
GRM=$D/pedigree_truth.sparseGRM.mtx
PH=$D/pheno_pcs.txt
PHB=$D/pheno_pcs_binnull.txt
if [ ! -s "$PHB" ] || [ "$PH" -nt "$PHB" ]; then
  $RS - <<EOF
p <- read.table("$PH", header = TRUE, stringsAsFactors = FALSE)
K <- mean(p\$y_binary)
p\$y_null_bin <- as.integer(p\$y_null > quantile(p\$y_null, 1 - K))
cat("y_null_bin: K =", K, " cases =", sum(p\$y_null_bin), "\n")
write.table(p, "$PHB", quote = FALSE, row.names = FALSE)
EOF
fi
fit() {
  local col=$1 tt=$2 ph=$3 extra=""
  [ "$tt" = quantitative ] && extra="--invNormalize=TRUE"
  echo "=== step 1: $col ($tt) ==="
  $RS extdata/step1_fitNULLGLMM.R \
    --plinkFile=$D/data_for_vr \
    --sparseGRMFile=$GRM --sparseGRMSampleIDFile=$GRM.sampleIDs.txt \
    --useSparseGRMtoFitNULL=TRUE --useSparseGRMforVarRatio=TRUE \
    --phenoFile=$ph --phenoCol=$col --sampleIDColinphenoFile=IID \
    --traitType=$tt $extra \
    --covarColList=age,sex,PC1,PC2,PC3,PC4 --qCovarColList=sex \
    --skipVarianceRatioEstimation=FALSE --isCateVarianceRatio=TRUE \
    --relatednessCutoff=0.125 --IsOverwriteVarianceRatioFile=TRUE \
    --nThreads=$NT --outputPrefix=$D/step1_$col > $D/s1_$col.flexrv.log 2>&1
  grep -m1 -i 'theta\|tau' $D/s1_$col.flexrv.log || true
  echo "$col: $(cat $D/step1_$col.varianceRatio.txt | wc -l) VR rows"
}
# The second-generation power phenotypes (make_power_pheno.py, 2026-09-14)
# live in pheno_pcs<suffix>.txt with their own y_null_bin<suffix> column, so a
# trait carrying a suffix is fitted straight from that file.
for t in ${@:-y_null y_null_bin y_quant y_binary}; do
  case $t in
    y_null)     fit y_null quantitative $PH ;;
    y_null_bin) fit y_null_bin binary $PHB ;;
    y_quant)    fit y_quant quantitative $PH ;;
    y_binary)   fit y_binary binary $PH ;;
    y_quant_*)  fit $t quantitative $D/pheno_pcs${t#y_quant}.txt ;;
    y_binary_*) fit $t binary $D/pheno_pcs${t#y_binary}.txt ;;
    y_null_bin_*) fit $t binary $D/pheno_pcs${t#y_null_bin}.txt ;;
    y_null_r*)  fit $t quantitative $D/pheno_pcs_null.txt ;;   # make_null_series.py replicates
    *) echo "unknown trait $t" >&2; exit 2 ;;
  esac
done
