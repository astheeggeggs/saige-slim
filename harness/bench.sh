#!/usr/bin/env bash
# Phase 4 benchmark harness. Times representative workloads and appends one
# TSV row per run to harness/bench/results.tsv (commit, workload, config,
# rep, elapsed_s, user_s, sys_s, maxRSS_MB).
#
#   bash harness/bench.sh <label> [reps]
#
# Workloads:
#   s1_sparse   step-1 binary fit, sparse GRM (golden config a)
#   s1_full     step-1 binary fit, FULL dense GRM (the heavy TBB-kernel path),
#               run at nThreads=1 and nThreads=4
#   s2_marker6k step-2 single-variant scan, ~5.8k markers (chunked config)
#   s2_group    step-2 SKAT-O group test (golden P2 config)
#
# Timing discipline: run on an otherwise-quiet machine; env pinned except the
# s1_full nThreads=4 config (which sets its own thread count via LEDGER-#7's
# fixed unconditional setThreadOptions).

set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="${1:?usage: bench.sh <label> [reps]}"
REPS="${2:-2}"
cd "$REPO/extdata"
source "$REPO/harness/env.sh"
mkdir -p "$REPO/harness/bench"
OUT="$REPO/harness/bench/results.tsv"
[ -f "$OUT" ] || echo -e "label\tcommit\tworkload\tconfig\trep\telapsed_s\tuser_s\tsys_s\tmaxRSS_MB" > "$OUT"
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)
SCRATCH=$(mktemp -d)

R2() { pixi run --manifest-path "$REPO/pixi.toml" Rscript "$@"; }

timed() { # timed <workload> <config> <rep> <cmd...>
  local wl="$1" cfg="$2" rep="$3"; shift 3
  local tf="$SCRATCH/time.$wl.$cfg.$rep"
  /usr/bin/time -l "$@" > "$SCRATCH/log.$wl.$cfg.$rep" 2> "$tf" || {
    echo "BENCH RUN FAILED: $wl/$cfg rep $rep"; tail -5 "$SCRATCH/log.$wl.$cfg.$rep"; tail -5 "$tf"; exit 1; }
  local real user sys rss
  real=$(awk '/real/{print $1}' "$tf" | head -1)
  user=$(awk '/real/{print $3}' "$tf" | head -1)
  sys=$(awk '/real/{print $5}' "$tf" | head -1)
  rss=$(awk '/maximum resident set size/{printf "%.1f", $1/1048576}' "$tf")
  echo -e "$LABEL\t$COMMIT\t$wl\t$cfg\t$rep\t$real\t$user\t$sys\t$rss" >> "$OUT"
  echo "  $wl/$cfg rep$rep: ${real}s elapsed, ${user}s user, ${rss}MB"
}

GOLD="$REPO/harness/goldens"
SG="--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

for rep in $(seq 1 "$REPS"); do
  echo "=== rep $rep"

  timed s1_sparse pinned "$rep" env RCPP_PARALLEL_NUM_THREADS=1 \
    pixi run --manifest-path "$REPO/pixi.toml" Rscript step1_fitNULLGLMM.R \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
    $SG --useSparseGRMtoFitNULL=TRUE \
    --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
    --phenoCol=y_binary --covarColList=x1,x2,a9,a10 --qCovarColList=a9,a10 \
    --sampleIDColinphenoFile=IID --traitType=binary \
    --IsOverwriteVarianceRatioFile=TRUE --isCateVarianceRatio=TRUE \
    --outputPrefix="$SCRATCH/s1_sparse"

  for nt in 1 4; do
    timed s1_full "nt$nt" "$rep" env RCPP_PARALLEL_NUM_THREADS=$nt \
      pixi run --manifest-path "$REPO/pixi.toml" Rscript step1_fitNULLGLMM.R \
      --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
      --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
      --phenoCol=y_binary --covarColList=x1,x2,a9,a10 --qCovarColList=a9,a10 \
      --sampleIDColinphenoFile=IID --traitType=binary --nThreads=$nt \
      --LOCO=FALSE --minMAFforGRM=0.01 \
      --IsOverwriteVarianceRatioFile=TRUE \
      --outputPrefix="$SCRATCH/s1_full_nt$nt"
  done

  timed s2_marker6k default "$rep" \
    pixi run --manifest-path "$REPO/pixi.toml" Rscript step2_SPAtests.R \
    --bedFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bed \
    --bimFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim \
    --famFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.fam \
    --AlleleOrder=alt-first --chrom=1 --minMAF=0 --minMAC=1 --LOCO=FALSE \
    --GMMATmodelFile="$GOLD/step1/binary_sparseGRM_a.rda" \
    --varianceRatioFile="$GOLD/step1/binary_sparseGRM_a.varianceRatio.txt" \
    --markers_per_chunk=1000 \
    --SAIGEOutputFile="$SCRATCH/s2_marker.txt"

  timed s2_group default "$rep" \
    pixi run --manifest-path "$REPO/pixi.toml" Rscript step2_SPAtests.R \
    --bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim \
    --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first \
    --chrom=1 --minMAF=0 --minMAC=0.5 --LOCO=FALSE \
    --GMMATmodelFile=./output/example_binary_sparseGRM.rda \
    --varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt \
    $SG --groupFile=./input/group_new_snpid.txt \
    --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 --is_fastTest=TRUE \
    --SAIGEOutputFile="$SCRATCH/s2_group.txt"
done

rm -rf "$SCRATCH"
echo "=== bench done; results in $OUT"
