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
#   s2_marker6k step-2 single-variant scan, chr1 only (~5.8k markers, ONE chunk)
#               -- ~90% process startup; kept for continuity, not for ratios
#   s2_marker_all step-2 single-variant scan, all 128,868 markers in 22 chunks
#               -- the workload a marker-kernel change can be measured on
#   s2_marker_all_sparse  the same scan through the SPARSE scoreTest + PCG
#               (--sparseGRMFile + the binary_sparseGRM_vr model, the only one
#               whose variance-ratio file has both sparse and null_noXadj rows).
#               s2_marker_all cannot measure that path: at default flags every
#               marker goes to scoreTestFast_noadjCov and the PCG never runs.
#               Added with the F8 work on getPCG1ofSigmaAndGtilde.
#   s2_group    step-2 SKAT-O group test (golden P2 config; 100 markers, 2 genes,
#               so also startup-dominated -- see harness/bench/NOTES.md)
#
# Columns beyond maxRSS_MB: startup_s / work_s / throughput (markers/s) are
# derived for s2_* runs from SAIGE's own per-chunk cumulative proc.time();
# pcg_iters is the summed PCG iteration count for s1_* runs, i.e. the GRM-pass
# budget. Rows written before 2026-08-31 carry NA in all four.
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
[ -f "$OUT" ] || printf 'label\tcommit\tworkload\tconfig\trep\telapsed_s\tuser_s\tsys_s\tmaxRSS_MB\tstartup_s\twork_s\tthroughput\tpcg_iters\n' > "$OUT"
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)
SCRATCH=$(mktemp -d)

R2() { pixi run --manifest-path "$REPO/pixi.toml" Rscript "$@"; }

# Split a step-2 run's wall clock into fixed overhead and marker work, using
# only what SAIGE already prints. R/SAIGE_SPATest_Marker.R:141 emits a
# CUMULATIVE proc.time() after every chunk, so with C >= 2 chunks:
#   work_per_chunk = (cum_last - cum_first) / (C - 1)
#   startup_s      = cum_first - work_per_chunk
#   work_s         = cum_last  - startup_s
# This matters because it is the ONLY way to see a marker-kernel change here:
# a whole-genome run at N=1000 is ~2.7 s of startup (R, library(SAIGE),
# ReadModel, bim/fam parse) around ~0.33 s per 5,858-marker chunk, so the
# single-chunk s2_marker6k workload is ~90% overhead and a 2x kernel speedup
# would move its elapsed_s by ~3% -- inside the noise NOTES.md documents.
# C == 1 cannot be decomposed; report NA rather than guess.
s2_breakdown() { # s2_breakdown <logfile> -> "startup work throughput"
  awk '
    /^[[:space:]]*[0-9.]+[[:space:]]+[0-9.]+[[:space:]]+[0-9.]+[[:space:]]*$/ {
      cum = $3; if (n++ == 0) first = cum; last = cum
    }
    / markers were tested/ { markers += $1 }
    END {
      if (n < 2) { print "NA\tNA\tNA"; exit }
      wpc = (last - first) / (n - 1)
      startup = first - wpc
      work = last - startup
      printf "%.3f\t%.3f\t%s\n", startup, work,
             (work > 0 ? sprintf("%.0f", markers / work) : "NA")
    }' "$1"
}

# Step 1 prints an iteration count per PCG solve. Summing them gives the GRM-pass
# budget directly, which is what decides whether blocking GetTrace's nrun probes
# is worth it (each PCG iteration is one full pass over the GRM markers).
s1_pcg_iters() { # s1_pcg_iters <logfile> -> total iterations, or NA
  awk '/iter from getPCG1ofSigmaAndVector/ { s += $NF; n++ }
       END { if (n) print s; else print "NA" }' "$1"
}

timed() { # timed <workload> <config> <rep> <cmd...>
  local wl="$1" cfg="$2" rep="$3"; shift 3
  local tf="$SCRATCH/time.$wl.$cfg.$rep" lf="$SCRATCH/log.$wl.$cfg.$rep"
  /usr/bin/time -l "$@" > "$lf" 2> "$tf" || {
    echo "BENCH RUN FAILED: $wl/$cfg rep $rep"; tail -5 "$lf"; tail -5 "$tf"; exit 1; }
  local real user sys rss startup work thru pcg
  real=$(awk '/real/{print $1}' "$tf" | head -1)
  user=$(awk '/real/{print $3}' "$tf" | head -1)
  sys=$(awk '/real/{print $5}' "$tf" | head -1)
  rss=$(awk '/maximum resident set size/{printf "%.1f", $1/1048576}' "$tf")
  case "$wl" in
    s2_*) IFS=$'\t' read -r startup work thru < <(s2_breakdown "$lf"); pcg="NA" ;;
    *)    startup="NA"; work="NA"; thru="NA"; pcg=$(s1_pcg_iters "$lf") ;;
  esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$LABEL" "$COMMIT" "$wl" "$cfg" "$rep" "$real" "$user" "$sys" "$rss" \
    "$startup" "$work" "$thru" "$pcg" >> "$OUT"
  # Report the number that can actually move: work_s for step 2, elapsed for step 1.
  if [ "$work" != "NA" ]; then
    echo "  $wl/$cfg rep$rep: ${real}s elapsed (${startup}s startup + ${work}s work), ${thru} markers/s, ${rss}MB"
  else
    echo "  $wl/$cfg rep$rep: ${real}s elapsed, ${user}s user, ${pcg} PCG iters, ${rss}MB"
  fi
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

  # The whole genome, not just chr1. This is the workload a step-2 change can
  # actually be measured on: 128,868 markers in 22 chunks, so the ~2.7 s fixed
  # prologue is ~30% rather than ~90% of the run, and s2_breakdown() can
  # separate it out. s2_marker6k is kept above only for continuity with the
  # phase4-* rows already in results.tsv -- do NOT quote it as a kernel ratio.
  timed s2_marker_all default "$rep" \
    pixi run --manifest-path "$REPO/pixi.toml" Rscript step2_SPAtests.R \
    --bedFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bed \
    --bimFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim \
    --famFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.fam \
    --AlleleOrder=alt-first --minMAF=0 --minMAC=1 --LOCO=FALSE \
    --GMMATmodelFile="$GOLD/step1/binary_sparseGRM_a.rda" \
    --varianceRatioFile="$GOLD/step1/binary_sparseGRM_a.varianceRatio.txt" \
    --markers_per_chunk=10000 \
    --SAIGEOutputFile="$SCRATCH/s2_marker_all.txt"

  # Same markers, same fixture, sparse-GRM branch: this is the only workload in
  # which getPCG1ofSigmaAndGtilde and the marker-path scoreTest are on the
  # clock. One PCG solve per marker (~1.3 iterations each at this N).
  timed s2_marker_all_sparse default "$rep" \
    pixi run --manifest-path "$REPO/pixi.toml" Rscript step2_SPAtests.R \
    --bedFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bed \
    --bimFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim \
    --famFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.fam \
    --AlleleOrder=alt-first --minMAF=0 --minMAC=1 --LOCO=FALSE \
    --GMMATmodelFile="$GOLD/step1/binary_sparseGRM_vr.rda" \
    --varianceRatioFile="$GOLD/step1/binary_sparseGRM_vr.varianceRatio.txt" \
    $SG --is_fastTest=FALSE \
    --markers_per_chunk=10000 \
    --SAIGEOutputFile="$SCRATCH/s2_marker_all_sparse.txt"

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
