#!/usr/bin/env bash
# PLINK1 (.bed) vs PGEN (.pgen) single-variant step-2 throughput, single-threaded.
#
# Fair-fight rules, each of which changes the answer if broken:
#  1. --nThreads=1 both arms. The PGEN reader refuses >1 anyway (v1), and
#     mclapply children interleave stdout, which destroys the chunk timings the
#     startup decomposition below depends on.
#  2. Both files pre-warmed into page cache before EVERY rep. The .bed is 6.4x
#     the .pgen on disk (2.63 GB vs 408 MB) on a 16 GB machine, so a cold run
#     measures I/O, not decode. Warm isolates decode; the size gap is reported
#     separately as its own axis.
#  3. Arms ALTERNATE rather than running in blocks, so machine drift hits both.
#  4. Startup is separated, not assumed small. SAIGE prints cumulative
#     proc.time() per chunk; regressing those on chunk index gives per-chunk work
#     as the slope and startup as the intercept. Reported alongside raw wall time
#     so the two can be reconciled rather than trusted.
set -uo pipefail
WT="$(cd "$(dirname "$0")/../.." && pwd)"; cd "$WT"
. harness/sim/00_env.sh
D="${SIMDIR:-$WT/harness/work/sim50k}"
OUT="$WT/harness/work/readerbench"; mkdir -p "$OUT"
REPS="${REPS:-3}"
TRAIT="${TRAIT:-y_quant}"
GRM="$D/pedigree_truth.sparseGRM.mtx"

common=(
  --chrom=20 --minMAF=0 --minMAC=20
  --GMMATmodelFile=$D/step1_$TRAIT.rda
  --varianceRatioFile=$D/step1_$TRAIT.varianceRatio.txt
  --sparseGRMFile=$GRM --sparseGRMSampleIDFile=$GRM.sampleIDs.txt
  --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_fastTest=TRUE
  --nThreads=1
)

warm() { cat "$D/chr20_sim.bed" > /dev/null 2>&1
         cat "$D/chr20_sim_pgen.pgen" > /dev/null 2>&1; }

run_arm() { # run_arm <bed|pgen> <rep>
  # Split, do NOT chain: `local a=$1 b="${a}x"` expands every argument BEFORE the
  # builtin assigns any of them, so ${a} is unbound under set -u. Same trap as
  # `export A=1 B=$A`.
  local arm=$1
  local rep=$2
  local log="$OUT/${arm}_${rep}.log"
  local res="$OUT/${arm}_${rep}.txt"
  local -a io
  if [ "$arm" = bed ]; then
    io=(--bedFile=$D/chr20_sim.bed --bimFile=$D/chr20_sim.bim --famFile=$D/chr20_sim.fam)
  else
    io=(--pgenFile=$D/chr20_sim_pgen.pgen --pvarFile=$D/chr20_sim_pgen.pvar --psamFile=$D/chr20_sim_pgen.psam)
  fi
  warm
  /usr/bin/time -p $RS extdata/step2_SPAtests.R "${io[@]}" "${common[@]}" \
      --SAIGEOutputFile="$res" > "$log" 2>&1
  local rc=$?
  # /usr/bin/time -p writes real/user/sys to stderr, which is inside $log
  local real user
  real=$(awk '/^real/{print $2}' "$log" | tail -1)
  user=$(awk '/^user/{print $2}' "$log" | tail -1)
  echo -e "$arm\t$rep\t$rc\t$real\t$user"
}

echo -e "arm\trep\trc\treal_s\tuser_s"
for r in $(seq 1 $REPS); do
  run_arm bed  $r
  run_arm pgen $r
done
