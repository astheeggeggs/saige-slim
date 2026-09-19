#!/bin/bash
# Region-level parallelism: does splitting the group file across PROCESSES beat
# threading inside one process? BURDEN only -- SKAT-O already has cross-region
# parallelism (step 5, 1c706bd2).
#
# USAGE:  bash bench.sh <fixture-dir> [reps]
#
# The fixture dir must hold: geno.{pgen,pvar,psam}, groupFile.txt,
# step1_y_binary.{rda,varianceRatio.txt}, pedigree_truth.sparseGRM.mtx{,.sampleIDs.txt}
# See README.md beside this script for how to build it (~28 MB).
#
# WHY THIS REFUSES TO RUN ON A BUSY BOX: the first version of this measurement
# (2026-09-19) was taken on a 4-core laptop at load ~3 with an editor and a
# browser running. The serial baseline got 0.84 of a core, the 4-thread arm 1.89
# and the 4-process arm 3.24 -- so nothing measured scaling, only how
# aggressively each arm grabbed contended CPU, and it reported the wrong winner.
# USER TIME is the load-insensitive number and is what this reports first.
set -uo pipefail
FIX=${1:?usage: bench.sh <fixture-dir> [reps]}
REPS=${2:-3}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT

ncores=$( sysctl -n hw.physicalcpu 2>/dev/null || nproc 2>/dev/null || echo 0 )
load=$(uptime | sed 's/.*load averages*: *//' | awk '{print $1}' | tr -d ,)
echo "cores=$ncores  1-min load=$load  reps=$REPS"
if [ "${BENCH_FORCE:-0}" != "1" ]; then
  busy=$(awk -v l="$load" 'BEGIN{print (l > 0.5) ? 1 : 0}')
  [ "$busy" = "1" ] && { echo "REFUSING: 1-min load $load > 0.5. A 4-way scaling test needs 4 idle cores."
                         echo "Close other work and retry, or set BENCH_FORCE=1 to override (results will not be a scaling measurement)."; exit 2; }
  [ "$ncores" -lt 4 ] && { echo "REFUSING: only $ncores physical cores; a 4-way test cannot scale here."; exit 2; }
fi

cd "$REPO" && source harness/env.sh && cd extdata
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
COMMON="--pgenFile=$FIX/geno.pgen --pvarFile=$FIX/geno.pvar --psamFile=$FIX/geno.psam \
 --chrom=20 --minMAF=0 --minMAC=0.5 \
 --GMMATmodelFile=$FIX/step1_y_binary.rda --varianceRatioFile=$FIX/step1_y_binary.varianceRatio.txt \
 --sparseGRMFile=$FIX/pedigree_truth.sparseGRM.mtx \
 --sparseGRMSampleIDFile=$FIX/pedigree_truth.sparseGRM.mtx.sampleIDs.txt \
 --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_output_moreDetails=TRUE \
 --is_fastTest=FALSE --is_single_in_groupTest=TRUE --r.corr=1 \
 --maxMAF_in_groupTest=0.0001,0.001,0.01 \
 --annotation_in_groupTest=high,medium,low,high:medium,high:medium:low"

ngenes=$(awk '$2=="var"' "$FIX/groupFile.txt" | wc -l | tr -d ' ')
awk -v W="$W" -v n="$ngenes" 'NR%2==1{g++} {q=int((g-1)*4/n)+1; if(q>4)q=4; print > (W"/q"q".txt")}' "$FIX/groupFile.txt"

# /usr/bin/time -p is POSIX: writes "real/user/sys" to stderr on macOS and Linux
run1() { /usr/bin/time -p $RS step2_SPAtests.R $COMMON "$@" 2>&1; }
field() { grep -E "^$1 " "$2" | tail -1 | awk '{print $2}'; }

for rep in $(seq 1 "$REPS"); do
  run1 --groupFile="$FIX/groupFile.txt" --nThreads=1 --SAIGEOutputFile=$W/s.txt > $W/s.log
  run1 --groupFile="$FIX/groupFile.txt" --nThreads=4 --SAIGEOutputFile=$W/t.txt > $W/t.log
  st=$(date +%s.%N)
  for q in 1 2 3 4; do run1 --groupFile=$W/q$q.txt --nThreads=1 --SAIGEOutputFile=$W/p$q.txt > $W/p$q.log & done
  wait; en=$(date +%s.%N)
  pu=$(for q in 1 2 3 4; do field user $W/p$q.log; done | awk '{s+=$1} END{print s}')
  printf "rep%s  serial  wall %6.2f user %6.2f | threads4 wall %6.2f user %6.2f | procs4 wall %6.2f user %6.2f\n" \
    "$rep" "$(field real $W/s.log)" "$(field user $W/s.log)" \
    "$(field real $W/t.log)" "$(field user $W/t.log)" \
    "$(echo "$en - $st" | bc)" "$pu"
  # correctness: the four quarters must reconstruct the whole-file result
  rs=$(tail -n +2 $W/s.txt 2>/dev/null | wc -l | tr -d ' ')
  cat $W/p?.txt 2>/dev/null | grep -v '^Region' | sort > $W/m.s
  tail -n +2 $W/s.txt 2>/dev/null | sort > $W/w.s
  rm_=$(wc -l < $W/m.s | tr -d ' ')
  if [ "$rs" -gt 0 ] && [ "$rs" = "$rm_" ] && cmp -s $W/w.s $W/m.s; then echo "      rows $rs, quarters reconstruct the whole run: IDENTICAL"
  else echo "      ROWS whole=$rs quarters=$rm_ -- MISMATCH or empty; the timings above are NOT valid"; fi
done
