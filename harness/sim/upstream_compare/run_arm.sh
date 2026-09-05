#!/usr/bin/env bash
# One SAIGE run of one build on the sim50k fixture, timed.
#
#   run_arm.sh <build> <workload> <trait> <rep> [nthreads]
#     build     strip | main      (main = upstream 365deeff, worktree $MAIN, library $MAIN/Rlib)
#     workload  step1 | s2single | s2group
#     trait     y_quant y_binary y_quant_pc y_binary_pc y_null
#     nthreads  default 1 (upstream's --nThreads is a no-op, LEDGER #30; >1 is strip-only)
#
# Env: MODELS=frozen|own (default frozen: both arms read the fixture's committed strip
# step-1 models; own: each arm reads the step-1 files ITS OWN step1 run wrote under
# harness/work/upstream_compare/<build>/models). RESULTS overrides the TSV path.
#
# The flags come from harness/sim/04_run_saige.sh unchanged; this script only points it
# at a build (RS, R_LIBS_USER, SAIGE_EXTDATA) and at prefixes outside the fixture dir.
# Every run records the SAIGE library the R process actually loaded (the null control
# from harness/bench/NOTES.md: an arm that silently loads the other build is the
# failure mode this catches).
set -uo pipefail
build=$1; wl=$2; trait=$3; rep=$4; nt=${5:-1}
REPO=/Users/dpalmer/Repositories/SAIGE
MAIN=${MAIN:-/Users/dpalmer/Repositories/SAIGE-main}
SIM=${SIM:-/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k}
WORK=$REPO/harness/work/upstream_compare/$build
RESULTS=${RESULTS:-$REPO/harness/work/upstream_compare/results.tsv}
MODELS=${MODELS:-frozen}
mkdir -p "$WORK/models" "$WORK/logs" "$WORK/out"

case $trait in
  y_binary|y_binary_pc) tt=binary ;;
  y_quant|y_quant_pc|y_null) tt=quantitative ;;
  *) echo "unknown trait $trait" >&2; exit 2 ;;
esac
case $build in
  strip) export RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript; TREE=$REPO ;;
  main)  export RS=$MAIN/.pixi/envs/default/lib/R/bin/Rscript; TREE=$MAIN
         export R_LIBS_USER=$MAIN/Rlib ;;
  *) echo "unknown build $build" >&2; exit 2 ;;
esac
export SAIGE_EXTDATA=$TREE/extdata
export SIMDIR=$SIM D=$SIM
export GRMFILE=$SIM/pedigree_truth.sparseGRM.mtx GRMIDS=$SIM/pedigree_truth.sparseGRM.mtx.sampleIDs.txt
export RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export NTHREADS=$nt GROUP_NT=$nt
export RELCUT=0.125
tag=${wl}_${trait}_r${rep}_nt${nt}_${MODELS}
log=$WORK/logs/$tag.log

case $MODELS in
  frozen) export MODELPREFIX=$SIM/step1_$trait ;;
  own)    export MODELPREFIX=$WORK/models/step1_$trait ;;
  *) echo "MODELS must be frozen|own" >&2; exit 2 ;;
esac
case $wl in
  step1)    cmd=(step1 "$trait" "$tt"); export S1OUTPREFIX=$WORK/models/step1_$trait
            warm=("$SIM/data_for_vr.bed") ;;
  s2single) cmd=(step2_single "$trait"); export S2OUT=$WORK/out/$tag.txt
            warm=("$SIM/chr20_sim.bed" "$MODELPREFIX.rda") ;;
  s2group)  cmd=(step2_group "$trait");  export S2OUT=$WORK/out/$tag.txt
            warm=("$SIM/chr20_sim.bed" "$MODELPREFIX.rda") ;;
  *) echo "unknown workload $wl" >&2; exit 2 ;;
esac
if [ "$wl" != step1 ] && [ ! -f "$MODELPREFIX.rda" ]; then
  echo "missing model $MODELPREFIX.rda" >&2; exit 3
fi
[ -n "${S2OUT:-}" ] && rm -f "$S2OUT" "$S2OUT".index

lib=$($RS -e 'cat(find.package("SAIGE"))' 2>/dev/null)
so=$(ls "$lib"/libs/SAIGE.so 2>/dev/null)
somd5=$(md5 -q "$so" 2>/dev/null)
for f in "${warm[@]}"; do cat "$f" > /dev/null 2>&1; done
{
  echo "# build=$build workload=$wl trait=$trait rep=$rep nthreads=$nt models=$MODELS"
  echo "# lib=$lib so_md5=$somd5 MODELPREFIX=$MODELPREFIX"
  echo "# start $(date '+%F %T') $(uptime) $(pmset -g therm | grep CPU_Speed_Limit | tr -s ' \t' ' ')"
} > "$log"
/usr/bin/time -l bash "$REPO/harness/sim/04_run_saige.sh" "${cmd[@]}" >> "$log" 2>&1
rc=$?
echo "# end $(date '+%F %T') $(uptime)" >> "$log"
read -r real user sys rss < <(awk '/ real .* user .* sys/{r=$1;u=$3;s=$5} /maximum resident set size/{m=$1} END{print r,u,s,m/1048576}' "$log")
[ -f "$RESULTS" ] || printf 'build\tworkload\ttrait\trep\tnthreads\tmodels\trc\treal_s\tuser_s\tsys_s\tmaxRSS_MB\tlib\tso_md5\tlog\n' > "$RESULTS"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%s\n' \
  "$build" "$wl" "$trait" "$rep" "$nt" "$MODELS" "$rc" "$real" "$user" "$sys" "$rss" "$lib" "$somd5" "$log" | tee -a "$RESULTS"
exit $rc
