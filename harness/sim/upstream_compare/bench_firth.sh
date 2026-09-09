#!/usr/bin/env bash
# Paired timing block, strip vs upstream main, with a Firth-on and a Firth-off arm.
#   bench_firth.sh <reps>
# Protocol as bench.sh (harness/bench/NOTES.md): the two BUILDS run ADJACENTLY within
# each (rep, workload, trait, variant) cell and the build order alternates with the rep,
# so a drift in machine load is differenced out of the strip-vs-main contrast. Refuses
# on a throttled or loaded machine unless FORCE=1.
# Step 2 reads the FROZEN fixture models so both arms test the identical null model.
# Firth-off is run on y_binary only: --is_Firth_beta is documented binary-only and the
# p-value phase checks that claim, so a quantitative Firth-off timing cell would be a
# second measurement of the same thing.
set -uo pipefail
reps=${1:-3}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=/Users/dpalmer/Repositories/SAIGE
export RESULTS=${RESULTS:-$REPO/harness/work/upstream_compare/results_20260909.tsv}
export MODELS=frozen
lim=$(pmset -g therm | awk '/CPU_Speed_Limit/{print $NF}')
load=$(uptime | sed 's/.*load averages*: *//' | awk '{print $1}' | tr -d ,)
echo "CPU_Speed_Limit=$lim load1=$load  $(date '+%F %T')"
if [ "${FORCE:-0}" != 1 ]; then
  if [ "$lim" != 100 ]; then echo "REFUSE: CPU throttled ($lim). FORCE=1 to override." >&2; exit 1; fi
  if awk -v l="$load" 'BEGIN{exit !(l > 2.5)}'; then echo "REFUSE: load $load > 2.5. FORCE=1 to override." >&2; exit 1; fi
fi
run() { # run <build> <wl> <trait> <rep> <variant>
  local v=$5 extra=""
  [ "$v" = firthoff ] && extra="--is_Firth_beta=FALSE"
  VARIANT=$v S2EXTRA=$extra bash "$HERE/run_arm.sh" "$1" "$2" "$3" "$4" 1
}
for r in $(seq 1 "$reps"); do
  if [ $((r % 2)) = 1 ]; then order="strip main"; else order="main strip"; fi
  echo "--- rep $r ($order) $(date '+%F %T') $(uptime) ---"
  for b in $order; do run "$b" step1 y_binary "$r" base; done
  for b in $order; do run "$b" step1 y_quant  "$r" base; done
  for wl in s2single s2group; do
    for t in y_binary y_quant; do
      for b in $order; do run "$b" "$wl" "$t" "$r" base; done
    done
    for b in $order; do run "$b" "$wl" y_binary "$r" firthoff; done
  done
done
echo "done $(date '+%F %T') $(uptime)"
