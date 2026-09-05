#!/usr/bin/env bash
# Paired timing block, strip vs upstream main, on sim50k.
#   bench.sh <reps> [workloads] [traits]
# Arms run ADJACENTLY within each (rep, workload, trait) and the arm ORDER alternates
# with the rep (odd: strip then main; even: main then strip) -- harness/bench/NOTES.md's
# protocol. Refuses on a throttled or loaded machine unless FORCE=1.
# Step 2 runs read the FROZEN fixture models (MODELS=frozen) so the two arms test
# the identical null model; the step1 workload writes each build's own models.
set -uo pipefail
reps=${1:-3}; wls=${2:-"step1 s2single s2group"}; traits=${3:-"y_binary y_quant"}
HERE=$(cd "$(dirname "$0")" && pwd)
lim=$(pmset -g therm | awk '/CPU_Speed_Limit/{print $NF}')
load=$(uptime | sed 's/.*load averages*: *//' | awk '{print $1}' | tr -d ,)
echo "CPU_Speed_Limit=$lim load1=$load  $(date '+%F %T')"
if [ "${FORCE:-0}" != 1 ]; then
  if [ "$lim" != 100 ]; then echo "REFUSE: CPU throttled ($lim). FORCE=1 to override." >&2; exit 1; fi
  if awk -v l="$load" 'BEGIN{exit !(l > 2.5)}'; then echo "REFUSE: load $load > 2.5. FORCE=1 to override." >&2; exit 1; fi
fi
for r in $(seq 1 "$reps"); do
  if [ $((r % 2)) = 1 ]; then order="strip main"; else order="main strip"; fi
  for wl in $wls; do
    for t in $traits; do
      for b in $order; do
        MODELS=${MODELS:-frozen} bash "$HERE/run_arm.sh" "$b" "$wl" "$t" "$r" "${NT:-1}"
      done
    done
  done
done
echo "done $(date '+%F %T') $(uptime)"
