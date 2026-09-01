#!/usr/bin/env bash
# step 1 + step 2 for every trait. Each stage is logged separately and a
# failure in one trait does not abort the rest.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
export D="${SIMDIR:-${D:-$REPO/harness/work/sim}}"
. harness/sim/00_env.sh
export NTHREADS="${NTHREADS:-6}"
export RELCUT="${RELCUT:-0.05}"
R=harness/sim/04_run_saige.sh
TRAITS=("y_quant quantitative" "y_binary binary" "y_quant_pc quantitative"
        "y_binary_pc binary" "y_null quantitative")

run() {  # run <logname> <args...>
  local log=$1; shift
  local t0=$SECONDS
  if bash $R "$@" > "$D/$log.log" 2>&1; then
    echo "  OK   $log  ($((SECONDS - t0))s)"
  else
    echo "  FAIL $log  ($((SECONDS - t0))s)  -- tail:"; tail -5 "$D/$log.log" | sed 's/^/       /'
  fi
}

echo "=== STEP 1 (relatednessCutoff=$RELCUT) ==="
for t in "${TRAITS[@]}"; do set -- $t; run "s1_$1" step1 "$1" "$2"; done
echo "=== STEP 2 group (SAIGE-GENE+) ==="
for t in "${TRAITS[@]}"; do set -- $t; run "s2g_$1" step2_group "$1"; done
echo "=== STEP 2 single (SAIGE-variant) ==="
for t in "${TRAITS[@]}"; do set -- $t; run "s2s_$1" step2_single "$1"; done
echo "=== pipeline done ==="
