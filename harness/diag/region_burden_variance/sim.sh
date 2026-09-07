#!/usr/bin/env bash
# sim.sh <label> <n> <beta_b> <full|s1> [parallel] [tau2]: n replicates of
# sim_rep.sh. tau2 (default: the model's) scales the simulated random effects.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="${DIAG_OUT:-$REPO/harness/work/diag-burden}"
export DIAG_OUT="$OUT" SAIGE_REPO="$REPO"
label=$1; n=$2; beta=$3; mode=$4; par=${5:-4}; tau2=${6:-}
mkdir -p "$OUT/sim/$label"
if [ -z "$tau2" ]; then prep="$OUT/sim/prep.rds"; else prep="$OUT/sim/prep_tau$tau2.rds"; fi
if [ ! -f "$prep" ]; then
  (cd "$REPO" && pixi run Rscript "$HERE/sim_prep.R" "$prep" $tau2) || exit 1
fi
export SIM_PREP="$prep"
seq 1 "$n" | xargs -P "$par" -I{} "$HERE/sim_rep.sh" "$label" {} "$beta" "$mode"
