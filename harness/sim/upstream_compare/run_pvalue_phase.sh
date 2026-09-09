#!/usr/bin/env bash
# P-value phase: one rep of every (build x workload x trait x firth) cell, MODELS=own
# (each build reads the step-1 model its OWN step 1 wrote -- the pipeline as a user runs
# it). Load-insensitive: outputs are deterministic, only the timings in the TSV are not.
# The paired timing block is bench_firth.sh, run separately on a quiet machine.
set -uo pipefail
REPO=/Users/dpalmer/Repositories/SAIGE
CMP=$REPO/harness/sim/upstream_compare
export RESULTS=$REPO/harness/work/upstream_compare/results_20260909.tsv
export MODELS=own
echo "=== step 1 (own models) $(date '+%F %T') ==="
for b in strip main; do
  for t in y_binary y_quant; do
    VARIANT=base bash "$CMP/run_arm.sh" "$b" step1 "$t" 0 1 || echo "FAIL step1 $b $t"
  done
done
for wl in s2single s2group; do
  for t in y_binary y_quant; do
    for b in strip main; do
      echo "=== $wl $t $b firth-on $(date '+%F %T') ==="
      VARIANT=base bash "$CMP/run_arm.sh" "$b" "$wl" "$t" 0 1 || echo "FAIL $wl $t $b base"
      echo "=== $wl $t $b firth-off $(date '+%F %T') ==="
      VARIANT=firthoff S2EXTRA=--is_Firth_beta=FALSE \
        bash "$CMP/run_arm.sh" "$b" "$wl" "$t" 0 1 || echo "FAIL $wl $t $b firthoff"
    done
  done
done
echo "=== done $(date '+%F %T') ==="
