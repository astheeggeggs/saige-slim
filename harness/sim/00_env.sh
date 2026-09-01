# Interpreters and paths for the synthetic-fixture pipeline.
#
# Replaces the SAIGE_LIB indirection, which pointed into a previous session's
# /private/tmp scratchpad (run_prep_chain.sh:6,17, 07_run_pipeline.sh:7,
# run_full_chain.sh:18) and no longer resolves. SAIGE is installed in the pixi
# env, so Rscript from that env finds it with no library juggling.
#
# The pixi env has NO python -- it is an R build environment. Simulation runs
# on miniconda3 (msprime 1.3.3 / stdpopsim 0.2.0 / tskit 0.5.8), verified.
# Calling the env's Rscript directly rather than via `pixi run` also avoids
# ~2.5-5 s of wrapper startup per invocation (same reasoning as
# harness/run_tests.sh:23-27).
REPO=/Users/dpalmer/Repositories/SAIGE
PY=/Users/dpalmer/miniconda3/bin/python3
P2=$REPO/.pixi/envs/default/bin/plink2
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
D=${D:-$REPO/harness/work/sim50k}

# Determinism for any SAIGE run (see harness/env.sh).
export RCPP_PARALLEL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Keep bench/step-2 scratch on this volume and visible to du, not /var/folders.
export TMPDIR=$REPO/harness/work/tmp
mkdir -p "$TMPDIR"

# Refuse to start a stage without headroom. Disk exhaustion mid-simulation
# wastes the most expensive stage in the pipeline.
require_free_gb() {
  local need=${1:-4} avail
  avail=$(df -g "$REPO" | awk 'NR==2{print $4}')
  if [ "$avail" -lt "$need" ]; then
    echo "ABORT: ${avail} GB free, need ${need} GB" >&2; exit 1
  fi
}
