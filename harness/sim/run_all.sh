#!/usr/bin/env bash
# Regenerate the whole fixture and run SAIGE steps 0/1/2. ~80 min, ~7 GB.
#
#   bash harness/sim/run_all.sh                      # full N=10,000 chr20
#   SIMDIR=... RIGHT=2000000 NEUR=600 NFAM=50 NSING=400 NGENES=40 \
#     bash harness/sim/run_all.sh                    # small smoke run
#
# RELCUT defaults to 0.125, chosen by measurement -- see README.md. Set
# RELCUT=0.05 to reproduce the diagnostic that shows why the recipe default
# is unusable on a single chromosome.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"
export D="${SIMDIR:-$REPO/harness/work/sim}"
. harness/sim/00_env.sh
mkdir -p "$D"
export RELCUT="${RELCUT:-0.125}" NRANDMARK="${NRANDMARK:-3000}" NTHREADS="${NTHREADS:-6}"
S=harness/sim
step() { echo; echo "### $*"; }

require_free_gb "${NEED_GB:-8}"

step "01 coalescent simulation"
$PY $S/01_simulate.py --out "$D/chr20.trees" \
    ${NEUR:+--n-eur $NEUR} ${RIGHT:+--right $RIGHT} 2>&1 | tee "$D/01_simulate.log"

step "02 genotypes, annotations, traits"
$PY $S/02_build_dataset.py --trees "$D/chr20.trees" --outdir "$D" \
    ${NSING:+--n-singletons $NSING} ${NFAM:+--n-families $NFAM} \
    ${NGENES:+--n-genes $NGENES} 2>&1 | tee "$D/02_build.log"

step "03 plink/plink2 sets, PCs"
bash $S/03_plink_prep.sh 2>&1 | tee "$D/03_plink.log"
$PY $S/04_add_pcs.py "$D" 4

step "04 step 0: sparse GRM (relatednessCutoff=$RELCUT)"
# createSparseGRM refuses to start if numRandomMarkerforSparseKin exceeds the
# markers passing its own MAF filter, so cap it against what pruning actually
# left (relevant at small --right scales; the full run has 3,055).
AVAIL=$(wc -l < "$D/data_for_grm.bim")
CAP=$(( AVAIL * 85 / 100 ))
if [ "$NRANDMARK" -gt "$CAP" ]; then
  echo "note: only $AVAIL LD-pruned markers; lowering numRandomMarkerforSparseKin $NRANDMARK -> $CAP"
  export NRANDMARK=$CAP
fi
if ! bash $S/04_run_saige.sh step0 > "$D/step0.log" 2>&1; then
  echo "step 0 FAILED:"; tail -12 "$D/step0.log"; exit 1
fi
$RS $S/05_check_grm.R \
  "$D/sparseGRM_relatednessCutoff_${RELCUT}_${NRANDMARK}_randomMarkersUsed.sparseGRM.mtx" \
  2>&1 | tee "$D/grm_check_$RELCUT.txt"

step "05 steps 1 and 2, all traits"
bash $S/07_run_pipeline.sh

step "06 results vs simulated truth"
$RS $S/06_analyse.R 2>&1 | tee "$D/RESULTS.txt"
echo; echo "### done -- outputs in $D"
