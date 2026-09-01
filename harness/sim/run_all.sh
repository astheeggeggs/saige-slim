#!/usr/bin/env bash
# Regenerate the whole fixture and run SAIGE steps 0/1/2.
#
#   bash harness/sim/run_all.sh                      # CONFIG=chr20 (default)
#   CONFIG=bench50k bash harness/sim/run_all.sh      # step-2 benchmark fixture
#   SIMDIR=... RIGHT=2000000 NEUR=600 NFAM=50 NSING=400 NGENES=40 \
#     bash harness/sim/run_all.sh                    # small smoke run
#
# ONE generator, two configurations -- they differ only in these variables, so
# a generator fix lands in both at once and they cannot drift.
#
#   chr20     correctness fixture. Whole chr20, N=10,000. ~80 min, ~7 GB.
#             This is the one whose controls are quoted in the plan.
#   bench50k  step-2 benchmark fixture. chr20:34-59 Mb, N=50,000. The N=1,000
#             bundled fixture INVERTS the ranking it is meant to measure --
#             writeOutfile_single (22% of the marker loop) and exceptions
#             (8.8%) are flat in N, while decode/impute/score/SPA all scale
#             with it. Window and thinning come from the measured pilot
#             (harness/sim/PILOT.md): mean_r 0.98x the chromosome mean,
#             ~19 s/Mb, 72,178 kept sites/Mb at 46,250 founders.
#
# Anything already set in the environment wins; a preset only fills in blanks.
#
# RELCUT defaults to 0.125, chosen by measurement -- see README.md. Set
# RELCUT=0.05 to reproduce the diagnostic that shows why the recipe default
# is unusable on a single chromosome.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO"

CONFIG="${CONFIG:-chr20}"
case "$CONFIG" in
  chr20)
    : "${SIMDIR:=$REPO/harness/work/sim}" "${NEED_GB:=8}"
    : "${NEUR:=9000}" "${NSING:=8000}" "${NFAM:=500}" "${NGENES:=550}"
    ;;
  bench50k)
    : "${SIMDIR:=$REPO/harness/work/sim50k}" "${NEED_GB:=14}"
    : "${LEFT:=34000000}" "${RIGHT:=59000000}"
    # 42,500 + 4*1,875 = 50,000 output samples from 42,500 + 2*1,875 = 46,250
    # founders, which is the founder count the pilot's site counts were taken at.
    : "${NEUR:=46250}" "${NSING:=42500}" "${NFAM:=1875}" "${NGENES:=210}"
    # Rare INTERGENIC variants only; genic and MAF >= 1% are never thinned, so
    # the group-test, GRM and variance-ratio pools are untouched. 0.070 targets
    # M ~ 200,000 (.bed ~ 2.50 GB) from the pilot's 72,178 kept sites/Mb.
    : "${THIN:=0.070}"
    ;;
  *) echo "unknown CONFIG=$CONFIG (want chr20 or bench50k)" >&2; exit 1 ;;
esac
export D="$SIMDIR"
. harness/sim/00_env.sh
mkdir -p "$D"
export RELCUT="${RELCUT:-0.125}" NRANDMARK="${NRANDMARK:-3000}" NTHREADS="${NTHREADS:-6}"
echo "CONFIG=$CONFIG  D=$D  N=$(( NSING + 4 * NFAM ))  founders=$(( NSING + 2 * NFAM ))"
echo "  window=${LEFT:-0}-${RIGHT:-end}  n_eur=$NEUR  n_genes=$NGENES" \
     "thin_intergenic=${THIN:-1.0}  relcut=$RELCUT  need=${NEED_GB}GB"
# The full run is ~80 min; DRYRUN resolves the configuration and stops, which is
# the only cheap way to check a preset before committing to it.
if [ -n "${DRYRUN:-}" ]; then echo "DRYRUN: stopping before stage 01"; exit 0; fi
S=harness/sim
step() { echo; echo "### $*"; }

require_free_gb "${NEED_GB:-8}"

step "01 coalescent simulation"
$PY $S/01_simulate.py --out "$D/chr20.trees" \
    ${NEUR:+--n-eur $NEUR} ${LEFT:+--left $LEFT} ${RIGHT:+--right $RIGHT} \
    2>&1 | tee "$D/01_simulate.log"

step "02 genotypes, annotations, traits"
$PY $S/02_build_dataset.py --trees "$D/chr20.trees" --outdir "$D" \
    ${NSING:+--n-singletons $NSING} ${NFAM:+--n-families $NFAM} \
    ${NGENES:+--n-genes $NGENES} ${THIN:+--thin-intergenic $THIN} \
    2>&1 | tee "$D/02_build.log"

step "03 plink/plink2 sets, PCs"
bash $S/03_plink_prep.sh 2>&1 | tee "$D/03_plink.log"
$PY $S/04_add_pcs.py "$D" 4

# A GRM estimated from ONE 25 Mb window is not usable, and the reason is
# genetic rather than a matter of counting markers: 25 Mb of chr20 is ~41 cM on
# the HapMap map (108.1 cM over 64.44 Mb, and the window was chosen at 0.98x the
# chromosome mean rate), so the number of effectively independent common markers
# is a few hundred HOWEVER MANY markers the window contains. That is the
# measured failure in the prior run -- true first-degree pairs sat at median
# 0.493 against an expected 0.5, but "unrelated" pairs sat at median 0.075 with
# the diagonal ranging 0.228-2.230, i.e. the estimator was fine and the input
# was not. Population structure is ruled out (02_build_dataset.py keeps EUR
# only). The fix is the synthetic independent scaffold -- founder frequencies
# drawn on [0.05,0.45] and gene-dropped through this same pedigree -- which is
# not built yet. Until it is, refuse rather than emit a GRM that looks fine and
# is not; <prefix>.pedigree_truth.txt holds the exact answer to check against.
if [ "$CONFIG" = "bench50k" ] && [ -z "${ALLOW_REGION_GRM:-}" ]; then
  echo
  echo "STOP: CONFIG=bench50k has no usable GRM input yet (25 Mb ~ 41 cM; see"
  echo "      the comment above and the GRM section of the plan). Stages 01-03"
  echo "      are complete in $D."
  echo "      Set ALLOW_REGION_GRM=1 to run step 0 anyway as a DIAGNOSTIC --"
  echo "      its output is not fit to fit a null model against."
  exit 0
fi

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
