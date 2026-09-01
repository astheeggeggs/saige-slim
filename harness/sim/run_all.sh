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
    # Trait scaling, per amendment A4. A fixed h2 gives every causal gene the
    # same power, which at N=50,000 saturates the positive control (p ~ 1e-25)
    # exactly where BURDEN/SKAT/SKAT-O most need to be told apart. The ladder
    # spreads genes across the power curve so the control becomes
    # Spearman(NCP, -log10 p). NCP 1 is undetectable and 60 is p ~ 1e-14, which
    # brackets the interesting range.
    : "${PREVALENCE:=0.01}" "${NCP_LADDER:=1,60}" "${H2_BINARY_MULT:=auto}"
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
echo "  prevalence=${PREVALENCE:-0.10}  ncp_ladder=${NCP_LADDER:-<none, fixed h2>}" \
     "h2_binary_mult=${H2_BINARY_MULT:-1.0}"
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
    ${PREVALENCE:+--prevalence $PREVALENCE} ${NCP_LADDER:+--ncp-ladder $NCP_LADDER} \
    ${H2_BINARY_MULT:+--h2-binary-mult $H2_BINARY_MULT} \
    ${H2:+--h2 $H2} ${H2_PC:+--h2-pc $H2_PC} ${SEED:+--seed $SEED} \
    2>&1 | tee "$D/02_build.log"

step "03 plink/plink2 sets, PCs"
bash $S/03_plink_prep.sh 2>&1 | tee "$D/03_plink.log"
$PY $S/04_add_pcs.py "$D" 4

step "03b pedigree-truth sparse GRM"
# The pedigree is known, so the relatedness matrix can be WRITTEN rather than
# estimated. Cheap, exact, and produced for both configs -- chr20 uses it as a
# reference to score step 0 against, bench50k consumes it directly.
$PY $S/03b_pedigree_grm.py --dir "$D" 2>&1 | tee "$D/03b_pedigree_grm.log"

# Which GRM steps 1-2 actually consume.
#
# bench50k uses the pedigree truth, because a GRM estimated from one 25 Mb
# window is expected to be unusable for a genetic rather than a counting reason:
# 25 Mb of chr20 is ~41 cM on the HapMap map (108.1 cM over 64.44 Mb, at 0.98x
# the chromosome mean rate), so the number of effectively INDEPENDENT common
# markers is a few hundred however many markers the window contains. That is an
# a priori argument -- no GRM has been built on a 25 Mb window here. What is
# measured is the full-chromosome case at 2.6x the genetic length, which already
# keeps 1,951 spurious pairs against 2,480 true ones. Population structure is
# ruled out (02_build_dataset.py keeps EUR only).
#
# Either way step 0 is now a CONTROL, not a dependency: it is scored against the
# truth rather than trusted. RUN_STEP0 defaults on for chr20, which consumes its
# output, and off for bench50k, which does not.
if [ "$CONFIG" = "bench50k" ]; then
  export GRMFILE="$D/pedigree_truth.sparseGRM.mtx"
  export GRMIDS="$GRMFILE.sampleIDs.txt"
  echo "steps 1-2 will consume the PEDIGREE-TRUTH GRM: $GRMFILE"
  : "${RUN_STEP0:=0}"
else
  : "${RUN_STEP0:=1}"
fi

if [ "$RUN_STEP0" != "1" ]; then
  echo
  echo "skipping step 0 (RUN_STEP0=$RUN_STEP0). Set RUN_STEP0=1 to run it as a"
  echo "control and score it against $D/pedigree_truth.sparseGRM.mtx."
else

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

fi   # RUN_STEP0

step "05 steps 1 and 2, all traits"
bash $S/07_run_pipeline.sh

step "06 results vs simulated truth"
$RS $S/06_analyse.R 2>&1 | tee "$D/RESULTS.txt"
echo; echo "### done -- outputs in $D"
