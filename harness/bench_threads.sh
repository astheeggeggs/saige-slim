#!/usr/bin/env bash
# S6 measurement: the in-process marker threading (nThreadsMarker), 1 vs N.
#
#   bash harness/bench_threads.sh <label> [reps] [threads] [mech]
#
# mech is `inproc` (default) or `fork`, selecting which of the two parallelism
# mechanisms the N-thread arm uses -- nThreadsMarker (S6, in-process TBB) or
# --nThreads (C2, mclapply shards). The 1-thread arm is the same serial run
# either way, which is what makes the two mechanisms comparable through it:
# each is measured against the SAME baseline rather than against each other in
# a single interleaved block. That matters for the retirement decision in
# docs/design/step2-threading-verdict.md, item 3.
#
# This is NOT harness/bench.sh's job. bench.sh times ONE configuration per row
# and appends; this one is a PAIRED A/B in the sense harness/bench/NOTES.md
# established as the only method that resolved anything here: the two arms run
# adjacent within each rep, and the arm ORDER is swapped between the two blocks
# so a cache or drift effect cannot masquerade as the treatment. The arms are
# the same installed library with a knob moved, so unlike the F8/S1 comparisons
# no second library is needed -- which removes one confounder rather than
# adding one.
#
# Three workloads, because the ER finding makes ONE workload unable to
# represent the change (docs/design/step2-threading-verdict.md, the amendment
# to the retirement criterion):
#   binary        s2_marker_all config. Binary + default max_MAC_for_ER=4, so
#                 every MAC<=4 marker is deferred and re-run on the main thread.
#   binary_sparse s2_marker_all_sparse config. Same deferral, but with the PCG
#                 and the marker-path scoreTest on the clock, so the threadable
#                 fraction of the run is much larger.
#   quant         the quantitative model over the same markers. Never enters
#                 the ER branch at all, so it is the upper bound: what the
#                 mechanism does when nothing is serialised.
# Reporting only the first would understate the mechanism; only the last would
# overstate it. Both are needed, which is the point.
#
# Byte-identity between the arms is checked on every pair -- a speedup measured
# on a run that stopped producing the right answer is not a speedup.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="${1:?usage: bench_threads.sh <label> [reps] [threads]}"
REPS="${2:-8}"
NT="${3:-4}"
MECH="${4:-inproc}"
case "$MECH" in inproc|fork) ;; *) echo "mech must be inproc or fork" >&2; exit 2 ;; esac
cd "$REPO/extdata"
source "$REPO/harness/env.sh"
mkdir -p "$REPO/harness/bench"
OUT="$REPO/harness/bench/threads.tsv"
[ -f "$OUT" ] || printf 'label\tcommit\tworkload\tmech\tthreads\tblock\trep\telapsed_s\tuser_s\tsys_s\tmaxRSS_MB\tstartup_s\twork_s\tthroughput\n' > "$OUT"
COMMIT=$(git -C "$REPO" rev-parse --short HEAD)
SCRATCH=$(mktemp -d)
trap 'rm -rf "$SCRATCH"' EXIT
GOLD="$REPO/harness/goldens"
SG="--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
BED="--bedFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bed --bimFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.bim --famFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr.fam --AlleleOrder=alt-first"
COMMON="--minMAF=0 --minMAC=1 --LOCO=FALSE --markers_per_chunk=10000"

# Copied verbatim from bench.sh: the ONLY way to see a marker-kernel change on
# this fixture, since a whole-genome run at N=1000 is ~2.7 s of fixed prologue.
# Threading moves work_s and cannot move startup_s, so work_s is the arm the
# statistics are run on and elapsed_s is reported beside it for honesty.
s2_breakdown() {
  awk '
    /^[[:space:]]*[0-9.]+[[:space:]]+[0-9.]+[[:space:]]+[0-9.]+[[:space:]]*$/ {
      cum = $3; if (n++ == 0) first = cum; last = cum
    }
    / markers were tested/ { markers += $1 }
    END {
      if (n < 2) { print "NA\tNA\tNA"; exit }
      wpc = (last - first) / (n - 1)
      startup = first - wpc
      work = last - startup
      printf "%.3f\t%.3f\t%s\n", startup, work,
             (work > 0 ? sprintf("%.0f", markers / work) : "NA")
    }' "$1"
}

run_arm() { # run_arm <workload> <threads> <block> <rep> -> writes a TSV row, echoes the output path
  local wl="$1" nt="$2" blk="$3" rep="$4"
  local tag="$wl.$MECH.t$nt.b$blk.r$rep"
  local tf="$SCRATCH/time.$tag" lf="$SCRATCH/log.$tag" of="$SCRATCH/out.$tag"
  local model args
  case "$wl" in
    binary)        model="binary_sparseGRM_a";  args="" ;;
    binary_sparse) model="binary_sparseGRM_vr"; args="$SG --is_fastTest=FALSE" ;;
    quant)         model="quant_sparseGRM";     args="" ;;
    *) echo "unknown workload $wl" >&2; exit 1 ;;
  esac
  # The 1-thread arm is identical under both mechanisms (both dispatch to the
  # same serial call), so only the N-thread arm differs.
  local mt=1 ft=1
  if [ "$nt" -gt 1 ]; then
    if [ "$MECH" = inproc ]; then mt=$nt; else ft=$nt; fi
  fi
  SAIGE_MARKER_THREADS="$mt" /usr/bin/time -l \
    pixi run --manifest-path "$REPO/pixi.toml" Rscript step2_SPAtests.R \
    $BED $COMMON $args --nThreads=$ft \
    --GMMATmodelFile="$GOLD/step1/$model.rda" \
    --varianceRatioFile="$GOLD/step1/$model.varianceRatio.txt" \
    --SAIGEOutputFile="$of" > "$lf" 2> "$tf" || {
      echo "BENCH RUN FAILED: $tag"; tail -20 "$lf"; tail -5 "$tf"; exit 1; }
  # A threaded arm must SAY it threaded. Without this the whole comparison can
  # be two serial runs -- the LEDGER #30 shape, and exactly what the standing
  # test in test-step2-threads-inprocess.R asserts for the same reason.
  if [ "$nt" -gt 1 ]; then
    if [ "$MECH" = inproc ]; then
      grep -q "will use $nt in-process threads" "$lf" || {
        echo "ARM DID NOT THREAD: $tag (no thread-count line)"; exit 1; }
    else
      grep -q "shards over $nt threads" "$lf" || {
        echo "ARM DID NOT FORK: $tag (no shard line)"; exit 1; }
    fi
  else
    grep -qE "in-process threads|shards over" "$lf" && {
      echo "SERIAL ARM PARALLELISED: $tag"; exit 1; } || true
  fi
  local real user sys rss startup work thru
  real=$(awk '/real/{print $1}' "$tf" | head -1)
  user=$(awk '/real/{print $3}' "$tf" | head -1)
  sys=$(awk '/real/{print $5}' "$tf" | head -1)
  rss=$(awk '/maximum resident set size/{printf "%.1f", $1/1048576}' "$tf")
  IFS=$'\t' read -r startup work thru < <(s2_breakdown "$lf")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$LABEL" "$COMMIT" "$wl" "$MECH" "$nt" "$blk" "$rep" "$real" "$user" "$sys" "$rss" \
    "$startup" "$work" "$thru" >> "$OUT"
  echo "  $wl/$MECH t=$nt b$blk r$rep: ${real}s (${startup}s startup + ${work}s work), ${thru} mk/s, ${rss}MB" >&2
  echo "$of"
}

echo "=== bench_threads: $LABEL @ $COMMIT, mech=$MECH, $REPS reps x 2 blocks, 1 vs $NT threads"
pmset -g therm 2>/dev/null | grep -E 'CPU_Speed_Limit' || true
uptime

for wl in binary binary_sparse quant; do
  echo "--- workload $wl"
  for blk in 1 2; do
    for rep in $(seq 1 "$REPS"); do
      # Block 1 runs the serial arm first, block 2 the threaded arm first. That
      # swap IS the ordering control; without it a warm-cache effect on the
      # second arm of every pair is indistinguishable from the treatment.
      if [ "$blk" = 1 ]; then
        a=$(run_arm "$wl" 1 "$blk" "$rep"); b=$(run_arm "$wl" "$NT" "$blk" "$rep")
      else
        b=$(run_arm "$wl" "$NT" "$blk" "$rep"); a=$(run_arm "$wl" 1 "$blk" "$rep")
      fi
      cmp -s "$a" "$b" || { echo "ARMS DISAGREE (bytes): $wl block $blk rep $rep"; exit 1; }
    done
  done
done

echo "=== done; rows in $OUT"
