#!/usr/bin/env bash
# Is the serial step-2 writer what caps --nThreads on a multi-model scan?
#
#   bash harness/bench_writer.sh <label> [reps]
#
# Measures, with NO source change, the share of mainMarkerInCPP's main-thread
# time that writeOutfile_single takes (w) at P = 1 and P = 5 models and
# T = 1/4/8 threads on the N = 50,000 fixture, plus paired wall/user timings.
# The question it answers is the docs/state/next.md item "format step-2 output
# rows on the workers": the ceiling that item can buy is 1/(1-w), so if w is
# small the item is closed by this number and not implemented.
#
# Two instruments, one run each per profiled cell and REPS timed runs per cell:
#   profile  harness/profile.sh (macOS `sample`, 1 ms, one call tree per
#            thread). harness/bench/writer_profile.py reads the tree.
#   time     /usr/bin/time -l, no profiler, for wall/user/sys/maxRSS.
# The T arms of one (P, trait) run adjacently within a rep and the order is
# reversed on even reps (the ordering control from harness/bench/NOTES.md).
#
# Controls the script enforces rather than trusts:
#   * refuses a throttled or loaded machine (as sim/upstream_compare/bench.sh);
#   * a threaded arm must PRINT its thread count and the P = 5 arm must PRINT
#     its model count (the LEDGER #30 shape: a knob that fails to arrive would
#     otherwise time two serial runs);
#   * every T = 4/8 output is cmp'd to the T = 1 output of the same cell, and
#     the P = 5 per-model outputs to the P = 1 outputs of the same trait, so a
#     number from a run that stopped producing the right bytes is never
#     recorded as a speed;
#   * every row carries the commit, whether the tree was dirty, and the md5 of
#     the SAIGE.so actually loaded, because the installed library is rebuilt by
#     whoever last ran build.sh and is not guaranteed to be HEAD.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="${1:?usage: bench_writer.sh <label> [reps]}"
REPS="${2:-2}"
SIM="${SIM:-/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k}"
THREADS="${THREADS:-1 4 8}"
PROFILE_AT="${PROFILE_AT:-1 8}"      # thread counts that also get a sampling profile
CELLS="${CELLS:-p1:y_binary p1:y_quant p5:mix p3:quant}"
RS="$REPO/.pixi/envs/default/lib/R/bin/Rscript"
SCRIPT="$REPO/extdata/step2_SPAtests.R"
WORK="$REPO/harness/work/writer"
OUT="$REPO/harness/bench/writer.tsv"
PROFDIR="$REPO/harness/bench/profiles"
mkdir -p "$WORK/out" "$WORK/logs" "$PROFDIR"
source "$REPO/harness/env.sh"          # pins the BLAS/TBB env; --nThreads overrides it explicitly

lim=$(pmset -g therm 2>/dev/null | awk '/CPU_Speed_Limit/{print $NF}')
load=$(uptime | sed 's/.*load averages*: *//' | awk '{print $1}' | tr -d ,)
echo "=== bench_writer: $LABEL  CPU_Speed_Limit=${lim:-?} load1=$load  $(date '+%F %T')"
if [ "${FORCE:-0}" != 1 ] && [ "${DRY:-0}" != 1 ]; then
  [ "${lim:-100}" = 100 ] || { echo "REFUSE: CPU throttled ($lim). FORCE=1 to override." >&2; exit 1; }
  awk -v l="$load" 'BEGIN{exit !(l > 2.5)}' && { echo "REFUSE: load $load > 2.5. FORCE=1 to override." >&2; exit 1; }
fi
for f in "$SIM/chr20_sim.bed" "$SIM/step1_y_binary.rda" "$SIM/pedigree_truth.sparseGRM.mtx"; do
  [ -f "$f" ] || { echo "missing fixture file $f" >&2; exit 1; }
done

COMMIT=$(git -C "$REPO" rev-parse --short HEAD)
DIRTY=$([ -z "$(git -C "$REPO" status --porcelain -- src R extdata)" ] && echo clean || echo dirty)
LIB=$("$RS" -e 'cat(find.package("SAIGE"))' 2>/dev/null)
SOMD5=$(md5 -q "$LIB/libs/SAIGE.so" 2>/dev/null)
echo "commit=$COMMIT ($DIRTY src/R/extdata)  SAIGE.so=$SOMD5"

# The P = 5 manifest: the five committed fixture models, absolute paths, the
# 3-column headerless shape step 1 writes (docs/state/multimodel-step2.md).
# p3:quant is the worst case for the writer: the per-model kernel is smallest
# there (no SPA/ER/Firth), so the P rows per marker weigh most against it.
MANIFEST="$WORK/manifest5.txt"; MANIFEST3="$WORK/manifest3q.txt"
: > "$MANIFEST"; : > "$MANIFEST3"
for t in y_quant y_binary y_quant_pc y_binary_pc y_null; do
  printf '%s\t%s\t%s\n' "$t" "$SIM/step1_$t.rda" "$SIM/step1_$t.varianceRatio.txt" >> "$MANIFEST"
done
for t in y_quant y_quant_pc y_null; do
  printf '%s\t%s\t%s\n' "$t" "$SIM/step1_$t.rda" "$SIM/step1_$t.varianceRatio.txt" >> "$MANIFEST3"
done

# Flags identical to step2_single in harness/sim/04_run_saige.sh, so these
# rows are comparable to sim/upstream_compare/RESULTS.md.
COMMON=(--bedFile="$SIM/chr20_sim.bed" --bimFile="$SIM/chr20_sim.bim" --famFile="$SIM/chr20_sim.fam"
        --chrom=20 --minMAF=0 --minMAC=20
        --sparseGRMFile="$SIM/pedigree_truth.sparseGRM.mtx"
        --sparseGRMSampleIDFile="$SIM/pedigree_truth.sparseGRM.mtx.sampleIDs.txt"
        --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_fastTest=TRUE)

model_args() { # model_args <cell> -> the model flags for that cell
  case "$1" in
    p1:*) local t="${1#p1:}"; echo "--GMMATmodelFile=$SIM/step1_$t.rda --varianceRatioFile=$SIM/step1_$t.varianceRatio.txt" ;;
    p5:*) echo "--GMMATmodel_varianceRatio_multiTraits_File=$MANIFEST" ;;
    p3:*) echo "--GMMATmodel_varianceRatio_multiTraits_File=$MANIFEST3" ;;
    *) echo "unknown cell $1" >&2; exit 2 ;;
  esac
}
outputs_of() { # outputs_of <cell> <prefix> -> the result files that run writes
  case "$1" in
    p1:*) echo "$2" ;;
    p5:*) local t; for t in y_quant y_binary y_quant_pc y_binary_pc y_null; do echo "${2}_$t"; done ;;
    p3:*) local t; for t in y_quant y_quant_pc y_null; do echo "${2}_$t"; done ;;
  esac
}

[ -f "$OUT" ] || printf 'label\tcommit\tdirty\tso_md5\tcell\tthreads\tmode\trep\treal_s\tuser_s\tsys_s\tmaxRSS_MB\trows\tlog\tprofile\n' > "$OUT"

check_log() { # check_log <log> <threads> <cell>: the run must SAY what it did
  local lf="$1" nt="$2" cell="$3"
  grep -q 'Analysis done' "$lf" || { echo "RUN DID NOT FINISH: $lf" >&2; tail -20 "$lf" >&2; return 1; }
  if [ "$nt" -gt 1 ]; then
    grep -q "will use $nt in-process threads" "$lf" || { echo "ARM DID NOT THREAD ($nt): $lf" >&2; return 1; }
  else
    grep -q 'in-process threads' "$lf" && { echo "SERIAL ARM PARALLELISED: $lf" >&2; return 1; }
  fi
  case "$cell" in
    p5:*) grep -q 'will test 5 models' "$lf" || { echo "MANIFEST DID NOT ARRIVE: $lf" >&2; return 1; } ;;
    p3:*) grep -q 'will test 3 models' "$lf" || { echo "MANIFEST DID NOT ARRIVE: $lf" >&2; return 1; } ;;
    *)    grep -q 'will test .* models' "$lf" && { echo "P=1 ARM RAN MULTI-MODEL: $lf" >&2; return 1; } ;;
  esac
  return 0
}

run_cell() { # run_cell <cell> <threads> <mode profile|time> <rep>
  local cell="$1" nt="$2" mode="$3" rep="$4"
  local tag="${cell//:/_}.t$nt.$mode.r$rep"
  local prefix="$WORK/out/$tag" lf="$WORK/logs/$tag.log" tf="$WORK/logs/$tag.time" prof=""
  local margs; margs=$(model_args "$cell")
  rm -f "$prefix" "$prefix".index "$prefix"_* 2>/dev/null
  if [ "${DRY:-0}" = 1 ]; then
    echo "DRY $tag: $RS $SCRIPT ${COMMON[*]} $margs --nThreads=$nt --SAIGEOutputFile=$prefix"; return 0
  fi
  cat "$SIM/chr20_sim.bed" > /dev/null   # warm the page cache, as run_arm.sh does
  if [ "$mode" = profile ]; then
    # profile.sh cd's into extdata and writes profiles/<label>.txt + .run.log
    prof="$PROFDIR/writer_$LABEL.$tag.txt"
    bash "$REPO/harness/profile.sh" "writer_$LABEL.$tag" "$SCRIPT" "${COMMON[@]}" $margs \
      --nThreads="$nt" --SAIGEOutputFile="$prefix" > "$tf" 2>&1 || { echo "PROFILE RUN FAILED: $tag" >&2; tail -5 "$tf" >&2; exit 1; }
    cp "$PROFDIR/writer_$LABEL.$tag.run.log" "$lf"
    : > "$tf"   # no /usr/bin/time under the profiler; the timed reps carry the seconds
  else
    /usr/bin/time -l "$RS" "$SCRIPT" "${COMMON[@]}" $margs --nThreads="$nt" \
      --SAIGEOutputFile="$prefix" > "$lf" 2> "$tf" || { echo "RUN FAILED: $tag" >&2; tail -20 "$lf" >&2; tail -5 "$tf" >&2; exit 1; }
  fi
  check_log "$lf" "$nt" "$cell" || exit 1
  local real user sys rss rows
  real=$(awk '/real/{print $1}' "$tf" | head -1); user=$(awk '/real/{print $3}' "$tf" | head -1)
  sys=$(awk '/real/{print $5}' "$tf" | head -1)
  rss=$(awk '/maximum resident set size/{printf "%.1f", $1/1048576}' "$tf")
  rows=$(awk '/ markers were tested/{s+=$1} END{print s+0}' "$lf")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$LABEL" "$COMMIT" "$DIRTY" "$SOMD5" "$cell" "$nt" "$mode" "$rep" "${real:-NA}" "${user:-NA}" "${sys:-NA}" \
    "${rss:-NA}" "$rows" "$lf" "${prof:-NA}" >> "$OUT"
  echo "  $tag: real=${real:-NA}s user=${user:-NA}s rss=${rss:-NA}MB rows=$rows" >&2
}

# Byte-identity: T > 1 against T = 1 of the same cell (rep 1, timed), and the
# P = 5 per-model files against the P = 1 files of the same trait.
cmp_outputs() {
  local cell="$1" nt="$2" mode="$3" rep="$4"
  [ "${DRY:-0}" = 1 ] && return 0
  local a="$WORK/out/${cell//:/_}.t$nt.$mode.r$rep" b="$WORK/out/${cell//:/_}.t1.time.r1" f
  if [ "$nt" != 1 ] || [ "$mode" != time ] || [ "$rep" != 1 ]; then
    for f in $(outputs_of "$cell" "$a"); do
      cmp -s "$f" "${f/$a/$b}" || { echo "ARMS DISAGREE (bytes): $f vs ${f/$a/$b}" >&2; exit 1; }
    done
  fi
  case "$cell" in p5:*|p3:*)
    for t in y_binary y_quant; do
      [ "$cell" = p3:quant ] && [ "$t" = y_binary ] && continue
      [ -f "$WORK/out/p1_$t.t1.time.r1" ] && { cmp -s "${a}_$t" "$WORK/out/p1_$t.t1.time.r1" \
        || { echo "P=5 DISAGREES WITH P=1 (bytes): ${a}_$t" >&2; exit 1; }; }
    done ;;
  esac
}

for cell in $CELLS; do
  echo "--- cell $cell"
  for rep in $(seq 1 "$REPS"); do
    order="$THREADS"; [ $((rep % 2)) = 0 ] && order=$(echo "$THREADS" | tr ' ' '\n' | tail -r | tr '\n' ' ')
    for nt in $order; do
      run_cell "$cell" "$nt" time "$rep"; cmp_outputs "$cell" "$nt" time "$rep"
    done
  done
  for nt in $PROFILE_AT; do
    run_cell "$cell" "$nt" profile 1; cmp_outputs "$cell" "$nt" profile 1
  done
done
echo "=== done $(date '+%F %T'); rows in $OUT; summarise with: python3 harness/bench/writer_profile.py harness/bench/profiles/writer_$LABEL.*.txt"
