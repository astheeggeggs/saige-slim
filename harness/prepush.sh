#!/usr/bin/env bash
# The pre-push gate: what CI runs, plus what it cannot.
#   bash harness/prepush.sh
#
# 1. `ledger.py lint`  -- every LEDGER.md row parses, no hand-written counts.
# 2. standard tier in BYTE compare mode     -- valid only on the golden platform.
# 3. standard tier in NUMERIC compare mode  -- what CI (ubuntu) runs.
# Both modes, because the trap fired twice: a cross-run comparison passed
# locally in byte mode and went red on CI in numeric mode, once for five hours
# while the session recorded the tier as green.
#
# It tests the INSTALLED package and does not build (harness/build.sh is ~4 min
# from clean objects). It REFUSES if a source file is newer than the installed
# SAIGE.so, because a green run against a stale install is exactly the silent
# pass this gate exists to prevent; rebuild, or set SAIGE_PREPUSH_ALLOW_STALE=1
# if you know the mtimes lie (a checkout that touched files without changing
# them). `.githooks/pre-push` runs this on every push to the `slim` remote once
# `git config core.hooksPath .githooks` is set.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"
mkdir -p harness/work

python3 harness/tools/ledger.py lint || exit 1

SO=".pixi/envs/default/lib/R/library/SAIGE/libs/SAIGE.so"
if [ ! -f "$SO" ]; then echo "prepush: SAIGE is not installed ($SO missing); run bash harness/build.sh" >&2; exit 1; fi
stale=$(find src R -type f \( -name '*.cpp' -o -name '*.hpp' -o -name '*.h' -o -name '*.R' \) -newer "$SO" | head -5)
if [ -n "$stale" ] && [ "${SAIGE_PREPUSH_ALLOW_STALE:-0}" != "1" ]; then
  echo "prepush: source newer than the installed package -- rebuild first (bash harness/build.sh)," >&2
  echo "         or SAIGE_PREPUSH_ALLOW_STALE=1 if the mtimes lie. Newer files:" >&2
  echo "$stale" | sed 's/^/           /' >&2
  exit 1
fi

# What the gate CANNOT be trusted to say about itself, and now does.
#
# 1. TIMING. Wall seconds per compare mode and for the whole gate, on stdout
#    and as the first line of each mode's log. Without it, "the gate is slow"
#    is a feeling: the standard tier was recorded at ~2.5 min in CLAUDE.md, at
#    4m20s on 2026-09-05, and nobody had a per-mode number. `uptime` and the
#    CPU idle at each mode's START go in too, because on this box the load
#    average swung between 2.5 and 23 inside one session and a timing without
#    its conditions is not a measurement.
# 2. WHAT IT ACTUALLY GATED. The old summary printed `git rev-parse HEAD`,
#    which names the last COMMIT, not the tree that was tested. On 2026-09-12
#    a run that gated the uncommitted TSAN annotation printed "OK at 8f4d5db7"
#    -- a commit that does not contain the change -- and that line was then
#    quoted as evidence. The tree is what is tested, so a dirty tree is now
#    said out loud, with the paths.
elapsed() { echo $(( $(date +%s) - $1 )); }
hms() { printf '%dm%02ds' $(( $1 / 60 )) $(( $1 % 60 )); }
# the idle percentage, by NAME rather than by position: `sed` on the tail of
# the line dragged in whatever trailed "idle" and printed "77.58 " with a space
cpu_idle() {
  top -l 1 -n 0 2>/dev/null | awk -F'[ ,]+' '/^CPU usage/ {
    for (i = 1; i <= NF; i++) if ($i == "idle") { gsub(/%/, "", $(i-1)); print $(i-1); exit }
  }'
}

# The two compare modes run CONCURRENTLY (2026-09-12). They are independent
# processes: their own testthat worker pool, their own log, read-only goldens,
# and no standard-tier test writes to a fixed path (audited -- every output goes
# through tempfile(); the only writers into extdata/output are the full-tier
# regeneration tests, which skip here).
#
# The worker budget is SPLIT, not doubled: each mode gets ceil(nproc/2), so the
# two together ask for the parallelism one mode used to ask for alone. The
# reason that is close to free is that a mode's wall time is bounded by its
# SLOWEST FILE (CLAUDE.md), not by its total CPU, so half the box sat idle for
# much of a serial mode; concurrency fills it with the other mode's long pole
# instead of queueing behind it. That is a hypothesis, and the measurement is
# in harness/bench/NOTES.md -- not an argument for taking it on faith.
#
# MEASURED AND OFF BY DEFAULT (2026-09-12). The reasoning above is sound and the
# implementation works; the speedup is ZERO. Paired alternating reps, warm box:
# concurrent 14m15s / 14m08s against serial 14m03s / 14m20s -- means identical
# to the second (851.5s each), and the two reps disagree on the sign. The
# predictor that was right is the saturation trace: a single mode already runs
# at 93.8% CPU (mean idle 6.2%, 92% of samples under 10% idle), so there was
# almost no headroom to recover and total CPU work is what binds.
#
# With wall time a tie, the tie-breakers all favour serial: TIME TO FIRST
# FAILURE is 7 minutes serial against 14 concurrent, there are half as many
# processes, and concurrency is exposed to any future test that writes to a
# shared path (it found one on its first run -- see 396e66b4).
#
# SAIGE_PREPUSH_CONCURRENT=1 opts in. It is kept rather than deleted because
# the null is a fact about THIS box (4 physical cores, a saturated tier) and a
# machine with more cores than the tier can use would see a real gain.
CONCURRENT="${SAIGE_PREPUSH_CONCURRENT:-0}"
NPROC=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)

run_tier() {   # everything into the mode's own log; nonzero return on failure
  local mode=$1 log="harness/work/prepush-$1.log"
  local t0 idle ld secs
  t0=$(date +%s)
  idle=$(cpu_idle)
  ld=$(uptime | sed 's/.*load averages*: //')
  {
    echo "prepush: mode=$mode start=$(date '+%F %T') load=[$ld] idle=${idle:-?}%"
    echo "prepush: mode=$mode workers=${TESTTHAT_CPUS:-default} concurrent=$CONCURRENT"
  } > "$log"
  if ! SAIGE_COMPARE_MODE="$mode" bash harness/run_tests.sh standard >> "$log" 2>&1; then
    secs=$(elapsed "$t0")
    echo "prepush: mode=$mode result=FAILED wall_seconds=$secs" >> "$log"
    return 1
  fi
  secs=$(elapsed "$t0")
  echo "prepush: mode=$mode result=OK wall_seconds=$secs" >> "$log"
  return 0
}

tier_secs() {
  sed -n "s/^prepush: mode=$1 result=.* wall_seconds=\([0-9]*\)/\1/p" \
    "harness/work/prepush-$1.log" | tail -1
}
tier_report() {
  local mode=$1 rc=$2 log="harness/work/prepush-$1.log"
  printf '   %-8s wall %-8s %s\n' "$mode" "$(hms "$(tier_secs "$mode")")" \
    "$(grep -E '^\[ FAIL' "$log" | tail -1)"
  if [ "$rc" != 0 ]; then
    echo "prepush: FAILED in $mode mode; last failures from $log:" >&2
    grep -nE 'FAIL|Error|Failure' "$log" | tail -20 >&2
  fi
}

GATE_T0=$(date +%s)
if [ "$CONCURRENT" = 1 ]; then
  HALF=$(( (NPROC + 1) / 2 ))
  echo "== standard tier, BOTH compare modes at once, $HALF testthat workers each"
  echo "   start $(date '+%H:%M:%S')  load [$(uptime | sed 's/.*load averages*: //')]  idle $(cpu_idle)%"
  ( export TESTTHAT_CPUS=$HALF; run_tier byte )    & pid_byte=$!
  ( export TESTTHAT_CPUS=$HALF; run_tier numeric ) & pid_numeric=$!
  wait "$pid_byte";    rc_byte=$?
  wait "$pid_numeric"; rc_numeric=$?
else
  echo "== standard tier, one mode then the other (SAIGE_PREPUSH_CONCURRENT=0)"
  run_tier byte;    rc_byte=$?
  run_tier numeric; rc_numeric=$?
fi
GATE_SECS=$(elapsed "$GATE_T0")

# Both modes are reported, and both RUN even when one fails. The byte/numeric
# split exists because a cross-run comparison passed in one mode and went red
# in the other -- so WHICH of the two failed is the answer the gate is for, and
# the old version exited on the first failure and never ran the second.
tier_report byte    "$rc_byte"
tier_report numeric "$rc_numeric"
if [ "$rc_byte" != 0 ] || [ "$rc_numeric" != 0 ]; then
  echo "prepush: FAILED (byte=$rc_byte numeric=$rc_numeric) after $(hms "$GATE_SECS")" >&2
  exit 1
fi

# The tree, not the commit: a dirty tree is what a green run then licensed.
DIRTY=$(git status --porcelain -- src R tests extdata harness | head -5)
if [ -n "$DIRTY" ]; then
  echo "prepush gate OK on a DIRTY TREE at $(git rev-parse --short HEAD) + uncommitted changes:"
  echo "$DIRTY" | sed 's/^/           /'
  echo "           (the commit named above does NOT contain what was just tested)"
else
  echo "prepush gate OK at $(git rev-parse --short HEAD) (clean tree)"
fi
if [ "$CONCURRENT" = 1 ]; then
  echo "prepush timing: byte $(hms "$(tier_secs byte)") | numeric $(hms "$(tier_secs numeric)") | gate $(hms "$GATE_SECS") wall (concurrent)"
else
  echo "prepush timing: byte $(hms "$(tier_secs byte)") + numeric $(hms "$(tier_secs numeric)") = $(hms "$GATE_SECS") total (serial)"
fi
