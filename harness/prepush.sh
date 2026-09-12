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

run_tier() {
  local mode=$1 log="harness/work/prepush-$1.log"
  local t0 idle ld secs
  t0=$(date +%s)
  idle=$(cpu_idle)
  ld=$(uptime | sed 's/.*load averages*: //')
  echo "== standard tier, SAIGE_COMPARE_MODE=$mode (log: $log)"
  echo "   start $(date '+%H:%M:%S')  load [$ld]  idle ${idle:-?}%"
  {
    echo "prepush: mode=$mode start=$(date '+%F %T') load=[$ld] idle=${idle:-?}%"
  } > "$log"
  if ! SAIGE_COMPARE_MODE="$mode" bash harness/run_tests.sh standard >> "$log" 2>&1; then
    secs=$(elapsed "$t0")
    echo "prepush: FAILED in $mode mode after $(hms "$secs"); last failures from $log:" >&2
    grep -nE 'FAIL|Error|Failure' "$log" | tail -20 >&2
    exit 1
  fi
  secs=$(elapsed "$t0")
  echo "prepush: mode=$mode wall=${secs}s ($(hms "$secs"))" >> "$log"
  printf '%s   %s wall %s\n' "$(grep -E '^\[ FAIL' "$log" | tail -1)" "$mode" "$(hms "$secs")"
  eval "SECS_$mode=$secs"
}

GATE_T0=$(date +%s)
run_tier byte
run_tier numeric
GATE_SECS=$(elapsed "$GATE_T0")

# The tree, not the commit: a dirty tree is what a green run then licensed.
DIRTY=$(git status --porcelain -- src R tests extdata harness | head -5)
if [ -n "$DIRTY" ]; then
  echo "prepush gate OK on a DIRTY TREE at $(git rev-parse --short HEAD) + uncommitted changes:"
  echo "$DIRTY" | sed 's/^/           /'
  echo "           (the commit named above does NOT contain what was just tested)"
else
  echo "prepush gate OK at $(git rev-parse --short HEAD) (clean tree)"
fi
echo "prepush timing: byte $(hms "${SECS_byte:-0}") + numeric $(hms "${SECS_numeric:-0}") = $(hms "$GATE_SECS") total"
