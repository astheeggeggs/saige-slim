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

run_tier() {
  local mode=$1 log="harness/work/prepush-$1.log"
  echo "== standard tier, SAIGE_COMPARE_MODE=$mode (log: $log)"
  if ! SAIGE_COMPARE_MODE="$mode" bash harness/run_tests.sh standard > "$log" 2>&1; then
    echo "prepush: FAILED in $mode mode; last failures from $log:" >&2
    grep -nE 'FAIL|Error|Failure' "$log" | tail -20 >&2
    exit 1
  fi
  grep -E '^\[ FAIL' "$log" | tail -1
}
run_tier byte
run_tier numeric
echo "prepush gate OK at $(git rev-parse --short HEAD)"
