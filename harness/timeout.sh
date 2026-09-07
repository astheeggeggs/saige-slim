#!/usr/bin/env bash
# A per-command timeout for the harness scripts. Sourced, not executed.
#
# macOS ships no timeout(1) and the pinned pixi env has no coreutils, so this
# is a bash watchdog. It exists because a sanitizer PLANT is licensed to
# misbehave arbitrarily: on 2026-09-07 a planted TSAN arm segfaulted, R's
# segfault handler printed its banner and then blocked on the interactive abort
# prompt, and harness/tsan.sh sat on it for 4h32m before it was killed by hand.
# A control that can hang forever is not a control.
#
#   source harness/timeout.sh
#   run_with_timeout <seconds> <logfile> <cmd> [args...]
#   # then read $RWT_STATUS: the command's exit status, or 124 if it timed out
#
# stdin is /dev/null so that R's abort prompt gets EOF rather than waiting, and
# any children the command forked (SkatoParallel's pool, for one) are killed
# with it -- killing only the parent leaves them alive and reparented.
run_with_timeout() {
  local secs="$1" log="$2"; shift 2
  "$@" > "$log" 2>&1 < /dev/null &
  local pid=$! waited=0
  while kill -0 "$pid" 2>/dev/null; do
    if [ "$waited" -ge "$secs" ]; then
      echo "run_with_timeout: no exit after ${secs}s; killing $pid" >> "$log"
      pkill -TERM -P "$pid" 2>/dev/null || true
      kill -TERM "$pid" 2>/dev/null || true
      sleep 2
      pkill -KILL -P "$pid" 2>/dev/null || true
      kill -KILL "$pid" 2>/dev/null || true
      wait "$pid" 2>/dev/null || true
      RWT_STATUS=124
      return 0
    fi
    sleep 1
    waited=$((waited + 1))
  done
  RWT_STATUS=0
  wait "$pid" || RWT_STATUS=$?   # `|| ...` so a non-zero exit is not fatal under set -e
  return 0
}
