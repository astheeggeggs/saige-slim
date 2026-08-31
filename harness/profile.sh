#!/usr/bin/env bash
# Sampling profile of one SAIGE workload.
#
#   bash harness/profile.sh <label> <script.R> [args...]
#
# Writes a call tree to harness/bench/profiles/<label>.txt.
#
# Why this exists: before 2026-08-31 the repo had NO profiling tooling of any
# kind (no Rprof, perf, gprof, -pg anywhere), and the harness/bench/prof*.txt
# files are NOT profiles -- they are copies of an ordinary step-2 result table.
# So every speedup decision was being made from reasoning rather than
# measurement. Pair this with `bash harness/bench.sh`, whose startup_s/work_s/
# throughput columns say WHETHER something moved; this says WHERE the time went.
#
# Symbolisation: the release build in src/Makevars keeps exported symbols, so
# top-level functions resolve, but static/inlined frames may not. For a sharper
# profile do a one-off instrumented build and REVERT it before any gate:
#   PKG_CXXFLAGS='-g1' bash harness/build.sh
# Never commit a Makevars change for this -- its release flags are load-bearing
# (Phase 4 batch 1) and a debug build would silently invalidate every bench row.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
LABEL="${1:?usage: profile.sh <label> <script.R> [args...]}"; shift
SCRIPT="${1:?usage: profile.sh <label> <script.R> [args...]}"; shift

source "$REPO/harness/env.sh"
OUTDIR="$REPO/harness/bench/profiles"
mkdir -p "$OUTDIR"
OUT="$OUTDIR/$LABEL.txt"
cd "$REPO/extdata"

case "$(uname -s)" in
  Darwin)
    command -v sample >/dev/null || { echo "profile.sh: 'sample' not found"; exit 1; }
    pixi run --manifest-path "$REPO/pixi.toml" Rscript "$SCRIPT" "$@" \
      > "$OUTDIR/$LABEL.run.log" 2>&1 &
    RPID=$!
    # `sample` needs the R process, not the pixi/shell wrapper that spawned it.
    TARGET="$RPID"
    for _ in $(seq 1 40); do
      KID=$(pgrep -P "$TARGET" 2>/dev/null | head -1) || true
      [ -n "${KID:-}" ] || break
      TARGET="$KID"
    done
    # 9999 s ceiling: sample exits on its own when the target does.
    sample "$TARGET" 9999 1 -f "$OUT" >/dev/null 2>&1 || true
    wait "$RPID" || { echo "profile.sh: workload FAILED, see $OUTDIR/$LABEL.run.log"; exit 1; }
    ;;
  Linux)
    command -v perf >/dev/null || { echo "profile.sh: 'perf' not found"; exit 1; }
    perf record -g --call-graph=dwarf -o "$OUTDIR/$LABEL.data" -- \
      pixi run --manifest-path "$REPO/pixi.toml" Rscript "$SCRIPT" "$@" \
      > "$OUTDIR/$LABEL.run.log" 2>&1
    perf report -i "$OUTDIR/$LABEL.data" --stdio --sort=overhead,symbol > "$OUT"
    ;;
  *) echo "profile.sh: unsupported platform $(uname -s)"; exit 1 ;;
esac

echo "profile written: $OUT"
