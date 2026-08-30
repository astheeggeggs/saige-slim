#!/usr/bin/env bash
# Compare a candidate run against the goldens.
#   bash harness/compare_goldens.sh <golden-dir> <candidate-dir>
#
# step0 + step2: byte-for-byte (step 2 is deterministic at defaults; p-values are
#   written pre-rounded %.6E). *.index files and logs are excluded.
# step1: toleranced via compare_step1.R (theta/VR at 1e-2 rel, rest at 1e-3).

set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
GOLD="$(cd "$1" && pwd)"; CAND="$(cd "$2" && pwd)"
fails=0

echo "### byte comparison: step0 + step2"
for sub in step0 step2; do
  while IFS= read -r -d '' gf; do
    rel="${gf#"$GOLD/"}"
    cf="$CAND/$rel"
    if [ ! -f "$cf" ]; then echo "MISSING in candidate: $rel"; fails=$((fails+1)); continue; fi
    if cmp -s "$gf" "$cf"; then echo "identical: $rel"
    else echo "DIFFERS:   $rel"; fails=$((fails+1)); fi
  done < <(find "$GOLD/$sub" -type f ! -name '*.index' ! -name '*.log' -print0 | sort -z)
  # files the candidate produced that the golden did not
  while IFS= read -r -d '' cf; do
    rel="${cf#"$CAND/"}"
    [ -f "$GOLD/$rel" ] || { echo "EXTRA in candidate: $rel"; fails=$((fails+1)); }
  done < <(find "$CAND/$sub" -type f ! -name '*.index' ! -name '*.log' -print0 | sort -z)
done

echo "### toleranced comparison: step1"
for grda in "$GOLD"/step1/*.rda; do
  base="$(basename "$grda")"
  prefix="${base%_noLOCO.rda}"; prefix="${prefix%.rda}"
  pixi run --manifest-path "$REPO/pixi.toml" Rscript "$REPO/harness/compare_step1.R" \
    "$GOLD/step1/$prefix" "$CAND/step1/$prefix" || fails=$((fails+1))
done

echo
if [ "$fails" -gt 0 ]; then echo "== HARNESS: $fails FAILURE(S)"; exit 1; fi
echo "== HARNESS: all comparisons PASS"
