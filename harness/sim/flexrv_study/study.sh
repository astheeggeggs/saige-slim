#!/bin/bash
# FlexRV simulation study, stages S2-S4: every (trait, arm, method, collapse)
# run through run.sh, P at a time. Methods that never read the score (unw,
# beta, skato) run once on the fixture's own group file (arm=plain); the
# synonymous arm runs flex192 on the low-only mask (run.sh sets CELL=low) and
# its unweighted twin. Idempotent: a run whose output exists is skipped, so an
# interrupted study resumes. usage: study.sh [P]   (env: OUT, D, G as run.sh)
set -uo pipefail
S=$(cd "$(dirname "$0")" && pwd)
P=${1:-3}
OUT=${OUT:-/Users/dpalmer/Repositories/SAIGE/harness/work/flexrv_study/results}
export OUT
TRAITS=${TRAITS:-"y_null y_null_bin y_quant y_binary"}
ARMS=${ARMS:-"informative realistic uninformative"}
FLEX=${FLEX:-"flex192 flexscore flexmaf flex1"}
PLAIN=${PLAIN:-"unw beta skato grid"}
COLL=${COLL:-"0 10"}
jobs=()
for t in $TRAITS; do
  for c in $COLL; do
    for a in $ARMS; do for m in $FLEX; do jobs+=("$t $a $m $c"); done; done
    for m in $PLAIN; do jobs+=("$t plain $m $c"); done
    jobs+=("$t synonymous flex192 $c" "$t synonymous unw $c")
  done
done
echo "${#jobs[@]} runs, $P at a time, results in $OUT"
printf '%s\n' "${jobs[@]}" | while read -r t a m c; do
  tag=${t}_${a}_${m}_c${c}_r1
  if [ -s "$OUT/$tag.txt" ] && grep -q 'exit=0' "$OUT/$tag.log" 2>/dev/null; then
    echo "skip $tag" >&2; continue
  fi
  echo "$t $a $m $c"
done | xargs -P "$P" -L 1 bash -c '"$0/run.sh" "$@" || echo "FAILED $*"' "$S"
echo "study done: $(ls $OUT/*.txt 2>/dev/null | wc -l) result files"
