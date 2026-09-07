#!/usr/bin/env bash
# The fixed runs: SKAT-O with each candidate fix installed (SAIGE_FIX), on the
# binary model at both category boundaries, the quantitative goldens (no-op
# check), and the sparse-ratio-2.0 file where the clamp already held (no-op).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="${DIAG_OUT:-$REPO/harness/work/diag-burden}"
export DIAG_OUT="$OUT"
run() { local fix=$1; shift; SAIGE_FIX="$fix" "$HERE/run.sh" "$@"; }
runq() { local fix=$1; shift; SAIGE_FIX="$fix" "$HERE/run_quant.sh" "$@"; }
{
  run 103  R103 0 "$OUT/vr/committed.txt" 10,20.5 20.5 &
  run 103A RA0  0 "$OUT/vr/committed.txt" 10,20.5 20.5 &
  wait
  run 103B RB0  0 "$OUT/vr/committed.txt" 10,20.5 20.5 &
  run 103A RA4  0 "$OUT/vr/committed.txt" 10,200 200 &
  wait
  run 103B RB4  0 "$OUT/vr/committed.txt" 10,200 200 &
  run 103A R2A  0 "$OUT/vr/sp2.txt" 10,20.5 20.5 &
  wait
  runq A RQA 0 &
  runq B RQB 0 &
  wait
}
G="$REPO/harness/goldens/step2"
cmp "$OUT/RQA/out.txt" "$G/plink_group_nofast_quant.txt" && echo "RQA (fix A, quantitative) == golden P37"
cmp "$OUT/RQB/out.txt" "$G/plink_group_nofast_quant.txt" && echo "RQB (fix B, quantitative) == golden P37"
