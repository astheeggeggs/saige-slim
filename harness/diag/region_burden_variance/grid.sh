#!/usr/bin/env bash
# The whole experiment: seven (SKAT-O, BURDEN) pairs on the binary model over
# four variance-ratio files and two MAC-category boundaries, the quantitative
# pair, the golden byte checks, and the analyses. ~4 min at -P 4.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
OUT="${DIAG_OUT:-$REPO/harness/work/diag-burden}"
export DIAG_OUT="$OUT" SAIGE_REPO="$REPO"
mkdir -p "$OUT/vr"
printf '1.0 sparse 1\n1.0 null 1\n1.0 sparse 2\n1.0 null 2\n' > "$OUT/vr/all1.txt"
printf '2.0 sparse 1\n1.0 null 1\n2.0 sparse 2\n1.0 null 2\n' > "$OUT/vr/sp2.txt"
printf '0.5 sparse 1\n1.0 null 1\n0.5 sparse 2\n1.0 null 2\n' > "$OUT/vr/sp05.txt"
printf '1.0 sparse 1\n0.5 null 1\n1.0 sparse 2\n0.5 null 2\n' > "$OUT/vr/null05.txt"
cp "$REPO/extdata/output/example_binary_sparseGRM.varianceRatio.txt" "$OUT/vr/committed.txt"
{
  echo "R0 0 $OUT/vr/committed.txt 10,20.5 20.5"
  echo "B0 1 $OUT/vr/committed.txt 10,20.5 20.5"
  echo "R1 0 $OUT/vr/all1.txt 10,20.5 20.5"
  echo "B1 1 $OUT/vr/all1.txt 10,20.5 20.5"
  echo "R2 0 $OUT/vr/sp2.txt 10,20.5 20.5"
  echo "B2 1 $OUT/vr/sp2.txt 10,20.5 20.5"
  echo "R3 0 $OUT/vr/sp05.txt 10,20.5 20.5"
  echo "B3 1 $OUT/vr/sp05.txt 10,20.5 20.5"
  echo "R4 0 $OUT/vr/committed.txt 10,200 200"
  echo "B4 1 $OUT/vr/committed.txt 10,200 200"
  echo "R5 0 $OUT/vr/sp2.txt 10,200 200"
  echo "B5 1 $OUT/vr/sp2.txt 10,200 200"
  echo "R6 0 $OUT/vr/null05.txt 10,20.5 20.5"
  echo "B6 1 $OUT/vr/null05.txt 10,20.5 20.5"
} | xargs -P 4 -L 1 "$HERE/run.sh"
"$HERE/run_quant.sh" RQ 0 &
"$HERE/run_quant.sh" BQ 1 &
wait
G="$REPO/harness/goldens/step2"
cmp "$OUT/R0/out.txt" "$G/plink_group_nofast.txt" && echo "R0 == golden P6"
cmp "$OUT/R0/out.txt.singleAssoc.txt" "$G/plink_group_nofast.txt.singleAssoc.txt" && echo "R0 singleAssoc == golden P6"
cmp "$OUT/B0/out.txt" "$G/plink_group_burden.txt" && echo "B0 == golden P5"
cmp "$OUT/RQ/out.txt" "$G/plink_group_nofast_quant.txt" && echo "RQ == golden P37"
cmp "$OUT/BQ/out.txt" "$G/plink_group_burden_quant.txt" && echo "BQ == golden P36"
cd "$REPO"
for k in 0 1 2 3 4 5 6; do
  pixi run Rscript "$HERE/analyze.R" "R$k" "B$k" > "$OUT/R$k/analysis.txt" 2>&1 || echo "analyze R$k failed"
done
pixi run Rscript "$HERE/summary.R"
pixi run Rscript "$HERE/numdiff.R"
pixi run Rscript "$HERE/pcg_check.R"
pixi run Rscript "$HERE/quant_check.R"
