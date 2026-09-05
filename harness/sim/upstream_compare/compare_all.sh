#!/usr/bin/env bash
# Run compare_pvalues.R + plot_pvalues.R over every (models, trait, kind) pair that has
# BOTH arms' outputs under harness/work/upstream_compare/{main,strip}/out.
#   compare_all.sh [rep] [nthreads]      (default rep 1, nthreads 1)
# Writes per-row TSVs to harness/work/upstream_compare/cmp/, a text summary to
# harness/sim/upstream_compare/cmp_summary.txt and PNGs to harness/sim/upstream_compare/plots/.
set -uo pipefail
rep=${1:-1}; nt=${2:-1}
REPO=/Users/dpalmer/Repositories/SAIGE
HERE=$REPO/harness/sim/upstream_compare
W=$REPO/harness/work/upstream_compare
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
mkdir -p "$W/cmp" "$HERE/plots"
SUM=$HERE/cmp_summary.txt; : > "$SUM"
traits="y_binary y_quant y_binary_pc y_quant_pc y_null"
for models in frozen own; do
  for kind in single gene gene_single; do
    files=()
    for t in $traits; do
      case $kind in
        single)      a=$W/main/out/s2single_${t}_r${rep}_nt${nt}_${models}.txt; b=$W/strip/out/s2single_${t}_r${rep}_nt${nt}_${models}.txt; k=single; pc=p.value ;;
        gene)        a=$W/main/out/s2group_${t}_r${rep}_nt${nt}_${models}.txt;  b=$W/strip/out/s2group_${t}_r${rep}_nt${nt}_${models}.txt;  k=gene;   pc=Pvalue ;;
        gene_single) a=$W/main/out/s2group_${t}_r${rep}_nt${nt}_${models}.txt.singleAssoc.txt; b=$W/strip/out/s2group_${t}_r${rep}_nt${nt}_${models}.txt.singleAssoc.txt; k=single; pc=p.value ;;
      esac
      [ -s "$a" ] && [ -s "$b" ] || continue
      out=$W/cmp/${kind}_${t}_${models}.tsv
      if [ "$models" = own ]; then vra=$W/main/models/step1_$t.varianceRatio.txt; vrb=$W/strip/models/step1_$t.varianceRatio.txt; else vra=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k/step1_$t.varianceRatio.txt; vrb=$vra; fi
      $RS "$HERE/compare_pvalues.R" "$k" "$a" "$b" "$out" "${kind} ${t} models=${models}" "$vra" "$vrb" | tee -a "$SUM"
      files+=("$out")
    done
    [ ${#files[@]} -gt 0 ] || continue
    title="strip HEAD vs upstream main 365deeff -- $kind tests, N=50,000, step-1 models: $models"
    $RS "$HERE/plot_pvalues.R" "$HERE/plots/${kind}_${models}.png" "$title" "$pc" "${files[@]}"
    if [ "$kind" = gene ]; then
      $RS "$HERE/plot_pvalues.R" "$HERE/plots/${kind}_${models}_burden.png" "$title (Pvalue_Burden)" Pvalue_Burden "${files[@]}"
      $RS "$HERE/plot_pvalues.R" "$HERE/plots/${kind}_${models}_skat.png" "$title (Pvalue_SKAT)" Pvalue_SKAT "${files[@]}"
    fi
  done
done
echo "summary: $SUM"
