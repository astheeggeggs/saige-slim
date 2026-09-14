#!/bin/bash
# Pooled-row calibration over the replicate nulls of make_null_series.py:
# per replicate, the unweighted burden (calibrated reference), the standard
# 6-cell grid, FlexRV 144 sets (informative score; under the null the arm is
# immaterial) and the MAF-only ladder, collapse 10. Then the ECDF at alpha
# 0.05 and 0.01 per replicate and pooled, with binomial 95% intervals.
# usage: null_series.sh [R] [P]   (env: OUT, D as run.sh)
set -uo pipefail
S=$(cd "$(dirname "$0")" && pwd)
R=${1:-5}; P=${2:-2}
export OUT=${OUT:-/Users/dpalmer/Repositories/SAIGE/harness/work/flexrv_study/results-nullseries}
export MINLOAD_REFUSE=${MINLOAD_REFUSE:-999}
D=${D:-/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k}
for r in $(seq 1 "$R"); do
  t=y_null_r$r
  [ -f "$D/step1_$t.rda" ] || { echo "no model $D/step1_$t.rda" >&2; exit 2; }
  for job in "plain unw" "plain grid" "informative flex192" "informative flexmaf"; do
    set -- $job
    tag=${t}_$1_$2_c10_r1
    if [ -s "$OUT/$tag.txt" ] && grep -q 'exit=0' "$OUT/$tag.log" 2>/dev/null; then continue; fi
    echo "$t $1 $2 10"
  done
done | xargs -P "$P" -L 1 bash -c '"$0/run.sh" "$@" || echo "FAILED $*"' "$S"
RS=/Users/dpalmer/Repositories/SAIGE/.pixi/envs/default/lib/R/bin/Rscript
$RS - "$OUT" "$R" <<'EOF'
args <- commandArgs(TRUE); OUT <- args[1]; R <- as.integer(args[2])
rd <- function(p) utils::read.table(p, header = TRUE, sep = "\t", colClasses = "character",
                                    quote = "", comment.char = "", na.strings = character(0))
gene_p <- function(f) {
  d <- rd(f)
  if (any(d$Group == "Cauchy")) { z <- d[d$Group == "Cauchy", ]
    as.numeric(if ("Pvalue" %in% colnames(z)) z$Pvalue else z$Pvalue_Burden)
  } else { z <- d[d$Group != "Cauchy", ]
    as.numeric(if ("Pvalue" %in% colnames(z)) z$Pvalue else z$Pvalue_Burden) }
}
cat("Pooled-row calibration over", R, "replicate nulls (210 genes each): ECDF at alpha, binomial 95% CI\n")
for (m in c("plain_unw", "plain_grid", "informative_flex192", "informative_flexmaf")) {
  all <- c()
  for (r in seq_len(R)) {
    f <- file.path(OUT, sprintf("y_null_r%d_%s_c10_r1.txt", r, m))
    if (!file.exists(f)) next
    p <- gene_p(f); all <- c(all, p)
    cat(sprintf("  %-20s r%d  alpha=0.05 obs=%.4f   alpha=0.01 obs=%.4f   n=%d\n", m, r,
                mean(p < 0.05, na.rm = TRUE), mean(p < 0.01, na.rm = TRUE), sum(!is.na(p))))
  }
  for (al in c(0.05, 0.01)) {
    k <- sum(all < al, na.rm = TRUE); n <- sum(!is.na(all)); ci <- binom.test(k, n)$conf.int
    flag <- if (al < ci[1] || al > ci[2]) " <-- OUTSIDE" else ""
    cat(sprintf("  %-20s POOLED alpha=%.2f obs=%.4f [%.4f, %.4f] n=%d%s\n", m, al, k / n, ci[1], ci[2], n, flag))
  }
}
EOF
