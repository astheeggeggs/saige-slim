#!/bin/bash
# Analyse saddlepoint-failure dumps INSIDE the SAIGE image (dsub): for each
# spaFail_<k>.tsv in the tarball, run the exported kernels on the vector the
# saddlepoint saw and print SUMMARY numbers only (no per-sample data leaves).
#   IN_DUMP  the .spaDump.tgz delivered by a run with SPA_DUMP=1
set -euo pipefail
: "${IN_DUMP:?}"
mkdir -p dumps && tar xzf "${IN_DUMP}" -C dumps
n=$(find dumps -name 'spaFail_*.tsv' | wc -l); echo "dumps: ${n}"
cat > check.R <<'R'
suppressMessages(library(SAIGE))
f <- commandArgs(TRUE)[1]
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]
d <- read.table(f, header = TRUE, skip = 1, sep = "\t")
q <- as.numeric(hdr[2]); qinv <- as.numeric(hdr[4]); m1 <- as.numeric(hdr[6]); why <- hdr[length(hdr)]
mu <- d$mu; g <- d$gtilde; tol <- .Machine$double.eps^0.25
gpos <- sum(g[g > 0]); gneg <- sum(g[g < 0]); gmax <- max(abs(g)); k <- sum(g > 0)
cat(sprintf("== %s | N=%d carriers(g>0)=%d gmax=%.3g | q=%.4g qinv=%.4g m1=%.4g gneg=%.4g gpos=%.4g | mu(carriers) min=%.2e max=%.2e | %s\n",
            basename(f), length(g), k, gmax, q, qinv, m1, gneg, gpos, min(mu[g > 0]), max(mu[g > 0]), why))
r1 <- getroot_K1_Binom(0, mu, g, q, tol, 1000L); r2 <- getroot_K1_Binom(0, mu, g, qinv, tol, 1000L)
cat(sprintf("   root(q)=%.4g conv=%s niter=%d | root(qinv)=%.4g conv=%s niter=%d | mirror attainable=%s\n",
            r1$root, r1$Isconverge, r1$niter, r2$root, r2$Isconverge, r2$niter, qinv > gneg && qinv < gpos))
if (r1$Isconverge) { s1 <- Get_Saddle_Prob_Binom(r1$root, mu, g, q, FALSE)
  cat(sprintf("   q-side saddle: isSaddle=%s pval=%.4g | K1(root1)-q=%.3g K2(root1)=%.3g Korg(root1)=%.4g\n",
              s1$isSaddle, s1$pval, K1_adj_Binom(r1$root, mu, g, q), K2_Binom(r1$root, mu, g), Korg_Binom(r1$root, mu, g))) }
side <- if (qinv < q) -1 else 1
b <- sapply(0:24, function(kk) { tt <- side * 2^kk / gmax; exp(Korg_Binom(tt, mu, g) - tt * qinv) })
cat(sprintf("   chernoff bound on the mirror tail: min=%.3g at k=%d (t=%.3g); K1-qinv at t=-2^k/gmax: %s\n",
            min(b, na.rm = TRUE), which.min(b) - 1, side * 2^(which.min(b) - 1) / gmax,
            paste(sprintf("%.2g", sapply(c(0, 5, 10, 15, 20, 24), function(kk) K1_adj_Binom(side * 2^kk / gmax, mu, g, qinv))), collapse = " ")))
p <- SPA_pval(mu, g, q, qinv, 1e-100, tol, FALSE, "binary", FALSE)
cat(sprintf("   composed p (this image): %.4g\n", p))
R
for f in $(find dumps -name 'spaFail_*.tsv' | sort -V); do Rscript check.R "$f"; done
echo done
