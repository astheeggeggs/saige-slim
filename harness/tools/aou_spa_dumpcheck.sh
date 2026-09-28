#!/bin/bash
# Analyse saddlepoint-failure dumps INSIDE the SAIGE image (dsub): for each
# spaFail_<k>.tsv in the tarball, run the exported kernels on the vector the
# saddlepoint saw and print SUMMARY numbers only (no per-sample data leaves).
#   IN_DUMP  the .spaDump.tgz delivered by a run with SPA_DUMP=1
set -euo pipefail
: "${IN_DUMP:?}"
mkdir -p dumps && tar xzf "${IN_DUMP}" -C dumps
n=$(find dumps -name 'spaFail_*.tsv' | wc -l); nb=$(find dumps -name 'spaBelow_*.tsv' | wc -l); echo "dumps: ${n} failed, ${nb} converged-below (SAIGE_SPA_DUMP_BELOW)"
cat > check.R <<'R'
suppressMessages(library(SAIGE))
f <- commandArgs(TRUE)[1]
l2 <- readLines(f, n = 2)
hdr <- strsplit(l2[1], "\t")[[1]]
# since 2026-09-28 every dump carries a second header line (kernel, roots, the
# Gaussian remainder, the reported p) and an `nz` carrier-flag column; a dump
# from an older image has neither, and both readers below accept it
kv <- function(line) { x <- strsplit(sub("^# ", "", line), "\t", fixed = TRUE)[[1]]; setNames(x[seq(2, length(x), 2)], x[seq(1, length(x), 2)]) }
h2 <- if (length(l2) > 1 && startsWith(l2[2], "# kernel")) kv(l2[2]) else NULL
d <- read.table(f, header = TRUE, skip = if (is.null(h2)) 1 else 2, sep = "\t")
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
cat(sprintf("   composed p (dense kernel, this image): %.4g\n", p))
if (!is.null(h2)) {
  cat(sprintf("   RUN: kernel=%s reported p=%s pNA=%s | root1=%s/%s root2=%s/%s qsaddle=%s bound=%s failReason=%s | NAmu=%s NAsigma=%s nNB=%s var2=%s varRatio=%s\n",
              h2["kernel"], h2["p"], h2["pNA"], h2["root1"], h2["niter1"], h2["root2"], h2["niter2"], h2["qsaddle"], h2["bound"], h2["failReason"],
              h2["NAmu"], h2["NAsigma"], h2["nNB"], h2["var2"], h2["varRatio"]))
  nz <- d$nz == 1
  if (!is.null(d$g)) {
    G <- d$g; Gc <- G[nz]
    cat(sprintf("   raw vector: max g=%.4g | carriers with g >= 1e-3*max: %d, >= 1e-6*max: %d, >= 1e-9*max: %d | sum g=%.4g | mu at max g=%.2e y=%d\n",
                max(Gc), sum(Gc >= 1e-3 * max(Gc)), sum(Gc >= 1e-6 * max(Gc)), sum(Gc >= 1e-9 * max(Gc)), sum(Gc), mu[nz][which.max(Gc)], d$y[nz][which.max(Gc)]))
  }
  top <- sort(g[nz], decreasing = TRUE)
  cat(sprintf("   carriers flagged=%d | largest g: %s | share of gpos in the largest: %.3f | mu of the largest: %.2e (y=%d)\n",
              sum(nz), paste(sprintf("%.3g", head(top, 6)), collapse = " "), top[1] / gpos, mu[nz][which.max(g[nz])], d$y[nz][which.max(g[nz])]))
  pNA <- as.numeric(h2["pNA"])
  if (identical(unname(h2["kernel"]), "fast")) {
    gNB <- g[nz]; muNB <- mu[nz]; gNA <- g[!nz]; muNA <- mu[!nz]
    NAmu <- as.numeric(h2["NAmu"]); NAsigma <- as.numeric(h2["NAsigma"])
    rf <- getroot_K1_fast_Binom(0, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol, 1000L)
    if (rf$Isconverge) {
      K <- Korg_fast_Binom(rf$root, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma)
      K2 <- K2_fast_Binom(rf$root, mu, g, gNA, gNB, muNA, muNB, NAmu, NAsigma)
      w <- sign(rf$root) * sqrt(2 * (rf$root * q - K)); v <- rf$root * sqrt(K2)
      cat(sprintf("   FAST replay: root(q)=%.6g niter=%d | K=%.8g t*q=%.8g t*q-K=%.6g | K2=%.4g w=%.4g v=%.4g Ztest=%.4g | K1-q at root=%.3g\n",
                  rf$root, rf$niter, K, rf$root * q, rf$root * q - K, K2, w, v, w + log(v / w) / w,
                  K1_adj_fast_Binom(rf$root, mu, g, q, gNA, gNB, muNA, muNB, NAmu, NAsigma)))
      # the same quantities from the rows alone (the dense CGF over every sample, no Gaussian remainder)
      cat(sprintf("   DENSE at the fast root: K=%.8g K2=%.4g | Gaussian remainder from the rows: NAmu=%.6g NAsigma=%.4g (line: %s %s)\n",
                  Korg_Binom(rf$root, mu, g), K2_Binom(rf$root, mu, g), sum(gNA * muNA), sum(muNA * (1 - muNA) * gNA^2), h2["NAmu"], h2["NAsigma"]))
    }
    pf <- SPA_binary_fast(mu, g, q, qinv, pNA, FALSE, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol)
    cat(sprintf("   composed p (fast kernel, this image): %.4g Isconverge=%s\n", pf$pvalue, pf$Isconverge))
  }
  # a reference for the tail: the carriers as exact Bernoullis, the rest as the
  # Gaussian remainder (the fast kernel's own model), 2e6 draws in chunks
  set.seed(1)
  gNB <- g[nz]; muNB <- mu[nz]; gNA <- g[!nz]; muNA <- mu[!nz]
  rm <- sum(gNA * muNA); rs <- sqrt(sum(muNA * (1 - muNA) * gNA^2))
  ndraw <- 2e6; chunk <- 2e5; up <- 0; lo <- 0
  for (b in seq_len(ndraw / chunk)) {
    Y <- matrix(rbinom(length(gNB) * chunk, 1, muNB), nrow = length(gNB))
    S <- as.vector(crossprod(Y, gNB)) + rnorm(chunk, rm, rs)
    up <- up + sum(S >= q); lo <- lo + sum(S <= qinv)
  }
  fmt <- function(k) if (k == 0) sprintf("< %.1e", 1 / ndraw) else sprintf("%.3g (+-%.1g)", k / ndraw, sqrt(k) / ndraw)
  cat(sprintf("   MONTE CARLO reference (%.0e draws): P(S>=q)=%s P(S<=qinv)=%s two-sided=%s\n", ndraw, fmt(up), fmt(lo), fmt(up + lo)))
}
R
for f in $(find dumps -name 'spaFail_*.tsv' -o -name 'spaBelow_*.tsv' | sort -V); do Rscript check.R "$f"; done
echo done
