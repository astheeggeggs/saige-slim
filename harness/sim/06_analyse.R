# Compare SAIGE step-2 output with the simulated truth.
#   negative control: the 494 non-causal genes (and the y_null trait) must give
#                     uniform p-values -- lambda ~ 1, KS not significant.
#   positive control: the 55 causal genes must rank at the top, and the h2=0.20
#                     traits must reach exome-wide significance.
args <- commandArgs(TRUE)
d <- Sys.getenv("D", "harness/work/sim")
cols <- if (length(args)) args else c("y_quant", "y_binary", "y_quant_pc", "y_binary_pc", "y_null")

truth <- read.delim(file.path(d, "truth_genes.tsv"), stringsAsFactors = FALSE)
vars <- NULL   # 2.3M rows: read once, not once per trait
lam <- function(p) {
  p <- p[is.finite(p) & p > 0]
  signif(median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1, lower.tail = FALSE), 4)
}

for (col in cols) {
  f <- file.path(d, paste0("step2_group_", col, ".txt"))
  if (!file.exists(f)) next
  g <- read.delim(f, stringsAsFactors = FALSE)
  cau <- g[g$Group == "Cauchy", c("Region", "Pvalue")]
  names(cau)[2] <- "P"
  m <- merge(cau, truth, by.x = "Region", by.y = "gene")
  nul <- m[m$causal == 0, ]
  pos <- m[m$causal == 1, ]
  ks <- suppressWarnings(ks.test(nul$P, "punif"))
  nt <- nrow(g[g$Group != "Cauchy", ])
  bonf <- 0.05 / nt
  cat(sprintf(
    "\n=== %s : gene-based (ACAT-O per region) ===\n%d regions (%d causal), %d sub-tests\n",
    col, nrow(m), nrow(pos), nt))
  # lambda from the MEDIAN is meaningless for the ACAT-O row: >25% of regions
  # return p exactly 1.0 (the Cauchy combination saturates in the upper tail).
  # What matters for association testing is the LOWER tail, so report that.
  tailcal <- function(p, a) sprintf("%.4f", mean(p < a))
  cat(sprintf("  NEGATIVE control (non-causal genes, n=%d): observed vs nominal tail\n", nrow(nul)))
  cat(sprintf("     P(p<0.05)=%s [0.05]  P(p<0.01)=%s [0.01]  P(p<0.001)=%s [0.001]  min p=%.2e\n",
              tailcal(nul$P, 0.05), tailcal(nul$P, 0.01), tailcal(nul$P, 0.001), min(nul$P)))
  cat(sprintf("     (median-based lambda=%s, saturated: %.0f%% of rows have p==1)\n",
              lam(nul$P), 100 * mean(nul$P >= 1)))
  cat(sprintf("  POSITIVE control (causal genes)   : lambda=%s  min p=%.2e  n<%.1e: %d\n",
              lam(pos$P), min(pos$P), bonf, sum(pos$P < bonf)))
  cat(sprintf("  causal genes in the top 20 by p   : %d / 20   (rank-sum p=%.2g)\n",
              sum(m$causal[order(m$P)][1:20]),
              suppressWarnings(wilcox.test(pos$P, nul$P, alternative = "less")$p.value)))
  o <- m[order(m$P), ][1:8, c("Region", "P", "causal", "h2", "NCP", "n_causal_var")]
  print(o, row.names = FALSE, digits = 3)

  # which annotation mask carries the signal in the causal genes?
  gg <- merge(g[g$Group != "Cauchy", ], truth, by.x = "Region", by.y = "gene")
  cs <- gg[gg$causal == 1, ]
  cat("  min p by annotation mask (causal genes):\n")
  print(tapply(cs$Pvalue, cs$Group, function(x) signif(min(x, na.rm = TRUE), 3)))
  cat("  sub-test lambda by annotation mask (non-causal genes; these are NOT saturated):\n")
  print(tapply(gg$Pvalue[gg$causal == 0], gg$Group[gg$causal == 0], lam))
  nn <- gg$Pvalue[gg$causal == 0]
  cat(sprintf("  sub-test tail (non-causal): P(p<0.05)=%.4f P(p<0.01)=%.4f P(p<0.001)=%.4f  n=%d\n",
              mean(nn < 0.05, na.rm = TRUE), mean(nn < 0.01, na.rm = TRUE),
              mean(nn < 0.001, na.rm = TRUE), sum(!is.na(nn))))

  fs <- file.path(d, paste0("step2_single_", col, ".txt"))
  if (file.exists(fs)) {
    if (file.info(fs)$size < 1e5) { cat("  (single-variant run still in progress)\n"); next }
    s <- read.delim(fs, stringsAsFactors = FALSE)
    if (is.null(vars)) {
      vars <- read.delim(gzfile(file.path(d, "variants.tsv.gz")), stringsAsFactors = FALSE)
      vars <- vars[, c("ID", "beta", "anno", "causal_gene")]
    }
    s <- merge(s, vars, by.x = "MarkerID", by.y = "ID")
    cat(sprintf("\n=== %s : single variant (%d markers with MAC>=20) ===\n", col, nrow(s)))
    cat(sprintf("  lambda_GC (all)      = %s\n", lam(s$p.value)))
    cat(sprintf("  lambda_GC (beta==0)  = %s   <- negative control\n",
                lam(s$p.value[s$beta == 0])))
    nz <- s[s$beta != 0, ]
    if (nrow(nz) > 2) {
      cat(sprintf("  causal markers tested: %d; cor(BETA_hat, beta_true) = %.3f (Pearson), %.3f (Spearman)\n",
                  nrow(nz), cor(nz$BETA, nz$beta), cor(nz$BETA, nz$beta, method = "spearman")))
      cat(sprintf("  min p among causal = %.2e; min p among null = %.2e\n",
                  min(nz$p.value), min(s$p.value[s$beta == 0])))
    }
  }
}
