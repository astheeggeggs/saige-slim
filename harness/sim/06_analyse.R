# Compare SAIGE step-2 output with the simulated truth.
#   negative control: the 494 non-causal genes (and the y_null trait) must give
#                     uniform p-values -- lambda ~ 1, KS not significant.
#   positive control: the 55 causal genes must rank at the top, and the h2=0.20
#                     traits must reach exome-wide significance.
args <- commandArgs(TRUE)
d <- Sys.getenv("D", "harness/work/sim")
cols <- if (length(args)) args else c("y_quant", "y_binary", "y_quant_pc", "y_binary_pc", "y_null")

truth <- read.delim(file.path(d, "truth_genes.tsv"), stringsAsFactors = FALSE)
# truth_genes.tsv carries NCP_quant/NCP_binary since the two arms are scaled
# separately (a binary trait at prevalence K buys only phi(t)^2/(K(1-K)) of the
# non-centrality). Older files have a single NCP column; accept both.
ncpcol <- function(col) {
  if ("NCP" %in% names(truth)) return("NCP")
  if (grepl("binary", col) && "NCP_binary" %in% names(truth)) "NCP_binary" else "NCP_quant"
}
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
  nc <- ncpcol(col)
  o <- m[order(m$P), ][1:8, c("Region", "P", "causal", "h2", nc, "n_causal_var")]
  print(o, row.names = FALSE, digits = 3)

  # The control the NCP ladder exists to enable. A hit count saturates once N is
  # large enough that every causal gene is found; a rank correlation between the
  # non-centrality each gene was BUILT with and the significance it actually
  # attains keeps discriminating, and it is the assertion amendment A4 asks for.
  # Predicted: strongly positive (rho -> 1 as the ladder spans the power curve),
  # and it must NOT be positive among non-causal genes, where NCP is 0 for all.
  if (nrow(pos) >= 5 && length(unique(pos[[nc]])) > 1) {
    ok <- is.finite(pos$P) & pos$P > 0
    rho <- suppressWarnings(cor(pos[[nc]][ok], -log10(pos$P[ok]), method = "spearman"))
    ct <- suppressWarnings(cor.test(pos[[nc]][ok], -log10(pos$P[ok]), method = "spearman"))
    cat(sprintf("  LADDER control: Spearman(%s, -log10 p) over %d causal genes = %.3f (p=%.3g)\n",
                nc, sum(ok), rho, ct$p.value))
    cat(sprintf("     %s range over causal genes: %.2f - %.2f\n",
                nc, min(pos[[nc]]), max(pos[[nc]])))
  } else {
    cat(sprintf("  LADDER control: skipped (%d causal genes, %d distinct %s)\n",
                nrow(pos), length(unique(pos[[nc]])), nc))
  }

  # A1: max_MAF = 1e-4 is STRUCTURALLY DEGENERATE. MAC <= 10 is the ultra-rare
  # collapse threshold (src/Main.cpp), and at 2N = 100,000 that is exactly
  # MAF <= 1e-4 -- so every marker in that stratum collapses into ONE
  # pseudo-marker and SKAT, BURDEN and SKAT-O become the same test. Measured
  # here: 1042/1042 cells identical to the last digit, 0 with >= 4 real markers.
  # Any BURDEN-vs-SKAT or SKAT-O comparison must therefore be restricted to
  # cells that actually contain several uncollapsed markers, or it is comparing
  # a number with itself and will report a spurious dead heat.
  sub <- g[g$Group != "Cauchy", ]
  cat("  A1 degeneracy by max_MAF (SKAT==BURDEN exactly / cells with >=4 real markers):\n")
  for (mm in sort(unique(sub$max_MAF))) {
    ss <- sub[sub$max_MAF == mm, ]
    cat(sprintf("     max_MAF=%-8s %4d cells  identical %4d (%5.1f%%)  usable %4d (%5.1f%%)\n",
                mm, nrow(ss), sum(ss$Pvalue_Burden == ss$Pvalue_SKAT),
                100 * mean(ss$Pvalue_Burden == ss$Pvalue_SKAT),
                sum(ss$Number_rare >= 4), 100 * mean(ss$Number_rare >= 4)))
  }
  usable <- sub[sub$Number_rare >= 4, ]
  cat(sprintf("  -> ranking controls below use %d of %d sub-test cells\n",
              nrow(usable), nrow(sub)))

  # DIRECTIONALITY control (plan Stage 7). BURDEN sums signed effects and SKAT
  # sums squares, so a gene whose causal effects mostly share a sign should be
  # found more easily by BURDEN, and a gene with mixed signs by SKAT. This is a
  # prediction about the RANKING of two tests on the same data, which is far
  # harder to satisfy by accident than "the causal genes come out significant",
  # and it is the reason the generator bothers to make class-2 effects 90%
  # same-sign.
  #
  # Directionality per gene is measured from the simulated betas as
  # |sum(beta)| / sum(|beta|): 1 = perfectly directional, 0 = perfectly
  # cancelling. Prediction: Spearman(dir, log10(P_SKAT / P_BURDEN)) > 0.
  # Restricted to cells with >= 4 uncollapsed markers, because below that the
  # two tests are literally the same number (see the A1 table above).
  if (is.null(vars)) {
    vf <- file.path(d, "variants.tsv.gz")
    if (file.exists(vf)) vars <<- read.delim(gzfile(vf), stringsAsFactors = FALSE)
  }
  if (!is.null(vars) && nrow(usable) > 10) {
    cv <- vars[vars$beta != 0 & vars$gene != ".", c("gene", "beta")]
    dir <- aggregate(beta ~ gene, cv, function(b) abs(sum(b)) / sum(abs(b)))
    names(dir)[2] <- "dirscore"
    nvar <- aggregate(beta ~ gene, cv, length); names(nvar)[2] <- "n_causal"
    dir <- merge(dir, nvar)
    dir <- dir[dir$n_causal >= 2, ]          # a 1-variant gene is trivially "directional"
    u <- merge(usable, dir, by.x = "Region", by.y = "gene")
    u <- u[is.finite(u$Pvalue_Burden) & is.finite(u$Pvalue_SKAT) &
             u$Pvalue_Burden > 0 & u$Pvalue_SKAT > 0, ]
    if (nrow(u) >= 8) {
      lr <- log10(u$Pvalue_SKAT) - log10(u$Pvalue_Burden)
      rho <- suppressWarnings(cor(u$dirscore, lr, method = "spearman"))
      ct <- suppressWarnings(cor.test(u$dirscore, lr, method = "spearman"))
      cat(sprintf("  DIRECTIONALITY control: Spearman(dir, log10 P_SKAT/P_BURDEN) = %.3f (p=%.3g) over %d cells in %d genes\n",
                  rho, ct$p.value, nrow(u), length(unique(u$Region))))
      cat(sprintf("     predicted > 0 (directional genes favour BURDEN). dir range %.2f-%.2f\n",
                  min(u$dirscore), max(u$dirscore)))
      hi <- u[u$dirscore >= 0.8, ]; lo <- u[u$dirscore <= 0.4, ]
      if (nrow(hi) && nrow(lo)) {
        cat(sprintf("     BURDEN beats SKAT in %5.1f%% of directional cells (dir>=0.8, n=%d) vs %5.1f%% of mixed (dir<=0.4, n=%d)\n",
                    100 * mean(hi$Pvalue_Burden < hi$Pvalue_SKAT), nrow(hi),
                    100 * mean(lo$Pvalue_Burden < lo$Pvalue_SKAT), nrow(lo)))
      }
    } else {
      cat(sprintf("  DIRECTIONALITY control: skipped (%d usable cells with >=2 causal variants)\n", nrow(u)))
    }
  }

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
