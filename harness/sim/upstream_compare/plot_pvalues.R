#!/usr/bin/env Rscript
# -log10 p, upstream main (x) vs strip (y), from compare_pvalues.R output.
#   Rscript plot_pvalues.R <out.png> <title> <pcol> <cmp1.tsv> [<cmp2.tsv> ...]
# One panel per input file (its basename minus .tsv is the panel title). Points are
# coloured by the TAG compare_pvalues.R assigned (none / tiny / AF>0.5 / UNPREDICTED);
# rows present in one arm only (from the .only_*.tsv siblings) are drawn as rug marks on
# that arm's axis; |dlog10 p| > 1 rows are labelled. Base graphics -- ggplot2 is not in
# the pinned env.
args <- commandArgs(TRUE)
out <- args[1]; title <- args[2]; pcol <- args[3]; files <- args[-(1:3)]
cols <- c(none = "#4d4d4d80", tiny = "#1f77b4", small = "#2ca02c", VR = "#17becf", "VR-mix" = "#9edae5", "SKATO=1" = "#bcbd22", Cauchy = "#8c564b", "Cauchy=1" = "#e377c2", "AF>0.5" = "#ff7f0e", UNPREDICTED = "#d62728")
n <- length(files); nc <- if (n <= 2) n else if (n <= 6) 3 else 4; nr <- ceiling(n / nc)
png(out, width = 5.2 * nc, height = 5.2 * nr + 0.6, units = "in", res = 130)
par(mfrow = c(nr, nc), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2.2, 0), pty = "s")
for (f in files) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  x <- d[[paste0("nlog10_", pcol, "_main")]]; y <- d[[paste0("nlog10_", pcol, "_strip")]]
  ok <- is.finite(x) & is.finite(y)
  lim <- c(0, max(c(x[ok], y[ok]), 1, na.rm = TRUE) * 1.04)
  tg <- d$tag; tg[!tg %in% names(cols)] <- "UNPREDICTED"
  o <- order(tg != "none")   # draw identical rows first so departures sit on top
  plot(x[o], y[o], xlim = lim, ylim = lim, pch = 16, cex = 0.55, col = cols[tg[o]],
       xlab = bquote(-log[10] ~ p ~ "(upstream main 365deeff)"), ylab = bquote(-log[10] ~ p ~ "(strip)"),
       main = sub("\\.tsv$", "", basename(f)), cex.main = 0.95)
  abline(0, 1, col = "grey40", lty = 2)
  for (side in c("main", "strip")) {
    of <- paste0(f, ".only_", side, ".tsv")
    if (file.exists(of)) {
      z <- read.delim(of, stringsAsFactors = FALSE); pc <- if (pcol %in% names(z)) pcol else names(z)[grep("^p", names(z), ignore.case = TRUE)[1]]
      v <- -log10(as.numeric(z[[pc]])); v <- v[is.finite(v)]
      if (length(v)) rug(v, side = if (side == "main") 1 else 2, col = "#9467bd", lwd = 1.2)
    }
  }
  big <- which(ok & abs(y - x) > 1); big <- big[order(-abs(y - x)[big])][seq_len(min(5, length(big)))]
  if (length(big)) text(x[big], y[big], labels = d$key[big], cex = 0.55, pos = 4, col = "#d62728")
  tb <- table(factor(tg, levels = names(cols)))
  legend("topleft", bty = "n", cex = 0.72, pch = 16, col = cols,
         legend = sprintf("%s (%d)", names(cols), as.integer(tb)))
  ok2 <- ok
  legend("bottomright", bty = "n", cex = 0.72, legend = c(
    sprintf("n = %d shared, Pearson r = %.4f", sum(ok2), if (sum(ok2) > 2) cor(x[ok2], y[ok2]) else NA),
    sprintf("identical p: %.1f%%   max |dlog10| = %.2g", 100 * mean((d[[paste0(pcol, "_main")]] == d[[paste0(pcol, "_strip")]])[ok2]), max(abs(y[ok2] - x[ok2])))))
}
mtext(title, outer = TRUE, cex = 1.1, font = 2)
dev.off()
cat("wrote", out, "\n")
