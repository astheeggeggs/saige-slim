# Upstream SAIGE (main, 365deeff) vs the strip branch, same fixture, same model,
# same GRM, single-variant, --minMAC=20.
#
# Parity is NOT expected: the strip branch carries deliberate fixes. For a
# QUANTITATIVE single-variant run the relevant ones are LEDGER #28 (variance
# mis-centred by scoreTestFast_noadjCov on the AF > 0.5 path) and #29 (markers at
# exactly --minMAC dropped by a float round trip in gate 1). #29 predicts the two
# builds may not even test the SAME MARKER SET, so check that before comparing
# any statistic -- a correlation computed over a silently different marker set
# would be meaningless.
A <- read.delim("/Users/dpalmer/Repositories/SAIGE-sim/harness/work/readerbench/bed_1.txt",
                stringsAsFactors = FALSE)          # strip
B <- read.delim("/Users/dpalmer/Repositories/SAIGE-main/bench/bed_1.txt",
                stringsAsFactors = FALSE)          # main
cat(sprintf("markers tested   strip %d   main %d   (difference %+d)\n",
            nrow(A), nrow(B), nrow(B) - nrow(A)))

only_s <- setdiff(A$MarkerID, B$MarkerID); only_m <- setdiff(B$MarkerID, A$MarkerID)
cat(sprintf("in strip only: %d | in main only: %d | shared: %d\n",
            length(only_s), length(only_m), length(intersect(A$MarkerID, B$MarkerID))))
if (length(only_s)) {
  x <- A[A$MarkerID %in% only_s, ]
  mac <- pmin(x$AC_Allele2, 2 * x$N - x$AC_Allele2)
  cat(sprintf("  strip-only markers: MAC range %g - %g   (LEDGER #29 predicts exactly minMAC=20)\n",
              min(mac), max(mac)))
  print(table(MAC = mac))
}
if (length(only_m)) {
  x <- B[B$MarkerID %in% only_m, ]
  mac <- pmin(x$AC_Allele2, 2 * x$N - x$AC_Allele2)
  cat(sprintf("  main-only markers: MAC range %g - %g\n", min(mac), max(mac)))
  print(table(MAC = mac))
}

k <- intersect(A$MarkerID, B$MarkerID)
a <- A[match(k, A$MarkerID), ]; b <- B[match(k, B$MarkerID), ]
cat(sprintf("\n== on the %d shared markers ==\n", length(k)))
pa <- a$p.value; pb <- b$p.value
ok <- is.finite(pa) & is.finite(pb) & pa > 0 & pb > 0
la <- -log10(pa[ok]); lb <- -log10(pb[ok])
cat(sprintf("cor(-log10 p)      Pearson %.10f   Spearman %.10f\n",
            cor(la, lb), cor(la, lb, method = "spearman")))
cat(sprintf("identical p-values : %d / %d (%.4f%%)\n",
            sum(pa[ok] == pb[ok]), sum(ok), 100 * mean(pa[ok] == pb[ok])))
d <- abs(la - lb)
cat(sprintf("max |dlog10 p|     %.3e   median over differing %.3e\n",
            max(d), if (any(d > 0)) median(d[d > 0]) else 0))
for (cl in c("BETA", "SE", "Tstat", "var", "AF_Allele2")) {
  if (!cl %in% names(a)) next
  x <- a[[cl]]; y <- b[[cl]]; o <- is.finite(x) & is.finite(y)
  cat(sprintf("  %-11s max |rel diff| %.3e   n differing %d\n", cl,
              max(abs(x[o] - y[o]) / pmax(abs(x[o]), abs(y[o]), .Machine$double.eps)),
              sum(x[o] != y[o])))
}
# Where do the differences sit? #28 is the AF > 0.5 path specifically.
hi <- a$AF_Allele2 > 0.5
cat(sprintf("\nmarkers with AF_Allele2 > 0.5: %d of %d (%.2f%%)\n", sum(hi), nrow(a), 100 * mean(hi)))
if (any(hi) && any(!hi)) {
  dl <- abs(-log10(pa) + log10(pb)); dl[!is.finite(dl)] <- 0
  cat(sprintf("  mean |dlog10 p|  AF<=0.5 %.3e   AF>0.5 %.3e\n", mean(dl[!hi]), mean(dl[hi])))
  cat(sprintf("  n differing      AF<=0.5 %d/%d      AF>0.5 %d/%d\n",
              sum(dl[!hi] > 0), sum(!hi), sum(dl[hi] > 0), sum(hi)))
}
