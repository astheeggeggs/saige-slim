#!/usr/bin/env Rscript
# Compare two step-2 result files, upstream main (A) vs strip (B), row by row.
#   Rscript compare_pvalues.R single <main.txt> <strip.txt> <out.tsv> [label]
#   Rscript compare_pvalues.R gene   <main.txt> <strip.txt> <out.tsv> [label] [vr_main.txt vr_strip.txt]
# With the two variance-ratio files, gene mode tags a cell VR when (SE_Burden_strip /
# SE_Burden_main)^2 equals 1 / (some row's strip/main ratio) to 1e-3: SE_Burden is 1/sqrt(var),
# so that quotient is the variance-ratio quotient the cell used (sparse or null row).
# Marker/gene SET differences are reported BEFORE any statistic (LEDGER #29 predicts the
# sets differ at MAC == --minMAC), then every shared row is written with its
# -log10 p in both arms, the difference, and a TAG naming the predicted cause:
#   AF>0.5   LEDGER #28 (variance mis-centred on the AF > 0.5 fast path)
#   tiny     |dlog10 p| < 1e-5 (rounding-level, e.g. #65 standardised step-1 fit)
#   small    |dlog10 p| < 0.05 (a variance-ratio difference between the arms' step-1 models:
#            LEDGER #53 decode fix moves upstream's category-1 ratio by ~0.14%; reconciled
#            below against the reported `var` column, whose per-category quotient must equal
#            the variance-ratio quotient if that is the whole story)
#   VR       the p-value movement is exactly what the two arms' reported `var` columns predict
#            under a chi-square rescale (|log10 p_pred - log10 p_strip| < 1e-3 with
#            chi2_strip = chi2_main * var_main / var_strip): a variance-ratio difference between
#            the arms' step-1 models (LEDGER #53 on binary; #97's sparse row on quantitative)
#   VR-mix   gene mode: the cell's quotient lies between the arms' variance-ratio row quotients
#            (the cell's markers span MAC categories, so its variance mixes rows)
#   SKATO=1  gene mode: either arm printed a SKAT-O p of exactly 1 (SKAT's saturation, #86)
#   Cauchy   gene mode: the per-gene Cauchy row, an aggregate of the tagged cells; Cauchy=1 when
#            either arm printed exactly 1 (SAIGE's CCT falls back to Bonferroni, LEDGER #86)
#   none     identical
#   UNPREDICTED anything else -- must be explained before numbers are quoted.
args <- commandArgs(TRUE)
kind <- args[1]; fa <- args[2]; fb <- args[3]; out <- args[4]; label <- if (length(args) >= 5) args[5] else ""
vrq <- NULL
if (length(args) >= 7 && file.exists(args[6]) && file.exists(args[7])) {
  va <- read.table(args[6]); vb <- read.table(args[7]); vrq <- vb$V1 / va$V1
  cat("variance-ratio quotients strip/main:", sprintf("%.5f", vrq), "\n")
}
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
unlink(paste0(out, c(".only_main.tsv", ".only_strip.tsv")))   # never leave a stale set-difference file
A <- rd(fa); B <- rd(fb)
cat(sprintf("== %s %s ==\nrows  main %d  strip %d\n", label, kind, nrow(A), nrow(B)))
if (kind == "single") {
  key <- function(d) d$MarkerID
  pcols <- c("p.value")
} else {
  key <- function(d) paste(d$Region, d$Group, d$max_MAF, sep = "|")
  pcols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT")
}
ka <- key(A); kb <- key(B)
stopifnot(!anyDuplicated(ka), !anyDuplicated(kb))
only_a <- setdiff(ka, kb); only_b <- setdiff(kb, ka)
cat(sprintf("only main %d | only strip %d | shared %d\n", length(only_a), length(only_b), length(intersect(ka, kb))))
mac_of <- function(d) if (all(c("AC_Allele2", "N") %in% names(d))) pmin(d$AC_Allele2, 2 * d$N - d$AC_Allele2) else
  if (all(c("AC_Allele2", "N_case", "N_ctrl") %in% names(d))) pmin(d$AC_Allele2, 2 * (d$N_case + d$N_ctrl) - d$AC_Allele2) else rep(NA_real_, nrow(d))
if (kind == "single") {
  if (length(only_b)) { m <- mac_of(B)[kb %in% only_b]; cat("strip-only MAC table (LEDGER #29 predicts all == minMAC):\n"); print(table(m)) }
  if (length(only_a)) { m <- mac_of(A)[ka %in% only_a]; cat("main-only MAC table:\n"); print(table(m)) }
}
k <- intersect(ka, kb); a <- A[match(k, ka), ]; b <- B[match(k, kb), ]
res <- data.frame(key = k, stringsAsFactors = FALSE)
if (kind == "single") {
  res$CHR <- a$CHR; res$POS <- a$POS; res$AF_Allele2 <- a$AF_Allele2; res$MAC <- mac_of(a)
  res$Is.SPA_main <- a$Is.SPA; res$Is.SPA_strip <- b$Is.SPA
} else {
  res$Region <- a$Region; res$Group <- a$Group; res$max_MAF <- a$max_MAF
  res$MAC <- a$MAC; res$Number_rare <- a$Number_rare
}
for (pc in pcols) {
  pa <- as.numeric(a[[pc]]); pb <- as.numeric(b[[pc]])
  la <- -log10(pa); lb <- -log10(pb)
  res[[paste0(pc, "_main")]] <- pa; res[[paste0(pc, "_strip")]] <- pb
  res[[paste0("nlog10_", pc, "_main")]] <- la; res[[paste0("nlog10_", pc, "_strip")]] <- lb
  d <- lb - la; d[is.na(pa) & is.na(pb)] <- 0
  res[[paste0("dlog10_", pc)]] <- d
  ok <- is.finite(la) & is.finite(lb)
  cat(sprintf("%-14s Pearson %.6f Spearman %.6f | identical %d/%d (%.2f%%) | max|d| %.3e | n(|d|>1) %d | NA main %d strip %d\n",
      pc, cor(la[ok], lb[ok]), cor(la[ok], lb[ok], method = "spearman"),
      sum(pa == pb, na.rm = TRUE), sum(ok), 100 * mean((pa == pb)[ok]),
      max(abs(d[ok])), sum(abs(d[ok]) > 1), sum(is.na(pa)), sum(is.na(pb))))
}
# Tag by the primary p column.
d <- abs(res[[paste0("dlog10_", pcols[1])]])
tag <- ifelse(d == 0 | is.na(d), "none", ifelse(d < 1e-5, "tiny", ifelse(d < 0.05, "small", "UNPREDICTED")))
if (kind == "single" && "var" %in% names(a)) {
  q <- as.numeric(b$var) / as.numeric(a$var)
  pa1 <- as.numeric(a[[pcols[1]]]); pb1 <- as.numeric(b[[pcols[1]]])
  chi_main <- qchisq(pa1, 1, lower.tail = FALSE)
  p_pred <- pchisq(chi_main / q, 1, lower.tail = FALSE)
  vr_ok <- is.finite(q) & abs(q - 1) > 1e-6 & is.finite(p_pred) & pb1 > 0 & abs(log10(p_pred) - log10(pb1)) < 1e-3
  res$var_quotient <- q; res$p_pred_from_var <- p_pred
  tag[d > 0 & !is.na(vr_ok) & vr_ok] <- "VR"
}
if (kind == "gene" && !is.null(vrq) && all(c("SE_Burden") %in% names(a))) {
  q <- (as.numeric(b$SE_Burden) / as.numeric(a$SE_Burden))^2   # = var_main / var_strip
  res$var_quotient <- 1 / q
  hit <- sapply(1 / q, function(x) any(abs(x / vrq - 1) < 1e-3))
  tag[d > 0 & !is.na(hit) & hit & abs(1 / q - 1) > 1e-6] <- "VR"
  # a cell whose markers span MAC categories mixes rows: its quotient lies between the extremes
  mix <- is.finite(q) & !hit & 1 / q > min(vrq) - 1e-3 & 1 / q < max(vrq) + 1e-3
  tag[d > 0 & mix & abs(1 / q - 1) > 1e-6] <- "VR-mix"
  cat("gene cells whose SE_Burden^2 quotient matches a variance-ratio row:", sum(hit, na.rm = TRUE), "; lies between the rows (mixed MAC categories):", sum(mix), "of", sum(is.finite(q)), "\n")
}
if (kind == "gene" && "Group" %in% names(a)) {
  # the Cauchy row aggregates the gene's cells, so its movement is theirs; when any cell's p is
  # exactly 1 SAIGE's CCT falls back to Bonferroni (LEDGER #86), which is what a printed 1 is
  # a SKAT-O p of exactly 1 is the point mass LEDGER #86 documents (SKAT's own saturation)
  one <- a$Group != "Cauchy" & (as.numeric(a[[pcols[1]]]) == 1 | as.numeric(b[[pcols[1]]]) == 1)
  tag[d > 0 & one & tag %in% c("small", "UNPREDICTED")] <- "SKATO=1"
  cau <- a$Group == "Cauchy"
  tag[d > 0 & cau] <- ifelse((as.numeric(a[[pcols[1]]]) == 1 | as.numeric(b[[pcols[1]]]) == 1)[d > 0 & cau], "Cauchy=1", "Cauchy")
}
if (kind == "single") tag[d > 0 & !is.na(res$AF_Allele2) & res$AF_Allele2 > 0.5] <- "AF>0.5"
res$tag <- tag
cat("tags: "); print(table(tag))
if (kind == "single") {
  hi <- res$AF_Allele2 > 0.5
  cat(sprintf("AF<=0.5: differing %d/%d   AF>0.5: differing %d/%d\n", sum(d[!hi] > 0), sum(!hi), sum(d[hi] > 0), sum(hi)))
  if ("var" %in% names(a)) {
    # Two-way reconciliation of the arms' variance ratios: on non-SPA rows the reported var is
    # var_score * VR, so var_strip / var_main per MAC category is the VR quotient.
    nsp <- !hi
    if ("Is.SPA" %in% names(a)) nsp <- nsp & tolower(a$Is.SPA) == "false" & tolower(b$Is.SPA) == "false"
    for (cat in c("MAC<=20.5", "MAC>20.5")) {
      sel <- nsp & (if (cat == "MAC<=20.5") res$MAC <= 20.5 else res$MAC > 20.5)
      if (any(sel)) { q <- as.numeric(b$var[sel]) / as.numeric(a$var[sel])
        cat(sprintf("  var_strip/var_main non-SPA AF<=0.5 %-9s n=%6d median %.7f  range [%.7f, %.7f]\n", cat, sum(sel), median(q), min(q), max(q))) }
    }
  }
  for (cl in c("BETA", "SE", "Tstat", "var", "AF_Allele2")) if (cl %in% names(a)) {
    x <- as.numeric(a[[cl]]); y <- as.numeric(b[[cl]]); o <- is.finite(x) & is.finite(y)
    cat(sprintf("  %-11s max|rel diff| %.3e  n differing %d\n", cl,
        max(abs(x[o] - y[o]) / pmax(abs(x[o]), abs(y[o]), .Machine$double.eps)), sum(x[o] != y[o])))
  }
}
if (length(only_a)) write.table(A[ka %in% only_a, ], paste0(out, ".only_main.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(only_b)) write.table(B[kb %in% only_b, ], paste0(out, ".only_strip.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n\n", out, nrow(res)))
