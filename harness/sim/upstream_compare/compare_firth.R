#!/usr/bin/env Rscript
# Firth on/off control, WITHIN one build: does --is_Firth_beta=FALSE move anything
# other than BETA and SE?
#   Rscript compare_firth.R <kind: single|gene> <base.txt> <firthoff.txt> <label> [out.tsv]
# Columns are read as CHARACTER so the comparison is on the printed bytes, not on a
# reparsed double: the claim under test (src/SAIGE_test.cpp:942-1017 -- Firth rewrites
# only t_Beta/t_seBeta, never pval) is a byte claim about the writer's output.
# Prints, for each column, the number of rows whose printed value differs; then the
# BETA/SE movement on the rows that moved.
args <- commandArgs(TRUE)
kind <- args[1]; fa <- args[2]; fb <- args[3]
label <- if (length(args) >= 4) args[4] else ""
out <- if (length(args) >= 5) args[5] else NA
# na.strings = character(0): a literal "NA" in the file must stay the STRING "NA", or the
# byte comparison below silently turns every NA cell into an NA logical and stops.
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE,
                             colClasses = "character", na.strings = character(0))
A <- rd(fa); B <- rd(fb)
cat(sprintf("== FIRTH CONTROL %s (%s) ==\nrows firth-on %d  firth-off %d\n", label, kind, nrow(A), nrow(B)))
stopifnot(identical(names(A), names(B)))
key <- if (kind == "single") function(d) d$MarkerID else
  function(d) paste(d$Region, d$Group, d$max_MAF, sep = "|")
ka <- key(A); kb <- key(B)
if (!identical(ka, kb)) {
  cat(sprintf("ROW SET / ORDER DIFFERS: only-on %d only-off %d\n",
              length(setdiff(ka, kb)), length(setdiff(kb, ka))))
  k <- intersect(ka, kb); A <- A[match(k, ka), ]; B <- B[match(k, kb), ]
} else cat("row set and order identical\n")
moved <- character(0)
for (cl in names(A)) {
  n <- sum(A[[cl]] != B[[cl]])
  if (n > 0) { moved <- c(moved, cl); cat(sprintf("  %-18s %7d / %d rows differ\n", cl, n, nrow(A))) }
}
if (length(moved) == 0) cat("  IDENTICAL: no column moved (whole file)\n")
expected <- c("BETA", "SE", "BETA_Burden", "SE_Burden")
unexpected <- setdiff(moved, expected)
cat(sprintf("columns moved: %s\nUNEXPECTED (anything but BETA/SE): %s\n",
            if (length(moved)) paste(moved, collapse = ",") else "(none)",
            if (length(unexpected)) paste(unexpected, collapse = ",") else "(none) -- control PASSES"))
for (cl in intersect(expected, moved)) {
  x <- suppressWarnings(as.numeric(A[[cl]])); y <- suppressWarnings(as.numeric(B[[cl]]))
  d <- which(A[[cl]] != B[[cl]] & is.finite(x) & is.finite(y))
  r <- y[d] / x[d]
  cat(sprintf("\n%s: %d rows moved; %s ratio (off/on) median %.4f  IQR [%.4f, %.4f]  range [%.4f, %.4f]\n",
              cl, length(d), cl, median(r), quantile(r, .25), quantile(r, .75), min(r), max(r)))
}
if (!is.na(out) && length(moved)) {
  x <- suppressWarnings(as.numeric(A$BETA)); y <- suppressWarnings(as.numeric(B$BETA))
  sx <- suppressWarnings(as.numeric(A$SE));  sy <- suppressWarnings(as.numeric(B$SE))
  res <- data.frame(key = key(A), BETA_firthon = x, BETA_firthoff = y,
                    SE_firthon = sx, SE_firthoff = sy,
                    p = suppressWarnings(as.numeric(A[[if (kind == "single") "p.value" else "Pvalue"]])),
                    firth_fired = A$BETA != B$BETA)
  write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("wrote %s\n", out))
}
cat("\n")
