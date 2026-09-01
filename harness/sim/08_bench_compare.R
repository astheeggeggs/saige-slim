# Do the two readers agree on the association results?
#
# Byte-identity is REPORTED here, not asserted: it is a claim about two different
# code paths, which the repo's standing rule says not to make in general (1-ULP
# differences reach the printed digits on some platforms and not others). The
# gate is: same markers, same order, statistics equal to float tolerance.
a <- "harness/work/readerbench/bed_1.txt"
b <- "harness/work/readerbench/pgen_1.txt"
A <- read.delim(a, stringsAsFactors = FALSE)
B <- read.delim(b, stringsAsFactors = FALSE)

cat(sprintf("rows      bed %d   pgen %d\n", nrow(A), nrow(B)))
cat(sprintf("columns   identical: %s\n", identical(names(A), names(B))))
stopifnot(identical(names(A), names(B)))

# Same marker set, in the same order? A reader that silently reordered or dropped
# markers would still produce plausible-looking numbers.
same_ids <- identical(A$MarkerID, B$MarkerID)
cat(sprintf("marker IDs identical and in the same order: %s\n", same_ids))
if (!same_ids) {
  cat(sprintf("  in bed only: %d | in pgen only: %d | shared: %d\n",
              length(setdiff(A$MarkerID, B$MarkerID)),
              length(setdiff(B$MarkerID, A$MarkerID)),
              length(intersect(A$MarkerID, B$MarkerID))))
  k <- intersect(A$MarkerID, B$MarkerID)
  A <- A[match(k, A$MarkerID), ]; B <- B[match(k, B$MarkerID), ]
}

num <- names(A)[vapply(A, is.numeric, TRUE)]
cat("\ncolumn                 max |abs diff|   max |rel diff|   n differing\n")
worst <- 0
for (cl in num) {
  x <- A[[cl]]; y <- B[[cl]]
  ok <- is.finite(x) & is.finite(y)
  ad <- abs(x[ok] - y[ok])
  rd <- ad / pmax(abs(x[ok]), abs(y[ok]), .Machine$double.eps)
  nd <- sum(ad > 0)
  worst <- max(worst, if (length(rd)) max(rd) else 0)
  cat(sprintf("%-22s %14.3e %16.3e %12d\n", cl,
              if (length(ad)) max(ad) else 0, if (length(rd)) max(rd) else 0, nd))
}
cat(sprintf("\nworst relative difference across all numeric columns: %.3e\n", worst))

# The p-value is what anyone actually reads off this table.
pa <- A$p.value; pb <- B$p.value
ok <- is.finite(pa) & is.finite(pb) & pa > 0 & pb > 0
cat(sprintf("p-values: max |log10 ratio| = %.3e over %d markers\n",
            max(abs(log10(pa[ok]) - log10(pb[ok]))), sum(ok)))
cat(sprintf("rank correlation of -log10 p: %.10f\n",
            cor(-log10(pa[ok]), -log10(pb[ok]), method = "spearman")))

bytes <- tools::md5sum(c(a, b))
cat(sprintf("\nmd5  bed  %s\nmd5  pgen %s\nbyte-identical: %s\n",
            bytes[1], bytes[2], bytes[1] == bytes[2]))
cat(sprintf("VERDICT: %s\n", if (worst < 1e-9 && same_ids)
    "same markers, same order, statistics agree to < 1e-9 relative" else "DIFFERENCES PRESENT - inspect above"))
