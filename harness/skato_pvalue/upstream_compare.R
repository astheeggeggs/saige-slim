S <- commandArgs(trailingOnly = TRUE)[1]
fm <- list.files(file.path(S, "upcmp"), "^out_main_[0-9]+[.]txt$", full.names = TRUE)[1]; fs <- list.files(file.path(S, "upcmp"), "^out_strip_[0-9]+[.]txt$", full.names = TRUE)[1]
m <- read.table(fm, header = TRUE, sep = "\t", stringsAsFactors = FALSE); s <- read.table(fs, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("rows upstream", nrow(m), " strip", nrow(s), "\n")
cols <- intersect(c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"), intersect(names(m), names(s)))
j <- merge(m[, c("Region", "Group", "max_MAF", cols)], s[, c("Region", "Group", "max_MAF", cols)], by = c("Region", "Group", "max_MAF"), suffixes = c(".up", ".st"))
j <- j[!grepl("Cauchy", j$Group), ]; cat("joined non-Cauchy cells:", nrow(j), "\n")
for (col in cols) { a <- j[[paste0(col, ".up")]]; b <- j[[paste0(col, ".st")]]; ok <- !is.na(a) & !is.na(b) & a > 0 & b > 0
  d <- abs(log10(a[ok]) - log10(b[ok])); low <- (a[ok] < 1e-5 | b[ok] < 1e-5)
  cat(sprintf("%s: n=%d identical=%d median|dlog10 p|=%.2e p90=%.3f max=%.3f (upstream %.3g vs strip %.3g); cells <1e-5 in either: %d, |dlog10|>0.1 among them: %d\n", col, sum(ok), sum(a[ok] == b[ok]), median(d), quantile(d, 0.9), max(d), a[ok][which.max(d)], b[ok][which.max(d)], sum(low), sum(d > 0.1 & low))) }
