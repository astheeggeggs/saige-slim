# LEDGER #158 regression guard, finished build: shipped answer vs SKAT's own
# kernel on IDENTICAL (Score, Phi, r.corr), paired per cell, over 32 replicate
# nulls x 3 trait arms.
#
# Controls, reported alongside every result (never a bare "0 misses"):
#   NEGATIVE  a planted miss must be counted (the detector fires)
#   POSITIVE  the threshold-sensitivity curve must show real straddles at some
#             threshold, or the detector is not seeing this data at all
suppressPackageStartupMessages({})
W <- commandArgs(TRUE)[1]
TH <- 0.05; TAIL <- 1e-5
fs <- sort(Sys.glob(file.path(W, "runs", "*.cells.rds")))
fs <- fs[!grepl("PILOT", fs)]
cat("arms:", length(fs), "\n\n")

arm_of <- function(f) {
  b <- sub("\\.cells\\.rds$", "", basename(f))
  if (grepl("^y_null", b)) "quantitative" else if (grepl("^y_nb01", b)) "binary K=0.01" else "binary K=0.002"
}
rows <- list(); allp <- list()
for (f in fs) {
  ce <- try(readRDS(f), silent = TRUE)
  if (inherits(ce, "try-error") || !length(ce)) { cat("SKIP unreadable:", basename(f), "\n"); next }
  p  <- vapply(ce, function(x) if (is.null(x$p)) NA_real_ else as.numeric(x$p), 0)
  pu <- vapply(ce, function(x) if (is.null(x$p_upstream)) NA_real_ else as.numeric(x$p_upstream), 0)
  m  <- vapply(ce, function(x) if (is.null(x$method)) NA_character_ else as.character(x$method), "")
  q  <- vapply(ce, function(x) length(x$Score), 0L)
  ok <- !is.na(p) & !is.na(pu) & p > 0 & pu > 0
  scr <- ok & m == "liu-screen"
  rows[[length(rows) + 1L]] <- data.frame(
    run = sub("\\.cells\\.rds$", "", basename(f)), arm = arm_of(f),
    cells = length(ce), paired = sum(ok),
    n_screen = sum(m == "liu-screen", na.rm = TRUE),
    n_davies = sum(m == "davies", na.rm = TRUE),
    n_single = sum(m == "single-rho", na.rm = TRUE),
    miss = sum(ok & p > TH & pu <= TH),
    miss_p075 = sum(ok & p > 0.075 & pu <= 0.075),
    miss_p10 = sum(ok & p > 0.10 & pu <= 0.10),
    margin = if (any(scr)) min(pu[scr]) / TH else NA_real_,
    maxratio = if (any(ok)) max(abs(log10(p[ok] / pu[ok]))) else NA_real_,
    tail_cells = sum(ok & pu < TAIL),
    tail_maxratio = if (any(ok & pu < TAIL)) max(abs(log10(p[ok & pu < TAIL] / pu[ok & pu < TAIL]))) else NA_real_,
    qmax = max(q), q_le10 = sum(q <= 10),
    stringsAsFactors = FALSE)
  allp[[length(allp) + 1L]] <- data.frame(p = p[ok], pu = pu[ok], m = m[ok], q = q[ok], arm = arm_of(f))
}
df <- do.call(rbind, rows); A <- do.call(rbind, allp)
cat("=== PER ARM (aggregated over replicates) ===\n")
agg <- do.call(rbind, lapply(split(df, df$arm), function(d) data.frame(
  arm = d$arm[1], runs = nrow(d), cells = sum(d$cells), paired = sum(d$paired),
  screened = sum(d$n_screen), davies = sum(d$n_davies), single = sum(d$n_single),
  MISSES = sum(d$miss), miss_at_0.075 = sum(d$miss_p075), miss_at_0.10 = sum(d$miss_p10),
  min_margin = signif(min(d$margin, na.rm = TRUE), 5),
  max_log10_ratio = signif(max(d$maxratio, na.rm = TRUE), 4),
  tail_cells = sum(d$tail_cells),
  tail_max_log10_ratio = if (all(is.na(d$tail_maxratio))) NA else signif(max(d$tail_maxratio, na.rm = TRUE), 4),
  stringsAsFactors = FALSE)))
print(agg, row.names = FALSE)

cat("\n=== TOTAL ===\n")
cat("paired cells:", sum(df$paired), " screened:", sum(df$n_screen),
    " davies:", sum(df$n_davies), " single-rho:", sum(df$n_single), "\n")
cat("SCREEN MISSES at 0.05:", sum(df$miss), "\n")
cat("tightest margin (min upstream p among screened / 0.05):", signif(min(df$margin, na.rm = TRUE), 5), "\n")
cat("max |log10(shipped/upstream)| over every paired cell:", signif(max(df$maxratio, na.rm = TRUE), 4),
    sprintf("(= factor %.4f)\n", 10^max(df$maxratio, na.rm = TRUE)))
cat("cells with upstream p < 1e-5 (the saddlepoint tail's regime):", sum(df$tail_cells), "\n")

cat("\n=== NEGATIVE CONTROL: the detector must count a planted miss ===\n")
cand <- which(A$pu <= TH & A$pu > TH * 0.5)
if (length(cand)) {
  Ap <- A; Ap$p[cand[1]] <- TH * 1.2
  cat("planted 1 miss (upstream p =", signif(A$pu[cand[1]], 4), ") -> counted:",
      sum(Ap$p > TH & Ap$pu <= TH), " (must be 1)\n")
} else cat("no candidate cell to plant on -- CONTROL DID NOT RUN\n")
cat("real miss count on the same pooled data:", sum(A$p > TH & A$pu <= TH), "\n")

cat("\n=== POSITIVE CONTROL: threshold sensitivity (real straddles exist) ===\n")
for (th in c(0.05, 0.075, 0.10, 0.15, 0.25)) {
  cat(sprintf("  threshold %.3f -> %d miss(es) over %d cells with upstream p <= th\n",
              th, sum(A$p > th & A$pu <= th), sum(A$pu <= th)))
}
cat("\n(If every row above is 0 the detector is blind to this data and the\n",
    " zero at 0.05 means nothing -- read it with the negative control.)\n", sep = "")

saveRDS(list(per_run = df, pooled = A), file.path(W, "guard_summary.rds"))
cat("\nwrote", file.path(W, "guard_summary.rds"), "\n")
