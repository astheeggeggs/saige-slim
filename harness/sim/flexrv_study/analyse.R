# FlexRV simulation study analysis (docs/design/flexrv-sim-study.md).
# Controls FIRST, each with its expectation printed beside it; then the
# tables. Reads results/<trait>_<arm>_<method>_c<collapse>_r1.txt written by
# run.sh and truth_genes.tsv from the fixture.
# usage: Rscript analyse.R [results dir] [fixture dir]
args <- commandArgs(TRUE)
R <- if (length(args) >= 1) args[1] else "/Users/dpalmer/Repositories/SAIGE/harness/work/flexrv_study/results"
D <- if (length(args) >= 2) args[2] else "/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k"
rd <- function(p) utils::read.table(p, header = TRUE, sep = "\t", colClasses = "character",
                                    quote = "", comment.char = "", na.strings = character(0))
truth <- rd(file.path(D, "truth_genes.tsv"))
truth$causal <- as.integer(truth$causal)
truth$directional <- as.integer(truth$directional)
truth$NCP_quant <- as.numeric(truth$NCP_quant)
truth$NCP_binary <- as.numeric(truth$NCP_binary)
if (!"link" %in% colnames(truth)) truth$link <- "."
files <- list.files(R, pattern = "_r1\\.txt$", full.names = TRUE)
parse_tag <- function(f) {
  t <- sub("_r1\\.txt$", "", basename(f))
  m <- regmatches(t, regexec("^(y_[a-z_]+?)_(informative|realistic|uninformative|synonymous|plain)_([a-z0-9]+)_c([0-9]+)$", t))[[1]]
  list(trait = m[2], arm = m[3], method = m[4], collapse = m[5])
}
# the gene-level p per file: the pooled Cauchy row (FlexRV) or the cell row
gene_p <- function(f) {
  d <- rd(f)
  if ("Weight" %in% colnames(d) || any(d$Group == "Cauchy")) {
    z <- d[d$Group == "Cauchy", ]
    p <- if ("Pvalue" %in% colnames(z)) z$Pvalue else z$Pvalue_Burden
    r <- d[d$Group != "Cauchy", ]
    # the oracle: the best single set per gene, Bonferroni over the sets
    pr <- if ("Pvalue" %in% colnames(r)) r$Pvalue else r$Pvalue_Burden
    orc <- tapply(as.numeric(pr), r$Region, function(x) min(1, min(x, na.rm = TRUE) * sum(!is.na(x))))
    data.frame(gene = z$Region, p = as.numeric(p), oracle = as.numeric(orc[z$Region]))
  } else {
    r <- d[d$Group != "Cauchy", ]
    p <- if ("Pvalue" %in% colnames(r)) r$Pvalue else r$Pvalue_Burden
    data.frame(gene = r$Region, p = as.numeric(p), oracle = NA_real_)
  }
}
res <- do.call(rbind, lapply(files, function(f) {
  tg <- parse_tag(f); g <- gene_p(f)
  cbind(g, trait = tg$trait, arm = tg$arm, method = tg$method, collapse = tg$collapse,
        stringsAsFactors = FALSE)
}))
res <- merge(res, truth[, c("gene", "causal", "directional", "NCP_quant", "NCP_binary", "link")], by = "gene")
res$p[!is.finite(res$p)] <- NA
cat("runs:", length(files), " rows:", nrow(res), "\n\n")

# ---- control 1: the null (y_null, y_null_bin; and non-causal genes on the signal traits)
binom_ci <- function(k, n) { if (n == 0) return(c(NA, NA)); b <- binom.test(k, n); b$conf.int }
ecdf_tab <- function(sub, label) {
  for (al in c(0.05, 0.01)) {
    k <- sum(sub$p < al, na.rm = TRUE); n <- sum(!is.na(sub$p))
    ci <- binom_ci(k, n)
    flag <- if (!is.na(ci[1]) && (al < ci[1] || al > ci[2])) " <-- OUTSIDE" else ""
    cat(sprintf("  %-55s alpha=%.2f  obs=%.4f  [%.4f, %.4f]  n=%d%s\n", label, al, k / n, ci[1], ci[2], n, flag))
  }
}
cat("CONTROL 1 -- null calibration: the ECDF at alpha with a binomial 95% CI; expected inside on every line\n")
for (tr in c("y_null", "y_null_bin")) for (cl in unique(res$collapse)) for (m in unique(res$method)) for (a in unique(res$arm)) {
  sub <- res[res$trait == tr & res$collapse == cl & res$method == m & res$arm == a, ]
  if (nrow(sub)) ecdf_tab(sub, paste(tr, "collapse", cl, m, a))
}
cat("\nCONTROL 1b -- the non-causal genes of the signal traits\n")
for (tr in c("y_quant", "y_binary")) for (cl in unique(res$collapse)) for (m in unique(res$method)) for (a in unique(res$arm)) {
  sub <- res[res$trait == tr & res$collapse == cl & res$method == m & res$arm == a & res$causal == 0, ]
  if (nrow(sub)) ecdf_tab(sub, paste(tr, "non-causal, collapse", cl, m, a))
}

# ---- control 2: the positive control (ladder) and the power comparison
cat("\nCONTROL 2 -- Spearman(NCP, -log10 p) over causal genes, expected > 0.5 for every method\n")
cat("   and detections at gene-level alpha = 0.05 / n_genes (Bonferroni) among causal genes\n")
pow <- list()
for (tr in c("y_quant", "y_binary")) {
  ncp <- if (tr == "y_quant") "NCP_quant" else "NCP_binary"
  ng <- length(unique(res$gene[res$trait == tr]))
  for (cl in unique(res$collapse)) for (a in unique(res$arm)) for (m in unique(res$method)) {
    sub <- res[res$trait == tr & res$collapse == cl & res$method == m & res$arm == a & res$causal == 1, ]
    if (!nrow(sub)) next
    rho <- suppressWarnings(cor(sub[[ncp]], -log10(sub$p), method = "spearman", use = "complete.obs"))
    det <- sum(sub$p < 0.05 / ng, na.rm = TRUE)
    det_dir <- sum(sub$p[sub$directional == 1] < 0.05 / ng, na.rm = TRUE)
    orc <- if (all(is.na(sub$oracle))) NA else sum(sub$oracle < 0.05 / ng, na.rm = TRUE)
    by_link <- sapply(c("step", "lof", "graded"), function(L) {
      z <- sub[sub$link == L, ]; if (!nrow(z)) "-" else sprintf("%d/%d", sum(z$p < 0.05 / ng, na.rm = TRUE), nrow(z))
    })
    cat(sprintf("  %-9s c%-2s %-13s %-9s rho=%5.2f  detected %2d/%d (directional %2d/%d)  oracle %s  by link step %s lof %s graded %s\n",
                tr, cl, a, m, rho, det, nrow(sub), det_dir, sum(sub$directional == 1),
                if (is.na(orc)) "-" else as.character(orc), by_link[1], by_link[2], by_link[3]))
    pow[[paste(tr, cl, a, m)]] <- c(det = det, n = nrow(sub))
  }
}

# ---- control 3: monotonicity in the score arm (FlexRV only)
cat("\nCONTROL 3 -- FlexRV detections must be monotone informative >= realistic >= uninformative,\n")
cat("   and on the uninformative arm must not beat the unweighted burden\n")
for (tr in c("y_quant", "y_binary")) for (cl in unique(res$collapse)) {
  g <- function(a, m) { x <- pow[[paste(tr, cl, a, m)]]; if (is.null(x)) NA else x[["det"]] }
  di <- g("informative", "flex192"); dr <- g("realistic", "flex192"); du <- g("uninformative", "flex192")
  dp <- g("plain", "unw")
  ok1 <- !is.na(di) && !is.na(dr) && !is.na(du) && di >= dr && dr >= du
  ok2 <- !is.na(du) && !is.na(dp) && du <= dp
  cat(sprintf("  %-9s c%-2s flex192: informative %s, realistic %s, uninformative %s, unw %s  monotone=%s  uninf<=unw=%s\n",
              tr, cl, di, dr, du, dp, ok1, ok2))
}

# ---- the collapse: 0 vs 10, same trait/arm/method, per-gene log10 p difference
cat("\nCOLLAPSE 0 vs 10 (flex192): median and max |log10 p_c0 - log10 p_c10| over genes\n")
for (tr in unique(res$trait)) for (a in unique(res$arm)) {
  s0 <- res[res$trait == tr & res$arm == a & res$method == "flex192" & res$collapse == "0", ]
  s10 <- res[res$trait == tr & res$arm == a & res$method == "flex192" & res$collapse == "10", ]
  if (!nrow(s0) || !nrow(s10)) next
  mm <- merge(s0, s10, by = "gene")
  dl <- abs(log10(mm$p.x) - log10(mm$p.y))
  cat(sprintf("  %-10s %-13s n=%d median=%.3f max=%.3f  (genes with p<1e-3 under either: %d)\n",
              tr, a, nrow(mm), median(dl, na.rm = TRUE), max(dl, na.rm = TRUE),
              sum(mm$p.x < 1e-3 | mm$p.y < 1e-3, na.rm = TRUE)))
}
