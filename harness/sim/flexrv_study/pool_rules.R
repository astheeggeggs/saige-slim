# LEDGER #142: what the pooled p_FlexRV would be under alternative pooling rules,
# measured on the sim50k study's 192-set outputs (harness/work/flexrv_study),
# without re-running SAIGE. Per region the 192 per-set burden p-values are
# re-pooled under:
#   all      the shipped rule: equal-weight CCT over every set row (192)
#   dedup    drop sets whose (p, beta, se) row is byte-identical to an earlier
#            set of the SAME score transform (LEDGER #142 candidate b)
#   apriori  drop the MAC rungs at or above the ingest bound 2*N*maxMAF before
#            pooling (computable from the flags and N alone; here bound = 100
#            so MAC100/200/500 go, 16 x 9 = 144 sets remain)
#   axis     weighted CCT: each score transform gets total weight 1/16 and its
#            DISTINCT MAF-axis sets share it equally (duplicates of unw count
#            once), so the score axis is not re-weighted by the duplicate count
# Reported: null-arm ECDF at 0.05 / 0.01 with a binomial CI, power-arm
# detections at Bonferroni 0.05/n_genes, and the per-rule rank agreement with
# the shipped rule. Usage: Rscript pool_rules.R [results_dir] [sim_dir]
suppressWarnings(library(SAIGE))
args <- commandArgs(TRUE)
R <- if (length(args) >= 1) args[1] else "/Users/dpalmer/Repositories/SAIGE/harness/work/flexrv_study/results"
D <- if (length(args) >= 2) args[2] else "/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k"
N <- 50000; maxMAF <- 0.001; bound <- 2 * N * maxMAF
rd <- function(p) utils::read.table(p, header = TRUE, sep = "\t", colClasses = "character",
                                    quote = "", comment.char = "", na.strings = character(0))
truth <- rd(file.path(D, "truth_genes.tsv")); truth$causal <- as.integer(truth$causal)
cct <- function(p, w = NULL) { p <- as.numeric(p); k <- !is.na(p); if (!any(k)) return(NA_real_)
  SAIGE:::CCT(pmax(p[k], .Machine$double.xmin), if (is.null(w)) NULL else w[k]) }
files <- list.files(R, pattern = "flex192_c(0|10)_r1\\.txt$", full.names = TRUE)
rung_of <- function(w) sub("^SIM_[^_]+_", "", w)          # MAF transform name
score_of <- function(w) sub("^SIM_([^_]+)_.*$", "\\1", w)  # score transform name
pool_region <- function(r) {
  sc <- score_of(r$Weight); mf <- rung_of(r$Weight)
  key <- paste(r$Pvalue_Burden, r$BETA_Burden, r$SE_Burden)
  dup <- ave(key, sc, FUN = function(k) duplicated(k)) == "TRUE"
  macrung <- suppressWarnings(as.numeric(sub("^MAC", "", mf)))
  apri <- !is.na(macrung) & macrung >= bound
  # axis weights: per score transform, 1/16 shared equally over its distinct sets
  w <- numeric(nrow(r)); for (s in unique(sc)) { i <- which(sc == s & !dup); w[i] <- (1 / 16) / length(i) }
  c(all = cct(r$Pvalue_Burden), dedup = cct(r$Pvalue_Burden[!dup]),
    apriori = cct(r$Pvalue_Burden[!apri]), axis = cct(r$Pvalue_Burden, w),
    n_dup = sum(dup), n_apri = sum(apri))
}
out <- list()
for (f in files) {
  d <- rd(f); z <- d[d$Group == "Cauchy", ]; r <- d[d$Group != "Cauchy", ]
  tag <- sub("_r1\\.txt$", "", basename(f))
  per <- do.call(rbind, lapply(split(r, r$Region), pool_region))
  per <- data.frame(gene = rownames(per), per, shipped = as.numeric(z$Pvalue_Burden[match(rownames(per), z$Region)]),
                    tag = tag, stringsAsFactors = FALSE)
  out[[tag]] <- per
}
res <- do.call(rbind, out); res <- merge(res, truth[, c("gene", "causal")], by = "gene")
# control: the recomputed 'all' must reproduce the shipped pooled row
cat("CONTROL: recomputed equal-weight CCT vs the file's Cauchy row -- max |dlog10 p| =",
    signif(max(abs(log10(res$all) - log10(res$shipped)), na.rm = TRUE), 3), "\n")
cat("duplicate sets per region: value-identical median", median(res$n_dup), " a-priori rungs >= bound:", unique(res$n_apri), "\n\n")
binom_ci <- function(k, n) if (n == 0) c(NA, NA) else binom.test(k, n)$conf.int
rules <- c("all", "dedup", "apriori", "axis")
cat("NULL ARMS -- ECDF at alpha per rule (binomial 95% CI)\n")
for (tag in grep("^y_null", names(out), value = TRUE)) for (al in c(0.05, 0.01)) {
  s <- res[res$tag == tag, ]
  line <- sprintf("  %-42s alpha=%.2f", tag, al)
  for (rl in rules) { k <- sum(s[[rl]] < al, na.rm = TRUE); n <- sum(!is.na(s[[rl]])); ci <- binom_ci(k, n)
    line <- paste0(line, sprintf("  %s=%.4f", rl, k / n)) }
  ci <- binom_ci(round(al * nrow(s)), nrow(s))
  cat(line, sprintf("  [CI %.4f, %.4f] n=%d\n", ci[1], ci[2], nrow(s)))
}
cat("\nPOWER ARMS -- detections at Bonferroni 0.05/n among causal genes, and Spearman rho vs 'all'\n")
for (tag in grep("^y_(quant|binary)_", names(out), value = TRUE)) {
  s <- res[res$tag == tag, ]; ng <- nrow(s); cz <- s[s$causal == 1, ]
  line <- sprintf("  %-42s causal %3d/%3d", tag, nrow(cz), ng)
  for (rl in rules) line <- paste0(line, sprintf("  %s=%2d", rl, sum(cz[[rl]] < 0.05 / ng, na.rm = TRUE)))
  rho <- sapply(rules[-1], function(rl) cor(s$all, s[[rl]], method = "spearman", use = "complete.obs"))
  cat(line, sprintf("  rho(dedup,apriori,axis)=%.3f,%.3f,%.3f\n", rho[1], rho[2], rho[3]))
}
cat("\nPER-GENE RATIO p_rule / p_all, all arms pooled: quantiles\n")
for (rl in rules[-1]) { q <- quantile(res[[rl]] / res$all, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)
  cat(sprintf("  %-8s %s\n", rl, paste(signif(q, 3), collapse = "  "))) }
