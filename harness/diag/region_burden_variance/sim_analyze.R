# sim_analyze.R: (2) null calibration of the three p-value columns under each
# candidate fix, (3) planted-effect recovery of BETA_Burden / SE_Burden, and
# the burden-variance identity across every replicate -- all offline from the
# per-replicate dumps and rep_stats, by the namespace functions in fixes.R.
suppressPackageStartupMessages(library(SAIGE))
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
source(file.path(repo, "harness", "diag", "region_burden_variance", "fixes.R"))
fx <- diag_fixes()
cell_fun <- diag_cell # bound at top level so the closures below see a defined name
options(width = 220, digits = 5)
specs <- c("103", "103A", "103B")
distinct <- c(1, 2, 3, 4, 7, 8, 9, 10) # the missense;lof;synonymous rows duplicate missense;lof
rep_dirs <- function(label) {
  sort(list.files(file.path(out_dir, "sim", label), pattern = "^rep_", full.names = TRUE))
}
cell_table <- function(label) {
  done <- Filter(function(rd) file.exists(file.path(rd, "phi_12.rds")), rep_dirs(label))
  cache <- file.path(out_dir, "sim", label, "cells.rds")
  if (file.exists(cache)) {
    cached <- readRDS(cache)
    if (identical(cached$n, length(done))) return(cached$tab)
  }
  rows <- list()
  for (rd in done) {
    if (!file.exists(file.path(rd, "phi_12.rds"))) next
    for (k in distinct) {
      d <- readRDS(file.path(rd, sprintf("phi_%02d.rds", k)))
      s <- readRDS(file.path(rd, sprintf("skat_%02d.rds", k)))
      o0 <- cell_fun(d, fx, "", s$r.corr)
      if (!identical(o0$Pvalue_SKATO, s$out$Pvalue_SKATO) ||
            !identical(o0$BETA_Burden, s$out$BETA_Burden)) {
        stop("offline pipeline does not reproduce the dumped output: ", rd, " cell ", k)
      }
      for (v in specs) {
        o <- cell_fun(d, fx, v, s$r.corr)
        rows[[length(rows) + 1]] <- data.frame(
          rep = basename(rd), cell = k, spec = v,
          p_skato = o$Pvalue_SKATO, p_burden = o$Pvalue_Burden, p_skat = o$Pvalue_SKAT,
          beta = o$BETA_Burden, se = o$SE_Burden, B_adj = o$B_adj, S = sum(d$Score),
          naive = sum(d$mu * (1 - d$mu) * d$g.sum^2), B_raw = sum(d$Phi)
        )
      }
    }
  }
  tab <- do.call(rbind, rows)
  saveRDS(list(n = length(done), tab = tab), cache)
  tab
}
stats_table <- function(label) {
  rows <- lapply(rep_dirs(label), function(rd) {
    f <- file.path(rd, "rep_stats.rds")
    if (!file.exists(f)) return(NULL)
    as.data.frame(readRDS(f))
  })
  do.call(rbind, rows)
}

# ---- (2) null calibration -------------------------------------------------
nul <- cell_table("null")
cat(sprintf("\n== null arm: %d replicates x %d distinct cells ==\n",
            length(unique(nul$rep)), length(distinct)))
fired <- abs(nul$B_adj[nul$spec == "103"] / nul$naive[nul$spec == "103"] - 1) < 1e-6
cat("fraction of cells where today's rescale fired (B_adj(103) == naive, r < 1):",
    round(mean(fired), 3), "\n")
# The single-pseudo-variant cells (MAC 5-15) carry exact-resampling and
# saddlepoint p-values whose null distribution is discrete, so a KS test
# rejects uniformity whatever the calibration; rejection rates and the mean of
# the chi-square-inverted p are the readable metrics, and the two cell types
# are reported separately. Within a replicate the cells are correlated, so the
# binomial SE quoted below is optimistic.
calib <- function(p) {
  p <- p[!is.na(p)]
  c(n = length(p), `rej@0.01` = mean(p < 0.01), `rej@0.05` = mean(p < 0.05),
    `rej@0.10` = mean(p < 0.10), mean_chi2 = mean(qchisq(p, 1, lower.tail = FALSE)))
}
multi <- c(2, 4, 8, 10)
for (col in c("p_burden", "p_skat", "p_skato")) {
  groups <- list(list("multi-marker cells", multi),
                 list("single-pseudo-variant cells", setdiff(distinct, multi)))
  for (grp in groups) {
    cat(sprintf("\n-- %s, %s --\n", col, grp[[1]]))
    sel <- nul$cell %in% grp[[2]]
    m <- t(sapply(specs, function(v) calib(nul[[col]][sel & nul$spec == v])))
    print(round(m, 4))
  }
}
n_p <- sum(nul$spec == "103" & nul$cell %in% multi)
cat(sprintf("binomial SE of a rejection rate at 0.05 with %d p-values: %.4f; of mean_chi2: %.3f\n",
            n_p, sqrt(0.05 * 0.95 / n_p), sqrt(2 / n_p)))

# ---- burden-variance identity from rep_stats (pure R, refit model) --------
# per arm: `null` at the model's tau2 (0.33), `nullx` at tau2 = 4 where the
# sparse-GRM and naive variances separate far enough for the test to have power
for (label in c("null", "nullx")) {
  st <- stats_table(label)
  if (is.null(st)) next
  cat(sprintf("\n== burden score variance, arm %s (n = %d), planted cell's burden ==\n",
              label, nrow(st)))
  v_emp <- var(st$S)
  cat(sprintf("Var_emp(S) = %.2f  (relative SE %.3f)\n", v_emp, sqrt(2 / (nrow(st) - 1))))
  for (cn in c("V_naive", "V_sparse", "V_sparse_vr", "V_burdenmode")) {
    cat(sprintf("  mean %-13s = %8.2f   Var_emp / mean = %.4f\n", cn, mean(st[[cn]]),
                v_emp / mean(st[[cn]])))
  }
  cat(sprintf("  mean tau2 %.3f   mean null ratio cat 2 %.4f   mean sparse ratio %.5f\n",
              mean(st$tau2), mean(st$vr_null2), mean(st$vr_sparse2)))
  z2 <- function(v) mean(st$S^2 / st[[v]])
  cat(sprintf("  S^2/V mean: naive %.4f  sparse %.4f  burdenmode %.4f  (1 if right; SE %.3f)\n",
              z2("V_naive"), z2("V_sparse"), z2("V_burdenmode"), sqrt(2 / nrow(st))))
}

# ---- (3) planted effect ---------------------------------------------------
pl <- cell_table("planted")
pl4 <- pl[pl$cell == 4, ]
beta_b <- 0.05
cat(sprintf("\n== planted arm: beta_b = %.3f on GENE1 missense;lof 0.01, %d replicates ==\n",
            beta_b, length(unique(pl4$rep))))
summ <- function(beta, se) {
  c(mean_BETA = mean(beta), `mean/beta_b` = mean(beta) / beta_b, sd_BETA = sd(beta),
    mean_SE = mean(se), `sd/mean_SE` = sd(beta) / mean(se),
    coverage95 = mean(abs(beta - beta_b) <= 1.96 * se),
    `SE of mean/beta_b` = sd(beta) / sqrt(length(beta)) / beta_b)
}
m <- t(sapply(specs, function(v) summ(pl4$beta[pl4$spec == v], pl4$se[pl4$spec == v])))
stp <- stats_table("planted")
bm <- stp$S / stp$V_burdenmode
m <- rbind(m, `BURDEN mode (offline, score SE)` = summ(bm, 1 / sqrt(stp$V_burdenmode)),
           `oracle S/V_sparse_exact` = summ(stp$S / stp$V_sparse, 1 / sqrt(stp$V_sparse)))
print(round(m, 4))
s_dump <- pl4$S[pl4$spec == "103"]
cat(sprintf("S from the region dump vs S from the refit model in R: ratio %.6f to %.6f\n",
            min(s_dump / stp$S), max(s_dump / stp$S)))
cat("fraction of planted replicates where the burden SPA fired (|z_naive| >= 2):",
    round(mean(abs(stp$S) / sqrt(stp$V_naive) >= 2), 3), "\n")
