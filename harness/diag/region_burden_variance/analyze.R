# analyze.R <skato-run> <burden-run>
# Reconciles, per cell, the region path's Phi before and after the R-side
# adjustment (SPA_ER_kernel_related_Phiadj_fast_new) against the BURDEN
# path's pseudo-variant test and an independent sparse-GRM variance solved
# from the model and the .mtx with Matrix. Reads the dumps drv.R wrote.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(TRUE)
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
skato_dir <- file.path(out_dir, args[1])
burden_dir <- file.path(out_dir, args[2])

# --- independent sparse Sigma: tau2 * GRM + diag(1 / (mu (1 - mu))) ---------
load(file.path(repo, "extdata/output/example_binary_sparseGRM.rda"))
mu <- modglmm$fitted.values
w_mu <- mu * (1 - mu)
tau2 <- modglmm$theta[2]
grm_name <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
grm_file <- file.path(repo, "extdata/output", grm_name)
grm <- readMM(grm_file)
ids <- readLines(paste0(grm_file, ".sampleIDs.txt"))
idx <- match(modglmm$sampleID, ids)
stopifnot(!anyNA(idx))
grm <- as(grm, "generalMatrix")[idx, idx]
sigma <- tau2 * grm + Diagonal(x = 1 / w_mu)
v_sparse <- function(g) as.numeric(crossprod(g, solve(sigma, g)))
v_null <- function(g) sum(w_mu * g^2)

read_cells <- function(path) {
  d <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
  d[d$Group != "Cauchy", ]
}
rt <- read_cells(file.path(skato_dir, "out.txt"))
bt <- read_cells(file.path(burden_dir, "out.txt"))
stopifnot(nrow(rt) == nrow(bt), all(rt$MAC == bt$MAC))

one_cell <- function(k) {
  a <- readRDS(file.path(skato_dir, sprintf("adj_%02d.rds", k)))
  p <- readRDS(file.path(skato_dir, sprintf("phi_%02d.rds", k)))
  s <- readRDS(file.path(skato_dir, sprintf("skat_%02d.rds", k)))
  phi <- as.matrix(a$Phi)
  score <- a$Score
  s_sum <- sum(score)
  # the function's own arithmetic, recomputed from its inputs
  var_s <- score^2 / qchisq(a$p.new, df = 1, lower.tail = FALSE, log.p = TRUE)
  d_fac <- var_s / diag(phi)
  var_q <- sum(t(t(phi * sqrt(d_fac))) * sqrt(d_fac))
  q_b <- qchisq(a$p_value_burden, df = 1, lower.tail = FALSE, log.p = TRUE)
  var_q2 <- s_sum^2 / q_b
  r_fac <- min(var_q / var_q2, 1)
  phi_adj <- as.matrix(a$out$val)
  g <- p$g.sum
  m1 <- sum(p$mu * g)
  naive <- v_null(g)
  data.frame(
    cell = paste(rt$Region[k], rt$Group[k], rt$max_MAF[k]), m = length(score),
    MAC = rt$MAC[k], S = s_sum, qsum_m1 = p$q.sum - m1,
    B_raw = sum(phi), A_raw = sum(diag(phi)), Dmax = max(d_fac),
    z_naive = (p$q.sum - m1) / sqrt(naive), VarQ = var_q, VarQ_2 = var_q2, r = r_fac,
    B_adj = sum(phi_adj), A_adj = sum(diag(phi_adj)),
    var1_naive = naive, Vsp_exact = v_sparse(g),
    BETA_R = s$out$BETA_Burden, BETA_R_fix = s_sum / sum(phi_adj), p_R = s$out$Pvalue_Burden,
    BETA_B = bt$BETA_Burden[k], var_B_implied = s_sum / bt$BETA_Burden[k],
    p_B = bt$Pvalue_Burden[k]
  )
}
d <- do.call(rbind, lapply(seq_len(nrow(rt)), one_cell))
options(width = 250, digits = 7)
cat("\n== identities ==\n")
print(data.frame(
  cell = d$cell, m = d$m, MAC = d$MAC,
  `S==qsum-m1` = d$qsum_m1 / d$S,
  `BETA_R==S/trPhiAdj` = d$BETA_R / (d$S / d$A_adj),
  Dmax = d$Dmax, z_naive = round(d$z_naive, 3), r = d$r,
  `B_adj/var1_naive` = d$B_adj / d$var1_naive,
  `B_adj/VarQ_2` = d$B_adj / d$VarQ_2,
  `B_raw/var1_naive` = d$B_raw / d$var1_naive,
  `B_raw/Vsp_exact` = d$B_raw / d$Vsp_exact,
  `Vsp/Vnull` = d$Vsp_exact / d$var1_naive,
  check.names = FALSE
), row.names = FALSE)
cat("\n== cross-mode ==\n")
print(data.frame(
  cell = d$cell, m = d$m, MAC = d$MAC,
  `BETA_Rfix/BETA_B` = d$BETA_R_fix / d$BETA_B,
  `varB_impl/B_raw` = d$var_B_implied / d$B_raw,
  `varB_impl/var1_naive` = d$var_B_implied / d$var1_naive,
  `varB_impl/Vsp_exact` = d$var_B_implied / d$Vsp_exact,
  p_R = d$p_R, p_B = d$p_B,
  p_naive = pchisq(d$S^2 / d$var1_naive, 1, lower.tail = FALSE),
  check.names = FALSE
), row.names = FALSE)
saveRDS(d, file.path(skato_dir, "analysis.rds"))
