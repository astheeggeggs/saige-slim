# Real per-rho SKAT-O spectra for cells with q = 4..10, exactly as
# skatoScreenKernel builds them: R_rho = (1-rho)I + rho 11', L = chol(pivot),
# lambda = Get_Lambda(L (Phi/2) L').
suppressPackageStartupMessages(library(SKAT))
ns <- asNamespace("SKAT"); GL <- get("Get_Lambda", ns)
W <- commandArgs(TRUE)[1]
fs <- Sys.glob(file.path(W, "runs", "*.cells.rds"))
fs <- fs[!grepl("PILOT", fs)][1:6]          # a sample of arms is enough for a condition-number census
rows <- list(); spec <- list()
for (f in fs) {
  ce <- readRDS(f)
  for (i in seq_along(ce)) {
    x <- ce[[i]]; q <- length(x$Score)
    if (q < 4 || q > 10 || length(x$r.corr) <= 1) next
    Phi2 <- x$Phi / 2; p.m <- nrow(x$Phi)
    r.all <- x$r.corr; r.all[r.all >= 0.999] <- 0.999
    for (rc in r.all) {
      RM <- diag(rep(1 - rc, p.m)) + matrix(rc, p.m, p.m)
      L <- chol(RM, pivot = TRUE)
      lam <- GL(L %*% (Phi2 %*% t(L)))
      lam <- lam[lam > 0]
      if (length(lam) < 2) next
      rows[[length(rows) + 1L]] <- data.frame(run = basename(f), cell = i, q = q, rho = rc,
                                              n_lam = length(lam), cond = max(lam) / min(lam))
      spec[[length(spec) + 1L]] <- lam
    }
  }
}
df <- do.call(rbind, rows)
cat("per-rho spectra collected:", nrow(df), "from", length(fs), "arms\n")
cat("n_lam (terms in the quadratic form):\n"); print(table(df$n_lam))
cat("\ncondition number quantiles:\n"); print(signif(quantile(df$cond, c(0,.25,.5,.75,.9,.99,1)), 4))
cat("\ncondition by n_lam (median):\n")
print(signif(tapply(df$cond, df$n_lam, median), 4))
saveRDS(list(meta = df, spectra = spec), file.path(dirname(commandArgs(TRUE)[2]), "spectra.rds"))
cat("\nsaved\n")
