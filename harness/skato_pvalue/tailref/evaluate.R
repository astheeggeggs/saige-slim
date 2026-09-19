# Evaluate the SHIPPED saddlepoint survival and SKAT's Davies at q-values that
# sit at target tail probabilities, on REAL small-q spectra.  Writes a CSV for
# the independent Ruben reference (Python) to score.
suppressPackageStartupMessages({library(SKAT); library(SAIGE)})
SP <- commandArgs(TRUE)[1]
d <- readRDS(file.path(SP, "spectra.rds")); meta <- d$meta; spec <- d$spectra
SPA <- SAIGE:::skatoSaddlepointSurvival
DAV <- get("SKAT_davies", asNamespace("SKAT"))

band <- cut(meta$cond, c(0, 10, 100, 1000, Inf), labels = c("le10", "10-100", "100-1e3", "gt1e3"))
set.seed(1)
pick <- unlist(lapply(split(seq_len(nrow(meta)), list(meta$n_lam, band), drop = TRUE),
                      function(ix) if (length(ix) <= 4) ix else sample(ix, 4)))
cat("spectra selected:", length(pick), "\n")
TARGETS <- c(1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
rows <- list()
for (j in pick) {
  lam <- spec[[j]]
  f <- function(lq) log(max(SPA(lam, exp(lq)), 1e-300)) - log(tg)
  for (tg in TARGETS) {
    qq <- try(exp(uniroot(f, c(log(1e-8), log(1e9)), tol = 1e-10)$root), silent = TRUE)
    if (inherits(qq, "try-error")) next
    ps <- SPA(lam, qq)
    dv <- DAV(qq, lam, lim = 1e4, acc = 1e-6)
    rows[[length(rows) + 1L]] <- data.frame(
      id = length(rows) + 1L, n_lam = length(lam), cond = max(lam) / min(lam),
      target = tg, q = qq, p_spa = ps, p_davies = dv$Qq, davies_ifault = dv$ifault,
      lambda = paste(sprintf("%.17g", lam), collapse = "|"), stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(SP, "eval.csv"), row.names = FALSE)
cat("rows:", nrow(out), "\n")
cat("davies ifault != 0:", sum(out$davies_ifault != 0), "\n")
cat("davies returning <= 0 (broken):", sum(out$p_davies <= 0), "\n")
