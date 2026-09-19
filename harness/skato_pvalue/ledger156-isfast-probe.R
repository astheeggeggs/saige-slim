suppressMessages(library(SKAT)); ns <- asNamespace("SKAT")
hits <- 0L; orig <- get("Get_Lambda_Approx", ns)
assignInNamespace("Get_Lambda_Approx", function(...) { hits <<- hits + 1L; orig(...) }, ns = ns)
set.seed(11); n <- 3000; q <- 2200          # ABOVE the p1 >= 2000 guard
G <- matrix(rbinom(n * q, 2, 0.01), n, q)
Phi <- crossprod(scale(G, scale = FALSE)) / n
Score <- as.numeric(rnorm(q) * 0.25)         # keep p off the underflow floor
grid <- c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 0.999)
probe <- function(label, rc) {
  for (fast in c(FALSE, TRUE)) {
    hits <<- 0L
    t0 <- proc.time()[[3]]
    p <- SKAT:::Met_SKAT_Get_Pvalue(Score, Phi, rc, "optimal.adj", NULL, isFast = fast)$p.value
    el <- proc.time()[[3]] - t0
    cat(sprintf("%-38s isFast=%-5s  Approx=%d  p=%.6e  %6.3f s\n", label, fast, hits, p, el))
  }
}
cat("kernel dimension q =", q, " (guard is p1 >= 2000)\n\n")
probe("r.corr=0  -> plain path, Get_PValue", 0)
probe("8-rho grid -> the real SKAT-O path", grid)
