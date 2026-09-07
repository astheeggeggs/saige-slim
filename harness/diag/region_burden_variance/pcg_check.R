# Replicates, in R, the sparse Sigma step 2 hands to C++ (the triplets of the
# symmetric .mtx as read by Matrix, i.e. ONE triangle) and the C++ PCG
# (getPCG1ofSigmaAndGtilde, verbatim: Jacobi, absolute ||r||^2 < 0.02, 100
# iterations), and compares both against the exact symmetric solve on the
# single-pseudo-variant cells of R0.
suppressPackageStartupMessages(library(Matrix))
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
load(file.path(repo, "extdata/output/example_binary_sparseGRM.rda"))
mu <- modglmm$fitted.values
w_mu <- mu * (1 - mu)
tau2 <- modglmm$theta[2]
grm_name <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
grm_file <- file.path(repo, "extdata/output", grm_name)
grm0 <- readMM(grm_file)
ids <- readLines(paste0(grm_file, ".sampleIDs.txt"))
idx <- match(modglmm$sampleID, ids)
cat("readMM class:", class(grm0), " nnz stored:", length(grm0@x), "\n")
sub <- grm0[idx, idx]
cat("after [idx, idx] class:", class(sub), " identity permutation:",
    identical(idx, seq_along(idx)), "\n")
# setSparseSigma_new (R/readInGLMM.R), step by step, relatednessCutoff = 0
x <- sub@x
i <- sub@i
j <- sub@j
drop <- which(x < 0)
if (length(drop) > 0) {
  x <- x[-drop]
  i <- i[-drop]
  j <- j[-drop]
}
x <- x * tau2
dg <- which(i == j)
x[dg] <- x[dg] + 1 / w_mu[i[dg] + 1]
n <- length(mu)
sigma_cpp <- sparseMatrix(i = i + 1, j = j + 1, x = x, dims = c(n, n))
cat("C++ Sigma: nnz", length(x), " symmetric:", isSymmetric(sigma_cpp),
    " off-diagonals below:", sum(i > j), " above:", sum(i < j), "\n")
sigma_sym <- forceSymmetric(sigma_cpp, uplo = if (sum(i > j) > 0) "L" else "U")
pcg <- function(a_mat, b, tol = 0.02, maxit = 100) {
  dinv <- 1 / diag(a_mat)
  x <- numeric(length(b))
  r <- b
  z <- dinv * r
  p <- z
  sumr2 <- sum(r * r)
  it <- 0
  while (sumr2 > tol && it < maxit) {
    it <- it + 1
    ap <- as.numeric(a_mat %*% p)
    alpha <- sum(r * z) / sum(p * ap)
    x <- x + alpha * p
    r1 <- r - alpha * ap
    z1 <- dinv * r1
    beta <- sum(z1 * r1) / sum(z * r)
    p <- z1 + beta * p
    z <- z1
    r <- r1
    sumr2 <- sum(r * r)
  }
  list(x = x, it = it)
}
quad <- function(a_mat, b) {
  s <- pcg(a_mat, b)
  list(v = sum(b * s$x), it = s$it)
}
d <- readRDS(file.path(out_dir, "R0", "analysis.rds"))
vr <- 0.999984072295751 # sparse, category 1, committed file
for (k in which(d$m == 1)) {
  p <- readRDS(file.path(out_dir, "R0", sprintf("phi_%02d.rds", k)))
  g <- p$g.sum
  mac <- d$MAC[k]
  w <- dbeta(mac / (2 * n), 1, 25)
  gu <- g / w # the unweighted residualised pseudo-marker the region path solves on
  tri_u <- quad(sigma_cpp, gu)
  tri_w <- quad(sigma_cpp, g)
  sym_u <- quad(sigma_sym, gu)
  sym_w <- quad(sigma_sym, g)
  cat(sprintf("\n%s  MAC=%d  w=%.4f  ||gu||^2=%.3f\n", d$cell[k], mac, w, sum(gu^2)))
  cat(sprintf("  region C++ var  (B_raw / vr / w^2)          = %.6f\n", d$B_raw[k] / vr / w^2))
  cat(sprintf("  BURDEN C++ var  (var_B / vr / w^2)          = %.6f\n",
              d$var_B_implied[k] / vr / w^2))
  cat(sprintf("  exact solve, symmetric Sigma                = %.6f\n",
              as.numeric(crossprod(gu, solve(sigma_sym, gu)))))
  cat(sprintf("  exact solve, one-triangle Sigma             = %.6f\n",
              as.numeric(crossprod(gu, solve(sigma_cpp, gu)))))
  cat(sprintf("  PCG replica, one-triangle: unweighted %.6f (%d it)  weighted %.6f (%d it)\n",
              tri_u$v, tri_u$it, tri_w$v / w^2, tri_w$it))
  cat(sprintf("  PCG replica, symmetric:    unweighted %.6f (%d it)  weighted %.6f (%d it)\n",
              sym_u$v, sym_u$it, sym_w$v / w^2, sym_w$it))
  cat(sprintf("  naive, no GRM                               = %.6f\n", d$var1_naive[k] / w^2))
}
