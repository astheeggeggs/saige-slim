# Region (per-marker) vs `--r.corr=1` (collapsed) burden variance against the EXACT model variance (2026-09-07)

Status: measurement only; no code, golden or ledger change. Companion to
`diag-2026-09-07-region-burden-variance.md` (reads its R0 dumps under
`harness/work/diag-burden/R0`). Question answered: after #103/#104 are fixed,
which of the two burden constructions is the better approximation of the model's
own burden variance `g_b' P_s g_b` (sparse-GRM model, `useSparseGRMtoFitNULL = TRUE`,
so the exact target is a dense symmetric solve at n = 1000)?

Controls (all pass): rebuilt weighted residualised burden vector equals the dumped
`g.sum` to 1e-14 on all 12 cells; the replica of the region path (per-row arm and
ratio, cross block = later chunk's P1 x earlier chunk's P2, one-triangle Jacobi
PCG with absolute tol 0.02) reproduces the C++ `1'Phi 1` to 5 digits; the burden
replica reproduces `var_B` to 5 digits; `Tstat[rs4]` matches `.singleAssoc`.

| cell | q | MAC pseudo | MAC cell | region construct, exact solves / exact | region as coded / exact | burden construct, exact / exact | burden as coded / exact |
|---|---|---|---|---|---|---|---|
| GENE1 lof 0.001 | 1 | 5 | 5 | 1.00000 | 0.99917 | 1.00000 | 1.00017 |
| GENE1 lof 0.01 | 4 | 45 | 98 | 1.00081 | 1.00520 | 0.99388 | 0.99388 |
| GENE1 missense;lof 0.001 | 1 | 10 | 10 | 1.00000 | 1.00410 | 1.00000 | 1.00516 |
| GENE1 missense;lof 0.01 | 5 | 54 | 119 | 1.00045 | 1.00637 | 1.00099 | 1.00098 |
| GENE2 lof 0.001 | 1 | 5 | 5 | 1.00000 | 1.00271 | 1.00000 | 1.00344 |
| GENE2 lof 0.01 | 2 | 9 | 22 | 0.99999 | 1.00884 | 0.99703 | 0.99703 |
| GENE2 missense;lof 0.001 | 1 | 15 | 15 | 1.00000 | 1.00429 | 1.00000 | 1.00481 |
| GENE2 missense;lof 0.01 | 4 | 41 | 77 | 0.99765 | 1.00486 | 0.99737 | 0.99737 |

(The synonymous cells duplicate the missense;lof rows.) Reading: with exact solves
the per-marker (paper) construction is within 0.25 % of exact, its residual being
the pseudo-marker's common-arm diagonal (`r_null x no-GRM variance`) when its MAC
exceeds 20.5. The collapsed construction is exact on the sparse arm and within
0.6 % on the common arm, where it applies the random-marker null ratio (0.9421)
to a vector whose true ratio is 0.941-0.948. As CODED the region path is the worse
of the two on multi-marker cells (+0.5 to +0.9 %) and that excess is entirely the
one-triangle PCG truncation (Evidence 4 of the companion note), not the
construction. The exact/naive ratios of the burden vectors (0.941-0.951) also
confirm the VR file's null rows. Both approximations' errors scale with how much
the GRM matters; here `1/W ~ 11` dominates `tau psi_s ~ 0.16`.

Also measured (scratch `proj_check.R`, same session): the projection term the
region `VarMat` omits relative to the paper's `P_s` is 0.0012-0.0020 % of the
variance on the sparse arm here and identically 0 on the no-GRM arm
(`X'W g~ = 0`, measured 1e-14).

## Scripts (verbatim; run from the repo root with `SC=<dir holding g100.raw>`; `g100.raw` is `plink2 --bfile extdata/input/genotype_100markers --export A`; SAIGE's counted allele is bim column 5, so every column is flipped)

### exact_check.R
```r
suppressPackageStartupMessages(library(Matrix))
repo <- getwd(); out_dir <- file.path(repo, "harness/work/diag-burden")
load(file.path(repo, "extdata/output/example_binary_sparseGRM.rda"))
mu <- as.numeric(modglmm$fitted.values); w_mu <- mu * (1 - mu); tau2 <- modglmm$theta[2]; X <- modglmm$X
grm_file <- file.path(repo, "extdata/output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
grm0 <- readMM(grm_file); ids <- readLines(paste0(grm_file, ".sampleIDs.txt")); idx <- match(modglmm$sampleID, ids)
sub <- grm0[idx, idx]; x <- sub@x; i <- sub@i; j <- sub@j
drop <- which(x < 0); if (length(drop)) { x <- x[-drop]; i <- i[-drop]; j <- j[-drop] }
x <- x * tau2; dg <- which(i == j); x[dg] <- x[dg] + 1 / w_mu[i[dg] + 1]; n <- length(mu)
sigma_sym <- as.matrix(forceSymmetric(sparseMatrix(i = i + 1, j = j + 1, x = x, dims = c(n, n)), uplo = "L"))
SX <- solve(sigma_sym, X); XSX_inv <- solve(crossprod(X, SX))
exact_var <- function(g) { Sg <- solve(sigma_sym, g); a <- sum(g * Sg); b <- as.numeric(t(g) %*% SX %*% XSX_inv %*% t(SX) %*% g); a - b }
d <- readRDS(file.path(out_dir, "R0", "analysis.rds"))
cat(sprintf("%-36s %3s %4s | %9s %9s %9s | %8s %8s %8s | %8s\n", "cell", "m", "MAC", "exact", "region", "burden", "reg/ex", "bur/ex", "naive/ex", "ex/naive"))
for (k in seq_along(d$cell)) {
  p <- readRDS(file.path(out_dir, "R0", sprintf("phi_%02d.rds", k))); g <- p$g.sum
  ex <- exact_var(g)
  cat(sprintf("%-36s %3d %4d | %9.5f %9.5f %9.5f | %8.5f %8.5f %8.5f | %8.5f\n", d$cell[k], d$m[k], d$MAC[k],
      ex, d$B_raw[k], d$var_B_implied[k], d$B_raw[k] / ex, d$var_B_implied[k] / ex, d$var1_naive[k] / ex, ex / d$var1_naive[k]))
}
```

### decompose.R
```r
suppressPackageStartupMessages(library(Matrix))
repo <- getwd(); out_dir <- file.path(repo, "harness/work/diag-burden"); SC <- Sys.getenv("SC")
load(file.path(repo, "extdata/output/example_binary_sparseGRM.rda"))
mu <- as.numeric(modglmm$fitted.values); w_mu <- mu * (1 - mu); tau2 <- modglmm$theta[2]; X <- modglmm$X
y <- as.numeric(modglmm$y); res <- y - mu; n <- length(mu)
grm_file <- file.path(repo, "extdata/output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
grm0 <- readMM(grm_file); ids <- readLines(paste0(grm_file, ".sampleIDs.txt")); idx <- match(modglmm$sampleID, ids)
sub <- grm0[idx, idx]; x <- sub@x; i <- sub@i; j <- sub@j
drop <- which(x < 0); if (length(drop)) { x <- x[-drop]; i <- i[-drop]; j <- j[-drop] }
x <- x * tau2; dg <- which(i == j); x[dg] <- x[dg] + 1 / w_mu[i[dg] + 1]
sigma_cpp <- sparseMatrix(i = i + 1, j = j + 1, x = x, dims = c(n, n))
sigma_sym <- as.matrix(forceSymmetric(sigma_cpp, uplo = "L"))
SX <- solve(sigma_sym, X); XSX_inv <- solve(crossprod(X, SX))
pcg <- function(a_mat, b, tol = 0.02, maxit = 100) {
  dinv <- 1 / diag(a_mat); xx <- numeric(length(b)); r <- b; z <- dinv * r; p <- z; sumr2 <- sum(r * r); it <- 0
  while (sumr2 > tol && it < maxit) {
    it <- it + 1; ap <- as.numeric(a_mat %*% p); alpha <- sum(r * z) / sum(p * ap); xx <- xx + alpha * p
    r1 <- r - alpha * ap; z1 <- dinv * r1; beta <- sum(z1 * r1) / sum(z * r); p <- z1 + beta * p; z <- z1; r <- r1; sumr2 <- sum(r * r)
  }
  xx
}
# genotypes
raw <- read.table(file.path(SC, "g100.raw"), header = TRUE, check.names = FALSE)
G <- as.matrix(raw[, 7:ncol(raw)]); counted <- sub(".*_", "", colnames(G)); colnames(G) <- sub("_[ACGT]$", "", colnames(G))
bim <- read.table(file.path(repo, "extdata/input/genotype_100markers.bim")); a2 <- setNames(bim$V6, bim$V2)
# SAIGE Allele2 (the counted allele) is bim column 5, per the singleAssoc rows (rs4: Allele1 A, Allele2 C, AC 20)
a1 <- setNames(bim$V5, bim$V2); flip <- counted != a1[colnames(G)]; G[, flip] <- 2 - G[, flip]; cat(sprintf("flipped %d columns to count bim column 5\n", sum(flip)))
G <- G[match(modglmm$sampleID, raw$IID), ]
AC <- colSums(G); AF <- AC / (2 * n); MAF <- pmin(AF, 1 - AF); MAC <- pmin(AC, 2 * n - AC)
stopifnot(all(AF <= 0.5))
cat(sprintf("control AC[rs4] = %d (expect 20)\n", AC["rs4"]))
XtWX_inv <- solve(crossprod(X, X * w_mu))
proj <- function(g) as.numeric(g - X %*% (XtWX_inv %*% crossprod(X, w_mu * g)))
Gt <- apply(G, 2, proj)
cat(sprintf("control Tstat[rs4] = %.6f (singleAssoc -0.249913)\n", sum(Gt[, "rs4"] * res)))
anno <- c(rep("lof", 10), rep("missense", 20), rep("lof", 20), rep("missense", 30), rep("lof", 20))
gene <- c(rep("GENE1", 50), rep("GENE2", 50)); names(anno) <- names(gene) <- colnames(G)
vr_sparse <- c(0.999984072295751, 0.999979229384575); vr_null <- c(0.947537894959354, 0.942146422285391)
cat_of <- function(m) ifelse(m <= 20.5, 1, 2)
d <- readRDS(file.path(out_dir, "R0", "analysis.rds"))
hdr <- sprintf("%-26s %2s %4s %4s | %8s %8s %8s | %8s %8s | %7s %7s %7s", "cell", "m", "MACc", "MACb",
               "regEx/ex", "regPCG/ex", "Braw/ex", "burEx/ex", "varB/ex", "gsum", "reg/Braw", "bur/varB")
cat(hdr, "\n")
for (k in seq_along(d$cell)) {
  parts <- strsplit(d$cell[k], " ")[[1]]; gn <- parts[1]; an <- strsplit(parts[2], ";")[[1]]; mm <- as.numeric(parts[3])
  sel <- gene == gn & anno %in% an & MAF <= mm & MAC >= 0.5
  rare <- which(sel & MAC > 10); ur <- which(sel & MAC <= 10)
  gc <- apply(G[, ur, drop = FALSE], 1, max); AFc <- mean(gc) / 2; MAFc <- min(AFc, 1 - AFc); MACc <- MAFc * 2 * n
  V <- cbind(Gt[, rare, drop = FALSE], UR = proj(gc))            # unweighted residualised vectors, pseudo LAST (its own chunk)
  macs <- c(MAC[rare], MACc); mafs <- c(MAF[rare], MAFc); u <- dbeta(mafs, 1, 25); q <- ncol(V)
  gsum <- as.numeric(V %*% u)
  p <- readRDS(file.path(out_dir, "R0", sprintf("phi_%02d.rds", k)))
  gsum_err <- max(abs(gsum - p$g.sum))
  # exact model variance of the burden vector
  S_gsum <- solve(sigma_sym, gsum); ex <- sum(gsum * S_gsum) - as.numeric(t(gsum) %*% SX %*% XSX_inv %*% t(SX) %*% gsum)
  # per-row arm and ratio as coded
  sparse_arm <- macs <= 20.5; r <- ifelse(sparse_arm, vr_sparse[cat_of(macs)], vr_null[cat_of(macs)])
  # P2 columns: exact and PCG replica; common arm -> W g
  P2ex <- solve(sigma_sym, V); P2pcg <- apply(V, 2, function(b) pcg(sigma_cpp, b)); P2W <- V * w_mu
  build <- function(P2s) {
    M <- crossprod(V, P2s)                                # M_ij = v_i' P2_j (sparse arm for column j)
    for (jj in which(!sparse_arm)) M[jj, jj] <- sum(V[, jj] * P2W[, jj])   # common-arm diagonal: no GRM
    # cross block: later chunk's P1 (pseudo row) times earlier chunk's P2 (rare, sparse arm), then mirrored
    if (q > 1) { cr <- as.numeric(crossprod(V[, q], P2s[, 1:(q - 1), drop = FALSE])); M[q, 1:(q - 1)] <- cr; M[1:(q - 1), q] <- cr }
    Phi <- M * (sqrt(r) %o% sqrt(r)) * (u %o% u); sum(Phi)
  }
  reg_ex <- build(P2ex); reg_pcg <- build(P2pcg)
  MACb <- d$MAC[k]; rb <- if (MACb <= 20.5) vr_sparse[cat_of(MACb)] else vr_null[cat_of(MACb)]
  bur_ex <- rb * (if (MACb <= 20.5) sum(gsum * S_gsum) else sum(gsum^2 * w_mu))
  bur_pcg <- rb * (if (MACb <= 20.5) sum(gsum * pcg(sigma_cpp, gsum)) else sum(gsum^2 * w_mu))
  cat(sprintf("%-26s %2d %4.0f %4d | %8.5f %8.5f %8.5f | %8.5f %8.5f | %7.0e %7.5f %7.5f\n", d$cell[k], q, MACc, MACb,
      reg_ex / ex, reg_pcg / ex, d$B_raw[k] / ex, bur_ex / ex, d$var_B_implied[k] / ex, gsum_err, reg_pcg / d$B_raw[k], bur_pcg / d$var_B_implied[k]))
}
```
