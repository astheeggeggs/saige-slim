# Characterization of LEDGER #106: step 2's sparse Sigma is ONE triangle, and
# the PCG that solves against it stops on an absolute residual after one or
# two iterations.
#
# setSparseSigma_new reads @i/@j/@x of the dsTMatrix that Matrix::readMM
# returns for the `symmetric` .mtx, which carries the diagonal and the entries
# BELOW it only, so the arma::sp_mat built from those triplets is lower
# triangular. getPCG1ofSigmaAndGtilde then iterates on that matrix until
# ||r||^2 < 0.02 with no ||b|| normalisation: an unweighted rare variant
# (||b||^2 ~ MAC) stops after ONE iteration.
#
# Pinned here from golden P27 (plink_marker_sparse.txt, every marker on the
# sparse arm): the printed `var` of rs4 equals the one-triangle one-iteration
# replica and NOT the exact symmetric solve. The replica is built from the
# same model, the same GRM through setSparseSigma_new itself, and rs4's
# genotypes read straight from the .bed. Tstat is reproduced first as the
# control that the model and the genotype are the right ones.
#
# LEDGER #106 FIXED: setSparseSigma_new expands the subset to a general
# triplet matrix (both triangles) and the PCG stops on a relative residual
# (||r||^2 <= 1e-8 ||b||^2). The assertions marked INVERTED pinned the defect;
# they now pin the fix: golden P27's var is the exact symmetric solve, and the
# one-triangle one-iteration replica is what it is NOT.

sigma_case <- function() {
  skip_if_no_goldens()
  model_rda <- file.path(saige_goldens, "step1", "binary_sparseGRM_vr.rda")
  vr_file <- file.path(saige_goldens, "step1", "binary_sparseGRM_vr.varianceRatio.txt")
  golden <- file.path(saige_goldens, "step2", "plink_marker_sparse.txt")
  skip_if_not(file.exists(model_rda) && file.exists(golden), "P27 inputs not present")
  e <- new.env()
  load(model_rda, envir = e)
  modglmm <- e$modglmm
  grm_prefix <- file.path(saige_extdata, "output",
                          "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  mu <- as.numeric(modglmm$fitted.values)
  W <- mu * (1 - mu)
  sig <- SAIGE:::setSparseSigma_new(grm_prefix, paste0(grm_prefix, ".sampleIDs.txt"), 0,
                                     modglmm$sampleID, modglmm$theta, W, "binary")
  # rs4: bim column 5 (SAIGE's Allele2 under --AlleleOrder=alt-first) counted.
  # The .bed holds 10,000 samples; the model's 1,000 are matched by IID, as
  # step 2 does.
  bed_prefix <- file.path(saige_extdata, "input", "genotype_100markers")
  fam_iid <- read.table(paste0(bed_prefix, ".fam"), stringsAsFactors = FALSE)$V2
  g <- read_bed_marker(bed_prefix, "rs4")[match(modglmm$sampleID, fam_iid)]
  stopifnot(!anyNA(g))
  X <- modglmm$X
  XtWX_inv <- solve(crossprod(X, X * W))
  gtilde <- as.numeric(g - X %*% (XtWX_inv %*% crossprod(X, W * g)))
  vr <- read.table(vr_file)
  vr_sparse1 <- vr$V1[vr$V2 == "sparse" & vr$V3 == 1]
  gold <- read.table(golden, header = TRUE, sep = "\t")
  list(sig = sig, mu = mu, gtilde = gtilde, vr = vr_sparse1, gold = gold[gold$MarkerID == "rs4", ],
       y = as.numeric(modglmm$y), g = g)
}

# Minimal SNP-major .bed reader for one marker: 00 = two copies of bim column
# 5, 10 = one, 11 = none, 01 = missing (none in this fixture).
read_bed_marker <- function(prefix, marker) {
  bim <- read.table(paste0(prefix, ".bim"), stringsAsFactors = FALSE)
  fam <- read.table(paste0(prefix, ".fam"), stringsAsFactors = FALSE)
  n <- nrow(fam)
  k <- which(bim$V2 == marker)
  stopifnot(length(k) == 1)
  bytes_per <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", 3)
  stopifnot(identical(as.integer(magic), c(108L, 27L, 1L)))
  seek(con, 3 + (k - 1) * bytes_per)
  raw_bytes <- as.integer(readBin(con, "raw", bytes_per))
  codes <- integer(bytes_per * 4)
  for (s in 0:3) codes[seq(s + 1, by = 4, length.out = bytes_per)] <- bitwAnd(bitwShiftR(raw_bytes, 2 * s), 3L)
  codes <- codes[seq_len(n)]
  stopifnot(!any(codes == 1L))
  c(2, NA, 1, 0)[codes + 1]
}

# The C++ PCG, verbatim (Jacobi, ||r||^2 against an absolute tol, at most
# maxit iterations), over the matrix `a`.
pcg_replica <- function(a, b, tol = 0.02, maxit = 100) {
  dinv <- 1 / Matrix::diag(a)
  x <- numeric(length(b)); r <- b; z <- dinv * r; p <- z; sumr2 <- sum(r * r); it <- 0
  while (sumr2 > tol && it < maxit) {
    it <- it + 1
    ap <- as.numeric(a %*% p); alpha <- sum(r * z) / sum(p * ap); x <- x + alpha * p
    r1 <- r - alpha * ap; z1 <- dinv * r1; beta <- sum(z1 * r1) / sum(z * r)
    p <- z1 + beta * p; z <- z1; r <- r1; sumr2 <- sum(r * r)
  }
  list(x = x, it = it)
}

test_that("the control holds: the replica's Tstat is golden P27's", {
  tc <- sigma_case()
  expect_equal(sum(tc$g), tc$gold$AC_Allele2)
  expect_equal(sum(tc$gtilde * (tc$y - tc$mu)), tc$gold$Tstat, tolerance = 1e-5)
})

test_that("setSparseSigma_new hands C++ both triangles of the symmetric GRM", {
  tc <- sigma_case()
  i <- tc$sig$locations[1, ]; j <- tc$sig$locations[2, ]
  below <- sum(i > j); above <- sum(i < j)
  expect_equal(sum(i == j), tc$sig$nSubj)
  # INVERTED (#106 fixed): both triangles, mirror images of each other.
  expect_gt(below, 0)
  expect_identical(above, below)
  lo <- data.frame(i = i[i > j], j = j[i > j], x = tc$sig$values[i > j])
  up <- data.frame(i = j[i < j], j = i[i < j], x = tc$sig$values[i < j])
  lo <- lo[order(lo$i, lo$j), ]; up <- up[order(up$i, up$j), ]
  expect_equal(lo$i, up$i); expect_equal(lo$j, up$j); expect_equal(lo$x, up$x)
})

test_that("golden P27's var is the exact symmetric solve, not the one-triangle one-iteration PCG", {
  tc <- sigma_case()
  n <- tc$sig$nSubj
  a_sym <- Matrix::sparseMatrix(i = tc$sig$locations[1, ] + 1, j = tc$sig$locations[2, ] + 1,
                                x = tc$sig$values, dims = c(n, n))
  expect_true(Matrix::isSymmetric(a_sym))
  a_tri <- Matrix::tril(a_sym)   # what C++ solved against before the fix
  rep <- pcg_replica(a_tri, tc$gtilde)
  var_replica <- tc$vr * sum(tc$gtilde * rep$x)
  var_exact <- tc$vr * as.numeric(crossprod(tc$gtilde, Matrix::solve(a_sym, tc$gtilde)))
  # The old rule still stops after one iteration on ||b||^2 ~ MAC = 20 ...
  expect_equal(rep$it, 1)
  # INVERTED (#106 fixed). The golden is the exact variance under the model's
  # Sigma, to its printed precision -- which is the tolerance: `var` is written
  # with 6 significant digits (rs4 prints 1.62837 for 1.6283745), so the golden
  # cannot express better than ~5e-6 relative and a tighter bound would be
  # measuring the writer's %g, not the solver.
  expect_equal(tc$gold$var, var_exact, tolerance = 1e-5)
  # INVERTED. ... and it is NOT the one-triangle one-iteration value any more:
  expect_gt(abs(tc$gold$var / var_replica - 1), 1e-4)
  # The code's own rule, replicated: relative 1e-8 on the symmetric matrix
  # converges in a few iterations and lands on the golden.
  new <- pcg_replica(a_sym, tc$gtilde, tol = 1e-8 * sum(tc$gtilde^2))
  expect_gt(new$it, 1)
  expect_equal(tc$vr * sum(tc$gtilde * new$x), tc$gold$var, tolerance = 1e-5)
})

test_that("a relative tolerance makes the variance independent of the vector's scale", {
  # Under the old absolute rule the same genotype times its Beta(1,25) weight
  # (~23 at MAF 0.01) took a second iteration and landed on a DIFFERENT
  # variance for one marker. INVERTED (#106 fixed): with the relative rule the
  # two scalings agree to ~1e-8; the old rule's dependence is kept below as
  # the record of what was wrong.
  tc <- sigma_case()
  n <- tc$sig$nSubj
  a_sym <- Matrix::sparseMatrix(i = tc$sig$locations[1, ] + 1, j = tc$sig$locations[2, ] + 1,
                                x = tc$sig$values, dims = c(n, n))
  w <- dbeta(tc$gold$AF_Allele2, 1, 25)
  rel <- function(b) pcg_replica(a_sym, b, tol = 1e-8 * sum(b^2))
  v_u <- sum(tc$gtilde * rel(tc$gtilde)$x)
  v_w <- sum(w * tc$gtilde * rel(w * tc$gtilde)$x) / w^2
  expect_equal(v_w / v_u, 1, tolerance = 1e-7)
  a_tri <- Matrix::tril(a_sym)
  old_u <- pcg_replica(a_tri, tc$gtilde); old_w <- pcg_replica(a_tri, w * tc$gtilde)
  expect_gt(old_w$it, old_u$it)
  expect_gt(abs(sum(w * tc$gtilde * old_w$x) / w^2 / sum(tc$gtilde * old_u$x) - 1), 1e-4)
})
