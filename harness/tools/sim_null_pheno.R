#!/usr/bin/env Rscript
# Simulate a NULL phenotype from a fitted SAIGE null model: the model's own fixed part, a
# random effect drawn from the SPARSE GRM at the model's tau2, residual noise at tau1 (binary:
# a Bernoulli draw on the logistic scale). Step 1 is then refitted on the simulated column and
# step 2 run on the REAL genotypes and the REAL group files, so every gene-level p-value is a
# null p-value on real genotype structure with no annotation confound -- the control the
# synonymous mask only approximates (synonymous sites sit under different selection, in
# different genes, with no AlphaMissense score).
#
# Why not naive. (1) Nothing N x N is ever formed: the GRM is read sparse, thinned at the
# cutoff step 1 used, cut into connected components, and each component is drawn on its own.
# (2) A thinned GRM is NOT positive definite in general -- zeroing the entries below the
# cutoff can leave a family block indefinite -- and the step-1 fit tolerates that only because
# Sigma adds the residual term. A sampler does not: a component up to --maxDenseComponent
# samples is eigendecomposed and its NEGATIVE eigenvalues are clipped to zero (the nearest
# PSD matrix in Frobenius norm), and every clip is counted and its lost mass reported; a
# larger component takes a sparse Cholesky with the smallest jitter that makes it succeed,
# reported likewise. (3) tau2 = 0 (most AoU fits) draws no random effect at all.
#
# The fixed part is what the model stored: under isCovariateOffset (the default) the covariate
# GLM's linear predictor `offset` plus the mixed model's intercept; otherwise X %*% coefficients.
# Checked on the bundled models: eta - fixed is the BLUP (sd 0.12 binary, 0.58 quantitative)
# and var(y - fixed) = 0.6795 against tau1 + tau2 = 0.6807 on quant_sparseGRM.
#
# Usage:
#   Rscript sim_null_pheno.R --modelFile=M.rda --sparseGRMFile=K.mtx --sparseGRMSampleIDFile=K.ids \
#     --phenoFile=pheno.txt --phenoCol=height --sampleIDCol=person_id --seed=1 \
#     --outPhenoFile=pheno_simnull.txt --reportFile=report.tsv [--relatednessCutoff=0.05]
#     [--tau2=<override>] [--maxDenseComponent=2000] [--mask=20]
# The output is the input phenotype file with phenoCol REPLACED by the simulated value for the
# model's samples and NA for every other row, so the production step-1 launcher runs on it
# unchanged. Counts of people in the report below --mask are written "<mask" (real data).
# Self-test: sim_null_pheno_selftest.sh.
suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(data.table)
})
opts <- list(
  make_option("--modelFile", type = "character"),
  make_option("--sparseGRMFile", type = "character"),
  make_option("--sparseGRMSampleIDFile", type = "character"),
  make_option("--phenoFile", type = "character"),
  make_option("--phenoCol", type = "character"), make_option("--sampleIDCol", type = "character"),
  make_option("--outPhenoFile", type = "character"),
  make_option("--reportFile", type = "character"),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--relatednessCutoff", type = "double", default = 0,
              help = "zero GRM entries below this, as step 1 did (amr: 0.125) [default %default]"),
  make_option("--tau2", type = "double", default = NA_real_, help = "override the model's tau2"),
  make_option("--maxDenseComponent", type = "integer", default = 2000L,
              help = paste("components up to this size are eigen-clipped;",
                           "larger take a jittered sparse Cholesky")),
  make_option("--mask", type = "integer", default = 20L,
              help = "mask people counts below this (real data)")
)
opt <- parse_args(OptionParser(option_list = opts))
for (k in c("modelFile", "sparseGRMFile", "sparseGRMSampleIDFile", "phenoFile", "phenoCol",
            "sampleIDCol", "outPhenoFile", "reportFile")) {
  if (is.null(opt[[k]])) stop("--", k, " is required")
}
showCount <- function(k) if (k > 0 && k < opt$mask) paste0("<", opt$mask) else as.character(k)
report <- new.env()
note <- function(key, value) {
  assign(key, as.character(value), envir = report)
  cat(sprintf("%-28s %s\n", key, value))
}
# ---- the model ------------------------------------------------------------------------------
e <- new.env()
load(opt$modelFile, envir = e)
m <- e$modglmm
if (is.null(m)) stop("no modglmm in ", opt$modelFile)
ids <- as.character(m$sampleID)
n <- length(ids)
trait <- as.character(m$traitType)
if (!trait %in% c("binary", "quantitative")) stop("traitType ", trait)
fixed <- if (isTRUE(m$isCovariateOffset)) {
  if (length(m$offset) != n) stop("isCovariateOffset but no offset of length N in the model")
  as.numeric(m$offset) + as.numeric(m$coefficients)[1]
} else {
  as.vector(m$X %*% m$coefficients)
}
if (length(fixed) != n || any(!is.finite(fixed))) stop("fixed part is not a finite length-N vector")
tau1 <- as.numeric(m$theta)[1]
tau2 <- as.numeric(m$theta)[2]
if (!is.na(opt$tau2)) tau2 <- opt$tau2
note("model", basename(opt$modelFile))
note("traitType", trait)
note("n_samples", n)
note("isCovariateOffset", isTRUE(m$isCovariateOffset))
note("tau1", tau1)
note("tau2", tau2)
note("seed", opt$seed)
if (trait == "binary") {
  note("observed_prevalence", sprintf("%.6f", mean(as.numeric(m$y))))
  note("fixed_part_prevalence", sprintf("%.6f", mean(plogis(fixed))))
} else {
  note("observed_var_resid", sprintf("%.6f", var(as.numeric(m$y) - fixed)))
  note("tau1_plus_tau2", sprintf("%.6f", tau1 + tau2))
}
# ---- the GRM, restricted to the model's samples, thinned as step 1 thinned it ----------------
set.seed(opt$seed)
b <- numeric(n)
diagK <- rep(1, n)
if (tau2 > 0) {
  K <- readMM(opt$sparseGRMFile)
  kid <- as.character(readLines(opt$sparseGRMSampleIDFile))
  kid <- kid[nzchar(kid)]
  if (length(kid) != nrow(K)) stop("GRM has ", nrow(K), " rows but ", length(kid), " sample IDs")
  pos <- match(ids, kid)
  if (anyNA(pos)) stop(sum(is.na(pos)), " model samples are absent from the GRM")
  K <- as(K[pos, pos, drop = FALSE], "TsparseMatrix")
  i <- K@i + 1L
  j <- K@j + 1L
  x <- K@x
  # step 1's thinning, `which(sparseGRM@x < cutoff)` zeroed: at the default 0 that also drops
  # every NEGATIVE entry, which is what the fit saw
  keep <- i == j | x >= opt$relatednessCutoff
  i <- i[keep]
  j <- j[keep]
  x <- x[keep]
  diagK <- numeric(n)
  diagK[i[i == j]] <- x[i == j]
  if (any(diagK <= 0)) stop(sum(diagK <= 0), " GRM diagonal entries are missing or non-positive")
  off <- i != j
  note("grm_offdiag_pairs_kept", sum(off))
  note("grm_mean_diag", sprintf("%.6f", mean(diagK)))
  # connected components by label propagation over the kept off-diagonal pairs (a sparse
  # symmetric matrix at a relatedness cutoff is families: the loop runs as many rounds as the
  # widest component's diameter)
  lab <- seq_len(n)
  if (any(off)) {
    ei <- i[off]
    ej <- j[off]
    repeat {
      a <- pmin(lab[ei], lab[ej])
      mi <- tapply(a, ei, min)
      mj <- tapply(a, ej, min)
      new <- lab
      ki <- as.integer(names(mi))
      new[ki] <- pmin(new[ki], as.integer(mi))
      kj <- as.integer(names(mj))
      new[kj] <- pmin(new[kj], as.integer(mj))
      if (identical(new, lab)) break
      lab <- new
    }
  }
  comp <- split(seq_len(n), lab)
  sizes <- lengths(comp)
  note("components", length(comp))
  note("components_singleton", showCount(sum(sizes == 1)))
  note("component_max_size", showCount(max(sizes)))
  note("components_above_dense_cap", sum(sizes > opt$maxDenseComponent))
  # one matrix for the multi-sample components, drawn block by block
  Koff <- sparseMatrix(i = c(i[off], j[off]), j = c(j[off], i[off]), x = c(x[off], x[off]),
                       dims = c(n, n))
  n_clipped <- 0L
  mass_clipped <- 0
  mass_total <- 0
  min_eig <- Inf
  n_jitter <- 0L
  max_jitter <- 0
  for (idx in comp) {
    s <- length(idx)
    if (s == 1L) {
      b[idx] <- sqrt(tau2 * diagK[idx]) * rnorm(1L)
      next
    }
    Kc <- as.matrix(Koff[idx, idx, drop = FALSE])
    diag(Kc) <- diagK[idx]
    z <- rnorm(s)
    if (s <= opt$maxDenseComponent) {
      ev <- eigen(Kc, symmetric = TRUE)
      min_eig <- min(min_eig, ev$values[s])
      mass_total <- mass_total + sum(abs(ev$values))
      if (ev$values[s] < 0) {
        n_clipped <- n_clipped + 1L
        mass_clipped <- mass_clipped + sum(abs(ev$values[ev$values < 0]))
        ev$values[ev$values < 0] <- 0
      }
      b[idx] <- sqrt(tau2) * as.vector(ev$vectors %*% (sqrt(ev$values) * z))
    } else {
      # too large to eigendecompose: sparse Cholesky, jittered until it succeeds
      Ks <- as(Koff[idx, idx, drop = FALSE], "CsparseMatrix")
      diag(Ks) <- diagK[idx]
      Ks <- forceSymmetric(Ks)
      jit <- 0
      ch <- NULL
      # CHOLMOD reports an indefinite matrix as a WARNING and returns a partial factor, so a
      # warning is a failure here, like an error
      for (f in c(0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.3, 0.5)) {
        jit <- f * mean(diagK[idx])
        ch <- tryCatch(Cholesky(Ks + Diagonal(s, jit), LDL = FALSE, perm = TRUE),
                       error = function(err) NULL, warning = function(w) NULL)
        if (!is.null(ch)) break
      }
      if (is.null(ch)) stop("component of ", s, " samples is not PD even at jitter 0.5 x diag")
      if (jit > 0) {
        n_jitter <- n_jitter + 1L
        max_jitter <- max(max_jitter, jit)
      }
      ex <- expand(ch)   # Ks + jit I = P' L L' P
      if (s <= 5000L) {   # the convention, checked where it is cheap (every self-test component)
        rec <- crossprod(ex$P, tcrossprod(ex$L)) %*% ex$P
        if (max(abs(rec - (Ks + Diagonal(s, jit)))) > 1e-8 * mean(diagK[idx])) {
          stop("Cholesky reconstruction failed")
        }
      }
      b[idx] <- sqrt(tau2) * as.vector(crossprod(ex$P, ex$L %*% z))
    }
  }
  note("components_eigen_clipped", n_clipped)
  note("min_eigenvalue", sprintf("%.6g", if (is.finite(min_eig)) min_eig else NA))
  note("clipped_mass_fraction",
       sprintf("%.3e", if (mass_total > 0) mass_clipped / mass_total else 0))
  note("components_jittered", n_jitter)
  note("max_jitter", sprintf("%.3e", max_jitter))
  # the pair diagnostic the self-test reads: over the kept pairs, mean(b_i b_j) / mean(tau2 k_ij)
  if (any(off)) {
    obs <- mean(b[i[off]] * b[j[off]])
    expct <- tau2 * mean(x[off])
    note("pair_cov_observed", sprintf("%.6f", obs))
    note("pair_cov_expected", sprintf("%.6f", expct))
  }
} else {
  note("components", "not formed (tau2 = 0: no random effect)")
}
note("var_b", sprintf("%.6f", var(b)))
# ---- the draw ------------------------------------------------------------------------------
eta <- fixed + b
y <- if (trait == "binary") rbinom(n, 1L, plogis(eta)) else eta + rnorm(n, sd = sqrt(tau1))
if (trait == "binary") {
  note("simulated_cases", showCount(sum(y == 1L)))
  note("simulated_prevalence", sprintf("%.6f", mean(y)))
} else {
  note("simulated_var_resid", sprintf("%.6f", var(y - fixed)))
}
# ---- the phenotype file, phenoCol replaced --------------------------------------------------
hdr <- readLines(opt$phenoFile, n = 1L)
sep <- if (grepl("\t", hdr, fixed = TRUE)) "\t" else if (grepl(",", hdr, fixed = TRUE)) "," else " "
ph <- fread(opt$phenoFile, sep = sep, colClasses = "character", na.strings = c("NA", ""))
if (!opt$sampleIDCol %in% names(ph)) stop("no column ", opt$sampleIDCol, " in ", opt$phenoFile)
if (!opt$phenoCol %in% names(ph)) stop("no column ", opt$phenoCol, " in ", opt$phenoFile)
pid <- as.character(ph[[opt$sampleIDCol]])
hit <- match(pid, ids)
if (sum(!is.na(hit)) != n) {
  stop("phenotype file holds ", sum(!is.na(hit)), " of the model's ", n, " samples")
}
val <- rep(NA_character_, nrow(ph))
ys <- y[hit[!is.na(hit)]]
val[!is.na(hit)] <- if (trait == "binary") as.character(ys) else sprintf("%.10g", ys)
set(ph, j = opt$phenoCol, value = val)
fwrite(ph, opt$outPhenoFile, sep = sep, na = "NA", quote = FALSE)
note("rows_written", nrow(ph))
note("rows_with_value", n)
keys <- ls(report)
tab <- data.table(field = keys, value = unlist(mget(keys, envir = report)))
fwrite(tab, opt$reportFile, sep = "\t")
cat("wrote", opt$outPhenoFile, "and", opt$reportFile, "\n")
