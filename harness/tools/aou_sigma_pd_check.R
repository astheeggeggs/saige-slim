#!/usr/bin/env Rscript
# Is the Sigma step 2 solves against positive definite? (LEDGER #183, 2026-10-28)
#
# Step 2 builds Sigma = tau1 I + tau2 K_thinned (quantitative; binary: 1/W on the diagonal),
# where K_thinned is the sparse GRM with every entry below --relatednessCutoff DELETED
# (R/SAIGE_SPATest_Region.R setSparseSigma_new). Deleting entries from a positive-definite
# matrix does not keep it positive definite, and the sparse-GRM variance of a test vector is
# g' Sigma^-1 g from a conjugate-gradient solve that checks only its residual -- so on an
# indefinite Sigma a vector in a negative-eigenvalue direction gets a NEGATIVE variance, which
# the LEDGER #44 guard prints as p 1 and a missing BETA, while a vector near a null direction
# gets a tiny variance and an inflated Z printed as a valid p. Found on an amr Height FlexRV
# cell whose variance walked from 10.7 through 0.39 to minus 2.9 and minus 7.9 across the
# weight sets with the solve converged.
#
# This script asks the question directly, with nothing N x N formed: the thinned GRM is cut
# into connected components (families), each component's Sigma block is eigendecomposed up to
# --maxDenseComponent samples, and the report counts the components whose smallest eigenvalue
# is at or below zero, the most negative one, and how many samples sit in such components
# (people counts below --mask written "<mask": real data). The whole sparse Sigma is also
# handed to CHOLMOD once: a failed Cholesky is the one-bit answer.
#
# Usage (inside the saige-slim image or any R with Matrix):
#   Rscript aou_sigma_pd_check.R --modelFile=M.rda --sparseGRMFile=K.mtx \
#     --sparseGRMSampleIDFile=K.ids --relatednessCutoff=0.125 --reportFile=report.tsv \
#     [--label=amr_Height] [--maxDenseComponent=2000] [--mask=20]
# or, ONCE PER GRM AND CUTOFF with no model (--grmOnly): Sigma = the thinned GRM itself. If
# that is positive semidefinite, tau1 I + tau2 K is positive definite for EVERY model fitted on
# it, so one check per cohort replaces one per phenotype; if it is not, the GRM is the thing
# to rebuild whatever the model. The per-model form remains for a diagnosis of one fit.
suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
})
opts <- list(
  make_option("--modelFile", type = "character"),
  make_option("--sparseGRMFile", type = "character"),
  make_option("--sparseGRMSampleIDFile", type = "character"),
  make_option("--relatednessCutoff", type = "double", default = 0,
              help = "entries below this are deleted, as step 2 does [default %default]"),
  make_option("--reportFile", type = "character"),
  make_option("--label", type = "character", default = "model"),
  make_option("--maxDenseComponent", type = "integer", default = 2000L),
  make_option("--mask", type = "integer", default = 20L),
  make_option("--plantIndefinite", action = "store_true", default = FALSE,
              help = "positive control: overwrite one family's block with an indefinite pattern"),
  make_option("--grmOnly", action = "store_true", default = FALSE,
              help = "no model: is the GRM thinned at the cutoff itself positive semidefinite?"),
  make_option("--modelDir", type = "character", default = NULL,
              help = paste("with --grmOnly: a directory of step-1 .rda models; each is judged",
                           "against the thinned GRM's block eigenvalues (exact for a quantitative",
                           "model, Weyl's bound for a binary one) without further eigen work"))
)
opt <- parse_args(OptionParser(option_list = opts))
for (k in c("sparseGRMFile", "sparseGRMSampleIDFile", "reportFile")) {
  if (is.null(opt[[k]])) stop("--", k, " is required")
}
if (!opt$grmOnly && is.null(opt$modelFile)) stop("--modelFile is required (or --grmOnly)")
report <- new.env()
report$rows <- list()
note <- function(k, v) {
  report$rows[[length(report$rows) + 1]] <-
    data.frame(label = opt$label, key = k, value = as.character(v))
  cat(sprintf("%-34s %s\n", k, as.character(v)))
}
showCount <- function(k) if (k > 0 && k < opt$mask) paste0("<", opt$mask) else as.character(k)

K <- readMM(opt$sparseGRMFile)
kid <- as.character(readLines(opt$sparseGRMSampleIDFile))
kid <- kid[nzchar(kid)]
if (length(kid) != nrow(K)) stop("GRM has ", nrow(K), " rows but ", length(kid), " sample IDs")
if (opt$grmOnly) {
  # Sigma = K thinned at the cutoff: tau1 = 0, tau2 = 1, every GRM sample
  ids <- kid
  n <- length(ids)
  tau1 <- 0
  tau2 <- 1
  binary <- FALSE
  note("mode", paste("grmOnly: Sigma = the GRM thinned at the cutoff (tau1 0, tau2 1);",
                     "PSD here means PD for every model"))
  note("samples", showCount(n))
  diagS <- rep(0, n)
} else {
  e <- new.env()
  load(opt$modelFile, envir = e)
  m <- get("modglmm", envir = e)
  ids <- as.character(m$sampleID)
  n <- length(ids)
  tau1 <- as.numeric(m$theta[1])
  tau2 <- as.numeric(m$theta[2])
  binary <- identical(m$traitType, "binary")
  note("trait_type", if (binary) "binary" else "quantitative")
  note("samples", showCount(n))
  note("tau1", sprintf("%.6g", tau1))
  note("tau2", sprintf("%.6g", tau2))
  # the covariate adjustment term is subtracted from the sparse variance only when step 1 stored
  # this matrix (skipVarianceRatioEstimation & useSparseGRMtoFitNULL); without it var is the raw
  # quadratic form, so a negative value is Sigma's alone
  note("has_Sigma_iXXSigma_iX", !is.null(m$Sigma_iXXSigma_iX))
  mu <- as.numeric(m$fitted.values)   # readInGLMM's mu; mu2 = mu (1 - mu) on a binary trait
  # 1/W on the diagonal, as setSparseSigma_new
  diagS <- if (binary) 1 / (mu * (1 - mu)) else rep(tau1, n)
  note("diag_min", sprintf("%.6g", min(diagS)))
}

if (tau2 <= 0) {
  note("verdict", "tau2 = 0: Sigma is diagonal, positive definite by construction")
} else {
  pos <- match(ids, kid)
  if (anyNA(pos)) stop(sum(is.na(pos)), " model samples are absent from the GRM")
  K <- as(K[pos, pos, drop = FALSE], "TsparseMatrix")
  i <- K@i + 1L
  j <- K@j + 1L
  x <- K@x
  # step 2's thinning: `which(sparseGRM@x < relatednessCutoff)` removed -- the DIAGONAL too if it
  # were below, which it never is at any cutoff in use
  keep <- x >= opt$relatednessCutoff
  i <- i[keep]
  j <- j[keep]
  x <- x[keep]
  off <- i != j
  if (opt$plantIndefinite) {
    # positive control: a star of four -- one hub related at r to three leaves that are
    # unrelated to each other. K's block has smallest eigenvalue 1 - r sqrt(3), so Sigma's block
    # is indefinite once tau2 (r sqrt(3) - 1) exceeds the largest diagonal entry among the four.
    # r is chosen FROM THE MODEL so that holds with a 10 % margin: the first version planted
    # r = 0.9 (K eigenvalue -0.559), which needs tau2 > 1.79 tau1 -- false of universal-saige's
    # e2e model (0.22) and of every bundled binary model (1/W diagonals of 4 and above), where
    # the plant read "positive definite" (2026-10-29). r above 1 is not a kinship, and need not
    # be: the control asks whether an indefinite Sigma is DETECTED, and any indefinite block
    # answers it. The four samples with the smallest diagonal are used, cut out of their
    # families first, so the block is exactly the star; ONE triangle is stored, as readMM's
    # symmetric triplets are, and mirrored with the rest below.
    ord <- order(diagS)
    four <- ord[1:4]
    dmax <- max(diagS[four])
    r <- max(0.9, 1.1 * (1.05 * dmax / tau2 + 1) / sqrt(3))
    # under --grmOnly there is no model to size r from: plant a block no model could sit on
    # (K eigenvalue -172, i.e. a tau2/tau1 bound of 0.006), so --modelDir must refuse every model
    if (opt$grmOnly) r <- 100
    a <- four[1L]
    leaves <- four[2:4]
    drop <- off & ((i %in% four) | (j %in% four))
    i <- c(i[!drop], rep(a, 3))
    j <- c(j[!drop], leaves)
    x <- c(x[!drop], rep(r, 3))
    off <- i != j
    note("planted_indefinite_block",
         sprintf(paste0("a star of 4 at %.4g: K block min eigenvalue 1 - %.4g sqrt(3) = %.4g, ",
                        "Sigma block min %.4g against a diagonal of at most %.4g"),
                 r, r, 1 - r * sqrt(3), dmax + tau2 * (1 - r * sqrt(3)), dmax))
  }
  note("grm_offdiag_pairs_kept", sum(off))
  diagK <- numeric(n)
  diagK[i[!off]] <- x[!off]
  if (any(diagK <= 0)) stop(sum(diagK <= 0), " GRM diagonal entries are missing or non-positive")
  # the scale the eigenvalues are read against: Sigma's diagonal, or K's own under --grmOnly
  diagRef <- if (opt$grmOnly) diagK else diagS
  # Sigma as a general sparse matrix (both triangles)
  S <- sparseMatrix(i = c(i, j[off]), j = c(j, i[off]), x = tau2 * c(x, x[off]), dims = c(n, n))
  S <- S + Diagonal(n, diagS)
  S <- forceSymmetric(S)
  # the one-bit answer: CHOLMOD reports an indefinite matrix as a WARNING with a partial factor
  chol_ok <- tryCatch({
    withCallingHandlers(Cholesky(S, LDL = FALSE, perm = TRUE),
                        warning = function(w) stop(conditionMessage(w)))
    TRUE
  }, error = function(err) FALSE)
  note("sparse_cholesky_succeeds", chol_ok)
  # components of the thinned GRM (label propagation over the kept pairs, as sim_null_pheno.R)
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
  multi <- comp[sizes > 1]
  note("components_multi", length(multi))
  note("largest_component", showCount(max(sizes)))
  note("components_above_dense_cap", sum(sizes > opt$maxDenseComponent))
  n_indef <- 0L
  n_people_indef <- 0L
  min_eig <- Inf
  worst_size <- NA
  n_eval <- 0L
  n_nearnull <- 0L   # smallest eigenvalue below 1e-3 of the diagonal scale: a near-null direction
  block_lo <- list()   # per evaluated block: smallest eigenvalue and samples (for --modelDir)
  for (idx in multi) {
    if (length(idx) > opt$maxDenseComponent) next
    n_eval <- n_eval + 1L
    ev <- eigen(as.matrix(S[idx, idx]), symmetric = TRUE, only.values = TRUE)$values
    lo <- min(ev)
    block_lo[[length(block_lo) + 1L]] <- list(lo = lo, idx = idx)
    if (lo <= 0) {
      n_indef <- n_indef + 1L
      n_people_indef <- n_people_indef + length(idx)
    }
    if (lo > 0 && lo < 1e-3 * mean(diagRef[idx])) n_nearnull <- n_nearnull + 1L
    if (lo < min_eig) {
      min_eig <- lo
      worst_size <- length(idx)
    }
  }
  note("components_evaluated", n_eval)
  note("components_indefinite", n_indef)
  note("samples_in_indefinite_components", showCount(n_people_indef))
  note("components_near_null", n_nearnull)
  note("min_eigenvalue_sigma", sprintf("%.6g", min_eig))
  note("min_eigenvalue_component_size", showCount(worst_size))
  note("min_eigenvalue_over_diag", sprintf("%.6g", min_eig / mean(diagRef)))
  if (opt$grmOnly && (!chol_ok || n_indef > 0)) {
    # The thinned GRM is indefinite, which is NOT by itself a refusal: Sigma = tau1 I + tau2 K
    # (or D + tau2 K, D = 1/W >= 4 on a binary trait) is positive definite whenever
    # tau2 |lambda_min(K)| < tau1 (resp. < min D). eur's thinned GRM at 0.05 carries blocks at
    # -0.145 (runbook 19.10) while every eur model's Sigma is positive definite. So: the bound,
    # and -- with --modelDir -- each model judged against the blocks it was evaluated on.
    if (n_indef > 0 && is.finite(min_eig)) {
      note("tau2_over_tau1_bound_quantitative", sprintf("%.6g", 1 / abs(min_eig)))
      note("tau2_bound_binary", sprintf("%.6g", 4 / abs(min_eig)))
    }
    # chol_ok FALSE with no evaluated block indefinite: the negative direction sits in a block
    # above the dense cap, so no bound can be stated and only a per-model Cholesky decides
    if (!chol_ok && n_indef == 0) {
      verdict <- "GRM indefinite in a block above --maxDenseComponent: bound unknown, run per model"
    } else if (is.null(opt$modelDir)) {
      verdict <- paste("GRM indefinite: positive definite for every model within the tau bound",
                       "(pass --modelDir)")
    } else {
      rdas <- list.files(opt$modelDir, pattern = "\\.rda$", full.names = TRUE)
      if (length(rdas) == 0) stop("--modelDir ", opt$modelDir, " holds no .rda")
      n_bad <- 0L
      for (f in rdas) {
        em <- new.env()
        load(f, envir = em)
        mm <- get("modglmm", envir = em)
        t1 <- as.numeric(mm$theta[1])
        t2 <- as.numeric(mm$theta[2])
        bin <- identical(mm$traitType, "binary")
        mids <- as.character(mm$sampleID)
        # the model's samples in GRM order; a sample the model drops leaves its block, which
        # can only raise the block's smallest eigenvalue (interlacing), so the full-block
        # eigenvalue is a conservative bound for the model's own block
        margin <- Inf
        if (t2 > 0) {
          if (bin) {
            mu <- as.numeric(mm$fitted.values)
            d <- 1 / (mu * (1 - mu))
            dmin_all <- min(d)
            dpos <- match(mids, kid)
            dmin_by_sample <- rep(NA_real_, n)
            dmin_by_sample[dpos[!is.na(dpos)]] <- d[!is.na(dpos)]
          }
          for (b in block_lo) {
            if (b$lo >= 0) next
            dmin <- if (bin) min(dmin_by_sample[b$idx], na.rm = TRUE) else t1
            if (!is.finite(dmin)) dmin <- if (bin) dmin_all else t1
            margin <- min(margin, dmin + t2 * b$lo)
          }
        }
        ok <- margin > 0
        if (!ok) n_bad <- n_bad + 1L
        note(paste0("model_", sub("\\.rda$", "", basename(f))),
             sprintf("%s tau1 %.4g tau2 %.4g min(Sigma block eigenvalue bound) %.4g: %s",
                     if (bin) "binary" else "quantitative", t1, t2, margin,
                     if (ok) "positive definite" else "NOT positive definite"))
      }
      note("models_checked", length(rdas))
      note("models_not_positive_definite", n_bad)
      verdict <- if (n_bad == 0) {
        "GRM indefinite: Sigma positive definite for every model in --modelDir"
      } else {
        "Sigma is NOT positive definite for some model in --modelDir"
      }
    }
  } else {
    verdict <- if (!chol_ok || n_indef > 0) {
      "Sigma is NOT positive definite"
    } else if (n_nearnull > 0) {
      "positive definite with near-null directions"
    } else {
      "positive definite"
    }
  }
  note("verdict", verdict)
}
out <- do.call(rbind, report$rows)
write.table(out, opt$reportFile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("done\n")
