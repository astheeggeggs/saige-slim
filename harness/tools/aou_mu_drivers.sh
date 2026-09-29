#!/bin/bash
# WHICH covariate puts a null model's fitted values in the tail? (LEDGER #174/#175)
#
# aou_mu_decades.sh says how many samples have a tiny flip probability. This
# says why, from the model alone: which covariates entered the fit, how many
# cases it has, and which covariate carries the tail samples' linear predictor
# away from everyone else's.
#
# HOW. A saved model's design matrix X is the QR-transformed one: orthogonal
# columns of mean square 1, still named after the covariates in the order they
# were given (R/SAIGE_fitGLMM_fast.R Covariate_Transform; a constant or
# collinear column was dropped before it and is simply absent). So the linear
# predictor eta = qlogis(fitted value) decomposes exactly and SEQUENTIALLY:
#   eta = sum_j gamma_j Q_j + remainder,   gamma_j = mean(Q_j * eta)
# where column j is covariate j with every earlier covariate projected out (the
# first covariate is itself, centred and scaled), and the remainder is what no
# covariate explains (the random effect). Per covariate the table gives
#   share_of_variance   gamma_j^2 / var(eta)
#   tail_minus_rest     mean of gamma_j Q_j over the TAIL samples (flip
#                       probability at or below TAIL) minus the mean over the
#                       rest, in log-odds. The column sums to the whole gap.
# A categorical covariate that separates a group shows as ONE column carrying
# the gap; a steep continuous effect shows as that covariate and its square.
#
# Runs INSIDE a dsub job on the SAIGE image (R only). Nothing individual is
# written: every number is a sum over the cohort or over the tail group, the
# tail group's size is written "<MASK" when it is from 1 to MASK-1 and its
# rows are then withheld, and the case count is masked the same way. Env:
#   MODELS   a directory, or a file pattern such as /path/eas_FemInf.rda
#   OUT_DIR  local dir dsub delocalizes (--output-recursive)
#   TAG      label for the output file, default models
#   TAIL     the flip probability that defines the tail, default 1e-6
#   MASK     default 20
#   RSCRIPT  default Rscript
# Output: <OUT_DIR>/mu_drivers.<TAG>.tsv.gz, and a summary per model on stdout.
set -euo pipefail
: "${MODELS:?}"; : "${OUT_DIR:?}"; TAG="${TAG:-models}"; TAIL="${TAIL:-1e-6}"; MASK="${MASK:-20}"; RSCRIPT="${RSCRIPT:-Rscript}"
mkdir -p "${OUT_DIR}"
if [[ -d ${MODELS} ]]; then root="${MODELS}"; pat='*.rda'; else root="$(dirname "${MODELS}")"; pat="$(basename "${MODELS}")"; fi
[[ -d ${root} ]] || { echo "ERROR: no such directory: ${root} (MODELS=${MODELS})"; exit 1; }
find "${root}" -type f -name "${pat}" -name '*.rda' ! -name '*_chr[0-9]*.rda' | sort > models.list
nfound=$(awk 'END{print NR}' models.list)
echo "models found: ${nfound} (directory ${root}, pattern ${pat})"
(( nfound > 0 )) || { echo "ERROR: no model file matched"; exit 1; }
cat > mu_drivers.R <<'RS'
args <- commandArgs(trailingOnly = TRUE)
files <- readLines(args[1]); out <- args[2]; tail_at <- as.numeric(args[3]); mask <- as.integer(args[4])
show <- function(k) ifelse(k > 0 & k < mask, paste0("<", mask), as.character(k))
num <- function(x) formatC(x, format = "g", digits = 4)
rows <- list()
for (f in files) {
  e <- new.env()
  ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
  if (!ok || is.null(e$modglmm)) { cat("skipped (no modglmm):", basename(f), "\n"); next }
  m <- e$modglmm
  if (!identical(as.character(m$traitType), "binary")) next
  mu <- as.numeric(m$fitted.values); y <- as.numeric(m$y); X <- as.matrix(m$X); n <- length(mu)
  if (nrow(X) != n) { cat("skipped (X has", nrow(X), "rows for", n, "fitted values):", basename(f), "\n"); next }
  # a design saved without column names still decomposes; its columns are numbered
  if (is.null(colnames(X))) colnames(X) <- c("column1_intercept", paste0("column", seq_len(ncol(X))[-1]))
  flip <- ifelse(y == 1, 1 - mu, mu)
  tail <- flip <= tail_at
  nt <- sum(tail)
  eta <- stats::qlogis(pmin(pmax(mu, 1e-300), 1 - 1e-16))
  # the orthogonal design, recomputed so that a model saved with the original
  # X decomposes the same way; each column signed like the covariate it is named after
  Q <- qr.Q(qr(X)) * sqrt(n)
  for (j in seq_len(ncol(Q))) if (sum(Q[, j] * X[, j]) < 0) Q[, j] <- -Q[, j]
  gam <- as.vector(crossprod(Q, eta)) / n
  part <- sweep(Q, 2, gam, `*`)
  rest <- eta - rowSums(part)
  v <- mean((eta - mean(eta))^2)
  share <- c(gam^2, mean((rest - mean(rest))^2)) / v
  share[1] <- 0                                   # the intercept column is a constant
  name <- c(colnames(X), "(no covariate: random effect)")
  gap <- if (nt >= mask && nt <= n - mask) {
    c(colMeans(part[tail, , drop = FALSE]) - colMeans(part[!tail, , drop = FALSE]), mean(rest[tail]) - mean(rest[!tail]))
  } else rep(NA_real_, length(name))
  gap[1] <- if (is.na(gap[1])) NA_real_ else 0
  cat("=====", basename(f), "\n")
  cat("  N", n, " cases", show(sum(y == 1)), " tail (flip probability at or below", format(tail_at), "):", show(nt), "\n")
  cat("  covariates in the model:", paste(colnames(X)[-1], collapse = " "), "\n")
  sexcols <- grep("sex", colnames(X), ignore.case = TRUE, value = TRUE)
  cat("  columns whose name contains 'sex':", if (length(sexcols)) paste(sexcols, collapse = " ") else "none", "\n")
  if (!anyNA(gap)) {
    total <- mean(eta[tail]) - mean(eta[!tail])
    o <- order(-abs(gap))[1:min(4, length(gap))]
    cat("  the tail's linear predictor sits", num(total), "log-odds from the rest's; carried by:\n")
    for (k in o) cat(sprintf("    %-32s %9s log-odds (%5.1f %% of the gap), %5.1f %% of the variance of the linear predictor\n",
                             name[k], num(gap[k]), 100 * gap[k] / total, 100 * share[k]))
  } else {
    cat("  tail statistics withheld: the tail or the rest is smaller than", mask, "\n")
  }
  for (k in seq_along(name)[-1]) {
    rows[[length(rows) + 1]] <- c(model = basename(f), N = n, cases = show(sum(y == 1)), n_tail = show(nt),
                                  covariate = name[k], order_in_model = k - 1,
                                  share_of_variance = num(share[k]),
                                  tail_minus_rest = if (is.na(gap[k])) "withheld" else num(gap[k]))
  }
}
if (length(rows) == 0) { cat("no binary model read\n"); quit(status = 1) }
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
utils::write.table(tab, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("binary models read:", length(unique(tab$model)), "\n")
RS
out="${OUT_DIR}/mu_drivers.${TAG}.tsv"
"${RSCRIPT}" mu_drivers.R models.list "${out}" "${TAIL}" "${MASK}"
gzip -f "${out}"
echo done
