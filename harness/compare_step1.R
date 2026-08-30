#!/usr/bin/env Rscript
# Toleranced comparison of two step-1 outputs (golden vs candidate).
# Usage: Rscript compare_step1.R <golden_prefix> <candidate_prefix>
# Compares <prefix>.rda (or <prefix>_noLOCO.rda) and <prefix>.varianceRatio.txt.
#
# Tolerances (see plan + LEDGER.md #6-#8 for why they cannot be tighter):
#   theta/tau, variance ratios : 1e-2 relative  (AI-REML tol=0.02, ratioCV design)
#   other numeric model fields : 1e-3 relative  (float PCG at tolPCG=1e-5)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
gold_prefix <- args[1]; cand_prefix <- args[2]

TOL_THETA <- 1e-2
TOL_OTHER <- 1e-3

find_rda <- function(prefix) {
  for (f in c(paste0(prefix, ".rda"), paste0(prefix, "_noLOCO.rda"))) {
    if (file.exists(f)) return(f)
  }
  stop("no .rda found for prefix: ", prefix)
}

load_modglmm <- function(path) {
  e <- new.env()
  load(path, envir = e)
  stopifnot(exists("modglmm", envir = e))
  get("modglmm", envir = e)
}

relerr <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) return(Inf)
  d <- abs(a - b) / (abs(a) + abs(b) + .Machine$double.eps)
  d[is.na(a) & is.na(b)] <- 0
  max(d, na.rm = TRUE)
}

fails <- 0L
report <- function(name, err, tol) {
  ok <- is.finite(err) && err <= tol
  cat(sprintf("%-28s max_rel_err=%.3e tol=%.0e %s\n", name, err, tol, if (ok) "PASS" else "FAIL"))
  if (!ok) fails <<- fails + 1L
}

g <- load_modglmm(find_rda(gold_prefix))
c_ <- load_modglmm(find_rda(cand_prefix))

# theta (variance components) at the loose tolerance
if (!is.null(g$theta)) report("theta", relerr(g$theta, c_$theta), TOL_THETA)

# every other shared numeric field at the tighter tolerance
skip <- c("theta")
shared <- intersect(names(g), names(c_))
for (nm in setdiff(shared, skip)) {
  gv <- g[[nm]]; cv <- c_[[nm]]
  if (is.numeric(gv) && is.numeric(cv) && length(gv) > 0) {
    report(nm, relerr(gv, cv), TOL_OTHER)
  } else if (is.list(gv) && nm == "obj.noK") {
    for (sub in intersect(names(gv), names(cv))) {
      if (is.numeric(gv[[sub]])) report(paste0("obj.noK$", sub), relerr(gv[[sub]], cv[[sub]]), TOL_OTHER)
    }
  }
}
only_g <- setdiff(names(g), names(c_)); only_c <- setdiff(names(c_), names(g))
if (length(only_g)) { cat("FIELDS ONLY IN GOLDEN:", paste(only_g, collapse=", "), "\n"); fails <- fails + 1L }
if (length(only_c)) { cat("FIELDS ONLY IN CANDIDATE:", paste(only_c, collapse=", "), "\n"); fails <- fails + 1L }

# variance ratio file: numeric col toleranced, tag columns exact
gvr_f <- paste0(gold_prefix, ".varianceRatio.txt")
cvr_f <- paste0(cand_prefix, ".varianceRatio.txt")
if (file.exists(gvr_f) && file.exists(cvr_f)) {
  gvr <- read.table(gvr_f, stringsAsFactors = FALSE)
  cvr <- read.table(cvr_f, stringsAsFactors = FALSE)
  if (nrow(gvr) != nrow(cvr) || !identical(gvr[, -1], cvr[, -1])) {
    cat("varianceRatio: row/tag structure differs FAIL\n"); fails <- fails + 1L
  } else {
    report("varianceRatio values", relerr(gvr[[1]], cvr[[1]]), TOL_THETA)
  }
} else if (file.exists(gvr_f) != file.exists(cvr_f)) {
  cat("varianceRatio file presence differs FAIL\n"); fails <- fails + 1L
}

if (fails > 0L) { cat(sprintf("== step1 compare: %d FAILURES (%s)\n", fails, basename(gold_prefix))); quit(status = 1) }
cat(sprintf("== step1 compare: all PASS (%s)\n", basename(gold_prefix)))
