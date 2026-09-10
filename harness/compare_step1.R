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

# Elementwise symmetric relative error. For the model's vector and matrix
# fields the denominator is floored at the field's scale (its largest
# magnitude): an entry a million-fold below that scale is determined only to
# the solve's working precision TIMES THE SCALE, so a floorless ratio there
# reads the last bits of a near-zero number. LEDGER #124's gate moved
# obj.noK$S_a entries sitting at 1e-8 of scale by 1e-9 and the floorless ratio
# called that 15% while every field's error against its scale was <= 2.4e-4.
# At the largest entry the floor is inactive (|a| + |b| >= scale there). theta
# and the variance ratios are compared WITHOUT the floor (floor = FALSE): their
# components are different quantities, and with tau1 fixed at 1 a floor at the
# vector's scale would let a 3% move of a binary tau2 through the 1e-2 gate
# (drilled 2026-09-10).
relerr <- function(a, b, floor = TRUE) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) return(Inf)
  scale <- if (floor) suppressWarnings(max(abs(c(a, b)), na.rm = TRUE)) else 0
  d <- abs(a - b) / pmax(abs(a) + abs(b), scale, .Machine$double.eps)
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
if (!is.null(g$theta)) report("theta", relerr(g$theta, c_$theta, floor = FALSE), TOL_THETA)

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
  } else if (is.list(gv) && nm == "LOCOResult") {
    # A LOCO fit's whole point lives in here: one refit per left-out
    # chromosome. Before this branch existed the loop skipped LOCOResult
    # silently -- it is a list and it is not obj.noK -- so a LOCO step-1
    # golden would have gated theta and the base fields while comparing not
    # one value of the thing it was added for.
    if (length(gv) != length(cv)) {
      cat("LOCOResult: length differs FAIL\n")
      fails <- fails + 1L
    } else {
      gflag <- vapply(gv, function(e) isTRUE(e$isLOCO), logical(1))
      cflag <- vapply(cv, function(e) isTRUE(e$isLOCO), logical(1))
      if (!identical(gflag, cflag)) {
        cat("LOCOResult: isLOCO pattern differs FAIL\n")
        fails <- fails + 1L
      }
      # Only the populated entries carry numbers; the rest are list(isLOCO =
      # FALSE) stubs for chromosomes the fixture has no markers on.
      for (j in which(gflag & cflag)) {
        ge <- gv[[j]]
        ce <- cv[[j]]
        for (sub in intersect(names(ge), names(ce))) {
          if (is.numeric(ge[[sub]]) && length(ge[[sub]]) > 0) {
            report(sprintf("LOCOResult[[%d]]$%s", j, sub),
                   relerr(ge[[sub]], ce[[sub]]), TOL_OTHER)
          } else if (is.list(ge[[sub]]) && sub == "obj.noK") {
            for (k in intersect(names(ge[[sub]]), names(ce[[sub]]))) {
              if (is.numeric(ge[[sub]][[k]])) {
                report(sprintf("LOCOResult[[%d]]$obj.noK$%s", j, k),
                       relerr(ge[[sub]][[k]], ce[[sub]][[k]]), TOL_OTHER)
              }
            }
          }
        }
      }
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
    report("varianceRatio values", relerr(gvr[[1]], cvr[[1]], floor = FALSE), TOL_THETA)
  }
} else if (file.exists(gvr_f) != file.exists(cvr_f)) {
  cat("varianceRatio file presence differs FAIL\n"); fails <- fails + 1L
}

if (fails > 0L) { cat(sprintf("== step1 compare: %d FAILURES (%s)\n", fails, basename(gold_prefix))); quit(status = 1) }
cat(sprintf("== step1 compare: all PASS (%s)\n", basename(gold_prefix)))
