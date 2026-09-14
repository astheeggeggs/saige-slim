# FlexRV (Schwartzentruber et al., bioRxiv 2025.11.04.686320) on the
# multi-weight region engine: the transform table, its names, and the R
# closed forms. The definitive table, with each form's derivation and the
# ambiguities the paper leaves open, is docs/design/flexrv-design.md; the
# C++ kernels that the pipeline actually evaluates are src/FlexRV.cpp. R never
# re-derives a weight the pipeline uses -- flexrvWeightsR() below exists so the
# unit tests can compute every transform a second way and reconcile.
#
# Kind codes, shared with src/FlexRV.hpp (FlexRV::ScoreKind / MafKind):
#   score  0 unw   1 lof   2 cbrt (param t)   3 sbeta (param b)
#   MAF    0 unw   1 mafbeta (param b)        2 mac (param t)

flexrvScoreKinds <- c(unw = 0L, lof = 1L, cbrt = 2L, sbeta = 3L)
flexrvMafKinds <- c(unw = 0L, mafbeta = 1L, mac = 2L)

# The paper's 16 score and 12 MAF transforms, as the spec strings the CLI
# defaults to (--flexRV_scoreTransforms / --flexRV_MAFTransforms).
flexrvDefaultScoreSpec <- "unw;lof;cbrt:0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9;sbeta:1.5,2,3,5,10,20"
flexrvDefaultMAFSpec <- "unw;mafbeta:2,5;mac:1,3,5,10,25,50,100,200,500"

# A spec string is ';'-separated families; a family is a bare name (unw, lof)
# or name:p1,p2,... (cbrt, sbeta, mafbeta, mac). Returns one row per
# transform, in spec order: kind code, parameter (NA for the bare ones) and
# the name suffix used in the set label (cbrt0.2, sbeta1.5, mafbeta2, MAC10).
parseFlexRVTransformSpec <- function(spec, axis = c("score", "maf")) {
  axis <- match.arg(axis)
  kinds <- if (axis == "score") flexrvScoreKinds else flexrvMafKinds
  bare <- if (axis == "score") c("unw", "lof") else "unw"
  flag <- if (axis == "score") "--flexRV_scoreTransforms" else "--flexRV_MAFTransforms"
  if (!is.character(spec) || length(spec) != 1 || is.na(spec) || !nzchar(trimws(spec))) {
    stop(flag, " must be a non-empty string such as '",
         if (axis == "score") flexrvDefaultScoreSpec else flexrvDefaultMAFSpec, "'\n")
  }
  families <- trimws(strsplit(spec, ";", fixed = TRUE)[[1]])
  families <- families[nzchar(families)]
  rows <- lapply(families, function(fam) {
    parts <- strsplit(fam, ":", fixed = TRUE)[[1]]
    name <- trimws(parts[1])
    if (!name %in% names(kinds)) {
      stop(flag, ": unknown transform '", name, "' in '", fam, "'; the ", axis,
           " transforms are ", paste(names(kinds), collapse = ", "), "\n")
    }
    if (name %in% bare) {
      if (length(parts) != 1) {
        stop(flag, ": '", name, "' takes no parameters (got '", fam, "')\n")
      }
      return(data.frame(kind = kinds[[name]], param = NA_real_, suffix = name,
                        stringsAsFactors = FALSE))
    }
    if (length(parts) != 2 || !nzchar(trimws(parts[2]))) {
      stop(flag, ": '", name, "' needs a comma-separated parameter list, e.g. '", name,
           ":", if (name == "cbrt") "0.2,0.5" else if (name == "mac") "1,3,5" else "2,5", "'\n")
    }
    vals <- suppressWarnings(as.numeric(strsplit(parts[2], ",", fixed = TRUE)[[1]]))
    if (length(vals) == 0 || anyNA(vals) || any(!is.finite(vals))) {
      stop(flag, ": non-numeric parameter in '", fam, "'\n")
    }
    if (anyDuplicated(vals) > 0) {
      stop(flag, ": duplicated parameter in '", fam, "'\n")
    }
    ok <- switch(name,
      cbrt = all(vals > 0 & vals < 1),
      sbeta = all(vals > 1),
      mafbeta = all(vals > 1),
      mac = all(vals >= 1)
    )
    if (!ok) {
      stop(flag, ": '", name, "' parameters must be ",
           switch(name,
             cbrt = "transition values strictly inside (0, 1)",
             sbeta = "beta shape parameters b > 1",
             mafbeta = "beta shape parameters b > 1",
             mac = "MAC thresholds >= 1"
           ), " (got '", fam, "')\n")
    }
    data.frame(kind = kinds[[name]], param = vals,
               suffix = paste0(if (name == "mac") "MAC" else name, as.character(vals)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) {
    stop(flag, " selects no transform\n")
  }
  if (anyDuplicated(out$suffix) > 0) {
    stop(flag, ": transform '", out$suffix[duplicated(out$suffix)][1], "' is listed twice\n")
  }
  rownames(out) <- NULL
  out
}

# The Kt x 4 transform table for C++ (scoreKind, scoreParam, mafKind,
# mafParam), score-major and MAF-minor: row k = c*nMAF + d (0-based) pairs
# score transform c with MAF transform d, the NORMATIVE order shared with
# src/Main.cpp. Bare transforms carry parameter 0 (never read). The set names
# are <scoreName>_<scoreSuffix>_<mafSuffix>; when maxMAF is not the paper's
# 0.001 the mafbeta suffixes carry '@<maxMAF>', because those transforms
# rescale to the inclusion boundary and the OUTPUT must say so.
#
# LEDGER #142 (decided 2026-09-13, option b): given the in-sample N, the MAC
# rungs that would weight EVERY marker 1 (a duplicate of `unw`) are dropped
# BEFORE the table is built; a rung that selects nothing (all zeros) is kept.
# Region ingest admits a marker only if its folded MAF is at or below maxMAF,
# so every admitted MAC is at most 2*N*maxMAF, and a rung t at or above that
# bound gives every marker weight 1 -- the `unw` MAF transform under another
# name, evaluated again and counted again in the pooled row. At the paper's N
# (~424k, maxMAF 0.001, bound ~847) no default rung is degenerate and the
# table is the full 16 x 12; at or below N = 250,000 rungs go one by one (t
# is degenerate for every N <= t / (2 * maxMAF), equality included). The ladder's values stay
# absolute (docs/design/flexrv-plan.md); what is dropped is only what cannot
# differ from `unw` at this N. The comparison is `>=` in doubles: the kernel
# is `MAC <= t` and a hard-call MAC is an integer, so a rung equal to the
# bound is degenerate too and no float edge can admit a marker above it.
# `dropped` names the rungs that went, for the caller's warning; N = NULL
# (the unit tests, the name builders) keeps the full ladder.
flexrvTransformTable <- function(scoreSpec, mafSpec, scoreName, maxMAF = 0.001, N = NULL) {
  sc <- parseFlexRVTransformSpec(scoreSpec, "score")
  mf <- parseFlexRVTransformSpec(mafSpec, "maf")
  if (!isTRUE(maxMAF > 0 & maxMAF <= 0.5)) {
    stop("--flexRV_maxMAF must be in (0, 0.5]\n")
  }
  dropped <- character(0)
  bound <- NA_real_
  if (!is.null(N)) {
    if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N < 1) {
      stop("flexrvTransformTable: N must be the in-sample sample count (got ",
           deparse(N), ")\n")
    }
    bound <- 2 * N * maxMAF
    degenerate <- mf$kind == flexrvMafKinds[["mac"]] & mf$param >= bound
    dropped <- mf$suffix[degenerate]
    if (all(degenerate)) {
      stop("--flexRV_MAFTransforms='", mafSpec, "': every MAF transform is a MAC rung at or ",
           "above 2 * N * maxMAF = ", bound, " (N = ", N, ", --flexRV_maxMAF=", maxMAF,
           "), which no admitted marker can exceed, so every set would be the unweighted ",
           "one and nothing is left to test. Add 'unw' or a rung below ", bound, "\n")
    }
    mf <- mf[!degenerate, , drop = FALSE]
  }
  mfSuffix <- mf$suffix
  if (maxMAF != 0.001) {
    isBeta <- mf$kind == flexrvMafKinds[["mafbeta"]]
    mfSuffix[isBeta] <- paste0(mfSuffix[isBeta], "@", as.character(maxMAF))
  }
  nS <- nrow(sc)
  nM <- nrow(mf)
  ci <- rep(seq_len(nS), each = nM)
  di <- rep(seq_len(nM), times = nS)
  param0 <- function(p) ifelse(is.na(p), 0, p)
  table <- cbind(sc$kind[ci], param0(sc$param[ci]), mf$kind[di], param0(mf$param[di]))
  dimnames(table) <- NULL
  storage.mode(table) <- "double"
  list(
    table = table,
    names = paste0(scoreName, "_", sc$suffix[ci], "_", mfSuffix[di]),
    scoreSuffix = sc$suffix,
    mafSuffix = mfSuffix,
    nScore = nS,
    nMAF = nM,
    usesLof = any(sc$kind == flexrvScoreKinds[["lof"]]),
    maxMAF = maxMAF,
    N = N,
    bound = bound,
    dropped = dropped
  )
}

# The R route to every weight, for reconciliation with flexrvWeightsInCPP()
# (the C++ kernels in src/FlexRV.cpp). Deliberately written from the PAPER's
# dbeta forms and the cubic-root formula verbatim, not from the closed forms
# the C++ uses, so that agreement is evidence about the derivation.
flexrvWeightsR <- function(table, score, isLof, MAF, MAC, maxMAF = 0.001) {
  q <- length(score)
  stopifnot(length(isLof) == q, length(MAF) == q, length(MAC) == q, ncol(table) == 4)
  scoreW <- function(kind, param) {
    switch(as.character(kind),
      "0" = rep(1, q),
      "1" = as.numeric(isLof),
      "2" = {
        t <- param
        (sign(score - t) * abs(score - t)^(1 / 3) + t^(1 / 3)) / (abs(1 - t)^(1 / 3) + t^(1 / 3))
      },
      "3" = dbeta(1 - abs(score), 1, param) / dbeta(0, 1, param),
      stop("flexrvWeightsR: unknown score kind ", kind)
    )
  }
  mafW <- function(kind, param) {
    switch(as.character(kind),
      "0" = rep(1, q),
      "1" = {
        x <- MAF / maxMAF
        ifelse(x < 1, dbeta(x, 1, param) / dbeta(0, 1, param), 0)
      },
      "2" = as.numeric(MAC <= param),
      stop("flexrvWeightsR: unknown MAF kind ", kind)
    )
  }
  W <- matrix(NA_real_, nrow = q, ncol = nrow(table))
  for (k in seq_len(nrow(table))) {
    W[, k] <- scoreW(table[k, 1], table[k, 2]) * mafW(table[k, 3], table[k, 4])
  }
  W
}

# The LoF label list for --flexRV_lofAnno: ',', ':' and ';' are
# interchangeable separators (one flat dimension, unlike
# --annotation_in_groupTest).
parseFlexRVLofAnno <- function(lofAnno) {
  if (is.null(lofAnno)) return(character(0))
  labs <- trimws(unlist(strsplit(lofAnno, "[,:;]")))
  labs <- labs[nzchar(labs)]
  if (length(labs) == 0) {
    stop("--flexRV_lofAnno must name at least one annotation label\n")
  }
  unique(labs)
}
