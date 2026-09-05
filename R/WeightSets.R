# The per-weight-set descriptor region_postprocess consumes (design D2, Phase
# 5): slot r describes C++ set r and aligns with the Weight labels handed to
# mainRegionInCPP. `kind` drives set r's R-side AnnoWeights --
#   custom   a group-file weight column (customCol indexes WEIGHT)
#   beta     dbeta(MAF, a, b); the rows of `betaPairs` fill these in order,
#            so the implicit single set of an unweighted default run is a
#            beta set holding the first weights.beta pair
#   noweight all ones
# -- and `force` (is_output_weight_column) puts a K = 1 run on the
# multi-weight output schema.
weightSetDescriptor <- function(kind, name, betaPairs = NULL, force = FALSE) {
  if (length(kind) != length(name) || !all(kind %in% c("custom", "beta", "noweight"))) {
    stop("weightSetDescriptor: kind and name must align, kind in custom/beta/noweight\n")
  }
  nBeta <- sum(kind == "beta")
  if (nBeta > 0 && (!is.matrix(betaPairs) || nrow(betaPairs) < nBeta || ncol(betaPairs) != 2)) {
    stop("weightSetDescriptor: ", nBeta, " beta set(s) need at least that many (a, b) rows\n")
  }
  betaA <- rep(NA_real_, length(kind))
  betaB <- rep(NA_real_, length(kind))
  betaA[kind == "beta"] <- betaPairs[seq_len(nBeta), 1]
  betaB[kind == "beta"] <- betaPairs[seq_len(nBeta), 2]
  list(
    kind = kind,
    name = name,
    # column index into the WEIGHT matrix/vector, valid where kind == "custom"
    customCol = cumsum(kind == "custom"),
    betaA = betaA,
    betaB = betaB,
    force = isTRUE(force)
  )
}
