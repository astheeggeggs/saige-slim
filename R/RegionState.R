# Phase G1 (2026-09-14), slice 1 of SAIGE.Region: the per-region state
# record. Everything the region's later stages read, built once so that a
# batch of regions can be held while the pool evaluates their cells (step
# 5). Moved verbatim from the scan loop, where it sat inline for ~100 lines;
# the stdout, the stop and the two warnings are the loop's, unchanged, and
# the goldens are byte-identical (every one takes this path).
#
# `region` is one entry of the parsed group file; `i` is the running count
# the loop keeps (the "Analyzing Region" line prints i - 1 of nRegions).
regionStateFor <- function(region, regionName, i, nRegions, nModels,
                           legacyNoWeight, nWeightSets, isMultiWeightSchema,
                           nTransformSets, flexLofAnno, scoreName) {
  SNP <- region$SNP
  if (legacyNoWeight) {
    # single-set is_no_weight run: the legacy all-ones vector (a
    # 1-column custom set for C++). When "No_Weights" is instead
    # APPENDED to other sets, t_includeNoWeights carries it and WEIGHT
    # keeps the group file's sets (or the c(0) sentinel).
    WEIGHT <- rep(1, length(SNP))
  } else if (nWeightSets <= 1) {
    # K = 0: region$WEIGHT is the c(0) sentinel and as.numeric() keeps
    # it, so the sum(WEIGHT) > 0 branches downstream still see 0.
    # K = 1: collapse the 1-column numeric matrix to the plain unnamed
    # vector the legacy path produced (as.numeric drops dim/dimnames).
    WEIGHT <- as.numeric(region$WEIGHT)
  } else {
    # K > 1: the q0 x K numeric matrix, colnames = the group file's weight
    # set names. It stays a matrix all the way through
    # run_mainRegionInCPP and region_postprocess, which both key
    # multi-weight mode off is.matrix(WEIGHT)/ncol(WEIGHT).
    WEIGHT <- region$WEIGHT
  }

  # Under the multi-weight schema ACAT-V hands custom weights to CCT as
  # Cauchy combination weights, which stop()s mid-run on a negative.
  # checkGroupFile already rejects negatives for K > 1 group FILES, but
  # a K <= 1 weighted file can still reach the schema through
  # composition (is_no_weight_in_groupTest / multiple weights.beta
  # pairs) or is_output_weight_column -- so guard here per region with a
  # clear message. Never triggers on legacy paths (single set, no
  # force): those keep the historical no-sign-check behavior.
  if (isMultiWeightSchema && any(WEIGHT < 0)) {
    stop(
      "Region ", regionName, ": negative marker weights cannot be ",
      "used with the multi-weight output schema (ACAT-V uses them as ",
      "Cauchy combination weights). Please use non-negative weights ",
      "in the group file.\n"
    )
  }

  print(paste0("Analyzing Region ", regionName, " (", i - 1, "/", nRegions, ")."))

  # FlexRV: the region's score and LoF flags (numeric(0) when no
  # transform set is composed; C++ then ignores them)
  SCORE <- if (nTransformSets > 0) as.numeric(region$SCORE) else numeric(0)
  # without a lof transform no label list was parsed and no ISLOF was
  # carried; C++ still wants one flag per marker (never read)
  ISLOF <- if (nTransformSets == 0) {
    numeric(0)
  } else if (length(flexLofAnno) > 0) {
    as.numeric(region$ISLOF)
  } else {
    rep(0, length(SCORE))
  }
  if (nTransformSets > 0 && length(SCORE) > 1) {
    # a property of the input, warned about rather than acted on: every
    # score transform of a constant column coincides (every set of a
    # score transform: 192 -> 12, or 112 -> 7 at N = 1,000),
    # and 0/1 are fixed points of every transform (16 -> at most 3)
    us <- unique(SCORE)
    if (length(us) == 1) {
      warning("Region ", regionName, ": the score '", scoreName, "' is constant (",
              us, ") over its ", length(SCORE), " markers, so every score transform ",
              "gives the same weights and the FlexRV score axis is degenerate\n",
              call. = FALSE, immediate. = TRUE)
    } else if (all(us %in% c(0, 1))) {
      warning("Region ", regionName, ": the score '", scoreName, "' takes only the ",
              "values 0 and 1, which every score transform maps to themselves, so ",
              "the nonlinear FlexRV score transforms coincide\n",
              call. = FALSE, immediate. = TRUE)
    }
  }
  list(
    regionName = regionName,
    i = i,
    SNP = SNP,
    WEIGHT = WEIGHT,
    SCORE = SCORE,
    ISLOF = ISLOF,
    annolistsub = region$annoVec,
    annoIndicatorMat = region$annoIndicatorMat,
    region = region,
    cctpvalK = rep(NA_real_, nModels),
    pvalRegionK = vector("list", nModels),
    resultDFListK = vector("list", nModels) # the per-set carry-forward cache, per model
  )
}
