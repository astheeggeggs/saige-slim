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
    # LEDGER #160: whether WEIGHT carries CUSTOM sets at all, recorded here
    # because only this function can tell. The three branches above produce a
    # plain numeric vector in two different senses -- a GENUINE one-column set
    # (legacyNoWeight's all-ones, or a K = 1 matrix collapsed by as.numeric)
    # and the K = 0 `c(0)` sentinel -- and at q0 == 1 those are the same 1 x 1
    # matrix by the time C++ sees them, which is the collision #160 fixes.
    hasCustomWeights = isTRUE(legacyNoWeight) || nWeightSets >= 1L,
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

# Phase G1 slice 2: the test plan r.corr selects. Returns r.corr as well as
# the method, because SKAT_Check_Method promotes 0 to the 7-element rho grid
# c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 1) that skat_fun, region_postprocess and
# the batch pipeline all consume; a plan that dropped it would move every
# SKAT-O golden or die. The two stdout lines are the loop's, unchanged.
regionTestPlan <- function(r.corr, is_single_in_groupTest) {
  if (r.corr == 0) {
    out.method <- SKAT:::SKAT_Check_Method(method = "optimal.adj", r.corr = 0)
    cat(
      "SKAT-O test will be performed. P-values for BURDEN and SKAT tests will also be output\n"
    )
    list(method = out.method$method, r.corr = out.method$r.corr,
         regionTestType = "SKAT-O", is_single_in_groupTest = TRUE)
  } else if (r.corr == 1) {
    cat("BURDEN test will be performed\n")
    # the BURDEN output file (written from Rcpp) is opened after the group
    # file check in SAIGE.Region: its header needs the number of weight sets K
    list(method = NULL, r.corr = r.corr,
         regionTestType = "BURDEN", is_single_in_groupTest = is_single_in_groupTest)
  } else {
    stop("r.corr needs to be either 1 (BURDEN test) or 0 (SKAT-O test)\n")
  }
}

# Phase G1 slice 3: the value SAIGE.Region returns, which the CLI prints.
# The sidecar summaries (#57, #127, the memo lines) are printed by the
# caller before this is assembled, as before.
regionDoneMessage <- function(OutputFiles, OutputFile, writeMarkerList,
                              is_single_in_groupTest) {
  message <- paste0(
    "Analysis done!",
    " The set-based tests results have been saved to '",
    paste(OutputFiles, collapse = "', '"),
    "'."
  )
  if (writeMarkerList) {
    message <- paste0(message, " The marker lists have been saved to '", OutputFile,
                      ".markerList.txt'.")
  }
  if (is_single_in_groupTest) {
    message <- paste0(message,
                      " The single-variant association tests results have been saved to '",
                      OutputFile, ".singleAssoc.txt'.")
  }
  message
}

# Phase G1 slice 4: the composed per-run weight-set list, in the NORMATIVE
# order shared with C++ (src/Main.cpp mainRegionInCPP):
#   custom sets (group-file weight lines, file order)
#   -> FlexRV transform sets (F3; score-major over the transform table)
#   -> explicit beta sets (weights.beta pairs in given order, named
#      "Beta_a_b"; composed ONLY when more than one pair is supplied -- a
#      single pair, default or not, keeps the legacy implicit-beta path)
#   -> "No_Weights" (all-ones; appended when is_no_weight_in_groupTest =
#      TRUE and the run has any other weight set)
# At exactly one composed set every consumer reduces to the legacy
# single-weight flow (same objects, same bytes). Moved verbatim from
# SAIGE.Region (design D5/D10); the three stdout lines and the two stops are
# unchanged. Pure apart from those: no file, no C++ state.
composeWeightSets <- function(weights_beta_mat, is_no_weight_in_groupTest, is_weight_included,
                              nTransformSets, nWeightSets, weightSetNames, flexNames,
                              is_output_weight_column) {
  nBetaPairs <- nrow(weights_beta_mat)
  useExplicitBeta <- nBetaPairs > 1
  betaSetNames <- if (useExplicitBeta) {
    paste0("Beta_", weights_beta_mat[, 1], "_", weights_beta_mat[, 2])
  } else {
    character(0)
  }
  # is_no_weight_in_groupTest composition:
  #  - no other weight set (unweighted group file, single weights.beta pair):
  #    the LEGACY path verbatim -- the single all-ones set replaces the
  #    implicit beta set, no new output columns
  #  - group-file weight lines (which used to hard-error here) or explicit
  #    beta sets: an all-ones "No_Weights" set is APPENDED as the last set
  includeNoWeightsSet <- is_no_weight_in_groupTest &&
    (is_weight_included || useExplicitBeta || nTransformSets > 0)
  legacyNoWeight <- is_no_weight_in_groupTest && !includeNoWeightsSet
  if (legacyNoWeight) {
    cat("No weights are used in the group test\n")
  }
  if (includeNoWeightsSet) {
    cat("is_no_weight_in_groupTest = TRUE: an all-ones 'No_Weights' weight set is appended",
        " to the group test\n")
  }
  # LEDGER #138 FIXED: a FlexRV run tests ONE cell with its own transform sets
  # and its pooled row is p_FlexRV, the Cauchy combination over exactly those.
  # Any other set composed beside them -- a group-file weight line, explicit
  # --weights.beta pairs, or the No_Weights set -- silently made the pooled
  # row a combination over 193+ rows. Refused here, beside the three FlexRV
  # cell guards in R/SAIGE_Test_main.R, because the weight lines are known
  # only once the group file is read. (A score line WITHOUT --flexRV_score is
  # ignored; a weight line WITH it is not, hence the explicit refusal.)
  if (nTransformSets > 0 && (nWeightSets > 0 || useExplicitBeta || includeNoWeightsSet)) {
    extra <- c(
      if (nWeightSets > 0) paste0("group-file weight set(s) ",
                                  paste(weightSetNames, collapse = ", ")),
      if (useExplicitBeta) paste0(nBetaPairs, " --weights.beta pairs"),
      if (includeNoWeightsSet) "--is_no_weight_in_groupTest=TRUE"
    )
    stop("A FlexRV run tests ONE cell with its ", nTransformSets, " transform sets, and its ",
         "pooled row is p_FlexRV, the Cauchy combination over exactly those; ",
         paste(extra, collapse = " and "), " would be composed beside them and the pooled ",
         "row would no longer be p_FlexRV. Drop the weight line(s) from the group file, or ",
         "run the custom-weight test separately without --flexRV_score\n")
  }
  composedNames <- c(
    if (legacyNoWeight) "No_Weights",
    if (nWeightSets > 0) weightSetNames,
    flexNames,
    betaSetNames,
    if (includeNoWeightsSet) "No_Weights"
  )
  nWeightSetsUsed <- max(1L, length(composedNames))
  isMultiWeightRun <- nWeightSetsUsed > 1
  if (isMultiWeightRun) {
    if (anyDuplicated(composedNames) > 0) {
      stop(
        "Duplicate weight set name '",
        composedNames[duplicated(composedNames)][1],
        "' in the composed weight-set list (group-file sets, then weights.beta pairs, ",
        "then No_Weights). Please rename the group-file weight set or drop the duplicate pair.\n"
      )
    }
    cat(
      nWeightSetsUsed, " weight sets are used in the group test: ",
      paste(composedNames, collapse = ", "), "\n"
    )
  }
  # Conditional analysis composes with the multi-weight schema (design D7,
  # C11): the conditioning-marker weight vector stays a SINGLE vector across
  # all sets (weights_for_condition, else pdf(beta(1,25)) of the conditioning
  # MAFs); what varies per set is the G1-side weight, applied in
  # region_postprocess from the raw exported G1tilde_P_G2tilde matrix.
  # Schema switch for the region output file: the Weight column (and, for
  # SKAT-O, Pvalue_ACATV/Pvalue_ACATO) is emitted when more than one weight
  # set is composed OR the user forces it at a single set.
  isMultiWeightSchema <- isMultiWeightRun || is_output_weight_column
  list(
    nBetaPairs = nBetaPairs,
    useExplicitBeta = useExplicitBeta,
    betaSetNames = betaSetNames,
    includeNoWeightsSet = includeNoWeightsSet,
    legacyNoWeight = legacyNoWeight,
    composedNames = composedNames,
    nWeightSetsUsed = nWeightSetsUsed,
    isMultiWeightRun = isMultiWeightRun,
    isMultiWeightSchema = isMultiWeightSchema
  )
}
