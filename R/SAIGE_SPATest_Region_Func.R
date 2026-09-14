#############################
#
# ReadGroupFile (its own error branch formatted an undefined `n2`, proof no
# input ever reached it) and check_close were deleted 2026-09-14 (dead-surface
# sweep C6); the group file is read by the C++ region loop.

# check Group file format and count the number of regions
# With detail = FALSE (the default) the return value is exactly
# list(nRegions, is_weight_included), the historical contract.
# With detail = TRUE three fields are appended: nWeightSets (K, the number of
# weight lines per region, defined by region 1), weightSetNames (character
# vector of length K; a bare 'weight' type is named "WEIGHT"). A fourth,
# scoreName, is appended ONLY when the file carries a FlexRV score line (a bare
# 'score' type is named "SCORE"); a file without one returns the historical
# four fields unchanged.
checkGroupFile <- function(groupFile, detail = FALSE) {
  cat("Start extracting marker-level information from 'groupFile' of", groupFile, "....\n")
  Check_File_Exist(groupFile, "RegionFile")

  # the 2nd field of a weight line is 'weight' (set name "WEIGHT"),
  # 'weight:NAME' or 'weight,NAME' (set name NAME); returns NULL for
  # non-weight lines
  parse_weight_type <- function(type) {
    type_parts <- strsplit(type, split = "[,:]")[[1]]
    if (length(type_parts) == 0 || type_parts[1] != "weight") {
      return(NULL)
    }
    if (length(type_parts) >= 2 && nzchar(type_parts[2])) {
      return(type_parts[2])
    }
    "WEIGHT"
  }

  # the 2nd field of a SCORE line is 'score' (name "SCORE"), 'score:NAME' or
  # 'score,NAME'; returns NULL for non-score lines. FlexRV (F1) reads one
  # per-marker pathogenicity score in [0, 1] per region and expands it into
  # weight sets in C++; nothing consumes it yet.
  #
  # POSITION IS FIXED: immediately after the anno line, before any weight
  # lines. The plan allowed "before or among" the weight lines; pinning it is
  # a deliberate tightening, because relaxing a position later is backward
  # compatible for existing files and tightening later would not be. It also
  # keeps every line offset in this parser and in SAIGE.getRegionList_new a
  # constant rather than a lookup into a per-region layout vector.
  parse_score_type <- function(type) {
    type_parts <- strsplit(type, split = "[,:]")[[1]]
    if (length(type_parts) == 0 || type_parts[1] != "score") {
      return(NULL)
    }
    if (length(type_parts) >= 2 && nzchar(type_parts[2])) {
      return(type_parts[2])
    }
    "SCORE"
  }

  # region 1 defines the format: a var line, an anno line, then K >= 0
  # consecutive weight lines; every region must follow the same layout
  gf <- file(groupFile, "r")
  marker_group_line <- readLines(gf, n = 1)
  marker_group_line <- readLines(gf, n = 1)
  is_weight_included <- FALSE
  weightSetNames <- character(0)
  scoreName <- character(0)
  # region 1's third line decides whether this file carries a score line
  peek <- readLines(gf, n = 1)
  if (length(peek) > 0) {
    peek_list <- strsplit(peek, split = "[\ \t]+")[[1]]
    if (length(peek_list) >= 3) {
      scoreName <- parse_score_type(peek_list[2])
      if (is.null(scoreName)) {
        scoreName <- character(0)
        pushBack(peek, gf)
      }
    } else {
      pushBack(peek, gf)
    }
  }
  repeat {
    marker_group_line <- readLines(gf, n = 1)
    if (length(marker_group_line) == 0) {
      break
    }
    marker_group_line_list <- strsplit(marker_group_line, split = "[\ \t]+")[[1]]
    if (length(marker_group_line_list) < 3) {
      break
    }
    if (!is.null(parse_score_type(marker_group_line_list[2]))) {
      stop("Error! A second 'score' line in region 1 of the group file, or a",
           " score line after a weight line: exactly one score line is allowed",
           " and it must come immediately after the anno line\n")
    }
    setName <- parse_weight_type(marker_group_line_list[2])
    if (is.null(setName)) {
      break
    }
    if (!is_weight_included) {
      is_weight_included <- TRUE
      print("weights are included for markers")
    }
    weightSetNames <- c(weightSetNames, setName)
  }
  close(gf)

  nWeightSets <- length(weightSetNames)
  if (anyDuplicated(weightSetNames) > 0) {
    stop("Error! Duplicate weight set name '", weightSetNames[duplicated(weightSetNames)][1], "' in region 1 of the group file\n")
  }
  if ("Cauchy" %in% weightSetNames) {
    stop("Error! 'Cauchy' is a reserved weight set name and cannot be used in the group file\n")
  }
  # one optional score line sits between anno and the weight lines
  nScoreSets <- length(scoreName)
  a <- 2 + nScoreSets + nWeightSets

  gf <- file(groupFile, "r")
  line <- 0
  nregion <- 0
  # every annotation label the file uses; collected only when a score line is
  # present, because that is when a caller (FlexRV's --flexRV_lofAnno) needs
  # to know whether a label matches anything in the WHOLE file before the
  # first chunk is read
  annoLabels <- character(0)
  while (TRUE) {
    marker_group_line <- readLines(gf, n = a)
    # line = line+a

    if (length(marker_group_line) == 0) {
      break
    } else {
      if (length(marker_group_line) < a) {
        marker_group_line_list <- strsplit(marker_group_line[1], split = "[\ \t]+")[[1]]
        if (length(marker_group_line_list) < 3) {
          stop("Error, group file line:", line - 2, " is incomplete.\n")
        }
        geneID <- marker_group_line_list[1]
        stop("Group file is incomplete for ", geneID, ".\n")
      }
    }


    marker_group_line_list <- strsplit(marker_group_line[1], split = "[\ \t]+")[[1]]
    line <- line + 1
    if (length(marker_group_line_list) < 3) {
      stop("Error, group file line:", line, " is incomplete.\n")
    }
    geneID <- marker_group_line_list[1]
    type <- marker_group_line_list[2]

    if (type != "var") {
      if (!is.null(parse_weight_type(type))) {
        stop("Error! Region ", geneID, " has more weight lines than region 1: every region must have the same ", nWeightSets, " weight line(s) in the same order\n")
      }
      stop("Error! No var is specified for ", geneID, "\n")
    }
    geneID0 <- geneID
    numMarkers <- length(marker_group_line_list) - 2

    marker_group_line_list <- strsplit(marker_group_line[2], split = "[\ \t]+")[[1]]
    line <- line + 1
    if (length(marker_group_line_list) < 3) {
      stop("Error, group file line:", line, " is incomplete.\n")
    }
    geneID <- marker_group_line_list[1]
    type <- marker_group_line_list[2]

    if (type != "anno") {
      stop("Error! No anno is specified for ", geneID, "\n")
    }
    if (geneID != geneID0) {
      stop("anno for ", geneID0, " is missing.\n")
    }
    numAnnos <- length(marker_group_line_list) - 2
    if (numAnnos != numMarkers) {
      stop("The length of annotations for markers in region ", geneID, " is not equal to the length of marker IDs\n")
    }
    if (nScoreSets > 0) {
      annoLabels <- union(annoLabels, marker_group_line_list[-c(1:2)])
    }

    if (nScoreSets > 0) {
      marker_group_line_list <- strsplit(marker_group_line[3], split = "[\ \t]+")[[1]]
      line <- line + 1
      if (length(marker_group_line_list) < 3) {
        stop("Error, group file line:", line, " is incomplete.\n")
      }
      geneID <- marker_group_line_list[1]
      thisScore <- parse_score_type(marker_group_line_list[2])
      if (is.null(thisScore)) {
        stop("Error! No score is specified for ", geneID0,
             ": region 1 has a score line, so every region must have one",
             " immediately after its anno line\n")
      }
      if (geneID != geneID0) {
        stop("score for ", geneID0, " is missing.\n")
      }
      if (thisScore != scoreName) {
        stop("Error! Score name '", thisScore, "' in region ", geneID0,
             " does not match score name '", scoreName,
             "' in region 1: every region must carry the same score\n")
      }
      numScores <- length(marker_group_line_list) - 2
      if (numScores != numMarkers) {
        stop("The length of scores for markers in region ", geneID,
             " is not equal to the length of marker IDs\n")
      }
      scoreValues <- suppressWarnings(as.numeric(marker_group_line_list[-c(1:2)]))
      if (anyNA(scoreValues)) {
        stop("Error! Non-numeric or NA score values in score '", scoreName,
             "' for region ", geneID, "\n")
      }
      # FlexRV's score transforms assume a normalised pathogenicity score;
      # cbrt/sbeta of a value outside [0, 1] is not defined by the method, so
      # reject at parse time with the region named rather than produce weights
      # nobody can interpret.
      if (any(scoreValues < 0 | scoreValues > 1)) {
        stop("Error! Score values outside [0, 1] in score '", scoreName,
             "' for region ", geneID, "\n")
      }
    }

    if (nWeightSets > 0) {
      for (k in 1:nWeightSets) {
        marker_group_line_list <- strsplit(marker_group_line[2 + nScoreSets + k],
                                           split = "[\ \t]+")[[1]]
        line <- line + 1
        if (length(marker_group_line_list) < 3) {
          stop("Error, group file line:", line, " is incomplete.\n")
        }
        geneID <- marker_group_line_list[1]
        type <- marker_group_line_list[2]
        setName <- parse_weight_type(type)
        if (is.null(setName)) {
          stop("Error! No weight is specified for ", geneID, "\n")
        }
        if (geneID != geneID0) {
          stop("weight for ", geneID0, " is missing.\n")
        }
        if (setName != weightSetNames[k]) {
          stop("Error! Weight set name '", setName, "' in region ", geneID0, " does not match weight set name '", weightSetNames[k], "' in region 1: weight lines must appear in the same order in every region\n")
        }
        numWeights <- length(marker_group_line_list) - 2
        if (numWeights != numMarkers) {
          stop("The length of weights for markers in region ", geneID, " is not equal to the length of marker IDs\n")
        }
        weightValues <- suppressWarnings(as.numeric(marker_group_line_list[-c(1:2)]))
        if (anyNA(weightValues)) {
          stop("Error! Non-numeric weight values in weight set '", weightSetNames[k], "' for region ", geneID, "\n")
        }
        # multi-weight mode only (K > 1): ACAT-V uses these values as Cauchy
        # combination weights, and CCT stop()s mid-run on any negative weight.
        # Fail fast at parse time with the file location instead. K <= 1 files
        # keep the historical behavior (weights enter only multiplicatively
        # and were never sign-checked upstream).
        if (nWeightSets > 1 && any(weightValues < 0)) {
          stop("Error! Negative weight values in weight set '", weightSetNames[k], "' for region ", geneID, ": weights must be non-negative in multi-weight mode\n")
        }
      }
    }
    nregion <- nregion + 1
  }
  # if(is.null(group_info_list[[geneID]])){
  #  group_info_list[[geneID]]<-list(geneID = geneID)
  # }
  # group_info_list[[geneID]][[type]] = marker_group_line_list[-c(1:2)]

  close(gf)
  res <- list(nRegions = nregion, is_weight_included = is_weight_included)
  if (detail) {
    res$nWeightSets <- nWeightSets
    res$weightSetNames <- weightSetNames
    # scoreName is appended ONLY when the file has a score line. Three tests
    # pin the historical detail = TRUE return as exactly
    # list(nRegions, is_weight_included, nWeightSets, weightSetNames), and the
    # FlexRV plan requires those to stay verbatim -- so a scoreless file must
    # keep returning four elements, not five with an empty fifth. Callers read
    # it as `region_list$scoreName`, which is NULL when absent, and
    # length(NULL) == 0 is exactly the arithmetic they want.
    if (nScoreSets > 0) {
      res$scoreName <- scoreName
      res$annoLabels <- annoLabels
    }
  }
  res
}

SPA_ER_kernel_related_Phiadj_fast_new <- function(p.new, Score, Phi, p.value_burden, regionTestType) {
  p.m <- length(Score)
  zscore.all_0 <- Score
  zscore.all_1 <- rep(0, p.m)
  if (regionTestType != "BURDEN") {
    VarS_org <- diag(as.matrix(Phi))
  } else {
    VarS_org <- as.vector(Phi)
  }
  stat.qtemp <- Score^2 / VarS_org


  idx_0 <- which(VarS_org > 0)
  # idx_p0<-which(p.new >0)
  idx_p0 <- which(!is.na(p.new))
  # idx_p1<-which(p.new <0)
  # if(length(idx_0) > 0){
  # 	zscore.all_1[idx_0]= qnorm(p.new[idx_0]/2, mean = 0, sd =sqrt(VarS_org[idx_0]),lower.tail = FALSE, log.p = FALSE)*sign(Score)
  # }
  VarS <- zscore.all_0^2 / 500

  if (length(idx_p0) > 0) {
    # VarS[idx_p0]= zscore.all_0[idx_p0]^2/qchisq(p.new[idx_p0], 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
    VarS[idx_p0] <- zscore.all_0[idx_p0]^2 / qchisq(p.new[idx_p0], df = 1, ncp = 0, lower.tail = FALSE, log.p = TRUE)
  }
  vars_inf <- which(VarS == Inf)
  if (regionTestType != "BURDEN") {
    if (length(vars_inf) > 0) {
      VarS[vars_inf] <- 0
      zscore.all_1[vars_inf] <- 0
      Phi[vars_inf, ] <- 0
      Phi[, vars_inf] <- 0
    }
  } else {
    if (length(vars_inf) > 0) {
      VarS[vars_inf] <- 0
      zscore.all_1[vars_inf] <- 0
      Phi[vars_inf] <- 0
    }
  }


  scaleFactor <- sqrt(VarS / VarS_org)
  ###################################
  # Burden test
  # print("SPA_ER_kernel_related_Phiadj_fast_new 5")
  zscore.all_1 <- zscore.all_0
  # print(as.matrix(Phi))
  # print(diag(as.vector(VarS/VarS_org)))
  VarStoorg <- as.vector(VarS / VarS_org)

  if (regionTestType != "BURDEN") {
    # t(sqrt(a)*t(b)*sqrt(a))
    G2_adj_n <- t(t(Phi * sqrt(VarStoorg)) * sqrt(VarStoorg))
    # v1=rep(1,nrow(G2_adj_n))
  } else {
    G2_adj_n <- Phi * VarStoorg
  }
  VarQ <- sum(G2_adj_n)
  Q_b <- sum(zscore.all_1)^2
  VarQ_2 <- Q_b / qchisq(p.value_burden, df = 1, ncp = 0, lower.tail = FALSE, log.p = TRUE)
  if (VarQ_2 == 0) {
    r <- 1
  } else {
    r <- VarQ / VarQ_2
  }
  r <- min(r, 1)

  if (regionTestType != "BURDEN") {
    Phi_ccadj <- as.matrix(G2_adj_n * 1 / r)
  } else {
    Phi_ccadj <- as.vector(G2_adj_n * 1 / r)
  }

  outlist <- list()
  outlist$val <- Phi_ccadj

  scaleFactor <- scaleFactor / sqrt(r)
  outlist$scaleFactor <- scaleFactor
  outlist$p.new <- p.new
  outlist
}


# LEDGER #143: this used to have an unsuffixed textual twin, get_newPhi_scaleFactor,
# 18 lines above -- the pre-#104 body without the standardisation below -- and
# that twin was the ONLY kernel the conditioning markers' scale factors went
# through (SAIGE_Test_main.R, the binary group-test --condition set-up). One
# rule, written twice, fixed once: the same shape as #44/#66 and #118/#134/#135.
# The twin is deleted; both call sites use this function.
get_newPhi_scaleFactor_traitType <- function(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType, traitType) {
  m1 <- sum(mu.a * g.sum)
  if (traitType == "binary") {
    var1 <- sum(mu.a * (1 - mu.a) * g.sum^2)
  }
  # LEDGER #104/#105 fix (fix A). var1 is the NAIVE burden variance -- no GRM,
  # no variance ratio -- and it is the variance of the saddlepoint's CGF over
  # (mu, g.sum). The burden statistic's variance under the null model is
  # 1'Phi 1 (GRM- and ratio-corrected, as C++ delivered it). Standardise the
  # statistic by that before the |z| < 2 test and the saddlepoint, exactly the
  # device the C++ marker SPA uses (q = Tstat / sqrt(var1 / var2) + m1 in
  # SAIGE_test.cpp): below the cutoff the implied VarQ_2 then equals 1'Phi 1
  # and the r = min(VarQ / VarQ_2, 1) rescale in
  # SPA_ER_kernel_related_Phiadj_fast_new is a no-op; above it, the rescale
  # carries only the case-control imbalance correction. Before this, VarQ_2
  # was var1 itself, so wherever the mixed model reduced the variance the
  # region path reported the naive burden test and inflated all of Phi by
  # var1 / 1'Phi 1 (every cell of golden P6, 5-6 %). With no GRM and unit
  # ratios 1'Phi 1 == var1 and the factor is exactly 1.
  var_b <- sum(Phi)
  q_adj <- m1 + (q.sum - m1) * sqrt(var1 / var_b)
  qinv <- -sign(q_adj - m1) * abs(q_adj - m1) + m1
  pval_noadj <- pchisq((q_adj - m1)^2 / var1, lower.tail = FALSE, df = 1, log.p = TRUE)
  isSPAConverge <- TRUE
  if (abs(q_adj - m1) / sqrt(var1) < 2) {
    p.value_burden <- pval_noadj
  } else {
    p.value_burden <- SPA_pval(mu.a, g.sum, q_adj, qinv, pval_noadj, .Machine$double.eps^0.25, TRUE, traitType, isSPAConverge)
  }

  re_phi <- SPA_ER_kernel_related_Phiadj_fast_new(p.new, Score, Phi, p.value_burden, regionTestType)
  re_phi
}


get_SKAT_pvalue <- function(Score, Phi, r.corr, regionTestType) {
  out_SKAT_List <- try(
    SKAT:::Met_SKAT_Get_Pvalue(
      Score = Score,
      Phi = Phi,
      r.corr = r.corr,
      method = "optimal.adj",
      Score.Resampling = NULL
    ),
    silent = TRUE
  )


  # LEDGER #103 fix: the burden information is the FULL sum of Phi,
  # 1'Phi 1, not its trace. The score for the burden coefficient b is
  # 1'Score and its variance is 1'Phi 1, so the estimator is
  # sum(Score) / sum(Phi) with SE 1 / sqrt(sum(Phi)); dividing by the trace
  # reported b times kappa_w = 1'Phi 1 / tr(Phi), the LD inflation of the
  # weighted score covariance (reconciled against a planted effect and
  # against the --r.corr=1 pseudo-variant path; see the ledger row). This
  # is also what SKAT's own r.corr = 1 branch uses (a <- as.matrix(sum(Phi))).
  # The p-value and Z were never affected: SE_Burden is back-derived from
  # the printed p, so it moves with BETA. (LEDGER #34, which this replaces,
  # was about diag() on a bare length-1 numeric; sum() has no such trap.)
  varSum_Phi <- sum(Phi)
  BETA_Burden <- sum(Score) / varSum_Phi
  # LEDGER #57 (the open half): the kernel's message used to die here -- the
  # error.code below is set and never read. It now rides out on the result
  # list as `error` (absent on success) and region_postprocess records it; see
  # R/SkatFailures.R. The NA rows themselves are unchanged.
  skat_error <- NULL
  if (inherits(out_SKAT_List, "try-error")) {
    Pvalue <- c(NA, NA, NA)
    error.code <- 2
    BETA_Burden <- NA
    SE_Burden <- NA
    cond <- attr(out_SKAT_List, "condition")
    skat_error <- if (is.null(cond)) trimws(as.character(out_SKAT_List)) else conditionMessage(cond)
    skat_error <- gsub("[\r\n]+", " ", trimws(skat_error))
  } else if (!any(c(0, 1) %in% out_SKAT_List$param$rho & !is.null(out_SKAT_List$p.value))) {
    # Pvalue = c(NA, NA, NA)
    Pvalue <- c(out_SKAT_List$p.value, out_SKAT_List$p.value, out_SKAT_List$p.value)
    # error.code = 3
    error.code <- 0
    # BETA_Burden = NA
    # SE_Burden = NA
    SE_Burden <- abs(BETA_Burden / qnorm((out_SKAT_List$p.value) / 2))
  } else {
    pos00 <- which(out_SKAT_List$param$rho == 0)
    pos01 <- which(out_SKAT_List$param$rho == 1)
    Pvalue <- c(
      out_SKAT_List$p.value, # SKAT-O
      out_SKAT_List$param$p.val.each[pos00], # SKAT
      out_SKAT_List$param$p.val.each[pos01]
    ) # Burden Test
    error.code <- 0
    SE_Burden <- abs(BETA_Burden / qnorm((out_SKAT_List$param$p.val.each[pos01]) / 2))
  }
  out <- list(Pvalue_SKATO = Pvalue[1], Pvalue_Burden = Pvalue[3], Pvalue_SKAT = Pvalue[2], BETA_Burden = BETA_Burden, SE_Burden = SE_Burden)
  if (!is.null(skat_error)) out$error <- skat_error
  out
}

# ACAT-V for one (annotation, max_MAF, weight set) cell: the Cauchy
# combination of the single-variant p-values of the markers in that cell,
# weighted by that weight set's marker weights (the same AnnoWeights the
# SKAT/Burden statistics use). Markers with a missing p-value are dropped, as
# in get_CCT_pvalue; markers with weight 0 contribute nothing to this set and
# are dropped too (an all-zero weight vector would make CCT's normalisation
# NaN and abort the run); NA is returned when nothing is left. p-values that
# underflowed double precision (exp(log p) == 0 for p < ~1e-308; CCT's
# small-p 1/p branch also overflows on denormals) are clamped to the smallest
# normal double, so extreme rows report a tiny nonzero ACAT-V instead of
# exactly 0. Only reached in multi-weight mode (K > 1) -- the single-weight
# schema has no ACAT-V column.
get_ACATV_pvalue <- function(pvalue, weights) {
  keep <- which(!is.na(pvalue) & weights > 0)
  if (length(keep) > 0) {
    acatvpval <- CCT(pmax(pvalue[keep], .Machine$double.xmin), weights[keep])
  } else {
    acatvpval <- NA
  }
  acatvpval
}

get_CCT_pvalue <- function(pvalue) {
  pvals <- pvalue
  notna <- which(!is.na(pvals))
  if (length(notna) > 0) {
    pvals <- pvals[!is.na(pvals)]
    cctpval <- CCT(pvals)
  } else {
    cctpval <- NA
  }
  cctpval
}
