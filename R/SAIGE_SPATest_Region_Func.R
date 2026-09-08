#############################
#
ReadGroupFile <- function(groupFile) {
  # groupFile = "./SAIGE_GENE_Test/group_multiSets_2.txt"
  group_info_list <- list()
  line <- 0
  gf <- file(groupFile, "r")
  while (TRUE) {
    line <- line + 1
    marker_group_line <- readLines(gf, n = 1)
    if (length(marker_group_line) == 0) {
      break
    }
    marker_group_line_list <- strsplit(marker_group_line, split = "[\ \t]+")[[1]]
    if (length(marker_group_line_list) < 3) {
      stop("Error, group file line:", line, "Each line should have geneID GroupName  \n")
    }

    geneID <- marker_group_line_list[1]
    type <- marker_group_line_list[2]

    # check type
    if (!(type %in% c("var", "anno"))) {
      stop("Error, group file, 2nd column should be either var or anno, line: ", line, " \n")
    }

    if (is.null(group_info_list[[geneID]])) {
      group_info_list[[geneID]] <- list(geneID = geneID)
    }
    group_info_list[[geneID]][[type]] <- marker_group_line_list[-c(1:2)]
  }
  close(gf)

  # check
  ngroup <- length(group_info_list)
  for (i in 1:ngroup) {
    n1 <- length(group_info_list[[i]])
    if (n1 < 3) {
      stop("Error: group file, each gene should have var and anno:", i, " \n")
    }

    var_g <- group_info_list[[i]][["var"]]
    anno_g <- group_info_list[[i]][["anno"]]
    n_var <- length(var_g)
    n_anno <- length(anno_g)

    if (n_var != n_anno) {
      stop("Error: group file, the length of var and anno are different:", i, n1, n2, " \n")
    }
    anno_name <- unique(anno_g)
    for (j in 1:length(anno_name)) {
      anno_name_1 <- anno_name[j]
      markerID <- var_g[anno_g == anno_name_1]
      group_info_list[[i]][[anno_name_1]] <- list(markerID = markerID)
    }
  }
  # print("group_info_list")
  # print(group_info_list)
  group_info_list
}


check_close <- function(groupFile) {
  if (isOpen(groupFile)) {
    close(groupFile)
  }
}


# check Group file format and count the number of regions
# With detail = FALSE (the default) the return value is exactly
# list(nRegions, is_weight_included), the historical contract.
# With detail = TRUE two fields are appended: nWeightSets (K, the number of
# weight lines per region, defined by region 1) and weightSetNames (character
# vector of length K; a bare 'weight' type is named "WEIGHT").
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

  # region 1 defines the format: a var line, an anno line, then K >= 0
  # consecutive weight lines; every region must follow the same layout
  gf <- file(groupFile, "r")
  marker_group_line <- readLines(gf, n = 1)
  marker_group_line <- readLines(gf, n = 1)
  is_weight_included <- FALSE
  weightSetNames <- character(0)
  repeat {
    marker_group_line <- readLines(gf, n = 1)
    if (length(marker_group_line) == 0) {
      break
    }
    marker_group_line_list <- strsplit(marker_group_line, split = "[\ \t]+")[[1]]
    if (length(marker_group_line_list) < 3) {
      break
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
  a <- 2 + nWeightSets

  gf <- file(groupFile, "r")
  line <- 0
  nregion <- 0
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

    if (nWeightSets > 0) {
      for (k in 1:nWeightSets) {
        marker_group_line_list <- strsplit(marker_group_line[2 + k], split = "[\ \t]+")[[1]]
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


get_newPhi_scaleFactor <- function(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType) {
  m1 <- sum(mu.a * g.sum)
  var1 <- sum(mu.a * (1 - mu.a) * g.sum^2)
  qinv <- -sign(q.sum - m1) * abs(q.sum - m1) + m1
  pval_noadj <- pchisq((q.sum - m1)^2 / var1, lower.tail = FALSE, df = 1, log.p = TRUE)
  isSPAConverge <- TRUE
  if (abs(q.sum - m1) / sqrt(var1) < 2) {
    p.value_burden <- pval_noadj
  } else {
    p.value_burden <- SPA_pval(mu.a, g.sum, q.sum, qinv, pval_noadj, .Machine$double.eps^0.25, TRUE, "binary", isSPAConverge)
  }
  # p.value_burden<-SPAtest:::Saddle_Prob(q.sum, mu=mu.a, g=g.sum, Cutoff=2,alpha=2.5*10^-6, log.p =TRUE)$p.value
  # cat("p.value_burden ", p.value_burden, "\n")
  re_phi <- SPA_ER_kernel_related_Phiadj_fast_new(p.new, Score, Phi, p.value_burden, regionTestType)
  re_phi
}


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
