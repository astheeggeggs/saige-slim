setSparseSigma <- function(sparseSigmaFile) {
  Check_File_Exist(sparseSigmaFile, "sparseSigmaFile")
  sparseSigma <- Matrix:::readMM(sparseSigmaFile)
  locations <- rbind(sparseSigma@i, sparseSigma@j)
  values <- sparseSigma@x
  nSubj <- dim(sparseSigma)[1]
  sigmaMatListR <- list(
    locations = locations,
    values = values,
    nSubj = nSubj
  )
  sigmaMatListR
}


setSparseSigma_new <- function(sparseGRMFile,
                               sparseGRMSampleIDFile,
                               relatednessCutoff,
                               sampleIDInModel,
                               tauVec,
                               W,
                               traitType) {
  Check_File_Exist(sparseGRMFile, "sparseGRMFile")
  Check_File_Exist(sparseGRMSampleIDFile, "sparseGRMSampleIDFile")

  sparseGRM <- Matrix:::readMM(sparseGRMFile)
  sparseGRMSampleID <- data.frame(
    data.table:::fread(
      sparseGRMSampleIDFile,
      header = FALSE,
      stringsAsFactors = FALSE,
      colClasses = c("character")
    )
  )
  colnames(sparseGRMSampleID) <- "sampleID"
  sparseGRMSampleID$IndexGRM <- c(1:nrow(sparseGRMSampleID))
  cat(
    "length(sparseGRMSampleID$IndexGRM): ",
    length(sparseGRMSampleID$IndexGRM),
    "\n"
  )
  cat("nrow(sparseGRMSampleID): ", nrow(sparseGRMSampleID), "\n")
  sampleInModel <- NULL
  sampleInModel$IID <- sampleIDInModel
  # rm(sampleIDInModel)
  sampleInModel <- data.frame(sampleInModel)
  sampleInModel$IndexInModel <- seq(1, length(sampleInModel$IID), by = 1)
  cat(
    nrow(sampleInModel),
    " samples have been used to fit the glmm null model\n"
  )
  mergeID <- merge(sampleInModel,
    sparseGRMSampleID,
    by.x = "IID",
    by.y = "sampleID"
  )
  mergeID <- mergeID[with(mergeID, order(IndexInModel)), ]
  print(dim(mergeID))
  print(head(mergeID))
  indexIDofGRM <- mergeID$IndexGRM
  sparseGRM <- sparseGRM[indexIDofGRM, indexIDofGRM]
  # LEDGER #106 fix: readMM returns the `symmetric` .mtx as a dsTMatrix whose
  # @i/@j/@x hold ONE triangle (Matrix 1.7.6 keeps it symmetric through the
  # subset), and C++ builds arma::sp_mat from exactly those triplets, so the
  # PCG used to solve against a lower-triangular Sigma. Expand to a general
  # triplet matrix first: both triangles, each diagonal entry once.
  sparseGRM <- as(as(sparseGRM, "generalMatrix"), "TsparseMatrix")
  if (length(indexIDofGRM) < nrow(sampleInModel)) {
    stop(
      nrow(sampleInModel) - length(indexIDofGRM),
      " samples were not found in the sparse GRM\n"
    )
  } else {
    print("Subsetting GRM")
  }
  removeIndex <- which(sparseGRM@x < relatednessCutoff)
  if (length(removeIndex) > 0) {
    cat(
      "Removing ",
      length(removeIndex),
      " elements in the sparse GRM < ",
      relatednessCutoff,
      ".\n"
    )
    sparseGRM@x <- sparseGRM@x[-removeIndex]
    sparseGRM@i <- sparseGRM@i[-removeIndex]
    sparseGRM@j <- sparseGRM@j[-removeIndex]
  }

  sparseGRM@x <- sparseGRM@x * tauVec[2]
  # sparseSigma = sparseGRM * tauVec[2]
  if (traitType == "binary") {
    sparseGRM@x[which(sparseGRM@i == sparseGRM@j)] <- sparseGRM@x[which(sparseGRM@i == sparseGRM@j)] + 1 / W
    # diag(sparseSigma) = W + diag(sparseSigma)
  } else if (traitType == "quantitative") {
    sparseGRM@x[which(sparseGRM@i == sparseGRM@j)] <- tauVec[1] + sparseGRM@x[which(sparseGRM@i == sparseGRM@j)]
  }
  # diag(sparseSigma) = tauVec[1] + diag(sparseSigma)

  locations <- rbind(sparseGRM@i, sparseGRM@j)
  values <- sparseGRM@x
  nSubj <- dim(sparseGRM)[1]
  sigmaMatListR <- list(
    locations = locations,
    values = values,
    nSubj = nSubj
  )
  sigmaMatListR
}


# Shared SKAT-path post-processing that runs after mainRegionInCPP in
# SAIGE.Region. Extracted verbatim from the two near-identical copies that
# used to live there (the primary pass and the fastTest re-run). The real
# differences between the two copies are parameterized by is_rerun:
#   - the re-run copy assigned annoMAFIndicatorMat inside the
#     is_single_in_groupTest block; the primary copy assigned it right after
#     that block, unconditionally
#   - the primary copy gated the core on
#     (sum(NumUltraRare_GroupVec) + sum(NumRare_GroupVec)) > 0 and
#     regionTestType != "BURDEN"; the re-run copy ran unguarded (it was
#     already nested inside those guards)
#   - only the primary copy set cctpval = 1 when the Cauchy combination is
#     skipped (single annotation and single max MAF)
# resultDF_carryover is the re-run copy's per-set carry-forward list
# (resultDFList; NULL on the first pass). The returned list carries
# pval.Region plus cctpval only when it was bound, so the caller can
# reproduce the exact frame state the inline blocks used to leave behind.
# (A scalar `resultDF` return field and the unpack that made `exists()` true
# for it were deleted 2026-09-14, C6: nothing read them.)
# LEDGER #132: the ACAT-V COMBINATION weight of a collapsed ultra-rare
# pseudo-marker, for a custom or transform set.
#
# It used to be a hardcoded 1, which is correct for the BURDEN (C++ has already
# baked w_r into the collapse max(w_r*G), LEDGER #17) but is not commensurable
# with the surviving markers' own weights when the SAME vector is reused as
# CCT's combination weights. CCT normalises by sum(weights)
# (R/CCT_modified.R), so a 1 sitting beside CADD-scale weights (0.5-10) and
# beside REVEL-scale weights (0.01-0.95) means the pseudo-marker's share of the
# combination is decided by the SCALE of the group file -- a uniform rescale,
# which is byte-exactly a no-op for Pvalue, Pvalue_Burden and Pvalue_SKAT,
# moved Pvalue_ACATV by up to 10.5%.
#
# The share is 1 / (sum(w) + 1), so as the scale falls ACAT-V tends to the
# pseudo-marker's own single-variant p and as it rises it tends to the Cauchy of
# the rare markers alone. Every scale is validly calibrated under the null; what
# the scale silently chooses is where between those two limits the statistic
# sits, which is a POWER choice the user never made.
#
# THE RULE (owner's decision, 2026-09-13): the pseudo-marker's combination
# weight is the MEAN of w_r over the markers that were actually folded into it
# IN THAT CELL.
#   - Per cell, not per set: annoMAFIndicatorMat marks a marker 2 in column jmr
#     when it was collapsed into that cell's pseudo-marker (src/Main.cpp:3025),
#     1 when it survived as a rare marker (:2955), and the pseudo-marker's own
#     row carries 1 in its own column (:3701) so it is in tempPos.
#   - NON-ZERO weights only, and that qualifier is load-bearing rather than
#     cosmetic: :3025 marks an ultra-rare marker 2 whatever its weight, but
#     :3052's `if(wUR != 0)` decides whether it is actually accumulated
#     (LEDGER #118). So the non-zero subset is exactly the set of markers the
#     collapse contains, not an approximation of it.
#   - MEAN rather than MAX: max matches the collapse operator, but it would make
#     the pseudo-marker the heaviest element of nearly every cell, and the mean
#     is what ACAT-V does with a very-rare burden component. Pending the paper
#     text.
#   - Beta sets need no equivalent: they already weight the pseudo-marker by
#     dbeta at its own collapsed MAF, over every output row.
#
# This is a KNOWING DIVERGENCE from SAIGE-QTL, recorded in the same shape as
# LEDGER #118's divergence from upstream SAIGE.
#
# APPLIED AT THE ACAT-V CALL SITE ONLY, and that separation is the whole
# subtlety. AnnoWeights is not just ACAT-V's combination vector: it also builds
# wStatVec (StatVec * w) and the weighted variance sub-matrix, i.e. the BURDEN
# and SKAT statistics. There the pseudo-marker's 1 is CORRECT, because C++ has
# already baked w_r into the collapse max(w_r*G) (#17) -- weighting its score
# again would count w_r twice. The first attempt at this fix changed the shared
# vector and did exactly that; test-step2-multiweight-behaviour.R caught it
# immediately, failing on Pvalue, Pvalue_Burden, Pvalue_SKAT, BETA_Burden and
# SE_Burden rather than on anything to do with ACAT-V. Hence two vectors:
# AnnoWeightsList for the statistics, AnnoWeightsACATVList for the combination.
#
# Fallback: a cell with no folded marker keeps 1. Its pseudo-marker row is
# zeroed in C++ (src/Main.cpp:3701 is not reached for an empty collapse) so the
# value is never read; 1 is the pre-#132 value and the safest thing to leave.
urCombinationWeight <- function(w, outList) {
  n <- outList$numofUR
  if (is.null(n) || n == 0) {
    return(numeric(0))
  }
  ind <- outList$annoMAFIndicatorMat
  if (is.null(ind) || ncol(ind) < n) {
    return(rep(1, n))
  }
  vapply(seq_len(n), function(jmr) {
    folded <- which(ind[, jmr] == 2)
    folded <- folded[folded <= length(w)]
    wf <- w[folded]
    wf <- wf[!is.na(wf) & wf != 0]
    if (length(wf) == 0) 1 else mean(wf)
  }, numeric(1))
}

region_postprocess <- function(outList,
                               mu,
                               WEIGHT,
                               traitType,
                               r.corr,
                               regionTestType,
                               isCondition,
                               annolistsub,
                               maxMAFlist,
                               regionName,
                               is_single_in_groupTest,
                               is_rerun,
                               resultDF_carryover = NULL,
                               weight_sets,
                               # LEDGER #126: the min-group-MAC gate, applied on
                               # THIS path too. The burden path enforces it in
                               # C++ (src/Main.cpp, the BURDEN emitter) on the
                               # same per-set MAC, so the two paths refuse the
                               # same cells.
                               minGroupMAC_in_BurdenTest = 0,
                               # 1-based index of the null model this call is
                               # post-processing, for the LEDGER #127 sidecar
                               model_index = 1L,
                               # FALSE on skatoCollect's drive: the parallel
                               # path runs this body twice (collect, then serve)
                               # and only the SERVE pass holds real p-values, so
                               # only it may record. The same discipline
                               # `skat_fun`'s stub answer already imposes on
                               # skatFailureRecord, made explicit here because
                               # the pooled record does not depend on an error
                               # field to stay silent.
                               record_pooled = TRUE,
                               # The SKAT-O evaluator, injected so the cell loop
                               # can be driven twice by skato_parallel(): once to
                               # COLLECT every (Score, Phi) it would test, once to
                               # SERVE the answers computed in parallel in
                               # between. Defaults to the real thing, which is
                               # the entire single-threaded path -- at
                               # nThreadsSKATO = 1 nothing below changes at all.
                               skat_fun = get_SKAT_pvalue) {
  # LEDGER #174: this model's pinned samples, for the burden rescale's saddlepoint
  pinnedSet <- if (traitType == "binary") getSpaPinnedModelInCPP(as.integer(model_index - 1L)) else NULL
  # ---- weight-set layout (jmr decode) --------------------------------------
  # The composed weight-set list arrives as `weight_sets`, a descriptor from
  # weightSetDescriptor() (R/WeightSets.R; built in SAIGE.Region, per-set
  # kind/name/params in the same custom -> beta -> No_Weights order as C++).
  # K > 1 is signalled by WEIGHT being the q0 x K matrix built in
  # SAIGE.getRegionList_new. The legacy shapes -- the c(0) sentinel of an
  # unweighted group file and the collapsed single-set numeric vector -- are
  # K = 1, and every expression below then reduces to the literal
  # single-weight code (same objects, same arithmetic, same bytes);
  # weight_sets$force (is_output_weight_column) switches a K = 1 run onto the
  # multi-weight SCHEMA while its single set keeps the K = 1 arithmetic.
  # C++ writes every (anno, max_MAF)-indexed vector weight-fastest with the
  # normative 0-based index jmr = j*q_maf*q_weight + m*q_weight + r (see
  # src/Main.cpp mainRegionInCPP). The 1-based R form of that index is
  #   jm  = (j - 1)*length(maxMAFlist) + m
  #   jmr = (jm - 1)*K + r
  # which at K = 1 (so r = 1) is jmr == jm, the legacy index.
  nWeightSets <- length(weight_sets$kind)
  ws_kind <- weight_sets$kind
  ws_customCol <- weight_sets$customCol
  ws_transformCol <- weight_sets$transformCol
  ws_betaA <- weight_sets$betaA
  ws_betaB <- weight_sets$betaB
  ws_force <- isTRUE(weight_sets$force)
  ws_names <- weight_sets$name
  isMultiWeight <- nWeightSets > 1 || ws_force
  weightNames <- if (isMultiWeight) {
    if (is.null(ws_names)) paste0("W", seq_len(nWeightSets)) else ws_names
  } else {
    NULL
  }
  # Carry-forward cache for the m > maxMAF0 branch, one slot per weight set;
  # at K = 1 this is the single legacy `resultDF` wrapped in a 1-element list
  # (exactly what the caller passed before).
  resultDFList <- if (!is.null(resultDF_carryover)) {
    resultDF_carryover
  } else {
    vector("list", nWeightSets)
  }
  pval.Region <- NULL
  if (is_single_in_groupTest) {
    # convert string p-value to log(p)
    pvalVec_new <- unlist(lapply(outList$pvalVec, convert_str_to_log))
    outList$pvalVec <- pvalVec_new ## log of p-values
    noNAIndices <- which(!is.na(outList$pvalVec))
    if (is_rerun) {
      annoMAFIndicatorMat <- outList$annoMAFIndicatorMat
    }
    if (isMultiWeight) {
      # per-set weights (design D2), by set kind:
      #  custom: the surviving rare markers carry set r's group-file weights
      #    and every collapsed ultra-rare pseudo-marker carries 1, because
      #    C++ already baked w_r into the ultra-rare collapse max(w_r*G).
      #  beta: dbeta(MAF, a_r, b_r) over EVERY output row -- the ultra-rare
      #    collapse is unweighted for beta sets, and C++ scaled only its
      #    genoSumMat contribution by pdf(beta_r, MAF_collapsed), so the
      #    pseudo-marker weights like any marker at its collapsed MAF (the
      #    literal legacy K = 1 flow). The implicit single set of an unweighted
      #    default run is one of these, holding the FIRST weights.beta pair --
      #    the pair C++ used for genoSumMat and the ultra-rare collapse
      #    (LEDGER #19 FIXED: it used to be a kind of its own hardcoded to
      #    dbeta(MAF, 1, 25) whatever --weights.beta said, a C++/R hybrid);
      #    reachable here only at K = 1 with is_output_weight_column.
      #  noweight: all ones.
      #  transform (FlexRV): custom-set semantics over the q0 x Kt weight
      #    matrix C++ computed per marker and exported as transformWeightMat
      #    (NaN on a marker that failed QC; such rows are dropped with the
      #    NA p-values below, exactly as a custom set's are)
      # outList$numofUR is q_anno*q_maf*K, so each vector is q0 + numofUR
      # long -- one entry per row of the marker-level output, as at K = 1
      # (outList$MAFVec has that same length). Rows of other sets are never
      # read for set r: each cell subsets by annoMAFIndicatorMat[, jmr].
      AnnoWeightsList <- lapply(seq_len(nWeightSets), function(r) {
        switch(ws_kind[r],
          custom = {
            w <- if (is.matrix(WEIGHT)) WEIGHT[, ws_customCol[r]] else as.numeric(WEIGHT)
            c(w, rep(1, outList$numofUR))
          },
          beta = dbeta(outList$MAFVec, ws_betaA[r], ws_betaB[r]),
          noweight = rep(1, length(outList$MAFVec)),
          transform = c(outList$transformWeightMat[, ws_transformCol[r]],
                        rep(1, outList$numofUR)),
          stop("region_postprocess: unknown weight-set kind '", ws_kind[r], "'\n")
        )
      })
      # LEDGER #132: the ACAT-V COMBINATION weights are a SEPARATE vector from
      # the statistics weights above. See urCombinationWeight for why the
      # pseudo-marker's entry differs between the two, and why getting that
      # wrong double-counts w_r in the burden.
      AnnoWeightsACATVList <- lapply(seq_len(nWeightSets), function(r) {
        aw <- AnnoWeightsList[[r]]
        if (!(ws_kind[r] %in% c("custom", "transform"))) {
          return(aw)
        }
        q0 <- length(aw) - outList$numofUR
        if (q0 < 0) {
          return(aw)
        }
        aw[seq_len(outList$numofUR) + q0] <-
          urCombinationWeight(aw[seq_len(q0)], outList)
        aw
      })
    } else if (sum(WEIGHT) > 0) {
      AnnoWeights <- c(WEIGHT, rep(1, outList$numofUR))
    } else if (identical(ws_kind, "transform")) {
      # LEDGER #119 FIXED. A FlexRV spec that selects ONE transform set
      # (e.g. unw x unw) is K = 1 on the legacy schema, with WEIGHT == 0
      # because the group file has no weight line; before this branch the
      # single set fell through to the implicit-beta arm below, whose guard
      # stopped the run. The set is a custom set here: its exported column,
      # 1 on every collapsed pseudo-marker.
      AnnoWeights <- c(outList$transformWeightMat[, 1], rep(1, outList$numofUR))
    }
  }
  if (!is_rerun) {
    annoMAFIndicatorMat <- outList$annoMAFIndicatorMat
  }

  run_core <- is_rerun
  if (!run_core) {
    if ((sum(outList$NumUltraRare_GroupVec) + sum(outList$NumRare_GroupVec)) > 0) {
      if (regionTestType != "BURDEN") {
        run_core <- TRUE
      }
    }
  }
  if (run_core) {
    if (traitType == "binary") {
      outList$gyVec <- outList$gyVec[noNAIndices]
    }

    if (isCondition) {
      if (nWeightSets == 1) {
        # the C++ baked adjustment exists only at q_weight == 1 (design D7);
        # at K > 1 mainRegionInCPP exports the G1-raw matrix plus the
        # VarInvMat_G2_cond_w0/Tstat_G2_cond_w0 pieces instead, and the
        # per-set adjustment is formed per cell below
        # LEDGER #34: drop = FALSE, else a one-marker region collapses this
        # to a scalar and the [tempPos, tempPos, drop = F] below raises
        # "incorrect number of dimensions".
        outList$VarMatAdjCond <- outList$VarMatAdjCond[noNAIndices, noNAIndices, drop = FALSE]
        outList$TstatAdjCond <- outList$TstatAdjCond[noNAIndices]
      }
      outList$G1tilde_P_G2tilde_Weighted_Mat <- outList$G1tilde_P_G2tilde_Weighted_Mat[noNAIndices, ,
        drop =
          FALSE
      ]
      weightMat_G2_G2 <- outList$G2_Weight_cond %*% t(outList$G2_Weight_cond)
    }

    ### Get annotation maf indicators

    ### 3. Adjust for saddlepoint approximation
    # notNAindice = which(!is.na(outList$TstatVec_flip))
    # print("notNAindice")
    # print(notNAindice)
    StatVec <- outList$TstatVec_flip[noNAIndices]
    VarSVec <- diag(outList$VarMat)
    VarSVec <- VarSVec[!is.na(VarSVec)]
    adjPVec <- outList$pvalVec[!is.na(outList$pvalVec)] ## adjPVec contains log(p)

    # varTestedIndices = which(apply(annoMAFIndicatorMat, 1, isContainValue, val=1))
    # print("varTestedIndices")
    # print(varTestedIndices)
    # varTestedIndices = which(rowSums(annoMAFIndicatorMat) > 0)
    # annoMAFIndicatorMat = annoMAFIndicatorMat[varTestedIndices, , drop=F]
    # MAFVec = outList$MAFVec[varTestedIndices]
    annoMAFIndicatorMat <- annoMAFIndicatorMat[noNAIndices, , drop = FALSE]
    MAFVec <- outList$MAFVec[noNAIndices]
    # AnnoWeights = dbeta(MAFVec,1,25)

    if (isMultiWeight) {
      AnnoWeightsList <- lapply(AnnoWeightsList, function(w) w[noNAIndices])
      AnnoWeightsACATVList <- lapply(AnnoWeightsACATVList, function(w) w[noNAIndices])
    } else if (sum(WEIGHT) > 0 || identical(ws_kind, "transform")) {
      # AnnoWeights = AnnoWeights[varTestedIndices]
      # (a single FlexRV transform set is a custom set here, LEDGER #119)
      AnnoWeights <- AnnoWeights[noNAIndices]
    } else {
      # LEDGER #19 FIXED. On the legacy single-set unweighted path the one
      # set is the implicit beta set, whose pair is the first weights.beta pair
      # (as C++'s already was; this line was beta(1, 25) whatever the flag
      # said). At the default pair the output bytes are unchanged.
      stopifnot(identical(ws_kind, "beta"))
      AnnoWeights <- dbeta(MAFVec, ws_betaA[1], ws_betaB[1])
    }

    if (isMultiWeight) {
      # P1Mat/P2Mat rows -- and hence VarMat -- are unweighted, so each set
      # applies its own weights here. Only the (cheap) score vectors are
      # materialised per set; the weighted variance sub-matrix is formed inside
      # the loop on the cell's markers, since
      #  (VarMat * w w')[pos, pos] == VarMat[pos, pos] * w[pos] w[pos]'
      # element for element (same two multiplications in the same order).
      wStatVecList <- lapply(AnnoWeightsList, function(w) StatVec * w)
    } else {
      weightMat <- AnnoWeights %*% t(AnnoWeights)


      if (isCondition) {
        weightMat_G1_G2 <- AnnoWeights %*% t(outList$G2_Weight_cond)
      }
      wStatVec <- StatVec * AnnoWeights

      wadjVarSMat <- outList$VarMat * weightMat

      if (isCondition) {
        wStatVec_cond <- wStatVec - outList$TstatAdjCond
        wadjVarSMat_cond <- wadjVarSMat - outList$VarMatAdjCond
      }
    }

    # gc()

    annoMAFIndVec <- c()
    # LEDGER #127: cells this region refused on the gate, counted where the
    # decision is made. Local to the call, so a region that takes the fastTest
    # escalation pass records its two passes separately.
    nRefusedMAC <- 0L
    for (j in 1:length(annolistsub)) {
      AnnoName <- annolistsub[j]
      maxMAF0 <- outList$q_maf_for_annoVec[j]
      # one flag per weight set; at K = 1 this is the legacy scalar isPolyRegion
      isPolyRegionVec <- rep(TRUE, nWeightSets)
      for (m in 1:length(maxMAFlist)) {
        jm <- (j - 1) * (length(maxMAFlist)) + m
        maxMAFName <- maxMAFlist[m]
        # rows are (anno x max_MAF x weight set), weight-innermost
        for (r in 1:nWeightSets) {
          jmr <- (jm - 1) * nWeightSets + r # K = 1, r = 1 -> jmr == jm
          if (isMultiWeight) {
            AnnoWeights <- AnnoWeightsList[[r]]
            # LEDGER #132: ACAT-V combines with its OWN weight vector
            AnnoWeightsACATV <- AnnoWeightsACATVList[[r]]
            wStatVec <- wStatVecList[[r]]
          }
          if (m <= maxMAF0) {
            tempPos <- which(annoMAFIndicatorMat[, jmr] == 1)
            # LEDGER #117. A weight of exactly 0 means the marker is not in this
            # set, so it must not reach the kernel. Left in, it zeroes its row
            # and column of the weighted Phi, and the binary-trait Phi
            # adjustment then computes Score^2 / diag(Phi) = 0/0: the NaN
            # propagates and the WHOLE CELL comes back NA for Pvalue,
            # Pvalue_Burden and Pvalue_SKAT through the #57 failure path. The
            # harm is not one NA row -- the failed set drops silently out of the
            # region's Cauchy combination, which then reports a surviving set's
            # p-value verbatim (#127 route (a)).
            # The gate is an identity, not a golden: zeroing a marker must equal
            # DELETING it from the group file, which the BURDEN path already
            # satisfies to every printed digit. Applied for both traits --
            # quantitative never enters the failing kernel, but "weight 0 means
            # not in the set" is the same statement there, and it matches the
            # per-set counters of #118, which are gated on the same predicate.
            tempPos <- tempPos[AnnoWeights[tempPos] != 0]
            # LEDGER #126. The C++ BURDEN emitter refuses a cell on
            #   indexNonZeroVec_arma.n_elem > 0 && MAC >= min_gourpmac_for_burdenonly
            # and this is the same predicate on the same quantity: length(tempPos)
            # is the cell's marker support and MAC_GroupVec[jmr] is the MAC C++
            # accumulated for it. A refused cell emits no row and clears the
            # carry-forward flag, exactly as `c.tested == false` does on the
            # burden path.
            cellMAC <- outList$MAC_GroupVec[jmr]
            gateOK <- length(cellMAC) == 1 && !is.na(cellMAC) &&
              cellMAC >= minGroupMAC_in_BurdenTest
            if (length(tempPos) > 0 && gateOK) {
              isPolyRegionVec[r] <- TRUE
              annoMAFIndVec <- c(annoMAFIndVec, jmr)
              Phi <- if (isMultiWeight) {
                outList$VarMat[tempPos, tempPos, drop = FALSE] *
                  (AnnoWeights[tempPos] %*% t(AnnoWeights[tempPos]))
              } else {
                wadjVarSMat[tempPos, tempPos, drop = FALSE]
              }
              Score <- wStatVec[tempPos]
              if (traitType == "binary") {
                p.new <- adjPVec[tempPos]
                g.sum <- outList$genoSumMat[, jmr]
                q.sum <- sum(outList$gyVec[tempPos] * AnnoWeights[tempPos])
                mu.a <- mu

                re_phi <- get_newPhi_scaleFactor_traitType(
                  q.sum,
                  mu.a,
                  g.sum,
                  p.new,
                  Score,
                  Phi,
                  regionTestType,
                  traitType,
                  pinned = pinnedSet
                )
                Phi <- re_phi$val
                # LEDGER #174: recorded on the pass that holds real p-values
                if (record_pooled) {
                  spaPinnedRecord("region_rescale", regionName, AnnoName, maxMAFName,
                                  if (isMultiWeight) weightNames[r] else "", model_index, re_phi$pinned)
                }
              }
              groupOutList <- skat_fun(Score, Phi, r.corr, regionTestType)
              # LEDGER #158: which path produced this cell's Pvalue (sidecar)
              skatoMethodRecord(regionName, AnnoName, maxMAFName,
                                if (isMultiWeight) weightNames[r] else NA_character_,
                                "marginal", groupOutList$skato_method,
                                groupOutList$skato_screen_p)
              if (!is.null(groupOutList$error)) {
                # LEDGER #57: the kernel failed and this cell's three p-values
                # are NA. Count it here, where the labels are in scope.
                skatFailureRecord(regionName, AnnoName, maxMAFName,
                                  if (isMultiWeight) weightNames[r] else NA_character_,
                                  "marginal", groupOutList$error)
              }
              if (isMultiWeight) {
                # additive schema (design D3): Pvalue keeps its SKAT-O meaning,
                # the Weight label sits right after max_MAF, and the two omnibus
                # columns sit right after Pvalue_SKAT. None of them exist at K = 1.
                # ACAT-V combines this cell's single-variant p-values with set r's
                # weights; ACAT-O is the Cauchy combination of (ACAT-V, Burden,
                # SKAT) for this row.
                # LEDGER #132: AnnoWeightsACATV, not AnnoWeights. The two differ
                # only in the collapsed pseudo-marker's entry -- 1 in the
                # statistics vector (C++ already baked w_r into the collapse,
                # #17, so weighting the score again would double-count it) and
                # the mean of the folded markers' non-zero weights here, where
                # a bare 1 is not commensurable with the surviving markers'
                # own weights and made ACAT-V depend on the group file's SCALE.
                Pvalue_ACATV <- get_ACATV_pvalue(exp(adjPVec[tempPos]),
                                                 AnnoWeightsACATV[tempPos])
                Pvalue_ACATO <- get_CCT_pvalue(c(
                  Pvalue_ACATV,
                  groupOutList$Pvalue_Burden,
                  groupOutList$Pvalue_SKAT
                ))
                resultDF <- data.frame(
                  Region = regionName,
                  Group = AnnoName,
                  max_MAF = maxMAFName,
                  Weight = weightNames[r],
                  Pvalue = groupOutList$Pvalue_SKATO,
                  Pvalue_Burden = groupOutList$Pvalue_Burden,
                  Pvalue_SKAT = groupOutList$Pvalue_SKAT,
                  Pvalue_ACATV = Pvalue_ACATV,
                  Pvalue_ACATO = Pvalue_ACATO,
                  BETA_Burden = groupOutList$BETA_Burden,
                  SE_Burden = groupOutList$SE_Burden
                )
              } else {
                resultDF <- data.frame(
                  Region = regionName,
                  Group = AnnoName,
                  max_MAF = maxMAFName,
                  Pvalue = groupOutList$Pvalue_SKATO,
                  Pvalue_Burden = groupOutList$Pvalue_Burden,
                  Pvalue_SKAT = groupOutList$Pvalue_SKAT,
                  BETA_Burden = groupOutList$BETA_Burden,
                  SE_Burden = groupOutList$SE_Burden
                )
              }

              if (isCondition) {
                if (traitType == "binary") {
                  G1tilde_P_G2tilde_Mat_sub <- outList$G1tilde_P_G2tilde_Weighted_Mat[tempPos, , drop = FALSE]
                  if (nWeightSets > 1) {
                    # design D7: at K > 1 the exported matrix is G1-raw, so
                    # apply set r's per-marker weights to its rows -- exactly
                    # the w0 the q_weight == 1 path bakes in C++ (custom sets:
                    # w_r for real markers and 1 for UR pseudo-markers, whose
                    # collapse already carries w_r; beta sets: dbeta(MAF))
                    G1tilde_P_G2tilde_Mat_sub <- G1tilde_P_G2tilde_Mat_sub * AnnoWeights[tempPos]
                  }
                  G1tilde_P_G2tilde_Mat_scaled <- t(t((
                    G1tilde_P_G2tilde_Mat_sub
                  ) * sqrt(as.vector(
                    re_phi$scaleFactor
                  ))) * sqrt(
                    as.vector(outList$scalefactor_G2_cond)
                  ))
                  # t(t(b * sqrt(a1)) * sqrt(a2))
                  adjCondTemp <- G1tilde_P_G2tilde_Mat_scaled %*% outList$VarInvMat_G2_cond_scaled
                  VarMatAdjCond <- adjCondTemp %*% t(G1tilde_P_G2tilde_Mat_scaled)
                  TstatAdjCond <- adjCondTemp %*% (outList$Tstat_G2_cond * outList$G2_Weight_cond)
                  # LEDGER #33 fix: subtract the FULL adjustment matrix, as the
                  # quantitative arms at :452 and :308 always did. The diag()
                  # here came from upstream 7c7e6da5 ("only for burden tests"),
                  # where re_phi$val is as.vector(); on this arm it is a q x q
                  # as.matrix(), so diag() EXTRACTED the diagonal and R recycled
                  # it column-major, subtracting d[i] across the whole of row i
                  # and leaving an asymmetric covariance. run_core (:230-237)
                  # keeps BURDEN out of this block entirely, so the vector case
                  # it was written for never reaches here.
                  Phi_cond <- re_phi$val - VarMatAdjCond
                  Score_cond <- Score - TstatAdjCond
                } else if (isMultiWeight) {
                  # quantitative trait under the multi-weight schema. The full
                  # wStatVec_cond/wadjVarSMat_cond matrices of the legacy arm
                  # below are per-set here, so the conditional adjustment is
                  # formed on the cell's rows only; subsetting the difference
                  # and differencing the subsets are the same elementwise
                  # operations.
                  if (nWeightSets > 1) {
                    # per-set adjustment from the G1-raw exports (design D7):
                    # the literal q_weight == 1 C++ arithmetic with set r's w0
                    # applied here instead of baked into the matrix
                    G1w_sub <- outList$G1tilde_P_G2tilde_Weighted_Mat[tempPos, , drop = FALSE] *
                      AnnoWeights[tempPos]
                    adjCondTemp <- G1w_sub %*% outList$VarInvMat_G2_cond_w0
                    VarMatAdjCond <- adjCondTemp %*% t(G1w_sub)
                    TstatAdjCond <- adjCondTemp %*% outList$Tstat_G2_cond_w0
                  } else {
                    # K = 1 with the forced Weight column: C++ baked the single
                    # set's w0 as always, so its adjustment applies unchanged
                    TstatAdjCond <- outList$TstatAdjCond[tempPos]
                    VarMatAdjCond <- outList$VarMatAdjCond[tempPos, tempPos, drop = FALSE]
                  }
                  Score_cond <- Score - TstatAdjCond
                  Phi_cond <- Phi - VarMatAdjCond
                } else {
                  Score_cond <- wStatVec_cond[tempPos]
                  # LEDGER #34: drop = FALSE, as the three sibling subsets
                  # already had. A single-marker cell otherwise degrades this
                  # to a bare scalar and get_SKAT_pvalue's trace sees it.
                  Phi_cond <- wadjVarSMat_cond[tempPos, tempPos, drop = FALSE]
                }

                groupOutList_cond <- skat_fun(
                  Score_cond,
                  Phi_cond,
                  r.corr,
                  regionTestType
                )
                skatoMethodRecord(regionName, AnnoName, maxMAFName,
                                  if (isMultiWeight) weightNames[r] else NA_character_,
                                  "conditional", groupOutList_cond$skato_method,
                                  groupOutList_cond$skato_screen_p)
                if (!is.null(groupOutList_cond$error)) {
                  skatFailureRecord(regionName, AnnoName, maxMAFName,
                                    if (isMultiWeight) weightNames[r] else NA_character_,
                                    "conditional", groupOutList_cond$error)
                }


                resultDF$Pvalue_cond <- groupOutList_cond$Pvalue_SKATO
                resultDF$Pvalue_Burden_cond <- groupOutList_cond$Pvalue_Burden
                resultDF$Pvalue_SKAT_cond <- groupOutList_cond$Pvalue_SKAT
                resultDF$BETA_Burden_cond <- groupOutList_cond$BETA_Burden
                resultDF$SE_Burden_cond <- groupOutList_cond$SE_Burden
              } # if(isCondition){
              pval.Region <- rbind.data.frame(pval.Region, resultDF)
              resultDFList[[r]] <- resultDF
            } else {
              # if(length(tempPos) > 0){
              if (length(tempPos) > 0) {
                # LEDGER #126/#127: refused by the gate rather than empty.
                nRefusedMAC <- nRefusedMAC + 1L
              }
              isPolyRegionVec[r] <- FALSE
            }
          } else {
            # if(m <= maxMAF0){
            if (isPolyRegionVec[r]) {
              annoMAFIndVec <- c(annoMAFIndVec, jmr)
              # this max_MAF cutoff adds no marker, so set r's last tested row is
              # re-emitted (its Weight label is already set r's own)
              resultDF <- resultDFList[[r]]
              resultDF$Region <- regionName
              resultDF$Group <- AnnoName
              resultDF$max_MAF <- maxMAFName
              pval.Region <- rbind.data.frame(pval.Region, resultDF)
              resultDFList[[r]] <- resultDF
              # pval.Region_temp = rbind.data.frame(pval.Region_temp, resultDF)
            }
          } # if(m <= maxMAF0){
        } # for(r in 1:nWeightSets){
      } # for(m in 1:length(maxMAFlist)){
    } # for(j in 1:length(annolist)){


    # (per-iteration gc() removed in Phase 4: it dominated the profile; R collects on demand)
    # }


    # if(regionTestType != "BURDEN"){

    if (length(annoMAFIndVec) > 0) {
      pval.Region$MAC <- outList$MAC_GroupVec[annoMAFIndVec]
      if (traitType == "binary") {
        pval.Region$MAC_case <- outList$MACCase_GroupVec[annoMAFIndVec]
        pval.Region$MAC_control <- outList$MACCtrl_GroupVec[annoMAFIndVec]
      }
      pval.Region$Number_rare <- outList$NumRare_GroupVec[annoMAFIndVec]
      pval.Region$Number_ultra_rare <- outList$NumUltraRare_GroupVec[annoMAFIndVec]
    }

    # LEDGER #127: this region's pooled-row provenance, taken BEFORE the pooled
    # row is appended (so `pval.Region` still holds only the cell rows). The
    # pooled statistic on this path is `Pvalue` -- SKAT-O's own -- and the rule
    # named here is the one CCT takes on exactly this vector.
    pooledCellRows <- if (is.null(pval.Region)) {
      numeric(0)
    } else if ("Pvalue" %in% colnames(pval.Region)) {
      as.numeric(pval.Region$Pvalue)
    } else {
      as.numeric(pval.Region$Pvalue_Burden)
    }
    pooledP <- NA_real_

    if (length(annolistsub) > 1 | length(maxMAFlist) > 1 | isMultiWeight) {
      if (is.null(pval.Region)) {
        # Degenerate corner: run_core was TRUE (the region has rare/ultra-rare
        # carriers) but no (anno, max_MAF, weight set) cell emitted a row --
        # e.g. every flagged marker's p-value came back NA, a weighted URV
        # collapse produced an all-zero pseudo-genotype, or (since LEDGER #126)
        # --minGroupMAC_in_BurdenTest refused every cell.
        #
        # LEDGER #128 FIXED: NO POOLED ROW IS BUILT HERE. The CCT values are
        # still computed, because `cctpval` is the fastTest escalation signal and
        # is returned, but nothing is appended to pval.Region.
        #
        # What this replaces: the legacy code built the row as a BARE UNNAMED
        # character vector and did `rbind(NULL, cctVec)`, giving a 1-row matrix
        # with no column names. ncol(NULL) is NULL, so the named construction in
        # the else-arm cannot run here -- that is why the legacy form existed.
        # The run then continued *inside* this region and died one region LATER,
        # in regionFinish's cross-region accumulator, where rbind dispatches to
        # rbind.data.frame and match.names stops with "names do not match
        # previous names". The old comment here claimed the run "continues
        # exactly as at HEAD", which was true within the region and false a
        # region later.
        #
        # Why not name the columns instead: the column set depends on
        # regionTestType, isCondition, isMultiWeight and traitType, so building
        # it without a pval.Region to read colnames() from means RE-DERIVING the
        # schema, and one wrong position silently puts an NA in the wrong column.
        # Writing no row needs no schema and matches what this path ALREADY does
        # for a region with one cell (the `cctpval <- 1` else-arm below, pinned
        # by test-step2-region-allrefused.R): on the SKAT-O path a region with no
        # tested cell contributes no row, and LEDGER #127's .pooledTests.txt
        # sidecar carries it with Rows = 0.
        #
        # LEDGER #116 (FIXED 2026-09-13) then harmonised this with the BURDEN
        # path: the NA Cauchy row for a region with no tested cell is built
        # AFTER this if/else, from R/RegionSchema.R, for every configuration.
        cctpval_Burden <- get_CCT_pvalue(pval.Region$Pvalue_Burden)
        if (regionTestType != "BURDEN") {
          cctpval <- get_CCT_pvalue(pval.Region$Pvalue)
          cctpval_SKAT <- get_CCT_pvalue(pval.Region$Pvalue_SKAT)
        }
        if (isCondition) {
          # Same overwrite-with-the-conditional-value order as the else-arm
          # (upstream behaviour, preserved): cctpval ends up conditional.
          cctpval_Burden <- get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)
          if (regionTestType != "BURDEN") {
            cctpval <- get_CCT_pvalue(pval.Region$Pvalue_cond)
            cctpval_SKAT <- get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)
          }
        }
        pooledP <- if (regionTestType != "BURDEN") cctpval else cctpval_Burden
      } else {
        # Exactly ONE pooled Cauchy row per region, combining every
        # (anno x max_MAF x weight set) row of the region. NOTE on the names:
        # rbind.data.frame IGNORES the names of an atomic-vector row and binds
        # purely positionally -- the names on cctRow document intent, they are
        # NOT load-bearing. The bytes are correct because names(cctRow) is
        # assigned from colnames(pval.Region) IN ORDER, so list order == column
        # order; never reorder these fields and never add one via $<- (a new
        # name appends at the END of the list, silently shifting every later
        # column). unlist() reproduces the legacy c() coercion (every field
        # becomes character, which promotes the whole frame -- the golden bytes
        # depend on that).
        cctRow <- as.list(rep(NA, ncol(pval.Region)))
        names(cctRow) <- colnames(pval.Region)
        cctRow$Region <- regionName
        cctRow$Group <- "Cauchy"
        cctpval_Burden <- get_CCT_pvalue(pval.Region$Pvalue_Burden)
        cctRow$Pvalue_Burden <- cctpval_Burden
        if (regionTestType != "BURDEN") {
          # print(pval.Region$Pvalue)
          cctpval <- get_CCT_pvalue(pval.Region$Pvalue)
          cctpval_SKAT <- get_CCT_pvalue(pval.Region$Pvalue_SKAT)
          cctRow$Pvalue <- cctpval
          cctRow$Pvalue_SKAT <- cctpval_SKAT
          if (isMultiWeight) {
            cctRow$Pvalue_ACATV <- get_CCT_pvalue(pval.Region$Pvalue_ACATV)
            cctRow$Pvalue_ACATO <- get_CCT_pvalue(pval.Region$Pvalue_ACATO)
          }
        }
        if (isCondition) {
          # NOTE (upstream behavior, preserved): these overwrite the unconditional
          # values, so the cctpval returned for the fastTest decision is the
          # conditional one
          cctpval_Burden <- get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)
          cctRow$Pvalue_Burden_cond <- cctpval_Burden

          if (regionTestType != "BURDEN") {
            cctpval <- get_CCT_pvalue(pval.Region$Pvalue_cond)
            cctpval_SKAT <- get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)
            cctRow$Pvalue_cond <- cctpval
            cctRow$Pvalue_SKAT_cond <- cctpval_SKAT
          }
        }

        pooledP <- if (regionTestType != "BURDEN") cctRow$Pvalue else cctRow$Pvalue_Burden
        pval.Region <- rbind(pval.Region, unlist(cctRow))
      } # if (is.null(pval.Region)) ... else
    } else {
      if (!is_rerun) {
        cctpval <- 1
      }
    }
    # LEDGER #116 FIXED: a region with NO tested cell -- every cell refused by
    # --minGroupMAC_in_BurdenTest, every flagged marker's p NA, or a 1 x 1
    # configuration whose one cell was refused -- contributes the one row the
    # BURDEN writer has always emitted for it, "Region Cauchy NA ...", instead
    # of no row at all. The schema comes from R/RegionSchema.R (the "why not
    # name the columns" objection above is answered there: the schema is now
    # stated once and asserted against every written row). The #127 sidecar
    # records the region with Rows = 0 below, as before.
    if (is.null(pval.Region) && regionTestType != "BURDEN") {
      pval.Region <- regionEmptyPooledRow(regionName, traitType, isCondition, isMultiWeight)
    }
    # LEDGER #127. One record per (region, model, pass). A region with one cell
    # writes no pooled row at all; it is still recorded, because a reader cannot
    # otherwise see that the region produced nothing.
    if (record_pooled) {
    pooledRecord(
      regionName, model_index,
      rows = length(pooledCellRows),
      refused = nRefusedMAC,
      combined = sum(!is.na(pooledCellRows)),
      rule = pooledRule(pooledCellRows),
      pvalue = if (length(pooledP) == 1) as.numeric(pooledP) else NA_real_
    )
    }
  }

  ret <- list(pval.Region = pval.Region)
  if (exists("cctpval", inherits = FALSE)) {
    ret$cctpval <- cctpval
  }
  if (any(!vapply(resultDFList, is.null, logical(1)))) {
    ret$resultDFList <- resultDFList
  }
  ret
}


# DosageCutoff_for_UltraRarePresence,
# method_to_CollapseUltraRare,
## working
SAIGE.Region <- function(mu,
                         OutputFile,
                         MACCutoff_to_CollapseUltraRare,
                         groupFile,
                         annolist,
                         maxMAFlist,
                         markers_per_chunk_in_groupTest,
                         genoType,
                         markerInfo,
                         traitType,
                         isImputation,
                         isCondition,
                         weight_cond,
                         groups_per_chunk,
                         r.corr,
                         isOverWriteOutput,
                         is_single_in_groupTest,
                         is_no_weight_in_groupTest,
                         is_output_markerList_in_groupTest,
                         chrom,
                         is_fastTest,
                         pval_cutoff_for_fastTest,
                         is_output_moreDetails,
                         weights_beta_mat = matrix(c(1, 25), nrow = 1, ncol = 2),
                         is_output_weight_column = FALSE,
                         # FlexRV (F3): NULL, or the list SPAGMMATtest builds
                         # from --flexRV_score/--flexRV_scoreTransforms/
                         # --flexRV_MAFTransforms/--flexRV_lofAnno/--flexRV_maxMAF
                         # (fields score, scoreTransforms, MAFTransforms,
                         # lofAnno, maxMAF). The transform sets it describes are
                         # composed right after the group file's custom sets.
                         flexRV = NULL,
                         # LEDGER #126: the min-group-MAC gate. C++ applies it on
                         # the BURDEN path; region_postprocess applies the same
                         # predicate on the SKAT-O path, which also reports
                         # Pvalue_Burden.
                         minGroupMAC_in_BurdenTest = 0,
                         # Workers for the SKAT-O post-processing only. The C++
                         # region loop stays single-threaded; this parallelises
                         # the 45% of a group run that Rprof puts in
                         # SKAT_davies. 1 means the serial path runs untouched.
                         nThreadsSKATO = 1L,
                         nThreads = 1L) {
  # Step 4 (2026-09-06): several null models per region. `mu` is a LIST of
  # fitted-value vectors (one per model), `OutputFile` and `traitType` vectors
  # of the same length; a single model arrives as the legacy scalars and is
  # wrapped here, so everything below is the one-model flow at nModels = 1.
  # One C++ call per region tests every model (they share a sample set, checked
  # by the caller); each model keeps its own result file, .index, marker list
  # and single-in-group file, byte-identical to its own single-model run.
  if (!is.list(mu)) {
    mu <- list(mu)
  }
  nModels <- length(mu)
  OutputFiles <- OutputFile
  traitTypes <- traitType
  if (length(OutputFiles) != nModels || length(traitTypes) != nModels) {
    stop("SAIGE.Region: mu, OutputFile and traitType must describe the same models\n")
  }
  OutputFile <- OutputFiles[1]
  traitType <- traitTypes[1]
  modelIdx <- seq_len(nModels) - 1L
  OutputFileIndexes <- paste0(OutputFiles, ".index")
  OutputFileIndex <- NULL
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- OutputFileIndexes[1]
  }

  # The worker pool outlives individual regions but must not outlive the scan.
  on.exit(stopSkatoCluster(), add = TRUE)

  # 2026-09-06: --nThreads also threads the C++ region marker loop (one marker
  # per worker inside each block, the accumulation serial), byte-identical to
  # one thread. A reader that cannot clone (the PLINK pread fallback) drops to
  # one thread with a printed reason, exactly as the marker path does. The
  # thread count is PRINTED so a knob that never arrived cannot pass as green
  # (LEDGER #30); the standing test asserts this line.
  nThreads <- as.integer(nThreads)
  if (nThreads > 1 && !isGenoReaderCloneableInCPP()) {
    warning("--nThreads=", nThreads, " requested but the genotype reader cannot be ",
            "cloned (PLINK .bed pread fallback); the region marker loop runs on 1 thread.")
    cat("Region marker loop runs on 1 thread: the genotype reader cannot be cloned.\n")
    nThreads <- 1L
  }
  cat("Region marker loop will use ", nThreads, " in-process threads.\n", sep = "")
  skatFailureReset() # LEDGER #57: per-run accounting, R/SkatFailures.R
  resetSpaFallbackRecordsInCPP() # LEDGER #164: per-run accounting, R/SpaFallbacks.R
  resetExactByWeightRecordsInCPP() # LEDGER #172: per-run accounting, R/ExactByWeight.R
  resetSpaPinnedRecordsInCPP() # LEDGER #174: per-run accounting, R/SpaPinned.R
  resetStretchGateRecordsInCPP() # LEDGER #176: per-run accounting, R/StretchGate.R
  skatoDedupeReset() # identical cells evaluated once: per-run counts, R/SkatoParallel.R
  burdenMemoReset() # identical BURDEN cells evaluated once: per-run counts, R/BurdenMemo.R
  pooledReset() # LEDGER #127: pooled-row provenance, R/PooledRow.R

  # LEDGER #25 FIXED: the index is written and validated with the SAME nEachChunk.
  # It records completed REGIONS one per line, so the chunk size it states is 1
  # on both sides (the writer used numberRegionsInChunk, still 0 when captured).
  # checkOutputFile is in 'Util.R'; one check per model, and the models' index
  # files must agree about where to restart (R/MultiModel.R).
  outChecks <- lapply(seq_len(nModels), function(k) {
    checkOutputFile(OutputFiles[k], OutputFileIndexes[k], "Region", 1, isOverWriteOutput)
  })
  outList <- outChecks[[1]]
  if (nModels > 1) {
    restartChunkForModels(outChecks, OutputFileIndexes)
    if (any(vapply(outChecks, function(x) !identical(x$Start, outList$Start) || !identical(x$End, outList$End), logical(1)))) {
      stop("The models' output files disagree about whether the analysis has started or finished; delete them and re-run, or set --is_overwrite_output=TRUE\n")
    }
  }

  # Phase G1 slice 5: the cross-region accumulator layer lives in ONE
  # environment, so the closures below write scan$x rather than `<<-` over
  # seven enclosing names (and no list-returning stage forces a rebind).
  scan <- new.env(parent = emptyenv())
  scan$indexChunk <- outList$indexChunk
  scan$Start <- outList$Start
  scan$End <- outList$End

  cat("Start ", scan$Start, "\n")
  cat("End ", scan$End, "\n")


  if (scan$End) {
    message <- paste0(
      "The analysis has been completed in earlier analysis. Results are saved in '",
      paste(OutputFiles, collapse = "', '"),
      "'. ",
      "If you want to change parameters and restart the analysis, please use another 'OutputFile'."
    )
    return(message)
  }

  isappend <- FALSE
  if (!scan$Start) {
    isappend <- TRUE
    # LEDGER #140 FIXED: S9 (#47) carried to the region path. The interrupted
    # run may have left rows on disk past the last chunk the index recorded
    # (the BURDEN rows sit in an unflushed C++ stream, the index hits the disk
    # at once, so the two could be out of step in EITHER direction). Cut each
    # model's results file, and its .singleAssoc.txt, back to the size the
    # index recorded for the last completed chunk -- indexChunk itself here:
    # the region loop below resumes AT indexChunk + 1 (LEDGER #25/#81), so the
    # target is the end of chunk indexChunk. truncateOutputToIndex is in
    # Util.R; an index without counts (an older version) warns and appends.
    for (k in seq_len(nModels)) {
      truncateOutputToIndex(OutputFiles[k], outChecks[[k]]$chunkBytes[as.character(scan$indexChunk)])
      singleOut <- paste0(OutputFiles[k], ".singleAssoc.txt")
      if (file.exists(singleOut)) {
        truncateOutputToIndex(singleOut,
                              outChecks[[k]]$chunkSingleBytes[as.character(scan$indexChunk)])
      }
    }
  }

  n <- length(mu[[1]]) # sample size


  # Phase G1 slice 2: r.corr -> test plan (R/RegionState.R). The plan's
  # r.corr is what every later consumer reads: SKAT_Check_Method promotes 0
  # to the 7-element rho grid.
  plan <- regionTestPlan(r.corr, is_single_in_groupTest)
  r.corr <- plan$r.corr
  regionTestType <- plan$regionTestType
  is_single_in_groupTest <- plan$is_single_in_groupTest


  # Start the SKAT-O pool BEFORE the first C++ region call: with --nThreads > 1
  # that call owns a TBB worker pool, and forking a multithreaded process is the
  # hazard (PgenReader is not fork-safe either), so the workers are forked while
  # the parent is still single-threaded. It used to be created lazily by the
  # first region that took the parallel branch; the announcement line is the
  # same and still printed once.
  if (regionTestType != "BURDEN" && canParallelSkato(nThreadsSKATO)) {
    startSkatoCluster(nThreadsSKATO)
  }

  if (is_single_in_groupTest) {
    cat(
      "is_single_in_groupTest = TRUE. Single-variant assoc tests results will be output\n"
    )
    for (k in seq_len(nModels)) {
      isOpenOutFile_singleinGroup <- openOutfile_singleinGroup(k - 1L, traitTypes[k], isImputation, isappend, is_output_moreDetails)
      if (!isOpenOutFile_singleinGroup) {
        stop(
          "Output file ",
          OutputFiles[k],
          ".singleAssoc.txt can't be opened\n"
        )
      }
    }
  } else {
    cat(
      "is_single_in_groupTest = FALSE. Single-variant assoc tests results will not be output\n"
    )
  }


  ## check group file
  region_list <- checkGroupFile(groupFile, detail = TRUE)
  nRegions <- region_list$nRegions
  is_weight_included <- region_list$is_weight_included
  nWeightSets <- region_list$nWeightSets
  weightSetNames <- region_list$weightSetNames
  # FlexRV (F1): one optional per-marker score line per region. Parsed and
  # carried into RegionData as SAIGESC_1; NOTHING consumes it yet -- the
  # transform sets that expand it are F3. character(0) for every group file
  # with no score line, which is every bundled one, so this is inert.
  scoreName <- region_list$scoreName

  # ---- FlexRV transform sets (F3) -------------------------------------------
  # Kt = nScore x nMAF sets expanded from the group file's score line inside
  # C++ (src/FlexRV.cpp), named <score>_<scoreT>_<mafT> in score-major order
  # (docs/design/flexrv-design.md). The LoF mask keys on annotation labels,
  # which are matched exactly (SAIGE has no annotation vocabulary); a label
  # that matches nothing anywhere in the file is an error here, because the
  # mask would otherwise silently zero a whole score transform's sets (12 of
  # 192, or 7 of 112 at N = 1,000 once LEDGER #142 has dropped the rungs).
  flexTable <- matrix(0, nrow = 0, ncol = 4)
  flexNames <- character(0)
  flexLofAnno <- character(0)
  flexMaxMAF <- 0.001
  nTransformSets <- 0L
  if (!is.null(flexRV)) {
    if (length(scoreName) == 0) {
      stop("--flexRV_score=", flexRV$score, " needs a 'score:", flexRV$score,
           "' line in every region of the group file, and ", groupFile, " has no score line\n")
    }
    if (!identical(scoreName, flexRV$score)) {
      stop("--flexRV_score=", flexRV$score, " but the group file's score line is named '",
           scoreName, "'\n")
    }
    flexMaxMAF <- flexRV$maxMAF
    # LEDGER #142: the in-sample N decides which MAC rungs can differ from
    # unw (every admitted MAC is at most 2 * N * maxMAF, so a rung at or above
    # that bound weights every marker 1); those rungs are dropped here, before
    # composition, and disclosed below. `n` is length(mu[[1]]) above, the same count mainRegionInCPP is
    # handed, so the R-side bound and the C++ ingest bound read one N; a
    # manifest's models share ONE sample set on the region path.
    flexTT <- flexrvTransformTable(flexRV$scoreTransforms, flexRV$MAFTransforms,
                                   scoreName, maxMAF = flexMaxMAF, N = n)
    if (length(flexTT$dropped) > 0) {
      warning("FlexRV: at N = ", n, " and --flexRV_maxMAF=", flexMaxMAF,
              " every admitted marker has MAC <= ", flexTT$bound, ", so the MAC ",
              if (length(flexTT$dropped) == 1) "rung " else "rungs ",
              paste(sub("^MAC", "", flexTT$dropped), collapse = ", "),
              " would weight every marker 1 and duplicate 'unw'; ",
              if (length(flexTT$dropped) == 1) "it is" else "they are",
              " dropped before the sets are composed (LEDGER #142): ",
              nrow(flexTT$table), " sets (", flexTT$nScore, " score x ", flexTT$nMAF,
              " MAF), not ", flexTT$nScore * (flexTT$nMAF + length(flexTT$dropped)), "\n",
              call. = FALSE, immediate. = TRUE)
    }
    flexTable <- flexTT$table
    flexNames <- flexTT$names
    nTransformSets <- nrow(flexTable)
    if (flexTT$usesLof) {
      flexLofAnno <- parseFlexRVLofAnno(flexRV$lofAnno)
      absent <- setdiff(flexLofAnno, region_list$annoLabels)
      if (length(absent) > 0) {
        # LEDGER #170 (2026-09-27): a warning, not an error. This is a
        # FILE-level check, and a group file subset to a few regions (the AoU
        # rerun wrapper writes one per chromosome) can lack the label although
        # the same regions ran inside the whole-chromosome file; the error
        # killed that chromosome's process while the job reported success and
        # its regions were silently absent. A region without the label is
        # handled already: its 'lof' sets keep no marker, have MAC 0 and are
        # refused by --minGroupMAC_in_BurdenTest, disclosed in .pooledTests.txt.
        present <- intersect(flexLofAnno, region_list$annoLabels)
        warning("--flexRV_lofAnno label(s) ", paste(shQuote(absent), collapse = ", "),
                " match no annotation in the group file, whose labels are ",
                paste(shQuote(region_list$annoLabels), collapse = ", "), ". ",
                if (length(present) > 0) {
                  paste0("The 'lof' score transform keys on the label(s) present: ",
                         paste(shQuote(present), collapse = ", "), ".")
                } else {
                  paste0("No named label is present, so in every region the 'lof' ",
                         "transform's sets keep no marker (MAC 0) and are refused by ",
                         "--minGroupMAC_in_BurdenTest, as they are in any region lacking ",
                         "the label. If the label is misspelled, fix --flexRV_lofAnno.")
                }, " (LEDGER #170)\n", call. = FALSE, immediate. = TRUE)
      }
    }
    if (any(flexNames %in% weightSetNames)) {
      stop("FlexRV set name '", flexNames[flexNames %in% weightSetNames][1],
           "' collides with a weight set in the group file\n")
    }
    if (flexMaxMAF != 0.001) {
      warning("--flexRV_maxMAF=", flexMaxMAF, " is not the paper's 0.001: the mafbeta ",
              "transforms rescale to (1 - MAF/", flexMaxMAF, ")^(b-1) and are named ",
              "'mafbeta<b>@", flexMaxMAF, "'; the MAC ladder is absolute and does not ",
              "rescale: at large N it may cover only the rare end of the tested range, ",
              "and any rung at or above 2 * N * maxMAF is dropped and named in the ",
              "FlexRV warning above (LEDGER #142)\n",
              call. = FALSE, immediate. = TRUE)
    }
    cat("FlexRV: ", nTransformSets, " transform sets (", flexTT$nScore, " score x ",
        flexTT$nMAF, " MAF) from score '", scoreName, "' at max MAF ", flexMaxMAF, "\n",
        sep = "")
  }

  # ---- weight-set composition (design D5/D10) ------------------------------
  # The composed per-run weight-set list, in the NORMATIVE order shared with
  # C++ (src/Main.cpp mainRegionInCPP):
  #   custom sets (group-file weight lines, file order)
  #   -> FlexRV transform sets (F3; score-major over the transform table)
  #   -> explicit beta sets (weights.beta pairs in given order, named
  #      "Beta_a_b"; composed ONLY when more than one pair is supplied -- a
  #      single pair, default or not, keeps the legacy implicit-beta path)
  #   -> "No_Weights" (all-ones; appended when is_no_weight_in_groupTest =
  #      TRUE and the run has any other weight set)
  # At exactly one composed set every branch below reduces to the legacy
  # single-weight flow (same objects, same bytes).
  # Phase G1 slice 4: pure composition (R/RegionState.R); the restart
  # schema guard, the openOutfile loop, nline_per_gene and
  # region_weight_names stay below -- they are I/O and C++ state, not
  # composition.
  ws <- composeWeightSets(weights_beta_mat, is_no_weight_in_groupTest, is_weight_included,
                          nTransformSets, nWeightSets, weightSetNames, flexNames,
                          is_output_weight_column)
  nBetaPairs <- ws$nBetaPairs
  useExplicitBeta <- ws$useExplicitBeta
  betaSetNames <- ws$betaSetNames
  includeNoWeightsSet <- ws$includeNoWeightsSet
  legacyNoWeight <- ws$legacyNoWeight
  composedNames <- ws$composedNames
  nWeightSetsUsed <- ws$nWeightSetsUsed
  isMultiWeightRun <- ws$isMultiWeightRun
  isMultiWeightSchema <- ws$isMultiWeightSchema

  # Restart guard: the multi-weight schema adds columns (Weight,
  # Pvalue_ACATV, Pvalue_ACATO), so appending to a file written with a
  # different weight configuration would interleave two schemas in one file.
  for (existingOut in OutputFiles[isappend & file.exists(OutputFiles)]) {
    if (file.info(existingOut)$size == 0) next
    existingHeader <- strsplit(readLines(existingOut, n = 1), "[\t]")[[1]]
    hasWeightCol <- "Weight" %in% existingHeader
    if (hasWeightCol != isMultiWeightSchema) {
      stop(
        "The existing 'OutputFile' of '",
        existingOut,
        "' was written ",
        if (hasWeightCol) "WITH" else "WITHOUT",
        " the multi-weight columns (Weight, Pvalue_ACATV, Pvalue_ACATO), but this run uses ",
        nWeightSetsUsed,
        " weight set(s). Please restart with the same group file and options, or use another 'OutputFile' (or --is_overwrite_output=TRUE).\n"
      )
    }
  }

  if (regionTestType == "BURDEN") {
    # output the result from Rcpp (moved after checkGroupFile: at K > 1 the
    # C++-written header gains the Weight column right after max_MAF; at
    # K = 1 it is byte-identical to the legacy header)
    cat("isappend ", isappend, "\n")
    for (k in seq_len(nModels)) {
      isOpenOutFile <- openOutfile(k - 1L, traitTypes[k], isappend, nWeightSetsUsed, is_output_weight_column)
      if (!isOpenOutFile) {
        stop("Output file ", OutputFiles[k], " can't be opened\n")
      }
    }
  }

  nline_per_gene <- 2 + length(scoreName) + nWeightSets

  # Per-set labels handed to mainRegionInCPP; consumed by C++ only for
  # per-set output labels when the Weight column is emitted, and by
  # region_postprocess as the Weight values. The implicit beta set of an
  # unweighted default run is labeled "Beta" (visible only with
  # is_output_weight_column = TRUE).
  region_weight_names <- if (length(composedNames) > 0) {
    composedNames
  } else {
    "Beta"
  }

  # Per-set descriptors for region_postprocess: slot r aligns with C++ set r
  # and with region_weight_names. `kind` drives set r's R-side AnnoWeights
  # (design D2):
  #   custom        -> c(group-file weights, 1s for the UR pseudo-markers)
  #                    (C++ baked w_r into the UR collapse max(w_r*G))
  #   beta          -> dbeta(MAFVec, a_r, b_r) over every output row,
  #                    UR pseudo-markers included (their collapse is
  #                    unweighted, the beta weight of the collapsed MAF
  #                    applies like any marker's)
  #                    The implicit single set of an unweighted default run
  #                    is a beta set holding the first weights.beta pair
  #                    (LEDGER #19 FIXED 2026-09-05; it was a kind of its own,
  #                    hardcoded to dbeta(MAFVec, 1, 25) whatever the flag said)
  #   noweight      -> all ones
  weight_sets_kind <- if (legacyNoWeight) {
    # legacy all-ones vector: enters C++ as a 1-column custom set
    "custom"
  } else {
    c(
      rep("custom", nWeightSets),
      rep("transform", nTransformSets),
      rep("beta", if (useExplicitBeta) nBetaPairs else 0L),
      if (includeNoWeightsSet) "noweight"
    )
  }
  if (length(weight_sets_kind) == 0) {
    # the implicit single beta set of an unweighted default run: the first
    # (only) weights.beta pair, labeled "Beta"
    weight_sets_kind <- "beta"
  }
  weight_sets <- weightSetDescriptor(weight_sets_kind, region_weight_names,
                                     betaPairs = weights_beta_mat,
                                     force = is_output_weight_column)
  stopifnot(length(weight_sets_kind) == nWeightSetsUsed)

  gf <- file(groupFile, "r")

  cat("indexChunk is ", scan$indexChunk, "\n")

  skipline <- scan$indexChunk * nline_per_gene
  if (scan$indexChunk > 0 & scan$indexChunk < nRegions) {
    for (k in 1:skipline) {
      marker_group_line_temp <- readLines(gf, n = 1)
      rm(marker_group_line_temp)
    }
  }

  # P1Mat/P2Mat used to be allocated here at markers_per_chunk x n and copied
  # into C++ on EVERY region call (an arma::mat& parameter is a copy). The
  # region loop now sizes its own scratch once per call and keeps the chunks in
  # memory (see RegionChunkStore in src/Main.cpp), so the arguments are 1 x 1
  # placeholders for both test types.
  P1Mat <- matrix(0, 1, 1)
  P2Mat <- matrix(0, 1, 1)

  chrom1 <- "FakeCHR"


  gc()
  num_region <- 0
  mth <- 0

  numberRegionsInChunk <- 0
  cat("indexChunk ", scan$indexChunk, "\n")
  cat("nRegions ", nRegions, "\n")
  scan$pval.Region.all <- vector("list", nModels)   # per model
  scan$Output_MarkerList.all <- NULL                  # model-independent content, written to every model's file
  scan$cth_chunk_to_output <- 1
  # LEDGER #144 FIXED: seeded FALSE regardless of a restart, so the first
  # resumed chunk with marker-list rows took the fresh-write branch and
  # truncated .markerList.txt to the resumed regions. Seeded from isappend
  # like the results-file flag beside it (LEDGER #59).
  scan$iswriteMarkerList <- rep(isappend, nModels)
  scan$iswriteRegion <- rep(isappend, nModels)
  # LEDGER #116 FIXED: the SKAT-O path writes its header when the run STARTS,
  # from the schema in R/RegionSchema.R, so a run in which no cell qualifies
  # leaves a parseable, visibly empty file rather than a 0-byte one at exit 0.
  # The BURDEN path's C++ writer has always done this (openOutfile above); on
  # a restart the header is already on disk. The first chunk that has rows
  # re-checks the schema against them (regionAssertSchema, in the output
  # block) and stops on any drift, so this cannot silently disagree with them.
  if (regionTestType != "BURDEN" && !isappend) {
    for (k in seq_len(nModels)) {
      cat(regionHeaderLine(traitTypes[k], isCondition, isMultiWeightSchema),
          file = OutputFiles[k])
      scan$iswriteRegion[k] <- TRUE
    }
  }
  i <- scan$indexChunk + 1
  nEachChunk <- 1 # LEDGER #25 FIXED: was numberRegionsInChunk, still 0 here

  # All three mainRegionInCPP call sites passed an identical argument list;
  # the only per-site difference (which sparse-GRM flag was set beforehand)
  # stays at each call site.
  #
  # C++ weight-set boundary (jmr extension; the normative flat index is
  # jmr = j*q_maf*q_weight + m*q_weight + r, weight fastest — see the
  # matching comment in src/Main.cpp mainRegionInCPP). The legacy WEIGHT
  # object is canonicalized into the K-set matrix form without touching any
  # numeric value:
  #   - c(0) sentinel (unweighted file) -> 0x0 matrix: C++ sees no custom
  #     sets and falls back to the implicit beta set parameterized by
  #     weights.beta via setAssocTest_GlobalVarsInCPP, exactly the legacy
  #     beta-weights path. It must carry no COLUMNS rather than one, because
  #     the 1x1 form it used to take collides with the q0 x 1 form below at
  #     q0 == 1 and C++ can only tell them apart by shape (LEDGER #160)
  #   - length-q0 numeric vector (weighted file, or is_no_weight's
  #     rep(1, q0)) -> q0 x 1 matrix: one custom set (all-ones still counts
  #     as "equal weights" for the per-set ER gate)
  #   - q0 x K numeric matrix (K > 1 weight lines in the group file) -> passed
  #     through unchanged: K custom sets, one per group-file weight line
  # t_beta_params carries the explicit beta pairs (only when more than one
  # weights.beta pair was supplied -- a single pair keeps the 0-matrix
  # sentinel and the implicit legacy beta path) and t_includeNoWeights
  # appends the all-ones "No_Weights" set; C++ composes the sets in the same
  # custom -> beta -> No_Weights order as weight_sets/region_weight_names.
  # hasCustom comes from the region state record, NOT from WEIGHT's shape:
  # a plain numeric vector is a GENUINE one-column set in two cases
  # (legacyNoWeight's all-ones, a K = 1 matrix collapsed by as.numeric) and
  # the K = 0 `c(0)` sentinel in a third, and at q0 == 1 all three reach C++
  # as the same 1 x 1 matrix. C++ decides `hasCustomWeights` on
  # `t_weight.n_rows == q0` ALONE since #146(b) dropped the `!is_zero()` half
  # that used to separate them, so it read the sentinel as one custom set and
  # stopped on the count it could not reconcile with R's (LEDGER #160).
  # Required argument, so a caller that forgets it is an error rather than a
  # silently weightless run.
  weight_matrix_for_CPP <- function(WEIGHT, hasCustom) {
    if (!isTRUE(hasCustom)) {
      # no custom sets: no COLUMNS, so q_weight_custom is 0 at every q0
      return(matrix(numeric(0), nrow = 0L, ncol = 0L))
    }
    if (is.matrix(WEIGHT)) {
      matrix(as.numeric(WEIGHT), nrow = nrow(WEIGHT), ncol = ncol(WEIGHT))
    } else {
      matrix(as.numeric(WEIGHT), ncol = 1)
    }
  }
  beta_params_for_CPP <- if (useExplicitBeta) {
    matrix(as.numeric(weights_beta_mat), nrow = nBetaPairs, ncol = 2)
  } else {
    matrix(0, nrow = 1, ncol = 2)
  }
  # Step 4 (2026-09-06): mainRegionInCPP tests a vector of null models per
  # call and returns one output list per model; a single-model run passes
  # index 0 and takes the first (only) element, byte-identical to before.
  # S is the region's state record (built in the loop below).
  run_mainRegionInCPP <- function(idx, S) {
    mainRegionInCPP(
      genoType,
      S$region$genoIndex,
      S$annoIndicatorMat,
      maxMAFlist,
      OutputFile,
      # the FIRST model of this call's subset: C++ checks it against that
      # model (an escalation call may start at a quantitative model while the
      # run's first is binary -- sim50k, five traits)
      traitTypes[idx[1] + 1L],
      n,
      P1Mat,
      P2Mat,
      regionTestType,
      isImputation,
      weight_matrix_for_CPP(S$WEIGHT, S$hasCustomWeights),
      beta_params_for_CPP,
      includeNoWeightsSet,
      region_weight_names,
      is_output_weight_column,
      S$SCORE,
      S$ISLOF,
      flexTable,
      flexMaxMAF,
      weight_cond,
      is_single_in_groupTest,
      is_output_markerList_in_groupTest,
      S$annolistsub,
      S$regionName,
      is_fastTest,
      is_output_moreDetails,
      nThreads,
      as.integer(idx)
    )
  }

  # ---- Step 5 (2026-09-06): the SKAT-O cells batched across regions -------
  # With the pool (--nThreads > 1 on a SKAT-O run) a region's post-processing
  # is split into stages so that the cells of a BATCH of regions -- a
  # group-file chunk, or fewer when the held C++ output lists pass
  # SAIGE_SKATO_BATCH_MB -- go to the pool in ONE call (R/SkatoParallel.R):
  #   A  per region, in order: C++ fast pass, stash the single-variant temp
  #      rows, COLLECT the first-pass cells                      (loop body)
  #   B  one pool call for the batch                          (flushSkatoBatch)
  #   C  per region, in order: SERVE, decide escalation, C++ accurate pass
  #      for the escalated models, append the stashed rows for the others,
  #      COLLECT the escalated models' cells                  (regionEscalate)
  #   D  one pool call for the escalated cells, then per region in order:
  #      SERVE them, marker list, result rows                  (regionFinish)
  # Serial (--nThreads=1) and BURDEN runs take A, C, D region by region with
  # one region_postprocess drive each: the same stages, no batch, no pool.
  useSkatoBatch <- regionTestType != "BURDEN" && canParallelSkato(nThreadsSKATO)
  scan$skatoBatch <- list()
  scan$skatoBatchBytes <- 0
  skatoBatchBudget <- skatoBatchBudgetBytes()

  # region_postprocess's argument list for model k of region S.
  postArgs <- function(S, k, is_rerun) {
    list(
      S$outLists[[k]],
      mu[[k]],
      S$WEIGHT,
      traitTypes[k],
      r.corr,
      regionTestType,
      isCondition,
      S$annolistsub,
      maxMAFlist,
      S$regionName,
      is_single_in_groupTest,
      is_rerun = is_rerun,
      # the per-set carry-forward cache of the first pass; at K = 1 this is
      # list(resultDF), exactly the legacy carry-over (NULL on the first pass)
      resultDF_carryover = if (is_rerun) S$resultDFListK[[k]] else NULL,
      weight_sets = weight_sets,
      minGroupMAC_in_BurdenTest = minGroupMAC_in_BurdenTest,
      model_index = k
    )
  }
  # Bind one model's post-processing result into the region's state.
  takePost <- function(S, k, res) {
    # `[k] <- list(...)` rather than `[[k]] <-`: assigning NULL into a list
    # DELETES the element and shifts the models after it, which bound one
    # model's rows to another's on sim50k
    S$pvalRegionK[k] <- list(res$pval.Region)
    if ("cctpval" %in% names(res)) {
      S$cctpvalK[k] <- res$cctpval
    }
    if ("resultDFList" %in% names(res)) {
      S$resultDFListK[[k]] <- res$resultDFList
    }
    S
  }
  # Stage C for one region: escalate the models whose Cauchy p crosses the
  # cutoff (their accurate pass writes their single-variant rows itself) and
  # append the stashed fast-pass rows for the others. Batched, the escalated
  # models' cells are collected here and served in regionFinish.
  regionEscalate <- function(S) {
    # the counts are genotype-only, identical for every model
    outList <- S$outLists[[1]]
    S$runCore <- (sum(outList$NumUltraRare_GroupVec) + sum(outList$NumRare_GroupVec)) > 0
    S$esc <- integer(0)
    S$jobs2 <- list()
    if (S$runCore && regionTestType != "BURDEN" && is_fastTest) {
      # The models whose Cauchy p crosses the cutoff are re-run on the
      # accurate pass in ONE call; the others keep the fast pass and
      # their single-in-group rows are copied over. One model: the
      # legacy escalate-or-copy.
      esc <- which(S$cctpvalK < pval_cutoff_for_fastTest)
      if (length(esc) > 0) {
        # if(cctpval < 0.1 & is_fastTest){
        cat("Non-fast test is performed\n")
        set_flagSparseGRM_cur_SAIGE(TRUE)

        rerun <- run_mainRegionInCPP(modelIdx[esc], S)
        for (j in seq_along(esc)) {
          S$outLists[[esc[j]]] <- rerun[[j]]
        }
        S$esc <- esc
        if (useSkatoBatch) {
          S$jobs2 <- lapply(esc, function(k) skatoCollect(postArgs(S, k, TRUE)))
        } else {
          for (k in esc) {
            S <- takePost(S, k, do.call(
              region_postprocess,
              c(postArgs(S, k, TRUE), list(skat_fun = skatoMemoKernel()))
            ))
          }
        }
      }
      for (k in setdiff(seq_len(nModels), esc)) {
        # cctpval < 0.05)
        if (!is.null(S$tempSingle[[k]])) {
          append_singleInGroup(k - 1L, S$tempSingle[[k]])
        }
      }
    } # if(is_fastTest){
    S
  }
  # Stage D for one region: the marker list and the result rows, into the
  # chunk accumulators of the enclosing frame.
  regionFinish <- function(S) {
    outList <- S$outLists[[1]]
    SNP <- S$SNP
    annolistsub <- S$annolistsub
    regionName <- S$regionName

    if (!S$runCore) {
      # LEDGER #116 FIXED. The region was analysed and NO marker survived the
      # ingest QC gate (src/Main.cpp; e.g. monomorphic in the analysed sample
      # at --minMAC=0.5, the row's bundled reproduction), so region_postprocess
      # never ran. This used to RETURN here, before the index bookkeeping
      # below: the region contributed no row AND left indexChunk/Start/End
      # stale, so a run whose last region was such recorded "chunk 0" and no
      # completion line. Now each model gets the one row a region with nothing
      # tested contributes -- the "Region Cauchy NA" row the BURDEN writer has
      # always emitted, disclosed in the #127 sidecar with Rows = 0 -- and the
      # bookkeeping runs for both paths.
      if (regionTestType != "BURDEN") {
        for (k in seq_len(nModels)) {
          S$pvalRegionK[k] <- list(regionEmptyPooledRow(regionName, traitTypes[k],
                                                        isCondition, isMultiWeightSchema))
          pooledRecord(regionName, k, rows = 0L, refused = 0L, combined = 0L,
                       rule = pooledRule(numeric(0)), pvalue = NA_real_)
        }
      }
    }

    Output_MarkerList <- NULL
    if (S$runCore && is_output_markerList_in_groupTest) {
      for (j in 1:length(annolistsub)) {
        AnnoName <- annolistsub[j]
        for (m in 1:length(maxMAFlist)) {
          jm <- (j - 1) * (length(maxMAFlist)) + m
          maxMAFName <- maxMAFlist[m]
          # the marker list is weight-independent (real markers belong to
          # the same (anno, max_MAF) cell in every weight set), so read
          # the first weight set's column: jmr with r = 1, which is jm
          # when there is a single set
          jmr1 <- (jm - 1) * nWeightSetsUsed + 1
          tempPos <- which(outList$annoMAFIndicatorMat[, jmr1] > 0)
          marker_rare_pos <- which(outList$markerIndcatorVec == 1)
          marker_ultrarare_pos <- which(outList$markerIndcatorVec == 2)
          if (length(tempPos) > 0) {
            if (length(marker_rare_pos) > 0) {
              markerind_b <- which(marker_rare_pos %in% tempPos)
              if (length(markerind_b) > 0) {
                markerind <- marker_rare_pos[markerind_b]
                SNPlist_rare <- paste(SNP[markerind], collapse = ",")
              } else {
                SNPlist_rare <- ""
              }
            } else {
              SNPlist_rare <- ""
            }
            if (length(marker_ultrarare_pos) > 0) {
              markerind_UR_b <- which(marker_ultrarare_pos %in% tempPos)
              if (length(markerind_UR_b) > 0) {
                markerindUR <- marker_ultrarare_pos[markerind_UR_b]
                SNPlist_Ultra_rare <- paste(SNP[markerindUR], collapse = ",")
              } else {
                SNPlist_Ultra_rare <- ""
              }
            } else {
              SNPlist_Ultra_rare <- ""
            }
            Output_MarkerList <- rbind(
              Output_MarkerList,
              c(
                regionName,
                AnnoName,
                maxMAFName,
                SNPlist_rare,
                SNPlist_Ultra_rare
              )
            )
          }
        }
      }
    }

    if (is_output_markerList_in_groupTest) {
      if (!is.null(Output_MarkerList)) {
        colnames(Output_MarkerList) <- c(
          "Region",
          "Group",
          "max_MAF",
          "Rare_Variants",
          "Ultra_Rare_Variants"
        )
        scan$Output_MarkerList.all <- rbind(scan$Output_MarkerList.all, Output_MarkerList)
      }
    } else {
      scan$Output_MarkerList.all <- NULL
    }

    # LEDGER #81 FIXED: i was advanced past this region when it was read,
    # so the index recorded one past the last region and the completion
    # line was never written (End compared i with nRegions).
    scan$indexChunk <- S$i - 1
    # ... and on a RESTART the first chunk written must APPEND, not truncate:
    # Start drove append = FALSE from cth_chunk_to_output alone, which is 1
    # on a restart too (unreachable before #25, so never seen).
    scan$Start <- (scan$cth_chunk_to_output == 1) && !isappend
    scan$End <- (scan$indexChunk == nRegions)
    if (regionTestType != "BURDEN") {
      for (k in seq_len(nModels)) {
        # `[k] <<- list(...)`, not `[[k]] <<-`: rbind(NULL, NULL) is NULL, and
        # assigning NULL into a list DELETES the element and shifts every model
        # after it, so the NEXT region's `[[k]]` reads out of bounds and the run
        # dies with "subscript out of bounds". Latent until LEDGER #126 gave the
        # SKAT-O path a way to produce a region with no rows at all -- the same
        # trap `takePost` above is commented for, met a second time.
        scan$pval.Region.all[k] <- list(rbind(scan$pval.Region.all[[k]], S$pvalRegionK[[k]]))
      }
    }
    invisible(NULL)
  }
  # Stages B, C and D for the held batch: one pool call for the first-pass
  # cells, then per region in order serve / escalate / collect, one pool
  # call for the escalated cells, then per region in order serve / finish.
  # The answers are consumed sequentially in the order the jobs were
  # collected (region, then model), which skatoServe's counters check.
  flushSkatoBatch <- function() {
    cl <- startSkatoCluster(nThreadsSKATO)
    flatten <- function(x) {
      out <- do.call(c, do.call(c, x))
      if (is.null(out)) list() else out
    }
    jobs <- flatten(lapply(scan$skatoBatch, `[[`, "jobs"))
    # Printed per batch so a test can assert the batching actually happened
    # (a knob that never arrived would otherwise pass as green, LEDGER #30),
    # with the distinct count so the same holds for the identical-cell memo.
    dd <- skatoDedupe(jobs)
    cat("SKAT-O batch: ", length(jobs), " cell(s) from ", length(scan$skatoBatch),
        " region(s) to the pool (", length(dd$uniq), " distinct)\n", sep = "")
    answers <- skatoEvaluate(cl, jobs, dd)
    pos <- 0L
    for (b in seq_along(scan$skatoBatch)) {
      S <- scan$skatoBatch[[b]]
      for (k in seq_len(nModels)) {
        nk <- length(S$jobs[[k]])
        S <- takePost(S, k, skatoServe(postArgs(S, k, FALSE), answers[pos + seq_len(nk)]))
        pos <- pos + nk
      }
      scan$skatoBatch[[b]] <- regionEscalate(S)
    }
    jobs <- flatten(lapply(scan$skatoBatch, `[[`, "jobs2"))
    dd <- skatoDedupe(jobs)
    if (length(jobs) > 0) {
      cat("SKAT-O batch: ", length(jobs), " escalated cell(s) to the pool (",
          length(dd$uniq), " distinct)\n", sep = "")
    }
    answers <- skatoEvaluate(cl, jobs, dd)
    pos <- 0L
    for (b in seq_along(scan$skatoBatch)) {
      S <- scan$skatoBatch[[b]]
      for (j in seq_along(S$esc)) {
        nk <- length(S$jobs2[[j]])
        S <- takePost(S, S$esc[j], skatoServe(postArgs(S, S$esc[j], TRUE), answers[pos + seq_len(nk)]))
        pos <- pos + nk
      }
      regionFinish(S)
    }
    scan$skatoBatch <- list()
    scan$skatoBatchBytes <- 0
    invisible(NULL)
  }
  while (i <= nRegions) {
    # for(i in (indexChunk+1):nRegions){
    if (mth == numberRegionsInChunk) {
      if (i + groups_per_chunk > nRegions) {
        nregions_ro_read <- nRegions - i + 1
      } else {
        nregions_ro_read <- groups_per_chunk
      }
      nlinetoread <- nregions_ro_read * nline_per_gene
      marker_group_line <- readLines(gf, n = nlinetoread)
      RegionList <- SAIGE.getRegionList_new(
        marker_group_line,
        nline_per_gene,
        annolist,
        markerInfo,
        chrom,
        nWeightSets,
        weightSetNames,
        scoreName,
        lofAnno = flexLofAnno
      )
      # print(RegionList)
      cat(
        "Read in ",
        nregions_ro_read,
        " region(s) from the group file.\n"
      )
      mth <- 0
      # numberRegionsInChunk = length(RegionList)
      numberRegionsInChunk <- nregions_ro_read
    }

    mth <- mth + 1
    if (!is.null(RegionList)) {
      pval.Region <- NULL
      region <- RegionList[[mth]]

      annolistsub <- region$annoVec
      regionName <- names(RegionList)[mth]
      i <- i + 1
      if (!is.null(region$SNP) & length(annolistsub) > 0) {
        # Phase G1 slice 1: the per-region state record (R/RegionState.R)
        S <- regionStateFor(region, regionName, i, nRegions, nModels,
                            legacyNoWeight, nWeightSets, isMultiWeightSchema,
                            nTransformSets, flexLofAnno, scoreName)
        rm(region)

        if (!is_fastTest) {
          set_flagSparseGRM_cur_SAIGE_org()
        } else {
          set_flagSparseGRM_cur_SAIGE(FALSE)
        }
        S$outLists <- run_mainRegionInCPP(modelIdx, S)   # one list per model

        if (regionTestType == "BURDEN" & is_fastTest) {
          # BURDEN decides in C++ whether the fast pass's answer is kept
          # (iswriteOutput); the models it did not write are re-run on the
          # accurate pass together. One model: exactly the legacy re-run.
          need <- which(vapply(S$outLists, function(o) !is.null(o$iswriteOutput) && !isTRUE(o$iswriteOutput), logical(1)))
          if (length(need) > 0) {
            set_flagSparseGRM_cur_SAIGE(TRUE)


            S$outLists[need] <- run_mainRegionInCPP(modelIdx[need], S)
          }
        }

        # The fast pass wrote each model's single-variant rows to that
        # model's temp file, which the NEXT fast pass truncates. Read them
        # now; regionEscalate appends them for the models the Cauchy screen
        # clears, whenever that decision is made (step 5 defers it past the
        # batch's pool call).
        if (regionTestType != "BURDEN" && is_fastTest && is_single_in_groupTest) {
          S$tempSingle <- lapply(seq_len(nModels), function(k) {
            if (isTRUE(S$outLists[[k]]$iswriteOutput)) NULL else read_singleInGroupTemp(k - 1L)
          })
        }

        if (useSkatoBatch) {
          S$jobs <- lapply(seq_len(nModels), function(k) skatoCollect(postArgs(S, k, FALSE)))
          scan$skatoBatch[[length(scan$skatoBatch) + 1L]] <- S
          scan$skatoBatchBytes <- scan$skatoBatchBytes + as.numeric(utils::object.size(S$outLists))
        } else {
          for (k in seq_len(nModels)) {
            S <- takePost(S, k, do.call(
              region_postprocess,
              c(postArgs(S, k, FALSE), list(skat_fun = skatoMemoKernel()))
            ))
          }
          S <- regionEscalate(S)
          regionFinish(S)
        }
        rm(S)
      } else {
        # if(!is.null(region)){
        cat(regionName, " is empty.\n")
      }
      if (useSkatoBatch && length(scan$skatoBatch) > 0 &&
          (mth == numberRegionsInChunk || scan$skatoBatchBytes > skatoBatchBudget)) {
        flushSkatoBatch()
      }


      # output
      if (mth == numberRegionsInChunk) {
        # Fix for LEDGER #140: flush the region output streams held in C++,
        # the BURDEN rows and .singleAssoc.txt, FIRST; then record sizes; then
        # write the index line carrying both: that order is what makes the
        # recorded counts mean "everything on disk belongs to chunks up to
        # indexChunk", so a restart can truncate to them; the marker path's
        # #47 fix, S9. The R-side rows, fwrite below, are on disk on return.
        flushOutfile_regionInCPP()
        # cat("n1 is ", n1, "\n")
        # cat("n2 is ", n2, "\n")
        for (k in seq_len(nModels)) {
        if (regionTestType != "BURDEN") {
          # LEDGER #116: the rows about to be written must carry exactly the
          # schema whose header the run wrote at start (R/RegionSchema.R)
          regionAssertSchema(scan$pval.Region.all[[k]], traitTypes[k], isCondition, isMultiWeightSchema)
          # The truncate-vs-append decision is untouched; only the HEADER moves,
          # to the first chunk that actually HAS rows (LEDGER #59).
          if (scan$Start) {
            if (!is.null(scan$pval.Region.all[[k]])) {
              fwrite(
                scan$pval.Region.all[[k]],
                OutputFiles[k],
                quote = FALSE,
                sep = "\t",
                append = FALSE,
                col.names = TRUE,
                row.names = FALSE,
                na = "NA"
              )
              scan$iswriteRegion[k] <- TRUE
            }
          } else {
            if (!is.null(scan$pval.Region.all[[k]])) {
              fwrite(
                scan$pval.Region.all[[k]],
                OutputFiles[k],
                quote = FALSE,
                sep = "\t",
                append = TRUE,
                col.names = !scan$iswriteRegion[k],
                row.names = FALSE,
                na = "NA"
              )
              scan$iswriteRegion[k] <- TRUE
            }
          }
        }
        if (is_output_markerList_in_groupTest) {
          markerListFile <- paste0(OutputFiles[k], ".markerList.txt")
          if (scan$Start) {
            if (!is.null(scan$Output_MarkerList.all)) {
              fwrite(
                scan$Output_MarkerList.all,
                markerListFile,
                quote = FALSE,
                sep = "\t",
                append = FALSE,
                col.names = TRUE,
                row.names = FALSE,
                na = "NA"
              )
              scan$iswriteMarkerList[k] <- TRUE
            }
          } else {
            if (!is.null(scan$Output_MarkerList.all)) {
              if (scan$iswriteMarkerList[k]) {
                fwrite(
                  scan$Output_MarkerList.all,
                  markerListFile,
                  quote = FALSE,
                  sep = "\t",
                  append = TRUE,
                  col.names = FALSE,
                  row.names = FALSE,
                  na = "NA"
                )
              } else {
                fwrite(
                  scan$Output_MarkerList.all,
                  markerListFile,
                  quote = FALSE,
                  sep = "\t",
                  append = FALSE,
                  col.names = TRUE,
                  row.names = FALSE,
                  na = "NA"
                )
              }
              scan$iswriteMarkerList[k] <- TRUE
            }
          }
        }
        # LEDGER #140: this model's line carries the sizes measured AFTER its
        # rows above were written. The index file itself is written by the
        # shared writeOutputFileIndex (E6, 2026-09-14) -- the marker path's
        # writer, whose singleBytes parameter had no caller until now -- so the
        # analysis-type and nEachChunk lines are built by the same code
        # checkOutputFile compares them against (the LEDGER #25 class). Its
        # "write to output" line now prints once per model, after the flush.
        regionBytes <- if (file.exists(OutputFiles[k])) file.size(OutputFiles[k]) else NA_real_
        singleOut <- paste0(OutputFiles[k], ".singleAssoc.txt")
        singleBytes <- if (file.exists(singleOut)) file.size(singleOut) else NA_real_
        writeOutputFileIndex(OutputFileIndexes[k], "Region", nEachChunk, scan$indexChunk,
                             scan$Start, scan$End, outputBytes = regionBytes, singleBytes = singleBytes)
        } # for k (models)
        # }#if(FALSE)

        scan$pval.Region.all <- vector("list", nModels)
        scan$Output_MarkerList.all <- NULL
        scan$cth_chunk_to_output <- scan$cth_chunk_to_output + 1
        # (per-iteration gc() removed in Phase 4: it dominated the profile; R collects on demand)
      }
    } else {
      # if(!is.null(RegionList)){
      cat("The chunk is empty\n")
      # mth = 0
      mth <- numberRegionsInChunk
      i <- i + numberRegionsInChunk
      pval.Region <- NULL
    }
  }

  # (a first "Analysis done! ... .markerInfo" message used to be assembled
  # here and overwritten unconditionally below; deleted, Phase G1 slice 3)

  # LEDGER #57: how many cells the SKAT kernel failed on (their rows are NA),
  # and a sidecar naming each with its message when that count is not zero.
  # With several models the count is over all of them and the sidecar goes
  # beside the FIRST model's file (the per-cell record has no model column).
  skatFailureSummary(OutputFile)
  spaFallbackSummary(OutputFile) # LEDGER #164
  exactByWeightSummary(OutputFile) # LEDGER #172
  exactByWeightAboveCapSummary(OutputFile) # LEDGER #176: observation only
  stretchGateSummary(OutputFile) # LEDGER #176: observation only
  spaPinnedSummary(OutputFile) # LEDGER #174
  skatoDedupeSummary()
  burdenMemoSummary()
  # LEDGER #127: how each region's pooled (Cauchy) row was formed -- cells
  # tested, cells the min-group-MAC gate refused, p-values that entered, and
  # which CCT rule produced the number -- with a sidecar beside the FIRST
  # model's file when any region has something to disclose.
  pooledSummary(OutputFile)
  # LEDGER #158: which path produced each cell's Pvalue, and the mode the run
  # took (always printed), with a sidecar when the screen is on.
  skatoMethodSummary(OutputFile)

  # Phase G1 slice 3 (R/RegionState.R)
  regionDoneMessage(OutputFiles, OutputFile, is_output_markerList_in_groupTest,
                    is_single_in_groupTest)
}


SAIGE.getRegionList_new <- function(marker_group_line,
                                    nline_per_gene,
                                    annoVec,
                                    # c("lof","lof;missense"
                                    markerInfo,
                                    chrom = "",
                                    nWeightSets = NULL,
                                    weightSetNames = NULL,
                                    scoreName = character(0),
                                    # FlexRV: annotation labels that mark a
                                    # marker as LoF; ISLOF is carried per region
                                    # only when this is non-empty
                                    lofAnno = character(0)) {
  if (is.null(nWeightSets)) {
    # legacy callers: infer the weight layout from nline_per_gene
    # (2 -> no weights; 3 -> one bare weight line, set name "WEIGHT")
    if (nline_per_gene == 2) {
      nWeightSets <- 0L
      weightSetNames <- character(0)
    } else if (nline_per_gene == 3) {
      nWeightSets <- 1L
      weightSetNames <- "WEIGHT"
    } else {
      stop("nWeightSets must be provided when nline_per_gene is not 2 or 3\n")
    }
  }
  # one optional FlexRV score line sits between anno and the weight lines
  nScoreSets <- length(scoreName)
  if (nline_per_gene != 2 + nScoreSets + nWeightSets) {
    # the historical message is kept verbatim when no score line is in play
    # (test-unit-regionlist.R pins it), and only mentions scores when they
    # are actually part of the layout
    if (nScoreSets == 0) {
      stop(
        "nline_per_gene (", nline_per_gene,
        ") must equal 2 + nWeightSets (nWeightSets = ", nWeightSets, ")\n"
      )
    }
    stop(
      "nline_per_gene (", nline_per_gene,
      ") must equal 2 + nScoreSets + nWeightSets (nScoreSets = ", nScoreSets,
      ", nWeightSets = ", nWeightSets, ")\n"
    )
  }
  if (nWeightSets > 0 & length(weightSetNames) != nWeightSets) {
    stop("weightSetNames must contain one name per weight set\n")
  }
  # weight columns live inside RegionData under sanitized internal names so
  # they follow every merge/filter/sort below; a user-facing set name such as
  # "POS", "CHROM", "ID", "SNP" or "genoIndex" must never collide with
  # markerInfo columns in the merges
  # paste0 recycles zero-length args to "" (recycle0=FALSE default), so the
  # K=0 case must be explicit or this becomes a length-1 "SAIGEWT_"
  weightColNames <- if (nWeightSets > 0) {
    paste0("SAIGEWT_", seq_len(nWeightSets))
  } else {
    character(0)
  }
  # the score rides the same sanitized-internal-name path as the weights, for
  # the same reason: it must survive every merge/filter/sort below without
  # colliding with a markerInfo column. Nothing reads it yet (FlexRV F1).
  scoreColNames <- if (nScoreSets > 0) "SAIGESC_1" else character(0)
  chrom_nochr <- gsub("CHR", "", chrom, ignore.case = TRUE)
  # read group file
  ngroup <- length(marker_group_line) / nline_per_gene
  # cat("ngroup is ", ngroup, "\n")
  RegionData <- NULL
  geneList <- c()
  for (i in 1:ngroup) {
    marker_group_line_list <- strsplit(marker_group_line[1 + (i - 1) * nline_per_gene],
      split =
        "[\ \t]+"
    )[[1]]
    gene <- marker_group_line_list[1]
    var <- marker_group_line_list[3:length(marker_group_line_list)]
    marker_group_line_list_anno <- strsplit(marker_group_line[2 + (i -
      1) * nline_per_gene], split = "[\ \t]+")[[1]]
    anno <- marker_group_line_list_anno[3:length(marker_group_line_list_anno)]
    score <- NULL
    if (nScoreSets > 0) {
      idxScore <- 3 + (i - 1) * nline_per_gene
      marker_group_line_list_score <- strsplit(marker_group_line[idxScore],
        split = "[\ \t]+"
      )[[1]]
      score <- cbind(
        marker_group_line_list_score[3:length(marker_group_line_list_score)]
      )
      colnames(score) <- scoreColNames
    }
    if (nWeightSets > 0) {
      weight <- NULL
      for (k in 1:nWeightSets) {
        idxWeight <- 2 + nScoreSets + k + (i - 1) * nline_per_gene
        marker_group_line_list_weight <- strsplit(marker_group_line[idxWeight],
          split = "[\ \t]+"
        )[[1]]
        weight <- cbind(
          weight,
          marker_group_line_list_weight[3:length(marker_group_line_list_weight)]
        )
      }
      colnames(weight) <- weightColNames
      RegionData <- rbind(RegionData, cbind(rep(gene, length(var)), var, anno, score, weight))
    } else {
      RegionData <- rbind(RegionData, cbind(rep(gene, length(var)), var, anno, score))
    }

    if (gene %in% geneList) {
      stop(gene, " is duplicated in the group File\n")
    } else {
      geneList <- c(geneList, gene)
    }
  }
  colnames(RegionData) <- c("REGION", "SNP", "ANNO", scoreColNames, weightColNames)
  RegionData <- as.data.frame(RegionData)
  setDT(RegionData)
  uRegion0 <- unique(RegionData$REGION)

  if (chrom != "" & is.null(markerInfo)) {
    RegionData[, c("chr") := tstrsplit(RegionData$SNP, ":")[[1]]]
    setkey(RegionData, "chr")
    RegionData <- RegionData[chr == chrom | chr == chrom_nochr]
    RegionData[, chr := NULL]
  }

  if (nrow(RegionData) != 0) {
    if (!is.null(markerInfo)) {
      setkey(RegionData, "SNP")
      RegionData <- merge(
        RegionData,
        markerInfo,
        by.x = "SNP",
        by.y = "ID",
        all.x = TRUE,
        sort = FALSE
      )

      if (!is.null(markerInfo$ID2)) {
        RegionData <- merge(
          RegionData,
          markerInfo,
          by.x = "SNP",
          by.y = "ID2",
          all.x = TRUE,
          sort = FALSE
        )
        # SNP REGION     ANNO CHROM.x POS.x genoIndex2.x  ID2 genoIndex_prev.x CHROM.y POS.y    ID genoIndex2.y genoIndex_prev.y
        setnames(RegionData, "genoIndex.x", "genoIndex")
        setnames(RegionData, "genoIndex.y", "genoIndex2")
        # setnames(RegionData, "genoIndex.y", "genoIndex")
        setnames(RegionData, "CHROM.x", "CHROM")
        setnames(RegionData, "CHROM.y", "CHROM2")
        setnames(RegionData, "POS.x", "POS")
        setnames(RegionData, "POS.y", "POS2")

        posNA <- which(is.na(RegionData$genoIndex) &
          !is.na(RegionData$genoIndex2))

        if (length(posNA) != 0) {
          RegionData$genoIndex[posNA] <- RegionData$genoIndex2[posNA]
          RegionData$CHROM[posNA] <- RegionData$CHROM2[posNA]
          RegionData$POS[posNA] <- RegionData$POS2[posNA]
        }
        # RegionData$genoIndex2 = NULL
      }
      posNA <- which(is.na(RegionData$genoIndex))

      nRegionMarkers <- nrow(RegionData)
      if (length(posNA) != 0) {
        RegionData <- RegionData[which(!is.na(RegionData$genoIndex))]
        cat(
          length(posNA),
          " markers in 'RegionFile' are not in 'GenoFile'.\n"
        )
      }
      # LEDGER #113: NOT ONE marker of the group file joined the genotype file.
      # Every region is then empty, the chunk loop has nothing to do and the
      # writer is reached with no header and no rows -- a zero-byte results file
      # at exit 0 with "Analysis done!", which nothing downstream can tell from
      # a run that legitimately found nothing. This is the emptiness guard #76
      # gave the --chrom filter, on the route #76 does not cover: there the
      # filter empties the marker set, here the ID JOIN does.
      #
      # The condition is deliberately "no marker matched AT ALL". Individual
      # empty regions are normal and are handled inside the loop below
      # (`if (length(posSNP) > 0)`), so they must not trip this.
      if (nrow(RegionData) == 0) {
        stop(
          "None of the ", nRegionMarkers, " markers in the group file were ",
          "found in the genotype file, so every region is empty. PLINK keys on ",
          "column 2 of the .bim and PGEN on the ID column of the .pvar, so a ",
          "group file written for another format -- chr:pos:ref:alt IDs, say -- ",
          "matches nothing here. Check that the group file's marker IDs are the ",
          "same identifiers the variant file uses.\n"
        )
      }
      setorderv(RegionData, col = c("CHROM", "POS"))
    }
  }


  if (nrow(RegionData) != 0) {
    # HeaderInRegionData = colnames(RegionData)
    HeaderInRegionData <- unique(RegionData$ANNO)
    RegionAnnoHeaderList <- list()
    if (length(annoVec) == 0) {
      stop("At least one annotation is required\n")
    }
    for (q in 1:length(annoVec)) {
      RegionAnnoHeaderList[[q]] <- strsplit(annoVec[q], ";")[[1]]
    }

    RegionList <- list()
    uRegion <- unique(RegionData$REGION)
    # RegionData = as.data.frame(RegionData)


    for (r in uRegion0) {
      # print(paste0("Analyzing region ",r,"...."))
      # print(RegionData$REGION)
      # print(r)
      # which(as.vector(RegionData$REGION) == r)

      posSNP <- which(RegionData$REGION == r)
      if (length(posSNP) > 0) {
        SNP <- RegionData$SNP[posSNP]

        if (nWeightSets > 0) {
          # character matrix, one column per weight set, rows aligned to SNP
          WEIGHT <- as.matrix(RegionData[posSNP, weightColNames, with = FALSE])
        }
        if (nScoreSets > 0) {
          # one per-marker score, rows aligned to SNP, subset with them below
          SCORE <- RegionData[posSNP, scoreColNames, with = FALSE][[1]]
        }
        if (length(lofAnno) > 0) {
          ISLOF <- as.numeric(RegionData$ANNO[posSNP] %in% lofAnno)
        }

        if (any(duplicated(SNP))) {
          stop(
            "Please check RegionFile: in region ",
            r,
            ": duplicated SNPs exist."
          )
        }

        if (!is.null(markerInfo)) {
          genoIndex <- as.numeric(RegionData$genoIndex[posSNP])
          chrom <- RegionData$CHROM[posSNP]
          # uchrom = unique(chrom)
          # if(length(uchrom) != 1)
          #  stop("In region ",r,", markers are from multiple chromosomes.")
        }

        annoIndicatorMat <- matrix(0, nrow = length(posSNP), ncol = length(annoVec))
        annoVecNew <- c()
        for (q in 1:length(annoVec)) {
          indiceVec <- which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]])
          if (length(indiceVec) > 0) {
            annoVecNew <- c(annoVecNew, annoVec[q])
            annoIndicatorMat[indiceVec, q] <- 1
            # annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
          }
        }

        RegionAnnoHeaderListNew <- list()
        if (length(annoVecNew) == 0) {
          warning(
            "No markers are found for at least one annotation, so region ",
            r,
            " is skipped\n"
          )
          # stop("At least one annotation is required\n")
        } else {
          if (length(annoVecNew) < length(annoVec)) {
            annoIndicatorMat <- matrix(0,
              nrow = length(posSNP),
              ncol = length(annoVecNew)
            )
            for (q in 1:length(annoVecNew)) {
              RegionAnnoHeaderListNew[[q]] <- strsplit(annoVecNew[q], ";")[[1]]
              indiceVec <- which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderListNew[[q]])
              # if(length(indiceVec) > 0){
              annoIndicatorMat[indiceVec, q] <- 1
              # annoIndicatorMat[which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]]), q] = 1
            }
          } else {
            annoVecNew <- annoVec
          }
        }

        annoIndicatorMat_rmind <- which(rowSums(annoIndicatorMat) == 0)
        if (length(annoIndicatorMat_rmind) > 0) {
          SNP <- SNP[-annoIndicatorMat_rmind]

          if (nWeightSets > 0) {
            WEIGHT <- WEIGHT[-annoIndicatorMat_rmind, , drop = FALSE]
          }
          if (nScoreSets > 0) {
            SCORE <- SCORE[-annoIndicatorMat_rmind]
          }
          if (length(lofAnno) > 0) {
            ISLOF <- ISLOF[-annoIndicatorMat_rmind]
          }

          annoIndicatorMat <- annoIndicatorMat[-annoIndicatorMat_rmind, ,
            drop =
              FALSE
          ]
          if (!is.null(markerInfo)) {
            genoIndex <- genoIndex[-annoIndicatorMat_rmind]
          }
        }

        if (nWeightSets == 0) {
          WEIGHT <- c(0)
        } else {
          # numeric matrix with the user-facing set names, in group-file
          # order; for a single set the caller collapses this back to the
          # legacy unnamed numeric vector
          WEIGHT <- matrix(as.numeric(WEIGHT),
            nrow = nrow(WEIGHT),
            ncol = nWeightSets,
            dimnames = list(NULL, weightSetNames)
          )
        }

        if (!is.null(markerInfo)) {
          RegionList[[r]] <- list(
            SNP = SNP,
            WEIGHT = WEIGHT,
            annoIndicatorMat = annoIndicatorMat,
            genoIndex = as.character(format(genoIndex, scientific = FALSE)),
            #                           chrom = uchrom,
            annoVec = annoVecNew
          )
          # appended only when the group file carries a score line, so the
          # entry shape is byte-for-byte what it was for every existing file
          if (nScoreSets > 0) {
            RegionList[[r]]$SCORE <- as.numeric(SCORE)
          }
          if (length(lofAnno) > 0) {
            RegionList[[r]]$ISLOF <- ISLOF
          }
        }
      } else {
        # if(length(posSNP) > 0){
        RegionList[[r]] <- list(SNP = NULL)
      }
    }
  } else {
    # if(nrow(RegionData) == 0){
    RegionList <- NULL
  }


  RegionList
}


mainRegionURV <- function(NullModelClass = "SAIGE_NULL_Model",
                          genoType,
                          genoIndex,
                          n) {
  if (NullModelClass == "SAIGE_NULL_Model") {
    obj.mainRegionURV <- mainRegionURVInCPP("SAIGE", genoType, genoIndex, n)
  }

  obj.mainRegionURV
}
