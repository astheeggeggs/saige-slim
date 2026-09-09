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
# resultDF_carryover reproduces the re-run copy's visibility of the primary
# pass's resultDF (pass list(resultDF), or NULL when it does not exist).
# The returned list carries pval.Region plus cctpval/resultDF only when they
# were bound, so the caller can reproduce the exact frame state the inline
# blocks used to leave behind.
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
                               # The SKAT-O evaluator, injected so the cell loop
                               # can be driven twice by skato_parallel(): once to
                               # COLLECT every (Score, Phi) it would test, once to
                               # SERVE the answers computed in parallel in
                               # between. Defaults to the real thing, which is
                               # the entire single-threaded path -- at
                               # nThreadsSKATO = 1 nothing below changes at all.
                               skat_fun = get_SKAT_pvalue) {
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
  if (!is.null(resultDF_carryover)) {
    resultDF <- resultDF_carryover[[1]]
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
          transform = c(outList$transformWeightMat[, ws_transformCol[r]], rep(1, outList$numofUR)),
          stop("region_postprocess: unknown weight-set kind '", ws_kind[r], "'\n")
        )
      })
    } else if (sum(WEIGHT) > 0) {
      AnnoWeights <- c(WEIGHT, rep(1, outList$numofUR))
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
    } else if (sum(WEIGHT) > 0) {
      # AnnoWeights = AnnoWeights[varTestedIndices]
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
            wStatVec <- wStatVecList[[r]]
          }
          if (m <= maxMAF0) {
            tempPos <- which(annoMAFIndicatorMat[, jmr] == 1)
            if (length(tempPos) > 0) {
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

                # re_phi = get_newPhi_scaleFactor(q.sum,
                re_phi <- get_newPhi_scaleFactor_traitType(
                  q.sum,
                  mu.a,
                  g.sum,
                  p.new,
                  Score,
                  Phi,
                  regionTestType,
                  traitType
                )
                Phi <- re_phi$val
              }
              groupOutList <- skat_fun(Score, Phi, r.corr, regionTestType)
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
                Pvalue_ACATV <- get_ACATV_pvalue(exp(adjPVec[tempPos]), AnnoWeights[tempPos])
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

    if (length(annolistsub) > 1 | length(maxMAFlist) > 1 | isMultiWeight) {
      if (is.null(pval.Region)) {
        # Degenerate corner: run_core was TRUE (the region has rare/ultra-rare
        # carriers) but no (anno, max_MAF, weight set) cell emitted a row --
        # e.g. every flagged marker's p-value came back NA, or a weighted URV
        # collapse produced an all-zero pseudo-genotype. ncol(NULL) is NULL, so
        # the named construction below would error; reproduce the legacy
        # construction verbatim instead (every get_CCT_pvalue(NULL) is NA and
        # rbind(NULL, cctVec) is a 1-row character matrix), so the run continues
        # exactly as at HEAD. Multi-weight mode reuses the same degenerate row.
        cctpval_Burden <- get_CCT_pvalue(pval.Region$Pvalue_Burden)
        if (regionTestType != "BURDEN") {
          cctpval <- get_CCT_pvalue(pval.Region$Pvalue)
          cctpval_SKAT <- get_CCT_pvalue(pval.Region$Pvalue_SKAT)
          cctVec <- c(
            regionName,
            "Cauchy",
            NA,
            cctpval,
            cctpval_Burden,
            cctpval_SKAT,
            NA,
            NA
          )
        } else {
          cctVec <- c(
            regionName,
            "Cauchy",
            NA,
            cctpval_Burden,
            NA,
            NA
          )
        }
        if (isCondition) {
          cctpval_Burden <- get_CCT_pvalue(pval.Region$Pvalue_Burden_cond)
          if (regionTestType != "BURDEN") {
            cctpval <- get_CCT_pvalue(pval.Region$Pvalue_cond)
            cctpval_SKAT <- get_CCT_pvalue(pval.Region$Pvalue_SKAT_cond)
            cctVec <- c(
              cctVec,
              cctpval,
              cctpval_Burden,
              cctpval_SKAT,
              NA,
              NA
            )
          } else {
            cctVec <- c(cctVec, cctpval_Burden, NA, NA)
          }
        }
        cctVec <- c(cctVec, NA)
        if (traitType == "binary") {
          cctVec <- c(cctVec, NA, NA)
        }
        cctVec <- c(cctVec, NA, NA)
        pval.Region <- rbind(pval.Region, cctVec)
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

        pval.Region <- rbind(pval.Region, unlist(cctRow))
      } # if (is.null(pval.Region)) ... else
    } else {
      if (!is_rerun) {
        cctpval <- 1
      }
    }
  }

  ret <- list(pval.Region = pval.Region)
  if (exists("cctpval", inherits = FALSE)) {
    ret$cctpval <- cctpval
  }
  if (exists("resultDF", inherits = FALSE)) {
    ret$resultDF <- resultDF
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
  skatoDedupeReset() # identical cells evaluated once: per-run counts, R/SkatoParallel.R

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

  indexChunk <- outList$indexChunk
  Start <- outList$Start
  End <- outList$End

  cat("Start ", Start, "\n")
  cat("End ", End, "\n")


  if (End) {
    message <- paste0(
      "The analysis has been completed in earlier analysis. Results are saved in '",
      paste(OutputFiles, collapse = "', '"),
      "'. ",
      "If you want to change parameters and restart the analysis, please use another 'OutputFile'."
    )
    return(message)
  }

  isappend <- FALSE
  if (!Start) {
    isappend <- TRUE
  }

  n <- length(mu[[1]]) # sample size


  if (r.corr == 0) {
    out.method <- SKAT:::SKAT_Check_Method(method = "optimal.adj", r.corr = 0)
    method <- out.method$method
    r.corr <- out.method$r.corr
    cat(
      "SKAT-O test will be performed. P-values for BURDEN and SKAT tests will also be output\n"
    )
    regionTestType <- "SKAT-O"
    is_single_in_groupTest <- TRUE
  } else if (r.corr == 1) {
    method <- NULL
    cat("BURDEN test will be performed\n")
    regionTestType <- "BURDEN"
    # the BURDEN output file (written from Rcpp) is opened after the group
    # file check below: its header needs the number of weight sets K
  } else {
    stop("r.corr needs to be either 1 (BURDEN test) or 0 (SKAT-O test)\n")
  }


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
  # mask would otherwise silently zero 12 of the 192 sets.
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
    flexTT <- flexrvTransformTable(flexRV$scoreTransforms, flexRV$MAFTransforms,
                                   scoreName, maxMAF = flexMaxMAF)
    flexTable <- flexTT$table
    flexNames <- flexTT$names
    nTransformSets <- nrow(flexTable)
    if (flexTT$usesLof) {
      flexLofAnno <- parseFlexRVLofAnno(flexRV$lofAnno)
      absent <- setdiff(flexLofAnno, region_list$annoLabels)
      if (length(absent) > 0) {
        stop("--flexRV_lofAnno label(s) ", paste(shQuote(absent), collapse = ", "),
             " match no annotation in the group file, whose labels are ",
             paste(shQuote(region_list$annoLabels), collapse = ", "),
             ". The LoF mask (the 'lof' score transform) would be all zeros.\n")
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
              "rescale, so it may cover only the rare end of the tested range\n",
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
    cat("is_no_weight_in_groupTest = TRUE: an all-ones 'No_Weights' weight set is appended to the group test\n")
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
        "' in the composed weight-set list (group-file sets, then weights.beta pairs, then No_Weights). Please rename the group-file weight set or drop the duplicate pair.\n"
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

  cat("indexChunk is ", indexChunk, "\n")

  skipline <- indexChunk * nline_per_gene
  if (indexChunk > 0 & indexChunk < nRegions) {
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
  cat("indexChunk ", indexChunk, "\n")
  cat("nRegions ", nRegions, "\n")
  pval.Region.all <- vector("list", nModels)   # per model
  OutList.all <- NULL
  Output_MarkerList.all <- NULL                  # model-independent content, written to every model's file
  cth_chunk_to_output <- 1
  iswriteMarkerList <- rep(FALSE, nModels)
  iswriteRegion <- rep(isappend, nModels)
  i <- indexChunk + 1
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
  #   - c(0) sentinel (unweighted file) -> 1x1 zero matrix: C++ sees no
  #     custom sets and falls back to the implicit beta set parameterized by
  #     weights.beta via setAssocTest_GlobalVarsInCPP, exactly the legacy
  #     beta-weights path
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
  weight_matrix_for_CPP <- function(WEIGHT) {
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
      weight_matrix_for_CPP(S$WEIGHT),
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
  skatoBatch <- list()
  skatoBatchBytes <- 0
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
      weight_sets = weight_sets
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
    if (!S$runCore) {
      return(invisible(NULL))
    }
    outList <- S$outLists[[1]]
    SNP <- S$SNP
    annolistsub <- S$annolistsub
    regionName <- S$regionName

    Output_MarkerList <- NULL
    if (is_output_markerList_in_groupTest) {
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
      colnames(Output_MarkerList) <- c(
        "Region",
        "Group",
        "max_MAF",
        "Rare_Variants",
        "Ultra_Rare_Variants"
      )
      Output_MarkerList.all <<- rbind(Output_MarkerList.all, Output_MarkerList)
    } else {
      Output_MarkerList.all <<- NULL
    }

    # LEDGER #81 FIXED: i was advanced past this region when it was read,
    # so the index recorded one past the last region and the completion
    # line was never written (End compared i with nRegions).
    indexChunk <<- S$i - 1
    # ... and on a RESTART the first chunk written must APPEND, not truncate:
    # Start drove append = FALSE from cth_chunk_to_output alone, which is 1
    # on a restart too (unreachable before #25, so never seen).
    Start <<- (cth_chunk_to_output == 1) && !isappend
    End <<- (indexChunk == nRegions)
    if (regionTestType != "BURDEN") {
      for (k in seq_len(nModels)) {
        pval.Region.all[[k]] <<- rbind(pval.Region.all[[k]], S$pvalRegionK[[k]])
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
    jobs <- flatten(lapply(skatoBatch, `[[`, "jobs"))
    # Printed per batch so a test can assert the batching actually happened
    # (a knob that never arrived would otherwise pass as green, LEDGER #30),
    # with the distinct count so the same holds for the identical-cell memo.
    dd <- skatoDedupe(jobs)
    cat("SKAT-O batch: ", length(jobs), " cell(s) from ", length(skatoBatch),
        " region(s) to the pool (", length(dd$uniq), " distinct)\n", sep = "")
    answers <- skatoEvaluate(cl, jobs, dd)
    pos <- 0L
    for (b in seq_along(skatoBatch)) {
      S <- skatoBatch[[b]]
      for (k in seq_len(nModels)) {
        nk <- length(S$jobs[[k]])
        S <- takePost(S, k, skatoServe(postArgs(S, k, FALSE), answers[pos + seq_len(nk)]))
        pos <- pos + nk
      }
      skatoBatch[[b]] <<- regionEscalate(S)
    }
    jobs <- flatten(lapply(skatoBatch, `[[`, "jobs2"))
    dd <- skatoDedupe(jobs)
    if (length(jobs) > 0) {
      cat("SKAT-O batch: ", length(jobs), " escalated cell(s) to the pool (",
          length(dd$uniq), " distinct)\n", sep = "")
    }
    answers <- skatoEvaluate(cl, jobs, dd)
    pos <- 0L
    for (b in seq_along(skatoBatch)) {
      S <- skatoBatch[[b]]
      for (j in seq_along(S$esc)) {
        nk <- length(S$jobs2[[j]])
        S <- takePost(S, S$esc[j], skatoServe(postArgs(S, S$esc[j], TRUE), answers[pos + seq_len(nk)]))
        pos <- pos + nk
      }
      regionFinish(S)
    }
    skatoBatch <<- list()
    skatoBatchBytes <<- 0
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

        # chrom = region$chrom
        print(paste0(
          "Analyzing Region ",
          regionName,
          " (",
          i - 1,
          "/",
          nRegions,
          ")."
        ))
        # Everything the region's later stages read, so that a batch of
        # regions can be held while the pool evaluates their cells (step 5).
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
          # score transform of a constant column coincides (192 sets -> 12),
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
        S <- list(
          regionName = regionName,
          i = i,
          SNP = SNP,
          WEIGHT = WEIGHT,
          SCORE = SCORE,
          ISLOF = ISLOF,
          annolistsub = annolistsub,
          annoIndicatorMat = region$annoIndicatorMat,
          region = region,
          cctpvalK = rep(NA_real_, nModels),
          pvalRegionK = vector("list", nModels),
          resultDFListK = vector("list", nModels) # the per-set carry-forward cache, per model
        )
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
          skatoBatch[[length(skatoBatch) + 1L]] <- S
          skatoBatchBytes <- skatoBatchBytes + as.numeric(utils::object.size(S$outLists))
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
      if (useSkatoBatch && length(skatoBatch) > 0 &&
          (mth == numberRegionsInChunk || skatoBatchBytes > skatoBatchBudget)) {
        flushSkatoBatch()
      }


      # output
      if (mth == numberRegionsInChunk) {
        message1 <- "This is the output index file for SAIGE package to record the end point in case users want to restart the analysis. Please do not modify this file."
        message2 <- "This is a Region level analysis."
        message3 <- paste("nEachChunk =", nEachChunk)
        message4 <- paste("Have completed the analysis of chunk", indexChunk)
        message5 <- "Have completed the analyses of all chunks."
        # n1 = length(Output)
        # n2 = length(OutputFile)
        cat("write to output\n")
        # cat("n1 is ", n1, "\n")
        # cat("n2 is ", n2, "\n")
        for (k in seq_len(nModels)) {
        if (regionTestType != "BURDEN") {
          # The truncate-vs-append decision is untouched; only the HEADER moves,
          # to the first chunk that actually HAS rows (LEDGER #59).
          if (Start) {
            if (!is.null(pval.Region.all[[k]])) {
              fwrite(
                pval.Region.all[[k]],
                OutputFiles[k],
                quote = FALSE,
                sep = "\t",
                append = FALSE,
                col.names = TRUE,
                row.names = FALSE,
                na = "NA"
              )
              iswriteRegion[k] <- TRUE
            }
          } else {
            if (!is.null(pval.Region.all[[k]])) {
              fwrite(
                pval.Region.all[[k]],
                OutputFiles[k],
                quote = FALSE,
                sep = "\t",
                append = TRUE,
                col.names = !iswriteRegion[k],
                row.names = FALSE,
                na = "NA"
              )
              iswriteRegion[k] <- TRUE
            }
          }
        }
        if (is_output_markerList_in_groupTest) {
          markerListFile <- paste0(OutputFiles[k], ".markerList.txt")
          if (Start) {
            if (!is.null(Output_MarkerList.all)) {
              fwrite(
                Output_MarkerList.all,
                markerListFile,
                quote = FALSE,
                sep = "\t",
                append = FALSE,
                col.names = TRUE,
                row.names = FALSE,
                na = "NA"
              )
              iswriteMarkerList[k] <- TRUE
            }
          } else {
            if (!is.null(Output_MarkerList.all)) {
              if (iswriteMarkerList[k]) {
                fwrite(
                  Output_MarkerList.all,
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
                  Output_MarkerList.all,
                  markerListFile,
                  quote = FALSE,
                  sep = "\t",
                  append = FALSE,
                  col.names = TRUE,
                  row.names = FALSE,
                  na = "NA"
                )
              }
              iswriteMarkerList[k] <- TRUE
            }
          }
        }
        if (Start) {
          write.table(
            c(message1, message2, message3),
            OutputFileIndexes[k],
            quote = FALSE,
            sep = "\t",
            append = FALSE,
            col.names = FALSE,
            row.names = FALSE
          )
        }
        write.table(
          message4,
          OutputFileIndexes[k],
          quote = FALSE,
          sep = "\t",
          append = TRUE,
          col.names = FALSE,
          row.names = FALSE
        )

        if (End) {
          write.table(
            message5,
            OutputFileIndexes[k],
            quote = FALSE,
            sep = "\t",
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE
          )
        }
        } # for k (models)
        # }#if(FALSE)

        pval.Region.all <- vector("list", nModels)
        OutList.all <- NULL
        Output_MarkerList.all <- NULL
        cth_chunk_to_output <- cth_chunk_to_output + 1
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

  message <- paste0(
    "Analysis done! The results have been saved to '",
    paste(OutputFiles, collapse = "', '"),
    "' and '",
    paste(paste0(OutputFiles, ".markerInfo"), collapse = "', '"),
    "'."
  )

  # LEDGER #57: how many cells the SKAT kernel failed on (their rows are NA),
  # and a sidecar naming each with its message when that count is not zero.
  # With several models the count is over all of them and the sidecar goes
  # beside the FIRST model's file (the per-cell record has no model column).
  skatFailureSummary(OutputFile)
  skatoDedupeSummary()

  message <- "Analysis done!"
  message <- paste0(
    message,
    " The set-based tests results have been saved to '",
    paste(OutputFiles, collapse = "', '"),
    "'."
  )
  if (is_output_markerList_in_groupTest) {
    message <- paste0(
      message,
      " The marker lists have been saved to '",
      OutputFile,
      ".markerList.txt'."
    )
  }
  if (is_single_in_groupTest) {
    message <- paste0(
      message,
      " The single-variant association tests results have been saved to '",
      OutputFile,
      ".singleAssoc.txt'."
    )
  }

  message
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
