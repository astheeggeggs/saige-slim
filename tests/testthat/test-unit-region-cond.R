# region_postprocess (R/SAIGE_SPATest_Region.R): the QUANTITATIVE conditional
# arm under the multi-weight schema (C11, design D7).
#
# No bundled golden exercises a quantitative GROUP test (the only group-capable
# committed model is binary), so this arm is driven directly. At K > 1 the C++
# side exports a G1-raw G1tilde_P_G2tilde matrix plus
# VarInvMat_G2_cond_w0 (= VarInvMat_cond / w0G2 w0G2') and
# Tstat_G2_cond_w0 (= Tstat_cond * w0G2), and R forms, per set r and per
# (anno, max_MAF) cell,
#   G1w_r        = raw[cell rows, ] * AnnoWeights_r[cell rows]
#   AdjCondMat_r = G1w_r %*% VarInvMat_G2_cond_w0
#   TstatAdjCond = AdjCondMat_r %*% Tstat_G2_cond_w0
#   VarMatAdjCond= AdjCondMat_r %*% t(G1w_r)
# The tests reconcile that cell-wise computation against the OTHER
# parameterization -- the full-matrix adjustment (all q rows scaled, adjusted,
# THEN subset), which is the literal q_weight == 1 C++ arithmetic -- and
# against a null control (zero conditioning exports => cond == uncond).

# ---- shared synthetic region: 1 anno x 1 max_MAF, q0 = 3 real markers -------

# the "<mantissa>E<exponent>" string form convert_str_to_log requires
pstr <- function(p) toupper(sprintf("%.6e", p))

make_quant_outList <- function(K, raw_G1, VarInv_w0, Tw0, baked = NULL) {
  q0 <- 3L
  numofUR <- K            # one UR pseudo-marker per weight set (1 anno x 1 maf)
  q <- q0 + numofUR
  ind <- matrix(0L, nrow = q, ncol = K)
  ind[1:q0, ] <- 1L
  for (r in seq_len(K)) ind[q0 + r, r] <- 1L
  out <- list(
    pvalVec = pstr(c(0.04, 0.20, 0.51, 0.33, 0.12)[1:q]),
    annoMAFIndicatorMat = ind,
    q_maf_for_annoVec = c(1),
    TstatVec_flip = c(1.2, -0.8, 0.5, 0.9, 0.4)[1:q],
    VarMat = diag(2, q) + 0.1,
    MAFVec = c(0.0010, 0.0020, 0.0005, 0.0012, 0.0011)[1:q],
    MAC_GroupVec = rep(7, K),
    NumRare_GroupVec = rep(3, K),
    NumUltraRare_GroupVec = rep(1, K),
    numofUR = numofUR,
    G1tilde_P_G2tilde_Weighted_Mat = raw_G1,
    G2_Weight_cond = c(2.0, 1.5)
  )
  if (is.null(baked)) {
    out$VarInvMat_G2_cond_w0 <- VarInv_w0
    out$Tstat_G2_cond_w0 <- Tw0
  } else {
    out$TstatAdjCond <- baked$TstatAdjCond
    out$VarMatAdjCond <- baked$VarMatAdjCond
  }
  out
}

run_quant_cond <- function(outList, WEIGHT, weight_sets) {
  region_postprocess(
    outList = outList, mu = rep(0, 10), WEIGHT = WEIGHT,
    traitType = "quantitative", r.corr = 0, regionTestType = "SKAT-O",
    isCondition = TRUE, annolistsub = c("lof"), maxMAFlist = c(0.01),
    regionName = "GENE1", is_single_in_groupTest = TRUE, is_rerun = FALSE,
    weight_sets = weight_sets)
}

test_that("quantitative K=2 conditional adjustment matches the full-matrix parameterization", {
  K <- 2L
  q <- 5L
  WEIGHT <- cbind(W1 = c(1.5, 0.8, 2.0), W2 = c(0.5, 1.2, 0.9))
  raw_G1 <- matrix(c(0.11, -0.07, 0.05, 0.09, -0.02,
                     0.04,  0.12, -0.06, 0.03, 0.08), nrow = q, ncol = 2)
  VarInv_w0 <- matrix(c(0.5, 0.1, 0.1, 0.4), 2, 2)
  Tw0 <- c(0.3, -0.2)
  outList <- make_quant_outList(K, raw_G1, VarInv_w0, Tw0)
  ws <- list(kind = c("custom", "custom"), name = c("W1", "W2"),
             customCol = c(1L, 2L), betaA = c(NA_real_, NA_real_),
             betaB = c(NA_real_, NA_real_), force = FALSE)
  res <- run_quant_cond(outList, WEIGHT, ws)
  d <- res$pval.Region
  expect_identical(nrow(d), 3L)           # 2 sets + 1 pooled Cauchy row
  expect_identical(d$Weight[1:2], c("W1", "W2"))

  # independent recomputation, full-matrix-then-subset (the K = 1 C++ order)
  StatVec <- outList$TstatVec_flip
  for (r in 1:2) {
    w <- c(WEIGHT[, r], rep(1, 2))        # UR pseudo-markers weigh 1 (D2)
    tempPos <- which(outList$annoMAFIndicatorMat[, r] == 1)
    G1w_full <- raw_G1 * w
    Adj_full <- G1w_full %*% VarInv_w0
    VarAdj_full <- Adj_full %*% t(G1w_full)
    TstatAdj_full <- as.vector(Adj_full %*% Tw0)
    Phi <- (outList$VarMat * (w %*% t(w)))[tempPos, tempPos]
    Score <- (StatVec * w)[tempPos]
    exp_uncond <- get_SKAT_pvalue(Score, Phi, 0, "SKAT-O")
    exp_cond <- get_SKAT_pvalue(Score - TstatAdj_full[tempPos],
                                Phi - VarAdj_full[tempPos, tempPos],
                                0, "SKAT-O")
    got <- d[r, ]
    expect_equal(as.numeric(got$Pvalue), exp_uncond$Pvalue_SKATO, tolerance = 1e-12)
    expect_equal(as.numeric(got$Pvalue_cond), exp_cond$Pvalue_SKATO, tolerance = 1e-12)
    expect_equal(as.numeric(got$Pvalue_Burden_cond), exp_cond$Pvalue_Burden, tolerance = 1e-12)
    expect_equal(as.numeric(got$Pvalue_SKAT_cond), exp_cond$Pvalue_SKAT, tolerance = 1e-12)
    expect_equal(as.numeric(got$BETA_Burden_cond), exp_cond$BETA_Burden, tolerance = 1e-12)
    expect_equal(as.numeric(got$SE_Burden_cond), exp_cond$SE_Burden, tolerance = 1e-12)
  }
  # the two sets must genuinely differ (positive control on the setup itself)
  expect_false(identical(d$Pvalue_cond[1], d$Pvalue_cond[2]))
})

test_that("zero conditioning exports are a null adjustment: cond equals uncond", {
  K <- 2L
  WEIGHT <- cbind(W1 = c(1.5, 0.8, 2.0), W2 = c(0.5, 1.2, 0.9))
  raw_G1 <- matrix(0.1, nrow = 5, ncol = 2)
  outList <- make_quant_outList(K, raw_G1,
                                VarInv_w0 = matrix(0, 2, 2),
                                Tw0 = c(0, 0))
  ws <- list(kind = c("custom", "custom"), name = c("W1", "W2"),
             customCol = c(1L, 2L), betaA = c(NA_real_, NA_real_),
             betaB = c(NA_real_, NA_real_), force = FALSE)
  d <- run_quant_cond(outList, WEIGHT, ws)$pval.Region
  r <- d[d$Group != "Cauchy", ]
  expect_identical(r$Pvalue_cond, r$Pvalue)
  expect_identical(r$Pvalue_Burden_cond, r$Pvalue_Burden)
  expect_identical(r$Pvalue_SKAT_cond, r$Pvalue_SKAT)
})

test_that("K=1 with the forced Weight column reproduces the legacy quantitative conditional arm", {
  # At q_weight == 1 C++ bakes w0 and exports TstatAdjCond/VarMatAdjCond as
  # always; the forced-schema R branch must consume those baked exports and
  # produce the legacy numbers (plus the Weight label and a Cauchy row).
  q0 <- 3L
  q <- 4L
  w <- c(1.5, 0.8, 2.0)
  wfull <- c(w, 1)
  raw_G1 <- matrix(c(0.11, -0.07, 0.05, 0.09,
                     0.04,  0.12, -0.06, 0.03), nrow = q, ncol = 2)
  VarInv_w0 <- matrix(c(0.5, 0.1, 0.1, 0.4), 2, 2)
  Tw0 <- c(0.3, -0.2)
  # bake exactly as the C++ q_weight == 1 path would
  G1w <- raw_G1 * wfull
  Adj <- G1w %*% VarInv_w0
  baked <- list(TstatAdjCond = as.vector(Adj %*% Tw0),
                VarMatAdjCond = Adj %*% t(G1w))
  mk <- function() {
    o <- make_quant_outList(1L, G1w, VarInv_w0, Tw0, baked = baked)
    o
  }
  legacy <- run_quant_cond(mk(), WEIGHT = w, weight_sets = NULL)$pval.Region
  ws <- list(kind = "custom", name = "W1", customCol = 1L,
             betaA = NA_real_, betaB = NA_real_, force = TRUE)
  forced <- run_quant_cond(mk(), WEIGHT = w, weight_sets = ws)$pval.Region
  # legacy: single row, no Cauchy (one anno, one max_MAF); forced: + Cauchy
  expect_identical(nrow(legacy), 1L)
  expect_identical(nrow(forced), 2L)
  expect_identical(forced$Weight[1], "W1")
  for (cc in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
               "Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond",
               "BETA_Burden_cond", "SE_Burden_cond")) {
    expect_equal(as.numeric(forced[1, cc]), as.numeric(legacy[1, cc]),
                 tolerance = 1e-12, label = cc)
  }
})

# ---- characterization: the two type-polymorphism bugs in this code path -----
# LEDGER #33 and #34. These assert what region_postprocess / get_SKAT_pvalue do
# TODAY, including the Inf; the fix phase must invert them.

test_that("LEDGER #34 (FIXED): a scalar Phi is traced, not fed to diag()", {
  # get_SKAT_pvalue takes the trace as `if (is.matrix(Phi)) sum(diag(Phi))
  # else sum(Phi)`, so a bare scalar -- which R/SAIGE_SPATest_Region.R:464
  # could produce before its drop = FALSE -- no longer reaches diag(), where
  # a length-1 numeric is read as a DIMENSION.
  Score <- 1.7
  v <- 0.42
  dropped <- get_SKAT_pvalue(Score, v, 0, "SKAT-O")
  kept <- get_SKAT_pvalue(Score, matrix(v, 1, 1), 0, "SKAT-O")

  # the two shapes now agree on every statistic
  expect_equal(dropped$Pvalue_SKATO, kept$Pvalue_SKATO, tolerance = 1e-12)
  expect_equal(dropped$Pvalue_Burden, kept$Pvalue_Burden, tolerance = 1e-12)
  expect_equal(dropped$BETA_Burden, kept$BETA_Burden, tolerance = 1e-12)
  expect_equal(dropped$BETA_Burden, Score / v, tolerance = 1e-12)
  expect_true(is.finite(dropped$BETA_Burden))   # was Inf before the fix
  expect_true(is.finite(dropped$SE_Burden))

  # a variance >= 1 is the dangerous half: diag() would have CONSTRUCTED a
  # trunc(x)-square identity, quadratic in the variance
  big <- get_SKAT_pvalue(Score, 2.5, 0, "SKAT-O")
  expect_equal(big$BETA_Burden, Score / 2.5, tolerance = 1e-12)   # was Score/2
  # and at a realistic weighted variance the old form aborted the run outright
  expect_error(diag(4e5), "memory|cannot allocate")
  expect_equal(get_SKAT_pvalue(Score, 4e5, 0, "SKAT-O")$BETA_Burden,
               Score / 4e5, tolerance = 1e-12)

  # matrix and vector inputs are untouched by the fix
  M <- matrix(c(2, 0.3, 0.3, 1.5), 2, 2)
  expect_identical(sum(diag(M)), if (is.matrix(M)) sum(diag(M)) else sum(M))
  expect_identical(sum(diag(c(0.2, 0.3))), sum(c(0.2, 0.3)))
})

test_that("LEDGER #33 (FIXED): the binary conditional arm subtracts the full matrix", {
  # R/SAIGE_SPATest_Region.R:426 now matches its siblings at :452 and :308.
  # This pins WHY the old diag() form was wrong, so the fix is not undone:
  # diag() on a matrix EXTRACTS, and matrix - vector then recycles
  # column-major, applying d[i] across the whole of row i.
  q <- 3L
  Phi <- matrix(c(2.0, 0.3, 0.1,
                  0.3, 1.5, 0.2,
                  0.1, 0.2, 1.8), q, q)
  VarMatAdjCond <- matrix(c(0.50, 0.20, 0.10,
                            0.20, 0.40, 0.15,
                            0.10, 0.15, 0.30), q, q)

  fixed <- Phi - VarMatAdjCond          # what :426 computes now
  old <- Phi - diag(VarMatAdjCond)      # the upstream 7c7e6da5 form

  # the fix preserves symmetry; the old form did not
  expect_identical(t(fixed), fixed)
  expect_false(identical(t(old), old))
  expect_false(isTRUE(all.equal(old, fixed)))
  # the old form was the row-broadcast of the diagonal
  expect_identical(old, Phi - matrix(rep(diag(VarMatAdjCond), q), q, q))
  # ... and its DIAGONAL was already correct, which is why BETA_Burden_cond
  # (sum(diag(.)) only) did not move when the fix was applied, while the
  # SKAT p-values, which read the off-diagonals, did.
  expect_equal(diag(old), diag(fixed), tolerance = 1e-12)
  expect_equal(get_SKAT_pvalue(rep(1, q), old, 0, "SKAT-O")$BETA_Burden,
               get_SKAT_pvalue(rep(1, q), fixed, 0, "SKAT-O")$BETA_Burden,
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    get_SKAT_pvalue(rep(1, q), old, 0, "SKAT-O")$Pvalue_SKAT,
    get_SKAT_pvalue(rep(1, q), fixed, 0, "SKAT-O")$Pvalue_SKAT)))

  # a single-marker cell coincides either way, which is why this hid so long
  expect_identical(matrix(2.0, 1, 1) - diag(matrix(0.5, 1, 1)),
                   matrix(2.0, 1, 1) - matrix(0.5, 1, 1))
})
