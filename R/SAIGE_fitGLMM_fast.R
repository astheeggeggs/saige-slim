# Run iterations to get converged alpha and eta
Get_Coef <- function(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, tolPCG, maxiter, verbose = FALSE, tol.coef = 0.1) {
  # tol.coef = 0.1
  # eta = eta0

  mu <- family$linkinv(eta0)
  mu.eta <- family$mu.eta(eta0)
  Y <- eta0 - offset + (y - mu) / mu.eta
  sqrtW <- mu.eta / sqrt(family$variance(mu))
  W <- sqrtW^2


  for (i in 1:maxiter) {
    re.coef <- getCoefficients(Y, X, W, tau, maxiter = maxiterPCG, tol = tolPCG)
    alpha <- re.coef$alpha
    eta <- re.coef$eta + offset

    if (verbose) {
      cat("Tau:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }

    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu) / mu.eta
    sqrtW <- mu.eta / sqrt(family$variance(mu))
    W <- sqrtW^2

    if (max(abs(alpha - alpha0) / (abs(alpha) + abs(alpha0) + tol.coef)) < tol.coef) {
      break
    }
    alpha0 <- alpha
  }

  re <- list(Y = Y, alpha = alpha, eta = eta, W = W, cov = re.coef$cov, sqrtW = sqrtW, Sigma_iY = re.coef$Sigma_iY, Sigma_iX = re.coef$Sigma_iX, mu = mu)
}


# Functon to get working vector and fixed & random coefficients when LOCO is TRUE
# getCoefficients_LOCO is used
# Run iterations to get converged alpha and eta
Get_Coef_LOCO <- function(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, tolPCG, maxiter, verbose = FALSE, tol.coef = 0.1) {
  # tol.coef = 0.1
  mu <- family$linkinv(eta0)
  mu.eta <- family$mu.eta(eta0)
  Y <- eta0 - offset + (y - mu) / mu.eta
  sqrtW <- mu.eta / sqrt(family$variance(mu))
  W <- sqrtW^2

  for (i in 1:maxiter) {
    re.coef <- getCoefficients_LOCO(Y, X, W, tau, maxiter = maxiterPCG, tol = tolPCG)
    alpha <- re.coef$alpha
    eta <- re.coef$eta + offset

    if (verbose) {
      cat("Tau:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }

    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    Y <- eta - offset + (y - mu) / mu.eta
    sqrtW <- mu.eta / sqrt(family$variance(mu))
    W <- sqrtW^2

    if (max(abs(alpha - alpha0) / (abs(alpha) + abs(alpha0) + tol.coef)) < tol.coef) break
    alpha0 <- alpha
  }

  re <- list(Y = Y, alpha = alpha, eta = eta, W = W, cov = re.coef$cov, sqrtW = sqrtW, Sigma_iY = re.coef$Sigma_iY, Sigma_iX = re.coef$Sigma_iX, mu = mu)
}


# Fits the null glmm for binary traits
glmmkin.ai_PCG_Rcpp_Binary <- function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau = c(0, 0), fixtau = c(0, 0), maxiter = 20, tol = 0.02, tolCoef = 1e-5, verbose = TRUE, nrun = 30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne, isLowMemLOCO, genoIn = NULL, ingestGates = NULL) {
  # Fits the null generalized linear mixed model for a binary trait
  # Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the glmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  # Returns:
  #  model output for the null glmm

  t_begin <- proc.time()
  print(t_begin)
  subSampleInGeno <- subPheno$IndexGeno
  if (verbose) {
    print("Start reading genotype plink file here")
  }


  if (is.null(genoIn)) genoIn <- genoInputPlink(bedFile, bimFile, famFile)
  # ensureGenoLoaded, not setgenoFromInput: with several phenotypes over the
  # same complete-case sample set this ingest is done ONCE and reused. It
  # reduces to exactly the old call whenever nothing is loaded, which is every
  # single-phenotype run. See R/SAIGE_fitGLMM_multi.R.
  ensureGenoLoaded(
    genoIn, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk,
    isDiagofKinSetAsOne, ingestGates
  )
  if (verbose) {
    print("Genotype reading is done")
  }

  if (LOCO) {
    MsubIndVec <- getQCdMarkerIndex()
    chrVec <- readVariantChromFromInput(genoIn)
    chrVec <- chrVec[which(MsubIndVec == TRUE)]
    updatechrList <- updateChrStartEndIndexVec(chrVec)
    LOCO <- updatechrList$LOCO
    chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
    chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
  }

  y <- fit0$y
  n <- length(y)
  X <- model.matrix(fit0)

  offset <- fit0$offset
  if (is.null(offset)) {
    offset <- rep(0, n)
  }

  family <- fit0$family
  eta <- fit0$linear.predictors

  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu) / mu.eta

  # mu = fit0$fitted.values
  # mu.eta = family$mu.eta(eta)
  # Y = eta - offset + (y - mu)/mu.eta
  alpha0 <- fit0$coef
  eta0 <- eta

  if (family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }

  # change, use 0.5 as a default value, and use Get_Coef before getAIScore
  q <- 1

  if (tauInit[fixtau == 0] == 0) {
    tau[fixtau == 0] <- 0
  } else {
    tau[fixtau == 0] <- tauInit[fixtau == 0]
  }
  cat("inital tau is ", tau, "\n")
  tau0 <- tau

  # LEDGER #125: the inner (alpha, b) pass converges to tolCoef, not to the outer
  # tol: at tol = 0.02 the REML score at a fixed tau was four times noisier than
  # the exact score, which is what the #122 guard was built to survive.
  re.coef <- Get_Coef(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, tol.coef = tolCoef)
  re <- getAIScore(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG = tolPCG, traceCVcutoff = traceCVcutoff)


  tau[2] <- max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace) / n)

  if (verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }

  guard <- tauGuardInit()   # LEDGER #122
  critPrev <- FALSE
  didConverge <- FALSE   # LEDGER #123: set by the two exits that mean it
  for (i in seq_len(maxiter)) {
    if (verbose) cat("\nIteration ", i, tau, ":\n")
    alpha0 <- re.coef$alpha
    tau0 <- tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 <- eta
    # use Get_Coef before getAIScore
    t_begin_Get_Coef <- proc.time()
    re.coef <- Get_Coef(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, tol.coef = tolCoef)
    t_end_Get_Coef <- proc.time()
    # cat("t_end_Get_Coef - t_begin_Get_Coef\n")
    cat("Updating fix effect coeffcieints took\n")
    print(t_end_Get_Coef - t_begin_Get_Coef)

    # LEDGER #122: the score terms come from getAIScore, the float Newton step
    # from tauProposalBinary (fitglmmaiRPCG's arithmetic, unchanged), and the
    # guard decides how much of that step is taken.
    re <- getAIScore(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG = tolPCG, traceCVcutoff = traceCVcutoff)
    score <- re$YPAPY - re$Trace
    prop <- tauProposalBinary(tau0[2], re$YPAPY, re$Trace, re$AI, tol)
    if (tau0[2] == 0 && prop[["step"]] == 0) {
      # zero is an iterate, not a decision: an arrival AT 0 from a positive
      # tau0 is evaluated again here. But once the iterate IS 0 and the
      # floored step from 0 is 0 again, another iteration has identical
      # inputs and cannot move it (LEDGER #161: with the exit keyed on
      # `score <= 0` a positive score whose Newton step was below tol --
      # floored to 0 by tauStepBinary's GMMAT rule -- ran out --maxiter at
      # (1, 0)). Two ways to be here, and the log says which. The post-loop
      # refit reads the inner pass's state, so hand it over as the step
      # below would have.
      cov <- re.coef$cov
      alpha <- re.coef$alpha
      eta <- re.coef$eta
      Y <- re.coef$Y
      mu <- re.coef$mu
      if (score <= 0) {
        cat("tau2 = 0 is the optimum: the REML score at 0 is ", score, "\n")
      } else {
        cat("tau2 = 0 by the floor: the REML score at 0 is ", score,
            " and the Newton step from 0, ", prop[["raw"]],
            ", is below tol = ", tol, ", so the estimate is set to 0 (LEDGER #161)\n")
      }
      didConverge <- TRUE
      break
    }
    sg <- tauGuardApply(guard, tau0[2], prop[["raw"]], prop[["step"]])
    guard <- sg$state
    fit <- list(tau = c(tau0[1], sg$p))
    t_end_fitglmmaiRPCG <- proc.time()
    # cat("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
    cat("Updating variance component estimate took\n")
    print(t_end_fitglmmaiRPCG - t_end_Get_Coef)

    tau <- as.numeric(fit$tau)
    cov <- re.coef$cov
    alpha <- re.coef$alpha
    eta <- re.coef$eta
    Y <- re.coef$Y
    mu <- re.coef$mu

    print(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol))
    cat("tau: ", tau, "\n")
    cat("tau0: ", tau0, "\n")


    # Convergence is judged on the UNDAMPED proposal (how far the score is from
    # zero); never on the first iteration, whose proposal comes from an inner
    # state one Get_Coef call away from the GLM start (a cold start next to the
    # answer "converged" there in one step); and, once the guard has damped a
    # step in this fit, on two consecutive iterations: a single pair within tol
    # is what let a cycle pass (#122). Keyed on a damped step, not on arming: a
    # fit with two reversals that converges on its own keeps today's criterion
    # and today's iteration count (the dense-GRM golden has exactly that).
    tauRaw <- c(tau0[1], prop[["raw"]])
    critNow <- max(abs(tauRaw - tau0) / (abs(tauRaw) + abs(tau0) + tol)) < tol
    if (critNow && i > 1 && (!guard$fired || critPrev)) {
      didConverge <- TRUE
      break
    }
    critPrev <- critNow

    if (max(tau) > tol^(-2)) {
      warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      i <- maxiter
      break
    }
  }

  if (verbose) cat("\nFinal ", tau, ":\n")

  # added these steps after tau is estimated 04-14-2018
  re.coef <- Get_Coef(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, tol.coef = tolCoef)
  cov <- re.coef$cov
  alpha <- re.coef$alpha
  eta <- re.coef$eta

  Y <- re.coef$Y
  mu <- re.coef$mu

  # LEDGER #123: the flag records the exit, not the iteration count (a fit that
  # met the criterion on its last allowed iteration was flagged unconverged).
  converged <- didConverge
  res <- y - mu

  mu2 <- mu * (1 - mu)
  if (!isCovariateOffset) {
    obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, X)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = X, traitType = "binary", isCovariateOffset = isCovariateOffset)
  } else {
    obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, Xorig)
    glmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = Xorig, traitType = "binary", isCovariateOffset = isCovariateOffset)
  }
  # LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme

  glmmResult$LOCO <- LOCO
  t_end_null <- proc.time()
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)


  if (!isLowMemLOCO & LOCO) {
    set_Diagof_StdGeno_LOCO()
    glmmResult$LOCOResult <- list()

    for (j in 1:22) {
      startIndex <- chromosomeStartIndexVec[j]
      endIndex <- chromosomeEndIndexVec[j]
      if (!is.na(startIndex) && !is.na(endIndex)) {
        cat("leave chromosome ", j, " out\n")
        setStartEndIndex(startIndex, endIndex, j - 1)
        t_begin_Get_Coef_LOCO <- proc.time()
        re.coef_LOCO <- Get_Coef_LOCO(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, tol.coef = tolCoef)
        t_end_Get_Coef_LOCO <- proc.time()
        # cat("t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO\n")
        cat("Updating fix effect coeffcieints and model residuals took\n")
        print(t_end_Get_Coef_LOCO - t_begin_Get_Coef_LOCO)
        cov <- re.coef_LOCO$cov
        alpha <- re.coef_LOCO$alpha
        eta <- re.coef_LOCO$eta
        Y <- re.coef_LOCO$Y
        mu <- re.coef_LOCO$mu
        mu2 <- mu * (1 - mu)

        res <- y - mu

        if (!is.null(out.transform) & is.null(fit0$offset)) {
          coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
        } else {
          coef.alpha <- alpha
        }


        if (!isCovariateOffset) {
          obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, X)
        } else {
          obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, Xorig)
        }
        glmmResult$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK, alpha0 = alpha)

        # obj.noK = ScoreTest_NULL_Model(mu, mu2, y, X)
        # obj.noK = ScoreTest_NULL_Model_binary(mu, y, X)
        # glmmResult$LOCOResult[[j]] = list(isLOCO = TRUE, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov)
      } else {
        glmmResult$LOCOResult[[j]] <- list(isLOCO = FALSE)
      }
    }
  }

  if (isLowMemLOCO & LOCO) {
    glmmResult$chromosomeStartIndexVec <- chromosomeStartIndexVec
    glmmResult$chromosomeEndIndexVec <- chromosomeEndIndexVec
  }
  glmmResult
}


# Fits the null glmm for a quantitative trait
glmmkin.ai_PCG_Rcpp_Quantitative <- function(bedFile, bimFile, famFile, Xorig, isCovariateOffset, fit0, tau = c(0, 0), fixtau = c(0, 0), maxiter = 20, tol = 0.02, tolCoef = 1e-5, verbose = TRUE, nrun = 30, tolPCG = 1e-5, maxiterPCG = 500, subPheno, indicatorGenoSamplesWithPheno, obj.noK, out.transform, tauInit, memoryChunk, LOCO, chromosomeStartIndexVec, chromosomeEndIndexVec, traceCVcutoff, isCovariateTransform, isDiagofKinSetAsOne, isLowMemLOCO, genoIn = NULL, ingestGates = NULL) {
  # Fits the null linear mixed model for a quantitative trait
  # Args:
  #  genofile: string. Plink file for the M1 markers to be used to construct the genetic relationship matrix
  #  fit0: glm model. linear model output (with no sample relatedness accounted for), family=gaussian(link = "identity")
  #  tau: vector for iniial values for the variance component parameter estimates
  #  fixtau: vector for fixed tau values
  #  maxiter: maximum iterations to fit the lmm model
  #  tol: tolerance for tau estimating to converge
  #  verbose: whether outputting messages in the process of model fitting
  #  nrun: integer. Number of random vectors used for trace estimation
  #  tolPCG: tolerance for PCG to converge
  #  maxiterPCG: maximum iterations for PCG to converge
  #  subPheno: data set with samples having non-missing phenotypes and non-missing genotypes (for M1 markers)
  #  obj.noK: model output from the SPAtest::ScoreTest_wSaddleApprox_NULL_Model
  #  out.transform: output from the function Covariate_Transform
  #  tauInit: vector for iniial values for the variance component parameter estimates
  #  memoryChunk: integer or float. The size (Gb) for each memory chunk
  #  LOCO:logical. Whether to apply the leave-one-chromosome-out (LOCO) option.
  #  chromosomeStartIndexVec: integer vector of length 22. Contains start indices for each chromosome, starting from 0
  #  chromosomeEndIndexVec: integer vector of length. Contains end indices for each chromosome
  #  traceCVcutoff: threshold for the coefficient of variation for trace estimation
  # Returns:
  #  model output for the null lmm

  t_begin <- proc.time()
  print(t_begin)
  subSampleInGeno <- subPheno$IndexGeno
  if (verbose) {
    print("Start reading genotype plink file here")
  }

  if (is.null(genoIn)) genoIn <- genoInputPlink(bedFile, bimFile, famFile)
  # ensureGenoLoaded, not setgenoFromInput: with several phenotypes over the
  # same complete-case sample set this ingest is done ONCE and reused. It
  # reduces to exactly the old call whenever nothing is loaded, which is every
  # single-phenotype run. See R/SAIGE_fitGLMM_multi.R.
  ensureGenoLoaded(
    genoIn, subSampleInGeno, indicatorGenoSamplesWithPheno, memoryChunk,
    isDiagofKinSetAsOne, ingestGates
  )
  if (verbose) {
    print("Genotype reading is done")
  }

  if (LOCO) {
    MsubIndVec <- getQCdMarkerIndex()
    chrVec <- readVariantChromFromInput(genoIn)
    chrVec <- chrVec[which(MsubIndVec == TRUE)]
    updatechrList <- updateChrStartEndIndexVec(chrVec)
    LOCO <- updatechrList$LOCO
    chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
    chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
  }


  if (verbose) {
    print("Genotype reading is done")
  }

  y <- fit0$y
  n <- length(y)
  offset <- fit0$offset
  if (is.null(offset)) offset <- rep(0, n)
  family <- fit0$family
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu) / mu.eta
  sqrtW <- mu.eta / sqrt(fit0$family$variance(mu))
  W <- sqrtW^2
  X <- model.matrix(fit0)

  alpha <- fit0$coef
  alpha0 <- alpha
  eta0 <- eta
  if (verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }

  if (family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }

  q <- 1
  if (sum(tauInit[fixtau == 0]) == 0) {
    # tau[fixtau == 0] = var(Y)/(q+1)
    tau[1] <- 1
    tau[2] <- 0
    # LEDGER #65 FIXED. The guard here used to stop when the variance of Y was
    # below 0.1, with a message that said WARNING, and depended only on the
    # phenotype's unit. fitNULLGLMM
    # now hands this function a phenotype standardised to var(Y) = 1, so the
    # (1, 0) start is always in scale; a direct caller is stopped only when
    # nothing can be fitted at all.
    if (!is.finite(var(Y)) || var(Y) <= 0) {
      stop("The phenotype has zero variance after covariate adjustment; nothing can be fitted\n")
    }
  } else {
    tau[fixtau == 0] <- tauInit[fixtau == 0]
  }

  tau0 <- tau
  cat("initial tau is ", tau, "\n")
  re.coef <- Get_Coef(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter)

  re <- getAIScore_q(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG = tolPCG, traceCVcutoff = traceCVcutoff)
  tau[2] <- max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace[2]) / n)
  tau[1] <- max(0, tau0[1] + tau0[1]^2 * (re$YPA0PY - re$Trace[1]) / n)

  if (verbose) {
    cat("Variance component estimates:\n")
    print(tau)
  }


  for (i in seq_len(maxiter)) {
    W <- sqrtW^2

    if (verbose) cat("\nIteration ", i, ":\n")
    alpha0 <- alpha
    tau0 <- tau
    eta0 <- eta
    #    cat("tau0: ", tau0,"\n")

    re.coef <- Get_Coef(y, X, tau, family, alpha0, eta0, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter)
    fit <- fitglmmaiRPCG_q(re.coef$Y, X, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov, nrun, maxiterPCG, tolPCG, tol = tol, traceCVcutoff = traceCVcutoff)

    #    cat("tau0_after_fit: ", tau0,"\n")
    #    print(fit)

    if ((fit$tau)[1] != 0) {
      tau <- as.numeric(fit$tau)
      cov <- re.coef$cov
      alpha <- re.coef$alpha
      eta <- re.coef$eta

      cat("cov: ", cov, "\n")

      if (verbose) {
        cat("Variance component estimates:\n")
        print(tau)
        cat("Fixed-effect coefficients:\n")
        print(alpha)
      }
      Y <- re.coef$Y
      mu <- re.coef$mu

      # mu = family$linkinv(eta)
      # mu.eta = family$mu.eta(eta)
      # Y = eta - offset + (y - mu)/mu.eta
      # sqrtW = mu.eta/sqrt(family$variance(mu))

      # if(tau[1]<=0){
      #  stop("ERROR! The first variance component parameter estimate is 0\n")
      # }

      if (tau[1] <= 0 | tau[2] <= 0) break


      if (max(abs(tau - tau0) / (abs(tau) + abs(tau0) + tol)) < tol) break

      #    if(2*max(max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol)), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
      if (max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
        i <- maxiter
        break
      }
    } else { # if(tau[1]!=0){
      break
    }
  }

  if (verbose) cat("\nFinal ", tau, ":\n")

  t_end_null <- proc.time()
  # print(t_end_null)
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)


  re.coef <- Get_Coef(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter)
  cov <- re.coef$cov
  alpha <- re.coef$alpha
  eta <- re.coef$eta
  Y <- re.coef$Y
  mu <- re.coef$mu

  converged <- ifelse(i < maxiter, TRUE, FALSE)
  res <- y - mu
  mu2 <- rep((1 / (tau[1])), length(res))

  # if(!is.null(out.transform) & is.null(fit0$offset)){
  #  coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  # }else{
  #  coef.alpha = alpha
  # }

  # coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
  if (!isCovariateOffset) {
    obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, X)
    lmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = X, traitType = "quantitative", isCovariateOffset = isCovariateOffset)
  } else {
    obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, Xorig)
    lmmResult <- list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, converged = converged, sampleID = subPheno$IID, obj.noK = obj.noK, y = y, X = Xorig, traitType = "quantitative", isCovariateOffset = isCovariateOffset)
  }
  # lmmResult = list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged, sampleID = subPheno$IID, Sigma_iy = Sigma_iy, Sigma_iX = Sigma_iX, obj.noK=obj.noK, obj.glm.null=fit0, traitType="quantitative")

  # LOCO: estimate fixed effect coefficients, random effects, and residuals for each chromoosme
  lmmResult$LOCO <- LOCO
  if (LOCO & !isLowMemLOCO) {
    lmmResult$LOCOResult <- list()
    set_Diagof_StdGeno_LOCO()
    for (j in 1:22) {
      cat("leave chromosome ", j, " out\n")
      startIndex <- chromosomeStartIndexVec[j]
      endIndex <- chromosomeEndIndexVec[j]
      if (!is.na(startIndex) && !is.na(endIndex)) {
        setStartEndIndex(startIndex, endIndex, j - 1)
        re.coef_LOCO <- Get_Coef_LOCO(y, X, tau, family, alpha, eta, offset, verbose = verbose, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter)

        # re.coef_LOCO=fitglmmaiRPCG_q_LOCO(Y, X, W, tau, nrun, maxiterPCG, tolPCG, tol, traceCVcutoff)
        cov <- re.coef_LOCO$cov
        alpha <- re.coef_LOCO$alpha
        eta <- re.coef_LOCO$eta
        Y <- re.coef_LOCO$Y
        mu <- re.coef_LOCO$mu

        res <- y - mu
        mu2 <- rep((1 / (tau[1])), length(res))
        if (!is.null(out.transform) & is.null(fit0$offset)) {
          # if(!is.null(out.transform)){
          coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
        } else {
          coef.alpha <- alpha
        }
        if (!isCovariateOffset) {
          obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, X)
        } else {
          obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, Xorig)
        }
        # coef.alpha<-Covariate_Transform_Back(alpha, out.transform$Param.transform)
        lmmResult$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK, alpha0 = alpha)
      } else {
        lmmResult$LOCOResult[[j]] <- list(isLOCO = FALSE)
      }
    }
  }

  if (isLowMemLOCO & LOCO) {
    lmmResult$chromosomeStartIndexVec <- chromosomeStartIndexVec
    lmmResult$chromosomeEndIndexVec <- chromosomeEndIndexVec
  }
  lmmResult
  # return(list(theta=tau, coefficients=coef.alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged, sampleID = subPheno$IID, Sigma_iy = Sigma_iy, Sigma_iX = Sigma_iX, obj.noK=obj.noK, obj.glm.null=fit0, traitType="quantitative"))
}


Saddle_Prob_q <- function(q, mu, g, tauVecNew) {
  m1 <- sum(mu * g)
  var2 <- sum(g^2)
  pval.noadj <- pchisq(((q - m1) / tauVecNew[1])^2 / var2, lower.tail = FALSE, df = 1)
  list(p.value = pval.noadj, p.value.NA = NA, Is.converge = NA, p1 = NA, p2 = NA)
}


# ScoreTest_wSaddleApprox_NULL_Model_q=function (formula, tau, data = NULL){
ScoreTest_wSaddleApprox_NULL_Model_q <- function(mu, y, X, tauVec) {
  V <- rep(1 / tauVec[1], length(y))
  res <- y - mu
  n1 <- length(res)
  XV <- t(X * V)
  XVX <- t(X) %*% (X * V)
  XVX_inv <- solve(XVX)
  XXVX_inv <- X1 %*% XVX_inv
  XVX_inv_XV <- XXVX_inv * V
  S_a <- colSums(X * res)

  re <- list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV)
  class(re) <- "SA_NULL"
  re
}


ScoreTest_NULL_Model <- function(mu, mu2, y, X) {
  V <- as.vector(mu2)
  res <- as.vector(y - mu)
  XV <- t(X * V)
  XVX <- t(X) %*% (t(XV))
  XVX_inv <- solve(XVX)
  XXVX_inv <- X %*% XVX_inv
  XVX_inv_XV <- XXVX_inv * V
  S_a <- colSums(X * res)
  re <- list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  class(re) <- "SA_NULL"
  re
}

solveSpMatrixUsingArma <- function(sparseGRMtest) {
  m4 <- gen_sp_v2(sparseGRMtest)
  m4
}


#' Fit the null logistic/linear mixed model and estimate the variance ratios by randomly selected variants
#'
#' @param plinkFile character. Path to plink file to be used for calculating elements of the genetic relationship matrix (GRM). minMAFforGRM can be used to specify the minimum MAF of markers in the plink file to be used for constructing GRM. Genetic markers are also randomly selected from the plink file to estimate the variance ratios
#' @param plink2File character. Prefix for plink2 pgen/pvar/psam files, used in
#'   place of plinkFile. Mutually exclusive with the PLINK input
#' @param pgenFile character. Path to plink2 file (pgen). If plink2File is
#'   specified, 'plink2File'.pgen will be used
#' @param pvarFile character. Path to plink2 file (pvar). If plink2File is
#'   specified, 'plink2File'.pvar will be used
#' @param psamFile character. Path to plink2 file (psam). If plink2File is
#'   specified, 'plink2File'.psam will be used
#' @param phenoFile character. Path to the phenotype file. The file can be either tab or space delimited. The phenotype file has a header and contains at least two columns. One column is for phentoype and the other column is for sample IDs. Additional columns can be included in the phenotype file for covariates in the null model. Please specify the names of the covariates using the argument covarColList and specify categorical covariates using the argument qCovarCol. All categorical covariates must also be included in covarColList.
#' @param phenoCol character. Column name for the phenotype in phenoFile e.g. "CAD"
#' @param traitType character. e.g. "binary" or "quantitative". By default, "binary"
#' @param invNormalize logical. Whether to perform the inverse normalization for the phentoype or not. e.g. TRUE or FALSE. By default, FALSE
#' @param covarColList vector of characters. Covariates to be used in the null model. e.g c("Sex", "Age")
#' @param qCovarCol vector of characters. Categorical covariates to be used in the null model. All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex").
#' @param sampleIDColinphenoFile character. Column name for the sample IDs in the phenotype file. By default "IID"; a column that does not exist stops the run.
#' @param tol numeric. The tolerance for fitting the null model to converge. By default, 0.02.
#' @param tolCoef numeric. The tolerance of the inner (alpha, b) pass of a binary fit at each
#'   variance-component iterate (LEDGER #125: it was tied to tol, which left the REML score
#'   noisier than its own signal). By default, 1e-5. Quantitative fits solve that pass directly
#'   and ignore it.
#' @param maxiter integer. The maximum number of iterations used to fit the null GLMMM. By default, 20.
#' @param tolPCG numeric. The tolerance for PCG to converge. By default, 1e-5.
#' @param maxiterPCG integer. The maximum number of iterations for PCG. By default, 500.
#' @param nThreads integer. Number of threads to be used. By default, 1
#' @param SPAcutoff numeric. The cutoff for the deviation of score test statistics from the mean in the unit of sd to perform SPA. By default, 2.
#' @param numMarkersForVarRatio integer (>0). Minimum number of markers to be used for estimating the variance ratio. By default, 30
#' @param skipModelFitting logical.  Whether to skip fitting the null model and only calculating the variance ratio, By default, FALSE. If TURE, the model file ".rda" is needed
#' @param memoryChunk integer or float. The size (Gb) for each memory chunk. By default, 2
#' @param tauInit vector of numbers. e.g. c(1,1), Initial values for tau. For binary traits, the first element will be always be set to 1. If the tauInit is 0,0, the second element will be 0.5 for binary traits and the initial tau vector for quantitative traits is 1,0. For quantitative traits the whole of step 1 runs on the phenotype divided by the sd of its covariate-adjusted values and the saved model and variance ratios are scaled back (LEDGER #65), so tauInit is given in phenotype units, the convergence and divergence rules are scale-free, and the old stop() on a phenotype variance below 0.1 is gone; only a constant phenotype stops
#' @param LOCO logical or NULL. Whether to apply the leave-one-chromosome-out
#'   (LOCO) option. LOCO requires a full GRM and is REFUSED with
#'   useSparseGRMtoFitNULL = TRUE, which has no per-chromosome structure to
#'   leave out. NULL, the default, resolves to FALSE for a sparse-GRM fit and
#'   TRUE otherwise
#' @param traceCVcutoff numeric. The threshold for coefficient of variation (CV) for the trace estimator to increase nrun. By default, 0.0025
#' @param ratioCVcutoff numeric. The threshold for coefficient of variation (CV) for the variance ratio estimate. If ratioCV > ratioCVcutoff. numMarkersForVarRatio will be increased by 10. By default, 0.001
#' @param outputPrefix character. Path to the output files with prefix.
#' @param outputPrefix_varRatio character. Path to the output variance ratio file with prefix. variace ratios will be output to outputPrefix_varRatio.varianceRatio.txt. If outputPrefix_varRatio is not specified, outputPrefix_varRatio will be the same as the outputPrefix. By default, ""
#' @param IsOverwriteVarianceRatioFile logical. Whether to overwrite the variance ratio file if the file exists. By default, FALSE
#' @param sparseGRMFile character. Path to the pre-calculated sparse GRM file. By default, ""
#' @param sparseGRMSampleIDFile character. Path to the sample ID file for the pre-calculated sparse GRM. No header is included. The order of sample IDs is corresponding to the order of samples in the sparse GRM. By default, ""
#' @param numRandomMarkerforSparseKin integer. number of randomly selected markers (MAF >= 0.01) to be used to identify related samples that are included in the sparse GRM. By default, 2000
#' @param relatednessCutoff float. Applied on LOAD of a sparse GRM (drop0 with this tolerance): entries at or below it are zeroed. By default 0, so the GRM is used exactly as step 0 wrote it; the pairs to keep are decided once, by step 0's own relatednessCutoff.
#' @param cateVarRatioIndexVec vector of integer 0 or 1. The length of cateVarRatioIndexVec is the number of MAC categories for variance ratio estimation. 1 indicates variance ratio in the MAC category is to be estimated, otherwise 0. By default, NULL. If NULL, variance ratios corresponding to all specified MAC categories will be estimated. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMinMACVecExclude vector of float. Exclusive lower bound of MAC for MAC categories; its minimum is also the MAC floor for markers entering the variance-ratio pool (mac >= min). The length equals to the number of MAC categories for variance ratio estimation. By default, c(10,20.5), the same as step 2 and as both CLI drivers. This argument is only activated when isCateVarianceRatio=TRUE
#' @param cateVarRatioMaxMACVecInclude vector of float. Higher bound of MAC for MAC categories. The length equals to the number of MAC categories for variance ratio estimation minus 1. By default, c(20.5). This argument is only activated when isCateVarianceRatio=TRUE
#' @param isCovariateTransform logical. Whether use qr transformation on non-genetic covariates. By default, TRUE
#' @param isDiagofKinSetAsOne logical. Whether to set the diagnal elements in GRM to be 1. By default, FALSE
#' @param useSparseGRMtoFitNULL logical. Whether to use sparse GRM to fit the null model. By default, FALSE
#' @param usePCGwithSparseGRM logical. Whether to use PCG for inverse of the sparse GRM, if useSparseGRMtoFitNULL is TRUE. By default, FALSE
#' @param useSparseGRMforVarRatio logical. Whether to use sparse GRM to estimate the variance Ratios. If TRUE, the variance ratios will be estimated using the full GRM (numerator) and the sparse GRM (denominator). By default, FALSE
#' @param minCovariateCount integer. If binary covariates have a count less than this, they will be excluded from the model to avoid convergence issues. By default, -1 (no covariates will be excluded)
#' @param minCaseCount integer. Binary traits only. Refuse to fit if the final sample set has fewer than this many cases. Below a few hundred cases the step-2 SKAT-O null tail is liberal (measured on this branch: 2.7x the nominal rate at p = 1e-4 with 100 cases, 1.5x at 1e-3 with 500), so a model fitted below the gate produces p-values that are wrong exactly where they would be reported. By default, 100; 0 disables the gate. A fit with fewer than 500 cases is warned about and still run.
#' @param minSampleCount integer. Refuse to fit if the final sample set has fewer than this many samples with a non-missing phenotype. Applies to both trait types. By default, 0 (no gate) -- the quantitative null is calibrated at small N, so this is a power decision rather than a correctness one.
#' @param minCovariateCellCount integer. Refuse to fit if any level of a categorical covariate has fewer than this many cases or controls (fewer than this many samples for a quantitative trait). A level with zero cases is complete separation: the logistic coefficient is infinite and, with isCovariateOffset = TRUE, is frozen into the offset, pinning fitted values at 0 or 1 and collapsing the step-2 score variance. Every covariate named in qCovarCol is checked, as is any covariate that is effectively binary. By default, 10 (the events-per-variable convention); 0 disables the gate. This composes with minCovariateCount, which EXCLUDES a low-count covariate rather than refusing: the exclusion runs first and this gate judges what survives.
#' @param separationCheck character. What to do when the covariate GLM shows separation -- it fails to converge, or a fitted probability reaches the 0/1 boundary. One of "error" (the default), "warn" or "none". Unlike minCovariateCellCount this also catches separation on a continuous covariate, such as an age-by-sex interaction. Binary traits only.
#' @param separationMaxFraction numeric. The fraction of samples whose fitted probability may sit at the 0/1 boundary before separationCheck treats the model as separated. A pinned sample contributes ~0 to the score variance instead of mu*(1 - mu), so it carries no information and a pinned fraction costs that same fraction of the effective sample size; a fraction of 1 is total collapse. Anything above zero but at or below this is reported as a warning and the run continues. By default, 0.02. This is a reasoned default, not a measured one -- no study here has established where p-values actually degrade as the fraction rises.
#' @param convergenceCheck character. What to do when the null fit does not converge -- the AI-REML loop ran out without meeting its criterion, so the variance components and the variance ratio derived from them did not settle. One of "error" (the default), "warn" or "none". Both trait types.
#' @param minMAFforGRM numeric. Minimum MAF for markers (in the Plink file) used for construcing the sparse GRM. By default, 0.01
#' @param includeNonautoMarkersforVarRatio logical. Whether to allow for non-autosomal markers for variance ratio. By default, FALSE
#' @param FemaleOnly logical. Whether to run Step 1 for females only. If TRUE, sexCol and FemaleCode need to be specified. By default, FALSE
#' @param MaleOnly logical. Whether to run Step 1 for males only. If TRUE, sexCol and MaleCode need to be specified. By default, FALSE
#' @param FemaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for females. By default, '1'
#' @param MaleCode character. Values in the column for sex (sexCol) in the phenotype file are used for males. By default, '0'
#' @param sexCol character. Coloumn name for sex in the phenotype file, e.g Sex. By default, ''
#' @param isCovariateOffset logical. TRUE (the default on both the CLI and this function): the covariate effects are estimated once by a fixed-effects GLM and entered into the mixed-model fit as an offset. FALSE: they are re-estimated inside the mixed model.
#' @return a file ended with .rda that contains the glmm model information, a file ended with .varianceRatio.txt that contains the variance ratio values, and a file ended with #markers.SPAOut.txt that contains the SPAGMMAT tests results for the markers used for estimating the variance ratio.
#' @param seed integer. Random seed for the run; also decouples variance-ratio marker selection from model-fitting convergence. By default, 1 (LEDGER #8)
#' @param bedFile character. Path to bed file. If plinkFile is specified, 'plinkFile'.bed is used
#' @param bimFile character. Path to bim file. If plinkFile is specified, 'plinkFile'.bim is used
#' @param famFile character. Path to fam file. If plinkFile is specified, 'plinkFile'.fam is used
#' @param isLowMemLOCO logical. Whether to write the model file by chromosome when LOCO=TRUE. If TRUE, memory use in Step 1 and Step 2 is lower. By default, FALSE
#' @param isCateVarianceRatio logical. Whether to estimate the variance ratio separately for each MAC category, corresponding to cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude. Combines freely with LOCO (nothing forces LOCO=FALSE; the earlier claim that it did was false). By default, FALSE
#' @param maxMissingRateforGRM numeric. Maximum missing rate of the markers used for the GRM. By default, 0.15
#' @param SampleIDIncludeFile character. Path to a file with one column of IDs of the samples to include in the null model fit
#' @param skipVarianceRatioEstimation logical. Whether to skip variance ratio estimation. By default, FALSE
#' @param nrun numeric. Number of runs in the trace estimation. By default, 30
#' @param releaseGenoAtEnd logical. Whether to free the in-memory genotype object when the fit finishes. By default, TRUE; set FALSE only when another call in the same process reuses it
#' @export
fitNULLGLMM <- function(plinkFile = "",
                        plink2File = "",
                        pgenFile = "",
                        pvarFile = "",
                        psamFile = "",
                        seed = 1,
                        bedFile = "",
                        bimFile = "",
                        famFile = "",
                        phenoFile = "",
                        phenoCol = "",
                        traitType = "binary",
                        invNormalize = FALSE,
                        covarColList = NULL,
                        qCovarCol = NULL,
                        sampleIDColinphenoFile = "IID", # LEDGER #151: was ""; the CLI has always passed "IID"
                        tol = 0.02,
                        tolCoef = 1e-5,
                        maxiter = 20,
                        tolPCG = 1e-5,
                        maxiterPCG = 500,
                        nThreads = 1,
                        SPAcutoff = 2,
                        numMarkersForVarRatio = 30,
                        skipModelFitting = FALSE,
                        memoryChunk = 2,
                        tauInit = c(0, 0),
                        LOCO = NULL,
                        isLowMemLOCO = FALSE,
                        traceCVcutoff = 0.0025,
                        ratioCVcutoff = 0.001,
                        outputPrefix = "",
                        outputPrefix_varRatio = "",
                        IsOverwriteVarianceRatioFile = FALSE,
                        sparseGRMFile = "",
                        sparseGRMSampleIDFile = "",
                        numRandomMarkerforSparseKin = 1000,
                        relatednessCutoff = 0, # LEDGER #151: was 0.125; the CLI has always passed 0
                        isCateVarianceRatio = FALSE,
                        cateVarRatioIndexVec = NULL,
                        cateVarRatioMinMACVecExclude = c(10, 20.5),
                        cateVarRatioMaxMACVecInclude = c(20.5),
                        isCovariateTransform = TRUE,
                        isDiagofKinSetAsOne = FALSE,
                        minCovariateCount = -1,
                        minCaseCount = 100,
                        minSampleCount = 0,
                        minCovariateCellCount = 10,
                        separationCheck = "error",
                        separationMaxFraction = 0.02,
                        convergenceCheck = "error",
                        minMAFforGRM = 0.01,
                        maxMissingRateforGRM = 0.15,
                        useSparseGRMtoFitNULL = FALSE,
                        usePCGwithSparseGRM = FALSE,
                        useSparseGRMforVarRatio = FALSE,
                        includeNonautoMarkersforVarRatio = FALSE,
                        sexCol = "",
                        FemaleCode = 1,
                        FemaleOnly = FALSE,
                        MaleCode = 0,
                        MaleOnly = FALSE,
                        SampleIDIncludeFile = "",
                        isCovariateOffset = TRUE, # LEDGER #151: was FALSE; the CLI has always passed TRUE
                        skipVarianceRatioEstimation = FALSE,
                        nrun = 30,
                        # FALSE keeps the loaded genotype store for the next
                        # phenotype in the same sample-set group; the multi
                        # driver is the only caller that passes FALSE, and it
                        # calls releaseGeno() itself at the end.
                        releaseGenoAtEnd = TRUE) {
  ## set up output files
  if (FemaleOnly & MaleOnly) {
    stop("Both FemaleOnly and MaleOnly are TRUE. Please specify only one of them as TRUE to run the sex-specific job\n")
  }

  if (FemaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_FemaleOnly")
    cat(
      "Female-specific model will be fitted. Samples coded as ",
      FemaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  } else if (MaleOnly) {
    outputPrefix <- paste0(outputPrefix, "_MaleOnly")
    cat(
      "Male-specific model will be fitted. Samples coded as ",
      MaleCode, " in the column ", sexCol, " in the phenotype file will be included\n"
    )
  }

  # LOCO is resolved HERE, before its first read at `if (LOCO & isLowMemLOCO)`
  # below. LOCO leaves one chromosome out of the GRM and a sparse GRM has no
  # per-chromosome structure to leave out, so the two are incompatible rather
  # than merely redundant: an EXPLICIT request is refused (LEDGER #112), and an
  # absent one (NULL) resolves to FALSE for a sparse-GRM fit and TRUE otherwise
  # -- which is exactly what the silent override produced, so no existing call
  # changes behaviour.
  if (useSparseGRMtoFitNULL && isTRUE(LOCO)) {
    stop(
      "LOCO = TRUE cannot be combined with useSparseGRMtoFitNULL = TRUE. LOCO leaves ",
      "one chromosome out of the GRM, and a sparse GRM has no per-chromosome structure ",
      "to leave out. Either omit LOCO (it then resolves to FALSE for a sparse-GRM fit), ",
      "pass LOCO = FALSE, or fit the null model with a full GRM ",
      "(useSparseGRMtoFitNULL = FALSE).\n"
    )
  }
  if (is.null(LOCO)) {
    LOCO <- !useSparseGRMtoFitNULL
  }

  modelOut <- paste0(outputPrefix, ".rda")

  # LEDGER #171 FIXED: the files created below are the writability check; if
  # this function exits by an error before it fills them, they are removed
  # (R/FitGates.R removeEmptyOutputsOnFailure). on.exit reads both variables
  # as they stand at exit, so the list grows as files are created and the
  # flag flips just before the return.
  createdOutputs <- character()
  fitSucceeded <- FALSE
  on.exit(removeEmptyOutputsOnFailure(createdOutputs, fitSucceeded), add = TRUE)

  if (skipModelFitting) {
    if (!file.exists(modelOut)) {
      stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
    }
  } else {
    if (LOCO & isLowMemLOCO) {
      modelOut <- paste(c(outputPrefix, "_noLOCO.rda"), collapse = "")
    }
    file.create(modelOut, showWarnings = TRUE)
    createdOutputs <- c(createdOutputs, modelOut)
  }

  # Genotype input: PLINK bed/bim/fam or PGEN pgen/pvar/psam, never both.
  # allowNone = TRUE because a sparse-GRM fit that skips variance-ratio
  # estimation needs no genotype file at all.
  genoIn <- checkGenoInputStep1(
    plinkFile = plinkFile, bedFile = bedFile, bimFile = bimFile, famFile = famFile,
    plink2File = plink2File, pgenFile = pgenFile, pvarFile = pvarFile, psamFile = psamFile,
    allowNone = TRUE
  )
  bimFile <- genoIn$bimFile
  bedFile <- genoIn$bedFile
  famFile <- genoIn$famFile
  pgenFile <- genoIn$pgenFile
  pvarFile <- genoIn$pvarFile
  psamFile <- genoIn$psamFile


  if (useSparseGRMtoFitNULL & genoIn$genoType == "") {
    cat("Sparse GRM is used to fit the null model and plink file is not specified, so variance ratios won't be estimated\n")
    skipVarianceRatioEstimation <- TRUE
  }

  if (!skipVarianceRatioEstimation) {
    SPAGMMATOut <- paste0(outputPrefix, "_", numMarkersForVarRatio, "markers.SAIGE.results.txt")
    # Check_OutputFile_Create(SPAGMMATOut)

    if (outputPrefix_varRatio == "") {
      outputPrefix_varRatio <- outputPrefix
    }
    varRatioFile <- paste0(outputPrefix_varRatio, ".varianceRatio.txt")

    if (!file.exists(varRatioFile)) {
      file.create(varRatioFile, showWarnings = TRUE)
      createdOutputs <- c(createdOutputs, varRatioFile) # LEDGER #171
    } else {
      if (!IsOverwriteVarianceRatioFile) {
        # LEDGER #136 (6): this is a fatal refusal and says so; it used to open
        # with "WARNING" and describe an overwrite it was not going to do
        stop(
          "ERROR: the variance ratio file ", varRatioFile,
          " already exists and --IsOverwriteVarianceRatioFile=FALSE, so the run stops ",
          "before fitting. Remove the file, use --outputPrefix_varRatio to write the ",
          "variance ratio(s) under a different prefix, or pass ",
          "--IsOverwriteVarianceRatioFile=TRUE to overwrite it\n"
        )
      } else {
        cat("The variance ratio file ", varRatioFile, " already exists. IsOverwriteVarianceRatioFile=TRUE so the file will be overwritten\n")
      }
    }
  } else {
    cat("Variance ratio estimation will be skipped\n.")
    useSparseGRMforVarRatio <- FALSE
  }


  if (useSparseGRMtoFitNULL) {
    # useSparseGRMforVarRatio = FALSE
    # LOCO was already resolved to FALSE above, where the incompatibility is
    # enforced; these announce the two options this branch takes over.
    nThreads <- 1
    # announced for every reader: the notice used to sit behind `bedFile != ""`,
    # so a --plink2File sparse fit overrode nThreads with no output at all
    cat("sparse GRM will be used to fit the NULL model and nThreads is set to 1\n")
    cat("LOCO is not applied: the null model is fitted with a sparse GRM\n")
  }


  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    if (!file.exists(sparseGRMFile)) {
      stop("sparseGRMFile ", sparseGRMFile, " does not exist!")
    }
    if (!file.exists(sparseGRMSampleIDFile)) {
      stop(
        "sparseGRMSampleIDFile ", sparseGRMSampleIDFile,
        " does not exist!"
      )
    }
  }


  # LEDGER #7 fix: always pin the TBB pool (see harness/env.sh history).
  RcppParallel:::setThreadOptions(numThreads = nThreads)
  cat(nThreads, " threads will be used ", "\n")

  # (LEDGER #90 FIXED: the sex-specific suffix block moved up, before the output
  # file names are built from outputPrefix -- it used to run after them and so
  # reached no file.)

  if (!useSparseGRMtoFitNULL | !skipVarianceRatioEstimation) {
    if (genoIn$genoType == "") {
      stop("ERROR! no genotype file was specified\n")
    }
    if (LOCO) {
      chrVec <- readVariantChromFromInput(genoIn)
      updatechrList <- updateChrStartEndIndexVec(chrVec)
      LOCO <- updatechrList$LOCO
      chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
      chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
    }
    if (!LOCO) {
      chromosomeStartIndexVec <- rep(NA, 22)
      chromosomeEndIndexVec <- rep(NA, 22)
    }

    sampleListwithGeno <- NULL
    sampleListwithGeno$IIDgeno <- readGenoSampleIDsFromInput(genoIn)
    sampleListwithGeno <- data.frame(sampleListwithGeno)
    sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
      by = 1
    )
    cat(nrow(sampleListwithGeno), " samples have genotypes\n")
  } else {
    sampleListwithGenov0 <- data.table:::fread(sparseGRMSampleIDFile,
      header = FALSE, , colClasses = c("character"), data.table = FALSE
    )
    colnames(sampleListwithGenov0) <- c("IIDgeno")
    sampleListwithGeno <- NULL
    sampleListwithGeno$IIDgeno <- sampleListwithGenov0$IIDgeno
    sampleListwithGeno <- data.frame(sampleListwithGeno)
    sampleListwithGeno$IndexGeno <- seq(1, nrow(sampleListwithGeno),
      by = 1
    )
    cat(nrow(sampleListwithGeno), " samples are in the sparse GRM\n")
  }


  pheno <- preparePhenotypeFrame(
    phenoFile, phenoCol, covarColList, qCovarCol, sampleIDColinphenoFile,
    sampleListwithGeno, traitType, invNormalize,
    sexCol = sexCol, FemaleOnly = FemaleOnly, FemaleCode = FemaleCode,
    MaleOnly = MaleOnly, MaleCode = MaleCode,
    SampleIDIncludeFile = SampleIDIncludeFile,
    sparseGRMSampleIDFile = sparseGRMSampleIDFile,
    restrictToSparseGRM = (useSparseGRMtoFitNULL & !skipVarianceRatioEstimation) | useSparseGRMforVarRatio
  )
  dataMerge_sort <- pheno$data
  formula.null <- pheno$formula.null
  covarColList <- pheno$covarColList
  hasCovariate <- pheno$hasCovariate
  indicatorGenoSamplesWithPheno <- pheno$indicatorGenoSamplesWithPheno
  # LEDGER #65 FIXED: a quantitative trait is fitted -- glm, mixed model and
  # variance ratios alike -- on the phenotype divided by the sd of its
  # covariate-adjusted values, and the saved model and variance-ratio file are
  # scaled back to phenotype units at the end (R/PhenoScale.R). tauInit is
  # given in phenotype units and converted here.
  yScale <- 1
  savedModelFiles <- character() # every .rda this run writes (LEDGER #98)
  if (traitType == "quantitative") {
    yScale <- quantPhenoScale(formula.null, dataMerge_sort, phenoCol)
    cat(
      "Quantitative trait: the null model and the variance ratios are estimated ",
      "on ", phenoCol, " divided by the sd of its covariate-adjusted values (",
      format(yScale), ") and scaled back to phenotype units before saving; ",
      "tau and coefficients printed during the fit are in the standardised ",
      "units (LEDGER #65)\n", sep = ""
    )
    dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] <-
      dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)] / yScale
    tauInit <- tauInit / yScale^2
  }
  print(head(dataMerge_sort))
  if (traitType == "binary" & (length(covarColList) > 0)) {
    out_checksep <- checkPerfectSep(formula.null,
      data = dataMerge_sort,
      minCovariateCount
    )
    covarColList <- covarColList[!(covarColList %in% out_checksep)]
    formula <- paste0(phenoCol, "~", paste0(quoteCovarNames(covarColList), collapse = "+"))
    formula.null <- as.formula(formula)


    if (length(covarColList) == 1 & traitType == "binary") {
      hasCovariate <- FALSE
    } else {
      hasCovariate <- TRUE
    }

    dataMerge_sort <- dataMerge_sort[, !(names(dataMerge_sort) %in%
      out_checksep)]
  }

  # The eligibility gates (R/FitGates.R). Placed here because dataMerge_sort is
  # final and any --minCovariateCount exclusion above has already been applied,
  # so the gates judge the model as it will actually be fitted.
  checkFitEligibility(
    data = dataMerge_sort, phenoCol = phenoCol, traitType = traitType,
    covarColList = covarColList, qCovarCol = qCovarCol,
    minCaseCount = minCaseCount, minSampleCount = minSampleCount,
    minCovariateCellCount = minCovariateCellCount
  )
  if (!hasCovariate) {
    print("No covariate is included so isCovariateOffset = FALSE")
    isCovariateOffset <- FALSE
  }

  if (isCovariateTransform & hasCovariate) {
    cat("qr transformation has been performed on covariates\n")
    out.transform <- Covariate_Transform(formula.null, data = dataMerge_sort)
    formulaNewList <- c(phenoCol, " ~ ", out.transform$Param.transform$X_name[1])
    if (length(out.transform$Param.transform$X_name) > 1) {
      for (i in c(2:length(out.transform$Param.transform$X_name))) {
        formulaNewList <- c(formulaNewList, "+", out.transform$Param.transform$X_name[i])
      }
    }
    formulaNewList <- paste0(formulaNewList, collapse = "")
    formulaNewList <- paste0(formulaNewList, "-1")
    formula.new <- as.formula(paste0(formulaNewList, collapse = ""))
    data.new <- data.frame(cbind(out.transform$Y, out.transform$X1))

    colnames(data.new) <- c(phenoCol, out.transform$Param.transform$X_name)
    cat("colnames(data.new) is ", colnames(data.new), "\n")
    cat(
      "out.transform$Param.transform$qrr: ", dim(out.transform$Param.transform$qrr),
      "\n"
    )
  } else {
    formula.new <- formula.null
    data.new <- dataMerge_sort
    out.transform <- NULL
  }

  # if(!hasCovariate){
  # 	covoffset = rep(0,nrow(data.new))
  #    }

  if (isCovariateOffset) {
    # E6 (2026-09-14): this design matrix and GLM used to be built for EVERY
    # run, ~15 lines above, and read only here -- so both branches below fitted
    # the identical GLM twice, and mmat (N x p+1) stayed live through the GRM
    # fit on the branch that never read it. One consequence to state rather
    # than discover: a malformed binary phenotype now meets the repo's own 0/1
    # check below before any glm() call, where before glm()'s message came first.
    mmat <- model.matrix(formula.new, data = data.new, na.action = NULL)
    if (traitType == "binary") {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = binomial
      )
    } else {
      modwitcov <- glm(formula.new,
        data = data.new,
        family = gaussian(link = "identity")
      )
    }
    if (traitType == "binary") {
      # These coefficients are about to be frozen into covoffset and never
      # re-estimated, so separation here is not recoverable downstream.
      checkSeparation(modwitcov, phenoCol, separationCheck,
                      separationMaxFraction)
    }
    covoffset <- mmat[, -1, drop = FALSE] %*% modwitcov$coefficients[-1]
    print("isCovariateOffset=TRUE, so fixed effects coefficients won't be estimated.")
    formula.new.withCov <- formula.new
    formula_nocov <- paste0(phenoCol, "~ 1")
    formula.new <- as.formula(formula_nocov)
    hasCovariate <- FALSE
  } else {
    covoffset <- rep(0, nrow(data.new))
  }
  data.new$covoffset <- covoffset

  if (useSparseGRMtoFitNULL | useSparseGRMforVarRatio) {
    sparseGRMtest <- getsubGRM(
      sparseGRMFile, sparseGRMSampleIDFile,
      relatednessCutoff, dataMerge_sort$IID
    )
    m4 <- gen_sp_v2(sparseGRMtest)
    cat("Setting up sparse GRM using ", sparseGRMFile, " and ", sparseGRMSampleIDFile, "\n")
    cat("Dimension of the sparse GRM is ", dim(m4), "\n")
    A <- summary(m4)
    locationMatinR <- rbind(A$i - 1, A$j - 1)
    valueVecinR <- A$x
    setupSparseGRM(dim(m4)[1], locationMatinR, valueVecinR)
    rm(sparseGRMtest)
  }

  if (!skipVarianceRatioEstimation) {
    isVarianceRatioinGeno <- TRUE
    if (isCateVarianceRatio) {
      minMAC_varRatio <- min(cateVarRatioMinMACVecExclude)
      maxMAC_varRatio <- max(cateVarRatioMaxMACVecInclude)
      # maxMAC_varRatio = max(max(cateVarRatioMaxMACVecInclude), minMACforGRM - 1)
      cat("Categorical variance ratios will be estimated. Please make sure there are at least 200 markers in each MAC category.\n")
    } else {
      minMAC_varRatio <- 20
      maxMAC_varRatio <- -1 # will randomly select markers from the plink file and leave them out when constructing GRM
    }
    setminMAC_VarianceRatio(minMAC_varRatio, maxMAC_varRatio, isVarianceRatioinGeno)
  }

  # set up parameters
  if (minMAFforGRM > 0) {
    cat(
      "Markers in the Plink file with MAF < ", minMAFforGRM,
      " will be removed before constructing GRM\n"
    )
  }
  if (maxMissingRateforGRM > 0) {
    cat("Markers in the Plink file with missing rate > ", maxMissingRateforGRM, " will be removed before constructing GRM\n")
  }

  setminMAFforGRM(minMAFforGRM)
  setmaxMissingRateforGRM(maxMissingRateforGRM)

  # The gates that decide WHICH markers the ingest keeps and which it stages for
  # the variance ratio. They are part of the ingest's identity: two phenotypes
  # may share a loaded store only if these agree as well as the sample set.
  ingestGates <- list(
    minMAFforGRM = minMAFforGRM,
    maxMissingRateforGRM = maxMissingRateforGRM,
    vr = if (skipVarianceRatioEstimation) {
      NULL
    } else {
      c(minMAC_varRatio, maxMAC_varRatio)
    }
  )


  # ---- the fit, the save and the variance ratios: ONE body for both traits.
  # The two traits used to be two 270-line copies of this block differing in
  # the glm family, the fitter, and the four lines marked "by trait" below
  # (the LOCO refit's mu2, offset and tol.coef, and the quantitative
  # rescaling); every other line was the same text (finding 2 of the
  # 2026-09-05 code-quality review).
  if (!traitType %in% c("binary", "quantitative")) {
    stop("traitType must be binary or quantitative")
  }
  isBinary <- traitType == "binary"
  family <- if (isBinary) binomial else gaussian(link = "identity")
  fitter <- if (isBinary) glmmkin.ai_PCG_Rcpp_Binary else glmmkin.ai_PCG_Rcpp_Quantitative
  if (isBinary) {
    cat(phenoCol, " is a binary trait\n")
    uniqPheno <- sort(unique(dataMerge_sort[, which(colnames(dataMerge_sort) == phenoCol)]))
    if (uniqPheno[1] != 0 | uniqPheno[2] != 1) {
      stop("ERROR! phenotype value needs to be 0 or 1 \n")
    }
    print("formula.new")
    print(formula.new)
    print("head(data.new)")
    print(head(data.new))
  } else {
    cat(phenoCol, " is a quantitative trait\n")
  }
  {
    if (!isCovariateOffset) {
      fit0 <- glm(formula.new, data = data.new, family = family)
      if (isBinary) {
        checkSeparation(fit0, phenoCol, separationCheck,
                        separationMaxFraction)
      }
      Xorig <- NULL
    } else {
      fit0orig <- glm(formula.new.withCov, data = data.new, family = family)
      Xorig <- model.matrix(fit0orig)
      rm(fit0orig)
      gc()
      fit0 <- glm(formula.new,
        data = data.new, offset = covoffset,
        family = family
      )
    }


    cat("glm:\n")
    print(fit0)
    if (isBinary) print(summary(fit0))

    obj.noK <- NULL
    if (!skipModelFitting) {
      setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
      if (useSparseGRMtoFitNULL) {
        cat("usePCGwithSparseGRM: ", usePCGwithSparseGRM, "\n")
        setisUsePCGwithSparseSigma(usePCGwithSparseGRM)
        resetSigmaFactoriser() # LEDGER #79: per-fit factorisation counters
      }
      cat("Start fitting the NULL GLMM\n")
      t_begin <- proc.time()
      print(t_begin)
      re_fit <- system.time(modglmm <- fitter(bedFile, bimFile, famFile, Xorig, isCovariateOffset,
        fit0,
        tau = c(0, 0), fixtau = c(0, 0), maxiter = maxiter,
        tol = tol, tolCoef = tolCoef, verbose = TRUE, nrun = nrun, tolPCG = tolPCG,
        maxiterPCG = maxiterPCG, subPheno = dataMerge_sort, indicatorGenoSamplesWithPheno = indicatorGenoSamplesWithPheno,
        obj.noK = obj.noK, out.transform = out.transform,
        tauInit = tauInit, memoryChunk = memoryChunk,
        LOCO = LOCO, chromosomeStartIndexVec = chromosomeStartIndexVec,
        chromosomeEndIndexVec = chromosomeEndIndexVec,
        traceCVcutoff = traceCVcutoff, isCovariateTransform = isCovariateTransform, genoIn = genoIn,
        ingestGates = ingestGates,
        isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        isLowMemLOCO = isLowMemLOCO
      ))
      # Refuse to SAVE a model whose variance components never settled: the
      # .rda would otherwise look like any other and be tested like any other.
      # Routed through the switch like every other gate -- SAIGE_FIT_GATES=0
      # must mean ALL of them, or the suite's degenerate fixtures would still
      # be refused here.
      if (fitGatesEnabled()) {
        checkConvergence(modglmm$converged, phenoCol, convergenceCheck)
      }
      if (useSparseGRMtoFitNULL) {
        # LEDGER #79 FIXED: the sparse Sigma factorisation is reused across
        # right-hand sides, so the fit does far fewer factorisations than solves.
        cnt <- getSigmaSolveCounts()
        cat("Sparse Sigma direct solves: ", cnt[2], " solves with ", cnt[1],
            " factorisations\n", sep = "")
      }
      modglmm$obj.glm.null$model <- data.frame(modglmm$obj.glm.null$model)

      for (x in names(modglmm$obj.glm.null)) {
        attr(modglmm$obj.glm.null[[x]], ".Environment") <- c()
      }
      # modglmm$offset = covoffset
      if (isCovariateOffset) {
        modglmm$offset <- covoffset
      } else {
        if (hasCovariate) {
          # modglmm$offset = covoffset
          # data.new = data.new[,-c(1,2), drop=F]
          # data.new = as.matrix(data.new[,-ncol(data.new), drop=F])
          data.new.X <- model.matrix(fit0)[, -1, drop = FALSE]
          modglmm$offset <- data.new.X %*% (as.vector(modglmm$coefficients[-1]))
          if (LOCO & !isLowMemLOCO) {
            for (j in 1:22) {
              if (modglmm$LOCOResult[[j]]$isLOCO) {
                modglmm$LOCOResult[[j]]$offset <- data.new.X %*% (as.vector(modglmm$LOCOResult[[j]]$alpha0[-1]))
              }
            }
          }
        } else {
          modglmm$offset <- covoffset
        }
      }

      if (skipVarianceRatioEstimation & useSparseGRMtoFitNULL) {
        family <- fit0$family
        eta <- modglmm$linear.predictors
        mu <- modglmm$fitted.values
        mu.eta <- family$mu.eta(eta)
        sqrtW <- mu.eta / sqrt(family$variance(mu))
        W <- sqrtW^2
        tauVecNew <- modglmm$theta
        Sigma_iX <- getSigma_X(W, tauVecNew, modglmm$X, maxiterPCG, tolPCG)
        Sigma_iXXSigma_iX <- Sigma_iX %*% (solve(t(modglmm$X) %*% Sigma_iX))
        modglmm$Sigma_iXXSigma_iX <- Sigma_iXXSigma_iX
      }

      modglmm$useSparseGRMtoFitNULL <- useSparseGRMtoFitNULL
      # save(modglmm, file = modelOut)
      tau <- modglmm$theta

      alpha0 <- modglmm$coefficients

      if (!is.null(out.transform) & is.null(fit0$offset)) {
        coef.alpha <- Covariate_Transform_Back(alpha0, out.transform$Param.transform)
        modglmm$coefficients <- coef.alpha
      }


      if (LOCO & isLowMemLOCO) {
        modglmm$LOCOResult <- NULL
        modglmm$LOCO <- FALSE
        chromosomeStartIndexVec <- modglmm$chromosomeStartIndexVec
        chromosomeEndIndexVec <- modglmm$chromosomeEndIndexVec
        modglmm$chromosomeStartIndexVec <- NULL
        modglmm$chromosomeEndIndexVec <- NULL
        modelOut <- paste(c(outputPrefix, "_noLOCO.rda"), collapse = "")
        save(modglmm, file = modelOut)
        savedModelFiles <- c(savedModelFiles, modelOut)
        modglmm$LOCO <- TRUE
        modglmm$Y <- NULL
        eta0 <- modglmm$linear.predictors
        modglmm$linear.predictors <- NULL
        modglmm$coefficients <- NULL
        modglmm$cov <- NULL
        modglmm$fitted.values <- NULL
        modglmm$residuals <- NULL
        modglmm$obj.noK <- NULL
        y <- fit0$y
        # by trait: the binary refit used the glm's offset (none unless
        # isCovariateOffset), the quantitative refit the model's covariate
        # offset built above. Preserved as found; see next.md.
        offset0 <- if (isBinary) fit0$offset else modglmm$offset
        if (is.null(offset0)) {
          offset0 <- rep(0, length(y))
        }


        modglmm$offset <- NULL
        gc()
        # save(modglmm, file = modelOut)
        set_Diagof_StdGeno_LOCO()
        modglmm$LOCOResult <- list()


        for (j in 1:22) {
          startIndex <- chromosomeStartIndexVec[j]
          endIndex <- chromosomeEndIndexVec[j]
          if (!is.na(startIndex) && !is.na(endIndex)) {
            cat("leave chromosome ", j, " out\n")
            setStartEndIndex(startIndex, endIndex, j - 1)
                        # by trait: the quantitative refit never passed tol.coef (default 0.1)
            re.coef_LOCO <- Get_Coef_LOCO(y, X = model.matrix(fit0), tau, family = fit0$family, alpha = alpha0, eta = eta0, offset = offset0, verbose = TRUE, maxiterPCG = maxiterPCG, tolPCG = tolPCG, maxiter = maxiter, tol.coef = if (isBinary) tolCoef else 0.1)
            cov <- re.coef_LOCO$cov
            alpha <- re.coef_LOCO$alpha
            eta <- re.coef_LOCO$eta
            Y <- re.coef_LOCO$Y
            mu <- re.coef_LOCO$mu
            res <- y - mu
            mu2 <- if (isBinary) mu * (1 - mu) else rep((1 / (tau[1])), length(res)) # by trait

            if (!is.null(out.transform) & is.null(fit0$offset)) {
              coef.alpha <- Covariate_Transform_Back(alpha, out.transform$Param.transform)
            } else {
              coef.alpha <- alpha
            }


            if (!isCovariateOffset) {
              obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, X = model.matrix(fit0))
            } else {
              obj.noK <- ScoreTest_NULL_Model(mu, mu2, y, Xorig)
            }

            modglmm$LOCOResult[[j]] <- list(isLOCO = TRUE, coefficients = coef.alpha, linear.predictors = eta, fitted.values = mu, Y = Y, residuals = res, cov = cov, obj.noK = obj.noK)

            if (!isCovariateOffset & hasCovariate) {
              data.new.X <- model.matrix(fit0)[, -1, drop = FALSE]
              modglmm$LOCOResult[[j]]$offset <- data.new.X %*% (as.vector(alpha[-1]))
            }
            modelOutbychr <- paste(c(outputPrefix, "_chr", j, ".rda"), collapse = "")
            if (j != 22) {
              for (j1 in (j + 1):22) {
                modglmm$LOCOResult[[j1]] <- list(NULL)
              }
            }
            save(modglmm, file = modelOutbychr)
            savedModelFiles <- c(savedModelFiles, modelOutbychr)
            modglmm$LOCOResult[[j]] <- list(NULL)
            gc()
          } else {
            modglmm$LOCOResult[[j]] <- list(isLOCO = FALSE)
          }
        }
        gc()
        # modelOut_nonauto = paste(c(outputPrefix,"_noLOCO.rda"), collapse="")
      } else {
        save(modglmm, file = modelOut)
        savedModelFiles <- c(savedModelFiles, modelOut)
      }

      gc()
      t_end <- proc.time()
      print(t_end)
      cat("t_end - t_begin, fitting the NULL model took\n")
      print(t_end - t_begin)
      cat("of which the null-model fit itself (glmmkin.ai_PCG_Rcpp) took\n")
      print(re_fit)
    } else {
      cat("Skip fitting the NULL GLMM\n")


      if (!file.exists(modelOut)) {
        stop("skipModelFitting=TRUE but ", modelOut, " does not exist\n")
      }

      load(modelOut)
      if (is.null(modglmm$LOCO)) {
        modglmm$LOCO <- FALSE
      }
      ensureGenoLoaded(
        genoIn, dataMerge_sort$IndexGeno, indicatorGenoSamplesWithPheno, memoryChunk,
        isDiagofKinSetAsOne, ingestGates
      )
      setisUseSparseSigmaforNullModelFitting(useSparseGRMtoFitNULL)
      if (useSparseGRMtoFitNULL) {
        setisUsePCGwithSparseSigma(usePCGwithSparseGRM)
      }
    }

    if (!skipVarianceRatioEstimation) {
      if (LOCO) {
        MsubIndVec <- getQCdMarkerIndex()
        # cat("MsubIndVec", MsubIndVec, "\n")
        print(length(MsubIndVec))
        chrVec <- readVariantChromFromInput(genoIn)
        print(length(chrVec))
        chrVec <- chrVec[which(MsubIndVec == TRUE)]
        updatechrList <- updateChrStartEndIndexVec(chrVec)
        LOCO <- updatechrList$LOCO
        chromosomeStartIndexVec <- updatechrList$chromosomeStartIndexVec
        chromosomeEndIndexVec <- updatechrList$chromosomeEndIndexVec
        set_Diagof_StdGeno_LOCO()
      }
      cat("Start estimating variance ratios\n")
      load(modelOut)
      if (skipModelFitting && !isBinary) {
        # the file is in phenotype units (LEDGER #65); the estimation runs standardised
        modglmm <- scaleQuantModel(modglmm, 1 / yScale)
      }
      extractVarianceRatio(
        obj.glmm.null = modglmm,
        seed = seed,
        obj.glm.null = fit0, maxiterPCG = maxiterPCG, varRatioScale = 1 / yScale^2,
        tolPCG = tolPCG, numMarkers = numMarkersForVarRatio, varRatioOutFile = varRatioFile,
        ratioCVcutoff = ratioCVcutoff, testOut = SPAGMMATOut,
        bedFile = bedFile, bimFile = bimFile, famFile = famFile, genoIn = genoIn, chromosomeStartIndexVec = chromosomeStartIndexVec,
        chromosomeEndIndexVec = chromosomeEndIndexVec,
        isCateVarianceRatio = isCateVarianceRatio, cateVarRatioIndexVec = cateVarRatioIndexVec,
        useSparseGRMforVarRatio = useSparseGRMforVarRatio, sparseGRMFile = sparseGRMFile,
        sparseGRMSampleIDFile = sparseGRMSampleIDFile,
        numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
        relatednessCutoff = relatednessCutoff, useSparseGRMtoFitNULL = useSparseGRMtoFitNULL,
        nThreads = nThreads, cateVarRatioMinMACVecExclude = cateVarRatioMinMACVecExclude,
        cateVarRatioMaxMACVecInclude = cateVarRatioMaxMACVecInclude,
        minMAFforGRM = minMAFforGRM, isDiagofKinSetAsOne = isDiagofKinSetAsOne,
        includeNonautoMarkersforVarRatio = includeNonautoMarkersforVarRatio
      )

      # print(gc(v=T))
      gc()
    } else {
      cat("Skip estimating variance ratios\n")
    }
    if (releaseGenoAtEnd) {
      releaseGeno()
    }
    if (!skipModelFitting && !isBinary) {
      # LEDGER #65 FIXED: everything above ran on the standardised phenotype;
      # the saved model(s) go back to phenotype units. The variance-ratio file
      # was already written in phenotype units (varRatioScale above). Exactly
      # the files THIS run wrote (LEDGER #98: a glob over _chr1..22 used to
      # rescale a previous run's leftovers under the same prefix).
      rescaleSavedQuantModels(savedModelFiles, yScale)
      cat("Saved model(s) scaled back to phenotype units (factor ", format(yScale), "): ",
          paste(basename(savedModelFiles), collapse = ", "), "\n", sep = "")
    }
    # The paths fitNULLGLMM ACTUALLY wrote, which the multi-phenotype driver
    # needs for the manifest: modelOut is reassigned under --isLowMemLOCO, and
    # outputPrefix is mutated in place by --FemaleOnly/--MaleOnly.
    fitSucceeded <- TRUE # LEDGER #171: everything above wrote what it meant to
    return(invisible(list(
      phenoCol = phenoCol, modelOut = modelOut,
      varRatioFile = if (skipVarianceRatioEstimation) NULL else varRatioFile,
      nSamples = nrow(dataMerge_sort)
    )))
  }
}


## suggested by Shawn 01-19-2018
# LEDGER #85 FIXED: the `traitType` formal was declared and never read (survival,
# which forced it, left in Phase 3a), and extractNglmm called without it -- legal
# only because the promise was never forced. Dropped from the signature.
Covariate_Transform <- function(formula, data) {
  X1 <- model.matrix(formula, data = data)
  #  X1=X1[,c(2:ncol(X1))] #remove intercept
  formula.frame <- model.frame(formula, data = data)
  Y <- model.response(formula.frame, type = "any")
  X_name <- colnames(X1)

  # First run linear regression to identify multi collinearity
  out.lm <- lm(Y ~ X1 - 1, data = data)
  #  out.lm<-lm(Y ~ X1, data=data)
  idx.na <- which(is.na(out.lm$coef))
  if (length(idx.na) > 0) {
    print(head(X1))
    X1 <- X1[, -idx.na, drop = FALSE] # LEDGER #36 FIXED: one surviving column stayed a matrix
    print(head(X1))
    cat("Warning: multi collinearity is detected in covariates! ", X_name[idx.na], " will be excluded in the model\n")
    X_name <- X_name[-idx.na]
    # cat("Warning: multi collinearity is detected in covariates! ", X_name[idx.na], " will be excluded in the model\n")
  }
  if (!(1 %in% idx.na)) {
    X_name[1] <- "minus1"
  }

  # QR decomposition
  Xqr <- qr(X1)
  X1_Q <- qr.Q(Xqr)
  qrr <- qr.R(Xqr)

  N <- nrow(X1)

  # Make square summation=N (so mean=1)
  X1_new <- X1_Q * sqrt(N)
  Param.transform <- list(qrr = qrr, N = N, X_name = X_name, idx.na = idx.na)
  re <- list(Y = Y, X1 = X1_new, Param.transform = Param.transform)
}

# In case to recover original scale coefficients
# X \beta = Q R \beta = (Q \sqrt(N)) ( R \beta / \sqrt(N))
# So coefficient from fit.new is the same as R \beta / \sqrt(N)
Covariate_Transform_Back <- function(coef, Param.transform) {
  # coef<-fit.new$coef; Param.transform=out.transform$Param.transform
  coef1 <- coef * sqrt(Param.transform$N)
  coef.org <- solve(Param.transform$qrr, coef1)

  names(coef.org) <- Param.transform$X_name
  coef.org
}


pcg <- function(A, b, M = NULL, maxiter = 1e+05, tol = 1e-06) {
  # A<-a; b<-c1[,1]; M<-NULL;maxiter = 1e+05; tol = 1e-06
  if (is.null(M)) {
    dA <- diag(A)
    dA[which(dA == 0)] <- 1e-04
    # print("dA")
    # print(dA)
    Minv <- 1 / dA
  } else {
    Minv <- solve(M)
  }
  print("R")
  print(Minv)
  x <- rep(0, length(b))
  r <- b
  if (is.null(M)) {
    z <- Minv * r
  } else {
    z <- Minv %*% r
  }
  p <- z
  iter <- 0
  sumr2 <- sum(r^2)
  while (sumr2 > tol & iter < maxiter) {
    iter <- iter + 1
    #    cat("iter is ", iter, "\n")
    Ap <- crossprod(p, A)[1, ]
    a <- as.numeric((t(r) %*% z) / (t(p) %*% Ap))
    x <- x + a * p
    r1 <- r - a * Ap

    if (is.null(M)) {
      z1 <- Minv * r1
    } else {
      z1 <- Minv %*% r1
    }


    bet <- as.numeric((t(z1) %*% r1) / (t(z) %*% r))

    p <- z1 + bet * p

    z <- z1
    r <- r1
    sumr2 <- sum(r^2)
  }
  if (iter >= maxiter) {
    x <- "pcg did not converge. You may increase maxiter number."
  }
  x
}


pcgSparse <- function(A, b, M = NULL, maxiter = 1e+05, tol = 1e-06) {
  # A<-a; b<-c1[,1]; M<-NULL;maxiter = 1e+05; tol = 1e-06
  if (is.null(M)) {
    dA <- diag(A)
    dA[which(dA == 0)] <- 1e-04
    # print("dA")
    # print(dA)
    Minv <- 1 / dA
  } else {
    Minv <- solve(M)
  }
  x <- rep(0, length(b))
  r <- b
  if (is.null(M)) {
    z <- Minv * r
  } else {
    z <- Minv %*% r
  }
  p <- z
  iter <- 0
  sumr2 <- sum(r^2)
  print("psparse0")
  psparse <- Matrix:::sparseMatrix(i = rep(1, length(p)), j = c(1:length(p)), x = as.vector(p))
  cat("nrow(psparse) ", nrow(psparse), "\n")
  cat("ncol(psparse) ", ncol(psparse), "\n")
  print(class(psparse))

  while (sumr2 > tol & iter < maxiter) {
    iter <- iter + 1
    # Ap = crossprod(p, A)[1,]
    print(class(psparse)[1])
    print(class(A)[1])
    print(nrow(psparse))
    print(ncol(psparse))
    print(nrow(A))
    print(ncol(A))

    # Ap = sparse_row_idx_mult(psparse, A)
    Ap <- psparse %*% A[1, ]
    print("psparse1")
    a <- as.numeric((t(r) %*% z) / (t(p) %*% Ap))
    x <- x + a * p
    r1 <- r - a * Ap

    if (is.null(M)) {
      z1 <- Minv * r1
    } else {
      z1 <- Minv %*% r1
    }


    bet <- as.numeric((t(z1) %*% r1) / (t(z) %*% r))

    p <- z1 + bet * p

    z <- z1
    r <- r1
    sumr2 <- sum(r^2)
  }
  if (iter >= maxiter) {
    x <- "pcg did not converge. You may increase maxiter number."
  }
  x
}


# https://gist.github.com/bobthecat/5024079
bigGRMPar <- function(x, nblocks = 10, verbose = TRUE, ncore = 1, relatednessCutoff = 0) {
  #  library(foreach)
  #  library(doParallel)
  # register cores
  doParallel:::registerDoParallel(ncore)

  NCOL <- ncol(x)
  NROW <- nrow(x)

  ## test if ncol(x) %% nblocks gives remainder 0
  #  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}

  #  if(NCOL %% nblocks == 0){
  ## split column numbers into 'nblocks' groups
  #  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

  SPLIT <- split(1:NCOL, ceiling(seq_along(1:NCOL) / (NCOL %/% nblocks)))
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  #  }else{
  #    NCOLNEW = NCOL - (NCOL%%nblocks + NCOL%/%nblocks)
  #    SPLIT <- split(1:NCOLNEW, rep(1:(nblocks-1), each = NCOLNEW/(nblocks-1)))
  #  }
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  `%dopar%` <- foreach::`%dopar%`
  results <- foreach:::foreach(i = 1:nrow(COMBS), .combine = "rbind") %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    GRM <- t(x[, G1]) %*% (x[, G2])

    if (sum(G1 != G2) == 0) {
      GRM[lower.tri(GRM, diag = FALSE)] <- 0
    }
    GRM <- GRM / NROW
    indice <- which(GRM >= relatednessCutoff, arr.ind = TRUE)
    cbind(G1[indice[, 1]], G2[indice[, 2]])
    #       resultsIVec =c(resultsIVec, G1[indice[,1]])
    #        resultsJVec =c(resultsJVec, G2[indice[,2]])
    # corMAT[G1, G2] <- COR
    # corMAT[G2, G1] <- t(COR)
    #       GRM <- NULL
    #       indice = NULL
  }
  # gc()
  results
}


bigGRMPar_new <- function(nblocks = 10, verbose = TRUE, ncore = 1, relatednessCutoff = 0) {
  doParallel:::registerDoParallel(ncore)

  NCOL <- getNColStdGenoMultiMarkersMat()
  NROW <- getNRowStdGenoMultiMarkersMat()
  cat("NCOL: ", NCOL, "\n")
  cat("NROW: ", NROW, "\n")

  nblocks <- NCOL %/% 10

  ## test if ncol(x) %% nblocks gives remainder 0
  #  if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}

  #  if(NCOL %% nblocks == 0){
  ## split column numbers into 'nblocks' groups
  #  SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))

  SPLIT <- split(1:NCOL, ceiling(seq_along(1:NCOL) / (NCOL %/% nblocks)))
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  #  }else{
  #    NCOLNEW = NCOL - (NCOL%%nblocks + NCOL%/%nblocks)
  #    SPLIT <- split(1:NCOLNEW, rep(1:(nblocks-1), each = NCOLNEW/(nblocks-1)))
  #  }
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  `%dopar%` <- foreach::`%dopar%`
  results <- foreach:::foreach(i = 1:nrow(COMBS), .combine = "rbind") %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]] - 1
    G2 <- SPLIT[[COMB[2]]] - 1
    # if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    G1M <- getColfromStdGenoMultiMarkersMat(G1)
    G2M <- getColfromStdGenoMultiMarkersMat(G2)
    # 	cat("G1M ", dim(G1M), "\n")
    # 	cat("G2M ", dim(G2M), "\n")
    GRM <- t(G1M) %*% (G2M)
    # 	cat("G1M: ", G1M, "\n")
    # 	cat("G2M: ", G2M, "\n")
    # 	cat("GRM: ", GRM, "\n")
    if (sum(G1 != G2) == 0) {
      GRM[lower.tri(GRM, diag = FALSE)] <- 0
    }
    GRM <- GRM / NROW
    indice <- which(GRM >= relatednessCutoff, arr.ind = TRUE)
    # 	cat("GRM[1:10,1:10]: ", GRM[1:10,1:10], "\n")
    # 	cat("relatednessCutoff: ", relatednessCutoff, "\n")
    # 	cat("indice: ", indice, "\n")
    cbind(G1[indice[, 1]], G2[indice[, 2]])
    #       resultsIVec =c(resultsIVec, G1[indice[,1]])
    #        resultsJVec =c(resultsJVec, G2[indice[,2]])
    # corMAT[G1, G2] <- COR
    # corMAT[G2, G1] <- t(COR)
    #       GRM <- NULL
    #       indice = NULL
  }
  # gc()
  results
}


# refineKinPar = function(iMat, relatednessCutoff, W, tauVecNew, nblocks = 10, verbose = TRUE, ncore= 1){
refineKinPar <- function(relatednessCutoff, W, tauVecNew, nblocks = 10, verbose = TRUE, ncore = 1) {
  # chunk iMat
  doParallel:::registerDoParallel(ncore)
  NROW <- nrow(iMat)
  SPLIT <- split(1:NROW, ceiling(seq_along(1:NROW) / (NROW %/% nblocks)))
  GRMvec <- rep(0, NROW)
  print(length(GRMvec))
  `%dopar%` <- foreach::`%dopar%`
  mMarkers <- gettotalMarker()
  for (j in 1:mMarkers) {
    cat("j is ", j, "\n")
    stdGeno <- Get_OneSNP_StdGeno(j - 1)
    results <- foreach:::foreach(i = 1:length(SPLIT)) %dopar% {
      G1 <- SPLIT[[i]]
      # print(G1)
      if (verbose) cat("Block", i, "\n")
      flush.console()
      for (m in G1) {
        GRMvec[m] <- GRMvec[m] + stdGeno[iMat[m, 1]] * stdGeno[iMat[m, 2]] / mMarkers
      }
    }
  }
  GRMvec
}


# createSparseKinParallel = function(markerIndexVec, nblocks, ncore, relatednessCutoff, W, tauVecNew){
createSparseKinParallel <- function(nblocks, ncore, relatednessCutoff) {
  # get MAT
  # MAT = Get_MultiMarkersBySample_StdGeno_Mat(markerIndexVec)
  # MAT = Get_MultiMarkersBySample_StdGeno_Mat()
  setRelatednessCutoff(relatednessCutoff)
  Get_MultiMarkersBySample_StdGeno_Mat()
  #  cat("dim(MAT) is ", dim(MAT), "\n")
  tp0 <- proc.time()
  #  indexVec = bigGRMPar(MAT, nblocks = nblocks, verbose = FALSE, ncore = nblocks, relatednessCutoff = relatednessCutoff)
  #  indexVec = bigGRMPar_new(nblocks = nblocks, verbose = TRUE, ncore = nblocks, relatednessCutoff = relatednessCutoff)

  # printComb(3) -- the R call was deleted 2026-09-14 (C6): it printed two junk
  # lines on every step-0 / internal-GRM run. The C++ printer stays (LEDGER #63).
  # indexVec = findIndiceRelatedSample()
  findIndiceRelatedSample()

  # print(indexVec)
  tp1 <- proc.time()
  cat("tp1 - tp0: ", tp1 - tp0, "\n")
  #  cat(indexVec)
  # sparseKinList = refineKin(indexVec-1, relatednessCutoff, W, tauVecNew)
  sparseKinList <- refineKin(relatednessCutoff)

  #  sparseKinList$kinValue = sparseKinList$kinValue * tauVecNew[2]

  Nval <- getNnomissingOut()

  sparseKinList$iIndex <- c(sparseKinList$iIndex, seq(1:Nval))
  sparseKinList$jIndex <- c(sparseKinList$jIndex, seq(1:Nval))
  #  diagKin = getDiagOfSigma(W, tauVecNew)
  #  diagKin = rep(1, Nval)
  diagKin <- get_DiagofKin()
  sparseKinList$kinValue <- c(sparseKinList$kinValue, diagKin)
  # sparseKinList = refineKin(indexVec, relatednessCutoff, W, tauVecNew)
  # GRMvec = refineKinPar(indexVec, relatednessCutoff = relatednessCutoff, W = W, tauVecNew = tauVecNew, nblocks = nblocks, verbose = TRUE, ncore= nblocks)
  # sparseKinList = shortenList(indexVec-1, GRMvec, relatednessCutoff, W, tauVecNew)
  tp2 <- proc.time()
  cat("tp2 - tp1: ", tp2 - tp1, "\n")
  sparseKinList
}


getSparseSigma <- function(bedFile,
                           bimFile,
                           famFile,
                           outputPrefix = "",
                           sparseGRMFile = NULL,
                           sparseGRMSampleIDFile = "",
                           numRandomMarkerforSparseKin = 1000,
                           relatednessCutoff = 0.125,
                           minMAFforGRM = 0,
                           nThreads = 1,
                           isDiagofKinSetAsOne = FALSE,
                           obj.glmm.null,
                           W,
                           tauVecNew,
                           genoIn = NULL) {
  if (is.null(genoIn)) genoIn <- genoInputPlink(bedFile, bimFile, famFile)
  cat("sparse GRM will be used\n")
  # cat("sparseGRMFile is ", sparseGRMFile, "\n")
  # sparseGRMFile = paste0(outputPrefix, ".sparseGRM.mtx")
  if (is.null(sparseGRMFile)) {
    cat("sparseGRMFile is not specified and the sparse GRM will be constructed\n")
    outputPrefix1 <- paste0(outputPrefix, "_allPlinksamples")
    sparseGRMList <- createSparseGRM(
      bedFile = genoIn$bedFile, bimFile = genoIn$bimFile, famFile = genoIn$famFile,
      pgenFile = genoIn$pgenFile, pvarFile = genoIn$pvarFile, psamFile = genoIn$psamFile,
      outputPrefix = outputPrefix1,
      numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
      relatednessCutoff = relatednessCutoff,
      isDiagofKinSetAsOne = isDiagofKinSetAsOne,
      nThreads = nThreads,
      minMAFforGRM = minMAFforGRM,
      isSetGeno = FALSE,
      isWritetoFiles = FALSE
    )
    sparseGRM <- sparseGRMList$sparseGRM
    sparseGRMSampleID <- data.frame(sampleID = sparseGRMList$sparseGRMSampleID)
    colnames(sparseGRMSampleID) <- c("sampleID")
    rm(sparseGRMList)
    # sparseGRMFile = paste0(outputPrefix1,"_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx")
    # cat("sparse GRM for all samples in the plink files is stored in ", sparseGRMFile, "\n")
    # sparseGRMSampleIDFile = paste0(outputPrefix1,"_relatednessCutoff_",relatednessCutoff,"_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
    # cat("sample IDs for the constructed sparse GRM are stored in ", sparseGRMSampleIDFile, "\n")
  } else { # if(sparseGRMFile=="")
    cat("sparse GRM has been specified\n")
    cat("read in sparse GRM from ", sparseGRMFile, "\n")
    sparseGRMLarge <- Matrix:::readMM(sparseGRMFile)
    # cat("sparseSigmaFile: ", sparseSigmaFile, "\n")
    if (sparseGRMSampleIDFile != "") {
      if (!file.exists(sparseGRMSampleIDFile)) {
        stop("ERROR! sparseSigmaSampleIDFile ", sparseGRMSampleIDFile, " does not exsit\n")
      } else {
        sparseGRMSampleID <- data.frame(data.table:::fread(sparseGRMSampleIDFile, header = FALSE, stringsAsFactors = FALSE, colClasses = c("character")))
        colnames(sparseGRMSampleID) <- c("sampleID")
      }
    } else { # end of if(sparseSigmaSampleIDFile != "")
      stop("ERROR! sparseSigmaSampleIDFile is not specified\n")
    }
    # }

    sparseGRMSampleID$IndexGRM <- c(1:nrow(sparseGRMSampleID))
    cat("length(sparseGRMSampleID$IndexGRM): ", length(sparseGRMSampleID$IndexGRM), "\n")
    cat("nrow(sparseGRMSampleID): ", nrow(sparseGRMSampleID), "\n")
    sampleInModel <- NULL
    sampleInModel$IID <- obj.glmm.null$sampleID
    sampleInModel <- data.frame(sampleInModel)
    sampleInModel$IndexInModel <- seq(1, length(sampleInModel$IID), by = 1)
    cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
    mergeID <- merge(sampleInModel, sparseGRMSampleID, by.x = "IID", by.y = "sampleID")
    mergeID <- mergeID[with(mergeID, order(IndexInModel)), ]
    print(dim(mergeID))
    print(head(mergeID))
    indexIDofGRM <- mergeID$IndexGRM
    # cat("indexIDofGRM = ", indexIDofGRM, "\n")
    # cat("Subset sparse GRM to be ", indexIDofSigma," by ", indexIDofSigma, "\n")
    sparseGRM <- sparseGRMLarge[indexIDofGRM, indexIDofGRM]
    rm(sparseGRMLarge)
  }
  # sparseGRMFile = paste0(outputPrefix,"_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseGRM.subset.mtx")
  # cat("write sparse GRM to ", sparseGRMFile ,"\n")
  # Matrix:::writeMM(sparseGRM, sparseGRMFile)
  Nval <- length(W)
  sparseSigma <- sparseGRM * tauVecNew[2]
  # diag(sparseSigma) = getDiagOfSigma(W, tauVecNew)
  diag(sparseSigma) <- diag(sparseSigma) + (1 / W) * tauVecNew[1]
  # sparseSigmaFile = paste0(outputPrefix, "_relatednessCutoff_",relatednessCutoff, "_", numRandomMarkerforSparseKin, "_randomMarkersUsed.sparseSigma.mtx")
  # cat("write sparse Sigma to ", sparseSigmaFile ,"\n")
  # Matrix:::writeMM(sparseSigma, sparseSigmaFile)
  sparseSigma
}


getsubGRM <- function(
  sparseGRMFile = NULL, sparseGRMSampleIDFile = "", relatednessCutoff,
  modelID = NULL
) {
  cat("extract sparse GRM\n")
  sparseGRMLarge <- Matrix:::readMM(sparseGRMFile)
  print(nnzero(sparseGRMLarge))
  cat(
    "set elements in the sparse GRM <= ", relatednessCutoff,
    " to zero\n"
  )
  sparseGRMLarge <- Matrix:::drop0(sparseGRMLarge, tol = relatednessCutoff)
  sparseGRMLarge <- sparseGRMLarge * 1
  print(nnzero(sparseGRMLarge))
  if (!file.exists(sparseGRMSampleIDFile)) {
    stop(
      "ERROR! sparseSigmaSampleIDFile ", sparseGRMSampleIDFile,
      " does not exist\n"
    )
  } else {
    sparseGRMSampleID <- data.frame(data.table:::fread(sparseGRMSampleIDFile,
      header = FALSE, stringsAsFactors = FALSE, colClasses = c("character")
    ))
    colnames(sparseGRMSampleID) <- c("sampleID")
    sparseGRMSampleID$IndexGRM <- seq(1, nrow(sparseGRMSampleID),
      by = 1
    )
    if (nrow(sparseGRMSampleID) != dim(sparseGRMLarge)[1] |
      nrow(sparseGRMSampleID) != dim(sparseGRMLarge)[2]) {
      stop(
        "ERROR! number of samples in the sparse GRM is not the same to the number of sample IDs in the specified sparseGRMSampleIDFile ",
        sparseGRMSampleIDFile, "\n"
      )
    } else {
      sampleInModel <- NULL
      sampleInModel$IID <- modelID
      sampleInModel <- data.frame(sampleInModel)
      sampleInModel$IndexInModel <- seq(1, length(sampleInModel$IID),
        by = 1
      )
      cat(nrow(sampleInModel), " samples have been used to fit the glmm null model\n")
      mergeID <- merge(sampleInModel, sparseGRMSampleID,
        by.x = "IID", by.y = "sampleID"
      )
      if (nrow(sampleInModel) > nrow(mergeID)) {
        stop(
          "ERROR: ", nrow(sampleInModel) - nrow(mergeID),
          "samples used for model fitting are not in the specified GRM\n"
        )
      } else {
        mergeID <- mergeID[with(mergeID, order(IndexInModel)), ]
        indexIDofGRM <- mergeID$IndexGRM
        sparseGRM <- sparseGRMLarge[indexIDofGRM, indexIDofGRM]
        rm(sparseGRMLarge)
        sparseGRM
      }
    }
  }
}


checkPerfectSep <- function(formula, data, minCovariateCount) {
  X1 <- model.matrix(formula, data = data)
  X_name <- colnames(X1)
  X1 <- as.matrix(X1[, -1])
  X_name <- X_name[-1]
  colnames(X1) <- X_name
  formula.frame <- model.frame(formula, data = data)
  Y <- model.response(formula.frame, type = "any")
  q <- length(X_name)
  colnamesDelete <- c()
  for (i in 1:q) {
    if (length(unique(X1[, i])) == 2) {
      sumTable <- table(Y, X1[, i])
      if (sum(sumTable < minCovariateCount) > 0) {
        colnamesDelete <- c(colnamesDelete, X_name[i])
        cat("less than ", minCovariateCount, " samples in a covariate  detected! ", X_name[i], " will be excluded in the model\n")
      }

      # if(sum(sumTable == 0) > 0){
      #    colnamesDelete = c(colnamesDelete, X_name[i])
      #    cat("perfect seperation is detected! ", X_name[i], " will be excluded in the model\n")
      #  }
    }
  }

  colnamesDelete
}

ScoreTest_NULL_Model_binary <- function(mu, y, X1) {
  V <- as.vector(mu * (1 - mu))
  res <- y - mu
  n1 <- length(res)
  XV <- t(X1 * V)
  XVX <- t(X1) %*% (X1 * V)
  XVX_inv <- solve(XVX)
  XXVX_inv <- X1 %*% XVX_inv
  S_a <- colSums(X1 * res)
  re <- list(XV = XV, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, XVX = XVX, S_a = S_a, XVX_inv_XV = XXVX_inv * V)
  class(re) <- "SA_NULL"
  re
}

extractVarianceRatio <- function(obj.glmm.null,
                                 seed = 1,
                                 obj.glm.null,
                                 maxiterPCG = 500,
                                 tolPCG = 0.01,
                                 numMarkers,
                                 varRatioOutFile,
                                 ratioCVcutoff,
                                 testOut,
                                 bedFile,
                                 bimFile,
                                 famFile,
                                 chromosomeStartIndexVec,
                                 chromosomeEndIndexVec,
                                 isCateVarianceRatio,
                                 cateVarRatioIndexVec,
                                 useSparseGRMforVarRatio,
                                 sparseGRMFile,
                                 sparseGRMSampleIDFile,
                                 numRandomMarkerforSparseKin,
                                 relatednessCutoff,
                                 useSparseGRMtoFitNULL,
                                 nThreads,
                                 cateVarRatioMinMACVecExclude,
                                 cateVarRatioMaxMACVecInclude,
                                 minMAFforGRM,
                                 isDiagofKinSetAsOne,
                                 includeNonautoMarkersforVarRatio,
                                 genoIn = NULL,
                                 # LEDGER #65: multiplies the null and null_noXadj
                                 # ratios before they are written. The estimation runs
                                 # on the standardised phenotype; those ratios'
                                 # numerator variance is in 1/y^2 units and their
                                 # denominator g'g carries no y, so 1/yScale^2 returns
                                 # them to phenotype units. The sparse ratio is NOT
                                 # scaled (LEDGER #97): both its sides go through
                                 # Sigma^-1, so it is scale-free. 1 for binary traits.
                                 varRatioScale = 1) {
  obj.noK <- obj.glmm.null$obj.noK
  if (is.null(genoIn)) genoIn <- genoInputPlink(bedFile, bimFile, famFile)
  if (file.exists(testOut)) {
    file.remove(testOut)
  }
  # Only column 1 (chromosome) of this table is ever used, at the CHR lookup
  # below; it is positionally aligned with the variant file.
  chrVecAll <- readVariantChromFromInput(genoIn)
  # if(sum(sapply(bimPlink[,1], is.numeric)) != nrow(bimPlink)){
  #  stop("ERROR: chromosome column in plink bim file is not numeric!\n")
  # }

  family <- obj.glm.null$family
  print(family)
  eta <- obj.glmm.null$linear.predictors
  mu <- obj.glmm.null$fitted.values
  X <- obj.glmm.null$X
  y <- obj.glm.null$y
  tauVecNew <- obj.glmm.null$theta

  mu.eta <- family$mu.eta(eta)
  sqrtW <- mu.eta / sqrt(obj.glm.null$family$variance(mu))
  W <- sqrtW^2 ## (mu*(1-mu) for binary)


  Sigma_iX_noLOCO <- getSigma_X(W, tauVecNew, X, maxiterPCG, tolPCG)


  ## randomize the marker orders to be tested
  # if(useSparseGRMtoFitNULL | useSparseGRMforVarRatio){
  if (useSparseGRMforVarRatio) {
    sparseSigma <- getSparseSigma(
      bedFile = bedFile, bimFile = bimFile, famFile = famFile, genoIn = genoIn,
      outputPrefix = varRatioOutFile,
      sparseGRMFile = sparseGRMFile,
      sparseGRMSampleIDFile = sparseGRMSampleIDFile,
      numRandomMarkerforSparseKin = numRandomMarkerforSparseKin,
      relatednessCutoff = relatednessCutoff,
      minMAFforGRM = minMAFforGRM,
      nThreads = nThreads,
      isDiagofKinSetAsOne = isDiagofKinSetAsOne,
      obj.glmm.null = obj.glmm.null,
      W = W, tauVecNew = tauVecNew
    )
  }
  mMarkers <- gettotalMarker()
  listOfMarkersForVarRatio <- list()
  MACvector <- getMACVec()
  isVarianceRatioinGeno <- getIsVarRatioGeno()

  if (isVarianceRatioinGeno) {
    MACvector_forVarRatio <- getMACVec_forVarRatio()
    Indexvector_forVarRatio <- getIndexVec_forVarRatio()
    cat("length(MACvector): ", length(MACvector), "\n")
    cat("length(MACvector_forVarRatio): ", length(MACvector_forVarRatio), "\n")

    if (length(MACvector_forVarRatio) > 0) {
      # MACdata = data.frame(MACvector = c(MACvector, MACvector_forVarRatio), geno_ind = c(rep(0, length(MACvector)), rep(1, length(MACvector_forVarRatio))), indexInGeno = c(seq(1,length(MACvector)), seq(1,length(MACvector_forVarRatio))))
      MACdata <- data.frame(MACvector = MACvector_forVarRatio, geno_ind = rep(1, length(MACvector_forVarRatio)), indexInGeno = seq(1, length(MACvector_forVarRatio)))
    } else {
      stop("No markers were found for variance ratio estimation. Please make sure there are at least 200 markers in each MAC category\n")
    }
  } else {
    MACdata <- data.frame(MACvector = MACvector, geno_ind = rep(0, length(MACvector)), indexInGeno = seq(1, length(MACvector)))
  }

  # LEDGER #8 fix: reseed here so variance-ratio marker selection depends
  # only on --seed, not on how many trace/AI-REML iterations reset and
  # advanced the global RNG stream (GetTrace calls set_seed(200) internally).
  set.seed(seed + 1)

  if (!isCateVarianceRatio) {
    cat("Only one variance ratio will be estimated using randomly selected markers with MAC >= 20\n")
    MACindex <- 1:nrow(MACdata)
    listOfMarkersForVarRatio[[1]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    cateVarRatioIndexVec <- c(1)
  } else {
    cat("Categorical variance ratios will be estimated.\n")

    if (is.null(cateVarRatioIndexVec)) {
      cateVarRatioIndexVec <- rep(1, length(cateVarRatioMinMACVecExclude))
    }
    numCate <- length(cateVarRatioIndexVec)
    for (i in seq_len(numCate - 1)) { # LEDGER #35 FIXED: 1:0 ran twice and indexed [[0]]
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[i] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[i])
      listOfMarkersForVarRatio[[i]] <- sample(MACindex, size = length(MACindex), replace = FALSE)
    }

    if (length(cateVarRatioMaxMACVecInclude) == (numCate - 1)) {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate])
    } else {
      MACindex <- which(MACdata$MACvector > cateVarRatioMinMACVecExclude[numCate] & MACdata$MACvector <= cateVarRatioMaxMACVecInclude[numCate])
    }

    listOfMarkersForVarRatio[[numCate]] <- sample(MACindex, size = length(MACindex), replace = FALSE)

    for (k in 1:length(cateVarRatioIndexVec)) {
      if (k <= length(cateVarRatioIndexVec) - 1) {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], "\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC <= ", cateVarRatioMaxMACVecInclude[k], " is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      } else {
        if (cateVarRatioIndexVec[k] == 1) {
          cat(cateVarRatioMinMACVecExclude[k], "< MAC\n")
          if (length(listOfMarkersForVarRatio[[k]]) < numMarkers) {
            stop("ERROR! number of genetic variants in ", cateVarRatioMinMACVecExclude[k], "< MAC  is lower than ", numMarkers, "\n", "Please include more markers in this MAC category in the plink file\n")
          }
        }
      }
    }
  } # if(!isCateVarianceRatio){

  freqVec <- getAlleleFreqVec()

  Nnomissing <- length(mu)
  varRatioTable <- NULL
  # LEDGER #83 FIXED. Whether a `sparse` row is emitted is a property of the
  # RUN (are sparse-GRM ratios being estimated at all?), not of whatever the
  # previous category left in varRatio_sparseGRM_vec -- which was NULL-ed only
  # inside the estimated arm, so a skipped category read the category before
  # it, or an unbound symbol at k = 1.
  emitSparseRow <- isTRUE(useSparseGRMforVarRatio)
  for (k in 1:length(listOfMarkersForVarRatio)) {
    varRatio_sparseGRM_vec <- NULL
    if (cateVarRatioIndexVec[k] == 1) {
      numMarkers0 <- numMarkers
      varRatio_NULL_vec <- NULL
      varRatio_NULL_noXadj_vec <- NULL
      indexInMarkerList <- 1
      numTestedMarker <- 0
      ratioCV <- ratioCVcutoff + 0.1

      while (ratioCV > ratioCVcutoff) {
        while (numTestedMarker < numMarkers0) {
          macdata_i <- listOfMarkersForVarRatio[[k]][indexInMarkerList]
          i <- (MACdata$indexInGeno)[macdata_i]
          genoInd <- (MACdata$geno_ind)[macdata_i]
          cat(i, "th marker in geno ", genoInd, "\n")
          cat("MAC: ", (MACdata$MACvector)[macdata_i], "\n")
          if (genoInd == 0) {
            G0 <- Get_OneSNP_Geno(i - 1)
          } else if (genoInd == 1) {
            G0 <- Get_OneSNP_Geno_forVarRatio(i - 1)
          }

          cat("G0", G0[1:10], "\n")
          CHR <- chrVecAll[Indexvector_forVarRatio[i] + 1]
          CHR <- gsub("(?i)^chr", "", CHR, perl = TRUE)
          CHR <- as.numeric(CHR)
          cat("CHR ", CHR, "\n")
          if (sum(G0) / (2 * Nnomissing) > 0.5) {
            G0 <- 2 - G0
          }
          NAset <- which(G0 == 0)
          AC <- sum(G0)

          indexInMarkerList <- indexInMarkerList + 1
          if (!is.na(CHR)) {
            if ((CHR >= 1 & CHR <= 22) | includeNonautoMarkersforVarRatio) {
              AF <- AC / (2 * Nnomissing)
              if (CHR >= 1 & CHR <= 22) {
                autoMarker <- TRUE
              } else {
                autoMarker <- FALSE
              }
              G <- G0 - obj.noK$XXVX_inv %*% (obj.noK$XV %*% G0) # G1 is X adjusted
              g <- G / sqrt(AC)
              q <- innerProduct(g, y)
              # eta = obj.glmm.null$linear.predictors
              # mu = obj.glmm.null$fitted.values
              # mu.eta = family$mu.eta(eta)
              # sqrtW = mu.eta/sqrt(obj.glm.null$family$variance(mu))
              # W = sqrtW^2
              G_noXadj <- as.vector(G0 - mean(G0))
              g_noXadj <- G_noXadj / sqrt(AC)


              Sigma_iG <- getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)

              Sigma_iX <- Sigma_iX_noLOCO

              var1a <- t(G) %*% Sigma_iG - t(G) %*% Sigma_iX %*% (solve(t(X) %*% Sigma_iX)) %*% t(X) %*% Sigma_iG
              var1 <- var1a / AC
              m1 <- innerProduct(mu, g)

              if (!useSparseGRMtoFitNULL) {
                if (useSparseGRMforVarRatio) {
                  pcginvSigma <- solve(sparseSigma, g, sparse = TRUE)
                  var2_a <- t(g) %*% pcginvSigma
                  var2sparseGRM <- var2_a[1, 1]
                  varRatio_sparseGRM_vec <- c(varRatio_sparseGRM_vec, var1 / var2sparseGRM)
                }
              } else {
                if (useSparseGRMforVarRatio) {
                  # varRatio_sparseGRM_vec = c(varRatio_sparseGRM_vec, 1)
                  # pcginvSigma = solve(sparseSigma, g, sparse=T)
                  pcginvSigma <- Sigma_iG / sqrt(AC)
                  var2_a <- t(g) %*% pcginvSigma
                  var2sparseGRM <- var2_a[1, 1]
                  # x=t(G)%*%Sigma_iG/AC
                  # cat(" x ", x, " var2 ", var2sparseGRM , "\n")
                  varRatio_sparseGRM_vec <- c(varRatio_sparseGRM_vec, var1 / var2sparseGRM)
                }
              }

              if (obj.glmm.null$traitType == "binary") {
                var2null <- innerProduct(mu * (1 - mu), g * g)
                var2null_noXadj <- innerProduct(mu * (1 - mu), g_noXadj * g_noXadj)
              } else if (obj.glmm.null$traitType == "quantitative") {
                var2null <- innerProduct(g, g)
                var2null_noXadj <- innerProduct(g_noXadj, g_noXadj)
              }

              varRatio_NULL_vec <- c(varRatio_NULL_vec, var1 / var2null)
              varRatio_NULL_noXadj_vec <- c(varRatio_NULL_noXadj_vec, var1 / var2null_noXadj)
              # indexInMarkerList = indexInMarkerList + 1
              numTestedMarker <- numTestedMarker + 1
            }
            # LEDGER #75 FIXED: this `else` used to advance indexInMarkerList a
            # SECOND time, on top of the unconditional increment above. The
            # commented-out increment inside the accepted branch is the record of
            # what happened: the increment was hoisted out to the top of the
            # iteration and these two copies were left behind. A skipped marker
            # therefore consumed the marker AFTER it, which was never read at
            # all, and a double step from the last position jumped straight over
            # the exhaustion guard below.
          }
          # LEDGER #75 FIXED: the same second increment sat in the branch taken
          # when the chromosome does not parse; it is gone for the same reason.


          # LEDGER #75: the exhaustion guard below is now an inequality rather
          # than an equality. With the two double increments gone
          # the cursor advances by exactly one per iteration, so equality is
          # sufficient and this changes nothing today -- it is here so that any
          # future path that advances by more than one stops at the end of the
          # list instead of reading an NA off it.
          if (indexInMarkerList - 1 >= length(listOfMarkersForVarRatio[[k]])) {
            numTestedMarker <- numMarkers0
          }
        } # end of while(numTestedMarker < numMarkers)

        ratioCV <- calCV(varRatio_NULL_vec)

        if (ratioCV > ratioCVcutoff) {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " > ", ratioCVcutoff, "\n")
          numMarkers0 <- numMarkers0 + 10
          cat("try ", numMarkers0, " markers\n")
        } else {
          cat("CV for variance ratio estimate using ", numMarkers0, " markers is ", ratioCV, " < ", ratioCVcutoff, "\n")
        }

        if (indexInMarkerList - 1 >= length(listOfMarkersForVarRatio[[k]])) {
          ratioCV <- ratioCVcutoff
          cat("no more markers are available in the MAC category ", k, "\n")
          print(indexInMarkerList - 1)
        }
      } # end of while(ratioCV > ratioCVcutoff)

      if (length(varRatio_sparseGRM_vec) > 0) {
        varRatio_sparse <- mean(varRatio_sparseGRM_vec)
        cat("varRatio_sparse", varRatio_sparse, "\n")
        # LEDGER #97 FIXED: the sparse ratio is scale-free (numerator AND
        # denominator go through Sigma^-1), so it is NOT multiplied by varRatioScale
        varRatioTable <- rbind(varRatioTable, c(varRatio_sparse, "sparse", k))
      }
      varRatio_null <- mean(varRatio_NULL_vec)
      varRatio_null_noXadj <- mean(varRatio_NULL_noXadj_vec)
      cat("varRatio_null", varRatio_null, "\n")
      cat("varRatio_null_noXadj", varRatio_null_noXadj, "\n")
      varRatioTable <- rbind(varRatioTable, c(varRatio_null * varRatioScale, "null", k))
      varRatioTable <- rbind(varRatioTable,
                             c(varRatio_null_noXadj * varRatioScale, "null_noXadj", k))
    } else { # if(cateVarRatioVec[k] == 1)
      varRatioTable <- rbind(varRatioTable, c(1 * varRatioScale, "null", k))
      varRatioTable <- rbind(varRatioTable, c(1 * varRatioScale, "null_noXadj", k))
      if (emitSparseRow) {
        varRatioTable <- rbind(varRatioTable, c(1, "sparse", k)) # LEDGER #97: scale-free
      }
    }
  } # for(k in 1:length(listOfMarkersForVarRatio)){
  write.table(varRatioTable, varRatioOutFile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  data <- read.table(varRatioOutFile, header = FALSE)
  print(data)
}
