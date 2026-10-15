# Code adpated from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value
#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{CCT} function takes in a numeric vector of p-values, a numeric
#' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
#' equal weights are assumed.
#' @return the aggregated p-value combining p-values from the vector \code{pvals}.
#' @examples pvalues <- c(2e-02, 4e-04, 0.2, 0.1, 0.8)
#' @examples CCT(pvals = pvalues)
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association 115}(529), 393-402.
#' (\href{https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485}{pub})
#' @export

CCT <- function(pvals, weights = NULL) {
  #### check if there is NA
  if (sum(is.na(pvals)) > 0) {
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals == 0) >= 1)
  is.one <- (sum(pvals == 1) >= 1)
  # if(is.zero && is.one){
  #  stop("Cannot have both 0 and 1 p-values!")
  # }
  if (is.zero) {
    return(0)
  }
  if (is.one) {
    # warning("There are p-values that are exactly 1!")
    return(min(1, (min(pvals)) * (length(pvals))))
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if (is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))
  } else if (length(weights) != length(pvals)) {
    stop("The length of weights should be the same as that of the p-values!")
  } else if (sum(weights < 0) > 0) {
    stop("All the weights must be positive!")
  } else {
    weights <- weights / sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
  } else {
    cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
    cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }

  #### check if the test statistic is very large.
  if (cct.stat > 1e+15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  pval
}

#' The Cauchy combination in log space: log10 of CCT()'s p-value, from log10
#' p-values, with equal weights. Written for --is_output_log10P (2026-10-14):
#' a cell whose p underflowed double precision reaches CCT() as 0 and the
#' combined p is then 0 by its `is.zero` rule, which is a statement about the
#' representation, not the evidence. Here every term is kept as a logarithm.
#' Same rules as CCT(): a p of exactly 1 anywhere gives min(1, k * min p); a
#' p below 1e-16 contributes 1/(pi p) (the tangent's asymptote) and is summed
#' as a log-sum-exp; the rest contribute tan((0.5 - p) pi) linearly; a
#' statistic above 1e15 is inverted as 1/(pi T), else the Cauchy CDF. NA
#' inputs are dropped; NA out when nothing is left. Agrees with log10(CCT(p))
#' wherever CCT(p) is representable (test-step2-log10p.R).
CCT_log10 <- function(log10p) {
  log10p <- log10p[!is.na(log10p)]
  k <- length(log10p)
  if (k == 0) {
    return(NA_real_)
  }
  if (any(log10p > 0)) {
    stop("All p-values must be between 0 and 1!")
  }
  # CCT()'s is.zero rule, first as there: a p of exactly 0 (log10 -Inf) makes
  # the combination 0. Without this the log-sum-exp below computes Inf - Inf =
  # NaN and the `if` on it errors, which on AoU (2026-10-15, GPT for ALT) killed
  # a whole chromosome from inside region_postprocess, outside the per-cell
  # kernel catch.
  if (any(is.infinite(log10p))) {
    return(-Inf)
  }
  # CCT()'s is.one rule: min(1, k * min p)
  if (any(log10p == 0)) {
    return(min(0, min(log10p) + log10(k)))
  }
  small <- log10p < -16
  lin <- if (any(!small)) sum(tan((0.5 - 10^log10p[!small]) * pi)) else 0
  # the Cauchy upper tail without cancellation: 0.5 - atan(t)/pi loses digits
  # for t above ~1e8 (measured 0.5 % at t ~ 1e13); atan(1/t)/pi is the same
  # quantity exactly for t > 0
  tail10 <- function(stat) {
    if (stat > 1e+15) return(-(log10(stat) + log10(pi)))
    if (stat > 0) return(log10(atan(1 / stat) / pi))
    log10(0.5 - atan(stat) / pi)
  }
  if (!any(small)) {
    return(tail10(lin / k))
  }
  big <- -log10p[small] - log10(pi)          # log10 of each 1/(pi p) term
  m <- max(big)
  log10S <- m + log10(sum(10^(big - m)))      # log-sum-exp, base 10
  if (lin != 0 && log10S < 300) {
    ratio <- lin / 10^log10S                  # the linear terms, relative to the large ones
    if (1 + ratio <= 0) {
      # cannot happen with a term at or above 1e16/pi against tangents bounded by
      # k * tan(0.5 pi - 1e-16) ~ k * 1e16: kept as a guard, not a path
      return(tail10((10^log10S + lin) / k))
    }
    log10S <- log10S + log10(1 + ratio)
  }
  log10T <- log10S - log10(k)
  if (log10T > 15) {
    return(-(log10T + log10(pi)))
  }
  tail10(10^log10T)
}
