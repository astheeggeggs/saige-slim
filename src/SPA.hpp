#ifndef SPA_HPP
#define SPA_HPP
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


#include "SPA_binary.hpp"

void SPA(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge);
// The same, over a caller-owned SPAWorkspace (see SPA_binary.hpp). getMarkerPval
// calls these; the exported SPA/SPA_fast above build a local workspace and
// delegate, so their signatures and results are unchanged.
void SPA_ws(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge, SPAWorkspace & ws);
void SPA_fast_ws(const arma::vec & g, double q, double qinv, double pval_noadj, bool logp, const arma::vec & gNB, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge, SPAWorkspace & ws);

void SPA_fast(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, bool logp, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge);

double SPA_pval(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, bool & isSPAConverge);

#endif // SPA_HPP
