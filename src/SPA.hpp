//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


void SPA(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge);

void SPA_fast(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, bool logp, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge);

double SPA_pval(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, bool & isSPAConverge);
