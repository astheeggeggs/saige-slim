//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double Korg_Binom(double t1, arma::vec & mu, arma::vec & g);
double K1_adj_Binom(double t1, arma::vec & mu, arma::vec & g, double q);
double K2_Binom(double t1, arma::vec & mu, arma::vec & g);
// gpos/gneg are the sums of the positive and negative entries of g. They depend
// only on g, but SPA()/SPA_pval()/SPA_fast() each call getroot twice with the
// SAME g (once for q, once for qinv), so computing them inside getroot did the
// two full-length find+elem+accu passes twice per SPA marker. The _pre variants
// take them precomputed; the exported wrappers below keep the old signature and
// compute them, so nothing R-visible changes.
void getGPosNeg_Binom(arma::vec & g, double & gpos, double & gneg);
Rcpp::List getroot_K1_Binom_pre(double init, arma::vec & mu, arma::vec & g, double q, double tol, double gpos, double gneg, int maxiter = 1000);
Rcpp::List getroot_K1_Binom(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter = 1000);
Rcpp::List Get_Saddle_Prob_Binom(double zeta,  arma::vec & mu, arma::vec & g, double q, bool logp=false);
Rcpp::List SPA_binary(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp = false);
double Korg_fast_Binom(double t1, arma::vec & mu, arma::vec & g,  arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);
double K1_adj_fast_Binom(double t1, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma);
double K2_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma);
Rcpp::List getroot_K1_fast_Binom_pre(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, double gpos, double gneg, int maxiter = 1000);
Rcpp::List getroot_K1_fast_Binom(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, int maxiter = 1000);
Rcpp::List Get_Saddle_Prob_fast_Binom(double zeta,  arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, bool logp=false);
Rcpp::List SPA_binary_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol);
