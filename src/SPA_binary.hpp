//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Plain-old-data results for the saddlepoint kernels (S1 in
// docs/design/step2-threading-verdict.md).
//
// The two root solves and two saddle evaluations that every SPA-corrected
// marker performs used to hand their results back in Rcpp::List. An Rcpp::List
// has no in-place named insert -- `result["root"] = root` reallocates the list
// and its names -- so each three-field result cost several R-heap allocations,
// and reading a field back cost an R API call. That is ~20 R allocations per
// SPA marker, each a potential GC over a heap holding the n x p model matrices,
// and it is why the SPA kernels cannot run on a worker thread: every Rcpp::List
// touch is an R API call.
//
// The `_core` functions below carry these two structs instead and use no R
// types at all. The exported wrappers keep their exact signatures and still
// return Rcpp::List, so RcppExports and everything R-visible are unchanged.
struct SPARoot {
        double root;
        int niter;
        bool Isconverge;
};

struct SPASaddle {
        double pval;
        bool isSaddle;
};

double Korg_Binom(double t1, arma::vec & mu, arma::vec & g);
double K1_adj_Binom(double t1, arma::vec & mu, arma::vec & g, double q);
double K2_Binom(double t1, arma::vec & mu, arma::vec & g);
// gpos/gneg are the sums of the positive and negative entries of g. They depend
// only on g, but SPA()/SPA_pval()/SPA_fast() each call getroot twice with the
// SAME g (once for q, once for qinv), so computing them inside getroot did the
// two full-length find+elem+accu passes twice per SPA marker. The _core variants
// take them precomputed; the exported wrappers below keep the old signature and
// compute them, so nothing R-visible changes.
void getGPosNeg_Binom(arma::vec & g, double & gpos, double & gneg);
SPARoot getroot_K1_Binom_core(double init, arma::vec & mu, arma::vec & g, double q, double tol, double gpos, double gneg, int maxiter = 1000);
Rcpp::List getroot_K1_Binom(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter = 1000);
SPASaddle Get_Saddle_Prob_Binom_core(double zeta,  arma::vec & mu, arma::vec & g, double q, bool logp);
Rcpp::List Get_Saddle_Prob_Binom(double zeta,  arma::vec & mu, arma::vec & g, double q, bool logp=false);
Rcpp::List SPA_binary(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp = false);
double Korg_fast_Binom(double t1, arma::vec & mu, arma::vec & g,  arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);
double K1_adj_fast_Binom(double t1, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma);
double K2_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma);
SPARoot getroot_K1_fast_Binom_core(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, double gpos, double gneg, int maxiter = 1000);
Rcpp::List getroot_K1_fast_Binom(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol, int maxiter = 1000);
SPASaddle Get_Saddle_Prob_fast_Binom_core(double zeta,  arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, bool logp);
Rcpp::List Get_Saddle_Prob_fast_Binom(double zeta,  arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, bool logp=false);
Rcpp::List SPA_binary_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB,  double NAmu, double NAsigma, double tol);
