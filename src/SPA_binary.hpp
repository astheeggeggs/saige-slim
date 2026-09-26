#ifndef SPA_BINARY_HPP
#define SPA_BINARY_HPP
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

// Scratch for the saddlepoint kernels. Every Korg/K1/K2 evaluation used to
// build three or four n-length Armadillo temporaries and free them on return:
// ~80 allocations per dense-SPA marker at N = 50,000, each above libmalloc's
// large threshold and therefore a mach_vm map/unmap pair that page-faults on
// first touch and serialises on the task's VM lock under --nThreads (measured
// 2026-09-05: sys time 19 s -> 61 s at 4 threads on the binary scan, with
// mach_vm_deallocate 6.8x and bzero 2.9x while pure-compute leaves grew 1.4x).
// The four buffers below are sized once and the _ws kernels write into strict
// non-owning views of their first m elements, so a fast kernel over the m
// non-zero-genotype samples and a dense one over all n share the same memory.
// Values are unchanged: every expression is the original verbatim, only its
// destination moved. The exported kernels keep their signatures and build a
// local workspace, so nothing R-visible changes.
struct SPAWorkspace {
        arma::vec a, b, c, d;
        // LEDGER #164: why the last SPA call did not converge (0 = it did):
        // 1 root(q) not found, 2 root(mirrored q) not found, 3 saddle(q)
        // invalid, 4 saddle(mirrored q) invalid. Set by SPA_ws/SPA_fast_ws;
        // read by getMarkerPval for the fallback sidecar.
        int failReason = 0;
        SPAWorkspace() {}
        explicit SPAWorkspace(arma::uword n) { ensure(n); }
        void ensure(arma::uword n) {
                if(a.n_elem < n){ a.set_size(n); b.set_size(n); c.set_size(n); d.set_size(n); }
        }
};
double Korg_Binom_ws(double t1, const arma::vec & mu, const arma::vec & g, SPAWorkspace & ws);
double K1_adj_Binom_ws(double t1, const arma::vec & mu, const arma::vec & g, double q, SPAWorkspace & ws);
double K2_Binom_ws(double t1, const arma::vec & mu, const arma::vec & g, SPAWorkspace & ws);
double Korg_fast_Binom_ws(double t1, const arma::vec & gNB, const arma::vec & muNB, double NAmu, double NAsigma, SPAWorkspace & ws);
double K1_adj_fast_Binom_ws(double t1, double q, const arma::vec & gNB, const arma::vec & muNB, double NAmu, double NAsigma, SPAWorkspace & ws);
double K2_fast_Binom_ws(double t1, const arma::vec & gNB, const arma::vec & muNB, double NAsigma, SPAWorkspace & ws);

double Korg_Binom(double t1, const arma::vec & mu, const arma::vec & g);
double K1_adj_Binom(double t1, const arma::vec & mu, const arma::vec & g, double q);
double K2_Binom(double t1, const arma::vec & mu, const arma::vec & g);
// gpos/gneg are the sums of the positive and negative entries of g. They depend
// only on g, but SPA()/SPA_pval()/SPA_fast() each call getroot twice with the
// SAME g (once for q, once for qinv), so computing them inside getroot did the
// two full-length find+elem+accu passes twice per SPA marker. The _core variants
// take them precomputed; the exported wrappers below keep the old signature and
// compute them, so nothing R-visible changes.
void getGPosNeg_Binom(const arma::vec & g, double & gpos, double & gneg);
SPARoot getroot_K1_Binom_core(double init, const arma::vec & mu, const arma::vec & g, double q, double tol, double gpos, double gneg, SPAWorkspace & ws, int maxiter = 1000);
Rcpp::List getroot_K1_Binom(double init, const arma::vec & mu, const arma::vec & g, double q, double tol, int maxiter = 1000);
SPASaddle Get_Saddle_Prob_Binom_core(double zeta,  const arma::vec & mu, const arma::vec & g, double q, bool logp, SPAWorkspace & ws);
Rcpp::List Get_Saddle_Prob_Binom(double zeta,  const arma::vec & mu, const arma::vec & g, double q, bool logp=false);
Rcpp::List SPA_binary(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp = false);
double Korg_fast_Binom(double t1, const arma::vec & mu, const arma::vec & g,  const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB, double NAmu, double NAsigma);
double K1_adj_fast_Binom(double t1, const arma::vec & mu, const arma::vec & g, double q, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma);
double K2_fast_Binom(double t1, const arma::vec & mu, const arma::vec & g, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma);
SPARoot getroot_K1_fast_Binom_core(double init, const arma::vec & g, double q, const arma::vec & gNB, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, double gpos, double gneg, SPAWorkspace & ws, int maxiter = 1000);
Rcpp::List getroot_K1_fast_Binom(double init, const arma::vec & mu, const arma::vec & g, double q, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, int maxiter = 1000);
SPASaddle Get_Saddle_Prob_fast_Binom_core(double zeta,  const arma::vec & g, double q, const arma::vec & gNB, const arma::vec & muNB,  double NAmu, double NAsigma, bool logp, SPAWorkspace & ws);
Rcpp::List Get_Saddle_Prob_fast_Binom(double zeta,  const arma::vec & mu, const arma::vec & g, double q, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, bool logp=false);
Rcpp::List SPA_binary_fast(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, bool logp, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, double tol);

#endif // SPA_BINARY_HPP
