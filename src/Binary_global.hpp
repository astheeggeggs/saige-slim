#ifndef GLOBAL_001 
#define GLOBAL_001 

#include <vector>
#include <stdexcept>

// S6: the ER exact test cannot run on a worker thread, and the reason is
// broader than the RNG.
//
// Two separate R dependencies live under SKATExactBin_Work:
//   1. SL_setseed, on the resampling arms only. It goes through
//      Rcpp::Environment/Rcpp::Function/GetRNGstate (the LEDGER #2 fix); S2
//      made it conditional on an arm actually drawing.
//   2. Rf_lchoose, on EVERY ER marker, via HyperGeo::lCombinations -- **REMOVED**.
//      It was never the arithmetic: R's stack check measures against the MAIN
//      thread's stack base, so a worker reads as impossibly deep and R aborts.
//      choose()/lchoose() recurse and carry that check; Rf_lgammafn, Rf_pnorm5,
//      Rf_qnorm5, Rf_pchisq and Rf_lbeta are all FINE off-thread, which is why
//      the SPA path always threaded and only ER did not. lchoose is now computed
//      from Rf_lbeta the way R computes it (HyperGeo::lchoose_threadsafe),
//      verified bit-identical over 2,009,004 argument pairs and on a worker
//      thread with the stack check left ON. So no R_CStackLimit games, and no
//      golden moves.
//
// What remains main-thread-only is therefore SL_setseed, on the RESAMPLING arms
// only, guarded where it lives in Binary_ComputeExact.cpp. The marker loop's
// predicate still defers every binary MAC <= max_MAC_for_ER marker, which is now
// more conservative than it needs to be; narrowing it is the next step and is
// gated on TSAN, not on this comment.
// This is the mechanism for that handback, and it is deliberately the
// AUTHORITY rather than a hint: the loop also has a cheap predicate (binary and
// MAC <= max_MAC_for_ER) for which markers take the branch, but a predicate can
// rot. With refuse() set, the ER path throws NeedsMainThread at the instant it
// reaches either R dependency, so a predicate that is wrong costs a re-run
// rather than an aborted session.
namespace SKATMainThread {

struct NeedsMainThread : public std::runtime_error {
  explicit NeedsMainThread(const char* t_what) : std::runtime_error(t_what) {}
};

void setRefuse(bool t_refuse);   // thread-local
bool refuse();

}  // namespace SKATMainThread

void * SL_calloc(size_t num, size_t size);
void SL_free(void * ptr);
void SL_Sample(int k, int n, std::vector<int> & y, std::vector<int> & x);


void SL_Binary_Boot1(int n, int ncase, std::vector<double> & pcase, std::vector<int> & buf1, std::vector<int> & buf2, std::vector<int> & z_one, int *err);
double SL_runif_double();
int SL_runif_INT(int max);
void  SL_setseed(int seed);
void  SL_out();
void SL_GetSample(int n, int k, std::vector<int> & y, std::vector<int> & x);
void SL_GetPermu(int n, std::vector<int> & y, std::vector<int> & x);

void SKAT_Exact(std::vector<int> & resarray, int nres, int * nres_k, double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon);

void SKATO_Exact(int * resarray, int nres, int * nres_k, double * Z0, double *Z1, double * r_corr, int n_r,
                 double * param, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon);

#endif 

