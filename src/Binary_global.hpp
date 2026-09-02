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
//   2. **Rf_lchoose**, on EVERY ER marker, via HyperGeo::lCombinations. This
//      one is not conditional on anything and was not in F4's inventory.
//      MEASURED, in a 20-line program with no SAIGE in it: calling
//      Rf_lchoose from a std::thread aborts the session with
//      "C stack usage ... is too close to the limit", for every argument pair
//      tried -- (1,0), (1,1), (5,2), (999,98) -- because R's stack check
//      measures against the MAIN thread's stack base. The same probe shows
//      Rf_lgammafn, Rf_pnorm5, Rf_qnorm5 and Rf_pchisq are all FINE off-thread,
//      which is why the SPA path threads happily and only ER does not.
//
// So the whole ER branch is main-thread-only, and the marker loop defers it.
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

