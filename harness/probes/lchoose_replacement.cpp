// Can Rf_lchoose be replaced, BIT FOR BIT, by a form that does not enter the
// recursive choose.c entry point (and so does not trip R's C-stack check on a
// worker thread)?
//
// On the domain HyperGeo::lCombinations actually uses -- integer n >= 0 and
// integer k, with k > n already short-circuited by the caller -- R's lchoose
// reduces to early returns plus
//     lfastchoose(n,k) = -log(n+1) - lbeta(n-k+1, k+1)
// and Rf_lbeta was measured safe off-thread and bit-identical there.
// This transcribes R's own branch structure and compares exhaustively.
#include <thread>
#include <vector>
#include <cmath>
#include <cstdint>
#include <cstring>
#define R_NO_REMAP
#define STRICT_R_HEADERS 1
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

static double lfastchoose_(double n, double k) {
  return -std::log(n + 1.) - Rf_lbeta(n - k + 1., k + 1.);
}

// Faithful to R's choose.c for INTEGER n >= 0. The two recursive arms of the
// original are unrolled: n < 0 cannot occur here, and the symmetry arm
// lchoose(n, n-k) lands on n-k in {0,1}, which are the k<2 early returns.
static double saige_lchoose(double n, double k) {
  if (k < 2) {
    if (k < 0) return R_NegInf;
    if (k == 0) return 0.;
    return std::log(std::fabs(n));          // k == 1
  }
  if (n < k) return R_NegInf;
  if (n - k < 2) {                          // symmetry: n-k is 0 or 1
    double k2 = n - k;
    if (k2 == 0) return 0.;
    return std::log(std::fabs(n));
  }
  return lfastchoose_(n, k);
}

extern "C" SEXP compare_exhaustive(SEXP nmaxS) {
  int nmax = Rf_asInteger(nmaxS);
  long long total = 0, diff = 0; int fn = -1, fk = -1;
  double fa = 0, fb = 0;
  for (int n = 0; n <= nmax; n++) {
    for (int k = -3; k <= n; k++) {
      double a = Rf_lchoose((double)n, (double)k);
      double b = saige_lchoose((double)n, (double)k);
      total++;
      // bit comparison, so NaN/Inf/-0.0 are all held to the same standard
      if (std::memcmp(&a, &b, sizeof(double)) != 0) {
        if (diff == 0) { fn = n; fk = k; fa = a; fb = b; }
        diff++;
      }
    }
  }
  SEXP r = PROTECT(Rf_allocVector(REALSXP, 6));
  REAL(r)[0]=(double)total; REAL(r)[1]=(double)diff;
  REAL(r)[2]=fn; REAL(r)[3]=fk; REAL(r)[4]=fa; REAL(r)[5]=fb;
  UNPROTECT(1); return r;
}

// The point of the exercise: run the replacement on a WORKER thread with R's
// C-stack check left ON (which is what aborts Rf_lchoose), and check both that
// it survives and that it agrees bit for bit with the main thread.
static std::vector<double> g_worker;
extern "C" SEXP offthread_check(SEXP nmaxS) {
  int nmax = Rf_asInteger(nmaxS);
  std::vector<double> main_vals;
  for (int n = 0; n <= nmax; n++)
    for (int k = -3; k <= n; k++) main_vals.push_back(saige_lchoose(n, k));
  g_worker.clear();
  std::thread t([nmax](){
    for (int n = 0; n <= nmax; n++)
      for (int k = -3; k <= n; k++) g_worker.push_back(saige_lchoose(n, k));
  });
  t.join();
  long long diff = 0;
  for (size_t i = 0; i < main_vals.size(); i++)
    if (std::memcmp(&main_vals[i], &g_worker[i], sizeof(double)) != 0) diff++;
  SEXP r = PROTECT(Rf_allocVector(REALSXP, 2));
  REAL(r)[0] = (double)main_vals.size(); REAL(r)[1] = (double)diff;
  UNPROTECT(1); return r;
}
