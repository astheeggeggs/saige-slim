#include <thread>
#include <cstdio>
#include <cstdint>
#define R_NO_REMAP
#define STRICT_R_HEADERS 1
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
// Rinterface.h does not expose this in a package build; it is an exported
// symbol of libR on unix, so declare it directly.
extern "C" uintptr_t R_CStackLimit;

static void run(double n, double k, double* out, int* ok) {
  *ok = 0;
  *out = Rf_lchoose(n, k);
  *ok = 1;
}

extern "C" SEXP probe_lchoose(SEXP disable) {
  uintptr_t saved = R_CStackLimit;
  if (Rf_asLogical(disable)) R_CStackLimit = (uintptr_t)-1;
  const double N[4] = {1,1,5,999}, K[4] = {0,1,2,98};
  SEXP res = PROTECT(Rf_allocVector(REALSXP, 8));
  for (int i = 0; i < 4; i++) {
    double main_val = Rf_lchoose(N[i], K[i]);
    double thr_val = -999; int ok = 0;
    std::thread t(run, N[i], K[i], &thr_val, &ok);
    t.join();
    REAL(res)[2*i] = main_val;
    REAL(res)[2*i+1] = ok ? thr_val : NA_REAL;
  }
  R_CStackLimit = saved;
  UNPROTECT(1);
  return res;
}
