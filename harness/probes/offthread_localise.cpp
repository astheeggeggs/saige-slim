// Localise WHICH call inside the ER path trips R's C-stack check on a worker
// thread. Each candidate is a separate entry point so it can be run in its own
// process: the failure mode is an R error raised on a non-main thread, which is
// not something to recover from and share a session with.
#include <thread>
#include <cstdint>
#define R_NO_REMAP
#define STRICT_R_HEADERS 1
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
extern "C" uintptr_t R_CStackLimit;

static double g_out;
#define ENTRY(name, expr)                                    \
  extern "C" SEXP name(SEXP disable) {                       \
    uintptr_t saved = R_CStackLimit;                         \
    if (Rf_asLogical(disable)) R_CStackLimit = (uintptr_t)-1;\
    double on_main = (expr);                                 \
    std::thread t([](){ g_out = (expr); });                  \
    t.join();                                                \
    R_CStackLimit = saved;                                   \
    SEXP r = PROTECT(Rf_allocVector(REALSXP, 2));            \
    REAL(r)[0] = on_main; REAL(r)[1] = g_out;                \
    UNPROTECT(1); return r;                                  \
  }

ENTRY(p_noop,      0.0)                       // negative control: empty worker
ENTRY(p_lgammafn,  Rf_lgammafn(100.0))        // S6 says safe
ENTRY(p_pnorm,     Rf_pnorm5(1.0,0.0,1.0,1,0))// S6 says safe
ENTRY(p_lbeta,     Rf_lbeta(3.0,4.0))
ENTRY(p_lchoose_k0,Rf_lchoose(1.0,0.0))       // early return, no math library call
ENTRY(p_lchoose,   Rf_lchoose(5.0,2.0))       // lfastchoose -> lbeta
ENTRY(p_choose,    Rf_choose(5.0,2.0))
