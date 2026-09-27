#define ARMA_USE_SUPERLU 1

// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



#include <sstream>
#include <memory>
#include "SAIGE_test.hpp"
#include "SPA.hpp"
#include "ER_binary_func.hpp"
#include "Binary_global.hpp"   // SKATMainThread::NeedsMainThread (LEDGER #164)
#include "UTIL.hpp"
#include <thread>         // std::this_thread::sleep_for
#include <atomic>
#include <mutex>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>       // getpid, for the dump file names        // getenv: the SAIGE_SPA_PLANT_FAIL test plant (LEDGER #164)
#include <chrono>         // std::chrono::seconds
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/normal.hpp>

// LEDGER #164: the fallback counter (records with cell names are kept in
// Main.cpp, where the names are in scope). Atomic because the region cell
// tests and the marker loop run on workers.
static std::atomic<unsigned long long> g_spaFallbackCount(0);
unsigned long long SAIGE::spaFallbackCount(){ return g_spaFallbackCount.load(); }
void SAIGE::spaFallbackCountReset(){ g_spaFallbackCount.store(0); }
// Test plant: SAIGE_SPA_PLANT_FAIL=1 declares every attempted saddlepoint
// non-converged, so the fallback and its sidecar can be exercised end to end
// on a fixture whose saddlepoints all converge (tests/testthat/
// test-step2-spa-fallback.R reconciles the planted run against a native
// all-ER run). Read once per process; never set it in a real run.
// LEDGER #165: the largest MAC at which the exact test is verified is 200
// (golden P33; it throws on this fixture at 607); --max_MAC_for_ER is capped
// there in R (checkArgs.R). The fallback's own ceiling is the carrier constant
// below, which is stricter.
// LEDGER #164/#165 (2026-09-26, after the first AoU rerun): the FALLBACK runs
// the exact test only in its enumeration regime, i.e. at most 13 carriers
// (2^13 = 8192 <= ExactMax 1e4). Above that the kernel resamples, and its
// draw (SL_Binary_Boot1) converts carriers to cases at rates p_i / mean(p)
// round by round, up to 500 rounds per draw: when the carriers' fitted
// probabilities are heterogeneous the low-probability ones convert at ~1% a
// round and a case-count target above the high-probability carriers costs
// hundreds of rounds x 2e6 draws -- hours on one cell (an AoU chromosome
// sat 40 minutes on a 209-carrier cell). Homogeneous probabilities (the
// sim fixture) finish in seconds, which is why the timing runs missed it.
static const double kFallbackCarrierCeiling = 13.0;
// LEDGER #164 (2026-09-26, found on the first AoU rerun): the ceiling counts
// CARRIERS, not the weighted allele count -- a FlexRV pseudo-marker under
// tiny weights has an AC near 0 and hundreds of carriers. The fallback's
// resampling budget is bounded for the same reason; the native ER branch
// keeps upstream's 2e6 (its bytes are golden).
static const int kFallbackResamples = 200000;
// Diagnostic dump (LEDGER #168 investigation): with SAIGE_SPA_DUMP_DIR set to a
// directory, every saddlepoint FAILURE writes the vector the saddlepoint saw --
// one row per sample: adjusted genotype gtilde, fitted probability mu, phenotype
// y -- under a header line carrying q, the mirrored q, m1, gneg, gpos, Tstat, var1
// and the failure reason, as <dir>/spaFail_<k>.tsv; the sidecar's Note names the
// file. Individual-level data: for use INSIDE the environment that holds it.
static const char* spaDumpDir(){
  static const char* d = std::getenv("SAIGE_SPA_DUMP_DIR");
  return (d != nullptr && d[0] != '\0') ? d : nullptr;
}
static std::atomic<unsigned long> g_spaDumpCounter(0);
static std::string spaDumpVector(const arma::vec & t_gtilde, const arma::vec & t_mu, const arma::vec & t_y,
                                 double t_q, double t_qinv, double t_m1, double t_Tstat, double t_var1, const std::string & t_why){
  const char* dir = spaDumpDir();
  if(dir == nullptr) return "";
  const unsigned long k = g_spaDumpCounter.fetch_add(1);
  // one SAIGE process per chromosome may share the directory (the AoU wrapper
  // does), so the name carries the pid: without it later processes overwrote
  // earlier ones' same-numbered dumps (24 fallbacks, 20 files, 2026-09-27)
  const std::string path = std::string(dir) + "/spaFail_" + std::to_string((long)getpid()) + "_" + std::to_string(k) + ".tsv";
  std::ofstream os(path);
  if(!os) return "";
  double gpos = 0.0, gneg = 0.0;
  for(arma::uword i = 0; i < t_gtilde.n_elem; i++){ if(t_gtilde(i) > 0) gpos += t_gtilde(i); else gneg += t_gtilde(i); }
  os.precision(17);
  os << "# q\t" << t_q << "\tqinv\t" << t_qinv << "\tm1\t" << t_m1 << "\tgneg\t" << gneg << "\tgpos\t" << gpos
     << "\tTstat\t" << t_Tstat << "\tvar1\t" << t_var1 << "\twhy\t" << t_why << "\n";
  os << "gtilde\tmu\ty\n";
  for(arma::uword i = 0; i < t_gtilde.n_elem; i++){ os << t_gtilde(i) << "\t" << t_mu(i) << "\t" << t_y(i) << "\n"; }
  return path;
}
static bool spaPlantFail(){
  static const bool plant = [](){ const char* v = std::getenv("SAIGE_SPA_PLANT_FAIL"); return v != nullptr && std::string(v) == "1"; }();
  return plant;
}
#include <boost/math/distributions/chi_squared.hpp>

// LEDGER #106 fix: the sparse-GRM PCG stops on a RELATIVE residual,
// ||r||^2 <= kSparsePCGRelTol * ||b||^2, so the solve no longer depends on the
// scale of the vector it is handed (on the old absolute 0.02 an unweighted rare
// variant, ||b||^2 ~ MAC, stopped after ONE iteration and the same vector times
// its beta weight after two). Measured against a converged Cholesky at
// N = 1,000 and N = 50,000: 2-3 iterations, variance error <= 1.2e-8, from one
// iteration and errors up to 0.3 %.
static const int kSparsePCGMaxIter = 100;
static const double kSparsePCGRelTol = 1e-8;

namespace SAIGE {

// See the declaration in SAIGE_test.hpp for why this is a named function.
double qnormTwoSidedLogP(double t_logp){
  return R::qnorm(t_logp - std::log(2.0), 0, 1, false, true);
}


SAIGEClass::SAIGEClass(
	arma::mat & t_XVX,
	arma::mat  t_XXVX_inv,
	arma::mat & t_XV,
	arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
	arma::mat & t_X,
	arma::vec &  t_S_a,
	arma::vec & t_res,
	arma::vec & t_mu2,
	arma::vec & t_mu,
	arma::vec & t_varRatio_sparse,
	arma::vec & t_varRatio_null,
	arma::vec & t_varRatio_null_noXadj,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
	double t_SPA_Cutoff,
	arma::vec & t_tauvec,
	std::string t_traitType,
	arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isFastTest,
	bool t_isnoadjCov,
	double t_pval_cutoff_for_fastTest,
	arma::umat & t_locationMat,
	arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
        std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
        double t_pCutoffforFirth,
        std::string t_firthSolver,
        arma::vec & t_offset,
	arma::vec & t_resout){


    m_XVX = t_XVX;
    m_XV = t_XV;
    m_XXVX_inv = t_XXVX_inv;
    m_XVX_inv_XV = t_XVX_inv_XV;
    m_Sigma_iXXSigma_iX = t_Sigma_iXXSigma_iX;
    m_isVarPsadj = false;
    if(m_Sigma_iXXSigma_iX.n_cols == 1 && m_Sigma_iXXSigma_iX.n_rows == 1){
	m_isVarPsadj = false;
    }else{
	m_isVarPsadj = true;
    }
    m_X = t_X;
    m_S_a = t_S_a;
    m_res = t_res;
    m_resout = t_resout;
    m_mu2 = t_mu2;
    // Cached here rather than per marker; same expression, same data.
    m_accu_mu2 = arma::accu(m_mu2);
    m_accu_res = arma::accu(m_res);
    m_mu = t_mu;
    m_varRatio_sparse = t_varRatio_sparse;
    m_varRatio_null = t_varRatio_null;
    m_varRatio_null_noXadj = t_varRatio_null_noXadj;
    m_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
    m_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
    m_tauvec = t_tauvec;
    m_traitType = t_traitType;
    m_y = t_y;

    m_case_indices = arma::find(m_y == 1);
    m_ctrl_indices = arma::find(m_y == 0);

    m_n = t_y.size();
    m_p = t_XV.n_rows;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_impute_method =  t_impute_method;
    m_isCondition = t_isCondition;
    m_condition_genoIndex = t_condition_genoIndex;
    if(m_isCondition){	
	        m_numMarker_cond = t_condition_genoIndex.size();      
    }else{
		m_numMarker_cond = 0;
    }	    
    //m_p = t_X.nrow();

    if(m_traitType == "binary"){
      //if(m_isOutputAFinCaseCtrl){
        m_case_indices = arma::find(m_y == 1);
        m_ctrl_indices = arma::find(m_y == 0);
	m_n_case = m_case_indices.n_elem;
	m_n_ctrl = m_ctrl_indices.n_elem;
	m_is_Firth_beta = t_is_Firth_beta;
	m_pCutoffforFirth = t_pCutoffforFirth;
	m_offset = t_offset;
	if(t_firthSolver == "qr"){
		m_firthSolver = 0;
	}else if(t_firthSolver == "closed"){
		m_firthSolver = 1;
	}else{
		Rcpp::stop("firth_solver must be 'qr' or 'closed', got '" + t_firthSolver + "'.");
	}
      //}
    }
    //m_dimNum = t_dimNum;
    m_flagSparseGRM = t_flagSparseGRM;
    m_isFastTest = t_isFastTest;
    m_isnoadjCov = t_isnoadjCov;
    m_pval_cutoff_for_fastTest = t_pval_cutoff_for_fastTest;
    if(t_dimNum != 0){
	//m_locationMat = t_locationMat;
    	//m_valueVec = t_valueVec;
        m_spSigmaMat = arma::sp_mat(t_locationMat, t_valueVec, t_dimNum, t_dimNum);
	m_diagSigma = arma::vec(m_spSigmaMat.diag());
	// Jacobi preconditioner for getPCG1ofSigmaAndGtilde. m_spSigmaMat and
	// m_diagSigma are written here and nowhere else, so the reciprocal is a
	// constant of the model; it used to be rebuilt (an n-vector allocation and
	// n divisions) on every PCG call, i.e. once per marker on the sparse-GRM
	// marker path and once per marker per region elsewhere.
	m_diagSigmaInv = 1.0 / m_diagSigma;
    }
}    


// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
// E3 (2026-09-14): the printed spelling of every step-2 p-value, written ONCE.
// It used to be pasted at six sites (LEDGER #42 was one defect fixed at three
// of them, #44 at three more), and the log-p branch -- an underflowed p printed
// as "%.1fE%d" from log p with the 9.95 carry -- is reached by no golden, so a
// missed edit at one of six would have moved no golden and failed no test.
// The arithmetic is the old lines verbatim, in the same order.
static std::string formatLogP(double logp){
  char buf[100];
  double log10p = logp/(log(10));
  int exponent = floor(log10p);
  double fraction = pow(10.0, log10p - exponent);
  if (fraction >= 9.95) {
    fraction = 1;
    exponent++;
  }
  sprintf(buf, "%.1fE%d", fraction, exponent);
  return std::string(buf);
}
// The four score-test sites: a non-zero p prints "%.6E"; a p that underflowed
// to 0 is re-derived as log p from the chi-square statistic, printed by
// formatLogP, and handed back IN PLACE as the log with t_islog = true.
static std::string formatPvalOrLog(double & t_pval, double t_stat, bool & t_islog){
  if (t_pval != 0){
    char buf[100];
    sprintf(buf, "%.6E", t_pval);
    t_islog = false;
    return std::string(buf);
  }
  double logp = R::pchisq(t_stat,1,false,true);
  std::string out = formatLogP(logp);
  t_pval = logp;
  t_islog = true;
  return out;
}

void SAIGEClass::set_seed(unsigned int seed){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

// LEDGER #44, second half. `var1 <= DBL_MIN` catches only a variance that has
// underflowed to (near) zero in absolute terms, which is what scoreTest does
// when a marker is exactly explained away. The two "fast" paths build var2 as a
// SUM OF SIGNED TERMS, so the same degeneracy shows up as CANCELLATION instead:
// on the bundled fixture rs20/rs22 give var2 ~ 2.7e-14 against terms of order
// 1e2, which is ~294 orders of magnitude ABOVE DBL_MIN. The old guard could
// never fire there, and BETA = Tstat/var was printed as a ratio of two noise
// quantities -- O(0.1), arbitrary, and not reproducible across platforms.
//
// So ask the question the absolute guard cannot: is var2 at the noise floor OF
// ITS OWN SUMMATION? The relative error of an n-term floating-point sum is
// bounded by about n*eps, so compare |var2| against that fraction of the sum of
// the term magnitudes. This is scale-free -- a rare variant has a small var2 AND
// small terms, so the ratio stays O(1) and it does not fire.
//
// MEASURED on the P30 configuration, all 70 markers: the two covariate-collinear
// markers give 1.7e-16 and 4.2e-16 (machine epsilon, i.e. total cancellation)
// and EVERY other marker gives >= 0.237. That is 14.8 orders of separation, and
// n*eps (2.2e-13 at n = 1000) sits ~500x above the cancellation and ~12 orders
// below the smallest real marker. The threshold is not tuned to the fixture; it
// is the arithmetic's own noise floor, and the fixture confirms the margin.
static inline bool varianceCollapsed(double t_var2, double t_termScale, int t_n)
{
  if(!(t_termScale > 0)) return false;
  return std::fabs(t_var2) <= double(t_n) * std::numeric_limits<double>::epsilon() * t_termScale;
}


// LEDGER #44, THIRD shape, surfaced by the #106 fix. `scoreTest`'s sparse arm
// used to report var == 0 EXACTLY on a covariate-collinear marker, so the
// absolute `var1 <= DBL_MIN` guard fired -- but that zero was an ACCIDENT of
// the old absolute PCG stopping rule: with gtilde ~ 1e-16 componentwise,
// ||b||^2 ~ 1e-29 was already below the old absolute tolerance 1e-8, so the
// solver returned x = 0 without taking a single step. Under the relative rule
// (#106) it converges and returns the true, genuinely tiny variance
// (6.8e-29 for rs20), ~279 orders of magnitude ABOVE DBL_MIN, and BETA =
// S/var1 became a ratio of two noise quantities: -4.4e13 with SE 1.2e14 and
// p = 0.72 on a marker the covariates explain exactly.
//
// `varianceCollapsed` cannot catch this one: there is no cancellation here.
// P2Vec and gtilde are both ~1e-16 componentwise, so var2 is a sum of ~1e-32
// terms that all have the same sign -- it is small because its INPUT is
// nothing, not because the sum destroyed itself. So ask the question one level
// up: did the covariate projection annihilate the genotype? gtilde is
// g - X (X'WX)^-1 X'W g, whose componentwise rounding error is bounded by
// about n*eps times the scale of g, so ||gtilde||^2 <= (n*eps)^2 ||g||^2 says
// "what is left of g after the projection is entirely roundoff". Scale-free by
// construction: it is a ratio of two norms of the same vector's scale.
//
// MEASURED on the P27 configuration, all 70 markers: rs20 and rs22 give
// 9.1e-31 and 7.5e-31; EVERY other marker gives >= 0.328. That is 29.5 orders
// of separation, with the threshold (n*eps)^2 = 4.9e-26 sitting ~5e4 above the
// collapse and ~25 orders below the smallest real marker.
static inline bool projectionCollapsed(const arma::vec & t_gtilde,
                                       const arma::vec & t_GVec, int t_n)
{
  double g2 = arma::dot(t_GVec, t_GVec);
  if(!(g2 > 0)) return false;
  double rel = double(t_n) * std::numeric_limits<double>::epsilon();
  return arma::dot(t_gtilde, t_gtilde) <= rel * rel * g2;
}

void SAIGEClass::scoreTest(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		     double& t_pval,
		     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
		     arma::vec & t_P2Vec,
		     double& t_gy, 
		     bool t_is_region,
		     arma::uvec & t_indexForNonZero,
		     const SAIGE::TestContext & t_ctx) const {
    arma::vec Sm, var2m;
    double S, var2;
    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y);
     }
    S = dot(t_gtilde, m_res);
    //std::cout << "S " << S << std::endl;
    S = S/m_tauvec[0];

    //std::cout << "S b " << S << std::endl;

    if(!t_ctx.useSparseGRM){
      t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0];  
      var2m = dot(t_P2Vec , t_gtilde);
       //std::cout << "!t_ctx.useSparseGRM" << std::endl;
    }else{
      //arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat();
      //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
      //double var2m_test = dot(t_P2Vec , t_gtilde);
      //std::cout << "var2m_test " << var2m_test << std::endl;
      t_P2Vec = getPCG1ofSigmaAndGtilde(t_gtilde, kSparsePCGMaxIter, kSparsePCGRelTol); 
      var2m = dot(t_P2Vec , t_gtilde);
      //std::cout << "var2m " << var2m << std::endl;
      if(m_isVarPsadj){
	var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;	
      }
             //std::cout << "m_flagSparseGRM_cur" << std::endl;
    }	      
    var2 = var2m(0,0);
    double var1 = var2 * t_ctx.varRatioVal;

    // LEDGER #44: the absolute guard catches the exactly-zero variance; the
    // projection test catches the marker whose covariate-adjusted genotype is
    // pure roundoff, whose variance is tiny but nowhere near DBL_MIN.
    const bool collapsed = (var1 <= std::numeric_limits<double>::min()) ||
                           projectionCollapsed(t_gtilde, t_GVec, m_n);
    double stat = S*S/var1;
    if (collapsed){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));
      }else{
          t_pval = 1;
	  stat = 0.0;
      }	      
    }
    t_pval_str = formatPvalOrLog(t_pval, stat, t_islogp);
    // LEDGER #44 FIXED. These two used to be computed unconditionally, even
    // under the guard above that already forces t_pval = 1 because the
    // covariate-adjusted variance has collapsed. At var1 == 0 exactly that gave
    // BETA = -inf, stat = inf and SE = inf/inf = nan; at var1 ~ 1e-14 a finite
    // but meaningless BETA with an SE of order 1e6. The p-value was right in
    // both cases -- only the effect sizes were junk -- so the fix is to say so
    // rather than to print a number nobody should read.
    //
    // NA_REAL, not a bare NaN, and the writer prints it as the literal "NA":
    // the spelling of a NaN is PLATFORM-DEPENDENT (glibc writes a sign-bit-set
    // NaN as "-nan" where macOS writes "nan"), which is what made these rows
    // non-byte-portable and cost two red CI runs. A row that says NA says the
    // same thing everywhere.
    if(collapsed){
      t_Beta = NA_REAL;
      t_seBeta = NA_REAL;
    }else{
      t_Beta = S/var1;
      t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    }
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     const SAIGE::TestContext & t_ctx) const {

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    arma::vec Z = A1.t() * g1;
    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    // LEDGER #152: var2, Bmu2 and g1tildemu2 are assigned only inside the
    // trait-type chain below and read after it. The chain is exhaustive on
    // this branch (binary, quantitative) and now says so with a terminal
    // else; the NaN initialisers make any future gap a NaN, not a stale
    // stack value.
    double var1, S, S1, S2;
    double var2 = std::numeric_limits<double>::quiet_NaN();
    double g1tildemu2 = std::numeric_limits<double>::quiet_NaN();
    arma::vec S_a2;
    double Bmu2 = std::numeric_limits<double>::quiet_NaN();
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    if(m_traitType == "binary"){
      mu21  = m_mu2.elem(t_indexForNonZero);
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }else{
      Rcpp::stop("scoreTest: unknown trait type '" + m_traitType + "' (LEDGER #152)");
    }

    // Term magnitudes for the cancellation test above. Free: all three are
    // already computed. g1tildemu2 is only assigned on the binary branch, so it
    // is folded in only there.
    double termScale = std::fabs(ZtXVXZ(0,0)) + std::fabs(Bmu2)
                     + (m_traitType == "binary" ? std::fabs(g1tildemu2) : std::fabs(dot(g1, g1)));
    bool collapsed = varianceCollapsed(var2, termScale, m_n);

    var1 = var2 * t_ctx.varRatioVal;
    S1 = dot(res1, g1_tilde);
    arma::mat res1X1_temp = (res1.t()) * X1;
    arma::vec res1X1 = res1X1_temp.t();
    S_a2 = m_S_a - res1X1;
    S2 = - arma::dot(S_a2,  Z);
    S = S1 + S2;
    S = S/m_tauvec[0];

    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){	    
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));

      }else{
          t_pval = 1;
	  stat = 0.0;
      }
    }
    t_pval_str = formatPvalOrLog(t_pval, stat, t_islogp);
    // LEDGER #44 FIXED -- rationale at scoreTest() above (E3: one copy, not three).
    if(var1 <= std::numeric_limits<double>::min() || collapsed){
      t_Beta = NA_REAL;
      t_seBeta = NA_REAL;
    }else{
      t_Beta = S/var1;
      t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    }
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast_noadjCov(arma::vec & t_GVec,
		     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     bool t_flip,
                     const SAIGE::TestContext & t_ctx) const {

    //arma::vec Sm, var2m;
      arma::vec g1 = t_GVec.elem(t_indexForNonZero);
      arma::vec m_mu21 = m_mu2.elem(t_indexForNonZero);
      arma::vec m_res1 = m_res.elem(t_indexForNonZero);
    // LEDGER #28: t_GVec is minor-coded (imputeGenoAndFlip mirrors G := 2-G when
    // altFreq > 0.5) but t_altFreq is restored to ALT orientation for reporting,
    // so centring on 2*t_altFreq is centring on the wrong allele whenever ALT is
    // major. Take the frequency of the allele t_GVec actually counts. When
    // t_flip is false this is t_altFreq itself -- no arithmetic, so unflipped
    // markers are bit-identical to the pre-fix code.
    double centreFreq = t_flip ? 1.0 - t_altFreq : t_altFreq;
    double S, var2;
    //getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);
    //getadjG(t_GVec, t_gtilde);


    //if(t_is_region && m_traitType == "binary"){
    //  t_gy = dot(t_gtilde, m_y);
    // }
     
    double var2_a = dot(m_mu21,pow(g1,2));
    double var2_b = dot(m_mu21, 2*2*centreFreq*g1);
    double var2_c = m_accu_mu2*pow(2*centreFreq, 2);
    var2 = var2_a - var2_b + var2_c;
    bool collapsed = varianceCollapsed(var2, std::fabs(var2_a) + std::fabs(var2_b)
                                             + std::fabs(var2_c), m_n);
    var2 = var2 *(m_tauvec[0]);
    //std::cout << "var2 " << var2 << std::endl;

    //arma::vec t_GVec_center = t_GVec-arma::mean(t_GVec);
    //var2 = dot(m_mu2, pow(t_GVec_center,2))*(m_tauvec[0]);  
    //std::cout << "var2 b " << var2 << std::endl;
    //S = dot(t_GVec_center, m_res);
    //std::cout << "S " << S << std::endl;
    S = dot(g1, m_res1)  - m_accu_res*(2*centreFreq);
    //std::cout << "S b " << S << std::endl; 
    S = S/m_tauvec[0];


    //var2 = var2m(0,0);
    double var1 = var2 * t_ctx.varRatioVal;
    //std::cout << "var1 " << var1 << std::endl;
    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));
      }else{
          t_pval = 1;
          stat = 0.0;
      }
    }
    t_pval_str = formatPvalOrLog(t_pval, stat, t_islogp);
    // LEDGER #44 FIXED -- rationale at scoreTest() above (E3: one copy, not three).
    if(var1 <= std::numeric_limits<double>::min() || collapsed){
      t_Beta = NA_REAL;
      t_seBeta = NA_REAL;
    }else{
      t_Beta = S/var1;
      t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    }
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}





void SAIGEClass::getadjG(arma::vec & t_GVec, arma::vec & g){
   g = m_XV * t_GVec;
   //t_GVec.t().print("t_GVec");
   //std::cout << "sum(t_GVec) " << arma::accu(t_GVec) << std::endl;
   //      m_XV.print("m_XV");
    g = t_GVec - m_XXVX_inv * g;
    //m_XXVX_inv.print("m_XXVX_inv");
    //g.t().print("g");
}


void SAIGEClass::getadjGFast(arma::vec & t_GVec, arma::vec & g, arma::uvec & iIndex) const
{

  // To increase computational efficiency when lots of GVec elements are 0
 arma::vec m_XVG(m_p, arma::fill::zeros);
  for(int i = 0; i < iIndex.n_elem; i++){
      m_XVG += m_XV.col(iIndex(i)) * t_GVec(iIndex(i));
  }
  g = t_GVec - m_XXVX_inv * m_XVG; 
}


void SAIGEClass::get_mu(arma::vec & t_mu){
    t_mu = m_mu;
}

void SAIGEClass::getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices){
     t_case_indices = m_case_indices;
     t_ctrl_indices = m_ctrl_indices;
  }

/*
void SAIGEClass::setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR) {
    m_locationMat = locationMatinR;
    m_valueVec = valueVecinR;
    m_dimNum = r;
}



arma::sp_mat SAIGEClass::gen_sp_SigmaMat() {
    arma::sp_mat resultMat(m_locationMat, m_valueVec, m_dimNum, m_dimNum);
    return resultMat;
}
*/

// revised
// need to add sparse Sigma version 
// This function only uses variance ratio and does not use sparse GRM
void SAIGEClass::getMarkerPval(arma::vec & t_GVec,
			       arma::uvec & iIndex,
			       arma::uvec & iIndexComVec,
                               double& t_Beta, 
                               double& t_seBeta,
			       std::string& t_pval,
			       std::string& t_pval_noSPA,	
                               double t_altFreq,
                               double& t_Tstat,
			       double& t_gy,
			       double& t_var1, 
			       bool & t_isSPAConverge, 
			       arma::vec & t_gtilde, 
			       bool & is_gtilde,
			       bool  is_region,
                               arma::vec & t_P2Vec,
			       bool t_isCondition,
			       double& t_Beta_c,
                           	double& t_seBeta_c,
			       std::string& t_pval_c,
                               std::string& t_pval_noSPA_c,	
                           	double& t_Tstat_c,
                           	double& t_varT_c,
			   	arma::rowvec & t_G1tilde_P_G2tilde, 
				bool & t_isFirth,
				bool & t_isFirthConverge, 
				bool t_isER,
				bool t_flip,
				const SAIGE::TestContext & t_ctx,
				KernelWorkspace * t_ws) const
{



  t_isFirth = false;
  // The kernel workspace: the caller's if it brought one (the marker path, one
  // per worker), else one local instance built the first time a kernel needs it.
  std::unique_ptr<KernelWorkspace> localWs;
  auto ws = [&]() -> KernelWorkspace & {
    if(t_ws != nullptr) return *t_ws;
    if(!localWs) localWs.reset(new KernelWorkspace());
    return *localWs;
  };
  ws().spaFallback = false;   // LEDGER #164: cleared at every entry
  ws().spaFallbackNote.clear();
  //arma::vec adjGVec = getadjGFast(t_GVec);
  std::string t_pval_str;
  double t_var2, t_SPApval;
  //iIndex = arma::find(t_GVec != 0);
  //arma::vec t_gtilde;
  bool isScoreFast = true;
  if(t_ctx.useSparseGRM){
    isScoreFast = false;
  }

 
  double pval_noadj, pval; //can be log or not raw
  double t_qval_Firth = std::numeric_limits<double>::quiet_NaN(); // LEDGER #152: assigned wherever t_isFirth is set, read under if(t_isFirth)
  bool ispvallog;

  //for test
  //arma::vec timeoutput3 = getTime();
  //if(!isScoreFast){
	//std::cout << "scoreTest " << std::endl;  
  //	is_gtilde = true;
 // 	scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex);
 // }else{
if(!t_ctx.useSparseGRM && t_ctx.noadjCov){
	is_gtilde = false;
        isScoreFast = true;
	scoreTestFast_noadjCov(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq,t_Tstat, t_var1, t_var2, t_flip, t_ctx);

}else if(t_ctx.useSparseGRM){
	is_gtilde = true;
        isScoreFast = false;
        scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex, t_ctx);
	//std::cout << "useSparseGRM " << t_ctx.useSparseGRM << std::endl;
}else{
	is_gtilde = false;
        isScoreFast = true;
        //std::cout << "scoreTestFast "  << std::endl;
        scoreTestFast(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_ctx);	
}


  double StdStat = std::abs(t_Tstat) / sqrt(t_var1);

  t_isSPAConverge = false;
  
 //arma::vec timeoutput3_a = getTime();
  // LEDGER #152: q/qinv are assigned on the binary branch and the SPA block
  // runs only for a non-quantitative trait; NAmu/NAsigma are assigned under
  // p_iIndexComVecSize >= 0.5, the predicate that also selects the fast
  // kernel that reads them. Flow-guarded, not compiler-provable; NaN so a
  // gap can never be a stale stack value.
  double m1, tol1, p_iIndexComVecSize;
  double q = std::numeric_limits<double>::quiet_NaN(), qinv = std::numeric_limits<double>::quiet_NaN(), NAmu = std::numeric_limits<double>::quiet_NaN(), NAsigma = std::numeric_limits<double>::quiet_NaN();


  //arma::uvec iIndexComVec = arma::find(t_GVec == 0);
  //arma::uvec iIndexVec = arma::find(t_GVec != 0);
 
  
  // 2026-09-16: the SIZE of the zero-genotype set is m_n minus the non-zero
  // count on every caller (the two index vectors partition the samples: the
  // marker path via imputeGenoAndFlip / imputeSparseAndFlip, the region path
  // via complementIndices), so the marker path stopped materialising the
  // vector and passes an empty one. This is that identity, not a new count.
  unsigned int iIndexComVecSize = (unsigned int)(m_n - (int)iIndex.n_elem);
  unsigned int iIndexSize = iIndex.n_elem;
  // Sized lazily. These two are only read on the SPA branches below, which
  // together take ~5% of binary markers and no quantitative one; sizing them
  // for every marker cost allocations of ~n doubles each, i.e. mmap/munmap
  // round trips per marker at biobank n, ~95% of them wasted.
  //
  // The marginal branch assigns them wholesale (gNB = t_gtilde(iIndex) etc.),
  // so it needs no explicit sizing. The CONDITIONAL branch is a different
  // predicate and can be reached without the marginal branch having run
  // (LEDGER #40); it computes its own saddlepoint context (see the block
  // around SPA_ws below).
  //
  // gNA/muNA -- the zero-genotype halves -- were materialised here too and
  // threaded through SPA_fast_ws and the two fast Binom kernels, which never
  // read them (inherited; upstream's bodies do not read them either). Deleted
  // 2026-09-14 (dead-surface sweep C4) from the three INTERNAL signatures
  // only; the Rcpp-exported wrappers keep theirs.
  //
  // LEDGER #40 history: 8f3ded00 fixed it with a prepSPAContext() lambda
  // armed at the top of the conditional block; b60141f9 (#84) made that block
  // compute its own CGF and removed the call, and the lambda plus its
  // 12-line "this is the fix" comment stayed here with no caller until the
  // 2026-09-14 sweep. The guard is the conditional block's own arithmetic, not
  // a lambda; goldens P12/P17 pin it.
  arma::vec gNB, muNB;
  double gmuNB;

//std::cout << "ok1" << std::endl;


if((StdStat > m_SPA_Cutoff || std::isnan(StdStat)) && m_traitType != "quantitative" && t_isER){
	t_isER = true;
}else{	
	t_isER = false;
}

//std::cout << "ok2" << std::endl;

  // LEDGER #164: the exact test as one callable, so the saddlepoint's
  // failure branch below can run it in place of the normal p. Body unchanged.
  auto runExactTest = [&](int t_nResampling){

    arma::mat Z_er(t_GVec.n_elem, 1);
    Z_er.col(0) = t_GVec;
    // LEDGER #167: the exact test ranks configurations by Q = (sum Z (y - pi))^2
    // and calls two configurations tied when |dQ| <= 1e-6 ABSOLUTE
    // (ComputeExact::Run, m_epsilon). Under weights of 1e-6 every
    // configuration is a tie and p is ~1. The test is scale-invariant in exact
    // arithmetic, so a vector whose largest entry is below 1 is scaled up to 1
    // first; a 0/1/2 genotype vector (largest entry 1 or 2) is untouched, which
    // keeps every existing exact-test row byte-identical.
    {
      double zmax = 0.0;
      for(arma::uword e = 0; e < iIndex.n_elem; e++){ zmax = std::max(zmax, std::fabs(t_GVec(iIndex(e)))); }
      if(zmax > 0.0 && zmax < 1.0){ Z_er.col(0) /= zmax; }
    }
    arma::vec res_er = m_res;
    arma::vec pi1_er = m_mu;
    arma::vec resout_er = m_resout;
    /*std::cout << "t_GVec(0) " << t_GVec(0) << std::endl;
    std::cout << "res_er(0) " << res_er(0) << std::endl;
    std::cout << "pi1_er(0) " << pi1_er(0) << std::endl;
    resout_er.print("resout_er");
    std::cout << "iIndex.n_elem" << iIndex.n_elem << std::endl;
    std::cout << "iIndexComVec.n_elem" << iIndexComVec.n_elem << std::endl;
    std::cout << "m_n_case " << m_n_case << std::endl;
        iIndex.print("iIndex");
*/	
    // 2026-09-16: the one consumer of the zero-genotype INDICES. The marker
    // path hands an empty vector (see iIndexComVecSize above), so this branch
    // -- main-thread only and reached for MAC <= max_MAC_for_ER markers --
    // builds the complement itself, in the increasing order every caller
    // built it in. A caller that still materialises it is used as is.
    arma::uvec iIndexComLocal;
    arma::uvec * iIndexComPtr = &iIndexComVec;
    if(iIndexComVec.n_elem == 0 && iIndex.n_elem < (arma::uword)m_n){
      iIndexComLocal.set_size((arma::uword)m_n - iIndex.n_elem);
      arma::uword* out = iIndexComLocal.memptr();
      const arma::uword* nz = iIndex.memptr();
      std::size_t kk = 0, o = 0;
      for(arma::uword i = 0; i < (arma::uword)m_n; i++){
        if(kk < iIndex.n_elem && nz[kk] == i){ kk++; }else{ out[o++] = i; }
      }
      iIndexComPtr = &iIndexComLocal;
    }
    double pval_ER =  SKATExactBin_Work(Z_er, res_er, pi1_er, m_n_case, iIndex, *iIndexComPtr, resout_er, t_nResampling, 1e+4, 1e-6, 1);
    char pValueBuf_ER[100];
    sprintf(pValueBuf_ER, "%.6E", pval_ER);
    std::string buffAsStdStr_ER = pValueBuf_ER;
    t_pval = pValueBuf_ER;
    
    //std::cout << "t_pval " << t_pval << std::endl;
    pval = pval_ER;
    boost::math::normal ns;
    double t_qval_ER;
    try{
      t_qval_ER = boost::math::quantile(ns, pval_ER/2);
      t_qval_ER = fabs(t_qval_ER);
      t_seBeta = fabs(t_Beta)/t_qval_ER;
      t_isSPAConverge = true;
    }catch (const std::overflow_error&) {
      t_qval_ER = std::numeric_limits<double>::infinity();
      t_seBeta = 0;
    }

    // `&` (not `&&`) and the parentheses are deliberate: <= binds tighter than
    // & so this is what the expression always parsed as; the parentheses only
    // make it explicit for -Wparentheses (Phase D, 2026-09-14). Do not change
    // the operator: both operands are evaluated today and that is the
    // preserved behaviour.
    if(m_is_Firth_beta & (pval <= m_pCutoffforFirth)){
	t_isFirth = true;
	t_qval_Firth = t_qval_ER;
    }
  };

bool spaFellBack = false;   // LEDGER #164
if(!t_isER){


  if(!std::isnan(StdStat) && (StdStat > m_SPA_Cutoff) && m_traitType != "quantitative"){

       if(!is_gtilde){
          t_gtilde.resize(m_n);
          getadjGFast(t_GVec, t_gtilde, iIndex);
	  is_gtilde = true;
       }
	//int t_gtilden = t_gtilde.n_elem;
        p_iIndexComVecSize = double(iIndexComVecSize)/m_n;
   	m1 = dot(m_mu, t_gtilde);

	if(p_iIndexComVecSize >= 0.5){
	gNB = t_gtilde(iIndex);
   	muNB = m_mu(iIndex);
  	gmuNB = dot(gNB,muNB);	 
   	NAmu= m1-gmuNB;

   }

   	if(m_traitType == "binary"){
                q = t_Tstat/sqrt(t_var1/t_var2) + m1;

                if((q-m1) > 0){
                        qinv = -1 * std::abs(q-m1) + m1;
                }else if ((q-m1) == 0){
                        qinv =  m1;
                }else{
                        qinv = std::abs(q-m1) + m1;
                }
		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % (1-muNB) % arma::pow(gNB,2));
		}
  }
    	//bool logp=false;
	double tol0 = std::numeric_limits<double>::epsilon();
	tol1 = std::pow(tol0, 0.25);
	if(p_iIndexComVecSize >= 0.5 && !t_ctx.useSparseGRM){
		//std::cout << "SPA_fast" << std::endl;
        	SPA_fast_ws(t_gtilde, q, qinv, pval_noadj, ispvallog, gNB, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge, ws().spa);
	
	}else{
		//std::cout << "SPA" << std::endl;
		SPA_ws(m_mu, t_gtilde, q, qinv, pval_noadj, tol1, ispvallog, m_traitType, t_SPApval, t_isSPAConverge, ws().spa);	
	}

    boost::math::normal ns;
    double t_qval;
    // LEDGER #64 FIXED: the score-test SE, kept so the rollback below can
    // restore it. The arm overwrites t_seBeta from the saddlepoint p and then
    // may discover that p is unusable (exactly 0 on the linear scale, where
    // R::qnorm returns Inf rather than throwing, so |Beta|/Inf = 0); the
    // p-value was rolled back to the score p and the SE was left at 0.
    const double seBeta_score = t_seBeta;

    if(t_isSPAConverge){
        try {
          // LEDGER #42 FIXED. The two-sided quantile wants p/2, and on the log
          // scale that is log(p) - log(2). Dividing the LOG by 2 is a square
          // root: it evaluated the quantile at sqrt(p) instead of p/2, which in
          // this regime inflates SE by ~sqrt(2) (measured 1.4188x at the very
          // edge, log p = -710, tending to 1.41421x deeper in). The !ispvallog
          // arm is spelled out rather than left polymorphic so the two cannot
          // drift again; it is the same call it always was.
          if(!ispvallog){
            t_qval = R::qnorm(t_SPApval/2, 0, 1, false, false);
          }else{
            t_qval = qnormTwoSidedLogP(t_SPApval);
          }
          t_qval = fabs(t_qval);
          t_seBeta = fabs(t_Beta)/t_qval;
        }catch (const std::overflow_error&) {
          t_qval = std::numeric_limits<double>::infinity();
	  t_isSPAConverge = false;
        }
    }


    if(!ispvallog && t_SPApval == 0){
          t_isSPAConverge = false;
          t_seBeta = seBeta_score; // LEDGER #64: the p goes back to the score p; so does the SE
    }


    if(spaPlantFail()){ t_isSPAConverge = false; }   // test plant only
    if(!t_isSPAConverge){
      // LEDGER #164 FIXED: the saddlepoint was attempted and failed (LEDGER
      // #163's root overflow, an unattainable mirrored q, a NaN step, the
      // iteration cap, or a p that underflowed). The normal-approximation p
      // is what the saddlepoint exists to replace, so it is NOT reported:
      // the exact test runs instead, at any MAC and for any weight set. It
      // may throw NeedsMainThread on a worker; the caller's deferral re-runs
      // this whole call on the main thread, where the flag is set again.
      // why the saddlepoint failed, for the sidecar
      std::string why;
      if(spaPlantFail()){ why = "planted"; }
      else if(!ispvallog && t_SPApval == 0){ why = "p underflowed to 0"; }
      else switch(ws().spa.failReason){
        case 1: why = "root(q) not found"; break;
        case 2: why = "root(mirrored q) not found"; break;
        case 3: why = "saddle(q) invalid"; break;
        case 4: why = "saddle(mirrored q) invalid"; break;
        default: why = "qnorm overflow or unknown"; break;
      }
      {
        const SPAWorkspace & sw = ws().spa;
        char dbuf[256];
        snprintf(dbuf, sizeof(dbuf), " [root1=%.4g/%d root2=%.4g/%d qsaddle=%s bound=%.3g]",
                 sw.diagRoot1, sw.diagNiter1, sw.diagRoot2, sw.diagNiter2,
                 sw.diagQSaddle == 1 ? "valid" : (sw.diagQSaddle == 0 ? "INVALID" : "n/a"), sw.diagBound);
        why += dbuf;
      }
      ws().spaFallbackNote = "spa: " + why;
      if(spaDumpDir() != nullptr){
        if(!is_gtilde){ getadjGFast(t_GVec, t_gtilde, iIndex); is_gtilde = true; }
        const std::string f = spaDumpVector(t_gtilde, m_mu, m_y, q, qinv, m1, t_Tstat, t_var1, why);
        if(!f.empty()) ws().spaFallbackNote += "; dump=" + f;
      }
      // carriers, not the (weighted) allele count: see kFallbackResamples
      const double carriers = (double)iIndex.n_elem;
      try {
        if(carriers > kFallbackCarrierCeiling){
          throw std::runtime_error(std::to_string((long long)carriers) +
                                   " carriers, above the fallback's exact-test ceiling " + std::to_string((long long)kFallbackCarrierCeiling) +
                                   " (enumeration only; LEDGER #164/#165); no usable p");
        }
        runExactTest(kFallbackResamples);
      } catch (const SKATMainThread::NeedsMainThread&) {
        throw;   // the caller's deferral re-runs this call on the main thread
      } catch (const Rcpp::LongjumpException&) {
        throw;   // an R-level unwind must reach Rcpp's handler; never swallow it
      } catch (const std::exception& e) {
        // The exact test itself failed (LEDGER #165: the native ER path throws
        // on this fixture above MAC ~600). Neither number is usable: the
        // normal p is what the saddlepoint exists to replace, and there is no
        // exact p. The p is NA, the row says why, the run goes on.
        t_pval = "NA"; pval = std::numeric_limits<double>::quiet_NaN();
        t_isSPAConverge = false;
        ws().spaFallbackNote += std::string("; exact test not run or failed: ") + e.what();
      } catch (...) {
        t_pval = "NA"; pval = std::numeric_limits<double>::quiet_NaN();
        t_isSPAConverge = false;
        ws().spaFallbackNote += "; exact test failed: non-standard exception (LEDGER #165)";
      }
      spaFellBack = true;
      ws().spaFallback = true;
      g_spaFallbackCount.fetch_add(1);
    }
  }
   //t_pval_noSPA = pval_noadj; 
   char pValueBuf_SPA[100];
   if(spaFellBack){
        // t_pval, pval, t_seBeta, t_isSPAConverge and the Firth decision were
        // set by runExactTest; nothing to format here.
   }else if(m_traitType!="quantitative"){
        if(t_isSPAConverge){
		if(!ispvallog){
		  sprintf(pValueBuf_SPA, "%.6E", t_SPApval);
		}else{
		  snprintf(pValueBuf_SPA, sizeof(pValueBuf_SPA), "%s", formatLogP(t_SPApval).c_str());
		}
		std::string buffAsStdStr_SPA = pValueBuf_SPA;
	        t_pval = buffAsStdStr_SPA;
		pval = t_SPApval;	
        }else{
                t_pval = t_pval_noSPA;
		pval = pval_noadj;
		//pval_log10 = pval_noadj_log10;
        }
     if(m_traitType=="binary"){
	if(!ispvallog){
		if(m_is_Firth_beta && pval <= m_pCutoffforFirth){
			t_isFirth = true;
			t_qval_Firth = R::qnorm(pval/2, 0, 1, false, false);
		}

	}else{
		if(m_is_Firth_beta && pval <= std::log(m_pCutoffforFirth)){
			t_isFirth = true;
			// LEDGER #42: log(p) - log(2), not log(p)/2 -- see the SE
			// back-calculation above. Same defect, same regime; this one
			// feeds the Firth SE instead of the score one.
			t_qval_Firth = qnormTwoSidedLogP(pval);
		}
	}
     }	

   }else{
        t_pval = t_pval_noSPA;
	pval = pval_noadj;
	//pval_log10 = pval_noadj_log10;
   }

}else{ //if(!t_isER){
    runExactTest(2000000);
}
   if(t_isFirth){
	if(!is_gtilde){
                getadjGFast(t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
	arma::vec init(2, arma::fill::zeros);
	if(m_firthSolver == 1){
		fast_logistf_fit_closed(t_gtilde, m_y, m_offset, true, init, 50, 15, 1e-5, 1e-5, t_Beta ,t_seBeta, t_isFirthConverge, ws());
	}else{
		fast_logistf_fit_simple(t_gtilde, m_y, m_offset, true, init, 50, 15, 15, 1e-5, 1e-5, 1e-5, t_Beta ,t_seBeta, t_isFirthConverge, ws());
	}
	//back calculates se based on beta from firth adjustion and the p-value that accounts for case-control imbalance
	t_seBeta = fabs(t_Beta)/fabs(t_qval_Firth);
   }
   
 //arma::vec timeoutput4 = getTime();
 //printTime(timeoutput3, timeoutput3_a, "Test Marker  ScoreTest");
//printTime(timeoutput3, timeoutput4, "Test Marker 3 to 4");
//printTime(timeoutput3_a, timeoutput4, "Test Marker SPA");

   //condition
   if(t_isCondition){
	if(!is_gtilde){
        	getadjGFast(t_GVec, t_gtilde, iIndex);
        	is_gtilde = true;
        }
        t_G1tilde_P_G2tilde = sqrt(t_ctx.varRatioVal) * t_gtilde.t() * m_P2Mat_cond;
        arma::vec t_Tstat_ctemp =  t_G1tilde_P_G2tilde * m_VarInvMat_cond * m_Tstat_cond;
	arma::mat tempgP2 = t_gtilde.t() * m_P2Mat_cond;

    	t_Tstat_c = t_Tstat - t_Tstat_ctemp(0);
    	arma::vec t_varT_ctemp = t_G1tilde_P_G2tilde * m_VarInvMat_cond * (t_G1tilde_P_G2tilde.t());
    	t_varT_c = t_var1 - t_varT_ctemp(0);

    // LEDGER #66 FIXED. The conditional variance is a DIFFERENCE, var1 minus
    // the part of it the conditioning markers explain, so a tested marker that
    // is collinear with the conditioning set cancels it to float noise. That is
    // the same geometry LEDGER #44 found in the two "fast" score tests, and it
    // defeats an absolute guard for the same reason: on the constructed high-LD
    // fixture a marker byte-identical to the conditioning marker gives
    // varT_c ~ 2.5e-14 against var1 ~ 28, which is ~294 orders of magnitude
    // ABOVE DBL_MIN. The guard below could never fire there, and the row was
    // published as p.value_c = 1 with an SE_c of order 1e6 and a BETA_c that is
    // the ratio of two cancelled quantities -- arbitrary, and not reproducible
    // across platforms.
    //
    // So ask #44's question, with #44's threshold and #44's helper: is varT_c at
    // the noise floor of its own two-term difference? The terms are var1 and
    // varT_ctemp(0), both of them sums over the n samples, so n*eps is the right
    // relative floor and nothing here is tuned to a fixture. MEASURED on the
    // constructed high-LD fixture (7 markers, n = 1000, so n*eps = 2.22e-13):
    // the markers collinear with the conditioning set give 4.89e-16 and
    // 7.19e-16 -- machine epsilon, i.e. total cancellation -- while the
    // r^2 = 0.9 partner gives 5.02e-02 and every other marker more. That is
    // 14 orders of separation, with the threshold ~450x above the cancellation
    // and ~11 orders below the smallest real value. Under the CLI default
    // (--is_noadjCov=TRUE) the smallest ratio on that same fixture is 4.1e-03,
    // so that path is untouched -- see LEDGER #82 for why the two differ.
    bool collapsed_c = varianceCollapsed(t_varT_c,
                                         std::fabs(t_var1) + std::fabs(t_varT_ctemp(0)),
                                         m_n);

    double S_c = t_Tstat_c;

    double stat_c = S_c*S_c/t_varT_c;

    double pval_noSPA_c;

     if (t_varT_c <= std::numeric_limits<double>::min() || collapsed_c){
        pval_noSPA_c = 1;
	stat_c = 0;
     }else{
       if(!std::isnan(stat_c) && std::isfinite(stat_c)){	     
        boost::math::chi_squared chisq_dist(1);
        pval_noSPA_c = boost::math::cdf(complement(chisq_dist, stat_c));
       }else{
        pval_noSPA_c = 1;
	stat_c = 0;
       }	       
     }

    t_pval_noSPA_c = formatPvalOrLog(pval_noSPA_c, stat_c, ispvallog);
    // LEDGER #66 FIXED, second half: these two used to be computed
    // unconditionally, PAST the guard above that has already forced
    // pval_noSPA_c = 1 and stat_c = 0 -- so the guard-fires row was published
    // with BETA_c = S_c/0 and SE_c = |BETA_c|/sqrt(0), which is exactly the
    // inf/nan pair upstream reports. NA_REAL rather than a bare NaN, for the
    // reason recorded at the marginal sites: a NaN's spelling is
    // platform-dependent ("nan" on macOS, "-nan" on glibc when the sign bit is
    // set) and the writers print an R_IsNA as the literal "NA" everywhere.
    if(t_varT_c <= std::numeric_limits<double>::min() || collapsed_c){
      t_Beta_c = NA_REAL;
      t_seBeta_c = NA_REAL;
    }else{
      t_Beta_c = S_c/t_varT_c;
      t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
    }
    t_Tstat_c = S_c;


    bool t_isSPAConverge_c;
    if(m_traitType != "quantitative" && stat_c > std::pow(m_SPA_Cutoff,2)){
        // LEDGER #84 FIXED (2026-09-05, owner decision: residualise). The
        // conditional score T_c = T - c'T2 with c = (g~'P G2~)(G2~'P G2~)^-1 is
        // exactly the score of the RESIDUALISED genotype g~_c = g~ - G2~ c, so
        // its saddlepoint must use g~_c's cumulant generating function -- its
        // own mean m1_c = sum(mu g~_c) and naive variance var2_c =
        // sum(mu(1-mu) g~_c^2) -- and not the marginal marker's, which is what
        // this block used to hand SPA_fast: anti-conservative, and growing with
        // the LD between the tested and the conditioning markers. g~_c is dense
        // (G2~ c touches every sample), so SPA_fast's zero-genotype shortcut
        // does not apply and the full SPA is used; the block runs only when the
        // conditional statistic clears the cutoff, so the extra n x k product
        // costs nothing measurable. Under --is_noadjCov=TRUE the marginal T is
        // the UNADJUSTED score while the conditioning terms are adjusted
        // (LEDGER #82), so there q_c is the mixture statistic saddlepointed
        // against g~_c's CGF; under FALSE the identity T_c = g~_c'res is exact.
        if(!is_gtilde){
                getadjGFast(t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
        arma::vec coef_c = (t_G1tilde_P_G2tilde * m_VarInvMat_cond).t();
        arma::vec gtilde_c = t_gtilde - m_G2tilde_cond * coef_c;
        double m1_c = arma::dot(m_mu, gtilde_c);
        double var2_c = arma::sum(m_mu % (1 - m_mu) % arma::square(gtilde_c));
        double tol1_c = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
	double SPApval_c;
	double q_c = std::numeric_limits<double>::quiet_NaN(), qinv_c = std::numeric_limits<double>::quiet_NaN(); // LEDGER #152
	if(m_traitType == "binary"){
                q_c = t_Tstat_c/sqrt(t_varT_c/var2_c) + m1_c;
                if((q_c-m1_c) > 0){
                        qinv_c = -1 * std::abs(q_c-m1_c) + m1_c;
                }else if ((q_c-m1_c) == 0){
                        qinv_c =  m1_c;
                }else{
                        qinv_c = std::abs(q_c-m1_c) + m1_c;
                }
        }else{
                // unreachable: the enclosing block is guarded by m_traitType != "quantitative"
                Rcpp::stop("conditional SPA: unknown trait type '" + m_traitType + "' (LEDGER #152)");
        }
        SPA_ws(m_mu, gtilde_c, q_c, qinv_c, pval_noSPA_c, tol1_c, ispvallog, m_traitType, SPApval_c, t_isSPAConverge_c, ws().spa);

	
        boost::math::normal ns;
        //double pval_SPA = t_SPApval;
        double t_qval_c;
     if(t_isSPAConverge_c){	
        try {
          if(!ispvallog){
           t_qval_c = boost::math::quantile(ns, SPApval_c/2);
          }else{
           // LEDGER #42: log(p) - log(2). The !ispvallog sibling one line up
           // takes the two-sided quantile at p/2 correctly, so before this fix
           // the same statistic was computed two different ways depending only
           // on whether the score p had underflowed.
	   t_qval_c = qnormTwoSidedLogP(SPApval_c);
          }
           t_qval_c = fabs(t_qval_c);
           t_seBeta_c = fabs(t_Beta_c)/t_qval_c;
        }catch (const std::overflow_error&) {
          t_qval_c = std::numeric_limits<double>::infinity();
          //t_seBeta_c = 0;
	  t_isSPAConverge_c = false;
        }
      }

	if(!ispvallog && SPApval_c == 0){
		t_isSPAConverge_c = false;	
	}


      char pValueBuf_SPA_c[100];
        if(t_isSPAConverge_c){
                if(!ispvallog){
                  //LEDGER #41: print the saddlepoint p itself. The /2
                  //belongs to the two-sided quantile at :929 that builds
                  //t_seBeta_c, not to the p-value; printing it here made
                  //every conditional SPA p-value 2x anti-conservative.
                  sprintf(pValueBuf_SPA_c, "%.6E", SPApval_c);
                }else{
                  snprintf(pValueBuf_SPA_c, sizeof(pValueBuf_SPA_c), "%s", formatLogP(SPApval_c).c_str());
                }
                std::string buffAsStdStr_SPA_c = pValueBuf_SPA_c;
                t_pval_c = buffAsStdStr_SPA_c;
        }else{
                t_pval_c = t_pval_noSPA_c;
        }

    }else{
	t_isSPAConverge_c = false;    
    	t_pval_c = t_pval_noSPA_c;	    
    }	    
 }


    gNB.clear();
    muNB.clear();   // was a second gNB.clear() -- the copy-paste typo the audit found (C4)


//std::cout << "ok3" << std::endl;

    if(is_region && !is_gtilde){
	getadjGFast(t_GVec, t_gtilde, iIndex);
	is_gtilde = true; 
    }

    if(is_region && isScoreFast){

      t_gy = dot(t_gtilde, m_y);
      if(!t_ctx.useSparseGRM){
        t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0];
      }else{
        //arma::sp_mat m_SigmaMat_sp = gen_sp_SigmaMat();
        //t_P2Vec = arma::spsolve(m_SigmaMat_sp, t_gtilde);
	t_P2Vec = getPCG1ofSigmaAndGtilde(t_gtilde, kSparsePCGMaxIter, kSparsePCGRelTol);
      }
    }
}


// The bin search that the assign* family performs, with no member write, so the
// marker path can own its result instead of parking it on the object between
// markers (S3). The arithmetic and the comparison order are the originals.
//
// LEDGER #32 is fixed here rather than at the call sites: when a MAC falls in no
// bin the old code returned false WITHOUT assigning, and every caller discarded
// that return value, so the marker silently reused the previous marker's ratio.
// Fixing it at the landmine means no caller can re-arm it. Coverage is
// (-inf, minExc[0]], each (minExc[i], maxInc[i]], and (maxInc.back(), inf), so a
// gap needs a specification whose bins leave a hole; the CLI validates only the
// lengths of the two vectors, which is what makes it reachable.
double SAIGEClass::selectVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj) const {
    arma::vec varRatio;
    if(issparseforVR){
        varRatio = m_varRatio_sparse;
    }else{
        if(!isnoXadj){
            varRatio = m_varRatio_null;
        }else{
            varRatio = m_varRatio_null_noXadj;
        }
    }

    double varRatioVal = 0.0;
    bool t_hasVarRatio = false;   // local since E1; see the header
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
                varRatioVal = varRatio(i);
                t_hasVarRatio = true;
        }
    }

    if(!t_hasVarRatio){
        if(MAC <= m_cateVarRatioMinMACVecExclude(0)){
                varRatioVal = varRatio(0);
                t_hasVarRatio = true;
        }
    }

    if(!t_hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
                varRatioVal = varRatio.back();
                t_hasVarRatio = true;
        }
    }

    if(!t_hasVarRatio){
        std::ostringstream msg;
        msg << "no categorical variance ratio covers MAC " << MAC
            << "; the bins from --cateVarRatioMinMACVecExclude / "
            << "--cateVarRatioMaxMACVecInclude leave a gap there. Bins (exclusive lower, "
            << "inclusive upper]:";
        for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++){
            msg << " (" << m_cateVarRatioMinMACVecExclude(i) << ", "
                << m_cateVarRatioMaxMACVecInclude(i) << "]";
        }
        throw std::runtime_error(msg.str());
    }

    return varRatioVal;
}


// Pure form of the two-argument assignSingleVarianceRatio. The arma::rowvec is
// the original's, deliberately: assigning a multi-element arma::vec to a rowvec
// throws "Mat::init(): requested size is not compatible with row vector layout"
// (the LEDGER #10 message), so this is safe only for a single-ratio model, which
// is exactly when its callers use it. Preserved rather than quietly widened.
SAIGE::TestContext SAIGEClass::contextForMAC(double MAC) const {
    SAIGE::TestContext ctx;
    ctx.useSparseGRM = (MAC > m_cateVarRatioMinMACVecExclude.back()) ? false : m_flagSparseGRM;
    ctx.noadjCov = false;
    if(isSingleVarianceRatio()){
        ctx.varRatioVal = selectSingleVarianceRatio(ctx.useSparseGRM, false);
    }else{
        ctx.varRatioVal = selectVarianceRatio(MAC, ctx.useSparseGRM, false);
    }
    return ctx;
}

double SAIGEClass::selectSingleVarianceRatio(bool issparseforVR, bool isnoXadj) const {
    arma::rowvec varRatio;
    if(issparseforVR){
        varRatio = m_varRatio_sparse;
    }else{
        if(isnoXadj){
            varRatio = m_varRatio_null_noXadj;
        }else{
            varRatio = m_varRatio_null;
        }
    }
    return varRatio(0);
}


// S12: currentContext() and the whole assign*/set_* family are gone with the
// three members they staged -- m_varRatioVal, m_flagSparseGRM_cur and
// m_isnoadjCov_cur. Both step-2 paths now build a TestContext where the choice
// is made and hand it to the test, so the model is read-only for the duration
// of a scan. The pure selectors above are what is left.




void SAIGEClass::assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_G2tilde_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      std::vector<std::string> & t_p_cond
      ){
	m_P2Mat_cond = t_P2Mat_cond;
	m_G2tilde_cond = t_G2tilde_cond;
	m_VarInvMat_cond = t_VarInvMat_cond;
	m_VarMat_cond = t_VarMat_cond;
	m_Tstat_cond = t_Tstat_cond;
	m_MAF_cond = t_MAF_cond;
	m_qsum_cond = t_qsum_cond;
	m_gsum_cond = t_gsum_cond;
	m_G2_Weight_cond = t_G2_Weight_cond;
	m_p_cond = t_p_cond;
}

void SAIGEClass::assignConditionFactors_scalefactor(
	arma::vec & t_scalefactor_G2_cond	
		){
	m_scalefactor_G2_cond = t_scalefactor_G2_cond;
	arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
	arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t(); 
	arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * m_VarMat_cond * scalefactor_G2_cond_Mat;
	arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
	m_VarInvMat_cond_scaled_weighted = arma::pinv(VarMat_cond_scaled_weighted);
	//m_VarInvMat_cond_region_binary = (1/scalefactor_G2_cond_Mat) * m_VarInvMat_cond	* (1/scalefactor_G2_cond_Mat);
	
}

void SAIGEClass::extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv) const {
	t_XV = m_XV;
	t_XXVX_inv = m_XXVX_inv;	
}



void SAIGEClass::fast_logistf_fit_simple(const arma::vec & t_g,
                const arma::vec & y,
                const arma::vec & offset,
                bool firth,
        arma::vec init,
        int maxit,
        int maxstep,
        int maxhs,
        double lconv,
        double gconv,
        double xconv,
        double & beta_G,
        double & sebeta_G, 
	bool & isfirthconverge,
	KernelWorkspace & w) const {
  // Workspace edition (2026-09-05). The original built, per Newton iteration,
  // wpi, W2, XW2, Q, R, Q%Q, h, ypih, xcol, XX_XW2 and XX_Fisher afresh -- ten
  // n-length or n x 2 allocations, each above libmalloc's large threshold at
  // biobank n -- plus the n x 2 design matrix per call. Every quantity below is
  // the same expression on the same operands, evaluated in the same order:
  // XX_XW2 was an element-for-element copy of XW2, so the Fisher matrix is the
  // same syrk over the same values; Q % Q * oneVec is, for k = 2, one addition
  // of two squares per row, which is what sum(square(Q), 1) computes; and
  // sqrt(pi % (1 - pi)) evaluates the same two elementwise operations the wpi
  // temporary did. Measured bit-identical on all 7,891 Firth markers of the
  // N = 50,000 binary scan and 23% faster per fit in isolation.
  isfirthconverge = false;	
  const arma::uword n = t_g.n_elem;
  const int k = 2;
  w.ensureFirth(n);
  arma::mat & x = w.x;
  x.col(0).ones();
  x.col(1) = t_g;
  arma::vec beta = init;
  int iter = 0;
  w.pi_0 = -x * beta - offset;
  w.pi_0 = arma::exp(w.pi_0) + 1;
  w.pi = 1/w.pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat XX_covs(k, k, arma::fill::zeros);
  arma::mat XX_Fisher;
  arma::vec U_star(2, arma::fill::zeros);
  arma::vec delta;
  while(iter <= maxit){
        beta_old = beta;
        w.W2 = arma::sqrt(w.pi % (1 - w.pi));
        for(int j = 0; j < k; j++){
                w.XW2.col(j) = x.col(j) % w.W2;
        }
        arma::qr_econ(w.Q, w.R, w.XW2);
        w.h = arma::sum(arma::square(w.Q), 1);
        if(firth){
                w.ypih = (y - w.pi) + (w.h % (0.5 - w.pi));
        }else{
                w.ypih = (y - w.pi);
        }
        U_star = x.t() * w.ypih;
        XX_Fisher = w.XW2.t() * w.XW2;
        bool isinv = arma::inv_sympd (XX_covs, XX_Fisher);
	if(!isinv){
                break;
        }
        delta = XX_covs * U_star;
        double mx = arma::max(arma::abs(delta))/maxstep;
        if(mx > 1){
                delta = delta/mx;
        }
        evals = evals + 1;
        iter = iter + 1;
        beta = beta + delta;
        w.pi_0 = -x * beta - offset;
        w.pi_0 = arma::exp(w.pi_0) + 1;
        w.pi = 1/w.pi_0;
        // LEDGER #26 FIXED (2026-09-05, owner decision: flag it). Reaching
        // maxit used to be reported as convergence, so the "N sucessfully
        // converged" line counted the markers that had drifted at maxstep
        // for every iteration (|BETA| = maxstep*maxit = 750 on the bundled
        // fixture's null markers). The flag is now true only when the step
        // and gradient criteria are met; at the cap the loop exits with it
        // false. BETA/SE are unchanged either way.
        if((arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv))){
		isfirthconverge = true;
                break;
        }
        if(iter == maxit){
                break;
        }
  }
        if(XX_covs.has_nan()){
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = sqrt(XX_covs(1,1));
        }
}

void SAIGEClass::fast_logistf_fit_closed(const arma::vec & t_g,
                const arma::vec & y,
                const arma::vec & offset,
                bool firth,
        arma::vec init,
        int maxit,
        int maxstep,
        double gconv,
        double xconv,
        double & beta_G,
        double & sebeta_G,
	bool & isfirthconverge,
	KernelWorkspace & w) const {
  // --firth_solver=closed. With k = 2 the Fisher matrix X'WX is three weighted
  // sums and the hat values are h_i = w_i (c00 + 2 c01 g_i + c11 g_i^2) with
  // C = (X'WX)^-1, so the QR of the n x 2 weighted design (12% of the binary
  // scan's CPU in the 2026-09-05 profile) and the design matrix itself are not
  // needed. Same penalised score, same step control, same stopping rule as
  // fast_logistf_fit_simple; the arithmetic differs in rounding only (measured
  // <= 8e-13 relative on BETA and SE at n = 50,000).
  isfirthconverge = false;
  const arma::uword n = t_g.n_elem;
  w.ensureFirth(n);
  arma::vec & wpi = w.W2;
  arma::vec beta = init;
  int iter = 0;
  w.pi_0 = -(beta(0) + beta(1) * t_g) - offset;
  w.pi_0 = arma::exp(w.pi_0) + 1;
  w.pi = 1/w.pi_0;
  arma::mat XX_covs(2, 2, arma::fill::zeros);
  arma::mat F(2, 2, arma::fill::zeros);
  arma::vec U_star(2, arma::fill::zeros);
  arma::vec delta;
  while(iter <= maxit){
        wpi = w.pi % (1 - w.pi);
        double f00 = arma::accu(wpi);
        double f01 = arma::dot(wpi, t_g);
        w.h = wpi % t_g;
        double f11 = arma::dot(w.h, t_g);
        F(0,0) = f00; F(0,1) = f01; F(1,0) = f01; F(1,1) = f11;
        bool isinv = arma::inv_sympd(XX_covs, F);
        if(!isinv){
                break;
        }
        if(firth){
                w.h = wpi % (XX_covs(0,0) + 2*XX_covs(0,1)*t_g + XX_covs(1,1)*arma::square(t_g));
                w.ypih = (y - w.pi) + (w.h % (0.5 - w.pi));
        }else{
                w.ypih = (y - w.pi);
        }
        U_star(0) = arma::accu(w.ypih);
        U_star(1) = arma::dot(t_g, w.ypih);
        delta = XX_covs * U_star;
        double mx = arma::max(arma::abs(delta))/maxstep;
        if(mx > 1){
                delta = delta/mx;
        }
        iter = iter + 1;
        beta = beta + delta;
        w.pi_0 = -(beta(0) + beta(1) * t_g) - offset;
        w.pi_0 = arma::exp(w.pi_0) + 1;
        w.pi = 1/w.pi_0;
        if((arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv))){
		isfirthconverge = true;
                break;
        }
        if(iter == maxit){
                break;
        }
  }
        if(XX_covs.has_nan()){
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = sqrt(XX_covs(1,1));
        }
}

arma::vec SAIGEClass::getPCG1ofSigmaAndGtilde(arma::vec& bVec, int maxiterPCG, double tolPCG) const {
    int Nnomissing = m_spSigmaMat.n_rows;
    arma::vec xVec(Nnomissing, arma::fill::zeros); // Initialize xVec to zeros
    arma::vec rVec = bVec; // Residual vector
    arma::vec zVec(Nnomissing);
    // Preconditioner: m_diagSigmaInv is 1.0 / m_diagSigma, built once with the
    // model instead of once per call (see the constructor).
    zVec = m_diagSigmaInv % rVec; // Apply preconditioner
    double sumr2 = arma::dot(rVec, rVec); // Initial residual norm, ||b||^2
    const double stopAt = tolPCG * sumr2;  // relative: tolPCG is a fraction of ||b||^2
    arma::vec pVec = zVec; // Search direction

    // Hoisted out of the loop: each of these was a fresh n-double allocation on
    // every iteration. Armadillo's operator= reuses the buffer when the
    // dimensions already match, so they are allocated once per call and the
    // expressions -- and therefore the arithmetic -- are untouched.
    arma::vec ApVec, r1Vec, z1Vec;

    int iter = 0;
    while (sumr2 > stopAt && iter < maxiterPCG) {
        iter++;
        ApVec = m_spSigmaMat * pVec; // Sparse matrix-vector multiplication
        double alpha = arma::dot(rVec, zVec) / arma::dot(pVec, ApVec); // Step size
        xVec = xVec + alpha * pVec;

        r1Vec = rVec - alpha * ApVec; // Update residual
        z1Vec = m_diagSigmaInv % r1Vec; // Apply preconditioner to new residual
        double beta = arma::dot(z1Vec, r1Vec) / arma::dot(zVec, rVec); // Update beta

        pVec = z1Vec + beta * pVec; // Update search direction
        // Was `zVec = z1Vec; rVec = r1Vec;` -- two n-double copies per
        // iteration. The old contents of z1Vec/r1Vec are dead (both are fully
        // reassigned at the top of the next iteration), so exchanging the
        // buffers leaves every value bit-identical and copies nothing.
        zVec.swap(z1Vec); // Update preconditioned residual
        rVec.swap(r1Vec); // Update residual
        sumr2 = arma::dot(rVec, rVec); // Update residual norm
    }

    if (iter >= maxiterPCG) {
        std::cout << "PCG in getPCG1ofSigmaAndGtilde did not converge. You may increase maxiter number." << std::endl;
    }

    //std::cout << "Iterations from getPCG1ofSigmaAndGtilde: " << iter << std::endl;
    return xVec;
}



}


// Test/reconciliation export (LEDGER #166/#167): the exact resampling test on
// one vector, exactly as getMarkerPval's runExactTest calls it (m = 1 column,
// test_type 1, ExactMax 1e4, epsilon 1e-6), plus the weight-independent
// configuration probabilities P(j cases among the carriers) it uses. The R
// reference in tests/testthat/test-unit-er.R takes probK from here and does
// the weight-dependent enumeration itself. Not used by the pipeline.
// [[Rcpp::export]]
Rcpp::List erExactTestInCPP(arma::vec t_Z, arma::vec t_y, arma::vec t_mu, int t_nResampling, bool t_normalise){
  if(t_Z.n_elem != t_y.n_elem || t_Z.n_elem != t_mu.n_elem){ Rcpp::stop("erExactTestInCPP: Z, y, mu must have one entry per sample"); }
  arma::uvec idx = arma::find(t_Z != 0.0);
  arma::uvec idxComp = arma::find(t_Z == 0.0);
  arma::vec res = t_y - t_mu;
  arma::mat Z(t_Z.n_elem, 1); Z.col(0) = t_Z;
  if(t_normalise){
    double zmax = 0.0;
    for(arma::uword e = 0; e < idx.n_elem; e++){ zmax = std::max(zmax, std::fabs(t_Z(idx(e)))); }
    if(zmax > 0.0 && zmax < 1.0){ Z.col(0) /= zmax; }
  }
  arma::mat resout;   // empty: one residual vector, as in getMarkerPval
  uint32_t ncase = (uint32_t)arma::accu(t_y);
  double p = SKATExactBin_Work(Z, res, t_mu, ncase, idx, idxComp, resout, t_nResampling, 1e+4, 1e-6, 1);
  const unsigned int k = idx.n_elem;
  std::vector<int> n_total_k(k + 1, 0), IsExactVec(k + 1, 1);
  std::vector<double> prob(k + 1, 0.0);
  int n_total = 0; bool isExact = false;
  SKATExactBin_ComputProb_New(idx, idxComp, t_mu, (uint32_t)t_Z.n_elem, ncase, t_nResampling, 1e+4, 1, 2, prob, IsExactVec, n_total_k, n_total, isExact);
  return Rcpp::List::create(Rcpp::Named("pval") = p, Rcpp::Named("probK") = prob,
                            Rcpp::Named("carriers") = (int)k, Rcpp::Named("exact") = isExact);
}
