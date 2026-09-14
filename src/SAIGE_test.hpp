
#ifndef SAIGE_HPP
#define SAIGE_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "SPA_binary.hpp"

namespace SAIGE{

// Per-worker scratch for the two n-length kernels the binary marker path runs
// after the score test: the saddlepoint kernels (SPAWorkspace, SPA_binary.hpp)
// and the Firth penalised fit, which used to allocate ~10 n-length temporaries
// per Newton iteration (9 iterations at 1% prevalence, 7,891 fits on the
// N = 50,000 binary scan, 56% of its CPU). Sized lazily on first use, so a scan
// that reaches neither (every quantitative scan) allocates nothing. Owned by
// Marker::MarkerScratch on the marker path, one per worker (S4); a caller that
// has none passes nullptr to getMarkerPval and gets one local instance per
// call, which is still one allocation set per marker rather than per iteration.
struct KernelWorkspace {
        SPAWorkspace spa;
        arma::mat x, XW2, Q, R;
        arma::vec pi_0, pi, W2, h, ypih;
        void ensureFirth(arma::uword n) {
                if(x.n_rows != n){
                        x.set_size(n, 2); XW2.set_size(n, 2);
                        pi_0.set_size(n); pi.set_size(n); W2.set_size(n); h.set_size(n); ypih.set_size(n);
                }
        }
};


// Everything the per-marker test needs that used to live in mutable members of
// SAIGEClass (S3 in docs/design/step2-threading-verdict.md). Both loops wrote
// three members before every call -- m_varRatioVal, m_flagSparseGRM_cur,
// m_isnoadjCov_cur -- and getMarkerPval read them back, which is a data race the
// moment two markers are tested at once, and (via m_varRatioVal) a genuine
// cross-marker read even when they are not: see LEDGER #32.
//
// S3 gave the marker path a context; S12 gave the region path one and DELETED
// the three members, so this struct is now the only place a per-marker choice
// is recorded. A test writes nothing to the model, on either path.
//
// One thing did not move here: the fastTest ESCALATION signal that R raises
// between region calls. It shared m_flagSparseGRM_cur with the per-marker
// sparse-GRM choice, which is what LEDGER #27 and #50 were, and it is now a
// separate step-2 variable in Main.cpp (g_regionEscalatedPass). If you find
// yourself wanting to put it in this struct, that is the collision coming back.
struct TestContext {
        double varRatioVal = 0.0;   // the variance ratio selected for this marker
        bool useSparseGRM = false;  // was m_flagSparseGRM_cur
        bool noadjCov = false;      // was the t_isnoadjCov argument / m_isnoadjCov_cur
};

class SAIGEClass
{
    private:
      arma::mat m_XVX;
      arma::mat m_XVX_inv_XV;
      arma::mat m_X;
      arma::mat m_Sigma_iXXSigma_iX;
      arma::vec m_res;
      arma::vec m_resout;
      arma::vec m_mu;
      arma::vec m_mu2;
      // accu(m_mu2) and accu(m_res): full-n reductions over model constants that
      // scoreTestFast_noadjCov (the default marker path) was recomputing for
      // every marker. Both vectors are set once and never mutated, so caching
      // the identical arma::accu call is bit-identical.
      double m_accu_mu2;
      double m_accu_res;
      arma::vec m_tauvec;
      arma::vec  m_S_a;
      //std::string m_traitType; 
      std::string m_impute_method;
      std::vector<uint32_t> m_condition_genoIndex;
      	

    public:
           std::string m_traitType;
      arma::mat m_XXVX_inv;
      arma::mat m_XV;
      int m_n, m_p; //MAIN Dimensions: sample size, number of covariates
      arma::vec m_varRatio_sparse;
      arma::vec m_varRatio_null;
      arma::vec m_varRatio_null_noXadj;
      arma::vec m_y;

      bool m_isOutputAFinCaseCtrl;
      bool m_isOutputNinCaseCtrl;
      bool m_isOutputHetHomCountsinCaseCtrl;
      arma::uvec m_case_indices;
      arma::uvec m_ctrl_indices;
      arma::uvec m_case_hom_indices;
      arma::uvec m_case_het_indices;
      arma::uvec m_ctrl_hom_indices;
      arma::uvec m_ctrl_het_indices;
      int m_n_case;
      int m_n_ctrl;
      arma::sp_mat m_SigmaMat_sp;
      bool m_flagSparseGRM;
      bool m_isFastTest;
      bool m_isnoadjCov;
      double m_pval_cutoff_for_fastTest; 
      double m_SPA_Cutoff;
      //arma::umat m_locationMat;
      //arma::vec m_valueVec;
      //int m_dimNum;	
      arma::vec m_cateVarRatioMinMACVecExclude; 
      arma::vec m_cateVarRatioMaxMACVecInclude;
      arma::mat m_P2Mat_cond;
      // LEDGER #84: the covariate-adjusted genotypes of the conditioning
      // markers, n x k, one column each, so the conditional saddlepoint can
      // form the residualised genotype g~ - G2~ c.
      arma::mat m_G2tilde_cond;
      int m_numMarker_cond;
      arma::mat m_VarInvMat_cond;
      arma::mat m_VarMat_cond;
      arma::vec m_Tstat_cond;
      arma::vec m_G2_Weight_cond;
      arma::vec m_MAF_cond;
      double  m_qsum_cond;
      arma::vec m_gsum_cond;
      std::vector<std::string> m_p_cond;
      arma::vec m_scalefactor_G2_cond;
      arma::mat m_VarInvMat_cond_scaled_weighted;
      //arma::mat m_VarInvMat_cond_region_binary;
      bool m_isCondition;
      bool m_is_Firth_beta;
      double m_pCutoffforFirth;
      // --firth_solver: 1 = "closed" (the default since 2026-09-06: 2x2 closed-form
      // Fisher matrix and hat values, no QR), 0 = "qr" (upstream's QR-based hat
      // values). Same estimator; the two differ in rounding only (<= 8e-13
      // relative at n = 50,000, never reaching the printed digits there).
      // A warm start at the null intercept was benchmarked and REJECTED: 118 of
      // 7,891 fits stopped converging and BETA moved up to 4x (NOTES.md).
      int m_firthSolver = 1;
     arma::vec  m_offset;	
      bool m_isVarPsadj;
      bool m_islog10p;
	arma::sp_mat m_spSigmaMat; // Declare the sparse matrix attribute
	arma::vec m_diagSigma;    // Precompute diagonal
	arma::vec m_diagSigmaInv; // 1/diag: the PCG's Jacobi preconditioner, built once


  ////////////////////// -------------------- functions ----------------
  //------------------ //////////////////////
  

  SAIGEClass(
        arma::mat & t_XVX,
        arma::mat  t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
	arma::mat & t_sigmainvX_XsignmainXtXinv,
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
	arma::vec & t_resout);

   void set_seed(unsigned int seed);

   void scoreTest(arma::vec & t_GVec,
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
		     const TestContext & t_ctx) const;


void scoreTestFast(arma::vec & t_GVec,
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
                     const TestContext & t_ctx) const;

// t_flip: true when imputeGenoAndFlip() mirrored this marker (ALT was the major
// allele), in which case t_GVec is minor-coded while t_altFreq has been restored
// to ALT orientation -- the centre must be the minor frequency. See LEDGER #28.
void scoreTestFast_noadjCov(arma::vec & t_GVec,
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
                     const TestContext & t_ctx) const;




     void get_mu(arma::vec & t_mu);

     void getadjG(arma::vec & t_GVec, arma::vec & g);
     void getadjGFast(arma::vec & t_GVec, arma::vec & g,  arma::uvec & iIndex) const;

     void getMarkerPval(arma::vec & t_GVec,
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
				const TestContext & t_ctx,
				KernelWorkspace * t_ws = nullptr) const;

    void getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices);


    //void setupSparseMat(int r, arma::umat & locationMatinR, arma::vec & valueVecinR);

    //arma::sp_mat gen_sp_SigmaMat();

    // Pure selectors: the bodies of the assign* family with no member write.
    // A MAC that falls in no bin throws rather than leaving the caller with
    // the previous marker's ratio (LEDGER #32). The former `bool &
    // t_hasVarRatio` out-parameter was written at every call site and read
    // by none of them -- only by the three guards INSIDE this function -- and
    // a frame-shared copy of it was the worker-vs-worker race TSAN caught on
    // 2026-09-06 (d36e7277 racing, e72e4186 clean); demoted to a local on
    // 2026-09-14 (audit E1), so that class cannot recur by construction.
    double selectVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj) const;
    double selectSingleVarianceRatio(bool issparseforVR, bool isnoXadj) const;
    // Single vs categorical, spelled ONCE (it was `== 1` at two sites and
    // `<= 1` at a third, which differ only at n_elem == 0 -- a `null`-less
    // ratio file, which Get_Variance_Ratio now refuses).
    bool isSingleVarianceRatio() const { return m_varRatio_null.n_elem <= 1; }
    // ONE derivation of "how is this marker / cell tested" for the four sites
    // that share it (the marker path's fastTest re-test, the region marker
    // prepare, the ultra-rare cell prepare and the BURDEN cell prepare; audit
    // E1, 2026-09-14): a common marker (MAC above the last category bound)
    // uses the non-sparse variance ratio, a rare one follows the model; the
    // covariate adjustment is on; the ratio is the single one or the MAC's
    // band. Drift between hand-written copies of exactly this rule is the
    // documented cause of LEDGER #118 -> #134 -> #135. Two sites deliberately
    // keep their own: testOneModel's first pass (fastTest forces the
    // non-sparse pass and noadjCov follows the model) and the
    // conditioning-marker setup (no MAC term at all -- routing it here would
    // ADD the gate and move golden P4's _c columns).
    SAIGE::TestContext contextForMAC(double MAC) const;

    void assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_G2tilde_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      std::vector<std::string> & t_p_cond);

     void assignConditionFactors_scalefactor(
        arma::vec & t_scalefactor_G2_cond);	


    void extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv) const;   // step 4: read by every model of a region call

    // The Firth fit over [1, t_g] with the workspace's buffers; t_ws.ensureFirth
    // is called inside. Same operations in the same order as the original,
    // proved bit-identical on the N = 50,000 binary scan (2026-09-05).
    void fast_logistf_fit_simple(const arma::vec & t_g,
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
	KernelWorkspace & t_ws) const;	

    // --firth_solver=closed: the same penalised Newton iteration with the 2x2
    // Fisher matrix and the hat values in closed form (h_i = w_i x_i' F^-1 x_i),
    // no design matrix, no QR. Not bit-identical to the QR path.
    void fast_logistf_fit_closed(const arma::vec & t_g,
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
	KernelWorkspace & t_ws) const;


     arma::vec getPCG1ofSigmaAndGtilde(arma::vec& bVec, int maxiterPCG, double tolPCG) const;

};

// Upper-tail two-sided normal quantile for a p-value given on the LOG scale --
// i.e. the |z| whose two-sided p is exp(t_logp).
//
// It exists as a named function because the same expression is needed at three
// places in getMarkerPval and all three had the SAME defect (LEDGER #42): they
// wrote qnorm(logp / 2, log.p = TRUE), and HALVING A LOG IS A SQUARE ROOT, so
// the quantile was evaluated at sqrt(p) rather than p/2 -- inflating SE by
// ~sqrt(2) wherever the score p had underflowed and the code switched to the
// log format. One function is one place to get it wrong.
//
// Exported for test through testQnormTwoSidedLogPInCPP (src/Main.cpp): the
// regime that reaches it needs |z| ~ 40, which no bundled fixture produces, so
// no golden can gate this and the change carries its own measurement.
double qnormTwoSidedLogP(double t_logp);
}
#endif
