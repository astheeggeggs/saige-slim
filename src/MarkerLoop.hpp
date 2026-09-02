#ifndef MARKERLOOP_HPP
#define MARKERLOOP_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <vector>

// S4: the step-2 marker loop's per-marker state, split into what belongs to the
// MARKER (one slot per index, written once, read by the writer) and what belongs
// to the WORKER (reusable buffers, overwritten every marker).
//
// Before this, both lived as ~30 q-length vectors and a handful of hoisted
// buffers in mainMarkerInCPP's frame. The slot vectors were fine except for
// std::vector<bool> isSPAConvergeVec, which is bit-packed: two markers at
// distinct indices in the same word are a data race, not independent writes.
// Giving each marker one MarkerResult makes every per-marker write land in its
// own object, so the race cannot be reintroduced by construction rather than by
// remembering to lock. MarkerScratch is the other half of the same statement:
// those buffers are shared across markers ON PURPOSE (the C3/F8 allocation
// work), so threading needs one instance per worker, and naming them as a unit
// is what makes that possible.
namespace Marker {

// One marker's output row. Defaults reproduce exactly what the q-length vectors
// were constructed with, so a marker that is never tested (QC gate, or a short
// read) leaves the same values behind as before; the writer skips it on
// pval == "NA" either way.
struct MarkerResult {
  std::string chr;
  std::string pos;
  std::string marker;
  std::string ref;
  std::string alt;

  double altCounts = 0.0;
  double altFreq = 0.0;
  double imputationInfo = 0.0;
  double missingRate = 0.0;

  double Beta = arma::datum::nan;
  double seBeta = arma::datum::nan;
  double Tstat = arma::datum::nan;
  double varT = arma::datum::nan;
  std::string pval = "NA";
  std::string pvalNA = "NA";
  bool isSPAConverge = false;

  double Beta_c = arma::datum::nan;
  double seBeta_c = arma::datum::nan;
  double Tstat_c = arma::datum::nan;
  double varT_c = arma::datum::nan;
  std::string pval_c = "NA";
  std::string pvalNA_c = "NA";

  double AF_case = 0.0;
  double AF_ctrl = 0.0;
  uint32_t N_case = 0;
  uint32_t N_ctrl = 0;
  double N_case_hom = 0.0;
  double N_ctrl_het = 0.0;
  double N_case_het = 0.0;
  double N_ctrl_hom = 0.0;
  uint32_t N = 0;

  // Tallied by the caller after the loop instead of incremented inside it: the
  // two Firth counters were the only cross-marker accumulators in the body.
  bool isFirth = false;
  bool isFirthConverge = false;
};

// The reusable buffers. Every one of these is fully overwritten before it is
// read on each marker (GVec and the index vectors by the reader,
// indexZero/NonZeroVec_arma by conv_to, P2Vec and G1tilde_P_G2tilde_Vec inside
// the score tests -- see the comments at their old declaration sites), which is
// why sharing them across markers is value-neutral. It is NOT thread-safe:
// one MarkerScratch per worker.
struct MarkerScratch {
  arma::vec GVec;
  arma::vec gtildeVec;
  arma::vec P2Vec;
  arma::rowvec G1tilde_P_G2tilde_Vec;

  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;
  std::vector<uint> indexForMissing;
  arma::uvec indexZeroVec_arma;
  arma::uvec indexNonZeroVec_arma;

  MarkerScratch(int t_n, int t_numMarkerCond)
    : GVec(t_n), gtildeVec(t_n), G1tilde_P_G2tilde_Vec(t_numMarkerCond)
  {}
};

}  // namespace Marker

#endif
