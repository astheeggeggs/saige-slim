// FlexRV weight kernels: see FlexRV.hpp and docs/design/flexrv-design.md.
#include "FlexRV.hpp"
#include <cmath>
#include <string>

namespace FlexRV {

double scoreWeight(int t_kind, double t_param, double t_s, bool t_isLof){
  switch(t_kind){
    case SCORE_UNW:
      return 1.0;
    case SCORE_LOF:
      return t_isLof ? 1.0 : 0.0;
    case SCORE_CBRT: {
      // w = (sign(s-t)*|s-t|^(1/3) + t^(1/3)) / (|1-t|^(1/3) + t^(1/3))
      const double t = t_param;
      const double d = t_s - t;
      const double sgn = (d > 0.0) - (d < 0.0);
      const double tc = std::cbrt(t);
      return (sgn * std::cbrt(std::fabs(d)) + tc) / (std::cbrt(std::fabs(1.0 - t)) + tc);
    }
    case SCORE_BETA:
      // dbeta(1-|s|, 1, b) / dbeta(0, 1, b) = |s|^(b-1)  (design D1)
      return std::pow(std::fabs(t_s), t_param - 1.0);
    default:
      throw std::invalid_argument("FlexRV::scoreWeight: unknown score transform kind " + std::to_string(t_kind));
  }
}

double mafWeight(int t_kind, double t_param, double t_maf, double t_mac, double t_maxMAF){
  switch(t_kind){
    case MAF_UNW:
      return 1.0;
    case MAF_BETA: {
      // dbeta(MAF/M, 1, b) / dbeta(0, 1, b) = (1 - MAF/M)^(b-1), 0 at and
      // above M (design D2, direction resolved in A1)
      const double x = t_maf / t_maxMAF;
      if(!(x < 1.0)) return 0.0;
      return std::pow(1.0 - x, t_param - 1.0);
    }
    case MAF_MAC:
      return (t_mac <= t_param) ? 1.0 : 0.0;
    default:
      throw std::invalid_argument("FlexRV::mafWeight: unknown MAF transform kind " + std::to_string(t_kind));
  }
}

double weight(const arma::mat & t_table, unsigned int t_row,
              double t_s, bool t_isLof, double t_maf, double t_mac, double t_maxMAF){
  return scoreWeight((int)t_table(t_row, 0), t_table(t_row, 1), t_s, t_isLof)
       * mafWeight((int)t_table(t_row, 2), t_table(t_row, 3), t_maf, t_mac, t_maxMAF);
}

bool isUnitRow(const arma::mat & t_table, unsigned int t_row){
  return (int)t_table(t_row, 0) == SCORE_UNW && (int)t_table(t_row, 2) == MAF_UNW;
}

void checkTable(const arma::mat & t_table){
  if(t_table.n_cols != 4){
    throw std::invalid_argument("FlexRV: the transform table must have 4 columns (scoreKind, scoreParam, mafKind, mafParam); got " + std::to_string(t_table.n_cols));
  }
  for(unsigned int k = 0; k < t_table.n_rows; k++){
    const int sk = (int)t_table(k, 0), mk = (int)t_table(k, 2);
    if(sk < SCORE_UNW || sk > SCORE_BETA){
      throw std::invalid_argument("FlexRV: row " + std::to_string(k + 1) + " has unknown score transform kind " + std::to_string(sk));
    }
    if(mk < MAF_UNW || mk > MAF_MAC){
      throw std::invalid_argument("FlexRV: row " + std::to_string(k + 1) + " has unknown MAF transform kind " + std::to_string(mk));
    }
    if(!std::isfinite(t_table(k, 1)) || !std::isfinite(t_table(k, 3))){
      throw std::invalid_argument("FlexRV: row " + std::to_string(k + 1) + " has a non-finite parameter");
    }
  }
}

}  // namespace FlexRV

// Test/reconciliation export: the q x Kt weight matrix for q markers, exactly
// the values mainRegionInCPP computes per marker (one call per marker there,
// on the same kernels). Not used by the pipeline.
// [[Rcpp::export]]
arma::mat flexrvWeightsInCPP(arma::mat t_table, arma::vec t_score, arma::vec t_isLof,
                             arma::vec t_MAF, arma::vec t_MAC, double t_maxMAF){
  FlexRV::checkTable(t_table);
  const unsigned int q = t_score.n_elem;
  if(t_isLof.n_elem != q || t_MAF.n_elem != q || t_MAC.n_elem != q){
    Rcpp::stop("flexrvWeightsInCPP: score, isLof, MAF and MAC must have one entry per marker");
  }
  if(!(t_maxMAF > 0.0)){
    Rcpp::stop("flexrvWeightsInCPP: maxMAF must be > 0");
  }
  arma::mat W(q, t_table.n_rows);
  for(unsigned int i = 0; i < q; i++){
    for(unsigned int k = 0; k < t_table.n_rows; k++){
      W(i, k) = FlexRV::weight(t_table, k, t_score(i), t_isLof(i) != 0.0, t_MAF(i), t_MAC(i), t_maxMAF);
    }
  }
  return W;
}
