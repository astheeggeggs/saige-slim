#ifndef FLEXRV_HPP
#define FLEXRV_HPP

// FlexRV weight kernels (docs/design/flexrv-design.md, step F2).
//
// A FlexRV weight set is one score transform times one MAF transform,
// evaluated per marker after QC on the marker's pathogenicity score s in
// [0, 1], its LoF indicator, and its IN-SAMPLE MAF and MAC. The transform
// table handed to mainRegionInCPP is Kt x 4: (scoreKind, scoreParam, mafKind,
// mafParam), one row per set, score-major and MAF-minor (k = c*nMAF + d),
// built in R by flexrvTransformTable() with the kind codes below. R and C++
// must agree on the codes; the R side names them in R/FlexRV.R.
//
// Every closed form here is the one pinned in the design doc, derived from
// the paper's dbeta ratios: score beta dbeta(1-s,1,b)/dbeta(0,1,b) = s^(b-1);
// MAF beta dbeta(MAF/M,1,b)/dbeta(0,1,b) = (1 - MAF/M)^(b-1), zero at and
// above M (the base would be negative). The cubic root uses std::cbrt, which
// can differ from pow(x, 1/3) by an ulp; the tests compare at 1e-12.

#include <RcppArmadillo.h>

namespace FlexRV {

enum ScoreKind { SCORE_UNW = 0, SCORE_LOF = 1, SCORE_CBRT = 2, SCORE_BETA = 3 };
enum MafKind   { MAF_UNW = 0, MAF_BETA = 1, MAF_MAC = 2 };

double scoreWeight(int t_kind, double t_param, double t_s, bool t_isLof);
double mafWeight(int t_kind, double t_param, double t_maf, double t_mac, double t_maxMAF);

// w = scoreWeight * mafWeight for one table row
double weight(const arma::mat & t_table, unsigned int t_row,
              double t_s, bool t_isLof, double t_maf, double t_mac, double t_maxMAF);

// true when the row is the unweighted x unweighted set (weights identically 1),
// the only transform set the ER exact test may take
bool isUnitRow(const arma::mat & t_table, unsigned int t_row);

// throws on an unknown kind code or a malformed table
void checkTable(const arma::mat & t_table);

}  // namespace FlexRV

#endif
