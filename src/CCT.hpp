#ifndef CCT_HPP
#define CCT_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

double CCT_cpp(arma::vec & pval);
// --is_output_log10P: log10 of CCT_cpp's p from log10 p-values (R: CCT_log10)
double CCT_log10_cpp(arma::vec & log10p);

#endif
