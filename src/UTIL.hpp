
#ifndef UTIL_HPP
#define UTIL_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <sys/time.h>

const static std::unordered_map<std::string,int> string_to_case{
   {"best_guess",1},
   {"mean",2},
   {"minor",3}
};

double getWeights(std::string t_kernel, 
                  double t_freq, 
                  arma::vec t_wBeta);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);


bool imputeGenoAndFlip(arma::vec& t_GVec,
                       double & t_altFreq,
		       double & t_altCount,
                       std::vector<uint32_t> &  t_indexForMissing,
                       std::string t_impute_method,
                       double t_dosage_zerod_cutoff,
                       double t_dosage_zerod_MAC_cutoff,
                       double & t_MAC,
		       std::vector<uint> & t_indexZero,
                       std::vector<uint> & t_indexNonZero,
                       bool t_buildIndexZero = true);

// The carrier-space twin of imputeGenoAndFlip (2026-09-16): the same
// decisions and the same arithmetic over the k listed samples of a
// Geno::SparseMarker whose common code is 0 or 2, so that altFreq, altCount,
// MAC, the flip and the non-zero index come out bit for bit as the dense call
// on the expanded vector would produce them. On return t_indexNonZero holds
// the analysis indices (increasing) of every sample whose post-impute value
// is non-zero and t_valNonZero those values, parallel.
bool imputeSparseAndFlip(uint32_t t_n,
                         uint32_t t_commonCode,
                         const std::vector<uint32_t> & t_idx,
                         const std::vector<uint8_t> & t_code,
                         double & t_altFreq,
                         double & t_altCount,
                         uint32_t t_nMissing,
                         std::string t_impute_method,
                         double t_dosage_zerod_cutoff,
                         double t_dosage_zerod_MAC_cutoff,
                         double & t_MAC,
                         std::vector<uint> & t_indexNonZero,
                         std::vector<double> & t_valNonZero);

arma::vec getTime();

void printTime(arma::vec t1, arma::vec t2, std::string message);

double getinvStd(double t_freq);

// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
// void set_seed(unsigned int seed) {
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(seed);  
// };

// 
// arma::vec nb(int n){
//   return(Rcpp::rbinom(n,1,0.5));
// }

arma::vec nb(unsigned int n);

double sum_arma1(arma::vec& X);

double add_logp(double p1, double p2);

#endif

