
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "UTIL.hpp"
#include <sys/time.h>

arma::vec nb(unsigned int n){
  return(Rcpp::rbinom(n,1,0.5));
}



double getWeights(std::string t_kernel, 
                  double t_freq, 
                  arma::vec t_wBeta)
{
  if(t_wBeta.size() != 2)
    Rcpp::stop("The size of argument t_wBeta should be 2.");
  
  double weights;
  if(t_kernel == "linear")
    weights = 1;
  
  if(t_kernel == "linear.weighted"){
    Rcpp::NumericVector freq = {t_freq};
    Rcpp::NumericVector temp = Rcpp::dbeta(freq, t_wBeta(0), t_wBeta(1));
    weights = temp(0);
  }
  
  return weights;
}

// imputeGeno (a second imputation rule with the vocabulary "bestguess"/"none",
// declared with 3 parameters and defined with 4, so neither half matched the
// other) was deleted 2026-09-14 (dead-surface sweep C5). imputeGenoAndFlip
// below is the live rule; anchors that cited it moved up by the deleted lines.
// used in Main.cpp::mainMarkerInCPP
bool imputeGenoAndFlip(arma::vec& t_GVec, 
                       double & t_altFreq,
		       double & t_altCount, 
                       std::vector<uint32_t> & t_indexForMissing,
                       std::string t_impute_method,
		       double t_dosage_zerod_cutoff,
		       double t_dosage_zerod_MAC_cutoff, 
		       double & t_MAC,
		       std::vector<uint> & t_indexZero,
		       std::vector<uint> & t_indexNonZero,
		       bool t_buildIndexZero)
{
  bool flip = false;
  t_indexNonZero.clear();
  t_indexZero.clear();
  int nMissing = t_indexForMissing.size();
  uint dosagesSize = t_GVec.size();
  double imputeG = 0;
  if(t_altFreq > 0.5){
    flip = true;
    t_GVec = 2 - t_GVec;
    t_altFreq = 1 - t_altFreq; 
  }

if(nMissing > 0){

switch(string_to_case.at(t_impute_method)) {
  case 1:
    imputeG = std::round(2 * t_altFreq);
    //std::cout << "t_impute_method " << t_impute_method << std::endl;
    break;
  case 2:
    imputeG = 2 * t_altFreq;
    //std::cout << "t_impute_method " << t_impute_method << std::endl;
    break;
  case 3:
    imputeG = 0;
    //std::cout << "t_impute_method " << t_impute_method << std::endl;
    break;
}


  for(int i = 0; i < nMissing; i++){
    uint32_t j = t_indexForMissing.at(i);
    t_GVec.at(j) = imputeG;
  }

 t_MAC = t_MAC + imputeG * nMissing;

}

  
  if(t_dosage_zerod_cutoff > 0){ 
    if(t_MAC <= t_dosage_zerod_MAC_cutoff){
      t_GVec.clean(t_dosage_zerod_cutoff);
    }	  
  } 

  //if(nMissing > 0 || t_dosage_zerod_cutoff > 0){
     t_altCount = arma::sum(t_GVec);
     t_altFreq = t_altCount / (2*dosagesSize);
     if(flip){
	t_altFreq = 1 - t_altFreq;
        t_altCount = 2*dosagesSize - t_altCount;	
	//t_altCount = 2 * t_altFreq * dosagesSize;
     }	     
  //}

 // 2026-09-16: the marker path passes t_buildIndexZero = false. Its only
 // consumers of the zero-genotype index were its SIZE (m_n minus the non-zero
 // count, an identity the test now uses directly) and the ER exact test, which
 // builds the complement itself when handed an empty one -- so the N push_backs
 // and the N-word copy that followed them were paid on every marker for
 // nothing. The region and conditional paths still build it.
 for(unsigned int i = 0; i < dosagesSize; i++){
   if(t_GVec(i) == 0){
     if(t_buildIndexZero) t_indexZero.push_back(i);
   }else{
     t_indexNonZero.push_back(i);
   }
 }

  return flip;
}


// The identities this rests on, each against the dense body above:
//  * the flip is decided on the pre-impute altFreq, and a difference list's
//    common code of 2 ALWAYS flips (hom-ALT holds at least 3N/4 of the samples
//    so altFreq >= 3/4) while a common code of 0 NEVER does (at most N/4
//    listed, so altFreq <= 1/3): after the fold every unlisted sample is an
//    exact 0. The caller admits no other common code (GenoReader.hpp).
//  * `2 - t_GVec` on the dense vector is 2.0 - code on each listed entry;
//    a missing entry is NaN there and imputeG here, both then overwritten /
//    set to imputeG, which is computed from the FLIPPED altFreq as the dense
//    body does.
//  * `.clean(cutoff)` zeroes |x| <= cutoff; an unlisted 0 is unchanged by it.
//  * arma::sum is arrayops::accumulate: two accumulators split by index
//    parity, added at the end (drilled 2026-09-12 for CarrierStore). The
//    unlisted zeros add nothing to either, so summing the listed values into
//    the accumulator their sample index selects, in increasing index order, is
//    the same sequence of additions.
//  * the non-zero index is built in increasing sample order, as the dense scan.
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
                         std::vector<double> & t_valNonZero)
{
  if(t_commonCode != 0 && t_commonCode != 2){
    Rcpp::stop("imputeSparseAndFlip: common genotype code " + std::to_string(t_commonCode) +
               " is not 0 or 2; the reader must have expanded this record (GenoReader.hpp)");
  }
  bool flip = false;
  t_indexNonZero.clear();
  t_valNonZero.clear();
  const uint dosagesSize = t_n;
  const int nMissing = (int)t_nMissing;
  double imputeG = 0;
  if(t_altFreq > 0.5){
    flip = true;
    t_altFreq = 1 - t_altFreq;
  }
  if(nMissing > 0){
    switch(string_to_case.at(t_impute_method)) {
      case 1: imputeG = std::round(2 * t_altFreq); break;
      case 2: imputeG = 2 * t_altFreq; break;
      case 3: imputeG = 0; break;
    }
    t_MAC = t_MAC + imputeG * nMissing;
  }
  const bool doClean = (t_dosage_zerod_cutoff > 0) && (t_MAC <= t_dosage_zerod_MAC_cutoff);
  double acc1 = 0, acc2 = 0;
  const std::size_t k = t_idx.size();
  for(std::size_t e = 0; e < k; e++){
    const uint8_t c = t_code[e];
    double v;
    if(c == 3){
      v = imputeG;
    }else{
      v = flip ? (2 - (double)c) : (double)c;
    }
    if(doClean && std::abs(v) <= t_dosage_zerod_cutoff){
      v = 0;
    }
    if(v != 0){
      const uint i = t_idx[e];
      t_indexNonZero.push_back(i);
      t_valNonZero.push_back(v);
      if((i % 2) == 0){ acc1 += v; }else{ acc2 += v; }
    }
  }
  t_altCount = acc1 + acc2;
  t_altFreq = t_altCount / (2*dosagesSize);
  if(flip){
    t_altFreq = 1 - t_altFreq;
    t_altCount = 2*dosagesSize - t_altCount;
  }
  return flip;
}


double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat)
{
  double innerProd = arma::accu(x1Mat % x2Mat);
  return(innerProd);
}

arma::vec getTime(){
  arma::vec Time(2, arma::fill::zeros);
  struct timeval time;
  Time(0) = 0;
  if(!gettimeofday(&time,NULL))
    Time(0) = (double)time.tv_sec + (double)time.tv_usec * .000001;
  Time(1) = (double)clock() / CLOCKS_PER_SEC;
  return Time;
}

void printTime(arma::vec t1, arma::vec t2, std::string message){
  double wallTime = t2(0) - t1(0);
  double cpuTime = t2(1) - t1(1);
  if(wallTime < 60){
    Rprintf ("It took %f seconds (%f CPU seconds) to %s.\n",
             wallTime, cpuTime, message.c_str());
  }else if(wallTime < 3600){
    Rprintf ("It took %f minutes (%f CPU minutes) to %s.\n",
             wallTime/60, cpuTime/60, message.c_str());
  }else{
    Rprintf ("It took %f hours (%f CPU hours) to %s.\n",
             wallTime/3600, cpuTime/3600, message.c_str());
  }
}

double getinvStd(double t_freq)
{
  double Std = sqrt(2 * t_freq * (1-t_freq));
  if(Std == 0)
    return 0;
  else
    return 1/Std;
}

// [[Rcpp::export]]
double sum_arma1(arma::vec& X) {
    double sum = 0;
    for (uint i = 0; i < X.n_elem; ++i) {
        if (arma::is_finite(X(i)))
            sum += X(i);
    }
    return sum;
}


// [[Rcpp::export]]
double add_logp(double p1, double p2)
{
        using namespace Rcpp;
        p1 = -std::abs(p1);
        p2 = -std::abs(p2);
        double maxp = std::max(p1,p2);
        double  minp = std::min(p1,p2);
        double result = maxp+std::log(1+std::exp(minp-maxp));
        return(result);
}
