// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>         // std::chrono::seconds
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/date_time.hpp>
#include <boost/math/distributions/cauchy.hpp>

#define _USE_MATH_DEFINES
#include <math.h> 


// [[Rcpp::export]]
double CCT_cpp(arma::vec & pval){
    double cauchyp, cctstat;
    uint np = pval.n_elem;    
    //check if there is NA	
    bool isna = pval.has_nan();
    if(isna){
        std::cout << "Cannot have NAs in the p-values!" << std::endl;
    }

    // check if all p-values are between 0 and 1

    arma::uvec between0and1Indice = arma::find(pval <= 1 && pval >= 0);
    if (between0and1Indice.n_elem != np){
        exit(EXIT_FAILURE); 
    }

    bool isnonzero = arma::all(pval);
    if(!isnonzero){
       cauchyp = 0;  
    }

    //if any of the p-value is one
    arma::uvec notoneIndice = arma::find(pval != 1.0);
    if(notoneIndice.n_elem < np){
	double minpval = pval.min();
	minpval = minpval * np;
	cauchyp = std::min(1.0, minpval);
    }else{	    
    double weight = 1.0/double(np);
    arma::vec weights(np, arma::fill::zeros);
    weights.replace(0, weight);

    arma::uvec issmallIndice = arma::find(pval < 1e-16);
    arma::uvec isNotSmallIndice = arma::find(pval >= 1e-16);
    if(issmallIndice.n_elem == 0){


	arma::vec tantemp = arma::tan((0.5 - pval) * M_PI);
	arma::vec weighttantemp = tantemp % weights;
        cctstat = arma::sum(weighttantemp);
    }else{
       arma::vec weightsub = weights.elem(issmallIndice);
       arma::vec pvalsub = pval.elem(issmallIndice);
       cctstat = arma::sum((weightsub/pvalsub)/M_PI);
       arma::vec pvalsub_b = pval.elem(isNotSmallIndice);
       arma::vec weightsub_b = weights.elem(isNotSmallIndice); 
       arma::vec tantemp_b = (0.5 - pvalsub_b) * M_PI;	
	tantemp_b = arma::tan(tantemp_b);
	arma::vec weighttantemp_b = tantemp_b % weightsub_b;
	cctstat = cctstat + arma::sum(weighttantemp_b);
    }

    if(cctstat > 1e+15){
       cauchyp = (1/cctstat)/M_PI;
    }else{
       boost::math::cauchy cauchy_dist(0,1);
       cauchyp = boost::math::cdf(complement(cauchy_dist, cctstat));
    }
    }//if(notoneIndice.n_elem < np){ 
    return(cauchyp);    
}




	


// --is_output_log10P (2026-10-14): the Cauchy combination in log space, the
// twin of R/CCT_modified.R's CCT_log10 (same rules, same branches): a p of 1
// anywhere gives min(1, k min p); terms below 1e-16 are 1/(pi p), summed as a
// log-sum-exp; the rest are tan((0.5 - p) pi) added linearly; a statistic above
// 1e15 is inverted as 1/(pi T), else the Cauchy CDF. Equal weights. NaN in ->
// NaN out (the caller drops NA cells first).
// [[Rcpp::export]]
double CCT_log10_cpp(arma::vec & log10p){
    const arma::uword k = log10p.n_elem;
    if(k == 0 || log10p.has_nan()){ return std::numeric_limits<double>::quiet_NaN(); }
    if(arma::any(log10p > 0)){ Rcpp::stop("All p-values must be between 0 and 1!"); }
    if(arma::any(log10p == 0)){ return std::min(0.0, log10p.min() + std::log10((double)k)); }
    double lin = 0.0;
    std::vector<double> big;
    for(arma::uword i = 0; i < k; i++){
        if(log10p(i) < -16.0){ big.push_back(-log10p(i) - std::log10(M_PI)); }
        else { lin += std::tan((0.5 - std::pow(10.0, log10p(i))) * M_PI); }
    }
    // the Cauchy upper tail without cancellation (R: tail10): atan(1/t)/pi
    // for t > 0 is exactly 0.5 - atan(t)/pi and keeps its digits at large t
    auto tail10 = [](double stat) -> double {
        if(stat > 1e+15){ return -(std::log10(stat) + std::log10(M_PI)); }
        if(stat > 0.0){ return std::log10(std::atan(1.0 / stat) / M_PI); }
        return std::log10(0.5 - std::atan(stat) / M_PI);
    };
    if(big.empty()){ return tail10(lin / (double)k); }
    double m = *std::max_element(big.begin(), big.end());
    double acc = 0.0;
    for(double b : big){ acc += std::pow(10.0, b - m); }
    double log10S = m + std::log10(acc);
    if(lin != 0.0 && log10S < 300.0){
        double ratio = lin / std::pow(10.0, log10S);
        if(1.0 + ratio <= 0.0){ return tail10((std::pow(10.0, log10S) + lin) / (double)k); }
        log10S += std::log10(1.0 + ratio);
    }
    double log10T = log10S - std::log10((double)k);
    if(log10T > 15.0){ return -(log10T + std::log10(M_PI)); }
    return tail10(std::pow(10.0, log10T));
}
