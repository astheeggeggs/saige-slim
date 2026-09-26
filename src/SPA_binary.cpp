// [[Rcpp::depends(BH)]]

//[[Rcpp::depends(RcppArmadillo)]]
#include <functional>
#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <sstream>
#include <time.h>
#include <Rcpp.h>
#include <stdint.h>
#include <boost/math/distributions/normal.hpp>
#include "SPA_binary.hpp"
#include "UTIL.hpp"
//boost::math::normal norm;




// [[Rcpp::export]]
double Korg_Binom(double t1, const arma::vec & mu, const arma::vec & g)
{
	SPAWorkspace ws(mu.n_elem);
	return Korg_Binom_ws(t1, mu, g, ws);
}

double Korg_Binom_ws(double t1, const arma::vec & mu, const arma::vec & g, SPAWorkspace & ws)
{
	const arma::uword n = mu.n_elem;
	ws.ensure(n);
	arma::vec temp(ws.a.memptr(), n, false, true);
	temp = arma::log(1 - mu + mu % arma::exp(g * t1));
        double out = arma::sum(temp);
	return(out);
}


// [[Rcpp::export]]
double K1_adj_Binom(double t1, const arma::vec & mu, const arma::vec & g, double q)
{
	SPAWorkspace ws(mu.n_elem);
	return K1_adj_Binom_ws(t1, mu, g, q, ws);
}

double K1_adj_Binom_ws(double t1, const arma::vec & mu, const arma::vec & g, double q, SPAWorkspace & ws)
{
	const arma::uword n = mu.n_elem;
	ws.ensure(n);
	arma::vec temp1(ws.a.memptr(), n, false, true);
	arma::vec temp2(ws.b.memptr(), n, false, true);
	arma::vec temp3(ws.c.memptr(), n, false, true);

	temp1 = (1 - mu) % arma::exp(-g * t1) + mu;
	temp2 = mu % g;
	temp3 = temp2/temp1;
	double out  = arma::sum(temp3) - q;
	return(out);
}


// [[Rcpp::export]]
double K2_Binom(double t1, const arma::vec & mu, const arma::vec & g)
{
	SPAWorkspace ws(mu.n_elem);
	return K2_Binom_ws(t1, mu, g, ws);
}

double K2_Binom_ws(double t1, const arma::vec & mu, const arma::vec & g, SPAWorkspace & ws)
{
	const arma::uword n = mu.n_elem;
	ws.ensure(n);
        arma::vec temp0(ws.a.memptr(), n, false, true);
        arma::vec temp1(ws.b.memptr(), n, false, true);
        arma::vec temp2(ws.c.memptr(), n, false, true);
        arma::vec temp3(ws.d.memptr(), n, false, true);
		
	temp0 = arma::exp(-g * t1);
        temp1 = (1 - mu) % temp0 + mu;
	temp1 = arma::pow(temp1, 2);
       	temp2 = arma::pow(g,2) % temp0;			
        temp2 = (1-mu) % mu % temp2;
        temp3 = temp2/temp1;
        double out = sum_arma1(temp3);

        return(out);
}





// The two full-length reductions over g, factored out so a caller that runs
// getroot twice on the same g pays for them once. The expressions are the
// originals, unchanged, so the values are bit-identical.
void getGPosNeg_Binom(const arma::vec & g, double & gpos, double & gneg){
	gpos = arma::accu( g.elem( find(g > 0) ) );
	gneg = arma::accu( g.elem( find(g < 0) ) );
}

// [[Rcpp::export]]
Rcpp::List getroot_K1_Binom(double init, const arma::vec & mu, const arma::vec & g, double q, double tol, int maxiter){
	double gpos, gneg;
	getGPosNeg_Binom(g, gpos, gneg);
	SPAWorkspace ws(mu.n_elem);
	SPARoot out = getroot_K1_Binom_core(init, mu, g, q, tol, gpos, gneg, ws, maxiter);
	Rcpp::List result;
	result["root"]	= out.root;
	result["niter"] = out.niter;
	result["Isconverge"] = out.Isconverge;
	return(result);
}


// LEDGER #163 (2026-09-26). The Newton iteration above starts at t = 0 with
// the step (q - K1(0)) / K2(0), and K2(0) is the score variance. When the
// carriers' fitted probabilities are tiny that step runs past the point where
// exp(-g t) underflows: K2 evaluates to exactly 0, the next iterate is -Inf,
// and 0 * Inf in the following K1 is NaN, which the #45 guard reads as
// non-convergence -- for a root that exists (q strictly inside (gneg, gpos),
// K1 monotone). This rescue runs ONLY after that iteration has failed, so
// every root the plain Newton finds keeps its bytes: it brackets the sign
// change of K1 - q by doubling from 0, then takes Newton steps that are
// replaced by bisection whenever they leave the bracket or K2 is unusable.
// The bracket cannot fail on an attainable q, because K1(t) -> gpos as
// t -> +Inf and -> gneg as t -> -Inf, monotonically.
struct SPARescue { double root; int niter; bool ok; };
static SPARescue rescueRootK1(const std::function<double(double)> & K1adj,
                              const std::function<double(double)> & K2,
                              double tol, int maxiter){
	SPARescue r; r.root = std::numeric_limits<double>::infinity(); r.niter = 0; r.ok = false;
	double f0 = K1adj(0.0);
	if(!std::isfinite(f0)){ return r; }
	if(f0 == 0.0){ r.root = 0.0; r.ok = true; return r; }
	const double dir = (f0 < 0.0) ? 1.0 : -1.0;   // K1 - q rises with t
	double lo = 0.0, flo = f0, hi = dir, fhi = K1adj(hi);
	int n = 1;
	while(std::isfinite(fhi) && fhi * flo > 0.0 && std::abs(hi) < 1048576.0){
		lo = hi; flo = fhi; hi = 2.0 * hi; fhi = K1adj(hi); n++;
	}
	if(!std::isfinite(fhi) || fhi * flo > 0.0){ return r; }
	if(lo > hi){ std::swap(lo, hi); std::swap(flo, fhi); }
	double t = 0.5 * (lo + hi);
	for(int it = 0; it < maxiter; it++){
		double ft = K1adj(t);
		if(!std::isfinite(ft)){ return r; }
		if(ft == 0.0){ r.root = t; r.niter = n + it; r.ok = true; return r; }
		if(ft * flo < 0.0){ hi = t; fhi = ft; } else { lo = t; flo = ft; }
		double k2 = K2(t);
		double tnew = t - ft / k2;
		if(!(k2 > 0.0) || !std::isfinite(tnew) || tnew <= lo || tnew >= hi){
			tnew = 0.5 * (lo + hi);   // bisection whenever Newton is unusable
		}
		if(std::abs(tnew - t) < tol || (hi - lo) < tol){
			r.root = tnew; r.niter = n + it + 1; r.ok = true; return r;
		}
		t = tnew;
	}
	return r;
}

SPARoot getroot_K1_Binom_core(double init, const arma::vec & mu, const arma::vec & g, double q, double tol, double gpos, double gneg, SPAWorkspace & ws, int maxiter){
	SPARoot result;
	double root;
	int niter;
	bool Isconverge;
        int rep;
	double K1_eval, K2_eval, t, tnew, newK1;
	double prevJump;
	//std::cout << "q is " << q << std::endl;
	//std::cout << "gpos is " << gpos << std::endl;
	//std::cout << "gneg is " << gneg << std::endl;
	if(q >= gpos || q <= gneg){
		// LEDGER #22 FIXED: an unattainable q (all-zero g, or q outside K1's
		// range) has no root; report NON-convergence with the infinite root, so
		// the caller takes its no-saddlepoint fallback deliberately rather than
		// after finding k1/k2 non-finite. Same composed p-value either way.
		root = std::numeric_limits<double>::infinity();
		niter = 0;
		Isconverge = false;
	} else{	
		t = init;
		K1_eval = K1_adj_Binom_ws(t,mu,g,q,ws);
		prevJump = std::numeric_limits<double>::infinity();
		int rep = 1;
		bool conv = true;
		while(rep <= maxiter){
			K2_eval = K2_Binom_ws(t,mu,g,ws);
			tnew = t-K1_eval/K2_eval;
			if(std::isnan(tnew)){ // LEDGER #45: NA_REAL is a NaN and == NaN is always false; ISNAN was meant
				conv = false;
				break;

			}
			
			if(std::abs(tnew-t)<tol){	
				conv = true;
				break;	
			}
			
			if(rep == maxiter)
                        {
                                conv = false;
                                break;
                        }

			newK1 = K1_adj_Binom_ws(tnew,mu,g,q,ws);
                        if((K1_eval * newK1) < 0) // LEDGER #46 FIXED: the fast kernel's test; sign(0) made an exact root look like a sign change
                        {
                                if(std::abs(tnew-t) > (prevJump-tol))
                                {
                                        tnew = t + (arma::sign(newK1-K1_eval))*prevJump/2;
                                        newK1 = K1_adj_Binom_ws(tnew,mu,g,q,ws);
                                        prevJump = prevJump/2;
                                } else {
                                        prevJump = std::abs(tnew-t);
                                }
                        }

			rep = rep + 1;
			t = tnew;
                        K1_eval = newK1;
		}
		if(!conv){ // LEDGER #163: rescue only after the plain Newton has failed
			SPARescue rs = rescueRootK1([&](double tt){ return K1_adj_Binom_ws(tt, mu, g, q, ws); },
				[&](double tt){ return K2_Binom_ws(tt, mu, g, ws); }, tol, maxiter);
			if(rs.ok){ t = rs.root; rep = rep + rs.niter; conv = true; }
		}
		root=t;
		niter=rep;
		Isconverge=conv;	
	}
	//std::cout << "here13" << std::endl;
        //std::cout << "root: " << root << std::endl;
        //std::cout << "niter: " << niter << std::endl;
	result.root	= root;
	result.niter = niter;
	result.Isconverge = Isconverge;
	return(result);
}



// [[Rcpp::export]]
Rcpp::List Get_Saddle_Prob_Binom(double zeta,  const arma::vec & mu, const arma::vec & g, double q, bool logp)
{
	SPAWorkspace ws(mu.n_elem);
	SPASaddle out = Get_Saddle_Prob_Binom_core(zeta, mu, g, q, logp, ws);
	Rcpp::List result;
	result["pval"] = out.pval;
	result["isSaddle"] = out.isSaddle;
	return(result);
}


SPASaddle Get_Saddle_Prob_Binom_core(double zeta,  const arma::vec & mu, const arma::vec & g, double q, bool logp, SPAWorkspace & ws)
{
	double k1 = Korg_Binom_ws(zeta, mu, g, ws);
	double k2 = K2_Binom_ws(zeta, mu, g, ws);
	double temp1, w, v, Ztest, pval;
	double negative_infinity = - std::numeric_limits<double>::infinity();

	temp1 = zeta * q - k1;
	SPASaddle result;
	bool isSaddle = false;
	//std::cout << "k1 " << k1 << std::endl;
	//std::cout << "k2 " << k2 << std::endl;
	//std::cout << "temp1 " << temp1 << std::endl;
	//std::cout << "zeta " << zeta << std::endl;
	//std::cout << "q " << q << std::endl;


        bool flagrun=false;
	if(std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0){
		 w = arma::sign(zeta) * std::sqrt(2 *temp1);
		 v = zeta *  std::sqrt(k2);
		 if(w != 0){
			flagrun = true;
		 } 	 
	}	

	//w = arma::sign(zeta) * std::sqrt(2 *temp1);
	//v = zeta *  std::sqrt(k2);

	//if(std::isfinite(k1) && std::isfinite(k2) && temp1 > 0 && w != 0 && (arma::sign(v) == arma::sign(w)))
	if(flagrun)
	{
		//temp1 = zeta * q - k1;
		Ztest = w + (1/w) * std::log(v/w);

		//std::cout << "Ztest: " << Ztest << std::endl;
		//std::cout << "w: " << w << std::endl;
		//std::cout << "v: " << v << std::endl;
   		//boost::math::normal norm_dist(0,1);

		double pval0;
	        if(Ztest > 0){
                        //pval0 = boost::math::cdf(complement(norm_dist,Ztest));
                        pval0 = R::pnorm(Ztest,0,1,false,logp);
                        pval = pval0;
                        //if(logp){
                        //      pval = std::log(pval);
                        //}
                }else {
                        //pval0 = boost::math::cdf(norm_dist,Ztest);
                        pval0 = R::pnorm(Ztest,0,1,true,logp);
                        pval = -pval0;

                        //      -Rcpp::pnorm( Ztest, mean = 0.0, sd = 1.0, lower = true, log = logp );
                }	

		isSaddle = true;
	} else {
			if(logp)
			{
				pval =  negative_infinity;
			}else{
				pval= 0;
			}
	}
	result.pval = pval;
	result.isSaddle = isSaddle;	
	return(result);
}



// [[Rcpp::export]]
Rcpp::List SPA_binary(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp){
	using namespace Rcpp;
	List result;
	double p1, p2, pval;
	bool Isconverge = true;
	Rcpp::List outuni1 = getroot_K1_Binom(0, mu, g, q, tol);
	Rcpp::List outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol);
	//double outuni1root = outuni1["root"];
	//double outuni2root = outuni2["root"];
	//bool Isconverge1 = outuni1["Isconverge"];
	//bool Isconverge2 = outuni2["Isconverge"];

	//std::cout << "outuni1root" << outuni1root << std::endl;
        //std::cout << "outuni2root" << outuni2root << std::endl;
        //std::cout << "Isconverge1" << Isconverge1 << std::endl;
        //std::cout << "Isconverge2" << Isconverge2 << std::endl;


	Rcpp::List getSaddle;
        Rcpp::List getSaddle2;
	if(outuni1["Isconverge"]  && outuni2["Isconverge"])
	{
		//std::cout << "q is " << q << " 3 qinv is " << qinv << std::endl;
		getSaddle = Get_Saddle_Prob_Binom(outuni1["root"], mu, g, q, logp);

		if(getSaddle["isSaddle"]){
			p1 = getSaddle["pval"];
		}else{	

			if(logp){
				p1 = pval_noadj-std::log(2);
			}else{
				p1 = pval_noadj/2;	
			}	
		}
		getSaddle2 = Get_Saddle_Prob_Binom(outuni2["root"], mu, g, qinv, logp);
		if(getSaddle2["isSaddle"]){
                        p2 = getSaddle2["pval"];
                }else{
			if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
		}
		
		if(logp)
		{
			pval = add_logp(p1,p2);
		} else {
			pval = std::abs(p1)+std::abs(p2);
		}
		Isconverge=true;
	}else {
			pval = pval_noadj;
			Isconverge=false;	
		}		
	result["pvalue"] = pval;
	result["Isconverge"] = Isconverge;
	return(result);
}


// [[Rcpp::export]]
double Korg_fast_Binom(double t1, const arma::vec & mu, const arma::vec & g,  const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB, double NAmu, double NAsigma)
{
	SPAWorkspace ws(gNB.n_elem);
	return Korg_fast_Binom_ws(t1, gNB, muNB, NAmu, NAsigma, ws);
}

double Korg_fast_Binom_ws(double t1, const arma::vec & gNB, const arma::vec & muNB, double NAmu, double NAsigma, SPAWorkspace & ws)
{
	const arma::uword m = gNB.n_elem;
	ws.ensure(m);
	arma::vec temp(ws.a.memptr(), m, false, true);
	temp = arma::log(1 - muNB + muNB % (arma::exp(gNB * t1)));
        double out = arma::sum(temp) + NAmu*t1+ 0.5*NAsigma*pow(t1,2);
	return(out);
}


// [[Rcpp::export]]
double K1_adj_fast_Binom(double t1, const arma::vec & mu, const arma::vec & g, double q, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma)
{
	SPAWorkspace ws(gNB.n_elem);
	return K1_adj_fast_Binom_ws(t1, q, gNB, muNB, NAmu, NAsigma, ws);
}

double K1_adj_fast_Binom_ws(double t1, double q, const arma::vec & gNB, const arma::vec & muNB, double NAmu, double NAsigma, SPAWorkspace & ws)
{
	const arma::uword m = gNB.n_elem;
	ws.ensure(m);
	arma::vec temp1(ws.a.memptr(), m, false, true);
	arma::vec temp2(ws.b.memptr(), m, false, true);
	double temp3;
	arma::vec temp4(ws.c.memptr(), m, false, true);

	temp1 = (1 - muNB) % arma::exp(-gNB * t1) + muNB;
	temp2 = muNB % gNB;
	temp3 = NAmu+NAsigma*t1;
	temp4 = temp2/temp1;
	double out  = arma::sum(temp4) + temp3 - q;
	return(out);
}



// [[Rcpp::export]]
double K2_fast_Binom(double t1, const arma::vec & mu, const arma::vec & g, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma)
{
	SPAWorkspace ws(gNB.n_elem);
	return K2_fast_Binom_ws(t1, gNB, muNB, NAsigma, ws);
}

double K2_fast_Binom_ws(double t1, const arma::vec & gNB, const arma::vec & muNB, double NAsigma, SPAWorkspace & ws)
{
	const arma::uword m = gNB.n_elem;
	ws.ensure(m);
        arma::vec temp0(ws.a.memptr(), m, false, true);
        arma::vec temp1(ws.b.memptr(), m, false, true);
        arma::vec temp2(ws.c.memptr(), m, false, true);
        arma::vec temp3(ws.d.memptr(), m, false, true);
		
	temp0 = arma::exp(-gNB * t1);
        temp1 = (1 - muNB) % temp0 + muNB;
	temp1 = pow(temp1, 2);
       	temp2 = arma::pow(gNB,2) % temp0;			
        temp2 = (1-muNB) % muNB % temp2;
        temp3 = temp2/temp1;
        double out = sum_arma1(temp3)+NAsigma;
        return(out);
}





// [[Rcpp::export]]
Rcpp::List getroot_K1_fast_Binom(double init, const arma::vec & mu, const arma::vec & g, double q, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, int maxiter){
	double gpos, gneg;
	getGPosNeg_Binom(g, gpos, gneg);
	SPAWorkspace ws(gNB.n_elem);
	SPARoot out = getroot_K1_fast_Binom_core(init, g, q, gNB, muNB, NAmu, NAsigma, tol, gpos, gneg, ws, maxiter);
	Rcpp::List result;
	result["root"]	= out.root;
	result["niter"] = out.niter;
	result["Isconverge"] = out.Isconverge;
	return(result);
}

SPARoot getroot_K1_fast_Binom_core(double init, const arma::vec & g, double q, const arma::vec & gNB, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, double gpos, double gneg, SPAWorkspace & ws, int maxiter){
	SPARoot result;
	double root;
	int niter;
	bool Isconverge;
        int rep;
	double K1_eval, K2_eval, t, tnew, newK1;
	double prevJump;
	if(q >= gpos || q <= gneg){
		// LEDGER #22 FIXED: an unattainable q (all-zero g, or q outside K1's
		// range) has no root; report NON-convergence with the infinite root, so
		// the caller takes its no-saddlepoint fallback deliberately rather than
		// after finding k1/k2 non-finite. Same composed p-value either way.
		root = std::numeric_limits<double>::infinity();
		niter = 0;
		Isconverge = false;
	} else{	
		t = init;
		K1_eval = K1_adj_fast_Binom_ws(t,q,gNB,muNB,NAmu, NAsigma,ws);
		prevJump = std::numeric_limits<double>::infinity();
		int rep = 1;
		bool conv = true;
		while(rep <= maxiter){
			K2_eval = K2_fast_Binom_ws(t,gNB,muNB,NAsigma,ws);
			tnew = t-K1_eval/K2_eval;
			if(std::isnan(tnew)){ // LEDGER #45: NA_REAL is a NaN and == NaN is always false; ISNAN was meant
				conv = false;
				break;

			}
			
			if(std::abs(tnew-t)<tol){	
				conv = true;
				break;	
			}
			
			if(rep == maxiter)
                        {
                                conv = false;
                                break;
                        }

			newK1 = K1_adj_fast_Binom_ws(tnew,q,gNB,muNB,NAmu, NAsigma,ws);
                        if((K1_eval * newK1) < 0)
                        {
                                if(std::abs(tnew-t) > (prevJump-tol))
                                {
                                        tnew = t + (arma::sign(newK1-K1_eval))*prevJump/2;
                                        newK1 = K1_adj_fast_Binom_ws(tnew,q,gNB,muNB,NAmu, NAsigma,ws);
                                        prevJump = prevJump/2;
                                } else {
                                        prevJump = std::abs(tnew-t);
                                }
                        }

			rep = rep + 1;
			t = tnew;
                        K1_eval = newK1;
		}
		if(!conv){ // LEDGER #163: rescue only after the plain Newton has failed
			SPARescue rs = rescueRootK1([&](double tt){ return K1_adj_fast_Binom_ws(tt, q, gNB, muNB, NAmu, NAsigma, ws); },
				[&](double tt){ return K2_fast_Binom_ws(tt, gNB, muNB, NAsigma, ws); }, tol, maxiter);
			if(rs.ok){ t = rs.root; rep = rep + rs.niter; conv = true; }
		}
		root=t;
		niter=rep;
		Isconverge=conv;	
	}
	result.root	= root;
	result.niter = niter;
	result.Isconverge = Isconverge;
	return(result);
}



// [[Rcpp::export]]
Rcpp::List Get_Saddle_Prob_fast_Binom(double zeta,  const arma::vec & mu, const arma::vec & g, double q, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, bool logp)
{
	SPAWorkspace ws(gNB.n_elem);
	SPASaddle out = Get_Saddle_Prob_fast_Binom_core(zeta, g, q, gNB, muNB, NAmu, NAsigma, logp, ws);
	Rcpp::List result;
	result["pval"] = out.pval;
	result["isSaddle"] = out.isSaddle;
	return(result);
}


SPASaddle Get_Saddle_Prob_fast_Binom_core(double zeta,  const arma::vec & g, double q, const arma::vec & gNB, const arma::vec & muNB,  double NAmu, double NAsigma, bool logp, SPAWorkspace & ws)
{
	double k1 = Korg_fast_Binom_ws(zeta, gNB, muNB, NAmu, NAsigma, ws);
	double k2 = K2_fast_Binom_ws(zeta, gNB, muNB, NAsigma, ws);
	double temp1, w, v, Ztest, pval;
	double negative_infinity = - std::numeric_limits<double>::infinity();
	temp1 = zeta * q - k1;
	SPASaddle result;
	bool isSaddle = false; 


	bool flagrun=false;
        if(std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0){
                 w = arma::sign(zeta) * std::sqrt(2 *temp1);
                 v = zeta *  std::sqrt(k2);
                 if(w != 0){
                        flagrun = true;
                 }
        }



	//if(std::isfinite(k1) && std::isfinite(k2) && temp1 > 0 )
	if(flagrun)
	{
		//temp1 = zeta * q - k1;
		//w = arma::sign(zeta) * std::sqrt(2 *temp1);
		//v = zeta *  std::sqrt(k2);
		Ztest = w + (1/w) * std::log(v/w);
	
		boost::math::normal norm_dist(0,1);
                double pval0;



		if(Ztest > 0){
			//pval0 = boost::math::cdf(complement(norm_dist,Ztest));
			pval0 = R::pnorm(Ztest,0,1,false,logp);
			pval=pval0;
			//if(logp){
			//	pval = std::log(pval);
			//}	
		}else {
			//pval0 = boost::math::cdf(norm_dist,Ztest);
			pval0 = R::pnorm(Ztest,0,1,true,logp);
			pval= -pval0;
			
			//	-Rcpp::pnorm( Ztest, mean = 0.0, sd = 1.0, lower = true, log = logp );
		}
		isSaddle = true;
	}else{
		if(logp)
		{
			pval =  negative_infinity;
		}else {
			pval=0;
		}
	}
	result.pval = pval;
        result.isSaddle = isSaddle;
        return(result);


}



// [[Rcpp::export]]
Rcpp::List SPA_binary_fast(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, bool logp, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, double tol){
	using namespace Rcpp;
	List result;
	double p1, p2, pval;
	bool Isconverge = true;
	Rcpp::List outuni1 = getroot_K1_fast_Binom(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
	//double qinv = -1*q;
	Rcpp::List outuni2 = getroot_K1_fast_Binom(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
	//std::cout << "outuni1root" << outuni1["root"] << std::endl;
	//std::cout << "outuni2root" << outuni2["root"] << std::endl;

	Rcpp::List getSaddle;
	Rcpp::List getSaddle2;
	if(outuni1["Isconverge"]  && outuni2["Isconverge"])
	{
		
		getSaddle  = Get_Saddle_Prob_fast_Binom(outuni1["root"], mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
		if(getSaddle["isSaddle"]){
			p1 = getSaddle["pval"];
		}else{	

		        if(logp){
				p1 = pval_noadj-std::log(2);
			}else{
				p1 = pval_noadj/2;	
			}	
		}

		getSaddle2 = Get_Saddle_Prob_fast_Binom(outuni2["root"], mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
		if(getSaddle2["isSaddle"]){
			p2 = getSaddle2["pval"];
		
		}else{
			if(logp){
				p2 = pval_noadj-std::log(2);
			}else{
				p2 = pval_noadj/2;	
			}	
		}

		std::cout << "p1  first " << p1 << "p2 " << p2 << std::endl;
		std::cout << "HEREHERE " << std::endl;
		if(logp){
			pval = add_logp(p1,p2);
		} else {
			pval = std::abs(p1)+std::abs(p2);
			//std::cout << "p1 " << p1 << "p2 " << p2 << std::endl; 
		}
		Isconverge=true;
	}else {
 			//std::cout << "Error_Converge" << std::endl;
			pval = pval_noadj;
			Isconverge=false;	
		}		
	result["pvalue"] = pval;
	result["Isconverge"] = Isconverge;
	return(result);
}

