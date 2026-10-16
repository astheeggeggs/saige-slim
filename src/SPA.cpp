// [[Rcpp::depends(BH)]]

//[[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <limits>
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
#include "SPA.hpp"
#include "UTIL.hpp"

// [[Rcpp::export]]
void SPA(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge){
        SPAWorkspace ws(mu.n_elem);
        SPA_ws(mu, g, q, qinv, pval_noadj, tol, logp, traitType, pval, isSPAConverge, ws);
}

void SPA_ws(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge, SPAWorkspace & ws){
        double p1, p2;
        bool Isconverge = true;
        ws.failReason = 0;
        if( traitType != "binary"){
                // Unreachable, and deliberately not an Rcpp::stop. Only the
                // binary saddlepoint kernels exist, and SAIGE_test.cpp gates the
                // whole SPA branch on traitType != "quantitative", which since
                // survival was removed in Phase 3a leaves "binary" as the only
                // value that can arrive here. Before S1 a non-binary trait fell
                // through to reading "Isconverge" out of a default-constructed
                // (empty) Rcpp::List, which throws from Rcpp's name lookup; this
                // keeps it a throw, but one that enters no R API and is safe to
                // reach from a worker thread.
                throw std::logic_error("SPA(): only the binary saddlepoint kernels exist");
        }
        // gpos/gneg depend only on g, and both roots use the same g
        double gpos, gneg;
        getGPosNeg_Binom(g, gpos, gneg);
        SPARoot outuni1 = getroot_K1_Binom_core(0, mu, g, q, tol, gpos, gneg, ws);
        SPARoot outuni2 = getroot_K1_Binom_core(0, mu, g, qinv, tol, gpos, gneg, ws);
        ws.diagRoot1 = outuni1.root; ws.diagNiter1 = outuni1.niter; ws.diagRoot2 = outuni2.root; ws.diagNiter2 = outuni2.niter; ws.diagBound = -1; ws.diagQSaddle = -1;

/*
        double outuni1root = outuni1["root"];
        double outuni2root = outuni2["root"];
        bool Isconverge1 = outuni1["Isconverge"];
        bool Isconverge2 = outuni2["Isconverge"];

        std::cout << "outuni1root" << outuni1root << std::endl;
        std::cout << "outuni2root" << outuni2root << std::endl;
        std::cout << "Isconverge1" << Isconverge1 << std::endl;
        std::cout << "Isconverge2" << Isconverge2 << std::endl;
*/

        // LEDGER #168 (2026-09-26): the MIRRORED q lies outside the attainable
        // range whenever the score is a large fraction of the carriers' mass
        // and the adjusted genotype's negative mass is small -- exactly the
        // tiny-mu cells (mu(1-mu)-weighted centring gives them a negative mass
        // of ~1e-2 against a mirrored q of ~-1). That tail has probability
        // EXACTLY 0 (S cannot go below gneg), which is how SPAtest treats it;
        // LEDGER #22's rule made it a non-convergence and the whole test fell
        // back. Now: root found on q, mirror unattainable => p = the q tail.
        // An unattainable q itself (every carrier a case: S at its maximum) is
        // still a failure, since that atom is the exact test's to compute.
        const bool mirrorUnattainable = (!outuni2.Isconverge && outuni2.niter == 0 && std::isinf(outuni2.root));
        if(outuni1.Isconverge && mirrorUnattainable){
                SPASaddle getSaddle1 = Get_Saddle_Prob_Binom_core(outuni1.root, mu, g, q, logp, ws);
                ws.diagQSaddle = getSaddle1.isSaddle ? 1 : 0;
                if(getSaddle1.isSaddle && std::isfinite(getSaddle1.pval)){   // LEDGER #181: a non-finite saddle is not a saddle
                        pval = logp ? getSaddle1.pval : std::abs(getSaddle1.pval);
                        isSPAConverge = true;
                        return;
                }
        }
        // LEDGER #168 (b): the mirrored q can be attainable by a hair -- then its
        // root exists but lies where K1 is within that hair of its limit, at
        // |t| of 1e6 and beyond, past any finder. That tail needs no root:
        // Chernoff gives P(S beyond qinv) <= exp(K(t) - t qinv) at every t on
        // the tail's side, and the bound at the far end of the bracket is
        // astronomically small. When it is below a thousandth of the q-side
        // tail, the two-sided p IS the q-side tail to that precision.
        if(outuni1.Isconverge && !outuni2.Isconverge && !mirrorUnattainable){
                SPASaddle sq = Get_Saddle_Prob_Binom_core(outuni1.root, mu, g, q, logp, ws);
                ws.diagQSaddle = sq.isSaddle ? 1 : 0;
                if(sq.isSaddle && std::isfinite(sq.pval)){   // LEDGER #181
                        const double p1 = logp ? std::exp(sq.pval) : std::abs(sq.pval);
                        double gmax = 0.0;
                        for(arma::uword e = 0; e < g.n_elem; e++){ gmax = std::max(gmax, std::fabs(g(e))); }
                        const double side = (qinv < q) ? -1.0 : 1.0;   // lower tail: t < 0; upper: t > 0
                        double bound = std::numeric_limits<double>::infinity();
                        for(int k = 0; k <= 20 && gmax > 0.0; k++){
                                const double tt = side * std::ldexp(1.0, k) / gmax;
                                const double b = std::exp(Korg_Binom_ws(tt, mu, g, ws) - tt * qinv);
                                if(std::isfinite(b) && b < bound) bound = b;
                        }
                        ws.diagBound = bound;
                        if(bound <= 1e-3 * p1){
                                pval = logp ? sq.pval : p1;
                                isSPAConverge = true;
                                return;
                        }
                }
        }
        if(outuni1.Isconverge  && outuni2.Isconverge)
        {
                SPASaddle getSaddle  = Get_Saddle_Prob_Binom_core(outuni1.root, mu, g, q, logp, ws);
                SPASaddle getSaddle2 = Get_Saddle_Prob_Binom_core(outuni2.root, mu, g, qinv, logp, ws);


                if(getSaddle.isSaddle){
                        p1 = getSaddle.pval;
        		//std::cout << "p1 " << p1 << std::endl;

                }else{
			Isconverge = false; ws.failReason = 3;
                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
                if(getSaddle2.isSaddle){
                        p2 = getSaddle2.pval;
        		//std::cout << "p2 " << p2 << std::endl;
                }else{
			Isconverge = false; ws.failReason = 4;
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
                //Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        // LEDGER #169: a q-side saddle refused in the #168 branches above is
                        // a saddle failure, not a missing mirrored root -- the label that hid
                        // the CGF overflow on the BenIntNeo cells
                        Isconverge=false; ws.failReason = !outuni1.Isconverge ? 1 : (ws.diagQSaddle == 0 ? 3 : 2);
                }
        // LEDGER #181: the both-roots sum can be non-finite too (a NaN saddle value
        // with isSaddle true); it is a non-convergence, never a p
        if(!std::isfinite(pval)){ pval = pval_noadj; Isconverge = false; ws.failReason = 2; }
        isSPAConverge = Isconverge;
	//std::cout << "pval " << pval << std::endl;
	//result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}


// [[Rcpp::export]]
void SPA_fast(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, bool logp, const arma::vec & gNA, const arma::vec & gNB, const arma::vec & muNA, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge){
        SPAWorkspace ws(gNB.n_elem);
        // mu/gNA/muNA are accepted for the exported signature and not forwarded (C4)
        SPA_fast_ws(g, q, qinv, pval_noadj, logp, gNB, muNB, NAmu, NAsigma, tol, traitType, pval, isSPAConverge, ws);
}

void SPA_fast_ws(const arma::vec & g, double q, double qinv, double pval_noadj, bool logp, const arma::vec & gNB, const arma::vec & muNB,  double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge, SPAWorkspace & ws){

        double p1, p2;
        bool Isconverge = true;
        ws.failReason = 0;
        if( traitType != "binary"){
                // Unreachable, and deliberately not an Rcpp::stop. Only the
                // binary saddlepoint kernels exist, and SAIGE_test.cpp gates the
                // whole SPA branch on traitType != "quantitative", which since
                // survival was removed in Phase 3a leaves "binary" as the only
                // value that can arrive here. Before S1 a non-binary trait fell
                // through to reading "Isconverge" out of a default-constructed
                // (empty) Rcpp::List, which throws from Rcpp's name lookup; this
                // keeps it a throw, but one that enters no R API and is safe to
                // reach from a worker thread.
                throw std::logic_error("SPA_fast(): only the binary saddlepoint kernels exist");
        }
        // gpos/gneg depend only on g, and both roots use the same g
        double gpos, gneg;
        getGPosNeg_Binom(g, gpos, gneg);
        SPARoot outuni1 = getroot_K1_fast_Binom_core(0, g, q, gNB, muNB,NAmu, NAsigma, tol, gpos, gneg, ws);
        //double qinv = -1*q;
        SPARoot outuni2 = getroot_K1_fast_Binom_core(0, g, qinv, gNB, muNB,NAmu, NAsigma, tol, gpos, gneg, ws);
        ws.diagRoot1 = outuni1.root; ws.diagNiter1 = outuni1.niter; ws.diagRoot2 = outuni2.root; ws.diagNiter2 = outuni2.niter; ws.diagBound = -1; ws.diagQSaddle = -1;

/*	double outuni1root = outuni1["root"];
        double outuni2root = outuni2["root"];
        bool Isconverge1 = outuni1["Isconverge"];
        bool Isconverge2 = outuni2["Isconverge"];

        std::cout << "outuni1root" << outuni1root << std::endl;
        std::cout << "outuni2root" << outuni2root << std::endl;
        std::cout << "Isconverge1" << Isconverge1 << std::endl;
        std::cout << "Isconverge2" << Isconverge2 << std::endl;
*/

        // LEDGER #168 (2026-09-26): the MIRRORED q lies outside the attainable
        // range whenever the score is a large fraction of the carriers' mass
        // and the adjusted genotype's negative mass is small -- exactly the
        // tiny-mu cells (mu(1-mu)-weighted centring gives them a negative mass
        // of ~1e-2 against a mirrored q of ~-1). That tail has probability
        // EXACTLY 0 (S cannot go below gneg), which is how SPAtest treats it;
        // LEDGER #22's rule made it a non-convergence and the whole test fell
        // back. Now: root found on q, mirror unattainable => p = the q tail.
        // An unattainable q itself (every carrier a case: S at its maximum) is
        // still a failure, since that atom is the exact test's to compute.
        const bool mirrorUnattainable = (!outuni2.Isconverge && outuni2.niter == 0 && std::isinf(outuni2.root));
        if(outuni1.Isconverge && mirrorUnattainable){
                SPASaddle getSaddle1 = Get_Saddle_Prob_fast_Binom_core(outuni1.root, g, q, gNB, muNB, NAmu, NAsigma, logp, ws);
                ws.diagQSaddle = getSaddle1.isSaddle ? 1 : 0;
                if(getSaddle1.isSaddle && std::isfinite(getSaddle1.pval)){   // LEDGER #181: a non-finite saddle is not a saddle
                        pval = logp ? getSaddle1.pval : std::abs(getSaddle1.pval);
                        isSPAConverge = true;
                        return;
                }
        }
        // LEDGER #168 (b): the mirrored q can be attainable by a hair -- then its
        // root exists but lies where K1 is within that hair of its limit, at
        // |t| of 1e6 and beyond, past any finder. That tail needs no root:
        // Chernoff gives P(S beyond qinv) <= exp(K(t) - t qinv) at every t on
        // the tail's side, and the bound at the far end of the bracket is
        // astronomically small. When it is below a thousandth of the q-side
        // tail, the two-sided p IS the q-side tail to that precision.
        if(outuni1.Isconverge && !outuni2.Isconverge && !mirrorUnattainable){
                SPASaddle sq = Get_Saddle_Prob_fast_Binom_core(outuni1.root, g, q, gNB, muNB, NAmu, NAsigma, logp, ws);
                ws.diagQSaddle = sq.isSaddle ? 1 : 0;
                if(sq.isSaddle && std::isfinite(sq.pval)){   // LEDGER #181
                        const double p1 = logp ? std::exp(sq.pval) : std::abs(sq.pval);
                        double gmax = 0.0;
                        for(arma::uword e = 0; e < g.n_elem; e++){ gmax = std::max(gmax, std::fabs(g(e))); }
                        const double side = (qinv < q) ? -1.0 : 1.0;   // lower tail: t < 0; upper: t > 0
                        double bound = std::numeric_limits<double>::infinity();
                        for(int k = 0; k <= 20 && gmax > 0.0; k++){
                                const double tt = side * std::ldexp(1.0, k) / gmax;
                                const double b = std::exp(Korg_fast_Binom_ws(tt, gNB, muNB, NAmu, NAsigma, ws) - tt * qinv);
                                if(std::isfinite(b) && b < bound) bound = b;
                        }
                        ws.diagBound = bound;
                        if(bound <= 1e-3 * p1){
                                pval = logp ? sq.pval : p1;
                                isSPAConverge = true;
                                return;
                        }
                }
        }
        if(outuni1.Isconverge  && outuni2.Isconverge)
        {
                SPASaddle getSaddle  = Get_Saddle_Prob_fast_Binom_core(outuni1.root, g, q, gNB, muNB,NAmu, NAsigma, logp, ws);
                SPASaddle getSaddle2 = Get_Saddle_Prob_fast_Binom_core(outuni2.root, g, qinv, gNB, muNB,NAmu, NAsigma, logp, ws);
                if(getSaddle.isSaddle){
                        p1 = getSaddle.pval;
                }else{
			Isconverge = false; ws.failReason = 3;

                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }

                if(getSaddle2.isSaddle){
                        p2 = getSaddle2.pval;

                }else{
			Isconverge = false; ws.failReason = 4;
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }

                if(logp){
                        pval = add_logp(p1,p2);
                }else {
                        pval = std::abs(p1)+std::abs(p2);
                }
                //Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        // LEDGER #169: a q-side saddle refused in the #168 branches above is
                        // a saddle failure, not a missing mirrored root -- the label that hid
                        // the CGF overflow on the BenIntNeo cells
                        Isconverge=false; ws.failReason = !outuni1.Isconverge ? 1 : (ws.diagQSaddle == 0 ? 3 : 2);
        }
        // LEDGER #181: the both-roots sum can be non-finite too (a NaN saddle value
        // with isSaddle true); it is a non-convergence, never a p
        if(!std::isfinite(pval)){ pval = pval_noadj; Isconverge = false; ws.failReason = 2; }
	isSPAConverge = Isconverge;
        //result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}




// [[Rcpp::export]]
double SPA_pval(const arma::vec & mu, const arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, bool & isSPAConverge){
       double pval;
        SPAWorkspace ws(mu.n_elem);
        double p1, p2;
        bool Isconverge = true;
        if( traitType != "binary"){
                // Unreachable, and deliberately not an Rcpp::stop. Only the
                // binary saddlepoint kernels exist, and SAIGE_test.cpp gates the
                // whole SPA branch on traitType != "quantitative", which since
                // survival was removed in Phase 3a leaves "binary" as the only
                // value that can arrive here. Before S1 a non-binary trait fell
                // through to reading "Isconverge" out of a default-constructed
                // (empty) Rcpp::List, which throws from Rcpp's name lookup; this
                // keeps it a throw, but one that enters no R API and is safe to
                // reach from a worker thread.
                throw std::logic_error("SPA_pval(): only the binary saddlepoint kernels exist");
        }
        // gpos/gneg depend only on g, and both roots use the same g
        double gpos, gneg;
        getGPosNeg_Binom(g, gpos, gneg);
        SPARoot outuni1 = getroot_K1_Binom_core(0, mu, g, q, tol, gpos, gneg, ws);
        SPARoot outuni2 = getroot_K1_Binom_core(0, mu, g, qinv, tol, gpos, gneg, ws);
        ws.diagRoot1 = outuni1.root; ws.diagNiter1 = outuni1.niter; ws.diagRoot2 = outuni2.root; ws.diagNiter2 = outuni2.niter; ws.diagBound = -1; ws.diagQSaddle = -1;


        double outuni1root = outuni1.root;
        double outuni2root = outuni2.root;
        bool Isconverge1 = outuni1.Isconverge;
        bool Isconverge2 = outuni2.Isconverge;

/*
        std::cout << "outuni1root" << outuni1root << std::endl;
        std::cout << "outuni2root" << outuni2root << std::endl;
        std::cout << "Isconverge1" << Isconverge1 << std::endl;
        std::cout << "Isconverge2" << Isconverge2 << std::endl;
*/

        // LEDGER #168 (2026-09-26): the MIRRORED q lies outside the attainable
        // range whenever the score is a large fraction of the carriers' mass
        // and the adjusted genotype's negative mass is small -- exactly the
        // tiny-mu cells (mu(1-mu)-weighted centring gives them a negative mass
        // of ~1e-2 against a mirrored q of ~-1). That tail has probability
        // EXACTLY 0 (S cannot go below gneg), which is how SPAtest treats it;
        // LEDGER #22's rule made it a non-convergence and the whole test fell
        // back. Now: root found on q, mirror unattainable => p = the q tail.
        // An unattainable q itself (every carrier a case: S at its maximum) is
        // still a failure, since that atom is the exact test's to compute.
        const bool mirrorUnattainable = (!outuni2.Isconverge && outuni2.niter == 0 && std::isinf(outuni2.root));
        if(outuni1.Isconverge && mirrorUnattainable){
                SPASaddle getSaddle1 = Get_Saddle_Prob_Binom_core(outuni1.root, mu, g, q, logp, ws);
                ws.diagQSaddle = getSaddle1.isSaddle ? 1 : 0;
                if(getSaddle1.isSaddle && std::isfinite(getSaddle1.pval)){   // LEDGER #181: a non-finite saddle is not a saddle
                        pval = logp ? getSaddle1.pval : std::abs(getSaddle1.pval);
                        isSPAConverge = true;
                        return pval;
                }
        }
        // LEDGER #168 (b): the mirrored q can be attainable by a hair -- then its
        // root exists but lies where K1 is within that hair of its limit, at
        // |t| of 1e6 and beyond, past any finder. That tail needs no root:
        // Chernoff gives P(S beyond qinv) <= exp(K(t) - t qinv) at every t on
        // the tail's side, and the bound at the far end of the bracket is
        // astronomically small. When it is below a thousandth of the q-side
        // tail, the two-sided p IS the q-side tail to that precision.
        if(outuni1.Isconverge && !outuni2.Isconverge && !mirrorUnattainable){
                SPASaddle sq = Get_Saddle_Prob_Binom_core(outuni1.root, mu, g, q, logp, ws);
                ws.diagQSaddle = sq.isSaddle ? 1 : 0;
                if(sq.isSaddle && std::isfinite(sq.pval)){   // LEDGER #181
                        const double p1 = logp ? std::exp(sq.pval) : std::abs(sq.pval);
                        double gmax = 0.0;
                        for(arma::uword e = 0; e < g.n_elem; e++){ gmax = std::max(gmax, std::fabs(g(e))); }
                        const double side = (qinv < q) ? -1.0 : 1.0;   // lower tail: t < 0; upper: t > 0
                        double bound = std::numeric_limits<double>::infinity();
                        for(int k = 0; k <= 20 && gmax > 0.0; k++){
                                const double tt = side * std::ldexp(1.0, k) / gmax;
                                const double b = std::exp(Korg_Binom_ws(tt, mu, g, ws) - tt * qinv);
                                if(std::isfinite(b) && b < bound) bound = b;
                        }
                        ws.diagBound = bound;
                        if(bound <= 1e-3 * p1){
                                pval = logp ? sq.pval : p1;
                                isSPAConverge = true;
                                return pval;
                        }
                }
        }
        if(outuni1.Isconverge  && outuni2.Isconverge)
        {
                SPASaddle getSaddle  = Get_Saddle_Prob_Binom_core(outuni1.root, mu, g, q, logp, ws);
                SPASaddle getSaddle2 = Get_Saddle_Prob_Binom_core(outuni2.root, mu, g, qinv, logp, ws);

                if(getSaddle.isSaddle){
                        p1 = getSaddle.pval;
                }else{
                        Isconverge = false;
                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
                if(getSaddle2.isSaddle){
                        p2 = getSaddle2.pval;
                }else{
                        Isconverge = false;
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }
		//std::cout << "p1 " << p1 << " p2 " << p2 << std::endl;
                if(logp)
                {
                        pval = add_logp(p1,p2);
                } else {
                        pval = std::abs(p1)+std::abs(p2);
                }
		//std::cout << "pval " << pval << std::endl;
                //Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        Isconverge=false;
                }
        // LEDGER #181: the both-roots sum can be non-finite too (a NaN saddle value
        // with isSaddle true); it is a non-convergence, never a p
        if(!std::isfinite(pval)){ pval = pval_noadj; Isconverge = false; }
        isSPAConverge = Isconverge;

	return(pval);
        //result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}

