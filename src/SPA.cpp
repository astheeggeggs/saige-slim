// [[Rcpp::depends(BH)]]

//[[Rcpp::depends(RcppArmadillo)]]
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

        if(outuni1.Isconverge  && outuni2.Isconverge)
        {
                SPASaddle getSaddle  = Get_Saddle_Prob_Binom_core(outuni1.root, mu, g, q, logp, ws);
                SPASaddle getSaddle2 = Get_Saddle_Prob_Binom_core(outuni2.root, mu, g, qinv, logp, ws);


                if(getSaddle.isSaddle){
                        p1 = getSaddle.pval;
        		//std::cout << "p1 " << p1 << std::endl;

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
        		//std::cout << "p2 " << p2 << std::endl;
                }else{
			Isconverge = false;
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
                        Isconverge=false;
                }
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

/*	double outuni1root = outuni1["root"];
        double outuni2root = outuni2["root"];
        bool Isconverge1 = outuni1["Isconverge"];
        bool Isconverge2 = outuni2["Isconverge"];

        std::cout << "outuni1root" << outuni1root << std::endl;
        std::cout << "outuni2root" << outuni2root << std::endl;
        std::cout << "Isconverge1" << Isconverge1 << std::endl;
        std::cout << "Isconverge2" << Isconverge2 << std::endl;
*/

        if(outuni1.Isconverge  && outuni2.Isconverge)
        {
                SPASaddle getSaddle  = Get_Saddle_Prob_fast_Binom_core(outuni1.root, g, q, gNB, muNB,NAmu, NAsigma, logp, ws);
                SPASaddle getSaddle2 = Get_Saddle_Prob_fast_Binom_core(outuni2.root, g, qinv, gNB, muNB,NAmu, NAsigma, logp, ws);
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

                if(logp){
                        pval = add_logp(p1,p2);
                }else {
                        pval = std::abs(p1)+std::abs(p2);
                }
                //Isconverge=true;
        }else {
                        //std::cout << "Error_Converge" << std::endl;
                        pval = pval_noadj;
                        Isconverge=false;
        }
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
        isSPAConverge = Isconverge;

	return(pval);
        //result["pvalue"] = pval;
        //result["Isconverge"] = Isconverge;
        //return(result);
}

