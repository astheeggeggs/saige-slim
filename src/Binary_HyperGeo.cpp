#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "Binary_HyperGeo.hpp"
#include "Binary_global.hpp"

#ifdef _STAND_ALONE_
    #define Rprintf printf
#else
    #include <R.h>
    #include <Rmath.h>
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/***************************************************

 Functions
 
 ***************************************************/

HyperGeo::~HyperGeo(){
    
    
    double * ptemp;
    for(int i=0;i<m_probtbl.size();i++){
        
        ptemp = m_probtbl[i];
        SL_free(ptemp);
    }
                        
}


int     HyperGeo::Run(int k, int ngroup, int ncase, int * group, double * weight){

    int i;
    m_ngroup = ngroup;
    m_ncase=ncase;
    m_k =k;

    for(i=0;i< ngroup;i++){
        m_group.push_back(group[i]);
        m_lweight.push_back(log(weight[i]));
    }
    
    for(i=0;i<=m_k;i++){
        m_kprob.push_back(0);
    }
    
    // generate prob table
    double * temp, lch;
    m_ref = 0;
    for(i=0;i<ngroup;i++){
        
        if(i < ngroup -1){
            temp = (double *) SL_calloc(m_group[i] + 1, sizeof(double));
            for(int j=0;j<=m_group[i];j++){
                
                lch = lCombinations(m_group[i], j);
                temp[j] = lch + m_lweight[i]*j;
                
                //Rprintf("[%d][%d][%e][%e]\n",i,j,lch, temp[j]);
            }
        
        } else {
            temp = (double *) SL_calloc(k+1, sizeof(double));
            
            for(int j=0;j<=k;j++){
                
                lch = lCombinations(m_group[i], m_ncase -j);
                temp[j] = lch + m_lweight[i]* ((double) m_ncase -j);
                
                m_ref = MAX(temp[j], m_ref);
                
                //Rprintf("[%d][%d][%e][%e]\n",i,j,lch, temp[j]);
            }
        }
        
        m_probtbl.push_back(temp);
    }
    
    Recursive(0, 0, 0);
    
    
    return 1;
}

double HyperGeo::GetLogProb(int idx, int i){

    //Rprintf("logprob[%d][%d][%e]\n",idx,i, m_probtbl[idx][i]);
    
    return m_probtbl[idx][i];

}

int HyperGeo::SaveProb(double lprob, int ncase_used){
    
    m_kprob[ncase_used] += exp(lprob - m_ref);
    //Rprintf("[%i][%e][%e][%e]\n",ncase_used, lprob, exp(lprob-m_ref), m_ref);
    
    return 1;

}

int     HyperGeo::Recursive(double lprob, int idx, int ncase_used){
    
	int i;
    double lprob1 ;
    
    int num = m_group[idx];
    
    if(idx == m_ngroup -1){
        lprob1 = GetLogProb(idx, ncase_used);
        SaveProb(lprob + lprob1, ncase_used);
        return 1;
    } 
    
	
    for(i=0; i<=num; i++) {
       
        if(ncase_used +i <=m_ncase){
            lprob1 = GetLogProb(idx, i);
            Recursive(lprob + lprob1, idx+1, ncase_used + i);
        }
    }
        
    return 0;
}

int     HyperGeo::Get_lprob(double * prob){
    
    double sum1=0;
    for(int i=0;i<=m_k;i++){
        sum1 += m_kprob[i];
    }

    for(int i=0;i<=m_k;i++){
        prob[i] = m_kprob[i] / sum1;
    }
    
    return 1;
    
}



int     HyperGeo::Print(){
    
    double sum1=0;
    double prob;
    for(int i=0;i<=m_k;i++){
        sum1 += m_kprob[i];
        
    }
    
    for(int i=0;i<=m_k;i++){
        prob = m_kprob[i] / sum1;
        Rprintf("%d:[%e][%e]\n",i, prob, m_kprob[i]);
    }
    
    return 1;
    
}


#ifndef _STAND_ALONE_
// Rf_lchoose cannot be called off the main thread, and the reason turned out
// NOT to be the arithmetic: R's C-stack check measures depth against the MAIN
// thread's stack base, so a worker's stack pointer reads as impossibly deep and
// R aborts the session. choose()/lchoose() are the two functions here that
// recurse and so carry that check; Rf_lgammafn, Rf_pnorm5 and -- the one that
// matters -- Rf_lbeta are all fine on a worker. The decisive measurement was
// lchoose(1, 0), which returns 0 through an early return that makes no
// math-library call at all and STILL aborts.
//
// So compute it the way R computes it, without entering that entry point. For
// the domain lCombinations can reach -- integer n >= 0, integer k, with k > n
// already short-circuited by the caller -- R's choose.c reduces to the early
// returns below plus lfastchoose(n,k) = -log(n+1) - lbeta(n-k+1, k+1). Both of
// R's recursive arms unroll here: n < 0 cannot occur, and the symmetry arm
// lchoose(n, n-k) lands on n-k in {0,1}, which are the k < 2 early returns.
//
// VERIFIED BIT FOR BIT, not by inspection: over n in [0,2000] and k in [-3,n],
// 2,009,004 argument pairs, ZERO differences from Rf_lchoose under memcmp (so
// NaN/Inf/-0.0 are held to the same standard) -- and the same 2,009,004 values
// computed on a worker thread with R's stack check LEFT ON, again with zero
// differences. harness/probes/lchoose_replacement.cpp reproduces both.
// Because it is bit-identical, no golden can move.
double lchoose_threadsafe(double n, double k)
{
    if (k < 2) {
        if (k < 0) return R_NegInf;
        if (k == 0) return 0.;
        return std::log(std::fabs(n));            // k == 1
    }
    if (n < k) return R_NegInf;
    if (n - k < 2) {                              // symmetry arm, n-k is 0 or 1
        if (n - k == 0) return 0.;
        return std::log(std::fabs(n));
    }
    return -std::log(n + 1.) - Rf_lbeta(n - k + 1., k + 1.);
}
#endif

double HyperGeo::lCombinations(int n, int k)
{
    if (k > n)
        return 0;
    double r = 0;
    
    #ifdef _STAND_ALONE_
    
    for (int d = 1; d <= k; ++d)
    {
        r += log(n--);
        r -= log(d);
    }
    #else
        // No main-thread guard here any more: lchoose_threadsafe is bit-identical
        // to Rf_lchoose and safe on a worker. The ER path's remaining R
        // dependency is SL_setseed, guarded where it lives, on the resampling
        // arms only (Binary_ComputeExact.cpp).
        r = lchoose_threadsafe (n, k);
    #endif
    
    return r;
}






