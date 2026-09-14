#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "Binary_ComputeExact.hpp"
#include "Binary_HyperGeo.hpp"
#include "Binary_global.hpp"

#ifdef _STAND_ALONE_
    #define Rprintf printf
#else
    #include <R.h>
    #include <Rmath.h>
#endif

/***************************************************
 
 Global functions
 
 ***************************************************/

void * SL_calloc(size_t num, size_t size){
	
	void * re =calloc(num, size);
	if(re == NULL){
		/*error("memory allocation error!");*/
        Rprintf("memory allocation error!");
        return(re) ;
	}
	
	return(re);
}


void SL_free(void * ptr){
	
	if(ptr != NULL){
		free(ptr);
        ptr=NULL;
	}
    
}

void SL_Sample(int k, int n, vector<int> & y, std::vector<int> & x){
    
    int i, j;
    
    for (i = 0; i < n; i++){
        x[i] = i;
    }
    for (i = 0; i < k; i++) {
        
        //temp = rand();
        //j = n * temp / (RAND_MAX+1);
        //j = temp % n;
        j = SL_runif_INT(n);
        y[i] = x[j] ;
        x[j] = x[--n];
    }
}


/*********************************************************
    Interface
 *********************************************************/

void SKAT_Exact(std::vector<int> & resarray, int nres, int * nres_k, double * Z0, double *Z1, int k, int m, int total, int * total_k, double *prob_k, double * odds, double * p1, int * IsExact, double * pval, double *pval_same, double *minP, int test_type, double epsilon)
{
    
    class ComputeExact exact;
    
    exact.Init(resarray, nres, nres_k, Z0, Z1, k, m, total, total_k, prob_k, odds, p1, IsExact, epsilon);
    exact.Run(test_type);
    //exact.PrintPval();
    exact.GetPvalues(pval, pval_same, prob_k, minP);
    
    
}


// Deleted 2026-09-14 (dead-surface sweep C1; do not re-import on an
// upstream merge): SKATO_Exact and the ComputeExactSKATO class it alone
// instantiated, SKAT_Permu and the Binary_Permu_SKAT class, GetProb (a
// verbatim twin of the live GetProb_new in ER_binary_func.cpp), the
// commented-out SKAT_ExactMC/RSKAT* wrappers, and Binary_ComputeExactMC.{cpp,h}
// -- git-tracked, #included here, never in Makevars' OBJECTS, and including a
// "Binary_global.h" that did not exist. SKAT_Exact above is the only entry
// point the shipped .so reaches (ER_binary_func.cpp).
