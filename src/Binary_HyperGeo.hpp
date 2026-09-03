#ifndef HyperGeo_001 
#define HyperGeo_001 



#include <vector>



using namespace std;


class HyperGeo {
    
public: 
    HyperGeo(){};
    ~HyperGeo();
    
    int     Run(int k, int ngroup, int ncase, int * group, double * weight);    
    int     Get_lprob(double * prob);
    int     Print();
    double  lCombinations(int n, int k);
    
protected:
    double  GetLogProb(int idx, int i);
    int     SaveProb(double lprob, int ncase_used);
    int     Recursive(double lprob, int idx, int ncase_used);

    int m_ngroup;
    int m_ncase;
    vector<int> m_group;   
    vector<double> m_lweight; 
    vector<double> m_kprob;
    int m_k;
    
    vector<double * > m_probtbl;
    
    double  m_ref;
   
    
};



// Bit-identical stand-in for Rf_lchoose that does not enter R's recursive
// choose.c entry point, and so is safe on a worker thread. Exposed only so
// tests/testthat/test-unit-lchoose.R can pin it against Rf_lchoose: the goldens
// execute it but are insensitive to small perturbations of it (a +1e-9 plant
// left the whole standard tier green), so they cannot serve as its regression
// test. Defined in Binary_HyperGeo.cpp.
double lchoose_threadsafe(double n, double k);

#endif 

