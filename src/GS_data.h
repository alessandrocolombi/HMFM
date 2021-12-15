#ifndef __GS_DATA_H
#define __GS_DATA_H

#include "include_headers.h"
#include "recurrent_traits.h"

struct GS_data{
    /* DATA */
    // single values
    unsigned int d; // number of groups
    unsigned int K; // number of allocated component ==> number of clusters
    unsigned int Mstar; // number of NON-allocated component
    unsigned int M; // total number of component
    double lambda; // M|lambda ~ Poi(lambda)
    double log_sum; // sum of log(U_j+1)*gamma_j : logarithm of 1/psi_prod
    // vectors
    std::vector<unsigned int> n_j; // number of elements for each group (dimension: d)
    std::vector<double> U; // auxiliary variable
    std::vector<double> gamma; // vector of d gamma, one for each group
    std::vector<double> mu; // vector of the mean for each component
    std::vector<double> sigma; // vector of the variance for each component
                               // N.B. sample::rnorm takes the s.d. as input ==> use sqrt(sigma[m])
    // matrix or vector of vectors
    GDFMM_Traits::MatRow S; // dxM matrix; allocated and NON-alloc together
    GDFMM_Traits::MatUnsCol N; // dxK matrix; only allocated components have n_jk>0
    //-----------------------------------------------------------//
    /* CONSTRUCTOR */
    /* METHODS */
    void initialize_S(unsigned int M);
    void initialize_tau(unsigned int M);
    void update_log_sum();
};

#endif