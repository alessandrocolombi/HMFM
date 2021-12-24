#ifndef __GS_DATA_H
#define __GS_DATA_H

#include "include_headers.h"
#include "recurrent_traits.h"
#include "Partition.h"

struct GS_data{
    /* DATA */
    // single values
    unsigned int d;
    unsigned int iterations;// number of groups
    unsigned int K; // number of allocated component ==> number of clusters
    unsigned int Mstar; // number of NON-allocated component
    unsigned int M; // total number of component
    double lambda; // M|lambda ~ Poi(lambda)
    double log_sum; // sum of log(U_j+1)*gamma_j : logarithm of 1/psi_prod
    // vectors
    std::vector<std::vector<double>> data; // our data, y_ji
    std::vector< std::vector<unsigned int>> Ctilde; //output partition ANDRE: DIMENSIONI?
    std::vector<unsigned int> n_j; // number of elements in  each group (dimension: d)
    std::vector<unsigned int> N_k; // number of elements in each cluster
    std::vector<double> U; // auxiliary variable
    std::vector<double> gamma; // vector of d gamma, one for each group
    std::vector<double> mu; // vector of the mean for each component
    std::vector<double> sigma; // vector of the variance for each component
                               // N.B. sample::rnorm takes the s.d. as input ==> use sqrt(sigma[m])
    std::string prior; //Which prior are we using for tau? noga or normal inverse gamma?
    // matrix or vector of vectors
    GDFMM_Traits::MatRow S; // dxM matrix; allocated and NON-alloc together
    GDFMM_Traits::MatUnsCol N; // dxK matrix; only allocated components have n_jk>0
    // Partition
    Partition p;
    //-----------------------------------------------------------//
    /* CONSTRUCTOR */
    /* METHODS */
    void initialize_S(unsigned int M);
    void initialize_tau(unsigned int M);
    void initialize_tau(unsigned int K);
    void update_log_sum();
};

#endif
