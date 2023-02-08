#ifndef __GS_DATA_H
#define __GS_DATA_H

#include "include_headers.h"
#include "recurrent_traits.h"
#include "GSL_wrappers.h"
#include "Individual.h"

struct GS_data{
    /* DATA */

    // single values
    unsigned int d; // number of levels
    unsigned int iterations;
    unsigned int K; // number of allocated component ==> number of clusters
    unsigned int Mstar; // number of NON-allocated component
    unsigned int M; // total number of component
    double lambda; // M|lambda ~ Poi(lambda)
    //double nu; // S_jm ~ gamma(gamma_j, nu), da togliere
    double log_sum; // sum of log(U_j+1)*gamma_j : logarithm of 1/psi_prod
    
    // vectors
    std::vector<std::vector<double>> data; // our data, y_ji
    std::vector<std::vector<Individual>> mv_data; // our multivariate data, now each y_ji is a generic individual whose values are defined in a specific class
    std::vector< std::vector<unsigned int>> Ctilde; //output partition
    std::vector<unsigned int> n_j; // number of elements in  each group (dimension: d)
    std::vector<unsigned int> N_k; // number of elements in each cluster
    std::vector<double> U; // auxiliary variable
    std::vector<double> gamma; // vector of d gamma, one for each group
    std::vector<double> mu; // vector of the mean for each component
    std::vector<double> sigma; // vector of the variance for each component // N.B. sample::rnorm takes the s.d. as input ==> use sqrt(sigma[m])
    std::vector<double> sum_cluster_elements; // vector of length K, each element contains the sum of the data within that cluster
    std::vector<double> squared_sum_cluster_elements; // vector of length K, each element contains the squared sum of the data within that cluster
    std::vector<std::vector<double>> log_prob_marginal_data; //same structure of data, it contains the marginal probabilities of each data in log-scale.

    // strings
    std::string prior; //Name of the prior for tau
    
    // matrix or vector of vectors
    GDFMM_Traits::MatRow S; // dxM matrix; allocated and NON-alloc together
    GDFMM_Traits::MatUnsCol N; // dxK matrix; only allocated components have n_jk>0

    //-----------------------------------------------------------//
    
    /* CONSTRUCTORS */

    // Constructor with default prior (Normal-InvGamma)
    GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin,
            const sample::GSL_RNG& gs_engine, unsigned int _Mstar0, double _Lambda0, double _mu0,
            double nu0, double sigma0, double _gamma0, 
            std::vector<double> _init_mean_clus, std::vector<double> _init_var_clus, 
            std::vector<unsigned int> _part_vec) : 
                        GS_data(dat, n_iter, burnin, thin, gs_engine, _Mstar0, _Lambda0, _mu0, nu0,
                                sigma0, _gamma0, _init_mean_clus, _init_var_clus, "Normal-InvGamma", _part_vec){}
    
    // Constructor with user defined prior
    GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin,
                const sample::GSL_RNG& gs_engine, unsigned int _Mstar0, double _Lambda0, double _mu0,
                double nu0, double sigma0, double _gamma0, 
                std::vector<double> _init_mean_clus, std::vector<double> _init_var_clus, 
                std::string P0_prior_name, std::vector<unsigned int> _part_vec);

    // Constructor for multivariate data
    GS_data(    const std::vector<std::vector<Individual>>& _dat, 
                const std::vector<unsigned int>& _n_j, const unsigned int _d,
                const sample::GSL_RNG& gs_engine, 
                unsigned int _Mstar0, double _Lambda0, double _mu0,
                double _nu0, double _sigma0, double _gamma0, 
                const std::vector<double>& _init_mean_clus, const std::vector<double>& _init_var_clus, 
                std::string P0_prior_name, const std::vector<unsigned int>& _part_vec);

    GS_data(){};
    ~GS_data(){};

    /* METHODS */
    // Initialize partition (Ctilde, N, N_k) when it is FIXED
    void initialize_Partition(const std::vector<unsigned int>& partition_vec);
    // Initialize partition (Ctilde, N, N_k) when M and K are RANDOM
    void initialize_Partition();
    void initialize_S(unsigned int M,  const sample::GSL_RNG& gs_engine);
    void allocate_S(unsigned int M);
    void initialize_tau(unsigned int M, const std::vector<double>& init_mean_clus, const std::vector<double>& init_var_clus,
                        double nu0, double mu0, double sigma0,
                        const sample::GSL_RNG& gs_engine);
    // update methods;
    void allocate_N(unsigned int K);
    void update_Ctilde(const std::vector< std::vector<unsigned int>> &C,
                            const std::vector<unsigned int> &clust_out);
    void update_log_sum();
    void allocate_tau(unsigned int M);

    // methods regarding the vectors to compute mean and variance in clusters
    void initialize_sums_in_clusters();
    double compute_var_in_cluster(const unsigned int& m) const;

    // compute log of marginal probabilities of each data, i.e log(m(x))
    void compute_log_prob_marginal_data(double nu_0, double sigma_0, double mu_0, double k_0);
};

#endif
