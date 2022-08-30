#include "GS_data.h"
#include <Rcpp.h>
#include <RcppEigen.h>

GS_data::GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin,
                const sample::GSL_RNG& gs_engine, unsigned int Mstar0, double Lambda0, double mu0,
                double nu0, double sigma0, double gamma0, std::string P0_prior_name, std::vector<unsigned int> part_vec, double _nu) :
                prior(P0_prior_name) {

    iterations = 0;
    lambda = Lambda0;
    Mstar = Mstar0;
    nu = _nu;

    // Read Data and extract d, n_j
    for (unsigned int j = 0; j < dat.rows(); ++j) {
        std::vector<double> v;
        for (unsigned int i = 0; i <dat.cols() ; ++i) {
            if(!std::isnan(dat(j,i))){ //na are not included. Not fine if data contains missing data
                v.push_back(dat(j,i));
            }
        }
        data.push_back(v);
    }
    d = dat.rows();
    // Rcpp::Rcout << "Data read, d is : " << d << std::endl;

    // Initialization of n_j
    n_j = std::vector<unsigned int>(d, 0);
    for (unsigned int j = 0; j < d; ++j) {
        for (unsigned int i = 0; i <dat.cols() ; ++i) {

            if(std::isnan(dat(j,i))){
                n_j[j] = i;
                break;
            }
            if(i == dat.cols()-1){
                n_j[j] = i+1;
            }
        }
    }
    // Rcpp::Rcout << "n_j Initialized : " << n_j[0] << " " << n_j[d-1]<< std::endl;
    
    // Initialization of partition data structures
    if( part_vec.empty() ){
        // Code should never arrive here
        K = 1;
        Mstar = Mstar0;
        M = K + Mstar;
        initialize_Partition();
        throw std::runtime_error("part_vec was found empty. This should no longer happen, why here? ");
    }
    else{
        initialize_Partition(part_vec);
    }
    
    // Initialization of gamma and U vector
    gamma = std::vector<double>(d, gamma0);
    // Rcpp::Rcout << "gamma vector Initialized "<< std::endl;
    U = std::vector<double>(d, 0.0);
    //Initialize log_sum
    update_log_sum();
    
    // Random Initialization of S and tau form the prior
    initialize_S(M, gs_engine); // inutile nel caso marginale ma non fa danni. NON va molto bene in ottica tener fisso S ad un valore iniziale!
    // Rcpp::Rcout << "S matrix Initialized "<< std::endl;
    
    initialize_tau(M, nu0, mu0, sigma0, gs_engine); // NON va molto bene in ottica tener fisso tau ad un valore iniziale!
    // Rcpp::Rcout << "tau Initialized "<< std::endl;
}


void GS_data::update_log_sum(){
    log_sum = 0.0;
    //double psi_S = 0.0; //just for testing
    //double psi_Sprime = 0.0; //just for testing
   for(size_t j=0; j<d; j++){
        //Rcpp::Rcout<<U[j]<<std::endl;
        //Rcpp::Rcout<<gamma[j]<<std::endl;
       //log_sum += log(U[j]+1.0)*gamma[j]; // this is the sum of log( psi_S(U) )
       //psi_S += log(U[j]+1.0)*gamma[j]; // this is the sum of log( psi_S(U) ) //just for testing
       
       log_sum += log(U[j]*1.0/nu + 1.0)*gamma[j]; // this is the sum of log( psi_S(U'/nu) )
       
       //log_sum += gamma[j]* ( log(nu) - log(U[j] + nu) ) ; // this is the sum of log( psi_S'(U') )
       //psi_Sprime += -gamma[j]* ( log(nu) - log(U[j] + nu) ) ; // this is the sum of log( psi_S'(U') ) //just for testing
    }
    // AL POSTO DEL FOR: log_sum = log( U.array() + 1).dot(gamma);
            //Rcpp::Rcout<<"old psi_S(U) = "<<psi_S<<std::endl;
            //Rcpp::Rcout<<"psi_S(U'/nu) = "<<log_sum<<std::endl;
            //Rcpp::Rcout<<"psi_S'(U') = "<< psi_Sprime <<std::endl;
}

// Initialize partition (Ctilde, N, N_k) when part_vect is not empty
void GS_data::initialize_Partition(const std::vector<unsigned int>& partition_vec){
    
    // Extract number of clusters from partition vec
    const auto max_it = std::max_element(partition_vec.cbegin(), partition_vec.cend()); //get the maximum 
    K =  *max_it + 1;
    M = K + Mstar;
    Rcpp::Rcout<<"initialize_Partition with non empty partition_vec"<<std::endl;
    Rcpp::Rcout<<"Watch out modification: Mstar is not set to zero but to Mstar0"<<std::endl;
    Rcpp::Rcout << " (K, Mstar, M) = ("<< K <<","<<Mstar<<","<<M<<")"<<std::endl;
    
    // Allocate Ctilde, N, N_k 
    Ctilde.clear();
    for(size_t j = 0; j < d; j++){
        std::vector<unsigned int> row_j(n_j[j], 0);
        Ctilde.push_back(row_j);
    }
    N = GDFMM_Traits::MatUnsCol(d, K);
    for (unsigned int j = 0; j <d ; ++j)
        for(unsigned int k = 0; k < K; ++k)
            N(j,k) = 0;
    N_k = std::vector<unsigned int>(K, 0);

    // Initialize Ctilde, N, N_k with info from partition_vec
    unsigned int j = 0;
    unsigned int i = 0;
    for(const unsigned int& clust_ji : partition_vec){
        Ctilde[j][i] = clust_ji;
        N(j, clust_ji)++;
        N_k[clust_ji]++;
        i++;
        // If we have just 
        if(i == n_j[j]){
            j++;
            i = 0;
        }
    }
}

// Initialize partition (Ctilde, N, N_k) when part_vect is empty
// This function is no longer use, isn't it?
void GS_data::initialize_Partition(){

    Rcpp::Rcout<<"initialize_Partition with EMPTY partition_vec"<<std::endl;
    // Initialize Ctilde so that it assign every observation to the same cluster
    Ctilde.clear();
    for(size_t j = 0; j < d; j++){
        std::vector<unsigned int> row_j(n_j[j], 1);
        Ctilde.push_back(row_j);
    }
    // Initialize Matrix N (recall that all observations are assigned to cluster 1)
    N = GDFMM_Traits::MatUnsCol(d, K);
    for (unsigned int j = 0; j <d ; ++j)
        N(j,0) = n_j[j];
    // Initialize N_k accordingly
    N_k = std::vector<unsigned int>(1, std::accumulate(n_j.begin(), n_j.end(), 0));
}

void GS_data::allocate_S(unsigned int M){
    //Rcpp::Rcout<<"Chiamo allocate_S"<<std::endl;
    S = GDFMM_Traits::MatRow::Constant(d, M, 1.0);

    // Old version
    //S = GDFMM_Traits::MatRow(d, M);
    //for (unsigned int j = 0; j <d ; ++j) {
        //for (unsigned int m = 0; m <M ; ++m) {
            //S(j,m) = 1;
        //}
    //}
}


// S influences the draw of the partition. Hence, this random initialization affects the first draw of the partition.
// Then it is sampled from the full conditional of set equal to 1.
void GS_data::initialize_S(unsigned int M, const sample::GSL_RNG& gs_engine){
  //sample::rgamma rgamma;
  sample::rgamma Gamma;
  S = GDFMM_Traits::MatRow(d, M);
  for (unsigned int j = 0; j <d ; ++j) {
    for (unsigned int m = 0; m <M ; ++m) {
      S(j,m)=Gamma(gs_engine, gamma[j],1);
    }
  }
}

void GS_data::initialize_tau(unsigned int M, double nu0, double mu0, double sigma0,
                             const sample::GSL_RNG& gs_engine){
    mu = std::vector<double>(M, 0.0);
    sigma = std::vector<double>(M, 1.0);

    sample::rgamma Gamma;
    sample::rnorm rnorm;

    for(unsigned m = 0; m < M; m++){
        sigma[m] =  1 / Gamma(gs_engine, nu0/2, 2 / (nu0*sigma0));
        mu[m] = rnorm(gs_engine, mu0, sqrt(sigma[m]));
  }
}

void GS_data::allocate_N(unsigned int K){
    //Rcpp::Rcout<<"Chiamo allocate_N"<<std::endl;
    N = GDFMM_Traits::MatUnsCol::Constant(d, K, 0);

    // Old Version
    //for (unsigned int j = 0; j <d ; ++j) { // Andre : penso che il ciclo for sia totalmente inutile, perchè
      //for (unsigned int k = 0; k<K ; ++k){ //         l'inizializzazione dovrebbe già avvenire a 0 --> VERIFICARE
        //N(j,k) = 0;}
    //}

    N_k = std::vector<unsigned int>(K, 0);
}

void GS_data::allocate_tau(unsigned int M){
    mu = std::vector<double>(M, 0.0);
    sigma = std::vector<double>(M, 1.0);
}

void GS_data::update_Ctilde(const std::vector< std::vector<unsigned int>> &C,
                            const std::vector<unsigned int> &clust_out) {
    // update of Ctilde, that imply the update of also N and N_k
    for (unsigned int m = 0; m < K; m++){
        for (unsigned int j = 0; j <d ; ++j) {
            for (unsigned int i = 0; i < n_j[j] ; ++i) {
                //Rcpp::Rcout<< C[j][i];
                //Rcpp::Rcout<< clust_out[m]<<std::endl;
                if(C[j][i] == clust_out[m]){
                    N(j,m)++;
                    Ctilde[j][i] = m;
                }
                // Rcpp::Rcout<< N(j,m);
            }
            N_k[m] = N_k[m]+N(j,m);
        }
    }
}
