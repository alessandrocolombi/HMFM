#include "GS_data.h"
#include <Rcpp.h>
#include <RcppEigen.h>



GS_data::GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin,
                const sample::GSL_RNG& gs_engine, unsigned int Mstar0, double Lambda0,
                double mu0, double nu0, double sigma0) {

    iterations = 0;
    K = 1; // Inizialmente tutte le osservazioni appartengono allo stesso gruppo
    Mstar = Mstar0; //Mstar inizializzata dopo
    lambda = Lambda0; //fixed

    M = K+Mstar;//M è somma di componenti Allocate e Non Allocate

    for (unsigned int j = 0; j < dat.rows(); ++j) {
        std::vector<double> v;
        for (unsigned int i = 0; i <dat.cols() ; ++i) {
            if(!std::isnan(dat(j,i))){
                v.push_back(dat(j,i));
            }
        }
        data.push_back(v);
    }
    d = dat.rows();
    Rcpp::Rcout << "Data read, d is : " << d << std::endl;

    // Initialization of n_j
    n_j = std::vector<unsigned int>(d, 0);    
    for (unsigned int j = 0; j < d; ++j) {
        for (unsigned int i = 0; i <dat.cols() ; ++i) {

            if(std::isnan(dat(j,i))){
                n_j[j] = i;
                break;

            }
            if(i == dat.cols()-1){n_j[j] = i+1;}
        }
        //Rcpp::Rcout<<n_j[j]<<std::endl;
    }
    Rcpp::Rcout << "n_j Initialized : " << n_j[0] << n_j[d-1]<< std::endl;
    // Initialization of partition data structures
    initialize_Partition(n_j);
    Rcpp::Rcout << "Partition Initialized "<< std::endl;
    // Initialization of gamma and U vector
    gamma = std::vector<double>(d, 1.0);
    Rcpp::Rcout << "gamma vector Initialized "<< std::endl;
    U = std::vector<double>(d, 0.0);
    Rcpp::Rcout << "U vector Initialized "<< std::endl;
    /*for (int l = 0; l <K ; ++l) {
        N_k.push_back(0);

    }*/
    //std::vector<double> U(d,0.0);
    // Random Initialization of S and tau form the prior
    initialize_S(M, gs_engine);
    Rcpp::Rcout << "S matrix Initialized "<< std::endl;
    initialize_tau(M, nu0, mu0, sigma0, gs_engine);
    Rcpp::Rcout << "tau Initialized "<< std::endl;
    //Rcpp::Rcout<< p<<std::endl;
}

void GS_data::update_log_sum(){
    log_sum = 0.0;
   for(size_t j=0; j<d; j++){
        //Rcpp::Rcout<<U[j]<<std::endl;
        //Rcpp::Rcout<<gamma[j]<<std::endl;
       log_sum += log(U[j]+1)*gamma[j];
    }
    // AL POSTO DEL FOR: log_sum = log( U.array() + 1).dot(gamma);
}

void GS_data::initialize_Partition(const std::vector<unsigned int>& n_j){
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
    S = GDFMM_Traits::MatRow(d, M);
    for (unsigned int j = 0; j <d ; ++j) {
        for (unsigned int m = 0; m <M ; ++m) {
            S(j,m) = 0;
        }
    }
}

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
        sigma[m] =  1 / Gamma(gs_engine, nu0, 2 / (nu0*sigma0)); // prendo
        mu[m] = rnorm(gs_engine, mu0, sqrt(sigma[m]));
  }
}

void GS_data::allocate_N(unsigned int K){
    N = GDFMM_Traits::MatUnsCol(d, K);
    for (unsigned int j = 0; j <d ; ++j) { // Andre : penso che il ciclo for sia totalmente inutile, perchè
      for (unsigned int k = 0; k<K ; ++k){ //         l'inizializzazione dovrebbe già avvenire a 0 --> VERIFICARE
        N(j,k) = 0;}
    }
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
