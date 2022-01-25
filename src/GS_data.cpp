#include "GS_data.h"
#include <Rcpp.h>
#include <RcppEigen.h>



GS_data::GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin,const sample::GSL_RNG& gs_engine) {

    iterations = 0;
    K = 1; // Inizialmente tutte le osservazioni appartengono allo stesso gruppo
    Mstar = 3;//Mstar inizializzata dopo
    M = 4;//da cambiare
    lambda = 2;//fixed

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
    //Rcpp::Rcout << "d is : " << d << std::endl;
    //std::cout<<d<<std::endl;

    for (unsigned int j = 0; j < d; ++j) {

        for (unsigned int i = 0; i <dat.cols() ; ++i) {

            if(std::isnan(dat(j,i))){
                n_j.push_back(i);
                break;

            }
            if(i == dat.cols()-1){n_j.push_back(i);}
        }
        //Rcpp::Rcout<<n_j[j]<<std::endl;
    }
    initialize_Ctilde(n_j);
    for (int k = 0; k < d; ++k) {
        gamma.push_back(1);
        U.push_back(0);
    }

    for (int l = 0; l <K ; ++l) {
        N_k.push_back(0);

    }

    //std::vector<double> U(d,0.0);


    initialize_S(M, gs_engine);
    initialize_tau(M, gs_engine);
    initialize_N(K);
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

void GS_data::initialize_Ctilde(const std::vector<unsigned int>& n_j){
    // Initialize Ctilde so that it assign every observation to the same cluster
    Ctilde.clear();
    for(size_t j = 0; j < d; j++){
        std::vector<unsigned int> row_j(n_j[j], 1);
        Ctilde.push_back(row_j);
    }
}

void GS_data::allocate_S(unsigned int M){
    S = GDFMM_Traits::MatRow(d, M);
    for (int i = 0; i <d ; ++i) {
        for (int j = 0; j <M ; ++j) {
            S(i,j)=0;
        }
    }
}

void GS_data::initialize_S(unsigned int M, const sample::GSL_RNG& gs_engine){
  //sample::rgamma rgamma;
  sample::rgamma Gamma;
  S = GDFMM_Traits::MatRow(d, M);
  for (int i = 0; i <d ; ++i) {
    for (int j = 0; j <M ; ++j) {
      S(i,j)=Gamma(gs_engine, gamma[j],1);;
    }
  }
}

void GS_data::initialize_N(unsigned int K){
    N = GDFMM_Traits::MatUnsCol(d, K);
    for (int i = 0; i <d ; ++i) {
      for (int j = 0; j<K ; ++j){
        N(i,j) = 0;}
    }
    N_k = std::vector<unsigned int>(K, 0);
}

void GS_data::initialize_tau(unsigned int M, const sample::GSL_RNG& gs_engine){
    mu = std::vector<double>(M, 0.0);
    sigma = std::vector<double>(M, 1.0);

    sample::rgamma Gamma;
    sample::rnorm rnorm;

    for( size_t m = 0; m < M; m++){
        sigma[m] =  1 / Gamma(gs_engine, 2.5, 2 / (2.5*40));
        mu[m] = rnorm(gs_engine, 60, sqrt(sigma[m]));
  }
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

