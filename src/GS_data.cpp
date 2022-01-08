#include "GS_data.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include "GSL_wrappers.h"

GS_data::GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin, const sample::GSL_RNG& gs_engine) {

    iterations=burnin + n_iter * thin;
    K=1;//da cambiare
    //Mstar inizializzata dopo
    M=4;//da cambiare
    lambda=2;//fixed

    for (unsigned int j = 0; j < dat.rows(); ++j) {
        std::vector<double> v;
        for (unsigned int i = 0; i <dat.cols() ; ++i) {
            if(!std::isnan(dat(j,i))){
                v.push_back(dat(j,i));
            }
        }
        data.push_back(v);
    }
    d=dat.rows();
    //std::cout<<d<<std::endl;
    for (unsigned int j = 0; j < dat.rows(); ++j) {
        for (unsigned int i = 0; i <dat.cols() ; ++i) {
            if(std::isnan(dat(j,i))){
                n_j.push_back(i);
                break;
            }
            if(i==dat.cols()-1){n_j.push_back(i);}
        }
    //std::cout<<n_j[j]<<std::endl;

    }

    std::vector<double> gamma(d,1.0);
    initialize_S(M, gs_engine);
    initialize_tau(M);
}

void GS_data::update_log_sum(){
    log_sum = 0.0;
    for(size_t j=0; j<d; j++){
        log_sum += log(U[j]+1)*gamma[j];
    }
    // AL POSTO DEL FOR: log_sum = log( U.array() + 1).dot(gamma);
}

void GS_data::initialize_S(unsigned int M, const sample::GSL_RNG& gs_engine){
  sample::rgamma rgamma;
  S = GDFMM_Traits::MatRow(d, M);
  for(size_t j=0; j<d; j++){
    for(size_t m=0; m<M; m++){
      S(j,m)=5;
    //S(j,m)=rgamma(gs_engine,gamma[j],1); dovrebbe essere così ma non funziona
    }
  }
}
void GS_data::initialize_N(unsigned int K){
  N = GDFMM_Traits::MatUnsCol(d, K);
}
void GS_data::initialize_tau(unsigned int M){
  mu = std::vector<double>(M, 57.0);
  sigma = std::vector<double>(M, 40.0);
}
