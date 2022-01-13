#include "GS_data.h"
#include <Rcpp.h>
#include <RcppEigen.h>


GS_data::GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin,const sample::GSL_RNG& gs_engine) {

    iterations=burnin + n_iter * thin;
    K=1;//da cambiare
    Mstar=3;//Mstar inizializzata dopo
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
        //Rcpp::Rcout<<n_j[j]<<std::endl;

    }
    for (int k = 0; k < d; ++k) {
        gamma.push_back(1.0);
        U.push_back(0.0);
    }
    for (int l = 0; l <K ; ++l) {
        N_k.push_back(0);

    }
    //std::vector<double> gamma(d,1.0);
    //std::vector<double> U(d,0.0);


    initialize_S(M);
    initialize_tau(M);
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

void GS_data::initialize_S(unsigned int M){
    //sample::rgamma rgamma;
    S = GDFMM_Traits::MatRow(d, M);
    for (int i = 0; i <d ; ++i) {
        for (int j = 0; j <M ; ++j) {
            S(i,j)=1;

        }

    }


}
void GS_data::initialize_N(unsigned int K){
    N = GDFMM_Traits::MatUnsCol(d, K);
    for (int i = 0; i <d ; ++i) {
        for (int j = 0; j <K ; ++j) {
            N(i,j)=1;

        }

    }



}
void GS_data::initialize_tau(unsigned int M){
  mu = std::vector<double>(M, 0.0);
  sigma = std::vector<double>(M, 0.0);
}
