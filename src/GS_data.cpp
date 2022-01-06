#include "GS_data.h"


GS_data::GS_data(Eigen::MatrixXd const &dat, unsigned int n_iter, unsigned int burnin, unsigned int thin) {
    d=data.size();
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
    };

    std::vector<unsigned int> n_j;
    for (unsigned int j = 0; j < dat.rows(); ++j) {

        for (unsigned int i = 0; i <dat.cols() ; ++i) {
            if(std::isnan(dat(j,i))){
                n_j.push_back(i);
            }
        }
    };

    std::vector<double> gamma(d,1.0);
    initialize_S(M);
    initialize_tau(M);
}

void GS_data::update_log_sum(){
    log_sum = 0.0;
    for(size_t j=0; j<d; j++){
        log_sum += log(U[j]+1)*gamma[j];
    }
    // AL POSTO DEL FOR: log_sum = log( U.array() + 1).dot(gamma);
}

void GS_data::initialize_S(unsigned int M){
  S = GDFMM_Traits::MatRow(d, M);
}
void GS_data::initialize_N(unsigned int K){
  N = GDFMM_Traits::MatUnsCol(d, K);
}
void GS_data::initialize_tau(unsigned int M){
  mu = std::vector<double>(M, 0.0);
  sigma = std::vector<double>(M, 0.0);
}
