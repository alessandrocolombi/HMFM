#include "GS_data.h"

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

void GS_data::update_log_sum(){
    log_sum = 0.0;
    for(size_t j=0; j<d; j++){
        log_sum += log(U[j]+1)*gamma[j];
    }
    // AL POSTO DEL FOR: log_sum = log( U.array() + 1).dot(gamma);
}
