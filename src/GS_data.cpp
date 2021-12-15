#include "GS_data.h"

void GS_data::initialize_S(unsigned int M){
    S = GDFMM_Traits::MatRow(d, M);
}
void GS_data::initialize_tau(unsigned int M){
    mu = std::vector<double>(M, 0.0);
    sigma = std::vector<double>(M, 0.0);
}