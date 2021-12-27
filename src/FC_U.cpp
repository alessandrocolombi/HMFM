//
// Created by dario on 15/12/2021.
//

#include "FC_U.h"

void FC_U::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) {

    const unsigned int& d = gs_data.d;
    const std::vector<unsigned int>& n_j = gs_data.n_j;
    const GDFMM_Traits::MatRow& S = gs_data.S;
    // T_j is computed for each group (T_j = sum of S_ji over i for each group j)
    Eigen::Vector3d T = S.rowwise().sum();
    // Sampler for new U
    sample::rgamma Gamma;

    for (unsigned j=0; j<d; j++) { // for loop per livelli
        gs_data.U[j] = Gamma(gs_engine, n_j[j], T[j]);
    }
    gs_data.update_log_sum();
}
