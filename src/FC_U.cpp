//
// Created by dario on 15/12/2021.
//

#include "FC_U.h"

void FC_S::update(GS_data& gs_data, sample::GSL_RNG gs_engine) {

    unsigned int d = gs_data.d;
    std::vector<unsigned int> n_j = gs_data.n_j;
    GDFMM_Traits::MatRow S = gs_data.S;


    vector<double> T = S.rowwise().sum();  //CONTROLLARE SE FA EFFETTIVAMENTE SOMMA PER RIGHE
    sample::rgamma Gamma;

    for (unsigned j=0; j<d; j++) { // for loop per livelli
        U[j] = Gamma(1, n_j[j], T[j]);
        gs_data.update_log_sum();
    }
}