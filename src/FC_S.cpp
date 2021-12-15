#include "FC_S.h"

void FC_S::update(GS_data& gs_data, sample::GSL_RNG gs_engine) {
    unsigned int d = gs_data.d;
    unsigned int K = gs_data.K;
    unsigned int Mstar = gs_data.Mstar;
    std::vector<double> U = gs_data.U;
    std::vector<double> gamma = gs_data.gamma;
    GDFMM_Traits::MatRow & S = gs_data.S;
    GDFMM_Traits::MatUnsCol N = gs_data.S;


    // Random sampler is created
    sample::rgamma Gamma;

    // Update routine
    for (j in 1:d) { //per ogni livello
        //S ALLOCATE
        for (k in 1:K) {//per ogni comp allocata
            S(j, k) = Gamma(1, N(j, k) + gamma[j], U[j] + 1);
        }

        //S NON ALLOCATE
        if (Mstar > 0) { // se c'è almeno una componente non allocata
            // S_na <- matrix(0, ncol = M_na, nrow = d)
            for (m_star in 1:Mstar) {
                S(j, m_star) = Gamma(1, gamma[j], U[j] + 1);
            }
        }
    }

}



}



