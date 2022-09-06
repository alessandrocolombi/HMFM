#include "FC_S.h"

void FC_S::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    const unsigned int& d = gs_data.d;
    const unsigned int& K = gs_data.K;
    const unsigned int& Mstar = gs_data.Mstar;
    //const double& nu = gs_data.nu;
    const std::vector<double>& U = gs_data.U;
    const std::vector<double>& gamma = gs_data.gamma;
    const GDFMM_Traits::MatUnsCol& N = gs_data.N;
    GDFMM_Traits::MatRow & S = gs_data.S; //non-const ref. because we modify it

    // Random sampler is created
    sample::rgamma Gamma;

    // Update routine
    // Rcpp::Rcout << "S: ";

    //Initialize S according to new M
    gs_data.allocate_S(gs_data.M);

    for (unsigned j=0; j<d; j++) { //per ogni livello
        //S ALLOCATE
        // Rcpp::Rcout << "[";
        for (unsigned k=0; k < K; k++) {//per ogni comp allocata

            S(j, k) = Gamma(gs_engine, N(j, k) + gamma[j], 1 /(U[j] + 1.0) ); //This is S' and U is U'
            /*
            NO ERRORE! LE U SALVATE SONO LE U', quindi questa estrazione è proprio sbagliata. 
            in pratica, cosi facendo si sta facendo uno scaling aggiuntivo
            S(j, k) = Gamma(gs_engine, N(j, k) + gamma[j], 1 /(U[j] + 1) );
            S(j, k) *= 1.0/nu; // compute S', S' is gamma(gamma_j, nu), S' = 1/nu * S
            */
            // Rcpp::Rcout << S(j,k)<< " ";
        }
        //S NON ALLOCATE
        if (Mstar > 0) { // se c'è almeno una componente non allocata
            for (unsigned mstar=0; mstar<Mstar; mstar++) {

                S(j, K + mstar) = Gamma(gs_engine, gamma[j], 1 /(U[j] + 1.0) ); //This is S' and U is U'
                /*
                Vecchio ordine
                S(j, K + mstar) = Gamma(gs_engine, gamma[j],  1 /(U[j] + 1) );
                S(j, K + mstar) *= 1.0/nu; // compute S', S' is gamma(gamma_j, nu), S' = 1/nu * S
                */
            }
        }
        // Rcpp::Rcout << "]";
    }

    // Rcpp::Rcout << std::endl;
}





