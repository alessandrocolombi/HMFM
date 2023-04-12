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
    const bool& UseData = gs_data.UseData; 

    // Random sampler is created
    sample::rgamma Gamma;

    // Update routine
    //Rcpp::Rcout<<"Dentro FC_S.update.cpp"<<std::endl;

    //Initialize S according to new M
    gs_data.allocate_S(gs_data.M);

    for (unsigned j=0; j<d; j++) { //per ogni livello
        //S ALLOCATE
        //Rcpp::Rcout<<"Allocate: K = "<<K<<std::endl;
        for (unsigned k=0; k < K; k++) {//per ogni comp allocata
                    //Rcpp::Rcout<<"N("<<j<<","<<k<<") = "<<N(j,k)<<std::endl;
                    //Rcpp::Rcout<<"gamma["<<j<<"] = "<<gamma[j]<<std::endl;
                    //Rcpp::Rcout<<"U["<<j<<"] = "<<U[j]<<std::endl;
            double shape_S{gamma[j]};
            //if(UseData){
                shape_S += (double)N(j,k);
            //}
            S(j, k) = Gamma(gs_engine, shape_S, 1 /(U[j] + 1.0) ); 
            //if(std::abs(S(j, k)) < 1e-12){
                //Rcpp::Rcout<<"S(j, k) :"<<std::endl<<S(j, k) <<std::endl;
                //Rcpp::Rcout<<"shape_S = "<<shape_S<<std::endl;
                //Rcpp::Rcout<<"U[j] = "<<U[j]<<std::endl;
                //throw std::runtime_error("Error in FC_S.cpp, S is too small ");
            //}
            // Rcpp::Rcout << S(j,k)<< " ";
        }
        //S NON ALLOCATE
        if (Mstar > 0) { // se c'è almeno una componente non allocata
            for (unsigned mstar=0; mstar<Mstar; mstar++) {
                S(j, K + mstar) = Gamma(gs_engine, gamma[j], 1 /(U[j] + 1.0) ); //This is S' and U is U'
            }
        }
        // Rcpp::Rcout << "]";
    }

    // Rcpp::Rcout << std::endl;
}





