#include "FC_gammaMarginal.h"
void FC_gammaMarginal::update(GS_data & gs_data, const sample::GSL_RNG & gs_engine){
    Rcpp::Rcout<<"Questo non e l'update di FC_gamma.CPP"<<std::endl;
/*
    // Samplers

    sample::rnorm rnorm;
    sample::runif runif;
    //sample::rgamma rgamma;

    // Data from G_stat
    std::vector<double> gamma = gs_data.gamma;
    const unsigned int & d = gs_data.d;
    double Lambda = gs_data.lambda;// Initialization of some variable
    double K = gs_data.K;
    double Mstar = gs_data.Mstar;
    double iter = gs_data.iterations;
    GDFMM_Traits::MatUnsCol N = gs_data.N;
    double acc=0; // <------- DA RIVEDERE
    double gamma_old;
    double lmedia;
    double ln_new;
    double gamma_new;
    double ln_acp;
    double ln_u;
    double ww_g;
    
    // Rcpp::Rcout<<"iter="<<iter<<std::endl;
    
    for (unsigned int j=0;j<d;j++){
        gamma_old= gamma[j];
        //Rcpp::Rcout<<"Gamma:"<<gamma[j]<<std::endl;
        lmedia = std::log(gamma_old);
        //cpp::Rcout<<"lmedia"<<lmedia<<std::endl;
        //cpp::Rcout<<"var"<<adapt_var_pop_gamma<<std::endl;
        // Update of Gamma via Adapting Metropolis Hastings
        // *computation of quantities is in logarithm for numerical reasons*
        ln_new = rnorm(gs_engine, lmedia, std::sqrt(adapt_var_pop_gamma));
        gamma_new = std::exp(ln_new);
        ln_acp = log_full_gamma(gamma_new, Lambda, K, Mstar, N.row(j)) - lmedia; //da rivedere il tipo
        ln_acp = ln_acp - (log_full_gamma(gamma_old, Lambda, K, Mstar, N.row(j)) - ln_new);
        ln_u= std::log(runif(gs_engine));
        
        if (ln_u  < ln_acp){
            gs_data.gamma[j] = gamma_new;
            acc = acc + 1;
        } else {
            gs_data.gamma[j] = gamma_old;
        }

        //std::string update_status = (ln_u  < ln_acp)? " updated" : "NOT updated";
        //Rcpp::Rcout << "gamma_" << j << gs_data.gamma[j] << update_status << std::endl;
        ww_g = pow(iter + 1,- hyp2);

        adapt_var_pop_gamma = adapt_var_pop_gamma *
                                     std::exp(ww_g *(std::exp(std::min(0.0, ln_acp)) -hyp1));

        if (adapt_var_pop_gamma < 1/pow(10, power)){
            adapt_var_pop_gamma = 1/pow(10, power);
        }
        if(adapt_var_pop_gamma > pow(10,power)){
            adapt_var_pop_gamma = pow(10,power);
        }

    }
    */
}
