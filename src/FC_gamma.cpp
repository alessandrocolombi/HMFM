//
// Created by ilari on 15/12/2021.
//

#include "FC_gamma.h"
void FC_gamma::update(GS_data & gs_data, const sample::GSL_RNG & gs_engine){
    // Samplers
    sample::rnorm rnorm;
    sample::runif runif;
    sample::rgamma rgamma;
    // Data from G_stat
    std::vector<double>& gamma = gs_data.gamma;
    const unsigned int & d = gs_data.d;
    // Initialization of some variable
    double acc=0; // <------- DA RIVEDERE
    double gamma_old;
    double lmedia;
    double ln_new;
    double gamma_new;
    double ln_acp;

    for (int j;j<d;j++){
        gamma_old= gamma[j];
        lmedia = std::log(gamma_old);

        // Update of Gamma via Adapting Metropolis Hastings
        // *computation of quantities is in logarithm for numerical reasons*
        double ln_new = rnorm(gs_engine, lmedia, sqrt(adapt_var_pop_gamma));  
        double gamma_new = std::exp(ln_new);
        
        double n_acp =log_full_gamma(gamma_new, Lambda,a1, b1 = b1, k = K, M_na = Mna, n_jk = N[j, ]) - lmedia;
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        // ARGOMENTI NON SI PARLANO CON HEADER FILE 
        double ln_acp = ln_acp - (log_full_gamma(gamma_old,  Lambda,a1 = a1, b1 = b1, k = K, M_na = Mna,n_jk = N[j, ]) - ln_new);

        double ln_u= std::log(runif(gs_engine));
        if (ln_u  < ln_acp){
            gamma[j] = gamma_new;
            acc = acc + 1;
        } else {
            gamma[j] = gamma_old;
        }

        double ww_g = pow(iter + 1,- hyp2); // <--- DA RIVEDERE, SECONDO ME ABBIAMO UN ww_g DIVERSO


        double adapt_var_pop_gamma = adapt_var_pop_gamma * 
                                     std::exp(ww_g *(std::exp(std::min(0, ln_acp)) -hyp1));
        
        if (adapt_var_pop_gamma < pow(10,50)){  //<--- ANDRE:
            adapt_var_pop_gamma = pow(10,-50);  //<--- DA RIVEDERE SECONDO ME
        }                                       //<--- DA COME HO CAPITO IO QUI SI
        else{adapt_var_pop_gamma = pow(10,50);  //<--- DOVREBBE FARE UNA COSA DIVERSA
        }

    }

    double  FC_gamma::log_full_gamma(double x, double Lambda, int a1, 
                                     int b1, int k, double M_na, double n_jk) const{
        // Gamma Distribution to compute new proposal
        std::gamma_distribution<double> distr(int alpha,int beta);
        // Computation of the output
        double out =  distr(x) + std::log(distr(x * (M_na + k)) -
                std::log(distr(x *(M_na + k) + std::sum(n_jk))) - k * std::log(distr(x)) +
                std::sum(std::log(d(n_jk + x)));
        return out
    }

}
