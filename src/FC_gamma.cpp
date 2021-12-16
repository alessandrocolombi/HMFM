//
// Created by ilari on 15/12/2021.
//

#include "FC_gamma.h"
void FC_gamma::update(GS_data& gs_data, sample::GSL_RNG gs_engine){
    sample::rnorm rnorm;
    sample::runif runif;
    sample::rgamma rgamma;
    std::vector<double> gamma=gs_data.gamma;
    unsigned int d = gs_data.d;
    double acc=0;
    double gamma_old;
    double lmedia;
    double ln_new;
    double gamma_new;
    double ln_acp;

    for (int j;j<d;j++) {
        gamma_old= gamma[j];
       lmedia = std::log(gamma_old);

        # Proposta
         ln_new = rnorm(lmedia, sqrt(adapt_var_pop_gamma));
         gamma_new = std::exp(ln_new);

        n_acp =log_full_gamma(gamma_new, Lambda,a1, b1 = b1, k = K, M_na = Mna, n_jk = N[j, ]) - lmedia;
        ln_acp = ln_acp - (log_full_gamma(gamma_old,  Lambda,a1 = a1, b1 = b1, k = K, M_na = Mna,n_jk = N[j, ]) - ln_new);

        double ln_u= std::log(runif());
        if (ln_u  < ln_acp){
            gamma[j] = gamma_new;
            acc = acc + 1;
        } else {
            gamma[j] = gamma_old
        }

        double ww_g = (iter + 1) ^ (- hyp2)


        adapt_var_pop_gamma<- adapt_var_pop_gamma * std::exp(ww_g *(std::exp(std::min(0, ln_acp)) -hyp1))

        if (adapt_var_pop_gamma< 10 ^ (-10)){
            adapt_var_pop_gamma<- 10 ^ (-10)
        }
        if (adapt_var_pop_gamma > 10 ^ (10)){
            adapt_var_pop_gamma <- 10 ^ 10
        }

    }

    double  FC_gamma::log_full_gamma(double x, double Lambda, int a1, int b1, int k, double M_na, double n_jk) {
        double out =dgamma(x, a1, b1, log = T) + lgamma(x * (M_na + k)) -
             std::log(rgamma(x *(M_na + k) + std::sum(n_jk)) - k * std::log(gamma(x) +
             std::sum(lgamma(n_jk + x))
        return out
    }
}