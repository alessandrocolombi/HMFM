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
    double K=gs_data.K;
    double Mna = gs_data.Mstar;
    double iter = gs_data.iterations;
    GDFMM_Traits::MatUnsCol N=gs_data.N;
    double acc=0; // <------- DA RIVEDERE
    double gamma_old;
    double lmedia;
    double ln_new;
    double gamma_new;
    double ln_acp=0;
   // Rcpp::Rcout<<"Step 0";
    for (int j;j<d;j++){
        gamma_old= gamma[j];
       // Rcpp::Rcout<<"Step 1";
        lmedia = std::log(gamma_old);
        //Rcpp::Rcout<<"Step 1";
        // Update of Gamma via Adapting Metropolis Hastings
        // *computation of quantities is in logarithm for numerical reasons*
        double ln_new = rnorm(gs_engine, lmedia, sqrt(adapt_var_pop_gamma));
        //Rcpp::Rcout<<"Step 1";
        double gamma_new = std::exp(ln_new);
       // Rcpp::Rcout<<"Step 1";
        double n_acp =log_full_gamma(gamma_new, Lambda, K, Mna, N.row(j)) - lmedia; //da rivedere il tipo
       // Rcpp::Rcout<<"Step 1";
        // ARGOMENTI NON SI PARLANO CON HEADER FILE
        double ln_acp = ln_acp - (log_full_gamma(gamma_old,  Lambda, K, Mna, N.row(j)) - ln_new);
       // Rcpp::Rcout<<"Step 1";
        double ln_u= std::log(runif(gs_engine));
        if (ln_u  < ln_acp){
            gamma[j] = gamma_new;
            acc = acc + 1;
        } else {
            gamma[j] = gamma_old;
        }

        double ww_g = pow(iter + 1,- hyp2); // <--- DA RIVEDERE, SECONDO ME ABBIAMO UN ww_g DIVERSO

       // Rcpp::Rcout<<"Step 2";
        double adapt_var_pop_gamma = adapt_var_pop_gamma *
                                     std::exp(ww_g *(exp(std::min(0.0, ln_acp)) -hyp1));

        if (adapt_var_pop_gamma < pow(10,50)){  //<--- ANDRE:
            adapt_var_pop_gamma = pow(10,-50);  //<--- DA RIVEDERE SECONDO ME
        }                                       //<--- DA COME HO CAPITO IO QUI SI
        else{adapt_var_pop_gamma = pow(10,50);  //<--- DOVREBBE FARE UNA COSA DIVERSA
        }

    }
}

double  FC_gamma::log_full_gamma(double x, double Lambda, int k, double M_na, GDFMM_Traits::MatUnsCol n_jk)const{
        // Gamma Distribution to compute new proposal
    sample::pdfgamma pdfgamma;
        // Computation of the output
        double out =  pdfgamma(x,alpha,beta) + lgamma(x * (M_na + k)) -
                lgamma(x *(M_na + k) + n_jk.sum()) -
                k * lgamma(x) + sumlgamma(x, n_jk);
    //Rcpp::Rcout<<"Step lfg";
        return out;
    }


double  FC_gamma::sumlgamma(double x, GDFMM_Traits::MatUnsCol n_jk) const{
  double sum=0;
  for(unsigned i=0; i<n_jk.size(); i++){
   sum+=lgamma(n_jk(i)+x);
  }
    //Rcpp::Rcout<<"Step slg";
  return sum;
}
