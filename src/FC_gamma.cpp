//
// Created by ilari on 15/12/2021.
//

#include "FC_gamma.h"
void FC_gamma::update(GS_data & gs_data, const sample::GSL_RNG & gs_engine){
    // Samplers

    sample::rnorm rnorm;
    sample::runif runif;
    //sample::rgamma rgamma;

    // Data from G_stat
    std::vector<double>& gamma = gs_data.gamma;
    const unsigned int & d = gs_data.d;
    double Lambda=gs_data.lambda;// Initialization of some variable
    double K=gs_data.K;
    double Mna = gs_data.Mstar;
    double iter = gs_data.iterations;
    GDFMM_Traits::MatUnsCol N=gs_data.N;
    double acc=0; // <------- DA RIVEDERE
    double gamma_old;
    double lmedia;
    double ln_new;
    double gamma_new;
    double ln_acp;
    double ln_u;
    double ww_g;
    //for (int i = 1; i < 100; ++i) {
        //for (auto j=0;j<d;j++){
       // Rcpp::Rcout<<"Gamma:"<<gamma[j]<<std::endl;
       // }
    //}
    for (auto j=0;j<d;j++){
        gamma_old= gamma[j];
        //Rcpp::Rcout<<"Gamma:"<<gamma[j]<<std::endl;
        lmedia = std::log(gamma_old);
      //cpp::Rcout<<"lmedia"<<lmedia<<std::endl;
       //cpp::Rcout<<"var"<<adapt_var_pop_gamma<<std::endl;
        // Update of Gamma via Adapting Metropolis Hastings
        // *computation of quantities is in logarithm for numerical reasons*
        ln_new = rnorm(gs_engine, lmedia, std::sqrt(adapt_var_pop_gamma));
        //Rcpp::Rcout<<"ln_new"<<ln_new<<std::endl;
        gamma_new = std::exp(ln_new);
       //Rcpp::Rcout<<gamma_new;
        ln_acp =log_full_gamma(gamma_new, Lambda, K, Mna, N.row(j)) - lmedia; //da rivedere il tipo
       //Rcpp::Rcout<<ln_acp<<"--"<<std::endl;
       //cpp::Rcout<<gamma_new<<"--"<<std::endl;
       //cpp::Rcout<<N.row(j)<<"--"<<std::endl;

         ln_acp = ln_acp - (log_full_gamma(gamma_old,  Lambda, K, Mna, N.row(j)) - ln_new);
        //Rcpp::Rcout<<ln_acp<<"--"<<std::endl;
        //Rcpp::Rcout<<log_full_gamma(gamma_old,  Lambda, K, Mna, N.row(j))<<std::endl;
       // Rcpp::Rcout<<"Step 1";
        ln_u= std::log(runif(gs_engine));
        if (ln_u  < ln_acp){
            gamma[j] = gamma_new;
            acc = acc + 1;
        } else {
            gamma[j] = gamma_old;
        }


        ww_g = pow(iter + 1,- hyp2);

       //Rcpp::Rcout<<"Step 2"<<std::endl;
       adapt_var_pop_gamma = adapt_var_pop_gamma *
                                     std::exp(ww_g *(std::exp(std::min(0.0, ln_acp)) -hyp1));

        if (adapt_var_pop_gamma < 0.01){
            adapt_var_pop_gamma = 0.01;

        } else {
            adapt_var_pop_gamma = pow(10,power);
        }

    }
}

double  FC_gamma::log_full_gamma(double x, double Lambda, int k, double M_na,const GDFMM_Traits::MatUnsCol & n_jk)const{
        // Gamma Distribution to compute new proposal
    sample::pdfgamma pdfgamma;
        // Computation of the output
        double out = std::log(pdfgamma(x,alpha,beta)) + lgamma(x * (M_na + k)) -lgamma(x *(M_na + k) + n_jk.sum()) -(k * lgamma(x)) + sumlgamma(x, n_jk);
    // Rcpp::Rcout<<"std::log(pdfgamma(x,alpha,beta))"<< std::log(pdfgamma(x,alpha,beta));
    //Rcpp::Rcout<<"sumlgamma"<<sumlgamma(x, n_jk);
        return out;
    }


double  FC_gamma::sumlgamma(double x, const GDFMM_Traits::MatUnsCol& n_jk) const{
  double sum=0;
  for(unsigned i=0; i<n_jk.size(); i++){
   sum+=lgamma(n_jk(i)+x);
  }
    //Rcpp::Rcout<<"Step slg";
  return sum;
}
