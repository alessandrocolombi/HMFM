#include "FC_gamma.h"


FC_gamma::FC_gamma(std::string na, double h1, double h2, double pow, unsigned int d, double adapt_var0, int a, int b, double _s_p, bool _keepfixed) : 
                    FullConditional(na,_keepfixed), hyp1(h1), hyp2(h2), power(pow), alpha(a), beta(b), s_p(_s_p)
            {
                adapt_var_pop_gamma.resize(d);
                std::fill(adapt_var_pop_gamma.begin(), adapt_var_pop_gamma.end(), adapt_var0);
            };

void FC_gamma::update(GS_data & gs_data, const sample::GSL_RNG & gs_engine){
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
        ln_new = rnorm(gs_engine, lmedia, std::sqrt(adapt_var_pop_gamma[j])); //ln_new = rnorm(gs_engine, lmedia, std::sqrt(adapt_var_pop_gamma));
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
        ww_g = pow(iter + 1, -hyp2);

        /*
        adapt_var_pop_gamma = adapt_var_pop_gamma *
                                     std::exp(ww_g *(std::exp(std::min(0.0, ln_acp)) -hyp1));

        if (adapt_var_pop_gamma < 1/pow(10, power)){
            adapt_var_pop_gamma = 1/pow(10, power);
        }
        if(adapt_var_pop_gamma > pow(10,power)){
            adapt_var_pop_gamma = pow(10,power);
        }
        */
        adapt_var_pop_gamma[j] *= std::exp(ww_g *(std::exp(std::min(0.0, ln_acp)) - hyp1));

        if (adapt_var_pop_gamma[j] < 1/pow(10, power)){
            adapt_var_pop_gamma[j] = 1/pow(10, power);
        }
        if(adapt_var_pop_gamma[j] > pow(10,power)){
            adapt_var_pop_gamma[j] = pow(10,power);
        }

    }
}

double FC_gamma::log_full_gamma(double gamma, double Lambda, unsigned int k,
                unsigned int M_star,const GDFMM_Traits::MatUnsCol & n_jk){

    // Computation of the output
    double out = l_dgamma(gamma, alpha, beta) + lgamma(gamma * (M_star + k)) - 
                      lgamma(gamma *(M_star + k) + n_jk.sum()) -(k * lgamma(gamma)) + sumlgamma(gamma, n_jk);
    //Rcpp::Rcout<<"std::log(pdfgamma(x,alpha,beta))"<< std::log(pdfgamma(x,alpha,beta));
    //Rcpp::Rcout<<"sumlgamma"<<sumlgamma(x, n_jk);
    return out;
    }


double FC_gamma::sumlgamma(double gamma, const GDFMM_Traits::MatUnsCol& n_jk) {
    double sum = 0.0;
    for(unsigned i=0; i<n_jk.size(); i++){
        sum += lgamma(n_jk(i)+gamma);
    }
    //Rcpp::Rcout<<"Step slg";
    return sum;
}


// This is the log of a gamma density with parameters (a,b) evaluated in gamma
double FC_gamma::l_dgamma(double gamma, double a, double b){
    return a*std::log(b) + (a-1)*std::log(gamma) - b*gamma - std::lgamma(a);
}