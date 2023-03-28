#include "FC_gamma.h"


FC_gamma::FC_gamma(std::string na, double h1, double h2, double pow, unsigned int d, double adapt_var0, double a, double b, double _s_p, bool _keepfixed) : 
                    FullConditional(na,_keepfixed), hyp1(h1), hyp2(h2), power(pow), alpha(a), beta(b), s_p(_s_p)
            {
                        //Rcpp::Rcout<<"Dentro al constructor di FC_gamma"<<std::endl;
                        //Rcpp::Rcout<<"a = "<<a<<std::endl;
                        //Rcpp::Rcout<<"b = "<<b<<std::endl;
                        //Rcpp::Rcout<<"alpha = "<<alpha<<std::endl;
                        //Rcpp::Rcout<<"beta = "<<beta<<std::endl;

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
    //double gamma_old;
    //double lmedia;
    //double ln_new;
    //double gamma_new;
    double ln_acp;
    double u_acc;
    double ww_g{ pow(iter + 1,-hyp2) };
    

    double log_gammaj_new{0.0};
    double gammaj_new{0.0};
    // Rcpp::Rcout<<"iter="<<iter<<std::endl;

    
    for (unsigned int j=0; j<d; j++){
                //Rcpp::Rcout<<"----------------------------------"<<std::endl;
                //Rcpp::Rcout<<"Dentro Adaptive for gamma, j = "<<j<<std::endl;
                //Rcpp::Rcout<<"alpha = "<<alpha<<std::endl;
                //Rcpp::Rcout<<"beta = "<<beta<<std::endl;
                //Rcpp::Rcout<<"gamma[j] = "<<gamma[j]<<std::endl;

        // ---------------------------------------------------
        // ADAPTIVE MH UPDATE for positive valued random vector
        // ---------------------------------------------------
        // Update of gamma_j via Adapting Metropolis Hastings - computation of quantities is in logarithm for numerical reasons
        //1) Sample proposed value in log scale
        log_gammaj_new = rnorm(gs_engine, std::log(gamma[j]), std::sqrt(adapt_var_pop_gamma[j]));
        gammaj_new = std::exp(log_gammaj_new);
                //Rcpp::Rcout<<"log_gammaj_new = "<<log_gammaj_new<<std::endl;
                //Rcpp::Rcout<<"gammaj_new = "<<gammaj_new<<std::endl;

        //2) Compute acceptance probability in log scale
        ln_acp = log_full_gamma(gammaj_new, Lambda, K, Mstar, N.row(j)) - 
                 log_full_gamma(gamma[j],   Lambda, K, Mstar, N.row(j)) +
                 log_gammaj_new  - std::log( gamma[j] );
                //Rcpp::Rcout<<"log_full new = "<<log_full_gamma(gammaj_new, Lambda, K, Mstar, N.row(j))<<std::endl;
                //Rcpp::Rcout<<"log_full old = "<<log_full_gamma(gamma[j], Lambda, K, Mstar, N.row(j))<<std::endl;
                //Rcpp::Rcout<<"log gamma new - log gamma old = "<<log_gammaj_new  - std::log( gamma[j] )<<std::endl;
                //Rcpp::Rcout<<"Prob accettazione = "<<std::exp(ln_acp)<<std::endl;

        // Checks
                 /*
        if(j == 0){
            double prior_ratio = (alpha - 1.0) * (log_gammaj_new  - std::log( gamma[j] )) - beta * (gammaj_new - gamma[j]);
            double likelihood_ratio = lgamma( gammaj_new * (double)(Mstar + K) ) - lgamma( gamma[j] * (double)(Mstar + K) ) +
                                    lgamma( gamma[j] * (double)(Mstar + K) + (double)N.row(j).sum() ) - lgamma( gammaj_new * (double)(Mstar + K) + (double)N.row(j).sum() ) +
                                    sumlgamma(gammaj_new, N.row(j)) - sumlgamma(gammaj_new, N.row(j)) + 
                                    K * (lgamma(gamma[j]) - lgamma(gammaj_new)) ;       
            double proposal_ratio = log_gammaj_new  - std::log( gamma[j] );
    
            Rcpp::Rcout<<"old value = "<<gamma[j]<<"; proposed value = "<<gammaj_new<<std::endl;
            Rcpp::Rcout<<"prior_ratio = "<<prior_ratio<<std::endl;
            Rcpp::Rcout<<"likelihood_ratio = "<<likelihood_ratio<<std::endl;
            Rcpp::Rcout<<"proposal_ratio = "<<proposal_ratio<<std::endl;
        }
        */
        //3) Acceptance rejection step
        u_acc = runif(gs_engine);
        
        if ( u_acc  < std::min(1.0, std::exp(ln_acp)) ){
                    //Rcpp::Rcout<<"Accetto"<<std::endl;
            gs_data.gamma[j] = gammaj_new;
        }

        adapt_var_pop_gamma[j] *=  std::exp(  ww_g *( std::exp(std::min(0.0, ln_acp)) - hyp1 
                                                    )  
                                            );
        
        //if(j == 0)
            //Rcpp::Rcout<<"adapt_var_pop_gamma[j]:"<<std::endl<<adapt_var_pop_gamma[j]<<std::endl;

        if (adapt_var_pop_gamma[j] < 1/pow(10, power)){
            adapt_var_pop_gamma[j] = 1/pow(10, power);
        }
        if(adapt_var_pop_gamma[j] > pow(10,power)){
            adapt_var_pop_gamma[j] = pow(10,power);
        }
        /*
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

        adapt_var_pop_gamma[j] *= std::exp(ww_g *(std::exp(std::min(0.0, ln_acp)) - hyp1));

        if (adapt_var_pop_gamma[j] < 1/pow(10, power)){
            adapt_var_pop_gamma[j] = 1/pow(10, power);
        }
        if(adapt_var_pop_gamma[j] > pow(10,power)){
            adapt_var_pop_gamma[j] = pow(10,power);
        }
        */

    }
}

double FC_gamma::log_full_gamma(double gamma, double Lambda, unsigned int k,
                                unsigned int M_star,const GDFMM_Traits::MatUnsCol & n_jk){
    if(n_jk.sum() == 0){
        Rcpp::Rcout<<"n_jk.sum() == 0"<<std::endl;
        throw std::runtime_error("Error in FC_gamma.cpp: all entries of n_jk are equal to 0 but at least one must be greater than 0 ");
    }
    if(n_jk.size() != k)
        throw std::runtime_error("Error in FC_gamma.cpp: n_jk should have length equal to k ");
    // Computation of the output
            //Rcpp::Rcout<<"********** Dentro log_full_gamma **********"<<std::endl;
            //Rcpp::Rcout<<"gamma = "<<gamma<<std::endl;
            //Rcpp::Rcout<<"Lambda = "<<Lambda<<std::endl;
            //Rcpp::Rcout<<"k = "<<k<<std::endl;
            //Rcpp::Rcout<<"M_star = "<<M_star<<std::endl;
            //Rcpp::Rcout<<"n_jk:"<<std::endl<<n_jk<<std::endl;
            //Rcpp::Rcout<<"l_dgamma("<<gamma<<","<<alpha<<", "<<beta<<" ) :"<<std::endl<<l_dgamma(gamma, alpha, beta)<<std::endl;
            //Rcpp::Rcout<<"lgamma(gamma * (M_star + k)):"<<std::endl<<lgamma(gamma * (M_star + k))<<std::endl;
            //Rcpp::Rcout<<"lgamma(gamma *(M_star + k) + n_jk.sum()):"<<std::endl<<lgamma(gamma *(M_star + k) + n_jk.sum())<<std::endl;
            //Rcpp::Rcout<<"(k * lgamma(gamma)):"<<std::endl<<(k * lgamma(gamma))<<std::endl;
            //Rcpp::Rcout<<"sumlgamma(gamma, n_jk):"<<std::endl<<sumlgamma(gamma, n_jk)<<std::endl;
    double out = (alpha - 1.0) * std::log(gamma) - 
                 beta * gamma +  //l_dgamma(gamma, alpha, beta) + 
                 lgamma( gamma * (double)(M_star + k) ) - 
                 lgamma( gamma * (double)(M_star + k) + (double)n_jk.sum() ) - 
                 k * lgamma(gamma) + 
                 sumlgamma(gamma, n_jk);
    //Rcpp::Rcout<<"std::log(pdfgamma(x,alpha,beta))"<< std::log(pdfgamma(x,alpha,beta));
    //Rcpp::Rcout<<"sumlgamma"<<sumlgamma(x, n_jk);
    //Rcpp::Rcout<<"********** Esco da log_full_gamma **********"<<std::endl;
    return out;
    }


double FC_gamma::sumlgamma(double gamma, const GDFMM_Traits::MatUnsCol& n_jk) {
    double sum = 0.0;
    for(unsigned i=0; i<n_jk.size(); i++){
        sum += lgamma( (double)n_jk(i) + gamma );
    }
    return sum;
}


// This is the log of a gamma density with parameters (a,b) evaluated in gamma
double FC_gamma::l_dgamma(double gamma, double a, double b){
            //Rcpp::Rcout<<"++++++++ Dentro l_dgamma ++++++++"<<std::endl;
            //Rcpp::Rcout<<"a*std::log(b) = "<<a*std::log(b)<<std::endl;
            //Rcpp::Rcout<<"(a-1.0)*std::log(gamma) = "<<(a-1.0)*std::log(gamma)<<std::endl;
            //Rcpp::Rcout<<"b*gamma = "<<b*gamma<<std::endl;
            //Rcpp::Rcout<<"std::lgamma(a) = "<<std::lgamma(a)<<std::endl;
            //Rcpp::Rcout<<"++++++++ Fuori l_dgamma ++++++++"<<std::endl;
    return a*std::log(b) + (a-1.0)*std::log(gamma) - b*gamma - std::lgamma(a);
}