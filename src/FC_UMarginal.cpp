#include <Rcpp.h>
#include <RcppEigen.h>

#include "FC_UMarginal.h"

FC_UMarginal::FC_UMarginal( std::string _na, bool _keepfixed,  
                            double _h1, double _h2, double _pow, unsigned int _d, double _adapt_var0):
                            FC_U(_na,_keepfixed),hyp1(_h1),hyp2(_h2),power(_pow)
                            {
                                adapt_var_proposal_U.resize(_d);
                                std::fill(adapt_var_proposal_U.begin(), adapt_var_proposal_U.end(), _adapt_var0);
                            };

void FC_UMarginal::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) {
    // Samplers
    sample::rnorm rnorm;
    sample::runif runif;

    // Data from GS_data
    std::vector<double>& U = gs_data.U; //U latent variables
    const std::vector<double>& gamma = gs_data.gamma; //gamma variables
    const unsigned int & d = gs_data.d; //number of levels
    const double& Lambda = gs_data.lambda; 
    const double& K = gs_data.K; // number of clusters
    const double& iter = gs_data.iterations; //iteration number
    const GDFMM_Traits::MatUnsCol& N = gs_data.N; // dxK matrix; 

    // Auxiliary variables
    double log_U_new; //log of proposed new value for U_j. U_new_j = exp(log_U_new)
    double U_new;     //proposed new value for U_j
    double ln_acp;    //log of acceptance probability
    double ln_u;      //log of uniform random variable used to decide if accept the move or not
    
    // Variables for Adaptive MH
    double ww_g{ pow(iter + 1, -hyp2) };
    
    // Rcpp::Rcout<<"iter="<<iter<<std::endl;
    
    for (unsigned int j = 0; j < d; j++){
        //Rcpp::Rcout<<"U:"<<U[j]<<std::endl;
        //cpp::Rcout<<"Adaptive variance var["<<j<<"] = "<<adapt_var_proposal_U[j]<<std::endl;

        // Update of U_j via Adapting Metropolis Hastings - computation of quantities is in logarithm for numerical reasons
        
        //1) Sample proposed value in log scale
        log_U_new = rnorm(gs_engine, std::log(U[j]), std::sqrt(adapt_var_proposal_U[j]));
        U_new = std::exp(log_U_new);
        
        //2) Compute acceptance probability in log scale
        ln_acp = log_FCU_marginal(U_new, Lambda, K, gamma[j], N.row(j).sum() ) - 
                 log_FCU_marginal(U[j],  Lambda, K, gamma[j], N.row(j).sum() ) +
                 log_U_new  - std::log(U[j]);

        //3) Acceptance rejection step
        ln_u = std::log(runif(gs_engine));
        
        if (ln_u  < ln_acp)
            gs_data.U[j] = U_new;

        //std::string update_status = (ln_u  < ln_acp)? " updated" : "NOT updated";
        //Rcpp::Rcout << "U_" << j << gs_data.U[j] << update_status << std::endl;

        adapt_var_proposal_U[j] *=  std::exp(  ww_g *( std::exp(std::min(0.0, ln_acp)) - hyp1 
                                                    )  
                                            );

        if (adapt_var_proposal_U[j] < 1/pow(10, power)){
            adapt_var_proposal_U[j] = 1/pow(10, power);
        }
        if(adapt_var_proposal_U[j] > pow(10,power)){
            adapt_var_proposal_U[j] = pow(10,power);
        }

    }
    
    //gs_data.update_log_sum();
    // Rcpp::Rcout<< "New log_sum : " << gs_data.log_sum <<std::endl;
}

double FC_UMarginal::log_FCU_marginal(const double& x, const double& Lambda, const unsigned int& K, const double& gamma_j, const unsigned int& n_j)const
{
    const double oneplusU{1+x};
    const double oneplusU_at_gamma{std::pow(oneplusU,gamma_j)};

    return(  (n_j - 1)*std::log(x) + Lambda/oneplusU_at_gamma + std::log( K + Lambda/oneplusU_at_gamma ) - oneplusU*(n_j + K*gamma_j)  );
}