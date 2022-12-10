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
    // Note: U_new_j = exp(log_U_new[j])
    std::vector<double> log_U_new(d, 1.0); //log of proposed new values for U_1,...,U_d. 
    std::vector<double> U_new(d, 1.0);     //proposed new value for U_1,...,U_d
    double ln_acp;    //log of acceptance probability
    double ln_u;      //log of uniform random variable used to decide if accept the move or not
    
    // Variables for Adaptive MH
    double ww_g{ pow(iter + 1, -hyp2) };
    
    // Rcpp::Rcout<<"iter="<<iter<<std::endl;
    
    /*
    // THIS IS THE OLD IMPLEMENTATION THAT IS WRONG
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
    */
    if(d == 1){
        //Rcpp::Rcout<<"U:"<<U[0]<<std::endl;
        //cpp::Rcout<<"Adaptive variance var["<<0<<"] = "<<adapt_var_proposal_U[0]<<std::endl;

        // Update of U_j via Adapting Metropolis Hastings - computation of quantities is in logarithm for numerical reasons
        
        //1) Sample proposed value in log scale
        log_U_new[0] = rnorm(gs_engine, std::log(U[0]), std::sqrt(adapt_var_proposal_U[0]));
        U_new[0] = std::exp(log_U_new[0]);
        
        //2) Compute acceptance probability in log scale
        /* 
        ln_acp = log_FCU_marginal(U_new[0], Lambda, K, gamma[0], N.row(0).sum() ) - 
                 log_FCU_marginal(U[0],  Lambda, K, gamma[0], N.row(0).sum() ) +
                 log_U_new[0]  - std::log(U[0]);
        */
                 
        ln_acp = log_FCU_marginal(U_new, Lambda, K, gamma, N ) - 
                 log_FCU_marginal(U,     Lambda, K, gamma, N ) +
                 log_U_new[0]  - std::log(U[0]);
        

        //3) Acceptance rejection step
        ln_u = std::log(runif(gs_engine));
        
        if (ln_u  < ln_acp)
            gs_data.U = U_new;

        //std::string update_status = (ln_u  < ln_acp)? " updated" : "NOT updated";
        //Rcpp::Rcout << "U_" << 0 << gs_data.U[0] << update_status << std::endl;

        adapt_var_proposal_U[0] *=  std::exp(  ww_g *( std::exp(std::min(0.0, ln_acp)) - hyp1 
                                                    )  
                                            );

        if (adapt_var_proposal_U[0] < 1/pow(10, power)){
            adapt_var_proposal_U[0] = 1/pow(10, power);
        }
        if(adapt_var_proposal_U[0] > pow(10,power)){
            adapt_var_proposal_U[0] = pow(10,power);
        }
    }
    else
        throw std::runtime_error("Error in FC_UMarginal: d>1 case not yet implemented ");


    // Update log_sum = sum_j( gamma_j * log(1+U_j) )
    gs_data.update_log_sum();
    
    //gs_data.update_log_sum();
    // Rcpp::Rcout<< "New log_sum : " << gs_data.log_sum <<std::endl;
}



// OLD FUNCTION THAT IS WRONG
/*
double FC_UMarginal::log_FCU_marginal(const double& x, const double& Lambda, const unsigned int& K, const double& gamma_j, const unsigned int& n_j)const
{
    const double oneplusU{1.0+x};
    const double oneplusU_at_gamma{std::pow(oneplusU,gamma_j)};

    return(  (n_j - 1)*std::log(x) + Lambda/oneplusU_at_gamma + std::log( K + Lambda/oneplusU_at_gamma ) - std::log(oneplusU)*(n_j + (double)K*gamma_j)  );
}
*/



// \pi(U_1, ... , U_d | rest ) \propto \exp{ \Lambda * \prod_j( \psi_j(U_j) ) }(K +  \Lambda * \prod_j( \psi_j(U_j) ))

double FC_UMarginal::log_FCU_marginal(const std::vector<double>& x, const double& Lambda, const unsigned int& K, const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N) const
{
    double res{0.0};
    
    double SumLog = 0.0;
    for(size_t j=0; j<x.size(); j++){
        double log_onePlusU{std::log(x[j] + 1.0)};
        SumLog += log_onePlusU*Gamma[j]; 

        res += ( (double)N.row(j).sum() - 1.0 )*std::log(x[j]) - 
               ( (double)N.row(j).sum() + (double)K*Gamma[j] )*log_onePlusU ; 
    }
    double Lambda_ProdPsi{Lambda*std::exp(-SumLog)};
    return( res +  Lambda_ProdPsi  +  std::log( (double)K + Lambda_ProdPsi  )  );

}

