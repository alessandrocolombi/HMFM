#include "FC_LambdaMarginal.h"

FC_LambdaMarginal::FC_LambdaMarginal(   std::string _na, double _a, double _b, bool _keepfixed, 
                                        double _h1, double _h2, double _pow, double _adapt_var0   ) : 
                                        FC_Lambda(_na,_a,_b,_keepfixed), hyp1(_h1), hyp2(_h2), power(_pow), adapt_var_proposal_Lambda(_adapt_var0){};

void FC_LambdaMarginal::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    
    throw std::runtime_error("Error: THIS UPDATE IS WRONG AND SHOULD NOT BE USED ");
    // Samplers
    sample::rnorm rnorm;
    sample::runif runif;

    // Data from GS_data
    double& Lambda = gs_data.lambda; 
    const std::vector<double>& U = gs_data.U; //U latent variables
    const std::vector<double>& gamma = gs_data.gamma; //gamma variables
    const double& K = gs_data.K; // number of clusters
    const double& iter = gs_data.iterations; //iteration number

    // Auxiliary variables
    double log_Lambda_new; //log of proposed new value for Lambda. Lambda_new = exp(log_Lambda_new)
    double Lambda_new;     //proposed new value for U_j
    double ln_acp;         //log of acceptance probability
    double ln_u;           //log of uniform random variable used to decide if accept the move or not
    
    // Variables for Adaptive MH
    double ww_g{ pow(iter + 1, -hyp2) };
    
    // Rcpp::Rcout<<"iter="<<iter<<std::endl;
    
    // Update of Lambda via Adapting Metropolis Hastings - computation of quantities is in logarithm for numerical reasons
        
    //1) Sample proposed value in log scale
    log_Lambda_new = rnorm(gs_engine, std::log(Lambda), std::sqrt(adapt_var_proposal_Lambda));
    Lambda_new = std::exp(log_Lambda_new);
        
    //2) Compute acceptance probability in log scale
    //double log_FCLambda_marginal(const double& x, const std::vector<double>& U, const std::vector<double>& gamma, const unsigned int& K) const;
    ln_acp = log_FCLambda_marginal(Lambda_new, U, gamma, K ) - 
             log_FCLambda_marginal(Lambda, U, gamma, K ) +
             log_Lambda_new  - std::log(Lambda);

    //3) Acceptance rejection step
    ln_u = std::log(runif(gs_engine));
        
    if (ln_u  < ln_acp)
        gs_data.lambda = Lambda_new;

    //std::string update_status = (ln_u  < ln_acp)? " updated" : "NOT updated";
    //Rcpp::Rcout << "Lambda" << gs_data.lambda << update_status << std::endl;

    adapt_var_proposal_Lambda *=  std::exp(  ww_g *( std::exp(std::min(0.0, ln_acp)) - hyp1 
                                             )  
                                    );

    if (adapt_var_proposal_Lambda < 1/pow(10, power)){
        adapt_var_proposal_Lambda = 1/pow(10, power);
    }
    if(adapt_var_proposal_Lambda > pow(10,power)){
        adapt_var_proposal_Lambda = pow(10,power);
    }

}

double FC_LambdaMarginal::log_FCLambda_marginal(const double& x, const std::vector<double>& U, const std::vector<double>& gamma, const unsigned int& K) const 
{
    /*
    const double sumPsi { std::inner_product(U.cbegin(), U.cend(), gamma.cbegin(), 0.0, std::plus<>(),
                                            [](const double& U_j, const double& gamma_j){return (1/std::pow(1.0+U_j,gamma_j));}
                                            ) 
                        }; //sum_{j=1}^d(Psi_j(U_j)) = sum_{j=1}^d( 1/(1+U_j)^gamma_j ) 
    */
    double sumPsi{0.0};                    
    double sum_log_K_plus_LambdaPsi{0.0};                    
    for(std::size_t j = 0; j < U.size(); j++){
        const double Psi_j{ 1.0/std::pow(1.0+U[j],gamma[j]) };
        sumPsi += Psi_j;
        sum_log_K_plus_LambdaPsi += std::log(K + x*Psi_j);
    }
    return(     std::log(x)*(a2 + U.size()*( (double)K-1.0) - 1.0) - 
                x*(b2 + U.size() - sumPsi) + 
                sum_log_K_plus_LambdaPsi  
          );
}