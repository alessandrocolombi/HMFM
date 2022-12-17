#include <Rcpp.h>
#include <RcppEigen.h>

#include "FC_UMarginal.h"

FC_UMarginal::FC_UMarginal( std::string _na, bool _keepfixed,  
                            double _h1, double _h2, double _pow, unsigned int _d, double _adapt_var0, double _s_p):
                            FC_U(_na,_keepfixed),hyp1(_h1),hyp2(_h2),power(_pow),s_p(_s_p)
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
    else{
        
        // ---------------------------------------------
        // MALA UPDATE for positive valued random vector
        // ---------------------------------------------
        // NB: This is ONLY if the matrix in the proposal is the identity. Otherwise the sampling from the multivariate normal must be different.
        //     Also the proposal ratio would be different

        // Compute log of current U vector
        std::vector<double> log_U(d,0.0); 
        std::transform(U.begin(),U.end(),log_U.begin(), [](const double& U_j){return(std::log(U_j));}); // compute log of each element
        
        // Compute MALA proposal values
        // NB: This is ONLY if the matrix in the proposal is the identity. Otherwise the sampling from the multivariate normal must be different 
        std::vector<double> mala_mean = grad_log_FCU_marginal(U, Lambda, K, gamma, N); // this is log_pi(U|rest) evaluated at U = exp(log(U))
        for(size_t j=0; j < d; j++){
            
            // compute mean
            mala_mean[j] = log_U[j] + 0.5*s_p*(mala_mean[j] * U[j] + 1.0); // adjust for mapping in the gradient and compute mala mean

            // draw value and compute its exponential
            log_U_new[j] = rnorm( gs_engine, mala_mean[j], std::sqrt(s_p) );
            U_new[j]     = std::exp(log_U_new[j]);
        }
        // MALA proposal is not symmetric. I also need to compute the mean of the inverse move
        double& sp_temp = s_p;
        std::vector<double> mala_mean_invmove = grad_log_FCU_marginal(U_new, Lambda, K, gamma, N); // this is log_pi(U|rest) evaluated at U = exp(log_U_new)
        std::transform(mala_mean_invmove.begin(), mala_mean_invmove.end(), log_U_new.cbegin(),mala_mean_invmove.begin(),
                        [&sp_temp](const double& x_grad, const double& log_x){
                                    return(  log_x + 0.5*sp_temp*( x_grad * std::exp(log_x) + 1.0 )  );});
        //for(size_t j=0; j < d; j++){
            //// compute mean of inverse move
            //mala_mean_invmove[j] = log_U_new[j] + 0.5*s_p*(mala_mean_invmove[j] * U_new[j] + 1.0); // adjust for mapping in the gradient and compute mala mean
        //}
        // Compute log acceptance probability
        ln_acp = log_FCU_marginal(U_new, Lambda, K, gamma, N ) - // target evaluated in exp( log_U_new ) = U_new
                 log_FCU_marginal(U,     Lambda, K, gamma, N ) + // target evaluated in exp( log_w )     = U
                 std::accumulate(log_U_new.cbegin(),log_U_new.cend(),0.0)  - // sum of log_U_new
                 std::accumulate(log_U.cbegin(),log_U.cend(),0.0)  -         // sum of log_U
                 1.0/(2*s_p)*std::inner_product(log_U.cbegin(),log_U.cend(),mala_mean_invmove.cbegin(),0.0,std::plus<>(),
                                                 [](const double& x1, const double& x2){return(x1-x2);}) + //proposal of inverse move
                 1.0/(2*s_p)*std::inner_product(log_U_new.cbegin(),log_U_new.cend(),mala_mean.cbegin(),0.0,std::plus<>(),
                                                 [](const double& x1, const double& x2){return(x1-x2);}); //proposal move
        

        //3) Acceptance rejection step
        ln_u = std::log(runif(gs_engine));
        
        if (ln_u  < ln_acp)
            gs_data.U = U_new;

        //throw std::runtime_error("Error in FC_UMarginal: d>1 case not yet implemented ");
    }


    // Update log_sum = sum_j( gamma_j * log(1+U_j) )
    gs_data.update_log_sum();
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
    double Lambda_ProdPsi{Lambda*std::exp(-SumLog)}; // this is Lambda*prod_psi
    return( res +  Lambda_ProdPsi  +  std::log( (double)K + Lambda_ProdPsi  )  );

}


std::vector<double> FC_UMarginal::grad_log_FCU_marginal(const std::vector<double>& x, const double& Lambda, 
                                                        const unsigned int& K, const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N) const
{
    // Naive implementation with double for-loop
    unsigned int d = x.size(); // gradient size

    // 1) Compute prod_psi
    double SumLog = 0.0; 
    for(size_t j=0; j<x.size(); j++){
        double log_onePlusU{std::log(x[j] + 1.0)}; // this is log(1+U_j)
        SumLog += log_onePlusU*Gamma[j]; 
    }
    double Lambda_ProdPsi{Lambda*std::exp(-SumLog)}; // this is Lambda*prod_psi

    // 2) Compute the gradient
    std::vector<double> grad_res(d,0.0); // inizialize the result
    for(size_t j=0; j<x.size(); j++){
        grad_res[j] = ( (double)N.row(j).sum() - 1.0 )/x[j] - // this is (n_j-1)/U_j
                      ( (double)N.row(j).sum() +(double)K*Gamma[j])/(1.0 + x[j] ) - // this is (n_j-Kgamma_j)/(1+U_j)
                      ( Lambda_ProdPsi*Gamma[j]/(1.0 + x[j]) )*( 1.0 + 1.0/((double)K + Lambda_ProdPsi)  ); // this is Lambda*gamma_j/(1+U_j)*prod_psi*(1 + 1/(K + Lambda*prod_psi))
    }

    return grad_res;
}
