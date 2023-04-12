#include "FC_gammaMarginal.h"
void FC_gammaMarginal::update(GS_data & gs_data, const sample::GSL_RNG & gs_engine){
    // Samplers
    sample::rnorm rnorm;
    sample::runif runif;

    // Data from GS_data
    std::vector<double>& gamma = gs_data.gamma; //gamma variables
    const std::vector<double>& U = gs_data.U; //U latent variables
    const unsigned int & d = gs_data.d; //number of levels
    const double& Lambda = gs_data.lambda; 
    //const double& log_sum = gs_data.log_sum; 
    const unsigned int& K = gs_data.K; // number of clusters
    const unsigned int& iter = gs_data.iterations; //iteration number
    const GDFMM_Traits::MatUnsCol& N = gs_data.N; // dxK matrix; 

    // Auxiliary variables
    // Note: gamma_new_j = exp(log_gamma_new[j])
    std::vector<double> log_gamma_new(d, 1.0); //log of proposed new values for gamma_1, ... , gamma_d. 
    std::vector<double> gamma_new(d, 1.0); //proposed new value for gamma_j
    double ln_acp;    //log of acceptance probability
    double ln_u;      //log of uniform random variable used to decide if accept the move or not
    
    // Variables for Adaptive MH
    double ww_g{ pow(iter + 1,-hyp2) };
    
    // Rcpp::Rcout<<"iter="<<iter<<std::endl;
    /*
    // OLD IMPLEMENTATION THAT IS WRONG!!
    for (unsigned int j = 0; j < d; j++){
        //Rcpp::Rcout<<"Gamma:"<<gamma[j]<<std::endl;
        //Rcpp::Rcout<<"Adaptive variamce var["<<j<<"] = "<<adapt_var_pop_gamma[j]<<std::endl;

        // Update of Gamma via Adapting Metropolis Hastings - computation of quantities is in logarithm for numerical reasons
        
        //1) Sample proposed value in log scale
        log_gamma_new = rnorm(gs_engine, std::log(gamma[j]), std::sqrt(adapt_var_pop_gamma[j]));
        gamma_new = std::exp(log_gamma_new);
        
        //2) Compute acceptance probability in log scale
        ln_acp = log_FCgamma_marginal(gamma_new, Lambda, K, U[j], N.row(j)) - 
                 log_FCgamma_marginal(gamma[j],  Lambda, K, U[j], N.row(j)) +
                 log_gamma_new  - std::log(gamma[j]);

        //3) Acceptance rejection step
        ln_u = std::log(runif(gs_engine));
        
        if (ln_u  < ln_acp)
            gs_data.gamma[j] = gamma_new;

        //std::string update_status = (ln_u  < ln_acp)? " updated" : "NOT updated";
        //Rcpp::Rcout << "gamma_" << j << gs_data.gamma[j] << update_status << std::endl;

        adapt_var_pop_gamma[j] *=  std::exp(  ww_g *( std::exp(std::min(0.0, ln_acp)) - hyp1 
                                                    )  
                                            );

        if (adapt_var_pop_gamma[j] < 1/pow(10, power)){
            adapt_var_pop_gamma[j] = 1/pow(10, power);
        }
        if(adapt_var_pop_gamma[j] > pow(10,power)){
            adapt_var_pop_gamma[j] = pow(10,power);
        }

    }
    */
    if(d == 0){ 

        // ---------------------------------------------------
        // ADAPTIVE MH UPDATE for positive valued random vector
        // ---------------------------------------------------
        // Update of Gamma via Adapting Metropolis Hastings - computation of quantities is in logarithm for numerical reasons
        // This implementation is only for d=1 case. For sake of code homogeneity, I set an if condition that can never be true
        // so that also in the case d=1 the MALA update is performed.
        
        //1) Sample proposed value in log scale
        log_gamma_new[0] = rnorm(gs_engine, std::log(gamma[0]), std::sqrt(adapt_var_pop_gamma[0]));
        gamma_new[0] = std::exp(log_gamma_new[0]);
        
        //2) Compute acceptance probability in log scale
        ln_acp = log_FCgamma_marginal(gamma_new, Lambda, K, U, N ) - 
                 log_FCgamma_marginal(gamma,     Lambda, K, U, N ) +
                 log_gamma_new[0]  - std::log(gamma[0]);

        //3) Acceptance rejection step
        ln_u = std::log(runif(gs_engine));
        
        if (ln_u  < ln_acp)
            gs_data.gamma[0] = gamma_new[0];

        //std::string update_status = (ln_u  < ln_acp)? " updated" : "NOT updated";
        //Rcpp::Rcout << "gamma_" << 0 << gs_data.gamma[0] << update_status << std::endl;

        adapt_var_pop_gamma[0] *=  std::exp(  ww_g *( std::exp(std::min(0.0, ln_acp)) - hyp1 
                                                    )  
                                            );

        if (adapt_var_pop_gamma[0] < 1/pow(10, power)){
            adapt_var_pop_gamma[0] = 1/pow(10, power);
        }
        if(adapt_var_pop_gamma[0] > pow(10,power)){
            adapt_var_pop_gamma[0] = pow(10,power);
        }
    }
    else{
        // ---------------------------------------------
        // MALA UPDATE for positive valued random vector
        // ---------------------------------------------
        // NB: This is ONLY if the matrix in the proposal is the identity. Otherwise the sampling from the multivariate normal must be different.
        //     Also the proposal ratio would be different

        // Compute log of current gamma vector
        std::vector<double> log_gamma(d,0.0); 
        std::transform(gamma.begin(),gamma.end(),log_gamma.begin(), 
                        [](const double& gamma_j){return(std::log(gamma_j));}
                      ); // compute log of each element

        // Compute MALA proposal values
        // NB: This is ONLY if the matrix in the proposal is the identity. Otherwise the sampling from the multivariate normal must be different 
        std::vector<double> mala_mean = grad_log_FCgamma_marginal(gamma, Lambda, K, U, N); // this is log_pi(gamma|rest) evaluated at gamma = exp(log_gamma)
        for(size_t j=0; j < d; j++){
            
            // compute mean
            mala_mean[j] = log_gamma[j] + 0.5*s_p*(mala_mean[j] * gamma[j] + 1.0); // adjust for mapping in the gradient and compute mala mean

            // draw value and compute its exponential
            log_gamma_new[j] = rnorm( gs_engine, mala_mean[j], std::sqrt(s_p) );
            gamma_new[j]     = std::exp(log_gamma_new[j]);
        }

        // MALA proposal is not symmetric. I also need to compute the mean of the inverse move
        const double& sp_temp = s_p;
        std::vector<double> mala_mean_invmove = grad_log_FCgamma_marginal(gamma_new, Lambda, K, U, N); // this is log_pi(gamma|rest) evaluated at gamma = exp(log_gamma_new)
        std::transform(mala_mean_invmove.begin(), mala_mean_invmove.end(), log_gamma_new.cbegin(),mala_mean_invmove.begin(),
                        [&sp_temp](const double& x_grad, const double& log_x){
                                    return(  log_x + 0.5*sp_temp*( x_grad * std::exp(log_x) + 1.0 )  );});
        
        // Compute log acceptance probability
        ln_acp = log_FCgamma_marginal(gamma_new, Lambda, K, U, N ) - // target evaluated in exp( log_gamma_new ) = gamma_new
                 log_FCgamma_marginal(gamma,     Lambda, K, U, N ) + // target evaluated in exp( log_gamma )     = gamma
                 std::accumulate(log_gamma_new.cbegin(),log_gamma_new.cend(),0.0)  - // sum of log_gamma_new
                 std::accumulate(log_gamma.cbegin(),log_gamma.cend(),0.0)  -         // sum of log_gamma
                 1.0/(2*s_p)*std::inner_product(log_gamma.cbegin(),log_gamma.cend(),mala_mean_invmove.cbegin(),0.0,std::plus<>(),
                                                 [](const double& x1, const double& x2){return( (x1-x2)*(x1-x2) );}) + //proposal of inverse move
                 1.0/(2*s_p)*std::inner_product(log_gamma_new.cbegin(),log_gamma_new.cend(),mala_mean.cbegin(),0.0,std::plus<>(),
                                                 [](const double& x1, const double& x2){return( (x1-x2)*(x1-x2) );;}); //proposal move
        
        //3) Acceptance rejection step
        ln_u = std::log(runif(gs_engine));
        
        if (ln_u  < ln_acp)
            gs_data.gamma = gamma_new;
        
    }
}

double FC_gammaMarginal::log_raising_factorial(const unsigned int& n, const double& a) const
{

    if(n==0)
        return 0.0;
    if(a<0)
        throw std::runtime_error("Error in my_log_raising_factorial, can not compute the raising factorial of a negative number in log scale"); 
    else if(a==0.0){
        return -std::numeric_limits<double>::infinity();
    }
    else{
        
        double val_max{std::log(a+n-1)};
        double res{1.0};
        if (n==1)
            return val_max;
        for(std::size_t i = 0; i <= n-2; ++i){
            res += std::log(a + (double)i) / val_max;
        }
        return val_max*res;

    }
}

/*
// OLD FUNCTION THAT IS WRONG
double FC_gammaMarginal::log_FCgamma_marginal(const double& x, const double& Lambda, const unsigned int& K, const double& U_j, const GDFMM_Traits::VecUnsRow& n_jk) const
{   
    const double oneplusU{1.0+U_j};
    const double Lambda_over_oneplusUatx{Lambda/std::pow(oneplusU,x)};
    //sum_{i=1}^K( log_raising_factorial(n_jk[i]), x )
    //const double sumPochammer{ std::accumulate(n_jk.cbegin(), n_jk.cend(), 0.0, [x](const unsigned int& n){return log_raising_factorial(n,x);})   }; //this can be done only from Eigen 3.4
    double sumPochammer{0.0};
    for(std::size_t i = 0; i < n_jk.size(); i++)
        sumPochammer += log_raising_factorial(n_jk[i], x);


    return( Lambda_over_oneplusUatx + 
            std::log( (double)K + Lambda_over_oneplusUatx ) - 
            (double)K * x*std::log(oneplusU) + 
            sumPochammer +
            (alpha - 1.0)*std::log(x) - 
            beta*x    
        );
}
*/

double FC_gammaMarginal::log_FCgamma_marginal(const std::vector<double>& x, const double& Lambda, const unsigned int& K,
                                              const std::vector<double>& U, const GDFMM_Traits::MatUnsCol& N ) const
{   

    double res{0.0};
    double SumLog = 0.0;

    for(size_t j=0; j<x.size(); j++){

        double gamma_log_onePlusU{x[j]*std::log(U[j] + 1.0)};
        SumLog += gamma_log_onePlusU; 

        //sum_{i=1}^K( log_raising_factorial(N(j,k), x[j] )
        double sumPochammer{0.0};
        for(std::size_t k = 0; k < N.cols(); k++)
            sumPochammer += log_raising_factorial(N(j,k), x[j]);

        double beta_gamma = beta * Lambda;
        res +=  (alpha - 1.0) * std::log(x[j]) - 
                beta_gamma * x[j] - 
                (double)K * gamma_log_onePlusU +
                sumPochammer ;
    }
    
    double Lambda_ProdPsi{Lambda*std::exp(-SumLog)};
    return ( res + std::log( (double)K + Lambda_ProdPsi ) + Lambda_ProdPsi );

}

// evaluate grad_log_pi(gamma_1,...,gamma_d = x_1,...,x_d | rest )
std::vector<double> FC_gammaMarginal::grad_log_FCgamma_marginal(const std::vector<double>& x, const double& Lambda, const unsigned int& K, 
                                                                const std::vector<double>& U, const GDFMM_Traits::MatUnsCol& N) const
{

    auto sum_diGamma = [](const double& gamma_j, const GDFMM_Traits::VecUnsCol& n_j)
    {
        double res = 0.0;
        for(size_t k=0; k<n_j.size(); k++){
            res += gsl_sf_psi( gamma_j + (double)n_j(k) );
        }
        res -= gsl_sf_psi(gamma_j);
        return res;
    };

    // Naive implementation with double for-loop
    unsigned int d = x.size(); // gradient size

    // 1) Compute prod_psi
    double SumLog = 0.0; 
    for(size_t j=0; j<x.size(); j++){
        double log_onePlusU{std::log(U[j] + 1.0)}; // this is log(1+U_j)
        SumLog += log_onePlusU*x[j]; 
    }
    double Lambda_ProdPsi{Lambda*std::exp(-SumLog)}; // this is Lambda*prod_psi

    // 2) Compute the gradient
    double beta_gamma = beta * Lambda;
    std::vector<double> grad_res(d,0.0); // inizialize the result
    for(size_t j=0; j<x.size(); j++){
        grad_res[j] = ( alpha - 1.0 )/(x[j]) - // this is (a_1-1)/(gamma_j)
                      beta_gamma - // this is b_1
                      (double)K * std::log(1 + U[j]) + // this is K*log(1+U_j)
                      sum_diGamma(x[j],N.row(j)) - // this is sum_k ( digamma(gamma_j+n_jk) - digamma(gamma_j) )
                      Lambda_ProdPsi * std::log(1.0 + U[j]) * ( 1.0 + 1.0/( (double)K + Lambda_ProdPsi ) );

    }

    return grad_res;
}