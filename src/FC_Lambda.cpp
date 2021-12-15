#include "FC_Lambda.h"

void FC_Lambda::update(GS_data& gs_data, sample::GSL_RNG gs_engine){
    // From param all needed variable are retrived
    unsigned int k = gs_data.k;
    unsigned int d = gs_data.d;
    std::vector<double> gamma = gs_data.gamma; // CHIEDERE A COLMBI SE ABBIAMO d gamma UGUALI O d gamma[j]
    std::vector<double> U = gs_data.U;

    // Random sampler are created
    sample::runif Unif;
    sample::rgamma Gamma;
    
    // Update routine
    double a2_star = static_cast<double>( d*(k-1) ) + a2; 

    double log_sum = 0.0;

    for(size_t j=0; j<d; j++){
        log_sum += log(U[j]+1)*gamma[j];
    }
    // Computation of the weight for the "first" gamma distr.
    double p0 = (a2_star)/((a2_star-k)+k*(b2+1)*exp(log_sum));
    // Select via extraction from a uniform which distribution extract from
    bool select_p0 = binary_decision(p0, gs_engine);
    
    if(select_p0)
        gs_data.lambda = Gamma(gs_engine, a2_star + 1, b2 + 1 - exp(- log_sum));
    else
        gs_data.lambda = Gamma(gs_engine, a2_star, b2 + 1 - exp(- log_sum));
}