#include "FC_Mstar.h"

void FC_Mstar::update(GS_data& gs_data, sample::GSL_RNG gs_engine){
    // From gs_data all needed variable are retrived
    unsigned int d = gs_data.d;
    unsigned int k = gs_data.k;
    double lambda = gs_data.lambda;
    std::vector<double> gamma = gs_data.gamma; // CHIEDERE A COLMBI SE ABBIAMO d gamma UGUALI o d gamma[j]
    std::vector<double> U = gs_data.U;

    // Random sampler is created
    sample::rpoisson Poisson;

    // Update routine
    double log_sum = 0.0;

    for(size_t j=0; j<d; j++){
        log_sum += log(U[j]+1)*gamma[j];
    }
    // Computation of the weight for the traslated Poisson ?? CHIEDERE A COLOMBI CHIARIMENTO 
    double p0 = lambda/( lambda + k*exp(-log_sum)); // CHIEDERE A COLMBI SE VA BENE COSI'
    // Select, via extraction from a uniform, which distribution sample from
    bool select_p0 = binary_decision(p0, gs_engine);

    if(select_p0)
        gs_data.Mstar = Poisson(gs_engine, lambda*exp(-log_sum)) + 1;
    else
        gs_data.Mstar = Poisson(gs_engine, lambda*exp(-log_sum));
    
    // Update M in the Gibbs Sampler
    gs_data.M = k + gs_data. Mstar;
}