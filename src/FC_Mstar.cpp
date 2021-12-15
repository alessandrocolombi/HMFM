#include "FC_Mstar.h"

void FC_Mstar::update(GS_data& gs_data, sample::GSL_RNG gs_engine){
    // From gs_data all needed variable are retrived
    unsigned int d = gs_data.d;
    unsigned int k = gs_data.K;
    double lambda = gs_data.lambda;
    std::vector<double> gamma = gs_data.gamma; // CHIEDERE A COLMBI SE ABBIAMO d gamma UGUALI o d gamma[j]
    std::vector<double> U = gs_data.U;
    double log_sum = gs_data.log_sum;

    // Random sampler is created
    sample::rpoisson Poisson;

    // Update routine

    // Computation of the weight for the traslated Poisson ?? CHIEDERE A COLOMBI CHIARIMENTO 
    double p0 = lambda/( lambda + k*exp(-log_sum)); // CHIEDERE A COLMBI SE VA BENE COSI'
    // Select, via extraction from a uniform, which distribution sample from
    bool select_p0 = binary_decision(p0, gs_engine);

    if(select_p0)
        gs_data.Mstar = Poisson(gs_engine, lambda*exp(-log_sum)) + 1;
    else
        gs_data.Mstar = Poisson(gs_engine, lambda*exp(-log_sum));
    
    // Update M in the Gibbs Sampler
    gs_data.M = k + gs_data.Mstar;
    // Initialize S according to new M
    gs_data.initialize_S(gs_data.M);
    // Initialize tau according to new M
    gs_data.initialize_tau(gs_data.M);
}