#include "FC_Mstar.h"

void FC_Mstar::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    // From gs_data all needed variable are retrived
    const unsigned int& d = gs_data.d;
    const unsigned int& k = gs_data.K;
    const double& lambda = gs_data.lambda;
    const std::vector<double>& gamma = gs_data.gamma;
    const std::vector<double>& U = gs_data.U;
    const double& log_sum = gs_data.log_sum;

    // Random sampler for the Poisson Distribution is created
    sample::rpoisson Poisson;

    // UPDATE ROUTINE
    double p0 = lambda/( lambda + k*exp(-log_sum)); 
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
    gs_data.allocate_tau(gs_data.M);
}
