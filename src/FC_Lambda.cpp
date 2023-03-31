#include "FC_Lambda.h"

void FC_Lambda::update(GS_data& gs_data, const sample::GSL_RNG& gs_engine){
    // From gs_data all needed variable are retrived
    const unsigned int& k = gs_data.K;
    const unsigned int& d = gs_data.d;
    const double& log_sum = gs_data.log_sum;
    const std::vector<double>& gamma = gs_data.gamma;

    //Rcpp::Rcout<<"log_sum:"<<std::endl<<log_sum<<std::endl;
    double a_Lambda{a2};
    double b_Lambda{b2};

    // NUOVA PRIOR: L(gamma,Lambda) = L(gamma|Lambda)*L(Lambda)
    a_Lambda += d * a_gamma;
    b_Lambda += b_gamma * std::accumulate(gamma.cbegin(), gamma.cend(), 0.0);


    // Random sampler is created
    sample::rgamma Gamma;

    // UPDATE ROUTINE
    double a2_star = static_cast<double>( (k-1) ) + a_Lambda;
    //OLD: double a2_star = static_cast<double>( d*(k-1) ) + a2; // by mistake, I added that ^d in the peppf formula

    // Computation of the weight for the "first" gamma distr.
    double p0 = (a2_star)/( (a2_star-1.0) +
                            (double)k * (b_Lambda+1.0)*exp(log_sum)
                           );
    // Select, via extraction from a uniform, which distribution sample from
    bool select_p0 = binary_decision(p0, gs_engine);

    if(select_p0)
        gs_data.lambda = Gamma(gs_engine, a2_star + 1.0, 1.0 /(b_Lambda + 1.0 - exp(-log_sum)) );
    else
        gs_data.lambda = Gamma(gs_engine, a2_star, 1.0 /(b_Lambda + 1.0 - exp(-log_sum)) );
}
