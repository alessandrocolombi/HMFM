#include "FC_Tau_A.h"

void FC_Tau_A::update(GS_data& gs_data, sample::GSL_RNG gs_engine){
  // From gs_data all needed variable are retrived
  unsigned int k = gs_data.K; // number of cluster
  unsigned int d = gs_data.d; // number of group
  std::vector<unsigned int> n_j = gs_data.n_j;//number of observation per group
  unsigned int M = gs_data.M;
  Eigen::Matrix S = gs_data.S;
  std::vector<double> mu=gs_data.mu;
  std::vector<double> sigma=gs_data.sigma;
  std::vector<std::vector<unsigned int>> probs;
  double probs_max;

  // Random sampler is created
  sample::rgamma Gamma;

  // Update routine
  double a2_star = static_cast<double>( d*(k-1) ) + a2;

  // Computation of the weight for the "first" gamma distr.
  double p0 = (a2_star)/((a2_star-k)+k*(b2+1)*exp(log_sum));
  // Select, via extraction from a uniform, which distribution sample from
  bool select_p0 = binary_decision(p0, gs_engine);

  if(select_p0)
    gs_data.lambda = Gamma(gs_engine, a2_star + 1, b2 + 1 - exp(-log_sum));
  else
    gs_data.lambda = Gamma(gs_engine, a2_star, b2 + 1 - exp(-log_sum));
}
