#ifndef __GDFMM_FC_POSTERIORMARGINAL_HPP__
#define __GDFMM_FC_POSTERIORMARGINAL_HPP__

#include "FullConditional.h"
#include "Partition.h"
#include "recurrent_traits.h"
#include <Rcpp.h>

class FC_PartitionMarginal : public Partition{
private: 
    double nu_0; 
    double sigma_0; 
    double mu_0; 
    double k_0; 
    double dnct(const double& x, double const & n0, double const & mu0, const double& gamma0) const;
    double log_dnct(const double& x, double const & n0, double const & mu0, const double& gamma0) const;
public:	
	//FC_PartitionMarginal(){};
    FC_PartitionMarginal(std::string _na, const unsigned int _d, const std::vector<unsigned int>& _n_j, bool _FixPart, double _nu_0, double _sigma_0, double _mu_0, double _k_0 );
    ~FC_PartitionMarginal(){};

    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);
};

#endif