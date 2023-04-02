#ifndef __GDFMM_FC_POSTERIOR_NEAL3_MV_HPP__
#define __GDFMM_FC_POSTERIOR_NEAL3_MV_HPP__

#include "FullConditional.h"
#include "FC_Partition_mv.h"
#include "recurrent_traits.h"
#include <Rcpp.h>

class Partition_Neal3_mv : public Partition_mv{
private: 
    double nu_0; 
    double sigma_0; 
    double mu_0; 
    double k_0; 
    double log_I(unsigned int N_ji, double Ybar_star_ji, double Vstar_ji, double mu1, double k1, double nu1, double sigma1 ) const;
    
    std::tuple<double,double,double,unsigned int> 
    compute_cluster_summaries(  const std::vector<unsigned int>& ind_i, 
                                const std::vector<unsigned int>& ind_j, 
                                const std::vector<std::vector<Individual>>& data,
                                const GDFMM_Traits::MatRow& beta )const;
    std::tuple<double,double,double,unsigned int>  
    compute_cluster_summaries(  const std::vector<unsigned int>& ind_i, 
                                const std::vector<unsigned int>& ind_j, 
                                const std::vector<std::vector<Individual>>& data )const;
public:	
	//Partition_Neal3_mv(){};
    Partition_Neal3_mv(std::string _na, const unsigned int _d, const std::vector<unsigned int>& _n_j, bool _FixPart, double _nu_0, double _sigma_0, double _mu_0, double _k_0 );
    ~Partition_Neal3_mv(){};

    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);
};

#endif