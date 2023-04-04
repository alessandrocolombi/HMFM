#ifndef GDFMM_FULLCONDITIONAL_H
#define GDFMM_FULLCONDITIONAL_H


#include "include_headers.h"
#include "recurrent_traits.h"
#include "GS_data.h"
#include "GSL_wrappers.h"
#include <Rcpp.h>

class FullConditional{
public :
    FullConditional(const std::string& _name, bool _keepfixed) : name(_name), keep_fixed(_keepfixed){};
    FullConditional() = default;
    std::string name;
    bool keep_fixed;
    virtual void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) = 0;
    void print() const;
    bool binary_decision(double p1, const sample::GSL_RNG& engine) const;
    double log_raising_factorial(const unsigned int& n, const double& a)const;
    std::tuple<double,double,double,unsigned int> 
    compute_cluster_summaries(  const GDFMM_Traits::uset_uiui& current_indicies,
                                const std::vector<std::vector<Individual>>& data,
                                const GDFMM_Traits::MatRow& beta )const;
    std::tuple<double,double,double,unsigned int>  
    compute_cluster_summaries(  const GDFMM_Traits::uset_uiui& current_indicies,
                                const std::vector<std::vector<Individual>>& data )const;

    virtual ~FullConditional() = default;
};
#endif
