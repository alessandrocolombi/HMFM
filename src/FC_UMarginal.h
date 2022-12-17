#ifndef GDFMM_FC_UMARGINAL_H
#define GDFMM_FC_UMARGINAL_H

#include "FullConditional.h"
#include "FC_U.h"
#include "recurrent_traits.h"

class FC_UMarginal : public FC_U{
private:
    /* MEMBERS */
    double hyp1 = 0.234; // adaptive tau_bar
    double hyp2 = 0.7;  
    double s_p  = 0.01; // mala parameter
    
    unsigned int power = 10;
    std::vector<double> adapt_var_proposal_U; //vector for variances to be adapted in MH steps
    
    // evaluate log_pi(U_1,...,U_d = x_1,...,x_d | rest )
    double log_FCU_marginal(const std::vector<double>& x, const double& Lambda, const unsigned int& K, 
                            const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N) const;
    // evaluate grad_log_pi(U_1,...,U_d = x_1,...,x_d | rest )
    std::vector<double> grad_log_FCU_marginal(const std::vector<double>& x, const double& Lambda, 
                                              const unsigned int& K, const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N) const;
public:
    //FC_UMarginal(bool _keepfixed):FC_U(_keepfixed){};
    FC_UMarginal(std::string _na, bool _keepfixed,  double _h1, double _h2, double _pow, unsigned int _d, double _adapt_var0, double _s_p);
    ~FC_UMarginal() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine);

};

#endif //GDFMM_FC_U_H
