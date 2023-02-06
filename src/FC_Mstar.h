#ifndef __FC_M_STAR_H__
#define __FC_M_STAR_H__

#include "FullConditional.h"

class FC_Mstar : public FullConditional{
private:
    bool Partition_fixed;
    unsigned int proposal;
    std::vector<int> support_proposal;
public:

    FC_Mstar(std::string na, unsigned int _proposalMstar, bool FixPart, bool _keepfixed);
    //FC_Mstar(unsigned int _proposalMstar, bool _keepfixed) : FullConditional("Mstar",_proposalMstar,_keepfixed){};
    ~FC_Mstar() {};
    void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
    double log_full_cond_Mstar( const unsigned int& m, const unsigned int& K, const double& Lambda, 
                                const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N)const;
    double log_prob_MH( const unsigned int& m, const unsigned int& m_new, 
                        const unsigned int& K, const double& Lambda,
                        const std::vector<double>& Gamma, const GDFMM_Traits::MatUnsCol& N)const;
};

#endif
