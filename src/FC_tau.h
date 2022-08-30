#ifndef GDFMM_FC_TAU_H
#define GDFMM_FC_TAU_H

#include "FullConditional.h"
#include "Partition.h"

    class FC_tau: public FullConditional {
    private:
        double nu_0; 
        double sigma_0; 
        double mu_0; 
        double k_0; 

    public:

        double var(double mean, const std::vector<unsigned int>& ind_i, const std::vector<unsigned int>& ind_j, const std::vector<std::vector<double>>& data);
        double mean(const std::vector<unsigned int>& ind_i, const std::vector<unsigned int>& ind_j, const std::vector<std::vector<double>>& data);
        FC_tau(bool _keepfixed) :FullConditional("tau", _keepfixed){};
        FC_tau(std::string na, double nu0, double sigma0, double mu0, double k0, bool _keepfixed):FullConditional(na,_keepfixed),nu_0(nu0)
                                , sigma_0(sigma0), mu_0(mu0), k_0(k0) {};
        ~FC_tau() {};
        void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
    };
#endif
