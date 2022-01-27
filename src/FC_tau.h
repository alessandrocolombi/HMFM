#ifndef GDFMM_FC_TAU_H
#define GDFMM_FC_TAU_H

#include "FullConditional.h"
#include "Partition.h"

    class FC_tau: public FullConditional {
    private:
        double nu_0; //= 2.5;
        double sigma_0; //= 40;//from data
        double mu_0; //= 57;//from data
        double k_0; // =8 mi sembra un po' sstrano

    public:

        double var(double mean,std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j, const std::vector<std::vector<double>>& data);
        double mean(std::vector<unsigned int> ind_i,std::vector<unsigned int>ind_j, const std::vector<std::vector<double>>& data);
        FC_tau(/* args */) {};
        FC_tau(std::string na, double nu0, double sigma0, double mu0, double k0):nu_0(nu0)
                                , sigma_0(sigma0), mu_0(mu0), k_0(k0){name = na;};
        ~FC_tau() {};
        void update(GS_data& gs_data, const sample::GSL_RNG& gs_engine) override;
    };
#endif
